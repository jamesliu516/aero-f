#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>


#include "Communicator.h"
#include "Domain.h"
#include "DistGeoState.h"
#include "DistVector.h"
#include "IoData.h"
#include "KirchhoffIntegrator.h"
#include "SubDomain.h"

#ifdef AEROACOUSTIC

#include "fftw3.h"
#include "gsl/gsl_sf.h"

#endif // AEROACOUSTIC

/////////////////
#define _UH_DEBUG_
/////////////////

#ifdef AEROACOUSTIC
#pragma message "Compiling KirchhoffIntegrator.C with Aeroacoustic support"
#endif // AEROACOUSTIC

KirchhoffIntegrator::KirchhoffIntegrator
(
  IoData &iod,
  Domain *domain
) :
  d_iod(iod), d_domain_p(domain), d_X_p((DistSVec<double, 3>*) 0),
  d_globalToNodeID(),
  d_SurfType(SPHERE), d_R(0.0)
{

#ifndef AEROACOUSTIC 
  domain->getCommunicator()->fprintf(stderr,"*** Error: AERO-F was compiled without Aeroacoustic capability."
                                            "  Rerun cmake with -DAEROACOUSTIC=ON to use this feature\n");
  exit(-1);
#endif
 
  d_X_p = new DistSVec<double, 3>(d_domain_p->getNodeDistInfo());
  
  DistVec<double> A(d_domain_p->getNodeDistInfo());
  
  //--- Get the mesh position
  DistGeoState geoState(d_iod, d_domain_p);
  // restart the geoState (positions of the mesh) 
  // At return X contains the latest position of the mesh.
  {
    char temp[1]; temp[0] = '\0';
    geoState.setup1(temp, d_X_p, &A);
  }
  
  if (d_iod.surfKI.d_surfaceType == KirchhoffData::SPHERICAL)
    d_SurfType = SPHERE;
  else if (d_iod.surfKI.d_surfaceType == KirchhoffData::CYLINDRICAL)
    d_SurfType = CYLINDER;
  
}


KirchhoffIntegrator::~KirchhoffIntegrator
(
)
{
  
  if (d_X_p)
    delete d_X_p;
  
}


void KirchhoffIntegrator::computeIntegral
(
)
{
#ifdef AEROACOUSTIC
  char prefix[192];
  sprintf(&prefix[0], "%s_", d_iod.input.strKPtraces);

  int *vecLen = new int[d_domain_p->getNumLocSub() + 1];
  vecLen[0] = 0;

  struct Data {
    int nid;
    double value;
  };

  //
  // Step 1 -- Read the pressure snapshots (snapshots in time)
  //
  
  int numSnapshots = 1;
  double Tf = 0.0;
  bool goodFile = true;
  
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    char filename[256];
    long pos = 0;
    sprintf(&filename[0], "%s%d", &prefix[0], iSub);
    ifstream pfile(&filename[0], ios::binary);
    pfile.seekg(pos, ios::beg);
    Data dtmp;
    if (!pfile.read((char*) &dtmp, sizeof(Data)))
    {
      fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
      goodFile = false;
      continue;
    }
    pos += sizeof(Data);
    if (iSub == 0)
      numSnapshots = dtmp.nid;
    vecLen[iSub+1] = (int) dtmp.value;
    pfile.close();
  }

  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    vecLen[iSub+1] += vecLen[iSub];

  std::vector<double> freqSamples(numSnapshots);
  std::vector<int> nodeID(vecLen[d_domain_p->getNumLocSub()]);

  std::vector< std::complex<double> > myBlock(2*numSnapshots*vecLen[d_domain_p->getNumLocSub()]);
  std::complex<double> *in = &myBlock[0];
  std::complex<double> *out = in + numSnapshots * vecLen[d_domain_p->getNumLocSub()];

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    char filename[256];
    long pos = 0;
    sprintf(&filename[0], "%s%d", &prefix[0], iSub);
    ifstream pfile(&filename[0], ios::binary);
    pfile.seekg(pos, ios::beg);
    Data dtmp;
    if (!pfile.read((char*) &dtmp, sizeof(Data)))
    {
      fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
    }
    pos += sizeof(Data);
    //---
    if (vecLen[iSub] == vecLen[iSub+1])
    {
      pfile.close();
      continue;
    }
    //---
    double ttt = 0.0;
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    Data *myData = new Data[myLen];
    for (int ii = 0; ii < numSnapshots; ++ii)
    {
      pfile.seekg(pos, ios::beg);
      if (!pfile.read((char*) &ttt, sizeof(double)))
      {
        fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
      }
      pos += sizeof(double);
      pfile.seekg(pos, ios::beg);
      pfile.read((char*) myData, sizeof(Data)*myLen);
      //---
      if (iSub == 0)
      {
        //
        // UH (08/2012)
        // Here we assume that the time step is constant
        //
        if (ii == 1)
        {
          Tf = ttt * numSnapshots;
        }
        freqSamples[ii] = (ii == 0) ? 0.0 : 2.0 * M_PI * ii / Tf;
        for (int jj = 0; jj < myLen; ++jj)
          nodeID[jj + vecLen[iSub]] = myData[jj].nid;
        // UH (08/2012)
        // myData[jj].nid is the local index of a node
        // It is local per subdomain (and per processor).
      }
      //---
      for (int jj = 0; jj < myLen; ++jj)
      {
        in[ii + (jj + vecLen[iSub])*numSnapshots] = std::complex<double>(myData[jj].value, 0.0);
      }
      pos += sizeof(Data) * myLen;
    } // for (int ii = 0; ii < numSnapshots; ++ii)

    free(myData);
    pfile.close();

    //
    //--- Step 2: Do the FFT
    //

    fftw_complex *tmpi = new fftw_complex[numSnapshots];
    fftw_complex *tmpo = new fftw_complex[numSnapshots];

    fftw_plan plan_forward = fftw_plan_dft_1d(numSnapshots, tmpi, tmpo, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int jj = 0; jj < myLen; ++jj)
    {
      memcpy(tmpi, in+(jj+vecLen[iSub])*numSnapshots, numSnapshots*sizeof(fftw_complex));
      fftw_execute(plan_forward);
      //
      // Scale the result
      //
      for (int ii = 0; ii < numSnapshots; ++ii)
      {
        tmpo[ii][0] = tmpo[ii][0] / numSnapshots;
        tmpo[ii][1] = tmpo[ii][1] / numSnapshots;
        out[jj + vecLen[iSub] * numSnapshots + ii * myLen] = std::complex<double>(tmpo[ii][0], tmpo[ii][1]);
      } // for (int ii = 0; ii < numSnapshots; ++ii)

    } // for (int jj = 0; jj < myLen; ++jj)

    delete[] tmpi;
    delete[] tmpo;

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)

  //
  // At this point, the time-snapshots have been converted
  // to frequency-snapshots (with nodal values for each frequency).
  //

  //
  // Create map of node numberings
  // (from global to the subdomain to a local numbering on the surface)
  //

  int totNumNodes = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    totNumNodes += (*d_X_p)(iSub).size();

  d_globalToNodeID.resize(totNumNodes);
  for (int jj = 0; jj < totNumNodes; ++jj)
    d_globalToNodeID[jj] = -1;

  totNumNodes = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    for (int jj = 0; jj < myLen; ++jj)
      d_globalToNodeID[totNumNodes + nodeID[jj + vecLen[iSub]]] = jj;
    totNumNodes += (*d_X_p)(iSub).size();
  }

  //
  // Compute the coefficients of p in series
  //
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the coefficients for the series of p ...\n";
  }

  int Nmax = 0;
  std::vector< std::complex<double> > coeff;
  
  switch (d_SurfType)
  {
    case SPHERE:
      convertToSHSeries(out, numSnapshots, vecLen, coeff, Nmax);
      break;
    default:
      std::cerr << "\n !!! The Kirchhoff Integral is not implemented for this surface !!!\n\n";
      break;
  }

  //
  // Get the coefficient for the normal derivative
  //
  
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the coefficients for dp/dn ...\n";
  }

  std::vector< std::complex<double> > coeff_dudn(coeff.size());
  switch (d_SurfType)
  {
    case SPHERE:
      getdpdnSHseries(&coeff[0], &coeff_dudn[0], &freqSamples[0], numSnapshots, Nmax);
      break;
    default:
      std::cerr << "\n !!! The Kirchhoff Integral is not implemented for this surface !!!\n\n";
      break;
  }

  //
  // Compute the Kirchhoff integral for all the frequencies 
  // and all the observation points.
  //
  
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the Kirchhoff integrals ...\n";
  }

  switch (d_SurfType)
  {
    case SPHERE:
      integrateOnSphere(out, vecLen, &coeff_dudn[0], Nmax, &freqSamples[0], numSnapshots);
      break;
    default:
      std::cerr << "\n !!! The Kirchhoff Integral is not implemented for this surface !!!\n\n";
      break;
  }


  //
  // Compute the far-field pattern for all the frequencies 
  //
  
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the far-field pattern ...\n";
  }
  
  switch (d_SurfType)
  {
    case SPHERE:
      ffpDataOnSphere(out, vecLen, &coeff_dudn[0], Nmax, &freqSamples[0], numSnapshots);
      break;
    default:
      std::cerr << "\n !!! The far-field pattern is not implemented for this surface !!!\n\n";
      break;
  }
  
#endif  // AEROACOUSTIC
}


void KirchhoffIntegrator::getL2NormSquare
(
 std::complex<double> *nodal,
 const int numFreq,
 const int *vecLen,
 double *l2Norm
 )
{

  
  // This function compute the L2-norm squared of the pressure 
  // on the surface.
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww);
  
  int sLen = numFreq;
  for (int ik = 0; ik < sLen; ++ik)
    l2Norm[ik] = 0.0;
  
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = nodal + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      if (faces[j].getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        x[jj] = XX[faces[j][jj]][0];
        y[jj] = XX[faces[j][jj]][1];
        z[jj] = XX[faces[j][jj]][2];
      }
      //
      double area = 0.0;
      area += pow((x[2] - x[0])*(y[1] - y[0]) - (y[2] - y[0])*(x[1] - x[0]), 2.0);
      area += pow((y[2] - y[0])*(z[1] - z[0]) - (z[2] - z[0])*(y[1] - y[0]), 2.0);
      area += pow((z[2] - z[0])*(x[1] - x[0]) - (x[2] - x[0])*(z[1] - z[0]), 2.0);
      area = 0.5 * sqrt(area);
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[faces[j][jj]]];
        }
      }
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xp = 0.0, yp = 0.0, zp = 0.0;
        xp = x[0]*xi[gp] + x[1]*xi[gp+ww.size()] + x[2]*xi[gp+2*ww.size()];
        yp = y[0]*xi[gp] + y[1]*xi[gp+ww.size()] + y[2]*xi[gp+2*ww.size()];
        zp = z[0]*xi[gp] + z[1]*xi[gp+ww.size()] + z[2]*xi[gp+2*ww.size()];
        //
        std::complex<double> pp;
        for (int ik = 0; ik < numFreq; ++ik)
        {
          pp = p[3*ik]*xi[gp];
          pp += p[1+3*ik]*xi[gp+ww.size()];
          pp += p[2+3*ik]*xi[gp+2*ww.size()];
          // 
#pragma omp critical
          {
            l2Norm[ik] += real(pp * conj(pp)) * ww[gp] * area;
          }
            
        } // for (int ik = 0; ik < numFreq; ++ik)
          
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(numFreq, l2Norm);

  
}


void KirchhoffIntegrator::convertToSHSeries
(
  std::complex<double> *nodal,
  const int numFreq, 
  const int *vecLen,
  std::vector< std::complex<double> > &series,
  int &nmax
)
{

  
  // This function computes the series coefficients for a P1 function
  // defined on the Kirchhoff surface (it should be a sphere).
  // The computation exploits the orthogonality of spherical harmonics.

  // --> UH (08/2012)
  // The next section is commented until the energy fraction is figured out.
  //
  //// 1. Compute the square of the L^2 norms
  //std::vector<double> l2NormSquare(numFreq);
  //this->getL2NormSquare(nodal, numFreq, vecLen, &l2NormSquare[0]);
    
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();

  std::vector<double> xi, ww;
  getQuadrature(xi, ww);

  nmax = d_iod.surfKI.d_nyquist;

  std::vector<double> yNorm((nmax+1)*(nmax+1));

  int sLen = numFreq * (nmax + 1) * (nmax + 1);
  series.resize(sLen);
  for (int ik = 0; ik < sLen; ++ik)
    series[ik] = std::complex<double>(0.0, 0.0);

  
#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif
  

  //
  // Get the radius of the surface
  //
  
  int count = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
        
    for (int j = 0; j < faces.size(); ++j)
    {
      if (faces[j].getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double r, t, phi;
      for (int jj = 0; jj < 3; ++jj)
      {
        cart2sph_x(XX[faces[j][jj]][0], XX[faces[j][jj]][1], XX[faces[j][jj]][2], r, t, phi);
        d_R += r;
        count += 1;
      }
      //
    }
    
  }
  d_R = d_R / count;
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(1, &d_R);
  d_R = d_R / MyCom_p->size();

  
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {

    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;

    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = nodal + numFreq * vecLen[iSub];

//////
#ifdef _UH_DEBUG_
    std::ofstream MyTrace("mytrace.txt");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      for (int i3 = 0; i3 < (*d_X_p)(iSub).size(); ++i3)
      {
        if (d_globalToNodeID[i3] < 0)
          continue;
        MyTrace << i3 << " " << d_globalToNodeID[i3] << " ";
        double ro = 0.0, to = 0.0, po = 0.0;
        cart2sph_x(XX[i3][0], XX[i3][1], XX[i3][2], ro, to, po);
        MyTrace << ro << " ";
        MyTrace << to << " ";
        MyTrace << po << " ";
        MyTrace << real(val[d_globalToNodeID[i3] + ik * myLen]) << " ";
        MyTrace << imag(val[d_globalToNodeID[i3] + ik * myLen]) << " ";
        MyTrace << std::endl;
      }
    }
    MyTrace.close();
#endif
//////

    for (int j = 0; j < faces.size(); ++j)
    {
      if (faces[j].getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        x[jj] = XX[faces[j][jj]][0];
        y[jj] = XX[faces[j][jj]][1];
        z[jj] = XX[faces[j][jj]][2];
      }
      double r[3], t[3], phi[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        cart2sph_x(x[jj], y[jj], z[jj], r[jj], t[jj], phi[jj]);
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[faces[j][jj]]];
        }
      }
      //
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {

        double xp = 0.0, yp = 0.0, zp = 0.0;
        xp = x[0]*xi[gp] + x[1]*xi[gp+ww.size()] + x[2]*xi[gp+2*ww.size()];
        yp = y[0]*xi[gp] + y[1]*xi[gp+ww.size()] + y[2]*xi[gp+2*ww.size()];
        zp = z[0]*xi[gp] + z[1]*xi[gp+ww.size()] + z[2]*xi[gp+2*ww.size()];
        double myrp = 0.0, mytp = 0.0, myphip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, mytp, myphip);
        double dMdr[3], dMds[3], tmpdot;
        tmpdot = xp * (x[1] - x[0]) + yp * (y[1] - y[0]) + zp * (z[1] - z[0]);
        dMdr[0] = d_R/myrp * (x[1] - x[0]);
        dMdr[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMdr[1] = d_R/myrp * (y[1] - y[0]);
        dMdr[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMdr[2] = d_R/myrp * (z[1] - z[0]);
        dMdr[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        tmpdot = xp * (x[2] - x[0]) + yp * (y[2] - y[0]) + zp * (z[2] - z[0]);
        dMds[0] = d_R/myrp * (x[2] - x[0]);
        dMds[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMds[1] = d_R/myrp * (y[2] - y[0]);
        dMds[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMds[2] = d_R/myrp * (z[2] - z[0]);
        dMds[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        double cross = 0.0;
        cross += std::pow(dMdr[1]*dMds[2] - dMdr[2]*dMds[1], 2.0);
        cross += std::pow(dMdr[2]*dMds[0] - dMdr[0]*dMds[2], 2.0);
        cross += std::pow(dMdr[0]*dMds[1] - dMdr[1]*dMds[0], 2.0);
        cross = sqrt(cross) * 0.5 / (d_R * d_R);
        
#ifdef _UH_DEBUG_
        surfaceCheck += d_R * d_R * cross * ww[gp];
#endif
        
        for (int nn = 0; nn <= nmax; ++nn)
        {
          std::vector< std::complex<double> > Yn;
          sphericalHarmonic(nn, mytp, myphip, Yn);
          //
          for (int mm = 0; mm < Yn.size(); ++mm)
          {
            yNorm[nn*nn + mm] += real( Yn[mm] * std::conj(Yn[mm]) * cross * ww[gp]) ;
            Yn[mm] = std::conj(Yn[mm]);
          }
          //
          std::complex<double> pp;
          for (int ik = 0; ik < numFreq; ++ik)
          {
            pp = p[3*ik]*xi[gp];
            pp += p[1+3*ik]*xi[gp+ww.size()];
            pp += p[2+3*ik]*xi[gp+2*ww.size()];
            //
#pragma omp critical
            {
              int shift = ik*(nmax+1)*(nmax+1) + nn*nn;
              for (int mm = 0; mm < Yn.size(); ++mm)
                series[shift + mm] += pp * Yn[mm] * cross * ww[gp];
            }

          } // for (int ik = 0; ik < numFreq; ++ik)

        } // for (int nn = 0; nn < nmax; ++nn)

      } // for (int gp = 0; gp < ww.size(); ++gp)

    } // for (int j = 0; j < faces.size(); ++j)

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)

  MyCom_p->globalSum(sLen, &series[0]);

#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  fprintf(stdout, " ... Spherical Surface = %4.3e (Error %4.3e)\n",
          surfaceCheck, std::abs(surfaceCheck - 4*M_PI*d_R*d_R)/(4*M_PI*d_R*d_R));
#endif
  
  //
  //--- Filter the coefficients
  //

  double errNorm = 0.0;
  for (int ik = 0; ik < yNorm.size(); ++ik)
  {
    errNorm = (errNorm > std::abs(yNorm[ik]-1.0)) ? errNorm : std::abs(yNorm[ik]-1.0);
  }
  errNorm *= 4.0;
  std::cout << " ... Filter level (relative) = " << scientific << errNorm << "\n";

  for (int ik = 0; ik < numFreq; ++ik)
  {
    //--- Maximum
    double maxEntry = 0.0;
    for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
    {
      if (std::abs(series[in + ik*(nmax+1)*(nmax+1)]) > maxEntry)
        maxEntry = std::abs(series[in + ik*(nmax+1)*(nmax+1)]);
    }
    //--- Filter
    for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
    {
      if (std::abs(series[in + ik*(nmax+1)*(nmax+1)]) < errNorm*maxEntry)
        series[in + ik*(nmax+1)*(nmax+1)] = std::complex<double>(0.0, 0.0);
    }
  }

////////////////////
#ifdef _UH_DEBUG_
  std::ofstream MyTrace("mytrace_2.txt");
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    SVec<double,3> XX = (*d_X_p)(iSub);
    for (int ik = 0; ik < numFreq; ++ik)
    {
      for (int i3 = 0; i3 < (*d_X_p)(iSub).size(); ++i3)
      {
        if (d_globalToNodeID[i3] < 0)
          continue;
        MyTrace << i3 << " " << d_globalToNodeID[i3] << " ";
        double ro = 0.0, to = 0.0, po = 0.0;
        cart2sph_x(XX[i3][0], XX[i3][1], XX[i3][2], ro, to, po);
        MyTrace << ro << " ";
        MyTrace << to << " ";
        MyTrace << po << " ";
        std::complex<double> pp = evaluateSHS(&series[0] + ik*(nmax+1)*(nmax+1), nmax, to, po);
        MyTrace << real(pp) << " " << imag(pp) << " ";
        MyTrace << std::endl;
      }
    }
  }
  MyTrace.close();

  std::ofstream MyFile("p_coeff.txt");
  for (int ik = 0; ik < numFreq; ++ik)
  {
    MyFile << ik << " 0.0 " << std::endl;
    MyFile << (nmax+1)*(nmax+1) << " 0 " << std::endl;
    for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
    {
      MyFile << scientific << real(series[in + ik*(nmax+1)*(nmax+1)]) << " ";
      MyFile << scientific << imag(series[in + ik*(nmax+1)*(nmax+1)]) << " ";
      MyFile << std::endl;
    }
  }
  MyFile.close();

#endif
///////////////////

}


void KirchhoffIntegrator::getdpdnSHseries
(
  std::complex<double> *coeff,
  std::complex<double> *coeff_dudn,
  const double *freqSamples,
  const int numFrequencies,
  const int nMax
) const
{
  
  int count = 0;
  for (int ik = 0; ik < numFrequencies; ++ik)
  {
    
    double kappa = freqSamples[ik];

    if ((ik == 0) || (kappa == 0.0))
    {
      for (int in = 0; in <= nMax; ++in)
      {
        double scalar = -(in + 1.0)/d_R;
        for (int jn = 0; jn <= 2*in; ++jn)
        {
          coeff_dudn[count] = coeff[count] * scalar;
          count += 1;
        } // for (int jn = 0; jn <= 2*in; ++jn)
      } // for (int in = 0; in <= nMax; ++in)
    }
    else
    {

      for (int in = 0; in <= nMax; ++in)
      {
        //
        double maxEntry = 0.0;
        for (int jn = 0; jn <= 2*in; ++jn)
          maxEntry = (maxEntry > std::abs(coeff[count+jn])) ? maxEntry : std::abs(coeff[count+jn]);
        if (maxEntry == 0.0)
        {
          count += 2*in+1;
          continue;
        }
        //
        std::complex<double> scalar = kappa * besselh_prime(in, kappa*d_R) / besselh(in, kappa*d_R);
        for (int jn = 0; jn <= 2*in; ++jn)
        {
          coeff_dudn[count] = coeff[count] * scalar;
          count += 1;
        } // for (int jn = 0; jn <= 2*in; ++jn)
        //
      } // for (int in = 0; in <= nMax; ++in)
    }
    
  } // for (int ik = 0; ik < numFrequencies; ++ik)
  

////////////////////
#ifdef _UH_DEBUG_
  std::ofstream MyFile("dpdn_coeff.txt");
  for (int ik = 0; ik < numFrequencies; ++ik)
  {
    MyFile << freqSamples[ik] << " 0.0 " << std::endl;
    MyFile << (nMax+1)*(nMax+1) << " 0 " << std::endl;
    for (int in = 0; in < (nMax+1)*(nMax+1); ++in)
    {
      MyFile << scientific << real(coeff_dudn[in + ik*(nMax+1)*(nMax+1)]) << " ";
      MyFile << scientific << imag(coeff_dudn[in + ik*(nMax+1)*(nMax+1)]) << " ";
      MyFile << std::endl;
    }
  }
  MyFile.close();
#endif
///////////////////

  
}


void KirchhoffIntegrator::integrateOnSphere
(
 std::complex<double> *pvalues,
 int *vecLen,
 std::complex<double> *dpdn_coeff,
 int nmax,
 const double *freqSamples,
 const int numFreq
)
{
  
  // This function computes the Kirchhoff integral.
  
  int numProbes = 0;
  for (int ii = 0; ii < Probes::MAXNODES; ++ii)
  {
    Probes& myProbes = d_iod.output.transient.probes;
    if (myProbes.myNodes[ii].locationX < -1.0e10)
      break;
    numProbes += 1;
  }

  std::vector<double> observations(3 * numProbes, 0.0);
  for (int ii = 0; ii < numProbes; ++ii)
  {
    Probes& myProbes = d_iod.output.transient.probes;
    //
    observations[3*ii] = myProbes.myNodes[ii].locationX;
    observations[3*ii+1] = myProbes.myNodes[ii].locationY;
    observations[3*ii+2] = myProbes.myNodes[ii].locationZ;
  }
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww);
  
  int intLen = numFreq * numProbes;
  std::vector< std::complex<double> > pnoise(intLen, std::complex<double>(0.0, 0.0));
  
#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = pvalues + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      if (faces[j].getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        x[jj] = XX[faces[j][jj]][0];
        y[jj] = XX[faces[j][jj]][1];
        z[jj] = XX[faces[j][jj]][2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[faces[j][jj]]];
        }
      }
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xp = 0.0, yp = 0.0, zp = 0.0;
        xp = x[0]*xi[gp] + x[1]*xi[gp+ww.size()] + x[2]*xi[gp+2*ww.size()];
        yp = y[0]*xi[gp] + y[1]*xi[gp+ww.size()] + y[2]*xi[gp+2*ww.size()];
        zp = z[0]*xi[gp] + z[1]*xi[gp+ww.size()] + z[2]*xi[gp+2*ww.size()];
        
        double myrp = 0.0, tp = 0.0, phip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, tp, phip);
        
        double dMdr[3], dMds[3], tmpdot;
        tmpdot = xp * (x[1] - x[0]) + yp * (y[1] - y[0]) + zp * (z[1] - z[0]);
        dMdr[0] = d_R/myrp * (x[1] - x[0]);
        dMdr[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMdr[1] = d_R/myrp * (y[1] - y[0]);
        dMdr[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMdr[2] = d_R/myrp * (z[1] - z[0]);
        dMdr[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        tmpdot = xp * (x[2] - x[0]) + yp * (y[2] - y[0]) + zp * (z[2] - z[0]);
        dMds[0] = d_R/myrp * (x[2] - x[0]);
        dMds[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMds[1] = d_R/myrp * (y[2] - y[0]);
        dMds[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMds[2] = d_R/myrp * (z[2] - z[0]);
        dMds[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        
        double cross = 0.0;
        cross += std::pow(dMdr[1]*dMds[2] - dMdr[2]*dMds[1], 2.0);
        cross += std::pow(dMdr[2]*dMds[0] - dMdr[0]*dMds[2], 2.0);
        cross += std::pow(dMdr[0]*dMds[1] - dMdr[1]*dMds[0], 2.0);
        cross = sqrt(cross) * 0.5 / (d_R * d_R);
        
#ifdef _UH_DEBUG_
        surfaceCheck += d_R * d_R * cross * ww[gp];
#endif
        //
        std::complex<double> pp;
        std::vector< std::complex<double> > dpdn(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          dpdn[ik] = evaluateSHS(dpdn_coeff + ik*(nmax+1)*(nmax+1), nmax, tp, phip);
        }
        //
        for (int id = 0; id < numProbes; ++id)
        {
          double rv[3];
          rv[0] = xp - observations[3*id];
          rv[1] = yp - observations[3*id+1];
          rv[2] = zp - observations[3*id+2];
          //
          double rho = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
          //
          double drho_dn = (rv[0] * xp + rv[1] * yp + rv[2] * zp ) / (myrp * rho);
          //
          for (int ik = 0; ik < numFreq; ++ik)
          {
            
            double kappa = freqSamples[ik];

            std::complex<double> Green = 0.25/(M_PI*rho) * exp(std::complex<double>(0.0, kappa*rho));
            std::complex<double> dGreendNu = Green * drho_dn * std::complex<double>(-1.0/rho, kappa);
            
            //
            pp = p[3*ik]*xi[gp];
            pp += p[1+3*ik]*xi[gp+ww.size()];
            pp += p[2+3*ik]*xi[gp+2*ww.size()];
            //
#pragma omp critical
            {
              pnoise[id + ik * numProbes] += (pp * dGreendNu - dpdn[ik] * Green) * cross * ww[gp] * d_R * d_R;
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)

        } // for (int id = 0; id < numProbes; ++id)
          
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(intLen, &pnoise[0]);
  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  std::cout << " Surface " << surfaceCheck << " " << std::abs(surfaceCheck - d_R * d_R * 4*M_PI)/(4*d_R * d_R * M_PI);
  std::cout << std::endl;
#endif

  if (MyCom_p->cpuNum() == 0)
  {
    
    char *myName = new char[strlen(d_iod.output.transient.probes.prefix) + strlen(d_iod.output.transient.probes.pressure) + 1];
    sprintf(myName, "%s%s", d_iod.output.transient.probes.prefix, d_iod.output.transient.probes.pressure);

    // Output PNOISE in ASCII file
    FILE *out = fopen(myName, "w");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      fprintf(out, "%d %e ", ik, freqSamples[ik]);
      for (int id = 0; id < numProbes; ++id)
      {
        fprintf(out, "%e ", std::real(pnoise[id + ik * numProbes]));
        fprintf(out, "%e ", std::imag(pnoise[id + ik * numProbes]));
        fprintf(out, "%e ", std::abs(pnoise[id + ik * numProbes]));
      }
      fprintf(out, "\n");
    }
    fclose(out);
    delete[] myName;
    
  } // if (MyCom_p->cpuNum() == 0)
  
  
}


void KirchhoffIntegrator::ffpDataOnSphere
(
 std::complex<double> *pvalues,
 int *vecLen,
 std::complex<double> *dpdn_coeff,
 const int nmax,
 const double *freqSamples,
 const int numFreq
 )
{
  
  // This function computes the far-field pattern for
  // (p, dp/dn) specified on a sphere.
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww);
  
  int numInc = d_iod.surfKI.d_angularIncrement;
  int numDir = (numInc/2 + 1) * numInc;
  
  std::vector<double> directions(3*numDir, 0.0);
  numDir = 0;
  for (int ib = 0; ib <= numInc/2; ++ib)
  {
    double beta = ib * M_PI * 2.0 / numInc - 0.5 * M_PI;
    double cb = cos(beta), sb = sin(beta);
    for (int ia = 0; ia < numInc; ++ia)
    {
      double alpha = ia * 2.0 * M_PI / numInc;
      directions[3*numDir] = cos(alpha) * cb;
      directions[3*numDir+1] = sin(alpha) * cb;
      directions[3*numDir+2] = sb;
      numDir += 1;
    }
  }
 
  int ffpLen = numDir * numFreq; 
  std::vector< std::complex<double> > ffp(ffpLen, std::complex<double>(0.0, 0.0));

#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif
  

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = pvalues + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      if (faces[j].getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        x[jj] = XX[faces[j][jj]][0];
        y[jj] = XX[faces[j][jj]][1];
        z[jj] = XX[faces[j][jj]][2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[faces[j][jj]]];
        }
      }
      //
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xp = 0.0, yp = 0.0, zp = 0.0;
        xp = x[0]*xi[gp] + x[1]*xi[gp+ww.size()] + x[2]*xi[gp+2*ww.size()];
        yp = y[0]*xi[gp] + y[1]*xi[gp+ww.size()] + y[2]*xi[gp+2*ww.size()];
        zp = z[0]*xi[gp] + z[1]*xi[gp+ww.size()] + z[2]*xi[gp+2*ww.size()];
        
        double myrp = 0.0, tp = 0.0, phip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, tp, phip);
        
        double dMdr[3], dMds[3], tmpdot;
        tmpdot = xp * (x[1] - x[0]) + yp * (y[1] - y[0]) + zp * (z[1] - z[0]);
        dMdr[0] = d_R/myrp * (x[1] - x[0]);
        dMdr[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMdr[1] = d_R/myrp * (y[1] - y[0]);
        dMdr[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMdr[2] = d_R/myrp * (z[1] - z[0]);
        dMdr[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        tmpdot = xp * (x[2] - x[0]) + yp * (y[2] - y[0]) + zp * (z[2] - z[0]);
        dMds[0] = d_R/myrp * (x[2] - x[0]);
        dMds[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
        dMds[1] = d_R/myrp * (y[2] - y[0]);
        dMds[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
        dMds[2] = d_R/myrp * (z[2] - z[0]);
        dMds[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
        
        double cross = 0.0;
        cross += std::pow(dMdr[1]*dMds[2] - dMdr[2]*dMds[1], 2.0);
        cross += std::pow(dMdr[2]*dMds[0] - dMdr[0]*dMds[2], 2.0);
        cross += std::pow(dMdr[0]*dMds[1] - dMdr[1]*dMds[0], 2.0);
        cross = sqrt(cross) * 0.5 / (d_R * d_R);

#ifdef _UH_DEBUG_
        surfaceCheck += d_R * d_R * cross * ww[gp];
#endif
        
       //
        std::complex<double> pp;
        std::vector< std::complex<double> > dpdn(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          dpdn[ik] = evaluateSHS(dpdn_coeff + ik*(nmax+1)*(nmax+1), nmax, tp, phip);
        }
        //
        for (int id = 0; id < numDir; ++id)
        {
          // 
          double DirDotX = xp * directions[3*id];
          DirDotX += yp * directions[3*id+1];
          DirDotX += zp * directions[3*id+2];
          //
          double DirDotNu = DirDotX / myrp;
          //
          for (int ik = 0; ik < numFreq; ++ik)
          {
            
            double kappa = freqSamples[ik];
            //
            pp = p[3*ik]*xi[gp];
            pp += p[1+3*ik]*xi[gp+ww.size()];
            pp += p[2+3*ik]*xi[gp+2*ww.size()];
            //
            std::complex<double> e_mi_k_x_d(cos(kappa*DirDotX), -sin(kappa*DirDotX));
            //
            std::complex<double> cTmp = e_mi_k_x_d * ww[gp] * d_R * d_R * cross;
#pragma omp critical
            {
              ffp[id + ik * numDir] += (dpdn[ik]
                                        + pp * std::complex<double>(0.0, kappa * DirDotNu)) * cTmp;
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)
          
        } // for (int id = 0; id < numDir; ++id)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(ffpLen, &ffp[0]);

  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  std::cout << " Surface " << surfaceCheck << " " << std::abs(surfaceCheck - d_R * d_R * 4*M_PI)/(4*d_R * d_R * M_PI);
  std::cout << std::endl;
#endif
  

  for (int ii = 0; ii < ffpLen; ++ii)
    ffp[ii] = ffp[ii] / (4.0 * M_PI);
  
  if (MyCom_p->cpuNum() == 0)
  {

    char *myName = new char[strlen(d_iod.output.transient.probes.prefix) + strlen(d_iod.output.transient.probes.farfieldpattern) + 1];
    sprintf(myName, "%s%s", d_iod.output.transient.probes.prefix, d_iod.output.transient.probes.farfieldpattern);

    // Output FFP in ASCII file
    FILE *out = fopen(myName, "w");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      fprintf(out, "%d %e ", ik, freqSamples[ik]);
      int iDir = 0;
      for (int ib = 0; ib <= numInc/2; ++ib)
      {
        double beta = ib * M_PI * 2.0 / numInc - 0.5 * M_PI;
        for (int ia = 0; ia < numInc; ++ia)
        {
          double alpha = ia * 2.0 * M_PI / numInc;
          fprintf(out, "%e %e ", alpha, beta);
          fprintf(out, "%e ", std::real(ffp[iDir + ik*numDir]));
          fprintf(out, "%e ", std::imag(ffp[iDir + ik*numDir]));
          fprintf(out, "%e ", std::abs(ffp[iDir + ik*numDir]));
          iDir += 1;
        }
      }
      fprintf(out, "\n");
    }
    fclose(out);
    delete[] myName;
  }
  
}


std::complex<double> KirchhoffIntegrator::evaluateSHS
(
 std::complex<double> *coeff,
 int nmax,
 double theta,
 double phi
 ) const
{
  
  std::complex<double> result(0.0, 0.0);
  
  int count = 0;
  std::vector< std::complex<double> > Yn(2*nmax + 1);
  for (int n = 0; n <= nmax; ++n)
  {
    //
    double maxEntry = 0.0;
    for (int in = 0; in <= 2*n; ++in)
      maxEntry = (maxEntry > std::abs(coeff[count+in])) ? maxEntry : std::abs(coeff[count+in]);
    if (maxEntry == 0.0)
    {
      count += 2*n+1;
      continue;
    }
    //
    sphericalHarmonic(n, theta, phi, Yn);
    //
    for (int m = -n; m <= n; ++m)
    {
      result += coeff[count] * Yn[m + n];
      count += 1;
    }
    //
  }
  
  return result;
  
}


void KirchhoffIntegrator::getQuadrature
(
  std::vector<double> &_xigauss,
  std::vector<double> &w
) const
{

  int rule = 12;

  if (w.size() != rule)
    w.resize(rule);

  if (_xigauss.size() != 3*rule)
    _xigauss.resize(3*rule);

  if (rule == 3)
  {

    w[0] = 1.0/3.0;
    w[1] = w[0];
    w[2] = w[0];

    _xigauss[0] = 4.0/6.0;
    _xigauss[1] = 1.0/6.0;
    _xigauss[2] = 1.0/6.0;

    _xigauss[3] = 1.0/6.0;
    _xigauss[4] = 4.0/6.0;
    _xigauss[5] = 1.0/6.0;

  }
  else if (rule == 12)
  {

    w[0] = 0.050844906370207;
    w[1] = 0.050844906370207;
    w[2] = 0.050844906370207;
    w[3] = 0.116786275726379;
    w[4] = 0.116786275726379;
    w[5] = 0.116786275726379;
    w[6] = 0.082851075618374; 
    w[7] = 0.082851075618374;
    w[8] = 0.082851075618374;
    w[9] = 0.082851075618374; 
    w[10] = 0.082851075618374;
    w[11] = 0.082851075618374;

    _xigauss[ 0] = 0.873821971016996;
    _xigauss[ 1] = 0.063089014491502;
    _xigauss[ 2] = 0.063089014491502;
    _xigauss[ 3] = 0.501426509658179;
    _xigauss[ 4] = 0.249286745170910;
    _xigauss[ 5] = 0.249286745170910;
    _xigauss[ 6] = 0.636502499121399;
    _xigauss[ 7] = 0.636502499121399;
    _xigauss[ 8] = 0.310352451033785;
    _xigauss[ 9] = 0.310352451033785;
    _xigauss[10] = 0.053145049844816;
    _xigauss[11] = 0.053145049844816;

    _xigauss[12] = 0.063089014491502;
    _xigauss[13] = 0.873821971016996;
    _xigauss[14] = 0.063089014491502;
    _xigauss[15] = 0.249286745170910;
    _xigauss[16] = 0.501426509658179;
    _xigauss[17] = 0.249286745170910;
    _xigauss[18] = 0.310352451033785;
    _xigauss[19] = 0.053145049844816;
    _xigauss[20] = 0.636502499121399;
    _xigauss[21] = 0.053145049844816;
    _xigauss[22] = 0.636502499121399;
    _xigauss[23] = 0.310352451033785;

  }

  for (int ig = 2*w.size(); ig < 3*w.size(); ++ig)
    _xigauss[ig] = 1.0 - _xigauss[ig-2*w.size()] - _xigauss[ig-w.size()];

}


std::complex<double> KirchhoffIntegrator::besselh
(
  const int n, 
  const double x
) const
{

  double hre, him;
#ifdef AEROACOUSTIC
  if (n >= 0)
  {
    hre = gsl_sf_bessel_jl(n, x);
    him = gsl_sf_bessel_yl(n, x);
  }
  else 
  {
    double scale = ((-n) % 2 == 0) ? 1.0 : -1.0;
    hre = scale * gsl_sf_bessel_jl(-n, x);
    him = scale * gsl_sf_bessel_yl(-n, x);
  }

  return std::complex<double>(hre, him);
#endif // AEROACOUSTIC
}


std::complex<double> KirchhoffIntegrator::besselh_prime
(
  const int n, 
  const double x
) const
{

  std::complex<double> h_nm1 = besselh(n-1, x);
  std::complex<double> h_np1 = besselh(n+1, x);

  std::complex<double> hprime;
  hprime = ((double) n)*h_nm1 - (n+1.0)*h_np1;
  hprime *= 1.0/(2.0*n+1.0);

  return hprime;

}


void KirchhoffIntegrator::sphericalHarmonic
(
  const int n,
  const double theta, 
  const double phi,
  std::vector< std::complex<double> > &Yn 
) const
{

//
// UH (08/30/2012)
// This routine has been verified against MATLAB.
//
#ifdef AEROACOUSTIC
  if (Yn.size() < 2*n+1)
    Yn.resize(2*n+1);

  double cost = cos(theta);

  std::vector<double> Pn(n+1, 0.0);
  for (int im = 0; im <= n; ++im)
    Pn[im] = gsl_sf_legendre_sphPlm(n, im, cost);

  for (int jj = 0; jj <= n; ++jj)
    Yn[n + jj] = Pn[jj] * std::complex<double>(cos(jj*phi), sin(jj*phi));

  double coeff = -1.0;
  for (int jj = 1; jj <= n; ++jj)
  {
    Yn[n - jj] = coeff * conj(Yn[n + jj]);
    coeff *= -1.0;
  }

#endif // AEROACOUSTIC
}


void KirchhoffIntegrator::cart2sph_x
(
  double xo, double yo, double zo,
  double &ro, double &to, double &po
) const
{

  ro = sqrt(xo*xo + yo*yo + zo*zo);
  to = acos(zo/ro);
  if ((abs(xo) < 1.0e-15 * ro) && (abs(yo) < 1.0e-15 * ro))
    po = 0.0;
  else if (abs(xo/ro) < 1.0e-15)
    po = (yo < 0.0) ? M_PI * 1.5 : M_PI * 0.5;
  else
    po = (yo < 0.0) ? 2*M_PI-acos(xo/(ro*sin(to))) : acos(xo/(ro*sin(to)));

}


void KirchhoffIntegrator::cart2sph_y
(
  double xo, double yo, double zo,
  double &ro, double &to, double &po
) const
{

  ro = sqrt(xo*xo + yo*yo + zo*zo);
  to = acos(zo/ro);
  if ((abs(xo) < 1.0e-15 * ro) && (abs(yo) < 1.0e-15 * ro))
    po = 0.0;
  else
    po = (xo > 0.0) ? asin(yo/(ro*sin(to))) : M_PI - asin(yo/(ro*sin(to)));

}


