#include <cstring>
#include <complex>
using std::complex;
#include <TsParameters.h>
#include <IoData.h>

#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif
#include <cstring>

#ifndef PI
#define PI 3.14159
#endif

//------------------------------------------------------------------------------

TsParameters::TsParameters(IoData &ioData)
{

  resType = ioData.ts.residual;

  cfllaw = ioData.ts.cfl.strategy;

  cfl0 = ioData.ts.cfl.cfl0;
  cflCoef1 = ioData.ts.cfl.cflCoef1;
  cflCoef2 = 0.0;
  cflMax = ioData.ts.cfl.cflMax;
  cflMin = ioData.ts.cfl.cflMin;
  dualtimecfl = ioData.ts.cfl.dualtimecfl;

  checksol = ioData.ts.cfl.checksol;
  checklinsolve = ioData.ts.cfl.checklinsolve;
  ser = ioData.ts.cfl.ser;
  angle_growth = ioData.ts.cfl.angle_growth;
  angle_zero = ioData.ts.cfl.angle_zero;
  dft_history = ioData.ts.cfl.dft_history;
  dft_freqcutoff = ioData.ts.cfl.dft_freqcutoff;
  dft_growth = ioData.ts.cfl.dft_growth;

  maxIts = ioData.ts.maxIts;
  eps = ioData.ts.eps;
  maxTime = ioData.ts.maxTime;

  output = new char[strlen(ioData.ts.output) + 1];
  sprintf(output, "%s", ioData.ts.output);

  cfl = cfl0;
  residual = 1.0;

  reshistory = new double[dft_history];
  for (int i=0; i<dft_history; i++) reshistory[i] = 1.0;
  dft = new complex<double>[dft_history];
  for (int i=0; i<dft_history; i++) dft[i]=0.0;

  unphysical = false;
  badlinsolve = false;

}

//------------------------------------------------------------------------------

TsParameters::~TsParameters()
{

  if (output) delete [] output;

}

//------------------------------------------------------------------------------

void TsParameters::computeCflNumber(int its, double res, double angle)
{
  //cfl = min( max( max(cflCoef1, cflCoef2*its), cfl0/pow(res,ser) ), cflMax );

  // First run automatic CFL checks
  if (unphysical){
    unphysical=false;
    badlinsolve=false;
    cfl *= 0.5;
    std::printf("Unphysicality detected. Reducing CFL number to %f.\n",cfl);
    if (cfl < cfl0/1000.) {std::printf("Could not resolve unphysicality by reducing CFL number. Aborting.\n"); exit(-1);}
    return;
  }
  if (badlinsolve){
    unphysical=false;
    badlinsolve=false;
    cfl *= 0.5;
    std::printf("Saturated linear solver detected. Reducing CFL number to %f.\n",cfl);
    if (cfl < cfl0/1000.) {std::printf("Linear solver does not converge for any feasible CFL number. Aborting.\n"); exit(-1);}
    return;
  }

  // Now run main code
  double cfl_res, cfl_dir, cfl_dft, hf_ratio;

  double cfl_prev = cfl;

  for (int i=dft_history-1; i>0; i--) reshistory[i]=reshistory[i-1];
  reshistory[0]=res;

  if (cfllaw == CFLData::RESIDUAL || cfllaw == CFLData::HYBRID){
    // compute residual strategy proposal
    cfl = cfl_prev;
    cfl *= pow(reshistory[1]/reshistory[0],ser);
    cfl_res = cfl;
    std::printf("CFL residual strategy: old residual: %e, new residual: %e, CFL proposal: %e\n",reshistory[1],reshistory[0],cfl_res);
  }
  if (cfllaw == CFLData::DIRECTION || cfllaw == CFLData::HYBRID){
    // compute direction strategy proposal
    cfl = cfl_prev;
    if (angle != -2.0) cfl *= pow(angle_growth,angle-angle_zero);
    cfl_dir = cfl;
    std::printf("CFL direction strategy: angle: %e, CFL proposal: %e\n",angle,cfl_dir);
  }
  if (cfllaw == CFLData::DFT || cfllaw == CFLData::HYBRID){
    // compute dft strategy proposal
    cfl = cfl_prev;

    double e_hf, e_ac, e_dc, e_total, hf_ratio;
    int cutofflow = dft_history/2 - (dft_freqcutoff-1)/2; 
    int cutoffhigh = (dft_history+1)/2 + (dft_freqcutoff-1)/2;

    for (int i=cutofflow; i<cutoffhigh+1; i++){
      dft[i]=0.0;
      for (int j=0; j<dft_history; j++) 
        dft[i] += 1/sqrt((double)dft_history)*complex<double>(reshistory[j],0.0)*exp(-2.0*complex<double>(0.0,1.0)*PI*(double)i*(double)j/(double)dft_history);
    }
    e_hf = 0.0;
    for (int i=cutofflow; i<cutoffhigh+1; i++)
      e_hf += norm(dft[i]);

    e_total = 0.0;
    for (int i=0; i<dft_history; i++) e_total += pow(reshistory[i],2);
    e_dc = 0.0;
    for (int i=0; i<dft_history; i++) e_dc += reshistory[i]/sqrt((double)dft_history);
    e_dc *= e_dc;
    e_ac = e_total-e_dc;

    hf_ratio = e_hf/e_ac;

    if (hf_ratio<0 || hf_ratio>1) {std::printf("Found invalid hf_ratio: %e, check for bugs in CFL Law\n",hf_ratio); exit(-1);}

    cfl *= pow(dft_growth, 1-2*hf_ratio);
    cfl_dft = cfl;
    std::printf("CFL DFT strategy: e_ac: %e, e_hf: %e, CFL proposal: %e\n",e_ac,e_hf,cfl_dft);
  }
  if (cfllaw == CFLData::HYBRID){
    // compute hybrid strategy cfl
    if (hf_ratio > 0.66) cfl = cfl_dft;
    else if (angle < angle_zero) cfl = cfl_dir;
    else cfl = max(cfl_dir, cfl_res);
    std::printf("CFL Hybrid strategy: dft_proposal: %e, direction proposal: %e, residual proposal: %e, chosen cfl: %e\n",cfl_dft,cfl_dir,cfl_res,cfl);
  }

  cfl = (min(max(cfl, cflCoef1),cflMax));
  std::printf("CFL number chosen: %e\n",cfl);
  return;
}

//------------------------------------------------------------------------------

// Included (MB)
void TsParameters::rstVar(IoData &ioData) {

  maxTime = ioData.ts.maxTime;

}

//------------------------------------------------------------------------------

