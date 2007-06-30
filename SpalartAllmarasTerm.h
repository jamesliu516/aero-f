#ifndef _SA_TERM_H_
#define _SA_TERM_H_

#include <IoData.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
using std::min;
#endif

//------------------------------------------------------------------------------
/*
@UNPUBLISHED{spalart-allmaras-92,
  author = "Spalart, P. R. and Allmaras, S. R.",
  title = "A One-Equation Turbulence Model for Aerodynamic Flows",
  month = jan,
  year = 1992,
  note = "{AIAA} paper 92-0439",
} 

Multigrid for the 2-D compressible Navier-Stokes equations
AIAA Paper 99-3336 
Steven R. Allmaras (Boeing Co., Seattle, WA)   
AIAA, Aerospace Sciences Meeting and Exhibit, 37th, Reno, NV, Jan. 11-14, 1999
*/

class SATerm {

  double oorey;
  
protected:

  double alpha;

  double cb1;
  double cb2;
  double cw1;
  double cw2;
  double cw3_pow6;
  double opcw3_pow;
  double cv1_pow3;
  double oocv2;
  double oosigma;
  double oovkcst2;

public:

  SATerm(IoData &);
  ~SATerm() {}

  double computeTurbulentViscosity(double *[4], double, double &);
  double computeTurbulentViscosity(double *, double);

  template<int neq, int shift>
  void computeJacobianVolumeTermSA(double [4][3], double [4], double [3][3], double, 
				   double, double *[4], double (*)[3][neq][neq],
				   double (*)[neq][neq]);

};

//------------------------------------------------------------------------------

inline
SATerm::SATerm(IoData &iod)
{

  oorey = 1.0 / iod.ref.reynolds_mu;
  alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.tc.prandtlTurbulent;

  cb1 = iod.eqs.tc.tm.sa.cb1;
  cb2 = iod.eqs.tc.tm.sa.cb2;
  cw2 = iod.eqs.tc.tm.sa.cw2;
  double cw3 = iod.eqs.tc.tm.sa.cw3;
  cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  double cv1 = iod.eqs.tc.tm.sa.cv1;
  cv1_pow3 = cv1*cv1*cv1;
  oocv2 = 1.0 / iod.eqs.tc.tm.sa.cv2;
  oosigma = 1.0 / iod.eqs.tc.tm.sa.sigma;
  oovkcst2 = 1.0 / (iod.eqs.tc.tm.sa.vkcst*iod.eqs.tc.tm.sa.vkcst);
  cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;

  cw1 /= iod.ref.reynolds_mu;
  oosigma /= iod.ref.reynolds_mu;

}

//------------------------------------------------------------------------------

inline
double SATerm::computeTurbulentViscosity(double *V[4], double mul, double &mutilde)
{

  mutilde = 0.25 * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		    V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  return mutilde*fv1;

}

//------------------------------------------------------------------------------

inline
double SATerm::computeTurbulentViscosity(double *V, double mul)
{

  double mutilde = V[0]*V[5];
  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  return mutilde*fv1;

}

//------------------------------------------------------------------------------

template<int neq, int shift>
void SATerm::computeJacobianVolumeTermSA(double dp1dxj[4][3], double d2w[4], 
					 double dudxj[3][3], double mul, double mutilde, 
					 double *V[4], double (*dRdU)[3][neq][neq], 
					 double (*dSdU)[neq][neq])
{

  const double sixth = 1.0/6.0;

  double mu5 = oosigma * (mul + mutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];
  double drdx = oosigma * 0.25 * dnutildedx;
  double drdy = oosigma * 0.25 * dnutildedy;
  double drdz = oosigma * 0.25 * dnutildedz;

  int k;
  for (k=0; k<4; ++k) {
    double nu = mu5 / V[k][0];
    dRdU[k][0][shift][shift] = drdx + nu * dp1dxj[k][0];
    dRdU[k][1][shift][shift] = drdy + nu * dp1dxj[k][1];
    dRdU[k][2][shift][shift] = drdz + nu * dp1dxj[k][2];
  }

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  if (d2wall < 1.e-15) {
    for (k=0; k<4; ++k)
      dSdU[k][shift][shift] = 0.0;
    return;
  }
  double chi = max(mutilde/mul, 0.001);
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);
  double fv2 = 1.0 + oocv2*chi;
  fv2 = 1.0 / (fv2*fv2*fv2);
  double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
  double ood2wall2 = 1.0 / (d2wall * d2wall);
  double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
  double oorho = 1.0 / rho;
  double zz = oorey * oovkcst2 * mutilde * oorho * ood2wall2;
  double s12 = dudxj[0][1] - dudxj[1][0];
  double s23 = dudxj[1][2] - dudxj[2][1];
  double s31 = dudxj[2][0] - dudxj[0][2];
  double s = sqrt(s12*s12 + s23*s23 + s31*s31);
  double Stilde = s*fv3 + zz*fv2;
  double rr = min(zz/Stilde, 2.0);
  double rr2 = rr*rr;
  double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
  double gg2 = gg*gg;
  double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);

  double chi2 = chi*chi;
  double dchi = 1.0 / mul;
  double coef1 = 1.0 / (chi3 + cv1_pow3);
  double dfv1 = 3.0*chi2*dchi*cv1_pow3 * coef1*coef1;
  double coef2 = 1.0 / (1.0 + chi*oocv2);
  double coef3 = coef2 * coef2;
  double dfv2 = -3.0*dchi*oocv2 * coef3*coef3;
  double dfv3 = ((dchi*fv1 + chi*dfv1)*(1.0 - fv2) - 
		 (1.0 + chi*fv1)*dfv2 - fv3*dchi) / chi;
  double dStilde = s*dfv3 + oorey*oovkcst2*oorho*ood2wall2 * (fv2 + mutilde*dfv2);
  double drr = oorey*oovkcst2*oorho*ood2wall2 * (Stilde - mutilde*dStilde) / (Stilde*Stilde);
  double dgg = (1.0 + cw2 * (6.0*rr2*rr2*rr - 1.0)) * drr;
  double dfw = pow(gg2*gg2*gg2 + cw3_pow6, 7.0*sixth);
  dfw = cw3_pow6 * opcw3_pow * dgg / dfw;

  double P = cb1 * Stilde;
  double D = cw1 * fw * oorho * mutilde * ood2wall2;
  double dP = cb1 * dStilde;
  double dD = cw1 * oorho * ood2wall2 * (fw + mutilde * dfw);
  double s00 = 0.25 * (max(D - P, 0.0) + max(dD - dP, 0.0) * mutilde);
  //s00 = 0.25 * (D - P + (dD - dP) * mutilde);
  double coef4 = oosigma * cb2 * rho * 2.0;

  for (k=0; k<4; ++k)
    dSdU[k][shift][shift] = coef4 / V[k][0] * 
      (dnutildedx*dp1dxj[k][0] + dnutildedy*dp1dxj[k][1] + dnutildedz*dp1dxj[k][2]) - s00;
  
}

//------------------------------------------------------------------------------

#endif
