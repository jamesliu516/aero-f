#include <FluxFcnDesc.h>
#include <LinkF77.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif


extern "C" {

// standard functions needed for computation of fluxes and jacobians


// PERFECT GAS
  void F77NAME(roeflux1)(const double&, const double&, const double&, double*, 
			 const double&, double*, double*, double*,
			 const double&, const double&, const double&,
                         const double&, const int&);
  void F77NAME(roeflux5)(const int&, const double&, const double&, const double&, double*,
			                   const double&, double*, double*, double*, double*, double*, 
                         const double&, const double&, const double&, const double&, const int&);
  void F77NAME(roeflux5prim)(const int&, const double&, const double&, const double&, double*,
			                   const double&, double*, double*, double*, double*, double*, 
                         const double&, const double&, const double&, const double&, 
                         const double&, const int&);
  void F77NAME(roeflux6)(const int&, const double&, const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*);
  void F77NAME(roejac2)(const int&, const double&, const double&, 
			const double&, double*,
			const double&, double*, double*, double*, double*,
			const double&, const double&, const double&,
			const double&, const int&);
  void F77NAME(roejac5)(const int&, const double&, const double&, const double&, double*, 
			const double&, double*, double*, double*);
  void F77NAME(roejac6)(const int&, const double&, const double&, const double&, double*, 
			const double&, double*, double*, double*,const int&,
                        const double&, const double&, const double&, const double&, const int&);
  void F77NAME(boundflux5)(const int&, const double&, double*, const double&, 
			   double*, double*, double*);
  void F77NAME(boundjac2)(const int&, const double&, double*, const double&, 
			  double*, double*, double*);
  void F77NAME(genbcfluxgas)(const int&, const double &, const double&, double*, 
			     const double&, double*, double*, double*);

// BAROTROPIC LIQUID			  
  void F77NAME(roeflux1water)(const double&, const double&, const double&,
			 const double&, const double&, double*,
                         const double&, double*, double*, double*);
//  void F77NAME(roejac2water)();   //for jacobian with turbulence in segregated solvers
  void F77NAME(roeflux5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*,
                         const double&, const double&, const double&, const double&, const int&);
  void F77NAME(roejac5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, 
			 const double&, const double&, const double&, 
			 const double&, const int&);
  void F77NAME(genbcfluxtait)(const int&, const double&, const double&, const double&, 
			 const double&,  double*, const double&, double*, double*, double*);
};

//------------------------------------------------------------------------------
/*
@ARTICLE{roe-81,
  author = "Roe, P. L.",
  title = "Approximate {R}iemann Solvers, Parameter Vectors
           and Difference Schemes",
  journal = j. comp. phys.,
  year = 1981,
  volume = 43,
  pages = "357--372",
} 
*/
void FluxFcnFDJacRoeEuler3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{

  F77NAME(roeflux5prim)(0, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux, betaRef, k1, cmach, irey, length, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeEuler3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux)
{
  F77NAME(roeflux5prim)(0, gamma, vfgam, vfp, normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, betaRef, k1, cmach, irey, length, prec);
}

//------------------------------------------------------------------------------

template<int dim>
inline
void roejacexact3D(int type, double gamma, VarFcn* varFcn, FluxFcn::Type typeJac, double* normal, 
		   double normalVel, double* VL, double* VR, double* jacL, double* jacR)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdVL[dim2], dfdVR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};
  
  if (varFcn->getType() == VarFcn::GAS){ // || varFcn->getType() == VarFcn::GASINGAS){
    F77NAME(roejac5)(type, gamma, varFcn->getGamma(), varFcn->getPressureConstant(), n, normalVel, VL, VR, dfdVL);
    n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
    F77NAME(roejac5)(type, gamma, varFcn->getGamma(), varFcn->getPressureConstant(), n, -normalVel, VR, VL, dfdVR);

  } else if(varFcn->getType() == VarFcn::LIQUID || varFcn->getType() == VarFcn::LIQUIDINLIQUID){
/*
     F77NAME(roejacwaterdissprim)(type, gamma, varFcn->getCv(), varFcn->getPrefWater(), varFcn->getAlphaWater(),
                                                      varFcn->getBetaWater(), n, normalVel, VL, VR, dfdVL);
     n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
     F77NAME(roejacwaterdissprim)(type, gamma, varFcn->getCv(), varFcn->getPrefWater(), varFcn->getAlphaWater(),
                                                      varFcn->getBetaWater(), n, -normalVel, VR, VL, dfdVR);
*/
  }
  
  
  int k;
  for (k=0; k<dim2; ++k)
    dfdVR[k] = -dfdVR[k];

  if (type == 1 || type == 2) {
    double f1;
    if(varFcn->getType() == VarFcn::GAS) // || varFcn->getType() == VarFcn::GASINGAS)
      F77NAME(roeflux1)(gamma, varFcn->getGamma(), varFcn->getPressureConstant(), normal, normalVel, VL, VR, &f1,1.0,1.0,0.0,0.0,0);
    else if (varFcn->getType() == VarFcn::LIQUID || varFcn->getType() == VarFcn::LIQUIDINLIQUID){
      //F77NAME(roeflux1water)(gamma, varFcn->getCv(), varFcn->getPrefWater(), varFcn->getAlphaWater(),
	//	      varFcn->getBetaWater(), normal, normalVel, VL, VR, &f1);
    }


    if (type == 1) {
      if (f1 >= 0.0) {
	dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
	dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
	dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
	dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
	dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
	dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
	dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
	dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
	dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
	dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
	dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
	dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
	dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
	dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
	dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
	dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
	dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
	dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
	dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

	dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
	dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
	dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
	dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
	dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
	dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
	dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
	dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
	dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
	dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
	dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
	dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
	dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
	dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

	dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
	dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
	dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
	dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
	dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
	dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
	dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
  }

  if (typeJac == FluxFcn::CONSERVATIVE) {
    varFcn->postMultiplyBydVdU(VL, dfdVL, jacL);
    varFcn->postMultiplyBydVdU(VR, dfdVR, jacR);
  }
  else {
    for (k=0; k<dim2; ++k) { 
      jacL[k] = dfdVL[k]; 
      jacR[k] = dfdVR[k]; 
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void roejacappr3Dgas(int type, double gamma, VarFcn* varFcn, double vfgam, double vfp, 
                  FluxFcn::Type typeJac, double* normal, double normalVel, double* VL, 
		  double* VR, double* jacL, double* jacR, double irey, double betaRef, 
		  double k1, double cmach, int prec, int flag)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;
  double dflag = double(flag);


  double dfdVL[dim2], dfdVR[dim2];
  double dfdUL[dim2], dfdUR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};


  F77NAME(roejac6)(type, gamma, vfgam, vfp, n, normalVel, VL, VR, dfdUL,1, betaRef, k1, cmach, irey, prec);
  F77NAME(roejac6)(type, gamma, vfgam, vfp, n, normalVel, VR, VL, dfdUR,2, betaRef, k1, cmach, irey, prec);

  if (type == 1 || type == 2) {
    varFcn->postMultiplyBydUdV(VL, dfdUL, dfdVL, dflag);
    varFcn->postMultiplyBydUdV(VR, dfdUR, dfdVR, dflag);

    double f1;
    F77NAME(roeflux1)(gamma, vfgam, vfp, n, normalVel, VL, VR, &f1, betaRef,k1,cmach,irey,prec);

    if (type == 1) {
      if (f1 >= 0.0) {
        dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
        dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
        dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
        dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    varFcn->postMultiplyBydVdU(VL, dfdVL, dfdUL, dflag);
    varFcn->postMultiplyBydVdU(VR, dfdVR, dfdUR, dflag);
  }

  int k;

  if (typeJac == FluxFcn::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) {
      jacL[k] = dfdUL[k];
      jacR[k] = dfdUR[k];
    }
  }
  else {
    varFcn->postMultiplyBydUdV(VL, dfdUL, jacL, dflag);
    varFcn->postMultiplyBydUdV(VR, dfdUR, jacR, dflag);
  }

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeEuler3D::computeJacobiansPerfectGas(double irey, double vfgam, double vfp, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, int flag)
{
	
  roejacappr3Dgas<5>(0, gamma, vf, vfgam, vfp, type, normal, normalVel, VL, VR, jacL, jacR, irey, betaRef, k1, cmach, prec, flag);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					double *VL, double *VR, double *flux)
{

  F77NAME(roeflux6)(0, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeEuler3D::computeJacobiansPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR)
{
	
  roejacexact3D<5>(0, gamma, vf, type, normal, normalVel, VL, VR, jacL, jacR);

}

//-----------------------------------------------------------------------
// VL and VR are the primitive state variables !!

void FluxFcnVanLeerEuler3D::evalFlux(double vfgam, double vfp, double *normal, double normalVel, 
				     double *V, double *f, int sign) 
{

  //gamma 
  double gam = vfgam;
  double invgam = 1.0/gam;
  double gam1 = gam - 1.0;
  double invgam1 = 1.0/gam1;

  // compute norm to remove area factor from computations
  double norm = sqrt( normal[0]*normal[0]
                    + normal[1]*normal[1]
                    + normal[2]*normal[2] );

  double invNorm = 1.0 / norm;
  	            	     
  // normalize the normal 
  double nVec[3];
  nVec[0] = normal[0]*invNorm;
  nVec[1] = normal[1]*invNorm;
  nVec[2] = normal[2]*invNorm;
  
  //normalize grid velocity
  double gridVel = normalVel * invNorm;

  //compute speed of sound
  double a = vf->computeSoundSpeed(V);
  
  // compute fluid velocity
  double fluidVel = V[1]*nVec[0] + V[2]*nVec[1] + V[3]*nVec[2];

  // compute total velocity --> fluid vel. - grid vel.
  double totVel = fluidVel - gridVel;

  //compute normal mach numbers
  double mach = totVel / a;

  // get conservative variables
  double U[5];
  vf->primitiveToConservative(V, U);
    
  double factor = mach*sign;
  if (factor >= 1.0)  {
    f[0] = norm * U[0] * totVel;

    f[1] = norm * ( U[1] * totVel
                  + V[4] * nVec[0] ); 

    f[2] = norm * ( U[2] * totVel
                  + V[4] * nVec[1] ); 

    f[3] = norm * ( U[3] * totVel
         	  + V[4] * nVec[2] ); 

    f[4] = norm * ( U[4] * totVel
 		  + V[4] * fluidVel );
  }  
  else if (factor > -1.0)  {
    // compute frequently used terms
    double h = 0.25 * norm * sign * V[0] * a * (mach+sign)*(mach+sign);  
    double f1 = (2*a*sign - totVel) * invgam;

    f[0] = h;

    f[1] = h * (V[1] + f1 * nVec[0]);

    f[2] = h * (V[2] + f1 * nVec[1]);

    f[3] = h * (V[3] + f1 * nVec[2]);


    f[4] = h * ( (2*a*a + 2*sign * gam1 * totVel * a 
                 - totVel*totVel*gam1)/(gam*gam - 1.0)
               + .5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
               + gridVel*f1 ) ;
  }  
  else {
    for (int i = 0; i < 5; i++)
      f[i] = 0;  
  }

}

//------------------------------------------------------------------------------

void FluxFcnVanLeerEuler3D::evalJac(double vfgam, double vfp, double *normal, double normalVel,
				    double *V, double *jac, int sign)

{

  //gamma 
  double gam = vfgam;
  double invgam = 1.0/gam;
  double gam1 = gam - 1.0;
  double invgam1 = 1.0/gam1;

  // compute norm to remove area factor from computations
  double norm = sqrt( normal[0]*normal[0]
                    + normal[1]*normal[1]
                    + normal[2]*normal[2] );

  double invNorm = 1.0 / norm;
 
  // normalize the normal
  double nVec[3];
  nVec[0] = normal[0]*invNorm;
  nVec[1] = normal[1]*invNorm;
  nVec[2] = normal[2]*invNorm;
 
  //normalize grid velocity
  double gridVel = normalVel * invNorm;

  //compute speed of sound
  double a = vf->computeSoundSpeed(V);

  // compute fluid velocity
  double fluidVel = V[1]*nVec[0] + V[2]*nVec[1] + V[3]*nVec[2];

  // compute total velocity --> fluid vel. - grid vel.
  double totVel = fluidVel - gridVel;

  //compute normal mach numbers
  double mach = totVel / a; 

  // get conservative variables
  double U[5];
  vf->primitiveToConservative(V, U);

  // compute frequently used terms
  double q2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
  // compute Energy term
  double E = U[4] / U[0];

  // jacs are computed in the following order: df0/du0, df1/du0, 
  // df2/du0 ... df0/du1, df2/du1, ... df4/du4
  double factor = mach*sign;
  if (factor >= 1.0)  {
    double pOverRho = V[4] / V[0];
    jac[0] = -gridVel;
    jac[1] = norm * nVec[0];
    jac[2] = norm * nVec[1];
    jac[3] = norm * nVec[2];
    jac[4] = 0;

    jac[5] = norm * (0.5*gam1*q2*nVec[0] - V[1]*fluidVel);
    jac[6] = norm * ( totVel + V[1]*nVec[0]*(2.0-gam) );
    jac[7] = norm * ( V[1] * nVec[1] - gam1 * V[2] * nVec[0] );
    jac[8] = norm * ( V[1] * nVec[2] - gam1 * V[3] * nVec[0] );
    jac[9] = norm * gam1 * nVec[0];

    jac[10] = norm * (0.5*gam1*q2*nVec[1] - V[2]*fluidVel);
    jac[11] = norm * ( V[2] * nVec[0] - gam1 * V[1] * nVec[1] );
    jac[12] = norm * ( totVel + V[2]*nVec[1]*(2.0-gam) );
    jac[13] = norm * ( V[2] * nVec[2] - gam1 * V[3] * nVec[1] );
    jac[14] = norm * gam1 * nVec[1];

    jac[15] = norm * (0.5*gam1*q2*nVec[2] - V[3]*fluidVel);
    jac[16] = norm * ( V[3] * nVec[0] - gam1 * V[1] * nVec[2] );
    jac[17] = norm * ( V[3] * nVec[1] - gam1 * V[2] * nVec[2] );
    jac[18] = norm * ( totVel + V[3]*nVec[2]*(2.0-gam) );
    jac[19] = norm * gam1 * nVec[2];

    jac[20] = norm*fluidVel * ( 0.5*+gam1*q2 - E - pOverRho );
    jac[21] = norm * ( nVec[0]*(E + pOverRho) - gam1*V[1]*fluidVel );
    jac[22] = norm * ( nVec[1]*(E + pOverRho) - gam1*V[2]*fluidVel );
    jac[23] = norm * ( nVec[2]*(E + pOverRho) - gam1*V[3]*fluidVel );
    jac[24] = norm * (fluidVel*gam - gridVel);    
  }
  else if (factor > -1.0)  {
    // compute frequently used terms
    double k3 = 1.0 - sign * mach;
    double k4 = 0.5*(1.0 + sign * mach);
    double coef = sign*V[0]*a*k4*k4;
    double f2 = invgam*a*sign;
    double invgamPlus1 = 1.0 / (gam + 1.0);
    double invRho = 1.0 / V[0];
    double invA = 1.0/a;
    double ggam1 = gam*gam1;
 
    double g[4];
    g[0] = V[1] + nVec[0]*f2*(1.0 + k3);
    g[1] = V[2] + nVec[1]*f2*(1.0 + k3);
    g[2] = V[3] + nVec[2]*f2*(1.0 + k3);
    g[3] = f2*gridVel*(1.0 + k3) + 0.5*q2 + a*a*invgam1 
	 - a*a*k3*k3*invgamPlus1; 

    double gGrad[4][5];
    //dg0/du
    gGrad[0][0] = invRho * 
		( nVec[0]*(gridVel*invgam + gam1*invA*sign*E) - g[0] );
    gGrad[0][1] = invRho * ( 1.0 - sign*invgam*nVec[0]
		           * (ggam1*V[1]*invA + sign*nVec[0]) );
    gGrad[0][2] = -sign*invgam*invRho*nVec[0]
        	* ( ggam1*invA*V[2] + sign*nVec[1] );
    gGrad[0][3] = -sign*invgam*invRho*nVec[0]
        	* ( ggam1*invA*V[3] + sign*nVec[2] );
    gGrad[0][4] = sign * invA * invRho * gam1 * nVec[0];

    // dg1/du
    gGrad[1][0] = invRho * ( nVec[1]*(gridVel*invgam + gam1*invA*sign*E)
                           - g[1] );
    gGrad[1][1] = -sign*invgam*invRho*nVec[1]
                * ( ggam1*invA*V[1] + sign*nVec[0] );
    gGrad[1][2] = invRho * ( 1.0 - sign*invgam*nVec[1]
                           * (ggam1*V[2]*invA + sign*nVec[1]) );
    gGrad[1][3] = -sign*invgam*invRho*nVec[1]
                * ( ggam1*invA*V[3] + sign*nVec[2] );
    gGrad[1][4] = sign * invA * invRho * gam1 * nVec[1];

    // dg2/du
    gGrad[2][0] = invRho * ( nVec[2]*(gridVel*invgam + gam1*invA*sign*E)
                           - g[2] );
    gGrad[2][1] = -sign*invgam*invRho*nVec[2]
                * ( ggam1*invA*V[1] + sign*nVec[0] );
    gGrad[2][2] = -sign*invgam*invRho*nVec[2]
                * ( ggam1*invA*V[2] + sign*nVec[1] );
    gGrad[2][3] = invRho * ( 1.0 - sign*invgam*nVec[2]
                           * (ggam1*V[3]*invA + sign*nVec[2]) );
    gGrad[2][4] = sign * invA * invRho * gam1 * nVec[2];

    //dg3/du
    gGrad[3][0] = invRho * ( (gam1*q2 - gam*E)
                           + k3*invgamPlus1 * ( 2*sign*mach*a*a 
                                              + ggam1*(q2-E) 
					      + 2*a*sign*gridVel )
    			   + gridVel * ( gridVel*invgam + sign*E*invA*gam1
                                       - f2*(1.0 + k3) ) );  

    gGrad[3][1] = invRho*( k3*invgamPlus1 * (ggam1*V[1]+2*a*sign*nVec[0])
                         - gridVel * sign * invgam
			 * (ggam1*invA*V[1] + sign*nVec[0]) - gam1*V[1] );
    gGrad[3][2] = invRho*( k3*invgamPlus1 * (ggam1*V[2]+2*a*sign*nVec[1])
                         - gridVel * sign * invgam
		  	 * (ggam1*invA*V[2] + sign*nVec[1]) - gam1*V[2] );
    gGrad[3][3] = invRho*( k3*invgamPlus1 * (ggam1*V[3]+2*a*sign*nVec[2])
                         - gridVel*sign*invgam*(ggam1*invA*V[3]+sign*nVec[2])
                         - gam1*V[3] );
    gGrad[3][4] = invRho*(gam - k3*ggam1*invgamPlus1 + gridVel*invA*sign*gam1); 

    //df/du0
    jac[0] = k4*(0.25*ggam1*sign*E*invA*k3 - gridVel);
    jac[1] = k4*(nVec[0] - 0.25*ggam1*sign*invA*k3*V[1]);
    jac[2] = k4*(nVec[1] - 0.25*ggam1*sign*invA*k3*V[2]);
    jac[3] = k4*(nVec[2] - 0.25*ggam1*sign*invA*k3*V[3]);
    jac[4] = 0.25*k4*k3*ggam1*sign*invA; 
    
    //df1/du through df4/du
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 5; j++)
        jac[5*(i+1) + j] = norm * (jac[j]*g[i] + coef*gGrad[i][j]);

    jac[0] *= norm;
    jac[1] *= norm;
    jac[2] *= norm;
    jac[3] *= norm;
    jac[4] *= norm;
     
  }
  else {
    for (int i = 0; i < 25; i++) 
      jac[i] = 0.0;
  }

}

//------------------------------------------------------------------------------
// VL and VR are the primitive state variables !!

void FluxFcnVanLeerEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel,
				    double *VL, double *VR, double *flux)
{

  //compute the split fluxes
  double fPlus[5], fMinus[5];
  evalFlux(vfgam, vfp, normal, normalVel, VL, fPlus, 1);
  evalFlux(vfgam, vfp, normal, normalVel, VR, fMinus, -1);

  for (int i = 0; i < 5; i++)
    flux[i] = fPlus[i] + fMinus[i];

} 

//------------------------------------------------------------------------------

void FluxFcnVanLeerEuler3D::computeJacobiansPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					     double *VL, double *VR, 
					     double *jacL, double *jacR)
{

  //compute the jacs for left and right states
  evalJac(vfgam, vfp, normal, normalVel, VL, jacL, 1);
  evalJac(vfgam, vfp, normal, normalVel, VR, jacR, -1);

} 

//------------------------------------------------------------------------------

void FluxFcnWallEuler3D::computePerfectGas(double *normal, double normalVel,
                                 double *V, double *Ub, double *flux)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;
 
}
//------------------------------------------------------------------------------

void FluxFcnGhidagliaEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux)
{

  F77NAME(genbcfluxgas)(0, vfgam, vfp, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnInflowEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux)
{
  F77NAME(boundflux5)(0, vfgam, normal, normalVel, V, Ub, flux);

}
//------------------------------------------------------------------------------

void FluxFcnOutflowEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux)
{

  F77NAME(boundflux5)(0, vfgam, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void influx3D(int type, VarFcn* varFcn, double* normal, 
	      double normalVel, double* V, double* Ub, double* flux, int flag)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double dflag = double(flag);

  double gam, gam1, invgam1, pstiff;
  if(flag==1){
    gam = varFcn->getGamma();
    pstiff = varFcn->getPressureConstant();
  }else{
    gam = varFcn->getGammabis();
    pstiff = varFcn->getPressureConstantbis();
  }
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double Vb[dim];
  varFcn->conservativeToPrimitive(Ub, Vb, dflag);
 
  double rho, u, v, w, p, nut, eps, k;
  
  if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] == 0.0){
    rho = Vb[0];
    u = V[1];
    v = V[2];
    w = V[3];
    p = V[4];
  }

  else if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] <= 0.0){
    rho = Vb[0];
    u = Vb[1];
    v = Vb[2];
    w = Vb[3];
    p = V[4];
    if (type == 1)
      nut = Vb[dimm1];
    else if (type == 2){
      k = Vb[dimm2];
      eps = Vb[dimm2];
    }
  }else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    p = Vb[4];
    if (type == 1)
      nut = V[dimm1];
    else if (type == 2){
      k = V[dimm2];
      eps = V[dimm2];
    }

  }
  double rhoun = rho * (u*n[0] + v*n[1] + w*n[2] - nVel);
  double q = u*u + v*v + w*w;

  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*u + p*n[0]);
  flux[2] = S * (rhoun*v + p*n[1]);
  flux[3] = S * (rhoun*w + p*n[2]);
  flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);

  if (type == 1) {
    flux[5] = S * rhoun*nut;
  }
  else if (type == 2) {
    flux[5] = S * rhoun*k;
    flux[6] = S * rhoun*eps;
  }

}

//------------------------------------------------------------------------------
template<int dim>
inline
void jacinflux3D(int type, VarFcn* varFcn, FluxFcn::Type typeJac, 
		 double* normal, double normalVel, double* V, double* Ub, double* jac, int flag)
{

  double _dfdV[dim][dim];
  double* dfdV = reinterpret_cast<double*>(_dfdV);
  int j;
  for (j=0; j<dim*dim; ++j)
    dfdV[j] = 0.0;

  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double dflag = double(flag);
                                                                                                    
  double gam, gam1, invgam1, pstiff;
  if(flag==1){
    gam = varFcn->getGamma();
    pstiff = varFcn->getPressureConstant();
  }else{
    gam = varFcn->getGammabis();
    pstiff = varFcn->getPressureConstantbis();
  }
  gam1 = gam - 1.0;
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double Vb[dim];
  varFcn->conservativeToPrimitive(Ub, Vb, dflag);
  double rho, u, v, w, p;

  if(V[1]*n[0]+V[2]*n[1]+V[3]*n[2]<= 0.0){
    u = Vb[1];
    v = Vb[2];
    w = Vb[3];

    _dfdV[1][4] = S * n[0];
    _dfdV[2][4] = S * n[1];
    _dfdV[3][4] = S * n[2];
    _dfdV[4][4] = S * (gam*invgam1 * (u*n[0] + v*n[1] + w*n[2] - nVel) + nVel);
  }else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    p = Vb[4];
                                                                                                                                                          
    double un = u*n[0] + v*n[1] + w*n[2] - nVel;
    double q = u*u + v*v + w*w;
    double rhoH = 0.5*rho*q + gam*invgam1*(p+pstiff);
                                                                                                                                                                                                     
    _dfdV[0][0] = S * un;
    _dfdV[0][1] = S * rho*n[0];
    _dfdV[0][2] = S * rho*n[1];
    _dfdV[0][3] = S * rho*n[2];
    _dfdV[1][0] = S * u*un;
    _dfdV[1][1] = S * (rho*un + rho*u*n[0]);
    _dfdV[1][2] = S * rho*u*n[1];
    _dfdV[1][3] = S * rho*u*n[2];
    _dfdV[2][0] = S * v*un;
    _dfdV[2][1] = S * rho*v*n[0];
    _dfdV[2][2] = S * (rho*un + rho*v*n[1]);
    _dfdV[2][3] = S * rho*v*n[2];
    _dfdV[3][0] = S * w*un;
    _dfdV[3][1] = S * rho*w*n[0];
    _dfdV[3][2] = S * rho*w*n[1];
    _dfdV[3][3] = S * (rho*un + rho*w*n[2]);
    _dfdV[4][0] = S * 0.5*q*un;
    _dfdV[4][1] = S * (rho*u*un + rhoH*n[0]);
    _dfdV[4][2] = S * (rho*v*un + rhoH*n[1]);
    _dfdV[4][3] = S * (rho*w*un + rhoH*n[2]);
                                                                                                                                                                                                     
    if (type == 1) {
      double nut = V[dimm1];
      _dfdV[dimm1][0] = S * nut*un;
      _dfdV[dimm1][1] = S * rho*nut*n[0];
      _dfdV[dimm1][2] = S * rho*nut*n[1];
      _dfdV[dimm1][3] = S * rho*nut*n[2];
      _dfdV[dimm1][dimm1] = S * rho*un;
    }
    else if (type == 2) {
      double k = V[dimm2];
      double eps = V[dimm1];
      _dfdV[dimm2][0] = S * k*un;
      _dfdV[dimm2][1] = S * rho*k*n[0];
      _dfdV[dimm2][2] = S * rho*k*n[1];
      _dfdV[dimm2][3] = S * rho*k*n[2];
      _dfdV[dimm2][dimm2] = S * rho*un;
      _dfdV[dimm1][0] = S * eps*un;
      _dfdV[dimm1][1] = S * rho*eps*n[0];
      _dfdV[dimm1][2] = S * rho*eps*n[1];
      _dfdV[dimm1][3] = S * rho*eps*n[2];
      _dfdV[dimm1][dimm1] = S * rho*un;
    }
  }


  if (typeJac == FluxFcn::CONSERVATIVE)
    varFcn->postMultiplyBydVdU(V, dfdV, jac, dflag);
  else
    for (j=0; j<dim*dim; ++j)
      jac[j] = dfdV[j]; 

}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, int flag)
{

  influx3D<5>(0, vf, normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowEuler3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, int flag)
{

  jacinflux3D<5>(0, vf, type, normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowEuler3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, int flag)
{

  influx3D<5>(0, vf, normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowEuler3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, int flag)
{

  jacinflux3D<5>(0, vf, type, normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
//turbulence

void FluxFcnFDJacRoeSA3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				  double *VL, double *VR, double *flux)
{

  F77NAME(roeflux5)(1, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux, betaRef, k1, cmach, irey, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeSA3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				    double *VL, double *VR, double *flux)
{

  F77NAME(roeflux5)(1, gamma, vfgam, vfp, normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, betaRef, k1, cmach, irey, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeSA3D::computeJacobiansPerfectGas(double irey, double vfgam, double vfp, double *normal, double normalVel, 
					     double *VL, double *VR, 
					     double *jacL, double *jacR, int flag)
{

  roejacappr3Dgas<6>(1, gamma, vf, vfgam, vfp, type, normal, normalVel, VL, VR, jacL, jacR, irey, betaRef, k1, cmach, prec, flag);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeSA3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{

  F77NAME(roeflux6)(1, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeSA3D::computeJacobiansPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					      double *VL, double *VR, 
					      double *jacL, double *jacR)
{

  roejacexact3D<6>(1, gamma, vf, type, normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnWallSA3D::computePerfectGas(double *normal, double normalVel, 
			      double *V, double *Ub, double *flux)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;
  flux[5] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnOutflowSA3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				 double *V, double *Ub, double *flux)
{

  F77NAME(boundflux5)(1, vfgam, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowSA3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					double *V, double *Ub, double *flux, int flag)
{

  influx3D<6>(1, vf, normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowSA3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						double *V, double *Ub, double *jacL, int flag)
{

  jacinflux3D<6>(1, vf, type, normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowSA3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					 double *V, double *Ub, double *flux, int flag)
{

  influx3D<6>(1, vf, normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowSA3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						 double *V, double *Ub, double *jacL, int flag)
{

  jacinflux3D<6>(1, vf, type, normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR

void FluxFcnRoeSAturb3D::computeJacobiansPerfectGas(double irey, double vfgam, double vfp, double *normal, double normalVel, 
					  double *VL, double *VR, 
					  double *jacL, double *jacR)
{

  F77NAME(roejac2)(1, gamma, vfgam, vfp, normal, normalVel, VL, VR, jacL, jacR, betaRef,k1,cmach,irey,prec);

}

//------------------------------------------------------------------------------

void FluxFcnWallSAturb3D::computeJacobianPerfectGas(double *normal, double normalVel, 
					  double *V, double *Ub, double *jacL)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnOutflowSAturb3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					     double *V, double *Ub, double *jacL)
{

  F77NAME(boundjac2)(1, vfgam, normal, normalVel, V, Ub, jacL);

}
//------------------------------------------------------------------------------

void FluxFcnInternalInflowSAturb3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, int flag)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowSAturb3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
						     double *V, double *Ub, double *jacL, int flag)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = V[0];
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (type == FluxFcn::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------

void FluxFcnFDJacRoeKE3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				  double *VL, double *VR, double *flux)
{

  F77NAME(roeflux5)(2, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux, betaRef, k1, cmach, irey, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeKE3D::computePerfectGas(double length, double irey, double vfgam, double vfp, double *normal, double normalVel, 
				    double *VL, double *VR, double *flux)
{

  F77NAME(roeflux5)(2, gamma, vfgam, vfp, normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, betaRef, k1, cmach, irey, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeKE3D::computeJacobiansPerfectGas(double irey, double vfgam, double vfp, double *normal, double normalVel, 
					     double *VL, double *VR, 
					     double *jacL, double *jacR, int flag)
{

  roejacappr3Dgas<7>(2, gamma, vf, vfgam, vfp, type, normal, normalVel, VL, VR, jacL, jacR, irey, betaRef, k1, cmach, prec, flag);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeKE3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{

  F77NAME(roeflux6)(2, gamma, vfgam, vfp, normal, normalVel, VL, VL, VR, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeKE3D::computeJacobiansPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					      double *VL, double *VR, 
					      double *jacL, double *jacR)
{

  roejacexact3D<7>(2, gamma, vf, type, normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnWallKE3D::computePerfectGas(double *normal, double normalVel, 
			      double *V, double *Ub, double *flux)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;
  flux[5] = 0.0;
  flux[6] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnOutflowKE3D::computePerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
				 double *V, double *Ub, double *flux)
{

  F77NAME(boundflux5)(2, vfgam, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR

void FluxFcnRoeKEturb3D::computeJacobiansPerfectGas(double irey, double vfgam, double vfp,
					  double *normal, double normalVel, 
					  double *VL, double *VR, 
					  double *jacL, double *jacR)
{

  F77NAME(roejac2)(1, gamma, vfgam, vfp, normal, normalVel, VL, VR, jacL, jacR, betaRef,k1,cmach,irey,prec);

}

//------------------------------------------------------------------------------

void FluxFcnWallKEturb3D::computeJacobianPerfectGas(double *normal, double normalVel, 
					  double *V, double *Ub, double *jacL)
{

  jacL[0] = 0.0;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = 0.0;

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnOutflowKEturb3D::computeJacobianPerfectGas(double vfgam, double vfp, double *normal, double normalVel, 
					     double *V, double *Ub, double *jacL)
{

  F77NAME(boundjac2)(2, vfgam, normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------BAROTROPIC LIQUID -----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//--------------------------------------------------------------------
/*
@ARTICLE{roe-81,
  author = "Roe, P. L.",
  title = "Approximate {R}iemann Solvers, Parameter Vectors
           and Difference Schemes",
  journal = j. comp. phys.,
  year = 1981,
  volume = 43,
  pages = "357--372",
} 
*/

//------------------------------------------------------------------------------

template<int dim>
inline
void roejacappr3Dwater(int type, double gamma, VarFcn* varFcn, double vfcv, double vfa, 
                  double vfb, double vfp, FluxFcn::Type typeJac, double* normal, 
		              double normalVel, double* VL, double* VR, double* jacL, 
		              double* jacR, double irey, double betaRef, double k1, 
									double cmach, int prec, int flag)
{

  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;
  double dflag = double(flag);


  double dfdUL[dim2], dfdUR[dim2];
  for (int kk = 0; kk<dim2; kk++) dfdUR[kk] = 0.0;
  double n[3] = {normal[0], normal[1], normal[2]};

  F77NAME(roejac5waterdissprec)(type, gamma, vfcv, vfp, vfa, vfb, n, normalVel, 
			                          VL, VR, dfdUL, betaRef, k1, cmach, irey, prec);
  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  F77NAME(roejac5waterdissprec)(type, gamma, vfcv, vfp, vfa, vfb, n, -normalVel,
		                          	VR, VL, dfdUR, betaRef, k1, cmach, irey, prec);

  for (int k=0; k<dim2; k++)
    dfdUR[k] = -dfdUR[k];
  
  int k;

  if (typeJac == FluxFcn::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) { 
      jacL[k] = dfdUL[k]; 
      jacR[k] = dfdUR[k]; 
    }
  }
  else {
    varFcn->postMultiplyBydUdV(VL, dfdUL, jacL, dflag);
    varFcn->postMultiplyBydUdV(VR, dfdUR, jacR, dflag);
  }

}

//------------------------------------------------------------------------------
void FluxFcnFDJacRoeEuler3D::computeBarotropicLiquid(double irey, double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{
    
   F77NAME(roeflux5waterdissprec)(0, gamma, vfCv, vfPr, vfa, vfb,
        normal, normalVel, VL, VL, VR, VR, flux, betaRef, k1, cmach, irey, prec);

}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeEuler3D::computeBarotropicLiquid(double irey, double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{
	 /*fprintf(stdout, "Vl[0] = %e\n", VL[0]);
	 fprintf(stdout, "Vl[1] = %e\n", VL[1]);
	 fprintf(stdout, "Vl[4] = %e\n", VL[4]);
	 fprintf(stdout, "Vr[0] = %e\n", VR[0]);
	 fprintf(stdout, "Vr[1] = %e\n", VR[1]);
	 fprintf(stdout, "Vr[4] = %e\n", VR[4]);
	 */
   F77NAME(roeflux5waterdissprec)(0, gamma, vfCv, vfPr, vfa, vfb,
        normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, betaRef, k1, cmach, irey, prec);
}

//------------------------------------------------------------------------------

void FluxFcnApprJacRoeEuler3D::computeJacobiansBarotropicLiquid(double irey, double vfCv, double vfPr,
						double vfa, double vfb, 
						double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, int flag)
{
  roejacappr3Dwater<5>(0, gamma, vf, vfCv, vfa, vfb, vfPr, type, normal, normalVel, VL, VR, jacL, jacR, irey, betaRef, k1, cmach, prec, flag);
}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux)
{
  fprintf(stderr, "no exact jacobian for Tait\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnExactJacRoeEuler3D::computeJacobiansBarotropicLiquid(double vfCv, double vfPr,
						 double vfa, double vfb, 
						 double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR)
{
  roejacexact3D<5>(0, gamma, vf, type, normal, normalVel, VL, VR, jacL, jacR);
}


//------------------------------------------------------------------------------

//NO VAN LEER

//------------------------------------------------------------------------------

void FluxFcnWallEuler3D::computeBarotropicLiquid(double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel,
            	                     double *V, double *Ub, double *flux)
{
  
  double P = vfPr + vfa*pow(V[0], vfb);
  flux[0] = 0.0;
  flux[1] = P * normal[0];
  flux[2] = P * normal[1];
  flux[3] = P * normal[2];
  flux[4] = P * normalVel;

}
//------------------------------------------------------------------------------

template<int dim>
inline
void flux3Dwater(int type, VarFcn *vf, double *normal, double normalVel, double *V, double *Ub, double *flux){
 
  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double VV[5];
  double machb = vf->computeMachNumber(Vb);
  double cb = vf->computeSoundSpeed(Vb);
  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2];//-nVel;
  double un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
  double c = vf->computeSoundSpeed(V);
  double rhoun, p, rhoE;
  
  if(un==0.0){                   // SLIP-WALL LIKE, if boundary velocity is in the same plane as the face of the boundary
    VV[0]=vf->getDensity(Vb);
    VV[1]=0.0; VV[2]=0.0; VV[3]=0.0; VV[4]=0.0;

    rhoun=0.0;
    p = vf->getPressure(VV);
    rhoE = 0.0;
  }else{
    if(un<0.0){ 		//INLET, as the normal is going outward.
      if(-un-c>0.0){ 			//SUPERSONIC
	VV[0]=vf->getDensity(Vb);
	VV[1]=Vb[1];
	VV[2]=Vb[2];
	VV[3]=Vb[3];
	VV[4]=vf->computeTemperature(Vb);
      }else{        			//SUBSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
      }
    }
    else{                  //OUTLET
      if(un-c>0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
      }else{                              //SUBSONIC
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
        VV[4] = vf->computeTemperature(V);
      }  
    }
    rhoun = VV[0] * ( VV[1]*n[0] + VV[2]*n[1] + VV[3]*n[2] - nVel );
    p = vf->getPressure(VV);
    rhoE = vf->computeRhoEnergy(VV);
  }
  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*VV[1] + p*n[0]);
  flux[2] = S * (rhoun*VV[2] + p*n[1]);
  flux[3] = S * (rhoun*VV[3] + p*n[2]);
  flux[4] = S * ((rhoE + p) * rhoun/VV[0] + p*nVel); 

}

//------------------------------------------------------------------------------



template<int dim>
inline
void jacflux3Dwater(int type, VarFcn *vf, FluxFcn::Type typeJac, double *normal,
		      double normalVel, double *V, double *Ub, double *jac){

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double VV[dim];
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2]; 
  double cb = vf->computeSoundSpeed(Vb);
  double un = V[1]*n[0]+V[2]*n[1]+V[3]*n[2];
  double c  = vf->computeSoundSpeed(V);

  double _dfdV[dim][dim];
  double *dfdV = reinterpret_cast<double *>(_dfdV);
  for(int k=0; k<dim*dim; ++k)
    dfdV[k] = 0.0;


  if(un == 0.0){                 //SLIP-WALL LIKE
    
    //nothing to do, the flux has only pressure terms
    //which depend on density only, which is the one 
    //at the boundary (not inside the volume)

  }else{  
    if(un < 0.0){ 		//INLET, as the normal is going outward.
      if(-un-c > 0.0){                  //SUPERSONIC
	
	//nothing to do: jacobian is null
	
      }else{                         //SUBSONIC
	
	//derivative wrt pressure
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	double unp = VV[1]*n[0] +VV[2]*n[1] + VV[3]*n[2] - nVel;
	_dfdV[0][0] = S * unp;
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[2][0] = S * (VV[2]*unp + cp2*n[1]);
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
      }
    }else{                  //OUTLET
      if(un-c > 0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	
	_dfdV[0][0] = S * unp;
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	_dfdV[2][0] = S * (VV[2]*unp +cp2*n[1]);
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
	
      }else{//subsonic outlet
	
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	
	
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
      }
    }
  }
  
  
  
  if (typeJac == FluxFcn::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<dim*dim; ++k)
      jac[k] = dfdV[k]; 
  
}

//------------------------------------------------------------------------------

void FluxFcnGhidagliaEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *V, double *Ub, double *flux)
{
  F77NAME(genbcfluxtait)(0, vfCv, vfPr, vfa, vfb, normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------

void FluxFcnInflowEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *V, double *Ub, double *flux)
{
    fprintf(stderr, "*** Error1: The computation of fluxes at the boundaries is not available for water\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnOutflowEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *V, double *Ub, double *flux)
{
    fprintf(stderr, "*** Error2: The computation of fluxes at the boundaries is not available for water\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *V, double *Ub, double *flux)
{
    flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);
}

//------------------------------------------------------------------------------

void FluxFcnInternalInflowEuler3D::computeJacobianBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                      double vfb, double *normal, double normalVel, 
				      double *V, double *Ub, double *jacL)
{
    jacflux3Dwater<5>(0,vf,type,normal,normalVel,V,Ub,jacL);
}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowEuler3D::computeBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
			             double *V, double *Ub, double *flux)
{   
    flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);
}

//------------------------------------------------------------------------------

void FluxFcnInternalOutflowEuler3D::computeJacobianBarotropicLiquid(double vfCv, double vfPr, double vfa,
                                     double vfb, double *normal, double normalVel, 
				     double *V, double *Ub, double *jacL)
{
    jacflux3Dwater<5>(0,vf,type,normal,normalVel,V,Ub,jacL);
}

//------------------------------------------------------------------------------
