#include <FluxFcnDescGasInLiquid.h>
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

//------------------------------------------------------------------------------

void FluxFcnGasInLiquidFDJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  fprintf(stderr, "Should never be called\n");
   exit(1);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnGasInLiquidApprJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
                       
  if ( flag == 1 ){
     computeBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnGasInLiquidApprJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
  if ( flag == 1 ) {
    computeJacobiansBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }else{
    computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag); 
  }                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnGasInLiquidExactJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-liquid simulations\n");
  exit(1);                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnGasInLiquidExactJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-liquid simulations\n");
  exit(1);
}

//-----------------------------------------------------------------------

void FluxFcnGasInLiquidWallEuler3D::compute(double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{
  
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);
  else
    computePerfectGas(normal, normalVel, VL, VR, flux);
     
}

//------------------------------------------------------------------------------

void FluxFcnGasInLiquidGhidagliaEuler3D::compute(double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                            normal, normalVel, V, Ub, flux);
  else
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------
                                                                                                                  
                                                                                                                  
void FluxFcnGasInLiquidInflowEuler3D::compute(double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
                                                                                                                  
  fprintf(stderr, "*** Error1: StegerWarming not available for Tait\n");
  exit(1);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnGasInLiquidOutflowEuler3D::compute(double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error1: StegerWarming not available for Tait\n");
  exit(1);
                                                                                                                  
}

//------------------------------------------------------------------------------
 
                                                                                                                  
void FluxFcnGasInLiquidInternalInflowEuler3D::compute(double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);
  else 
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);
}

//------------------------------------------------------------------------------

void FluxFcnGasInLiquidInternalInflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);
  else
    computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);
}

//------------------------------------------------------------------------------

                                                                                                                  
void FluxFcnGasInLiquidInternalOutflowEuler3D::compute(double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub,  flux);
  else 
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);
}

//------------------------------------------------------------------------------

void FluxFcnGasInLiquidInternalOutflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);
  else
    computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);
}
