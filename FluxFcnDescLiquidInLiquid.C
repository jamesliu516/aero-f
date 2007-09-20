#include <FluxFcnDescLiquidInLiquid.h>
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

void FluxFcnLiquidInLiquidFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

   fprintf(stderr, "Should never be called\n");
   exit(1);
}

//------------------------------------------------------------------------------
void FluxFcnLiquidInLiquidApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);
  
  else
    computeBarotropicLiquid(irey, vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, VL, VR, flux);
  
}

//------------------------------------------------------------------------------
void FluxFcnLiquidInLiquidApprJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
  for(int j=0; j<25; j++){
    jacL[j] = 0.0;
    jacR[j] = 0.0;
  }


  if ( flag == 1 ){ 
    computeJacobiansBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }else
    computeJacobiansBarotropicLiquid(irey, vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, VL, VR, jacL, jacR, flag);
  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnLiquidInLiquidExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{

  fprintf(stderr, "no exact jacobian for liquid-liquid simulations\n");
  exit(1);
}

//------------------------------------------------------------------------------
void FluxFcnLiquidInLiquidExactJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{

  fprintf(stderr, "no exact jacobian for liquid-liquid simulations\n");
  exit(1);
}

//-----------------------------------------------------------------------
void FluxFcnLiquidInLiquidWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1)
    computeBarotropicLiquid(vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);
  else
    computeBarotropicLiquid(vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, VL, VR, flux);
  
}

//------------------------------------------------------------------------------

void FluxFcnLiquidInLiquidGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                            normal, normalVel, V, Ub, flux);
  else
    computeBarotropicLiquid(vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(),
                            normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------
void FluxFcnLiquidInLiquidInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error2: StegerWarming not available for Tait\n");
  exit(1);
}

//------------------------------------------------------------------------------
void FluxFcnLiquidInLiquidOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error2: StegerWarming not available for Tait\n");
  exit(1);
}

//------------------------------------------------------------------------------
 
void FluxFcnLiquidInLiquidInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);
  else 
    computeBarotropicLiquid(vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------

void FluxFcnLiquidInLiquidInternalInflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);
  else
    computeJacobianBarotropicLiquid(vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, V, Ub, jacL);
}

//------------------------------------------------------------------------------

void FluxFcnLiquidInLiquidInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub,  flux);
  else 
    computeBarotropicLiquid(vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, V, Ub, flux);
                                                                                                                  
}

//------------------------------------------------------------------------------

void FluxFcnLiquidInLiquidInternalOutflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);
  else
    computeJacobianBarotropicLiquid(vf->getCvbis(), vf->getPrefWaterbis(), vf->getAlphaWaterbis(), vf->getBetaWaterbis(), normal, normalVel, V, Ub, jacL);
}
