#include <FluxFcnDescWaterCompressible.h>
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

void FluxFcnWaterCompressibleFDJacRoeEuler3D::compute(double irey, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux, int flag)
{

   fprintf(stderr, "Should never be called\n");
   exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleApprJacRoeEuler3D::compute(double irey, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux, int flag)
{

  computeBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);
 
}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleApprJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, int flag)
{
  computeJacobiansBarotropicLiquid(irey, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, jacL, jacR, flag);
}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleExactJacRoeEuler3D::compute(double irey, double *normal, double normalVel, 
					double *VL, double *VR, double *flux, int flag)
{

  fprintf(stderr, "no exact jacobian for Tait\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleExactJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR, int flag)
{

  fprintf(stderr, "no exact jacobian for Tait\n");
  exit(1);

}

//------------------------------------------------------------------------------

//NO VAN LEER

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleWallEuler3D::compute(double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{

  computeBarotropicLiquid(vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleGhidagliaEuler3D::compute(double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, int flag)
{

  computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                          normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleInflowEuler3D::compute(double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, int flag)
{
    fprintf(stderr, "*** Error1: StegerWarming not available for Tait\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleOutflowEuler3D::compute(double irey, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux, int flag)
{
    fprintf(stderr, "*** Error2: StegerWarming not available for Tait\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleInternalInflowEuler3D::compute(double irey, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, int flag)
{

  computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleInternalInflowEuler3D::computeJacobian(double irey, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleInternalOutflowEuler3D::compute(double irey, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, int flag)
{

  computeBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnWaterCompressibleInternalOutflowEuler3D::computeJacobian(double irey, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, int flag)
{
 
  computeJacobianBarotropicLiquid(vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

