#include <FluxFcnDescJWLInGas.h>
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

void FluxFcnJWLInGasFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  fprintf(stderr, "*** Error : FluxFcnJWLInGasFDJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computeJWL(length, irey, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
  if ( flag == 1 ) {
    computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }else{
    computeJacobiansJWL(length, irey, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }                                                                                                                  
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-gas simulations\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-gas simulations\n");
  exit(1);
}

//-----------------------------------------------------------------------

void FluxFcnJWLInGasVanLeerEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{
  fprintf(stderr, "the Van Leer fluxes are not implemented for a two-phase flow simulation\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasVanLeerEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
  fprintf(stderr, "the Van Leer fluxes are not implemented for a two-phase flow simulation\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(normal, normalVel, VL, VR, flux);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
   if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);
  else 
    computeJWL(vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(), normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasInflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasOutflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------
 
void FluxFcnJWLInGasInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasInternalInflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasInternalInflowEuler3D::computeJacobian not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasInternalOutflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInGasInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{
  fprintf(stderr, "*** Error : FluxFcnJWLInGasInternalOutflowEuler3D::computeJacobian not implemented\n");
  exit(1);
}
