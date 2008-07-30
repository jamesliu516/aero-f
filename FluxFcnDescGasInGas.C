#include <FluxFcnDescGasInGas.h>
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

void FluxFcnGasInGasFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  fprintf(stdout, "Should never be called\n");
  exit(1);
  if ( flag == 1 )
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  else 
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnGasInGasApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasApprJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{


  if ( flag == 1 ) {
    computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }else{
    computeJacobiansPerfectGas(irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }                                                                                                                  
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-gas simulations\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasExactJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{
  fprintf(stderr, "no exact jacobian for gas-gas simulations\n");
  exit(1);
}

//-----------------------------------------------------------------------

void FluxFcnGasInGasFDJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);
  }

}

//-----------------------------------------------------------------------

void FluxFcnGasInGasApprJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);
  }
}

//-----------------------------------------------------------------------

void FluxFcnGasInGasApprJacHLLEEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
  if ( flag == 1 ) {
    computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }else{
    computeJacobiansPerfectGas(irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, jacL, jacR, flag);
  }
}


//-----------------------------------------------------------------------

void FluxFcnGasInGasFDJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);
  }

}

//-----------------------------------------------------------------------

void FluxFcnGasInGasApprJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  if ( flag == 1 ){
    computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
  else{
    computePerfectGas(length, irey, vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, VL, VR, flux);
  }
}

//-----------------------------------------------------------------------

void FluxFcnGasInGasVanLeerEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{
  fprintf(stderr, "the Van Leer fluxes are not implemented for a two-phase flow simulation\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasVanLeerEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
  fprintf(stderr, "the Van Leer fluxes are not implemented for a two-phase flow simulation\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(normal, normalVel, VL, VR, flux);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
   if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);
  else 
    computePerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, flux);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{
   if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);
  else 
    computePerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnGasInGasOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);
  else 
    computePerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
 
void FluxFcnGasInGasInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{

  if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);
  else 
    computePerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnGasInGasInternalInflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);
  else
    computeJacobianPerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, jacL, flag);
}

//------------------------------------------------------------------------------

void FluxFcnGasInGasInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{
  if ( flag == 1 )
    computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub,  flux, flag);
  else 
    computePerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnGasInGasInternalOutflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{
  if ( flag == 1 )
    computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);
  else
    computeJacobianPerfectGas(vf->getGammabis(), vf->getPressureConstantbis(), normal, normalVel, V, Ub, jacL, flag);

}
