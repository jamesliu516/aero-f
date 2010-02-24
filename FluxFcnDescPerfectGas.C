#include <FluxFcnDescPerfectGas.h>
//#include <FluxFcnDesc.C>
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

void FluxFcnPerfectGasFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasFDJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasApprJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasExactJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//-----------------------------------------------------------------------

void FluxFcnPerfectGasFDJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
}


//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLEEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
  
  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasVanLeerEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasVanLeerEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasVanLeerEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasWallEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnPerfectGasWallEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jac);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{

   computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{

   computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnPerfectGasInflowEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

   computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnPerfectGasOutflowEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
 
                                                                                                                  
void FluxFcnPerfectGasInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasInternalInflowEuler3D::computeDerivative(double vfgam, double vfp, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

                                                                                                                  
void FluxFcnPerfectGasInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasInternalOutflowEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
//turbulence
                                                                                                                  
                                                                                                                  
void FluxFcnPerfectGasFDJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasFDJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{
                                                                                                                  
  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  
                                                                                                                  
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasApprJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasExactJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                              double *VL, double *VR,
                                              double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacHLLESA3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLESA3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);


}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLESA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallSA3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasWallSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnPerfectGasWallSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jac);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasOutflowSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasInternalInflowSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                             double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasInternalOutflowSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                             double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                 double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnPerfectGasRoeSAturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnPerfectGasOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasFDJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasApprJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(irey, dIrey, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnPerfectGasExactJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VR, dVR, 0.0*dMach, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                              double *VL, double *VR,
                                              double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacHLLEKE3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLEKE3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacHLLEKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallKE3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnPerfectGasWallKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                   double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnPerfectGasWallKE3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jac);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasOutflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnPerfectGasOutflowKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel, double *V,
                                      double *Ub, double *dUb, double *flux, double *dFlux, int flag)
{

  computeDerivativeOfPerfectGas(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnPerfectGasRoeKEturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(length, irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnPerfectGasOutflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
