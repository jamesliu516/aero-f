#include <FluxFcnDescPerfectGas.h>
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


/* Note:
** flag is a dummy index for this class.
** It is only used for multiphase flows.
** It does not need to be specified for single phase flows.
*/

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, int flag)
{
  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeEuler3D::compute(double irey, double *normal, double normalVel,
                                        double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                                 double *VL, double *VR,
                                                 double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//-----------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasVanLeerEuler3D::compute(double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasVanLeerEuler3D::computeJacobians(double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallEuler3D::compute(double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasGhidagliaEuler3D::compute(double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{

   computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInflowEuler3D::compute(double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, int flag)
{

   computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasOutflowEuler3D::compute(double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
 
                                                                                                                  
void FluxFcnPerfectGasInternalInflowEuler3D::compute(double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

                                                                                                                  
void FluxFcnPerfectGasInternalOutflowEuler3D::compute(double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowEuler3D::computeJacobian(double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
//turbulence
                                                                                                                  
                                                                                                                  
                                                                                                                  
void FluxFcnPerfectGasFDJacRoeSA3D::compute(double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeSA3D::compute(double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{
                                                                                                                  
  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeSA3D::computeJacobians(double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeSA3D::compute(double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeSA3D::computeJacobians(double irey, double *normal, double normalVel,
                                              double *VL, double *VR,
                                              double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasWallSA3D::compute(double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasOutflowSA3D::compute(double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSA3D::compute(double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSA3D::computeJacobian(double irey, double *normal, double normalVel,
                                                double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSA3D::compute(double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSA3D::computeJacobian(double irey, double *normal, double normalVel,
                                                 double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnPerfectGasRoeSAturb3D::computeJacobians(double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallSAturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnPerfectGasOutflowSAturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalInflowSAturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasInternalOutflowSAturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL, flag);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasFDJacRoeKE3D::compute(double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasApprJacRoeKE3D::compute(double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasApprJacRoeKE3D::computeJacobians(double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, int flag)
{
                                                                                                                  
  computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, flag);
                                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeKE3D::compute(double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasExactJacRoeKE3D::computeJacobians(double irey, double *normal, double normalVel,
                                              double *VL, double *VR,
                                              double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasWallKE3D::compute(double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnPerfectGasOutflowKE3D::compute(double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, int flag)
{

  computePerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnPerfectGasRoeKEturb3D::computeJacobians(double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, int flag)
{

  computeJacobiansPerfectGas(irey, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR);

}

//------------------------------------------------------------------------------

void FluxFcnPerfectGasWallKEturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnPerfectGasOutflowKEturb3D::computeJacobian(double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, int flag)
{

  computeJacobianPerfectGas(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
