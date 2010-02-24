#include <FluxFcnDescJWL.h>
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

void FluxFcnJWLFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux, int flag)
{

   fprintf(stderr, "FluxFcnJWLFDJacRoeEuler3D::compute not implemented\n");
   exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux, int flag)
{

  computeJWL(length, irey, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), 
             vf->getR2r(), normal, normalVel, VL, VR, flux);
 
}

//------------------------------------------------------------------------------

void FluxFcnJWLApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, int flag)
{
  computeJacobiansJWL(length, irey, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(),
             vf->getR2r(), normal, normalVel, VL, VR, jacL, jacR, flag);
}

//------------------------------------------------------------------------------

void FluxFcnJWLExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					double *VL, double *VR, double *flux, int flag)
{

  fprintf(stderr, "FluxFcnJWLExactJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR, int flag)
{

  fprintf(stderr, "FluxFcnJWLExactJacRoeEuler3D::computeJacobians not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *VL, double *VR, double *flux, int flag)
{

  computeJWL(normal, normalVel, VL, VR, flux);

}

//------------------------------------------------------------------------------

void FluxFcnJWLGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, int flag)
{

  computeJWL(vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(),
             vf->getR2r(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnJWLInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, int flag)
{
    fprintf(stderr, "*** Error1: StegerWarming not available for JWL\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux, int flag)
{
    fprintf(stderr, "*** Error2: StegerWarming not available for JWL\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJWLInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, int flag)
{

    fprintf(stderr, "*** Error: FluxFcnJWLInternalInflowEuler3D::compute not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, int flag)
{

    fprintf(stderr, "*** Error: FluxFcnJWLInternalInflowEuler3D::computeJacobian not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, int flag)
{

    fprintf(stderr, "*** Error: FluxFcnJWLInternalOutflowEuler3D::compute not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJWLInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, int flag)
{
 
    fprintf(stderr, "*** Error: FluxFcnJWLInternalOutflowEuler3D::computeJacobian not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

