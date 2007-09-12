#include <BcFcn.h>
#include <BcDef.h>
#include <IoData.h>

#include <stdlib.h>
#include <stdio.h>

//------------------------------------------------------------------------------

void BcFcn::applyToSolutionVector(int t, double *v, double *u)
{

  fprintf(stderr, "*** Error: applyToSolutionVector function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToResidualTerm(int t, double *v, double *u, double *f)
{

  fprintf(stderr, "*** Error: applyToResidualTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, float *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, double *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, float *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, double *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------


void BcFcn::zeroDiagonalTerm(int t, float *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::zeroDiagonalTerm(int t, double *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::zeroDiagonalTerm(int t, bcomp *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

inline
void BcFcnNS::template_applyToSolutionVectorTerm(int type, double *Vwall, double *U)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    U[1] = U[0] * Vwall[1];
    U[2] = U[0] * Vwall[2];
    U[3] = U[0] * Vwall[3];

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
      U[4] = U[0] * ( Vwall[4] + 0.5*(Vwall[1]*Vwall[1]+Vwall[2]*Vwall[2]+Vwall[3]*Vwall[3]) );
  }

}

//------------------------------------------------------------------------------

inline
void BcFcnNS::template_applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    F[1] = - Vwall[1] * U[0] + U[1];
    F[2] = - Vwall[2] * U[0] + U[2];
    F[3] = - Vwall[3] * U[0] + U[3];

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
      F[4] = (U[4] - U[0] * Vwall[4]) * U[0] - 0.5 * (U[1]*U[1] + U[2]*U[2] + U[3]*U[3]);
  }
  
#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    F[0] = - Vwall[0] + U[0];
    F[1] = - Vwall[0]*Vwall[1] + U[1];
    F[2] = - Vwall[0]*Vwall[2] + U[2];
    F[3] = - Vwall[0]*Vwall[3] + U[3];
  }
#endif

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToSolutionVector(int type, double *Vwall, double *U)
{

  template_applyToSolutionVectorTerm(type, Vwall, U);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  template_applyToResidualTerm(type, Vwall, U, F);  

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm<float,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm<double,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm<bcomp,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, float *A)
{

  template_zeroDiagonalTerm<float,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, double *A)
{

  template_zeroDiagonalTerm<double,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, bcomp *A)
{

  template_zeroDiagonalTerm<bcomp,5>(type, A);

}

void BcFcnNS::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm<float,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm<double,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm<bcomp,5>(type, A);

}

//------------------------------------------------------------------------------

BcFcnSA::BcFcnSA(IoData& iod)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = true;
  else
    wallFcn = false;

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToSolutionVector(int type, double *Vwall, double *U)
{
  if (!wallFcn)
    BcFcnNS::template_applyToSolutionVectorTerm(type, Vwall, U);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    U[5] = U[0] * Vwall[5];

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (!wallFcn)
    BcFcnNS::template_applyToResidualTerm(type, Vwall, U, F);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    F[5] = - U[0] * Vwall[5] + U[5];

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToSolutionVector(int type, double *Vwall, double *U)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    U[5] = Vwall[5];


}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    F[5] = - Vwall[5] * U[0] + U[5];
    F[6] = - Vwall[6] * U[0] + U[6];
  }

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------
