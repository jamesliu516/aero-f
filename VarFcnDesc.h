#ifndef _VAR_FCN_DESC_H_
#define _VAR_FCN_DESC_H_

#include <VarFcn.h>

//------------------------------------------------------------------------------
class VarFcnPerfectGasEuler3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasEuler3D(IoData &);
  ~VarFcnPerfectGasEuler3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0);
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0);

};

//------------------------------------------------------------------------------

inline
VarFcnPerfectGasEuler3D::VarFcnPerfectGasEuler3D(IoData &iod) : VarFcnPerfectGas(iod)
{
  pname = new const char*[5];

  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasEuler3D::conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV, double phi)
{
  conservativeToPrimitiveDerivativeGasEuler(gam, Pstiff, dPstiff, U, dU, V, dV);
}

//------------------------------------------------------------------------------

inline
int VarFcnPerfectGasEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
  return VerificationGasEuler(glob, pmin, gam, Pstiff, U, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasEuler3D::primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU, double phi)
{
  primitiveToConservativeDerivativeGasEuler(gam, invgam1, Pstiff, dPstiff, V, dV, U, dU);
}

//------------------------------------------------------------------------------    

inline
void VarFcnPerfectGasEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUGasEuler(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdUGasEuler(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUTGasEuler(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdUTGasEuler(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  preMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydVdUGasEuler(gam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                  double *Vinter, double *V, double phi)
{
  extrapolatePrimitiveGasEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  extrapolateCharacteristicGasEuler(gam,Pstiff,n,un,c,Vb,dV);
}

//------------------------------------------------------------------------------

class VarFcnPerfectGasSA3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasSA3D(IoData &);
  ~VarFcnPerfectGasSA3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  double getTurbulentNuTilde(double *V) { return V[5]; }

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0);
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0);

};

//------------------------------------------------------------------------------

inline
VarFcnPerfectGasSA3D::VarFcnPerfectGasSA3D(IoData &iod) : VarFcnPerfectGas(iod)
{
  pname = new const char*[6];

  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
  pname[5] = "nut";
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  conservativeToPrimitiveGasSA(gam1, U, V);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasSA3D::conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV, double phi)
{
   conservativeToPrimitiveDerivativeGasSA(gam1, U, dU, V, dV);
}

//------------------------------------------------------------------------------

inline
int VarFcnPerfectGasSA3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  conservativeToPrimitiveGasSA(gam1, U, V);
  return VerificationGasSA(glob, pmin, gam1, U, V);
}

//------------------------------------------------------------------------------

// Modified (MB)
inline
void VarFcnPerfectGasSA3D::primitiveToConservative(double *V, double *U, double phi)
{
  primitiveToConservativeGasSA(invgam1, V, U);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasSA3D::primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU, double phi)
{
  primitiveToConservativeDerivativeGasSA(invgam1, V, dV, U, dU);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
// Modified (MB)
  multiplyBydVdUGasSA(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  preMultiplyBydUdVGasSA(invgam1, V, mat, res);
}
 
//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydVdUGasSA(gam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydUdVGasSA(invgam1, V, mat, res);
}

//------------------------------------------------------------------------------

class VarFcnPerfectGasKE3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasKE3D(IoData &);
  ~VarFcnPerfectGasKE3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  double getTurbulentKineticEnergy(double *V) { return V[5]; }
  double getTurbulentDissipationRate(double *V) { return V[6]; }

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0);
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0);

};

//------------------------------------------------------------------------------

inline
VarFcnPerfectGasKE3D::VarFcnPerfectGasKE3D(IoData &iod) : VarFcnPerfectGas(iod)
{
  pname = new const char*[7];

  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
  pname[5] = "k";
  pname[6] = "epsilon";
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  conservativeToPrimitiveGasKE(gam1, U, V);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasKE3D::conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV, double phi)
{
   conservativeToPrimitiveDerivativeGasKE(gam1, U, dU, V, dV);
}

//------------------------------------------------------------------------------

inline
int VarFcnPerfectGasKE3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  conservativeToPrimitiveGasKE(gam1, U, V);
  return VerificationGasKE(glob, pmin, gam1, U, V);
}

//------------------------------------------------------------------------------

// Modified (MB)
inline
void VarFcnPerfectGasKE3D::primitiveToConservative(double *V, double *U, double phi)
{
  primitiveToConservativeGasKE(invgam1, V, U);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void VarFcnPerfectGasKE3D::primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU, double phi)
{
  primitiveToConservativeDerivativeGasKE(invgam1, V, dV, U, dU);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
// Modified (MB)
  multiplyBydVdUGasKE(gam1, V, vec, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  preMultiplyBydUdVGasKE(invgam1, V, mat, res);
}

//------------------------------------------------------------------------------
inline
void VarFcnPerfectGasKE3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydVdUGasKE(gam1, V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydUdVGasKE(invgam1, V, mat, res);
}

//-------------------------------------------------------------------------------
class VarFcnWaterCompressibleEuler3D : public VarFcnWaterCompressible {
public:
  VarFcnWaterCompressibleEuler3D(IoData &);
  ~VarFcnWaterCompressibleEuler3D() {}
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);
  void computeNewPrimitive(double *, double *);
  void computeOldPrimitive(double *, double *);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//--------------------------------------------------------------------------
inline
VarFcnWaterCompressibleEuler3D::VarFcnWaterCompressibleEuler3D(IoData &iod) : VarFcnWaterCompressible(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnWaterCompressibleEuler3D::computeNewPrimitive(double *U, double *V)
{
  computeNewPrimitiveLiquidEuler(Pref_water,alpha_water,beta_water,U,V);
}
//------------------------------------------------------------------------------
inline
void VarFcnWaterCompressibleEuler3D::computeOldPrimitive(double *U, double *V)
{
  computeOldPrimitiveLiquidEuler(Pref_water,alpha_water,beta_water,U,V);
}
//------------------------------------------------------------------------------

inline
void VarFcnWaterCompressibleEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  conservativeToPrimitiveLiquidEuler(invCv, U, V);
}

//------------------------------------------------------------------------------

inline
int VarFcnWaterCompressibleEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  conservativeToPrimitiveLiquidEuler(invCv, U, V);
  return VerificationLiquidEuler(glob, pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
}

//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  primitiveToConservativeLiquidEuler(Cv, V, U);
}
                                                                                         
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdULiquidEuler(invCv, V, vec, res);
}

//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdULiquidEuler(invCv, V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnWaterCompressibleEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  preMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydVdULiquidEuler(invCv, V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------


inline
void VarFcnWaterCompressibleEuler3D::extrapolateBoundaryPrimitive(double un, double c,
                                                 double *Vb, double *Vinter, double *V,
                                                 double phi)
{
  extrapolatePrimitiveLiquidEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnWaterCompressibleEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  extrapolateCharacteristicLiquidEuler(Cv,Pref_water,alpha_water,beta_water,
					 n,un,c,Vb,dV);
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
class VarFcnJWLEuler3D : public VarFcnJWL {
public:
  VarFcnJWLEuler3D(IoData &);
  ~VarFcnJWLEuler3D() {}
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//--------------------------------------------------------------------------
inline
VarFcnJWLEuler3D::VarFcnJWLEuler3D(IoData &iod) : VarFcnJWL(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------

inline
void VarFcnJWLEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  conservativeToPrimitiveJWLEuler(omega, computeFrho(U), U, V);
}

//------------------------------------------------------------------------------

inline
int VarFcnJWLEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  conservativeToPrimitiveJWLEuler(omega, computeFrho(U), U, V);
  return VerificationJWLEuler(glob, pmin, invomega, computeFrho(U), U, V);
}

//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  primitiveToConservativeJWLEuler(invomega, computeFrho(V), V, U);
}
                                                                                         
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUJWLEuler(omega, computeFrhop(V), V, vec, res);
}

//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdUJWLEuler(omega, computeFrhop(V), V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUTJWLEuler(omega, computeFrhop(V), V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnJWLEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  multiplyBydVdUTJWLEuler(omega, computeFrhop(V), V, vec, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  preMultiplyBydUdVJWLEuler(invomega, computeFrhop(V), V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydVdUJWLEuler(omega, computeFrhop(V), V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V), V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnJWLEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V), V, mat, res);
}
                                                                                            
//------------------------------------------------------------------------------


inline
void VarFcnJWLEuler3D::extrapolateBoundaryPrimitive(double un, double c,
                                                 double *Vb, double *Vinter, double *V,
                                                 double phi)
{
  extrapolatePrimitiveJWLEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnJWLEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  extrapolateCharacteristicJWLEuler(un,c,n,Vb,dV);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

class VarFcnGasInGasEuler3D : public VarFcnGasInGas {

public:

  VarFcnGasInGasEuler3D(IoData &);
  ~VarFcnGasInGasEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0);
                                                                                            
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//------------------------------------------------------------------------------
                                                                                           
inline
VarFcnGasInGasEuler3D::VarFcnGasInGasEuler3D(IoData &iod) : VarFcnGasInGas(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  if (phi>=0.0)
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
  else  conservativeToPrimitiveGasEuler(gamp, Pstiffp, U, V);

}
//------------------------------------------------------------------------------
inline
int VarFcnGasInGasEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
    return VerificationGasEuler(glob, pmin, gam, Pstiff, U, V);
  }else{
    conservativeToPrimitiveGasEuler(gamp, Pstiffp, U, V);
    return VerificationGasEuler(glob, pminp, gamp, Pstiffp, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  if (phi>=0.0)
    primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
  else  primitiveToConservativeGasEuler(gamp, invgamp1, Pstiffp, V, U);
}
//------------------------------------------------------------------------------    
inline
bool VarFcnGasInGasEuler3D::updatePhaseChange(double *V, double *U, double phi,
                            double phin, double *Riemann, double weight)
{

  //nature of fluid at this node has not changed over time
  if(phi*phin > 0.0){
    if(phi>=0.0)
      primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    else primitiveToConservativeGasEuler(gamp, invgamp1, Pstiffp, V, U);
    return false;

  //nature of fluid at this node has changed over time
  }else{
    if(weight<=0.0){
      fprintf(stdout, "negative weight in updatePhaseChange\n");
      exit(1);
    }
    for(int k=0; k<5; k++)
    //for(int k=0; k<1; k++)
      V[k] = Riemann[k]/weight;

    // from Fluid2 to Fluid1, ie Gas To Gas
    if(phi>=0.0 && phin<0.0)
      primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);

    // from Fluid1 to Fluid2, ie Gas To Gas
    else primitiveToConservativeGasEuler(gamp, invgamp1, Pstiffp, V, U);

    return true;
  }

}
//------------------------------------------------------------------------------    
inline
void VarFcnGasInGasEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUGasEuler(gamp1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUGasEuler(gamp1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUTGasEuler(gamp1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUTGasEuler(gamp1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    preMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  preMultiplyBydUdVGasEuler(invgamp1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydVdUGasEuler(gam1, V, mat, res);
  else  postMultiplyBydVdUGasEuler(gamp1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  postMultiplyBydUdVGasEuler(invgamp1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  postMultiplyBydUdVGasEuler(invgamp1, V, mat, res);
}

//----------------------------------------------------------------------------------

inline
void VarFcnGasInGasEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                double *Vinter, double *V, double phi)
{
  extrapolatePrimitiveGasEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnGasInGasEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  if(phi>=0.0)
    extrapolateCharacteristicGasEuler(gam,Pstiff,n,un,c,Vb,dV);
  else
    extrapolateCharacteristicGasEuler(gamp,Pstiffp,n,un,c,Vb,dV);
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------                                                                                           
class VarFcnLiquidInLiquidEuler3D : public VarFcnLiquidInLiquid {

public:
  VarFcnLiquidInLiquidEuler3D(IoData &);
  ~VarFcnLiquidInLiquidEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int  conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0);

  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//------------------------------------------------------------------------------
                                                                                           
inline
VarFcnLiquidInLiquidEuler3D::VarFcnLiquidInLiquidEuler3D(IoData &iod) : VarFcnLiquidInLiquid(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  if (phi>=0.0)
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
  else  conservativeToPrimitiveLiquidEuler(invCvbis, U, V);

}
//------------------------------------------------------------------------------
inline
int VarFcnLiquidInLiquidEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
    return VerificationLiquidEuler(glob, pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
  }else{
    conservativeToPrimitiveLiquidEuler(invCvbis, U, V);
    return VerificationLiquidEuler(glob, pminp, alpha_waterbis, beta_waterbis, Pref_waterbis, invCvbis, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  if (phi>=0.0)
    primitiveToConservativeLiquidEuler(Cv, V, U);
  else  primitiveToConservativeLiquidEuler(Cvbis, V, U);
}
//------------------------------------------------------------------------------    
inline
bool VarFcnLiquidInLiquidEuler3D::updatePhaseChange(double *V, double *U, double phi,
                                  double phin, double *Riemann, double weight)
{

  //nature of fluid at this node has not changed over time
  if(phi*phin > 0.0){
    if(phi>=0.0)
      primitiveToConservativeLiquidEuler(Cv, V, U);
    else primitiveToConservativeLiquidEuler(Cvbis, V, U);

  //nature of fluid at this node has changed over time
  }else{
    for(int k=0; k<5; k++)
      V[k] = Riemann[k]/weight;

    // from Fluid2 to Fluid1, ie Liquid to Liquid
    if(phi>=0.0 && phin<0.0)
      primitiveToConservativeLiquidEuler(Cv, V, U);

    // from Fluid1 to Fluid2, ie Liquid to Liquid
    else primitiveToConservativeLiquidEuler(Cvbis, V, U);

  }

}
//------------------------------------------------------------------------------    
inline
void VarFcnLiquidInLiquidEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdULiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdULiquidEuler(invCvbis, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdULiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdULiquidEuler(invCv, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUTLiquidEuler(invCvbis, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUTLiquidEuler(invCvbis, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    preMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  preMultiplyBydUdVLiquidEuler(Cvbis, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydVdULiquidEuler(invCv, V, mat, res);
  else  postMultiplyBydVdULiquidEuler(invCvbis, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  postMultiplyBydUdVLiquidEuler(Cvbis, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  postMultiplyBydUdVLiquidEuler(Cvbis, V, mat, res);
}

//----------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                double *Vinter, double *V, double phi)
{
  extrapolatePrimitiveLiquidEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnLiquidInLiquidEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  if(phi>=0.0)
    extrapolateCharacteristicLiquidEuler(Cv,Pref_water,alpha_water,beta_water,n,un,c,Vb,dV);
  else
    extrapolateCharacteristicLiquidEuler(Cvbis,Pref_waterbis,alpha_waterbis,beta_waterbis,n,un,c,Vb,dV);
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------                                                                                           

class VarFcnGasInLiquidEuler3D : public VarFcnGasInLiquid {

public:
  VarFcnGasInLiquidEuler3D(IoData &);
  ~VarFcnGasInLiquidEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int  conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0);

  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//------------------------------------------------------------------------------
                                                                                           
inline
VarFcnGasInLiquidEuler3D::VarFcnGasInLiquidEuler3D(IoData &iod) : VarFcnGasInLiquid(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature_or_pressure";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  if (phi>=0.0)
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
  else  conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);

}
//------------------------------------------------------------------------------
inline
int VarFcnGasInLiquidEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
    return VerificationLiquidEuler(glob, pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
  }else{
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
    return VerificationGasEuler(glob, pminp, gam, Pstiff, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::primitiveToConservative(double *V, double *U, double phi)
{

  if (phi>=0.0)
    primitiveToConservativeLiquidEuler(Cv, V, U);
  else  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);

}
//------------------------------------------------------------------------------    
inline
bool VarFcnGasInLiquidEuler3D::updatePhaseChange(double *V, double *U, double phi,
                               double phin, double *Riemann, double weight)
{

  //nature of fluid at this node has not changed over time
  if(phi*phin > 0.0){
    if(phi>=0.0)
      primitiveToConservativeLiquidEuler(Cv, V, U);
    else primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    return false;

  //nature of fluid at this node has changed over time
  }else{

    if(weight<=0.0){
      fprintf(stdout, "negative weight in updatePhaseChange\n");
      exit(1);
    }
    // from Fluid2 to Fluid1, ie Gas To Liquid
    if(phi>=0.0 && phin<0.0){
      fprintf(stdout, "V1 = %e %e %e\n", V[0],V[1],V[4]);
      for(int k=0; k<5; k++)
        V[k] = Riemann[k]/weight;
      fprintf(stdout, "V2 = %e %e %e\n", V[0],V[1],V[4]);
      primitiveToConservativeLiquidEuler(Cv,V,U);

    // from Fluid1 to Fluid2, ie Liquid To Gas
    }else{
      V[4] = getPressure(V,phin);
      primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);

      //for(int k=0; k<5; k++)
      //  V[k] = Riemann[k]/weight;
      //primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    }
    return true;
  }

}
//------------------------------------------------------------------------------
/*
inline
void VarFcnGasInLiquidEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf, int node)
{
  if (!phi1) {
    if (phi>=0.0)
      primitiveToConservativeLiquidEuler(Cv, V, U);
    else  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
  }

  else {
    if (phi * (*phi1) < 0.0) { // the fluid at this node has changed
      if (phi >= 0.0 && (*phi1) < 0.0) { // gas to liquid

        if(vgf[0]<0.0) fprintf(stdout, "*** Warning: interface jumped more than 2 nodes in 1 time step\n");
        if(node_change) {
          fprintf(stdout, "local node %d changed from gas to liquid\n", node);
          fprintf(stdout, "gaspressure = %e and density = %e", V[4], V[0]);
        }


        V[0]   = vgf[0];
        V[1]   = vgf[1];
        V[2]   = vgf[2];
        V[3]   = vgf[3];
        V[4]   = vgf[4];

        if(node_change) fprintf(stdout, " become %e and %e\n", getPressure(V, phi),V[0]);

        primitiveToConservativeLiquidEuler(Cv, V, U);
      }
      if (phi < 0.0 && (*phi1) >= 0.0) { // liquid to gas
        if(node_change) fprintf(stdout, "node changed from liquid to gas\n");
        double P_g  = getPressure(V, *phi1); 
        V[4] = P_g;
        primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
      }
    }
    else {  // the fluid at this node did not change
      if (phi>=0.0)
        primitiveToConservativeLiquidEuler(Cv, V, U);
      else  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    }
  }
}
*/
//------------------------------------------------------------------------------    
inline
void VarFcnGasInLiquidEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdULiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUGasEuler(gam1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdULiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUGasEuler(gam1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUTGasEuler(gam1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTLiquidEuler(invCv, V, vec, res);
  else  multiplyBydVdUTGasEuler(gam1, V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    preMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  preMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydVdULiquidEuler(invCv, V, mat, res);
  else  postMultiplyBydVdUGasEuler(gam1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVLiquidEuler(Cv, V, mat, res);
  else  postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
}

//----------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                double *Vinter, double *V, double phi)
{
  if (phi>=0.0)
    extrapolatePrimitiveLiquidEuler(un, c, Vb, Vinter, V);
  else
    extrapolatePrimitiveGasEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnGasInLiquidEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  if (phi>=0.0)
    extrapolateCharacteristicLiquidEuler(Cv,Pref_water,alpha_water,beta_water,n,un,c,Vb,dV);
  else
    extrapolateCharacteristicGasEuler(gam,Pstiff,n,un,c,Vb,dV);
}
//------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

class VarFcnJWLInGasEuler3D : public VarFcnJWLInGas {

public:

  VarFcnJWLInGasEuler3D(IoData &);
  ~VarFcnJWLInGasEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0);
                                                                                            
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//------------------------------------------------------------------------------
                                                                                           
inline
VarFcnJWLInGasEuler3D::VarFcnJWLInGasEuler3D(IoData &iod) : VarFcnJWLInGas(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  if (phi>=0.0)
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
  else  conservativeToPrimitiveJWLEuler(omega, computeFrho(U,phi), U, V);

}
//------------------------------------------------------------------------------
inline
int VarFcnJWLInGasEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
    return VerificationGasEuler(glob, pmin, gam, Pstiff, U, V);
  }else{
    conservativeToPrimitiveJWLEuler(omega, computeFrho(U,phi), U, V);
    return VerificationJWLEuler(glob, pminp, invomega, computeFrho(U,phi), U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  if (phi>=0.0)
    primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
  else  primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);
}
//------------------------------------------------------------------------------    
inline
bool VarFcnJWLInGasEuler3D::updatePhaseChange(double *V, double *U, double phi,
                            double phin, double *Riemann, double weight)
{

  //nature of fluid at this node has not changed over time
  if(phi*phin > 0.0){
    if(phi>=0.0)
      primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    else primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);

  //nature of fluid at this node has changed over time
  }else{
    if(weight<=0.0){
      fprintf(stdout, "negative weight in updatePhaseChange\n");
      exit(1);
    }
    for(int k=0; k<5; k++)
    //for(int k=0; k<1; k++)
      V[k] = Riemann[k]/weight;

    // from Fluid2 to Fluid1, ie JWL To Gas
    if(phi>=0.0 && phin<0.0)
      primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);

    // from Fluid1 to Fluid2, ie Gas To JWL
    else primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);

  }

}
//------------------------------------------------------------------------------    
inline
void VarFcnJWLInGasEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUTJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTGasEuler(gam1, V, vec, res);
  else  multiplyBydVdUTJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    preMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  preMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydVdUGasEuler(gam1, V, mat, res);
  else  postMultiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInGasEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVGasEuler(invgam1, V, mat, res);
  else  postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
}

//----------------------------------------------------------------------------------

inline
void VarFcnJWLInGasEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                double *Vinter, double *V, double phi)
{
  if(phi>=0.0)
    extrapolatePrimitiveGasEuler(un, c, Vb, Vinter, V);
  else
    extrapolatePrimitiveJWLEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnJWLInGasEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  if(phi>=0.0)
    extrapolateCharacteristicGasEuler(gam,Pstiff,n,un,c,Vb,dV);
  else
    extrapolateCharacteristicJWLEuler(un,c,n,Vb,dV);
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

class VarFcnJWLInJWLEuler3D : public VarFcnJWLInJWL {

public:

  VarFcnJWLInJWLEuler3D(IoData &);
  ~VarFcnJWLInJWLEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  int conservativeToPrimitiveVerification(int, double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0);

  bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0);
                                                                                            
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0);
  void multiplyBydVdUT(double *, double *, double *, double = 0.0);
  void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0);
                                                                                            
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
                                                                                            
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0);

  void extrapolateBoundaryPrimitive(double, double, double*, double*, double*, double = 0.0);
  void extrapolateBoundaryCharacteristic(double*, double, double, double*, double*, double = 0.0);

// Included (MB)
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) {}
  void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) {}

};
                                                                                            
//------------------------------------------------------------------------------
                                                                                           
inline
VarFcnJWLInJWLEuler3D::VarFcnJWLInJWLEuler3D(IoData &iod) : VarFcnJWLInJWL(iod)
{
                                                                                            
  pname = new const char*[5];
                                                                                            
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
                                                                                            
}
                                                                                            
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::conservativeToPrimitive(double *U, double *V, double phi)
{
  if (phi>=0.0)
    conservativeToPrimitiveJWLEuler(omega, computeFrho(U,phi), U, V);
  else  conservativeToPrimitiveJWLEuler(omegap, computeFrho(U,phi), U, V);

}
//------------------------------------------------------------------------------
inline
int VarFcnJWLInJWLEuler3D::conservativeToPrimitiveVerification(int glob, double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveJWLEuler(omega, computeFrho(U,phi), U, V);
    return VerificationJWLEuler(glob, pmin, invomega, computeFrho(U,phi), U, V);
  }else{
    conservativeToPrimitiveJWLEuler(omegap, computeFrho(U,phi), U, V);
    return VerificationJWLEuler(glob, pminp, invomegap, computeFrho(U,phi), U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::primitiveToConservative(double *V, double *U, double phi)
{
  if (phi>=0.0)
    primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);
  else  primitiveToConservativeJWLEuler(invomegap, computeFrho(V,phi), V, U);
}
//------------------------------------------------------------------------------    
inline
bool VarFcnJWLInJWLEuler3D::updatePhaseChange(double *V, double *U, double phi,
                            double phin, double *Riemann, double weight)
{

  //nature of fluid at this node has not changed over time
  if(phi*phin > 0.0){
    if(phi>=0.0)
      primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);
    else primitiveToConservativeJWLEuler(invomegap, computeFrho(V,phi), V, U);

  //nature of fluid at this node has changed over time
  }else{
    if(weight<=0.0){
      fprintf(stdout, "negative weight in updatePhaseChange\n");
      exit(1);
    }
    for(int k=0; k<5; k++)
    //for(int k=0; k<1; k++)
      V[k] = Riemann[k]/weight;

    // from Fluid2 to Fluid1, ie JWL To Gas
    if(phi>=0.0 && phin<0.0)
      primitiveToConservativeJWLEuler(invomega, computeFrho(V,phi), V, U);

    // from Fluid1 to Fluid2, ie Gas To JWL
    else primitiveToConservativeJWLEuler(invomegap, computeFrho(V,phi), V, U);

  }

}
//------------------------------------------------------------------------------    
inline
void VarFcnJWLInJWLEuler3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
  else  multiplyBydVdUJWLEuler(omegap, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::multiplyBydVdU(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
  else  multiplyBydVdUJWLEuler(omegap, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::multiplyBydVdUT(double *V, double *vec, double *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
  else  multiplyBydVdUTJWLEuler(omegap, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, double phi)
{
  if (phi>=0.0)
    multiplyBydVdUTJWLEuler(omega, computeFrhop(V,phi), V, vec, res);
  else  multiplyBydVdUTJWLEuler(omegap, computeFrhop(V,phi), V, vec, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::preMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    preMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
  else  preMultiplyBydUdVJWLEuler(invomegap, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::postMultiplyBydVdU(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydVdUJWLEuler(omega, computeFrhop(V,phi), V, mat, res);
  else  postMultiplyBydVdUJWLEuler(omegap, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::postMultiplyBydUdV(double *V, double *mat, double *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
  else  postMultiplyBydUdVJWLEuler(invomegap, computeFrhop(V,phi), V, mat, res);
}
//------------------------------------------------------------------------------
inline
void VarFcnJWLInJWLEuler3D::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, double phi)
{
  if (phi>=0.0)
    postMultiplyBydUdVJWLEuler(invomega, computeFrhop(V,phi), V, mat, res);
  else  postMultiplyBydUdVJWLEuler(invomegap, computeFrhop(V,phi), V, mat, res);
}

//------------------------------------------------------------------------------

inline
void VarFcnJWLInJWLEuler3D::extrapolateBoundaryPrimitive(double un, double c, double *Vb,
                                                double *Vinter, double *V, double phi)
{
  extrapolatePrimitiveJWLEuler(un, c, Vb, Vinter, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnJWLInJWLEuler3D::extrapolateBoundaryCharacteristic(double n[3], double un,
                                double c, double *Vb, double *dV, double phi)
{
  extrapolateCharacteristicJWLEuler(un,c,n,Vb,dV);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#endif
