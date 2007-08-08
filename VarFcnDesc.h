#ifndef _VAR_FCN_DESC_H_
#define _VAR_FCN_DESC_H_

#include <VarFcn.h>

//------------------------------------------------------------------------------
class VarFcnPerfectGasEuler3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasEuler3D(IoData &);
  ~VarFcnPerfectGasEuler3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);

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

inline
void VarFcnPerfectGasEuler3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
  VerificationGasEuler(pmin, gam, Pstiff, U, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
   primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
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
//------------------------------------------------------------------------------

class VarFcnPerfectGasSA3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasSA3D(IoData &);
  ~VarFcnPerfectGasSA3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  double getTurbulentNuTilde(double *V) { return V[5]; }
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

inline
void VarFcnPerfectGasSA3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  conservativeToPrimitiveGasSA(gam1, U, V);
  VerificationGasSA(pmin, gam1, U, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
  primitiveToConservativeGasSA(U, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasSA3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUGasSA(V, vec, res);
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
//------------------------------------------------------------------------------

class VarFcnPerfectGasKE3D : public VarFcnPerfectGas {

public:
  VarFcnPerfectGasKE3D(IoData &);
  ~VarFcnPerfectGasKE3D() {}

  void conservativeToPrimitive(double *, double *, double = 0.0);
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
  void multiplyBydVdU(double *, double *, double *, double = 0.0);
  void preMultiplyBydUdV(double *, double *, double *, double = 0.0);
  void postMultiplyBydVdU(double *, double *, double *, double = 0.0);
  void postMultiplyBydUdV(double *, double *, double *, double = 0.0);
  double getTurbulentKineticEnergy(double *V) { return V[5]; }
  double getTurbulentDissipationRate(double *V) { return V[6]; }
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

inline
void VarFcnPerfectGasKE3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
   conservativeToPrimitiveGasKE(gam1, U, V);
   VerificationGasKE(pmin, gam1, U, V);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
  primitiveToConservativeGasKE(V, U);
}

//------------------------------------------------------------------------------

inline
void VarFcnPerfectGasKE3D::multiplyBydVdU(double *V, double *vec, double *res, double phi)
{
  multiplyBydVdUGasKE(V, vec, res);
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
//-------------------------------------------------------------------------------
class VarFcnWaterCompressibleEuler3D : public VarFcnWaterCompressible {
public:
  VarFcnWaterCompressibleEuler3D(IoData &);
  ~VarFcnWaterCompressibleEuler3D() {}
  void conservativeToPrimitive(double *, double *, double = 0.0);
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
                                                                                            
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
void VarFcnWaterCompressibleEuler3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  conservativeToPrimitiveLiquidEuler(invCv, U, V);
  VerificationLiquidEuler(pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
}

//------------------------------------------------------------------------------
                                                                                            
inline
void VarFcnWaterCompressibleEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
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

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------                                                                                           

class VarFcnGasInGasEuler3D : public VarFcnGasInGas {

public:

  VarFcnGasInGasEuler3D(IoData &);
  ~VarFcnGasInGasEuler3D() {}
  
  void conservativeToPrimitive(double *, double *, double = 0.0);
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
                                                                                            
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
void VarFcnGasInGasEuler3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
    VerificationGasEuler(pmin, gam, Pstiff, U, V);
  }else{
    conservativeToPrimitiveGasEuler(gamp, Pstiffp, U, V);
    VerificationGasEuler(pminp, gamp, Pstiffp, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnGasInGasEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
  if (phi>=0.0)
    primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
  else  primitiveToConservativeGasEuler(gamp, invgamp1, Pstiffp, V, U);
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
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
                                                                                            
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
void VarFcnLiquidInLiquidEuler3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
    VerificationLiquidEuler(pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
  }else{
    conservativeToPrimitiveLiquidEuler(invCvbis, U, V);
    VerificationLiquidEuler(pminp, alpha_waterbis, beta_waterbis, Pref_waterbis, invCvbis, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnLiquidInLiquidEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
  if (phi>=0.0)
    primitiveToConservativeLiquidEuler(Cv, V, U);
  else  primitiveToConservativeLiquidEuler(Cvbis, V, U);
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
  void conservativeToPrimitiveVerification(double *, double *, double = 0.0);
  void primitiveToConservative(double *, double *, double = 0.0, double * = 0, double * = 0);
                                                                                            
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
void VarFcnGasInLiquidEuler3D::conservativeToPrimitiveVerification(double *U, double *V, double phi)
{
  if (phi>=0.0){
    conservativeToPrimitiveLiquidEuler(invCv, U, V);
    VerificationLiquidEuler(pmin, alpha_water, beta_water, Pref_water, invCv, U, V);
  }else{
    conservativeToPrimitiveGasEuler(gam, Pstiff, U, V);
    VerificationGasEuler(pminp, gam, Pstiff, U, V);
  }

}
//------------------------------------------------------------------------------
inline
void VarFcnGasInLiquidEuler3D::primitiveToConservative(double *V, double *U, double phi, double *phi1, double *vgf)
{
  if (!phi1) {
    if (phi>=0.0)
      primitiveToConservativeLiquidEuler(Cv, V, U);
    else  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
  }
  else {
    if (phi * (*phi1) < 0.0) { // the fluid at this node has changed
      if (phi >= 0.0 && (*phi1) < 0.0) { // gas to liquid
        V[0]   = vgf[0];
        V[1]   = vgf[1];
        V[2]   = vgf[2];
        V[3]   = vgf[3];
        V[4]   = vgf[4];
        primitiveToConservativeLiquidEuler(Cv, V, U);
      }
      if (phi < 0.0 && (*phi1) >= 0.0) { // liquid to gas
        double P_g  = getPressure(V, *phi1); 
        V[4] = P_g;
        primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
      }
    }
    else {
      if (phi>=0.0)
        primitiveToConservativeLiquidEuler(Cv, V, U);
      else  primitiveToConservativeGasEuler(gam, invgam1, Pstiff, V, U);
    }
  }
}
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
//------------------------------------------------------------------------------                                                                                     
#endif

