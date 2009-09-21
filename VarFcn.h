#ifndef _VAR_FCN_H_
#define _VAR_FCN_H_

#include <DistVector.h>
#include <Vector3D.h>

#include <math.h>
#include <complex>
typedef std::complex<double> bcomp;
#include <iostream>

using std::cout;
using std::endl;

#include <IoData.h>

template<int dim> class DistExactRiemannSolver;

//--------------------------------------------------------------------------
class VarFcn {

public:
  
  const char** pname;
  enum Type{ GAS = 0, LIQUID = 1, JWL = 2, GASINGAS = 3, GASINLIQUID = 4, LIQUIDINLIQUID = 5, JWLINGAS = 6, JWLINJWL = 7} type;
  enum SubType { NONE = 0, IDEAL = 1, STIFFENED = 2} subType;
  Vec3D meshVel;

  double pmin;
  double pminp;
  bool verif_clipping;
  bool node_change;

  double gravity[3];
  double gravity_norm;
  double depth;

  VarFcn(IoData &iod); 
  VarFcn() { meshVel = 0.0; fprintf(stderr, "MeshVel initialized: %e %e %e\n", meshVel[0], meshVel[1], meshVel[2]);}
  ~VarFcn() {}

  double *getMeshVel()  { return meshVel.v; }

  template<int dim>
  void conservativeToPrimitive(SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0);
  template<int dim>
  void conservativeToPrimitive(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0);
  template<int dim>
  void primitiveToConservative(SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0);
  template<int dim>
  void primitiveToConservative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0);

  template<int dim>
  void updatePhaseChange(DistSVec<double,dim> &V, DistSVec<double,dim> &U,
          DistVec<double> &Phi, DistVec<double> &Phin,
          DistSVec<double,dim> *Vgf, DistVec<double> *Vgfweight,
          DistExactRiemannSolver<dim> *riemann);
  template<int dim>
  void updatePhaseChange(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &, 
			DistVec<double> &, DistSVec<double,dim> *, DistVec<double> * = 0);

  void setMeshVel(Vec3D &v)  { meshVel = v; }

protected:
  //GAS_EULER
  void conservativeToPrimitiveGasEuler(double, double, double *, double *);
  void primitiveToConservativeGasEuler(double, double, double, double *, double *);
  void primitiveToCharacteristicVariationsGasEuler(double *, double, double, double *,
                                                   double *, double *);
  void characteristicToPrimitiveVariationsGasEuler(double *, double, double, double *,
                                                   double *, double *);
  void computedVdUGasEuler(double, double *, double *);
  void computedUdVGasEuler(double, double *, double *);
  void multiplyBydVdUGasEuler(double, double *, double *, double *);
  void multiplyBydVdUGasEuler(double, double *, bcomp *, bcomp *);
  void multiplyBydVdUTGasEuler(double, double *, double *, double *) ;
  void multiplyBydVdUTGasEuler(double, double *, bcomp *, bcomp *);
  void postMultiplyBydUdVGasEuler(double, double *, double *, double *);
  void postMultiplyBydVdUGasEuler(double, double *, double *, double *);
  void postMultiplyBydUdVGasEuler(double, double *, bcomp *, bcomp *);
  void preMultiplyBydUdVGasEuler(double, double *, double *, double *);
  void extrapolatePrimitiveGasEuler(double, double, double*, double*, double*);
  void extrapolatePrimitiveGasEulerPUT(double, double, double*, double*, double*);
  void extrapolateCharacteristicGasEuler(double ,double , double*,
                                         double, double, double*, double*);
  int VerificationGasEuler(int, double, double, double, double*, double*);

// Included (MB)
  void conservativeToPrimitiveDerivativeGasEuler(double, double, double, double *, double *, double *, double *);
  void primitiveToConservativeDerivativeGasEuler(double, double, double, double, double *, double *, double *, double *);

  //GAS_SA
  void conservativeToPrimitiveGasSA(double, double *, double *);
// Modified (MB)
  void primitiveToConservativeGasSA(double, double *, double *);

  void computedUdVGasSA(double, double *, double *);
  void computedVdUGasSA(double, double *, double *);

// Modified (MB)
  void multiplyBydVdUGasSA(double, double *, double *, double *);

  void postMultiplyBydUdVGasSA(double, double *, double *, double *);
  void postMultiplyBydVdUGasSA(double, double *, double *, double *);
  void preMultiplyBydUdVGasSA(double, double *, double *, double *);
  int VerificationGasSA(int, double, double , double *, double *);

// Included (MB)
  void conservativeToPrimitiveDerivativeGasSA(double, double *, double *, double *, double *);
  void primitiveToConservativeDerivativeGasSA(double, double *, double *, double *, double *);

  //GAS_KE
  void conservativeToPrimitiveGasKE(double, double *, double *);
// Modified (MB)
  void primitiveToConservativeGasKE(double, double *, double *);

  void computedUdVGasKE(double, double *, double *);
  void computedVdUGasKE(double, double *, double *);

// Modified (MB)
  void multiplyBydVdUGasKE(double, double *, double *, double *);

  void postMultiplyBydUdVGasKE(double, double *, double *, double *);
  void postMultiplyBydVdUGasKE(double, double *, double *, double *);
  void preMultiplyBydUdVGasKE(double, double *, double *, double *);
  int VerificationGasKE(int, double, double , double *, double *);

// Included (MB)
  void conservativeToPrimitiveDerivativeGasKE(double, double *, double *, double *, double *);
  void primitiveToConservativeDerivativeGasKE(double, double *, double *, double *, double *);

  //LIQUID_EULER
  void conservativeToPrimitiveLiquidEuler(double, double *, double *);
  void primitiveToConservativeLiquidEuler(double, double *, double *);
  void primitiveToCharacteristicVariationsLiquidEuler(double *, double, double, double, double, 
                                                      double *, double *, double *);
  void characteristicToPrimitiveVariationsLiquidEuler(double *, double, double, double, double,
                                                      double *, double *, double *);
  void computedVdULiquidEuler(double, double *, double *);
  void computedUdVLiquidEuler(double, double *, double *);
  void multiplyBydVdULiquidEuler(double, double *, double *, double *);
  void multiplyBydVdULiquidEuler(double, double *, bcomp *, bcomp *);
  void multiplyBydVdUTLiquidEuler(double, double *, double *, double *) ;
  void multiplyBydVdUTLiquidEuler(double, double *, bcomp *, bcomp *);
  void postMultiplyBydUdVLiquidEuler(double, double *, double *, double *);
  void postMultiplyBydVdULiquidEuler(double, double *, double *, double *);
  void postMultiplyBydUdVLiquidEuler(double, double *, bcomp *, bcomp *);
  void preMultiplyBydUdVLiquidEuler(double, double *, double *, double *);
  void extrapolatePrimitiveLiquidEuler(double, double, double*, double*, double*);
  void extrapolateCharacteristicLiquidEuler(double, double, double, double, 
					double*, double, double, double*, double*);

  void computeNewPrimitiveLiquidEuler(double, double, double, double *, double *);
  void computeOldPrimitiveLiquidEuler(double, double, double, double *, double *);
  int VerificationLiquidEuler(int, double, double, double, double, double, double *, double *);

  //JWL EOS
  void conservativeToPrimitiveJWLEuler(double, double, double *, double *);
  void primitiveToConservativeJWLEuler(double, double, double *, double *);
  void primitiveToCharacteristicVariationsJWLEuler(double, double *, double *,
                                                   double *, double *);
  void characteristicToPrimitiveVariationsJWLEuler(double, double *, double *,
                                                   double *, double *);
  void computedVdUJWLEuler(double, double, double *, double *);
  void computedUdVJWLEuler(double, double, double *, double *);
  void multiplyBydVdUJWLEuler(double, double, double *, double *, double *);
  void multiplyBydVdUJWLEuler(double, double, double *, bcomp *, bcomp *);
  void multiplyBydVdUTJWLEuler(double, double, double *, double *, double *) ;
  void multiplyBydVdUTJWLEuler(double, double, double *, bcomp *, bcomp *);
  void preMultiplyBydUdVJWLEuler(double, double, double *, double *, double *);
  void postMultiplyBydUdVJWLEuler(double, double, double *, double *, double *);
  void postMultiplyBydVdUJWLEuler(double, double, double *, double *, double *);
  void postMultiplyBydUdVJWLEuler(double, double, double *, bcomp *, bcomp *);
  void extrapolatePrimitiveJWLEuler(double, double, double*, double*, double*);
  void extrapolateCharacteristicJWLEuler(double ,double , double*,
                                         double*, double*);
  int VerificationJWLEuler(int, double, double, double, double*, double*);

// Included (MB)
  void conservativeToPrimitiveDerivativeJWLEuler(double, double, double, double *, double *, double *, double *);
  void primitiveToConservativeDerivativeJWLEuler(double, double, double, double, double *, double *, double *, double *);
public:

  template<int dim>
  void computeNewPrimitive(SVec<double,dim> &, SVec<double,dim> &);
  virtual void computeNewPrimitive(double *, double *){}
  virtual void computeOldPrimitive(double *, double *){}

  virtual void conservativeToPrimitive(double *, double *, double = 0.0) = 0; 
  virtual int  conservativeToPrimitiveVerification(int, double *, double *, double = 0.0) = 0; 
  virtual void primitiveToConservative(double *, double *, double = 0.0) = 0;
  virtual bool updatePhaseChange(double *, double *, double, double, double *, double = 1.0){};
  bool doVerification(){ return !(pmin<0 && pminp<0); }
  
  virtual void multiplyBydVdU(double *, double *, double *, double = 0.0) {
    cout << "ERROR: multiplyBydVdU Function in VarFcn" << endl; }
  virtual void multiplyBydVdU(double *, bcomp *, bcomp *, double = 0.0)  {
    cout << "ERROR: multiplyBydVdU Function in VarFcn for bcomp" << endl; }
  virtual void multiplyBydVdUT(double *, double *, double *, double = 0.0) {
    cout << "ERROR: multiplyBydVdUT in VarFcn for type double" << endl;}
  virtual void multiplyBydVdUT(double *, bcomp *, bcomp *, double = 0.0) {
    cout << "ERROR: multiplyBydVdUT in VarFcn for type bcomp" << endl;}
  virtual void postMultiplyBydUdV(double *, double *, double *, double = 0.0)  {
    cout << "ERROR: postMultiplyBydUdV Function in VarFcn" << endl; }
  virtual void postMultiplyBydVdU(double *, double *, double *, double = 0.0)  {
    cout << "ERROR: postMultiplyBydVdU Function in VarFcn" << endl; }
  virtual void postMultiplyBydUdV(double *, bcomp *, bcomp *, double = 0.0)  {
    cout << "ERROR: postMultiplyBydVdU Function in VarFcn for becomp" << endl; }
  virtual void preMultiplyBydUdV(double *, double *, double *, double = 0.0) = 0;
  virtual void extrapolateBoundaryPrimitive(double, double, double *, double *,
                                   double *, double = 0.0) {}
  virtual void extrapolateBoundaryCharacteristic(double*, double, double, double*,
                                   double*, double = 0.0) {}

  virtual int getType() { return type; }
  virtual int getSubType() { return subType; }
  
  virtual double getGamma() {
        fprintf(stderr, "*** Warning:  getGamma Function not defined\n");
        return 0.0; }
  virtual double getGammabis() {
        fprintf(stderr, "*** Warning:  getGammabis Function not defined\n");
        return 0.0; }
  virtual double getGamma1() {
        fprintf(stderr, "*** Warning:  getGamma1 Function not defined\n");
        return 0.0; }
  virtual double getGammabis1() {
        fprintf(stderr, "*** Warning:  getGammabis1 Function not defined\n");
        return 0.0; }
  virtual double getPressureConstant() {
        fprintf(stderr, "*** Warning:  getPressureConstant Function not defined\n");
        return 0.0; }
  virtual double getPressureConstantbis() {
        fprintf(stderr, "*** Warning:  getPressureConstantbis Function not defined\n");
        return 0.0; }
  virtual double getCv() {
        fprintf(stderr, "*** Warning:  getCv Function not defined\n");
        return 0.0; }
  virtual double getCvbis() {
        fprintf(stderr, "*** Warning:  getCvbis Function not defined\n");
        return 0.0; }
  virtual double getAlphaWater() {
        fprintf(stderr, "*** Warning: getAlphaWater Function not defined\n");
        return 0.0;}
  virtual double getAlphaWaterbis() {
        fprintf(stderr, "*** Warning: getAlphaWaterbis Function not defined\n");
        return 0.0;}
  virtual double getBetaWater() {
        fprintf(stderr, "*** Warning: getBetaWater Function not defined\n");
        return 0.0;}
  virtual double getBetaWaterbis() {
        fprintf(stderr, "*** Warning: getBetaWaterbis Function not defined\n");
        return 0.0;}
  virtual double getPrefWater() {
        fprintf(stderr, "*** Warning: getPrefWater Function not defined\n");
        return 0.0;}
  virtual double getPrefWaterbis() {
        fprintf(stderr, "*** Warning: getPrefWaterbis Function not defined\n");
        return 0.0;}
  virtual double getOmega() {
        fprintf(stderr, "*** Warning: getOmega Function not defined\n");
        return 0.0;}
  virtual double getOmegap1() {
        fprintf(stderr, "*** Warning: getOmegap1 Function not defined\n");
        return 0.0;}
  virtual double getInvOmega() {
        fprintf(stderr, "*** Warning: getInvOmega Function not defined\n");
        return 0.0;}
  virtual double getInvOmegap1() {
        fprintf(stderr, "*** Warning: getInvOmegap1 Function not defined\n");
        return 0.0;}
  virtual double getA1() {
        fprintf(stderr, "*** Warning: getA1 Function not defined\n");
        return 0.0;}
  virtual double getA2() {
        fprintf(stderr, "*** Warning: getA2 Function not defined\n");
        return 0.0;}
  virtual double getR1() {
        fprintf(stderr, "*** Warning: getR1 Function not defined\n");
        return 0.0;}
  virtual double getR2() {
        fprintf(stderr, "*** Warning: getR2 Function not defined\n");
        return 0.0;}
  virtual double getRhoref() {
        fprintf(stderr, "*** Warning: getRhoref Function not defined\n");
        return 0.0;}
  virtual double getR1r() {
        fprintf(stderr, "*** Warning: getR1r Function not defined\n");
        return 0.0;}
  virtual double getR2r() {
        fprintf(stderr, "*** Warning: getR2r Function not defined\n");
        return 0.0;}
  virtual double getOmegabis() {
        fprintf(stderr, "*** Warning: getOmegabis Function not defined\n");
        return 0.0;}
  virtual double getOmegap1bis() {
        fprintf(stderr, "*** Warning: getOmegap1bis Function not defined\n");
        return 0.0;}
  virtual double getInvOmegabis() {
        fprintf(stderr, "*** Warning: getInvOmegabis Function not defined\n");
        return 0.0;}
  virtual double getInvOmegap1bis() {
        fprintf(stderr, "*** Warning: getInvOmegap1bis Function not defined\n");
        return 0.0;}
  virtual double getA1bis() {
        fprintf(stderr, "*** Warning: getA1bis Function not defined\n");
        return 0.0;}
  virtual double getA2bis() {
        fprintf(stderr, "*** Warning: getA2bis Function not defined\n");
        return 0.0;}
  virtual double getR1bis() {
        fprintf(stderr, "*** Warning: getR1bis Function not defined\n");
        return 0.0;}
  virtual double getR2bis() {
        fprintf(stderr, "*** Warning: getR2bis Function not defined\n");
        return 0.0;}
  virtual double getRhorefbis() {
        fprintf(stderr, "*** Warning: getRhorefbis Function not defined\n");
        return 0.0;}
  virtual double getR1rbis() {
        fprintf(stderr, "*** Warning: getR1rbis Function not defined\n");
        return 0.0;}
  virtual double getR2rbis() {
        fprintf(stderr, "*** Warning: getR2rbis Function not defined\n");
        return 0.0;}
        
  virtual Vec3D getVelocity(double *V) { return Vec3D(V[1], V[2], V[3]); }
  virtual double getDensity(double *V) { return V[0]; }
  virtual double getVelocityX(double *V) {return V[1];}
  virtual double getVelocityY(double *V) {return V[2];}
  virtual double getVelocityZ(double *V) {return V[3];}
  virtual double computeFrho(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeFrhop(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeFrho(double rho, double phi = 0.0) {return 0.0;}
  virtual double computeFrhop(double rho, double phi = 0.0) {return 0.0;}
  virtual double getPressure(double *V, double phi = 0.0) {return 0.0;}
  virtual double checkPressure(double *V, double phi = 0.0) {return 0.0;}
  virtual void   setDensity(double *V, double *Vorig){}
  virtual void   setPressure(double *V, double *Vorig, double phi=0.0){}
  virtual double computeTemperature(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeRhoEnergy(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeRhoEpsilon(double *V, double phi = 0.0) {return 0.0;} //this function computes the internal energy (=rho*e-0.5*rho*u^2)
  virtual double computeU2(double *V) { return V[1]*V[1]+V[2]*V[2]+V[3]*V[3]; }
  virtual double computeWtU2(double *V) { 
          double v1 = V[1] - meshVel[0]; double v2 = V[2] - meshVel[1]; double v3 = V[3] - meshVel[2];
          return v1*v1+v2*v2+v3*v3;
  }

  virtual double computeMachNumber(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeWtMachNumber(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeSoundSpeed(double *V, double phi = 0.0) {return 0.0;}
  virtual double computeSoundSpeed(double density, double entropy, double phi = 0.0) {return 0.0;}
  virtual double computeEntropy(double density, double pressure, double phi = 0.0){ return 0.0; }
  virtual double computeIsentropicPressure(double entropy, double density, double phi = 0.0){ return 0.0; }
  virtual double computeTotalPressure(double machr, double* V, double phi = 0.0) {return 0.0;}
                                                                     
  virtual double getTurbulentNuTilde(double *V) { return 0.0; }
  virtual double getTurbulentKineticEnergy(double *V) { return 0.0; }
  virtual double getTurbulentDissipationRate(double *V) { return 0.0; }

// Included (MB)
  template<int dim>
  void conservativeToPrimitiveDerivative(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0);
  template<int dim>
  void conservativeToPrimitiveDerivative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0);
  template<int dim>
  void primitiveToConservativeDerivative(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0);
  template<int dim>
  void primitiveToConservativeDerivative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0);

  virtual void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) = 0;
  virtual void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0) = 0;

  virtual Vec3D getDerivativeOfVelocity(double *dV) { return Vec3D(dV[1], dV[2], dV[3]); }
  virtual double computeDerivativeOfTemperature(double *V, double *dV, double phi = 0.0) { return 0.0; }
  virtual double computeDerivativeOfMachNumber(double *V, double *dV, double dMach, double phi = 0.0) { return 0.0; }
  virtual double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach, double phi = 0.0) { return 0.0; }
  virtual double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dmach, double phi = 0.0) { return 0.0; }

  virtual double getVelocityNorm(double *V) { return 0.0; }
  virtual double getDerivativeOfVelocityNorm(double *V, double *dV) { return 0.0; }
  virtual void setDensity(double rho, double *V) { V[0] = rho; }
  virtual void setPressure(double p, double *V) { V[4] = p; }
  virtual double computeDensity(double p, double temp) { return 0.0; }
  virtual double computePressure(double rho, double temp) { return 0.0; }
  virtual double getTotalEnergy(double *U) { return U[4] / U[0]; }
  virtual double getDerivativeOfPressureConstant() { return 0.0; }
  virtual void rstVar(IoData &iod) {}
  virtual void rV(IoData &iod) {  pmin  = iod.eqs.fluidModel.pmin;
                                  pminp = iod.eqs.fluidModel2.pmin; }

  virtual double gravity_value() const { return gravity_norm; }
  virtual double hydrostaticPressure(double rho, double *X) { 
    return rho*(gravity_norm*depth + 
                gravity[0]*X[0]+gravity[1]*X[1]+gravity[2]*X[2]); }
  virtual double hydrodynamicPressure(double *V, double *X, double phi = 0.0)
  { //check that getPressure does not call the virtual fcn of the present
    //class which always returns 0.0!!!!
    double localP = getPressure(V,phi);
    double localStaticP = hydrostaticPressure(V[0],X);
    return localP - localStaticP; }
  virtual double DerivativeHydrostaticPressure(double rho, double drho,
                                               double *X, double *dX)
  {
    return drho*(gravity_norm*depth +
                 gravity[0]*X[0]+gravity[1]*X[1]+gravity[2]*X[2])
         + rho*(gravity[0]*dX[0]+gravity[1]*dX[1]+gravity[2]*dX[2]);
  }
  virtual double computePressureCoefficient(double *V, double pinfty, 
                                            double mach, bool dimFlag, double phi = 0.0) { 
     fprintf(stderr, "... ERROR: computePressureCoefficient is only implemented for Perfect Gas Simulations. Exiting.\n"); 
     exit(-1);
  }
                                                                    
};
//------------------------------------------------------------------------------

class VarFcnPerfectGas : public VarFcn {

protected:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;  

// Included (MB)
  double dPstiff;  

public:
  VarFcnPerfectGas(IoData &);
  ~VarFcnPerfectGas() {}
                            
  double getGamma()  {return gam;}
  double getGamma1()  {return gam1;}
  double getPressureConstant() {return Pstiff;}

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) { return V[4]+Pstiff; }
  void   setDensity(double *V, double *Vorig)  { V[0] = Vorig[0]; }
  void   setPressure(double *V, double *Vorig, double phi=0.0) { V[4] = Vorig[4]; }
  double computeTemperature(double *V, double phi = 0.0 ) {
    if (isnan(1.0/V[0])) {
      fprintf(stderr, "ERROR*** computeTemp\n");
      exit(1);
    }
    return invgam1 * (V[4]+gam*Pstiff) / V[0];
  }
  double computeRhoEnergy(double *V, double phi = 0.0) {
    return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V, double phi = 0.0) { return invgam1 * (V[4]+gam*Pstiff); }

  double computeMachNumber(double *V, double phi = 0.0) {
    return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)));
  }

  double computeWtMachNumber(double *V, double phi = 0.0) {
         double v1 = V[1] - meshVel[0];
         double v2 = V[2] - meshVel[1];
         double v3 = V[3] - meshVel[2];
         double m = sqrt((v1*v1+v2*v2+v3*v3) * V[0] / (gam * (V[4]+Pstiff)));
         if (isnan(m))  {
           fprintf(stderr, "Nan: %e %e %e\n", meshVel[0], meshVel[1], meshVel[2]);
           exit(-1);
         }
         return m;
         //return sqrt((v1*v1+v2*v2+v3*v3) * V[0] / (gam * (V[4]+Pstiff)));
  }
  double computeSoundSpeed(double *V, double phi = 0.0) { return sqrt(gam * (V[4]+Pstiff) / V[0]); }
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) 
    { return sqrt(gam * entropy*pow(density,gam-1.0)); }
  double computeEntropy(double density, double pressure, double phi = 0.0)
    { return (pressure+Pstiff)/pow(density,gam); }
  double computeIsentropicPressure(double entropy, double density, double phi = 0.0)
    { return entropy*pow(density,gam)-Pstiff; }
  double computeTotalPressure(double machr, double* V, double phi = 0.0) {
    double mach = computeMachNumber(V);
    double machr2 = machr*machr;
    double popr = V[4]*gam*machr2;
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    //double opmachr = 1.0 + 0.5*gam1*machr2;
    return V[4]*pow(opmach, gam*invgam1);
  }

  double computePressureCoefficient(double *V, double pinfty, double mach, bool dimFlag, double phi = 0.0) { 
    if (dimFlag)
      return 2.0 * (V[4] - pinfty);
    else 
      return 2.0 * (V[4] - 1.0/(gam*mach*mach));
  }

// Included (MB)
  double getVelocityNorm(double *V) { return sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]); }
  double getDerivativeOfVelocityNorm(double *V, double *dV) { return (V[1]*dV[1]+V[2]*dV[2]+V[3]*dV[3])/sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]); }
  double computeDerivativeOfTemperature(double *V, double *dV, double phi = 0.0) {return invgam1 * ( dV[4]*V[0] - V[4]*dV[0] ) / ( V[0]*V[0] ); }
  double computeDerivativeOfMachNumber(double *V, double *dV, double dMach, double phi = 0.0) {
    return 1/(2.0*sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)))) * ( ( (2.0*(V[1]*dV[1] + V[2]*dV[2] + V[3]*dV[3]) * V[0] + (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * dV[0]) * (V[4]+Pstiff) - (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] * (dV[4] + dPstiff*dMach) ) / ( (V[4]+Pstiff) * (V[4]+Pstiff) ) );
  }
  double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach, double phi = 0.0) { return 1.0/( 2.0*sqrt(gam * (V[4]+Pstiff) / V[0]) ) * gam * ( (dV[4]+dPstiff*dMach) * V[0] - (V[4]+Pstiff) * dV[0] ) / ( V[0] * V[0] ); }
  double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dMach, double phi = 0.0) {
    double mach = computeMachNumber(V);
    double dmach = computeDerivativeOfMachNumber(V, dV, dMach);
    double machr2 = machr*machr;
    double dmachr2 = 2.0*machr*dmachr;
    double popr = V[4]*gam*machr2;
    double dpopr = dV[4]*gam*machr2 + V[4]*gam*dmachr2;
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    double dopmach = gam1*mach*dmach;

//    double opmachr = 1.0 + 0.5*gam1*machr2;
//    double dopmachr = 0.5*gam1*dmachr2;
//    return dpopr*pow(opmach/opmachr, gam*invgam1) + popr*gam*invgam1*pow(opmach/opmachr, (gam*invgam1-1))*((dopmach*opmachr - opmach*dopmachr)/(opmachr*opmachr));

    return dV[4]*pow(opmach, gam*invgam1) + V[4]*gam*invgam1*pow(opmach, (gam*invgam1-1))*dopmach;
  }
  double computeDensity(double p, double temp) { return invgam1 * p / temp; }
  double computePressure(double rho, double temp) { return temp * rho * gam1; }
  double getDerivativeOfPressureConstant() {return dPstiff;}
  void rstVar(IoData &iod) { dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach)); rV(iod);}

};

//-----------------------------------------------------------------------------------

inline
VarFcnPerfectGas::VarFcnPerfectGas(IoData &iod) : VarFcn(iod) {

  type = GAS;
  if (iod.eqs.fluidModel.gasModel.type == GasModelData::IDEAL)
    subType == IDEAL;
  else if(iod.eqs.fluidModel.gasModel.type == GasModelData::STIFFENED)
    subType == STIFFENED;
  else{
    subType == NONE;
    fprintf(stderr, "no subType defined for the VarFcnPerfectGas\n");
  }
  gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;

// Included (MB)
  dPstiff = 0.0;//iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));

}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

class VarFcnWaterCompressible : public VarFcn {

protected:
  double Cv;
  double invCv;
  double alpha_water;
  double beta_water;
  double Pref_water;

public:
  VarFcnWaterCompressible(IoData &);
  ~VarFcnWaterCompressible() {}
    
  double getCv()  {return Cv;}
  double getAlphaWater()  {return alpha_water;}
  double getBetaWater()  {return beta_water;}
  double getPrefWater()  {return Pref_water;}

  double getPressure(double *V, double phi = 0.0) { return Pref_water + alpha_water * pow(V[0], beta_water); }
  double checkPressure(double *V, double phi = 0.0) { return Pref_water + alpha_water * pow(V[0], beta_water); }
  void   setDensity(double *V, double *Vorig)  { V[0] = Vorig[0];}
  void   setPressure(double *V, double *Vorig, double phi=0.0) { V[0] = Vorig[0];}
  double computeTemperature(double *V, double phi = 0.0) { return V[4]; }
  double computeRhoEnergy(double *V, double phi = 0.0) {
    return V[0] * Cv * V[4] + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V, double phi = 0.0) { return V[0] * Cv * V[4]; }
  double computeSoundSpeed(double *V, double phi = 0.0) { return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0)); }
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) 
    { return sqrt(alpha_water * beta_water * pow(density, beta_water - 1.0)); }
  double computeMachNumber(double *V, double phi = 0.0) {
        double c=sqrt(alpha_water * beta_water *pow(V[0], beta_water - 1.0));
        double u=sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
        return u/c;
  }
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
  return 0.0;
  }

// Included (MB)
  void rstVar(IoData &iod) { rV(iod);}

};

//--------------------------------------------------------------------------------------

inline
VarFcnWaterCompressible::VarFcnWaterCompressible(IoData &iod)
: VarFcn(iod) {

  type = LIQUID;
  subType = NONE;

  Cv=1.0;
  invCv=1.0/Cv;
  alpha_water = iod.eqs.fluidModel.liquidModel.alpha;
  beta_water  = iod.eqs.fluidModel.liquidModel.beta;
  Pref_water  = iod.eqs.fluidModel.liquidModel.Pref;

}

//-----------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

class VarFcnJWL: public VarFcn {

protected:
  double omega;
  double invomega;
  double omegap1;
  double invomegap1;
  double A1, A2, R1, R2, rhoref;
  double R1r, R2r;

public:
  VarFcnJWL(IoData &);
  ~VarFcnJWL() {}
    
  double getOmega()       {return omega;}
  double getOmegap1()     {return omegap1;}
  double getInvOmega()    {return invomega;}
  double getInvOmegap1()  {return invomegap1;}
  double getA1()          {return A1;}
  double getA2()          {return A2;}
  double getR1()          {return R1;}
  double getR2()          {return R2;}
  double getRhoref()      {return rhoref;}
  double getR1r()         {return R1r;}
  double getR2r()         {return R2r;}
  
  double computeFrho(double *V, double phi = 0.0) {
    return A1*(1.0-omega*V[0]/R1r)*exp(-R1r/V[0])
         + A2*(1.0-omega*V[0]/R2r)*exp(-R2r/V[0]);
  }
  double computeFrhop(double *V, double phi = 0.0) {
    double rho2 = V[0]*V[0];
    return A1*(-omega/R1r+(1.0-omega*V[0]/R1r)*R1r/rho2)*exp(-R1r/V[0])
         + A2*(-omega/R2r+(1.0-omega*V[0]/R2r)*R2r/rho2)*exp(-R2r/V[0]);
  }
  double computeFrho(double rho, double phi = 0.0) {
    return A1*(1.0-omega*rho/R1r)*exp(-R1r/rho)
         + A2*(1.0-omega*rho/R2r)*exp(-R2r/rho);
  }
  double computeFrhop(double rho, double phi = 0.0) {
    double rho2 = rho*rho;
    return A1*(-omega/R1r+(1.0-omega*rho/R1r)*R1r/rho2)*exp(-R1r/rho)
         + A2*(-omega/R2r+(1.0-omega*rho/R2r)*R2r/rho2)*exp(-R2r/rho);
  }
  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) { return V[4] - (computeFrho(V) - computeFrhop(V)*V[0])*invomegap1;}
  void   setDensity(double *V, double *Vorig)  { V[0] = Vorig[0];}
  void   setPressure(double *V, double *Vorig, double phi=0.0) { V[4] = Vorig[4];}
  double computeTemperature(double *V, double phi = 0.0) { return invomega * (V[4]-computeFrho(V)) / V[0]; }
  double computeRhoEnergy(double *V, double phi = 0.0) {
    return invomega * (V[4]-computeFrho(V)) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V, double phi = 0.0) { return invomega * (V[4]-computeFrho(V)); }
  double computeSoundSpeed(double *V, double phi = 0.0) { return sqrt((omegap1*V[4] - computeFrho(V) + V[0]*computeFrhop(V))/V[0]); }
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) 
    { double pressure = entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
      return sqrt((omegap1*pressure - computeFrho(density) + density*computeFrhop(density))/density); }
  double computeEntropy(double density, double pressure, double phi = 0.0)
    { return (pressure - A1*exp(-R1r/density) - A2*exp(-R2r/density))/pow(density,omegap1); }
  double computeIsentropicPressure(double entropy, double density, double phi = 0.0)
    { return entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density); }
  double computeMachNumber(double *V, double phi = 0.0) {
    double c=computeSoundSpeed(V);
    double u=sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]));
    return u/c;
  }
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    fprintf(stderr, "*** Error: computeTotalPressure not yet implemented for VarFcnJWL\n");
    return 0.0;
  }

// Included (MB)
  void rstVar(IoData &iod) { rV(iod);}

};

//--------------------------------------------------------------------------------------

inline
VarFcnJWL::VarFcnJWL(IoData &iod)
: VarFcn(iod) {

  type = JWL;
  subType == NONE;

  omega        = iod.eqs.fluidModel.jwlModel.omega;
  omegap1      = iod.eqs.fluidModel.jwlModel.omega + 1.0;
  invomega     = 1.0/iod.eqs.fluidModel.jwlModel.omega;
  invomegap1   = 1.0/(iod.eqs.fluidModel.jwlModel.omega+1.0);

  A1           = iod.eqs.fluidModel.jwlModel.A1;
  A2           = iod.eqs.fluidModel.jwlModel.A2;
  R1           = iod.eqs.fluidModel.jwlModel.R1;
  R2           = iod.eqs.fluidModel.jwlModel.R2;
  rhoref       = iod.eqs.fluidModel.jwlModel.rhoref;
  R1r          = R1*rhoref;
  R2r          = R2*rhoref;

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

class VarFcnGasInGas : public VarFcn {

 protected:
  double gam;
  double gamp;    //gamma prime
  double gam1;
  double gamp1;
  double invgam1;
  double invgamp1;
  double Pstiff;
  double Pstiffp;

 public:
  VarFcnGasInGas(IoData &);
  ~VarFcnGasInGas() {}

  double getGamma() {return gam;}
  double getGamma1() {return gam1;}
  double getGammabis() {return gamp;}
  double getGammabis1() {return gamp1;}
  double getPressureConstant() {return Pstiff;}
  double getPressureConstantbis() {return Pstiffp;}

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[4]+Pstiff;
    else          return V[4]+Pstiffp; }
  void setDensity(double *V, double *Vorig){ V[0] = Vorig[0];}
  void setPressure(double *V, double *Vorig, double phi=0.0){ V[4] = Vorig[4];}
  
  double computeTemperature(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1*(V[4]+gam*Pstiff)/V[0];
    else         return invgamp1*(V[4]+gamp*Pstiffp)/V[0]; }
    
  double computeRhoEnergy(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    else         return invgamp1 * (V[4]+gamp*Pstiffp) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);}
    
  double computeRhoEpsilon(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1*(V[4]+gam*Pstiff); 
    else         return invgamp1*(V[4]+gamp*Pstiffp); }

  double computeSoundSpeed(double *V, double phi = 0.0) { 
    if (phi>=0.0){
      if(gam * (V[4]+Pstiff) / V[0]<0.0) 
        fprintf(stdout, "gam = %e, Pres = %e, pstiff = %e, rho = %e\n", gam,V[4],Pstiff,V[0]);
      return sqrt(gam * (V[4]+Pstiff) / V[0]); 
    }
    else{
      if(gamp * (V[4]+Pstiffp) / V[0]<0.0) 
        fprintf(stdout, "gamp = %e, Pres = %e, pstiffp = %e, rho = %e\n", gamp,V[4],Pstiffp,V[0]);
      return sqrt(gamp * (V[4]+Pstiffp) / V[0]);
    }
  }

  double computeSoundSpeed(double density, double entropy, double phi = 0.0) { 
    if (phi>=0.0)
      return sqrt(gam * entropy*pow(density,gam-1.0));
    else
      return sqrt(gamp * entropy*pow(density,gamp-1.0));
      //return sqrt(gamp * (V[4]+Pstiffp) / V[0]);
  }

  double computeEntropy(double density, double pressure, double phi = 0.0){
    if(phi>=0.0)
      return (pressure+Pstiff)/pow(density,gam);
    else
      return (pressure+Pstiffp)/pow(density,gamp);
  }

  double computeIsentropicPressure(double entropy, double density, double phi = 0.0){
    if(phi>=0.0) return entropy*pow(density,gam)-Pstiff;
    return entropy*pow(density,gamp)-Pstiffp;}

  double computeMachNumber(double *V, double phi = 0.0) { 
    if (phi>=0.0) return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)));
    else         return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gamp * (V[4]+Pstiffp))); }
    
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    double mach, machr2, popr, opmach, opmachr;
    mach = computeMachNumber(V, phi);
    machr2 = machr*machr;
    if (phi>=0.0) {
      popr = V[4]*gam*machr2;
      opmach = 1.0 + 0.5*gam1*mach*mach;
      //opmachr = 1.0 + 0.5*gam1*machr2;
      return V[4]*pow(opmach, gam*invgam1);
    }else{
      popr = V[4]*gamp*machr2;
      opmach = 1.0 + 0.5*gamp1*mach*mach;
      //opmachr = 1.0 + 0.5*gamp1*machr2;
      return V[4]*pow(opmach, gamp*invgamp1);
    }
  }

};
//--------------------------------------------------------------------------------------

inline
VarFcnGasInGas::VarFcnGasInGas(IoData &iod)
: VarFcn(iod) {

  type = GASINGAS;
  subType = NONE;

  gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;

  gamp = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
  gamp1 = gamp - 1.0;
  invgamp1 = 1.0/gamp1;
  Pstiffp = iod.eqs.fluidModel2.gasModel.pressureConstant;

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

class VarFcnLiquidInLiquid : public VarFcn {

protected:
  double Cv;
  double invCv;
  double alpha_water;
  double beta_water;
  double Pref_water;
  
  double Cvbis;
  double invCvbis;
  double alpha_waterbis;
  double beta_waterbis;
  double Pref_waterbis;

public:
  VarFcnLiquidInLiquid(IoData &);
  ~VarFcnLiquidInLiquid() {}
    
  double getCv()  {return Cv;}
  double getAlphaWater()  {return alpha_water;}
  double getBetaWater()  {return beta_water;}
  double getPrefWater()  {return Pref_water;}
  double getCvbis()  {return Cvbis;}
  double getAlphaWaterbis()  {return alpha_waterbis;}
  double getBetaWaterbis()  {return beta_waterbis;}
  double getPrefWaterbis()  {return Pref_waterbis;}

  double getPressure(double *V, double phi = 0.0) { 
    if (phi>=0.0) return Pref_water + alpha_water * pow(V[0], beta_water);
    else          return Pref_waterbis + alpha_waterbis * pow(V[0], beta_waterbis); }
  double checkPressure(double *V, double phi = 0.0) { 
    return 1.0;
    if (phi>=0.0) return Pref_water + alpha_water * pow(V[0], beta_water);
    else          return Pref_waterbis + alpha_waterbis * pow(V[0], beta_waterbis); }
  void   setDensity(double *V, double *Vorig)  { V[0] = Vorig[0];}
  void   setPressure(double *V, double *Vorig, double phi=0.0) { V[0] = Vorig[0];}

  double computeTemperature(double *V, double phi = 0.0) { return V[4]; }
  double computeRhoEnergy(double *V, double phi = 0.0) {
    if (phi>=0.0)
      return V[0] * Cv * V[4] + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    return V[0] * Cvbis * V[4] + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  
  double computeRhoEpsilon(double *V, double phi = 0.0) { 
    if (phi>=0.0) return V[0] * Cv * V[4];
    return  V[0] * Cvbis * V[4];}
    
  double computeSoundSpeed(double *V, double phi = 0.0) { 
    if (phi>=0.0) return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0));
    return  sqrt(alpha_waterbis * beta_waterbis * pow(V[0], beta_waterbis - 1.0)); }
    
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) { 
    if (phi>=0.0) return sqrt(alpha_water * beta_water * pow(density, beta_water - 1.0));
    return  sqrt(alpha_waterbis * beta_waterbis * pow(density, beta_waterbis - 1.0)); }
    
  double computeMachNumber(double *V, double phi = 0.0) {
    double c, u;
    if (phi>=0.0) c=sqrt(alpha_water * beta_water *pow(V[0], beta_water - 1.0));
    else        c=sqrt(alpha_waterbis * beta_waterbis *pow(V[0], beta_waterbis - 1.0));
    u=sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
    return u/c;
  }
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
  return 0.0;
  }

// Included (MB)
  void rstVar(IoData &iod) { rV(iod);}

};

//--------------------------------------------------------------------------------------

inline
VarFcnLiquidInLiquid::VarFcnLiquidInLiquid(IoData &iod)
: VarFcn(iod) {

  type = LIQUIDINLIQUID;
  subType = NONE;

  Cv=1.0;
  invCv=1.0/Cv;
  alpha_water = iod.eqs.fluidModel.liquidModel.alpha;
  beta_water  = iod.eqs.fluidModel.liquidModel.beta;
  Pref_water  = iod.eqs.fluidModel.liquidModel.Pref;
  
  Cvbis=1.0;
  invCvbis=1.0/Cvbis;
  alpha_waterbis = iod.eqs.fluidModel2.liquidModel.alpha;
  beta_waterbis  = iod.eqs.fluidModel2.liquidModel.beta;
  Pref_waterbis  = iod.eqs.fluidModel2.liquidModel.Pref;

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class VarFcnGasInLiquid : public VarFcn {

protected:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;

  double Cv;
  double invCv;
  double alpha_water;
  double beta_water;
  double Pref_water;

public:
  VarFcnGasInLiquid(IoData &);
  ~VarFcnGasInLiquid() {}
                                         
  double getGamma() {return gam;}
  double getGamma1() {return gam1;}
  double getPressureConstant() {return Pstiff;}
  double getCv()  {return Cv;}
  double getAlphaWater()  {return alpha_water;}
  double getBetaWater()  {return beta_water;}
  double getPrefWater()  {return Pref_water;}

  double getPressure(double *V, double phi = 0.0) { 
    if (phi>=0.0) 
      return Pref_water + alpha_water * pow(V[0], beta_water); 
    else
      return V[4];
  }
  double checkPressure(double *V, double phi = 0.0) { 
    if (phi>=0.0) 
      return Pref_water + alpha_water * pow(V[0], beta_water); 
    else
      return V[4]+Pstiff;
  }
  void   setDensity(double *V, double *Vorig)  { V[0] = Vorig[0];}
  void   setPressure(double *V, double *Vorig, double phi=0.0) {
    if (phi>=0.0) V[0] = Vorig[0];
    else          V[4] = Vorig[4]; }
    
  double computeTemperature(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[4];
    return invgam1*(V[4]+gam*Pstiff)/V[0]; }
    
  double computeRhoEnergy(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[0] * Cv * V[4] + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);}
    
  double computeRhoEpsilon(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[0] * Cv * V[4];
    return invgam1*(V[4]+gam*Pstiff);}
    
  double computeSoundSpeed(double *V, double phi = 0.0) {
    if (phi>=0.0){
      if (alpha_water * beta_water * pow(V[0], beta_water - 1.0) < 0.0) 
        fprintf(stdout, "c2_water = %e\n", alpha_water * beta_water * pow(V[0], beta_water - 1.0));
      return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0));
    }
    if (gam * (V[4]+Pstiff) / V[0] < 0.0) 
      fprintf(stdout, "c2_air = %e - P = %e - rho = %e\n", gam * (V[4]+Pstiff) / V[0], V[4], V[0]);
    return sqrt(gam * (V[4]+Pstiff) / V[0]);}
    
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) {
    if (phi>=0.0)
      return sqrt(alpha_water * beta_water * pow(density, beta_water - 1.0));
    return sqrt(gam * entropy*pow(density,gam-1.0));}

  double computeEntropy(double density, double pressure, double phi = 0.0){
    if(phi>=0.0) return 0.0;
    return (pressure+Pstiff)/pow(density,gam);}

  double computeIsentropicPressure(double entropy, double density, double phi = 0.0){
    if(phi>=0.0) return Pref_water + alpha_water * pow(density, beta_water);
    return entropy*pow(density,gam)-Pstiff;
  }

  double computeMachNumber(double *V, double phi = 0.0) {
    double c, u;
    if (phi>=0.0) c=sqrt(alpha_water * beta_water *pow(V[0], beta_water - 1.0));
    else c = sqrt(gam*(V[4]+Pstiff)/V[0]);
    u = sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
    return u/c;
  }
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    return 0.0;
  }

};

//--------------------------------------------------------------------------------------

inline
VarFcnGasInLiquid::VarFcnGasInLiquid(IoData &iod)
: VarFcn(iod) {

  type = GASINLIQUID; 
  subType = NONE;

  if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
     iod.eqs.fluidModel2.fluid == FluidModelData::GAS){
    gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    gam1 = gam -1.0;
    invgam1 = 1.0/gam1;
    Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant;

    Cv=1.0;
    invCv=1.0/Cv;

    alpha_water = iod.eqs.fluidModel.liquidModel.alpha;
    beta_water  = iod.eqs.fluidModel.liquidModel.beta;
    Pref_water  = iod.eqs.fluidModel.liquidModel.Pref;

  }else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
           iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID){
    gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    gam1 = gam -1.0;
    invgam1 = 1.0/gam1;
    Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;

    Cv=1.0;
    invCv=1.0/Cv;

    alpha_water = iod.eqs.fluidModel2.liquidModel.alpha;
    beta_water  = iod.eqs.fluidModel2.liquidModel.beta;
    Pref_water  = iod.eqs.fluidModel2.liquidModel.Pref;

  }

}

//------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

class VarFcnJWLInGas : public VarFcn {

 protected:
  double gam, gam1, invgam1;
  double Pstiff;

  double omega, invomega, omegap1, invomegap1;
  double A1,A2,R1,R2,rhoref;
  double R1r,R2r;

 public:
  VarFcnJWLInGas(IoData &);
  ~VarFcnJWLInGas() {}

  double getGamma() {return gam;}
  double getGamma1() {return gam1;}
  double getPressureConstant() {return Pstiff;}

  double getOmega()       {return omega;}
  double getOmegap1()     {return omegap1;}
  double getInvOmega()    {return invomega;}
  double getInvOmegap1()  {return invomegap1;}
  double getA1()          {return A1;}
  double getA2()          {return A2;}
  double getR1()          {return R1;}
  double getR2()          {return R2;}
  double getRhoref()      {return rhoref;}
  double getR1r()         {return R1r;}
  double getR2r()         {return R2r;}

  double computeFrho(double *V, double phi = 0.0) {
    return A1*(1.0-omega*V[0]/R1r)*exp(-R1r/V[0])
         + A2*(1.0-omega*V[0]/R2r)*exp(-R2r/V[0]);
  }
  double computeFrhop(double *V, double phi = 0.0) {
    double rho2 = V[0]*V[0];
    return A1*(-omega/R1r+(1.0-omega*V[0]/R1r)*R1r/rho2)*exp(-R1r/V[0])
         + A2*(-omega/R2r+(1.0-omega*V[0]/R2r)*R2r/rho2)*exp(-R2r/V[0]);
  }
  double computeFrho(double rho, double phi = 0.0) {
    return A1*(1.0-omega*rho/R1r)*exp(-R1r/rho)
         + A2*(1.0-omega*rho/R2r)*exp(-R2r/rho);
  }
  double computeFrhop(double rho, double phi = 0.0) {
    double rho2 = rho*rho;
    return A1*(-omega/R1r+(1.0-omega*rho/R1r)*R1r/rho2)*exp(-R1r/rho)
         + A2*(-omega/R2r+(1.0-omega*rho/R2r)*R2r/rho2)*exp(-R2r/rho);
  }

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[4]+Pstiff;
    else          return V[4] - (computeFrho(V,phi) - computeFrhop(V,phi)*V[0])*invomegap1; }
  void setDensity(double *V, double *Vorig){ V[0] = Vorig[0];}
  void setPressure(double *V, double *Vorig, double phi=0.0){ V[4] = Vorig[4];}
  
  double computeTemperature(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1*(V[4]+gam*Pstiff)/V[0];
    else         return invomega*(V[4]-computeFrho(V,phi))/V[0]; }
    
  double computeRhoEnergy(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    else         return invomega * (V[4]-computeFrho(V,phi)) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
    
  double computeRhoEpsilon(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invgam1*(V[4]+gam*Pstiff); 
    else         return invomega * (V[4]-computeFrho(V,phi)); }
     
  double computeSoundSpeed(double *V, double phi = 0.0) { 
    if (phi>=0.0){
      if(gam * (V[4]+Pstiff) / V[0]<0.0) 
        fprintf(stdout, "gam = %e, Pres = %e, pstiff = %e, rho = %e\n", gam,V[4],Pstiff,V[0]);
      return sqrt(gam * (V[4]+Pstiff) / V[0]); 
    }
    else{
      if ((omegap1*V[4] - computeFrho(V,phi) + V[0]*computeFrhop(V,phi))/V[0] < 0.0){
      fprintf(stdout, "computeSoundSpeed^2 = %e\n", (omegap1*V[4] - computeFrho(V,phi) + V[0]*computeFrhop(V,phi))/V[0]);
      fprintf(stdout, "computeFrho = %e abd compureFrhop = %e\n", computeFrho(V,phi), computeFrhop(V,phi));
      fprintf(stdout, "V = %e %e %e\n", V[0],V[1],V[4]);
      //return 0.0;
      }
      return sqrt((omegap1*V[4] - computeFrho(V,phi) + V[0]*computeFrhop(V,phi))/V[0]);
    }
  }

  double computeSoundSpeed(double density, double entropy, double phi = 0.0) { 
    if (phi>=0.0){
      //fprintf(stderr, "gas side=%e\n", gam * entropy*pow(density,gam-1.0));
      return sqrt(gam * entropy*pow(density,gam-1.0));
    }else{
      double pressure = entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
      //fprintf(stderr, "jwl side=%e (pressure = %e)\n", (omegap1*pressure - computeFrho(density,phi) + density*computeFrhop(density,phi))/density, pressure);
      return sqrt((omegap1*pressure - computeFrho(density,phi) + density*computeFrhop(density,phi))/density);
    }
  }

  double computeEntropy(double density, double pressure, double phi = 0.0){
    if(phi>=0.0) return (pressure+Pstiff)/pow(density,gam);
    return (pressure - A1*exp(-R1r/density) - A2*exp(-R2r/density))/pow(density,omegap1);
  }

  double computeIsentropicPressure(double entropy, double density, double phi = 0.0){
    if(phi>=0.0) return entropy*pow(density,gam)-Pstiff;
    return entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
  }
    
  double computeMachNumber(double *V, double phi = 0.0) { 
    if (phi>=0.0) return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)));
    else {/*fprintf(stderr, "computeMachNumber = %e\n", computeSoundSpeed(V,phi));*/        return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]))/computeSoundSpeed(V,phi); }
  }
    
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    fprintf(stderr, "*** Error: computeTotalPressure not yet implemented for VarFcnJWL\n");
    return 0.0;
  }

};
//--------------------------------------------------------------------------------------

inline
VarFcnJWLInGas::VarFcnJWLInGas(IoData &iod)
: VarFcn(iod) {

// Assumes that fluidModel1 is the PGas
//          and fluidModel2 is the JWL
  type = JWLINGAS;
  subType = NONE;

  gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;


  omega        = iod.eqs.fluidModel2.jwlModel.omega;
  omegap1      = omega + 1.0;
  invomega     = 1.0/omega;
  invomegap1   = 1.0/omegap1;

  A1           = iod.eqs.fluidModel2.jwlModel.A1;
  A2           = iod.eqs.fluidModel2.jwlModel.A2;
  R1           = iod.eqs.fluidModel2.jwlModel.R1;
  R2           = iod.eqs.fluidModel2.jwlModel.R2;
  rhoref       = iod.eqs.fluidModel2.jwlModel.rhoref;
  R1r          = R1*rhoref;
  R2r          = R2*rhoref;

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

class VarFcnJWLInJWL : public VarFcn {

 protected:
  double omega, invomega, omegap1, invomegap1;
  double A1,A2,R1,R2,rhoref;
  double R1r,R2r;

  double omegap, invomegap, omegap1p, invomegap1p;
  double A1p,A2p,R1p,R2p,rhorefp;
  double R1rp,R2rp;

 public:
  VarFcnJWLInJWL(IoData &);
  ~VarFcnJWLInJWL() {}

  double getOmega(double phi)   {
    if(phi>=0.0) return omega;
    else         return omegap; }
  double getOmegap1(double phi) {
    if(phi>-0.0) return omegap1;
    else         return omegap1p; }

  double getOmega()       {return omega;}
  double getOmegap1()     {return omegap1;}
  double getInvOmega()    {return invomega;}
  double getInvOmegap1()  {return invomegap1;}
  double getA1()          {return A1;}
  double getA2()          {return A2;}
  double getR1()          {return R1;}
  double getR2()          {return R2;}
  double getRhoref()      {return rhoref;}
  double getR1r()         {return R1r;}
  double getR2r()         {return R2r;}

  double getOmegabis()       {return omegap;}
  double getOmegap1bis()     {return omegap1p;}
  double getInvOmegabis()    {return invomegap;}
  double getInvOmegap1bis()  {return invomegap1p;}
  double getA1bis()          {return A1p;}
  double getA2bis()          {return A2p;}
  double getR1bis()          {return R1p;}
  double getR2bis()          {return R2p;}
  double getRhorefbis()      {return rhorefp;}
  double getR1rbis()         {return R1rp;}
  double getR2rbis()         {return R2rp;}

  double computeFrho(double *V, double phi = 0.0) {
    if(phi>=0.0)
      return A1*(1.0-omega*V[0]/R1r)*exp(-R1r/V[0])
           + A2*(1.0-omega*V[0]/R2r)*exp(-R2r/V[0]);
    else
      return A1p*(1.0-omegap*V[0]/R1rp)*exp(-R1rp/V[0])
           + A2p*(1.0-omegap*V[0]/R2rp)*exp(-R2rp/V[0]);
  }
  double computeFrhop(double *V, double phi = 0.0) {
    double rho2 = V[0]*V[0];
    if(phi>=0.0)
      return A1*(-omega/R1r+(1.0-omega*V[0]/R1r)*R1r/rho2)*exp(-R1r/V[0])
           + A2*(-omega/R2r+(1.0-omega*V[0]/R2r)*R2r/rho2)*exp(-R2r/V[0]);
    else
      return A1p*(-omegap/R1rp+(1.0-omegap*V[0]/R1rp)*R1rp/rho2)*exp(-R1rp/V[0])
           + A2p*(-omegap/R2rp+(1.0-omegap*V[0]/R2rp)*R2rp/rho2)*exp(-R2rp/V[0]);
  }
  double computeFrho(double rho, double phi = 0.0) {
    if(phi>=0.0)
      return A1*(1.0-omega*rho/R1r)*exp(-R1r/rho)
           + A2*(1.0-omega*rho/R2r)*exp(-R2r/rho);
    else
      return A1p*(1.0-omegap*rho/R1rp)*exp(-R1rp/rho)
           + A2p*(1.0-omegap*rho/R2rp)*exp(-R2rp/rho);
  }
  double computeFrhop(double rho, double phi = 0.0) {
    double rho2 = rho*rho;
    if(phi>=0.0)
      return A1*(-omega/R1r+(1.0-omega*rho/R1r)*R1r/rho2)*exp(-R1r/rho)
           + A2*(-omega/R2r+(1.0-omega*rho/R2r)*R2r/rho2)*exp(-R2r/rho);
    else
      return A1p*(-omegap/R1rp+(1.0-omegap*rho/R1rp)*R1rp/rho2)*exp(-R1rp/rho)
           + A2p*(-omegap/R2rp+(1.0-omegap*rho/R2rp)*R2rp/rho2)*exp(-R2rp/rho);
  }

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) {
    if (phi>=0.0) return V[4] - (computeFrho(V,phi) - computeFrhop(V,phi)*V[0])*invomegap1;
    else          return V[4] - (computeFrho(V,phi) - computeFrhop(V,phi)*V[0])*invomegap1p; }
  void setDensity(double *V, double *Vorig){ V[0] = Vorig[0];}
  void setPressure(double *V, double *Vorig, double phi=0.0){ V[4] = Vorig[4];}
  
  double computeTemperature(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invomega *(V[4]-computeFrho(V,phi))/V[0]; 
    else          return invomegap*(V[4]-computeFrho(V,phi))/V[0]; }
    
  double computeRhoEnergy(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invomega  * (V[4]-computeFrho(V,phi)) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    else          return invomegap * (V[4]-computeFrho(V,phi)) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
    
  double computeRhoEpsilon(double *V, double phi = 0.0) { 
    if (phi>=0.0) return invomega * (V[4]-computeFrho(V,phi));
    else          return invomegap* (V[4]-computeFrho(V,phi)); }
     
  double computeSoundSpeed(double *V, double phi = 0.0) { 
    if (phi>=0.0)
      return sqrt((omegap1 *V[4] - computeFrho(V,phi) + V[0]*computeFrhop(V,phi))/V[0]);
    else
      return sqrt((omegap1p*V[4] - computeFrho(V,phi) + V[0]*computeFrhop(V,phi))/V[0]);
  }
    
  double computeSoundSpeed(double density, double entropy, double phi = 0.0) { 
    if (phi>=0.0){
      double pressure = entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
      return sqrt((omegap1 *pressure - computeFrho(density,phi) + density*computeFrhop(density,phi))/density);
    }else{
      double pressure = entropy*pow(density,omegap1p)+A1p*exp(-R1rp/density)+A2p*exp(-R2rp/density);
      return sqrt((omegap1p*pressure - computeFrho(density,phi) + density*computeFrhop(density,phi))/density);
    }
  }

  double computeEntropy(double density, double pressure, double phi = 0.0){
    if(phi>=0.0) return (pressure - A1*exp(-R1r/density) - A2*exp(-R2r/density))/pow(density,omegap1);
    return (pressure - A1p*exp(-R1rp/density) - A2p*exp(-R2rp/density))/pow(density,omegap1p);
  }

  double computeIsentropicPressure(double entropy, double density, double phi = 0.0){
    if(phi>=0.0) return entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
    return entropy*pow(density,omegap1p)+A1p*exp(-R1rp/density)+A2p*exp(-R2rp/density);
  }
    
  double computeMachNumber(double *V, double phi = 0.0) { 
    return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]))/computeSoundSpeed(V,phi);
  }
    
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    fprintf(stderr, "*** Error: computeTotalPressure not yet implemented for VarFcnJWLInJWL\n");
    return 0.0;
  }

};
//--------------------------------------------------------------------------------------

inline
VarFcnJWLInJWL::VarFcnJWLInJWL(IoData &iod)
: VarFcn(iod) {

  type = JWLINJWL;
  subType = NONE;

  omega        = iod.eqs.fluidModel.jwlModel.omega;
  omegap1      = omega + 1.0;
  invomega     = 1.0/omega;
  invomegap1   = 1.0/omegap1;

  A1           = iod.eqs.fluidModel.jwlModel.A1;
  A2           = iod.eqs.fluidModel.jwlModel.A2;
  R1           = iod.eqs.fluidModel.jwlModel.R1;
  R2           = iod.eqs.fluidModel.jwlModel.R2;
  rhoref       = iod.eqs.fluidModel.jwlModel.rhoref;
  R1r          = R1*rhoref;
  R2r          = R2*rhoref;

  omegap       = iod.eqs.fluidModel2.jwlModel.omega;
  omegap1p     = omegap + 1.0;
  invomegap    = 1.0/omegap;
  invomegap1p  = 1.0/omegap1p;

  A1p          = iod.eqs.fluidModel2.jwlModel.A1;
  A2p          = iod.eqs.fluidModel2.jwlModel.A2;
  R1p          = iod.eqs.fluidModel2.jwlModel.R1;
  R2p          = iod.eqs.fluidModel2.jwlModel.R2;
  rhorefp      = iod.eqs.fluidModel2.jwlModel.rhoref;
  R1rp         = R1p*rhorefp;
  R2rp         = R2p*rhorefp;

}

//-----------------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitive(SVec<double,dim> &U, SVec<double,dim> &V, Vec<double> *Phi)
{
  if (Phi){
    for (int i=0; i<U.size(); ++i){
      conservativeToPrimitive(U[i], V[i], (*Phi)[i]);
    }
  }else{
    for (int i=0; i<U.size(); ++i)
      conservativeToPrimitive(U[i], V[i]);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::conservativeToPrimitiveDerivative(SVec<double,dim> &U, SVec<double,dim> &dU, SVec<double,dim> &V, SVec<double,dim> &dV, Vec<double> *Phi)
{

  if (Phi){
    for (int i=0; i<U.size(); ++i){
      conservativeToPrimitiveDerivative(U[i], dU[i], V[i], dV[i], (*Phi)[i]);
    }
  }else{
    for (int i=0; i<U.size(); ++i){
      conservativeToPrimitiveDerivative(U[i], dU[i], V[i], dV[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitive(DistSVec<double,dim> &U, DistSVec<double,dim> &V, DistVec<double> *Phi)
{
 
  int numLocSub = U.numLocSub();
 
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        conservativeToPrimitive(u[i], v[i], phi[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        conservativeToPrimitive(u[i], v[i]);
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::conservativeToPrimitiveDerivative(DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistVec<double> *Phi)
{

  int numLocSub = U.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double (*dv)[dim] = dV.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        conservativeToPrimitiveDerivative(u[i], du[i], v[i], dv[i], phi[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        conservativeToPrimitiveDerivative(u[i], du[i], v[i], dv[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(SVec<double,dim> &V, SVec<double,dim> &U, Vec<double> *Phi)
{
  if (Phi){
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(V[i], U[i], (*Phi)[i]);
  }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(V[i], U[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(DistSVec<double,dim> &V, DistSVec<double,dim> &U, DistVec<double> *Phi)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservative(v[i], u[i], phi[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservative(v[i], u[i]);
    }
  }

}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::primitiveToConservativeDerivative(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<double> *Phi)
{

  if (Phi){
    for (int i=0; i<U.size(); ++i)
      primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i], (*Phi)[i]);
  }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i]);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::primitiveToConservativeDerivative(DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistVec<double> *Phi)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double (*dv)[dim] = dV.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i], phi[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::updatePhaseChange(DistSVec<double,dim> &V, DistSVec<double,dim> &U,
          DistVec<double> &Phi, DistVec<double> &Phin,
          DistSVec<double,dim> *Vgf, DistVec<double> *Vgfweight,
          DistExactRiemannSolver<dim> *riemann)
{
  if (riemann->DoUpdatePhase())
    // the solution of the riemann problem is used to replace values of a node
    // that changed nature (fluid1 to fluid2 or vice versa)
    // **** GFMPAR-like ****
    updatePhaseChange(V, U, Phi, Phin, riemann->getRiemannUpdate(), riemann->getRiemannWeight());

  else if (Vgf && Vgfweight)
    // an extrapolation is used to replace values of a node
    // that changed nature (fluid1 to fluid2 or vice versa)
    // **** GFMPAR-variation ****
    updatePhaseChange(V, U, Phi, Phin, Vgf, Vgfweight);

  else
    // no solution of the riemann problem was computed and we just use
    // the values that we have to convert back to conservative variables
    // **** GFMP-like ****
    primitiveToConservative(V, U, &Phi);

}

//------------------------------------------------------------------------------
template<int dim>
void VarFcn::updatePhaseChange(DistSVec<double,dim> &V, DistSVec<double,dim> &U,
                               DistVec<double> &Phi, DistVec<double> &Phin,
                               DistSVec<double,dim> *Riemann, DistVec<double> *weight)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double *phi = Phi.subData(iSub);
    double *phin = Phin.subData(iSub);
    double (*r)[dim] = (*Riemann).subData(iSub);
    double *w = 0;
    bool change = false;
    if(weight){
      w = (*weight).subData(iSub);
      for ( int i=0; i<U.subSize(iSub); ++i){
        change = updatePhaseChange(v[i],u[i],phi[i],phin[i],r[i],w[i]);
      }
    }else{
      for ( int i=0; i<U.subSize(iSub); ++i){
        change = updatePhaseChange(v[i],u[i],phi[i],phin[i],r[i]);
      }
    }
  }


}

//------------------------------------------------------------------------------
template<int dim>
void VarFcn::computeNewPrimitive(SVec<double,dim> &U, SVec<double,dim> &V)
{
  for (int i=0; i<U.size(); ++i)
    computeNewPrimitive(U[i], V[i]);
}

//------------------------------------------------------------------------------


#endif
