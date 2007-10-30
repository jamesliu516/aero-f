#ifndef _VAR_FCN_H_
#define _VAR_FCN_H_

#include <DistVector.h>
#include <Vector3D.h>

#include <math.h>
#include <complex.h>
typedef complex<double> bcomp;

#include <IoData.h>

//--------------------------------------------------------------------------
class VarFcn {

public:
  
  const char** pname;
  enum Type{ GAS = 0, LIQUID = 1, GASINGAS = 2, GASINLIQUID = 3, LIQUIDINLIQUID = 4} type;
  enum SubType { NONE = 0, IDEAL = 1, STIFFENED = 2} subType;
  Vec3D meshVel;

  double pmin;
  double pminp;

  double gravity;
  double ngravity[3];
  
  VarFcn(IoData &iod); 
  VarFcn() { meshVel = 0.0; fprintf(stderr, "MeshVel initialized: %e %e %e\n", meshVel[0], meshVel[1], meshVel[2]);}
  ~VarFcn() {}

  double *getMeshVel()  { return meshVel.v; }

  template<int dim>
  void conservativeToPrimitive(SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0);
  template<int dim>
  void conservativeToPrimitive(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0);
  template<int dim>
  void primitiveToConservative(SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0, Vec<double> * = 0, SVec<double, dim> * = 0);
  template<int dim>
  void primitiveToConservative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0, DistVec<double> * = 0, DistSVec<double,dim> * = 0);

  template<int dim>
  void mergeFluidVectors(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, Vec<double> *);
  template<int dim>
  void mergeFluidVectors(DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> *);
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
  void VerificationGasEuler(double, double, double, double*, double*);

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
  void VerificationGasSA(double, double , double *, double *);

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
  void VerificationGasKE(double, double , double *, double *);

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
  void VerificationLiquidEuler(double, double, double, double, double, double *, double *);
public:

  template<int dim>
  void computeNewPrimitive(SVec<double,dim> &, SVec<double,dim> &);
  virtual void computeNewPrimitive(double *, double *){}
  virtual void computeOldPrimitive(double *, double *){}

  virtual void conservativeToPrimitive(double *, double *, double = 0.0) = 0; 
  virtual void conservativeToPrimitiveVerification(double *, double *, double = 0.0) = 0; 
  virtual void primitiveToConservative(double *, double *, double = 0.0, double* = 0, double* = 0) = 0;
  bool doVerification(){ return !(pmin<0 && pminp<0); }
  void mergeFluidVectors(double *V1, double *V2, double *V, double phi) 
  { 
   if (phi >= 0.0) V  = V1;
   else            V  = V2;
  };
  
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
        fprintf(stderr, "*** Warning:  getCv Function not defined\n");
        return 0.0; }
  virtual double getAlphaWater() {
        fprintf(stderr, "*** Warning: getAlphaWater Function not defined\n");
        return 0.0;}
  virtual double getAlphaWaterbis() {
        fprintf(stderr, "*** Warning: getAlphaWater Function not defined\n");
        return 0.0;}
  virtual double getBetaWater() {
        fprintf(stderr, "*** Warning: getBetaWater Function not defined\n");
        return 0.0;}
  virtual double getBetaWaterbis() {
        fprintf(stderr, "*** Warning: getBetaWater Function not defined\n");
        return 0.0;}
  virtual double getPrefWater() {
        fprintf(stderr, "*** Warning: getPrefWater Function not defined\n");
        return 0.0;}
  virtual double getPrefWaterbis() {
        fprintf(stderr, "*** Warning: getPrefWater Function not defined\n");
        return 0.0;}
                                                               
  virtual Vec3D getVelocity(double *V) { return Vec3D(V[1], V[2], V[3]); }
  virtual double getDensity(double *V) { return V[0]; }
  virtual double getVelocityX(double *V) {return V[1];}
  virtual double getVelocityY(double *V) {return V[2];}
  virtual double getVelocityZ(double *V) {return V[3];}
  virtual double getPressure(double *V, double phi = 0.0) {return 0.0;}
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
  void primitiveToConservativeDerivative(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0, Vec<double> * = 0, SVec<double, dim> * = 0);
  template<int dim>
  void primitiveToConservativeDerivative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0, DistVec<double> * = 0, DistSVec<double,dim> * = 0);

  virtual void conservativeToPrimitiveDerivative(double *, double *, double *, double *, double = 0.0) = 0;
  virtual void primitiveToConservativeDerivative(double *, double *, double *, double *, double = 0.0, double* = 0, double* = 0) = 0;

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
  virtual double getDerivativeOfBetaRef() { return 0.0; }
  virtual double getDerivativeOfPressureConstant() { return 0.0; }
  virtual void rstVar(IoData &iod) {}
  virtual void rV(IoData &iod) {  pmin  = iod.eqs.fluidModel.pmin;
                                  pminp = iod.eqs.volumes.fluidModel2.pmin; }
                                                                    
};
//------------------------------------------------------------------------------

class VarFcnPerfectGas : public VarFcn {

protected:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;  

// Included (MB)
  double k2;
  double dPstiff;  

public:
  VarFcnPerfectGas(IoData &);
  ~VarFcnPerfectGas() {}
                            
  double getGamma()  {return gam;}
  double getGamma1()  {return gam1;}
  double getPressureConstant() {return Pstiff;}

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double computeTemperature(double *V, double phi = 0.0 ) {
        if (isnan(1.0/V[0])) {
                fprintf(stderr, "ERROR*** computeTemp\n");
                exit(1);
        }
        return invgam1 * (V[4]+gam*Pstiff) / V[0]; }
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
  double computeTotalPressure(double machr, double* V, double phi = 0.0) {
    double mach = computeMachNumber(V);
    double machr2 = machr*machr;
    double popr = V[4]*gam*machr2;
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    //double opmachr = 1.0 + 0.5*gam1*machr2;
    return V[4]*pow(opmach, gam*invgam1);
  }

// Included (MB)
  double getVelocityNorm(double *V) { return sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]); }
  double getDerivativeOfVelocityNorm(double *V, double *dV) { return (V[1]*dV[1]+V[2]*dV[2]+V[3]*dV[3])/sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]); }
  double computeDerivativeOfTemperature(double *V, double *dV, double phi = 0.0) { return invgam1 * ( dV[4]*V[0] - V[4]*dV[0] ) / ( V[0]*V[0] ); }
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
  double getDerivativeOfBetaRef() { return k2; }
  double getDerivativeOfPressureConstant() {return dPstiff;}
  void rstVar(IoData &iod) { dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach)); rV(iod);}

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
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;

// Included (MB)
  k2 = iod.prec.k2;
  dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));

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
  double computeTemperature(double *V, double phi = 0.0) { return V[4]; }
  double computeRhoEnergy(double *V, double phi = 0.0) {
    return V[0] * Cv * V[4] + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V, double phi = 0.0) { return V[0] * Cv * V[4]; }
  double computeSoundSpeed(double *V, double phi = 0.0) { return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0)); }
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

// Included (MB)
  double dPstiff;  

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
    if (phi>=0.0) return sqrt(gam * (V[4]+Pstiff) / V[0]); 
    else         return sqrt(gamp * (V[4]+Pstiffp) / V[0]);}
    
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

// Included (MB)
  double getDerivativeOfPressureConstant() {return dPstiff;}
  void rstVar(IoData &iod) { dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));  rV(iod);}

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
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;

  gamp = iod.eqs.volumes.fluidModel2.gasModel.specificHeatRatio;
  gamp1 = gamp - 1.0;
  invgamp1 = 1.0/gamp1;
  Pstiffp = iod.eqs.volumes.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;

// Included (MB)
  //dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));
  dPstiff = 0.0;

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
  alpha_waterbis = iod.eqs.volumes.fluidModel2.liquidModel.alpha;
  beta_waterbis  = iod.eqs.volumes.fluidModel2.liquidModel.beta;
  Pref_waterbis  = iod.eqs.volumes.fluidModel2.liquidModel.Pref;

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

// Included (MB)
  double dPstiff;  

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
    if (phi>=0.0) return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0));
    return sqrt(gam * (V[4]+Pstiff) / V[0]);}
    
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

// Included (MB)
  double getDerivativeOfPressureConstant() {return dPstiff;}
  void rstVar(IoData &iod) { dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));  rV(iod);}

};

//--------------------------------------------------------------------------------------

inline
VarFcnGasInLiquid::VarFcnGasInLiquid(IoData &iod)
: VarFcn(iod) {

  type = GASINLIQUID; 
  subType = NONE;

  //gam  = 1.0  +iod.eqs.fluidModel2.gasModel.idealGasConstant/iod.eqs.fluidModel.liquidModel.Cv;
  gam = iod.eqs.volumes.fluidModel2.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;

  Cv=1.0;
  invCv=1.0/Cv;

  alpha_water = iod.eqs.fluidModel.liquidModel.alpha;
  beta_water  = iod.eqs.fluidModel.liquidModel.beta;
  Pref_water  = iod.eqs.fluidModel.liquidModel.Pref;

// Included (MB)
  dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/(iod.bc.inlet.pressure*iod.ref.rv.pressure)*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));

}

//------------------------------------------------------------------------------
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
      for (int i=0; i<U.subSize(iSub); ++i){
        if (u[i][0]<0.0 || isnan(u[i][0])) fprintf(stderr, "in subD %d, rho[%d] is %f\n", iSub, i, u[i][0]);
        conservativeToPrimitive(u[i], v[i], phi[i]);
      }
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
void VarFcn::mergeFluidVectors(SVec<double, dim> &V1, SVec<double,dim> &V2, SVec<double,dim> &V, Vec<double> *Phi)
{
  if (Phi){
    for (int i=0; i<V.size(); ++i) {
      if ((*Phi)[i] >= 0.0) V[i]  = V1[i];
      else                  V[i]  = V2[i];
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::mergeFluidVectors(DistSVec<double,dim> &V1, DistSVec<double,dim> &V2, DistSVec<double,dim> &V, DistVec<double> *Phi)
{

  int numLocSub = V.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*v1)[dim] = V1.subData(iSub);
    double (*v2)[dim] = V2.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<V.subSize(iSub); ++i)
        mergeFluidVectors(v1[i], v2[i], v[i], phi[i]);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(SVec<double,dim> &V, SVec<double,dim> &U, Vec<double> *Phi, Vec<double> *Phi1, SVec<double,dim> *Vgf)
{

  if (!Phi1) {
    if (Phi){
      for (int i=0; i<U.size(); ++i)
        primitiveToConservative(V[i], U[i], (*Phi)[i]);
    }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(V[i], U[i]);
    }
  }
  else {
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(V[i], U[i], (*Phi)[i], &((*Phi1)[i]), (*Vgf)[i]);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::primitiveToConservativeDerivative(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<double> *Phi, Vec<double> *Phi1, SVec<double,dim> *Vgf)
{

  if (!Phi1) {
    if (Phi){
      for (int i=0; i<U.size(); ++i)
        primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i], (*Phi)[i]);
    }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i]);
    }
  }
  else {
    for (int i=0; i<U.size(); ++i)
      primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i], (*Phi)[i], &((*Phi1)[i]), (*Vgf)[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(DistSVec<double,dim> &V, DistSVec<double,dim> &U, DistVec<double> *Phi, DistVec<double> *Phi1, DistSVec<double,dim> *Vgf)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if (!Phi1) {
      if (Phi){
        double *phi = (*Phi).subData(iSub);
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(v[i], u[i], phi[i]);
      }else{
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(v[i], u[i]);
      }
    }
    else {
        double *phi  = (*Phi).subData(iSub);
        double *phi1 = (*Phi1).subData(iSub);
        double (*vgf)[dim]  = (*Vgf).subData(iSub);
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(v[i], u[i], phi[i], &(phi1[i]), vgf[i] );
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void VarFcn::primitiveToConservativeDerivative(DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistVec<double> *Phi, DistVec<double> *Phi1, DistSVec<double,dim> *Vgf)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double (*dv)[dim] = dV.subData(iSub);
    if (!Phi1) {
      if (Phi){
        double *phi = (*Phi).subData(iSub);
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i], phi[i]);
      }else{
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i]);
      }
    }
    else {
        double *phi  = (*Phi).subData(iSub);
        double *phi1 = (*Phi1).subData(iSub);
        double (*vgf)[dim]  = (*Vgf).subData(iSub);
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(v[i], dv[i], u[i], du[i], phi[i], &(phi1[i]), vgf[i] );
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
void VarFcn::computeNewPrimitive(SVec<double,dim> &U, SVec<double,dim> &V){
 
  for (int i=0; i<U.size(); ++i)
    computeNewPrimitive(U[i], V[i]);
}

//------------------------------------------------------------------------------


#endif
