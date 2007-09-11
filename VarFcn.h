#ifndef _VAR_FCN_H_
#define _VAR_FCN_H_

#include <DistVector.h>
#include <Vector3D.h>
//#include <DistExactRiemannSolver.h>

#include <math.h>
#include <complex.h>
typedef complex<double> bcomp;

#include <IoData.h>

template<int dim> class DistExactRiemannSolver;

//--------------------------------------------------------------------------
class VarFcn {

public:
  
  const char** pname;
  enum Type{ GAS = 0, LIQUID = 1, GASINGAS = 2, GASINLIQUID = 3, LIQUIDINLIQUID = 4} type;
  enum SubType { NONE = 0, IDEAL = 1, STIFFENED = 2} subType;
  Vec3D meshVel;

  double pmin;
  double pminp;
  bool verif_clipping;
  bool node_change;

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
  /*template<int dim>
  void primitiveToConservative(SVec<double,dim> &, SVec<double,dim> &, Vec<double> * = 0, Vec<double> * = 0, SVec<double, dim> * = 0);
  template<int dim>
  void primitiveToConservative(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> * = 0, DistVec<double> * = 0, DistSVec<double,dim> * = 0);
	*/
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
  //GAS-EULER
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

  //GAS_SA
  void conservativeToPrimitiveGasSA(double, double *, double *);
  void primitiveToConservativeGasSA(double *, double *);
  void computedUdVGasSA(double, double *, double *);
  void computedVdUGasSA(double, double *, double *);
  void multiplyBydVdUGasSA(double *, double *, double *);
  void postMultiplyBydUdVGasSA(double, double *, double *, double *);
  void postMultiplyBydVdUGasSA(double, double *, double *, double *);
  void preMultiplyBydUdVGasSA(double, double *, double *, double *);
  int VerificationGasSA(int, double, double , double *, double *);

  //GAS_KE
  void conservativeToPrimitiveGasKE(double, double *, double *);
  void primitiveToConservativeGasKE(double *, double *);
  void computedUdVGasKE(double, double *, double *);
  void computedVdUGasKE(double, double *, double *);
  void multiplyBydVdUGasKE(double *, double *, double *);
  void postMultiplyBydUdVGasKE(double, double *, double *, double *);
  void postMultiplyBydVdUGasKE(double, double *, double *, double *);
  void preMultiplyBydUdVGasKE(double, double *, double *, double *);
  int VerificationGasKE(int, double, double , double *, double *);

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
public:

  template<int dim>
  void computeNewPrimitive(SVec<double,dim> &, SVec<double,dim> &);
  virtual void computeNewPrimitive(double *, double *){}
  virtual void computeOldPrimitive(double *, double *){}

  virtual void conservativeToPrimitive(double *, double *, double = 0.0) = 0; 
  virtual int  conservativeToPrimitiveVerification(int, double *, double *, double = 0.0) = 0; 
  //virtual void primitiveToConservative(double *, double *, double = 0.0, double* = 0, double* = 0, int = 0) = 0;
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
                                                               
  virtual Vec3D getVelocity(double *V) { return Vec3D(V[1], V[2], V[3]); }
  virtual double getDensity(double *V) { return V[0]; }
  virtual double getVelocityX(double *V) {return V[1];}
  virtual double getVelocityY(double *V) {return V[2];}
  virtual double getVelocityZ(double *V) {return V[3];}
  virtual double getPressure(double *V, double phi = 0.0) {return 0.0;}
  virtual double checkPressure(double *V, double phi = 0.0) {return 0.0;}
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
                                                                    
};
//------------------------------------------------------------------------------

class VarFcnPerfectGas : public VarFcn {

protected:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;  

public:
  VarFcnPerfectGas(IoData &);
  ~VarFcnPerfectGas() {}
                            
  double getGamma()  {return gam;}
  double getGamma1()  {return gam1;}
  double getPressureConstant() {return Pstiff;}

  double getPressure(double *V, double phi = 0.0) { return V[4]; }
  double checkPressure(double *V, double phi = 0.0) { return V[4]+Pstiff; }
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
    double opmachr = 1.0 + 0.5*gam1*machr2;
    return popr*pow(opmach/opmachr, gam*invgam1);
  }

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
  double checkPressure(double *V, double phi = 0.0) { return 1.0; /*alpha_water * pow(V[0], beta_water);*/ }
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
      if(gam * (V[4]+Pstiff) / V[0]<0.0) fprintf(stdout, "gam = %e, Pres = %e, pstiff = %e, rho = %e\n", gam,V[4],Pstiff,V[0]);
      return sqrt(gam * (V[4]+Pstiff) / V[0]); 
    }
    else{
      if(gamp * (V[4]+Pstiffp) / V[0]<0.0) fprintf(stdout, "gamp = %e, Pres = %e, pstiffp = %e, rho = %e\n", gamp,V[4],Pstiffp,V[0]);
      return sqrt(gamp * (V[4]+Pstiffp) / V[0]);
    }
  }
    
  double computeMachNumber(double *V, double phi = 0.0) { 
    if (phi>=0.0) return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)));
    else         return sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gamp * (V[4]+Pstiffp))); }
    
  double computeTotalPressure(double machr, double *V, double phi = 0.0){
    double mach, machr2, popr, opmach, opmachr;
    if (phi>=0.0) {
      mach = computeMachNumber(V, phi);
      machr2 = machr*machr;
      popr = V[4]*gam*machr2;
      opmach = 1.0 + 0.5*gam1*mach*mach;
      opmachr = 1.0 + 0.5*gam1*machr2;
      return popr*pow(opmach/opmachr, gam*invgam1);
    }else{
      mach = computeMachNumber(V, phi);
      machr2 = machr*machr;
      popr = V[4]*gamp*machr2;
      opmach = 1.0 + 0.5*gamp1*mach*mach;
      opmachr = 1.0 + 0.5*gamp1*machr2;
      return popr*pow(opmach/opmachr, gamp*invgamp1);
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
  Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;

  gamp = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
  gamp1 = gamp - 1.0;
  invgamp1 = 1.0/gamp1;
  Pstiffp = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;

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
      return 1.0; //Pref_water + alpha_water * pow(V[0], beta_water); 
    else
      return V[4]+Pstiff;
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
    if (phi>=0.0){
      if (alpha_water * beta_water * pow(V[0], beta_water - 1.0) < 0.0) fprintf(stdout, "c2_water = %e\n", alpha_water * beta_water * pow(V[0], beta_water - 1.0));
      return sqrt(alpha_water * beta_water * pow(V[0], beta_water - 1.0));
    }
    if (gam * (V[4]+Pstiff) / V[0] < 0.0) fprintf(stdout, "c2_air = %e - P = %e - rho = %e\n", gam * (V[4]+Pstiff) / V[0], V[4], V[0]);
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

};

//--------------------------------------------------------------------------------------

inline
VarFcnGasInLiquid::VarFcnGasInLiquid(IoData &iod)
: VarFcn(iod) {

  type = GASINLIQUID; 
  subType = NONE;

	if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
     iod.eqs.fluidModel2.fluid == FluidModelData::GAS){
    //gam  = 1.0  +iod.eqs.fluidModel2.gasModel.idealGasConstant/iod.eqs.fluidModel.liquidModel.Cv;
    gam = iod.eqs.fluidModel2.gasModel.specificHeatRatio;
    gam1 = gam -1.0;
    invgam1 = 1.0/gam1;
    Pstiff = iod.eqs.fluidModel2.gasModel.pressureConstant/iod.ref.rv.pressure;

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
    Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.ref.rv.pressure;

    Cv=1.0;
    invCv=1.0/Cv;

    alpha_water = iod.eqs.fluidModel2.liquidModel.alpha;
    beta_water  = iod.eqs.fluidModel2.liquidModel.beta;
    Pref_water  = iod.eqs.fluidModel2.liquidModel.Pref;

	}

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
/*
template<int dim>
void VarFcn::primitiveToConservative(SVec<double,dim> &U, SVec<double,dim> &V, Vec<double> *Phi, Vec<double> *Phi1, SVec<double,dim> *Vgf)
{

  if (!Phi1) {
    if (Phi){
      for (int i=0; i<U.size(); ++i)
        primitiveToConservative(U[i], V[i], (*Phi)[i]);
    }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(U[i], V[i]);
    }
  }
  else {
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(U[i], V[i], (*Phi)[i], &((*Phi1)[i]), (*Vgf)[i], i);
  }
}
*/
//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(SVec<double,dim> &U, SVec<double,dim> &V, Vec<double> *Phi)
{
  if (Phi){
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(U[i], V[i], (*Phi)[i]);
  }else{
    for (int i=0; i<U.size(); ++i)
      primitiveToConservative(U[i], V[i]);
  }
}

//------------------------------------------------------------------------------
/*
template<int dim>
void VarFcn::primitiveToConservative(DistSVec<double,dim> &U, DistSVec<double,dim> &V, DistVec<double> *Phi, DistVec<double> *Phi1, DistSVec<double,dim> *Vgf)
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
          primitiveToConservative(u[i], v[i], phi[i]);
      }else{
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(u[i], v[i]);
      }
    }
    else {
        if(node_change) fprintf(stdout, "node change in locsubD = %d\n", iSub);
        double *phi  = (*Phi).subData(iSub);
        double *phi1 = (*Phi1).subData(iSub);
        double (*vgf)[dim]  = (*Vgf).subData(iSub);
        for (int i=0; i<U.subSize(iSub); ++i)
          primitiveToConservative(u[i], v[i], phi[i], &(phi1[i]), vgf[i], i);
    }
  }

}
*/
//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(DistSVec<double,dim> &U, DistSVec<double,dim> &V, DistVec<double> *Phi)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if (Phi){
      double *phi = (*Phi).subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservative(u[i], v[i], phi[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        primitiveToConservative(u[i], v[i]);
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
    if(weight){
      w = (*weight).subData(iSub);
		  bool change = false;
      for ( int i=0; i<U.subSize(iSub); ++i){
			  change = updatePhaseChange(v[i],u[i],phi[i],phin[i],r[i],w[i]);
	  	}
    }else{
      bool change = false;
      for ( int i=0; i<U.subSize(iSub); ++i)
        change = updatePhaseChange(v[i],u[i],phi[i],phin[i],r[i]);
    }
	}


}

//------------------------------------------------------------------------------
template<int dim>
void VarFcn::computeNewPrimitive(SVec<double,dim> &U, SVec<double,dim> &V){
 
  for (int i=0; i<U.size(); ++i)
    computeNewPrimitive(U[i], V[i]);
}



#endif
