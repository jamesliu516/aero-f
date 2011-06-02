#ifndef _VISCO_FCN_H_
#define _VISCO_FCN_H_

#include <IoData.h>

#include <cmath>


const double twothird = 2.0/3.0;
//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	new behavior of Lame coefficients (lambda and mu)
//------------------------------------------------------------------------------

class ViscoFcn {

protected:
  double ooreynolds_mu;
  double reynolds_lambda;
// Included (MB)
  double dReMach;
  double dRe_muMach;
  double dRe_lambdaMach;

public:

  ViscoFcn(IoData &iod) { ooreynolds_mu = 1.0/iod.ref.reynolds_mu; reynolds_lambda = iod.ref.reynolds_lambda;
// Included (MB)
  dRe_muMach = iod.ref.dRe_mudMach;
  dRe_lambdaMach = iod.ref.dRe_lambdadMach;
  }
  ~ViscoFcn() {}

  virtual double compute_mu(double) = 0;
  virtual double compute_lambda(double, double) = 0;

// Included (MB)
  virtual double compute_muDerivative(double, double, double) = 0;
  virtual double compute_lambdaDerivative(double, double, double) = 0;
  virtual void rstVar(IoData &) = 0;

};

//------------------------------------------------------------------------------

class ConstantViscoFcn : public ViscoFcn {

public:

  ConstantViscoFcn(IoData &iod):ViscoFcn(iod) {}
  ~ConstantViscoFcn() {}

  double compute_mu(double T) { return 1.0; }
  double compute_lambda(double T, double mu) { return -twothird*mu*reynolds_lambda; }
//  double compute_lambda(double T, double mu) { return mu;}

// Included (MB)
  double compute_muDerivative(double T, double dT, double dMach) { return 0.0; }
  double compute_lambdaDerivative(double mu, double dmu, double dMach) { return -twothird*dmu*reynolds_lambda-twothird*mu*dRe_lambdaMach*dMach; }

  void rstVar(IoData &iod) {  
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    reynolds_lambda = iod.ref.reynolds_lambda;
    dRe_muMach = iod.ref.dRe_mudMach;
    dRe_lambdaMach = iod.ref.dRe_lambdadMach;
  }

};

//------------------------------------------------------------------------------

class SutherlandViscoFcn : public ViscoFcn {

  double alpha;
  double Ts;

// Included (MB)
  double dalpha;

public:

  SutherlandViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach; 
    Ts = iod.eqs.viscosityModel.sutherlandReferenceTemperature / iod.ref.temperature;

// Included (MB)
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }
  ~SutherlandViscoFcn() {}

  double compute_mu(double Tadim)
  {
    double T = alpha * Tadim;
    return T * sqrt(T) * (1.0 + Ts) / (T + Ts); 
  }
  
  double compute_lambda(double Tadim, double mu) { return -twothird*mu*reynolds_lambda; }
//  double compute_lambda(double T, double mu) { return mu;}

// Included (MB)
  double compute_muDerivative(double Tadim, double dTadim, double dMach)
  {
    double dAlpha = dalpha*dMach;
    double T = alpha * Tadim;
    double dT = dAlpha * Tadim + alpha * dTadim;
    return ( ( 1.5*(1.0 + Ts)*sqrt(T)*dT*(T + Ts) - T*sqrt(T)*(1.0 + Ts)*dT ) / ( (T + Ts)*(T + Ts) ) );
  }
  double compute_lambdaDerivative(double mu, double dmu, double dMach) { return -twothird*dmu*reynolds_lambda-twothird*mu*dRe_lambdaMach*dMach; }

  void rstVar(IoData &iod)
  {
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    reynolds_lambda = iod.ref.reynolds_lambda;
    dRe_muMach = iod.ref.dRe_mudMach;
    dRe_lambdaMach = iod.ref.dRe_lambdadMach;

    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach; 
    Ts = iod.eqs.viscosityModel.sutherlandReferenceTemperature / iod.ref.temperature;
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }

};

//------------------------------------------------------------------------------

class PrandtlViscoFcn : public ViscoFcn {

  double alpha;

// Included (MB)
  double dalpha;

public:

  PrandtlViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach;

// Included (MB)
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }
  ~PrandtlViscoFcn() {}

  double compute_mu(double Tadim) { return alpha * Tadim; }
  double compute_lambda(double Tadim, double mu) { return -twothird*mu*reynolds_lambda;}
//  double compute_lambda(double T, double mu) { return mu;}

// Included (MB)
  double compute_muDerivative(double Tadim, double dTadim, double dMach)
  {
    double dAlpha = dalpha*dMach;
    return (dAlpha * Tadim + alpha * dTadim);
  }
  double compute_lambdaDerivative(double mu, double dmu, double dMach) { return -twothird*dmu*reynolds_lambda-twothird*mu*dRe_lambdaMach*dMach; }

  void rstVar(IoData &iod)
  {
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    reynolds_lambda = iod.ref.reynolds_lambda;
    dRe_muMach = iod.ref.dRe_mudMach;
    dRe_lambdaMach = iod.ref.dRe_lambdadMach;

    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach;
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }

};

//------------------------------------------------------------------------------

/*class WaterViscoFcn : public ViscoFcn {

  double ...

public:

  WaterViscoFcn(IoData &iod) : ViscoFcn(iod){
  }
  ~WaterViscoFcn() {}

  double compute_mu() {}
  double compute_lambda() {}


};
*/
//------------------------------------------------------------------------------

#endif
