#ifndef _VISCO_FCN_H_
#define _VISCO_FCN_H_

#include <IoData.h>

#include <math.h>


const double twothird = 2.0/3.0;
//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	new behavior of Lame coefficients (lambda and mu)
//------------------------------------------------------------------------------

class ViscoFcn {

protected:
  double ooreynolds_mu;
  double reynolds_lambda;

public:

  ViscoFcn(IoData &iod) { ooreynolds_mu = 1.0/iod.ref.reynolds_mu; reynolds_lambda = iod.ref.reynolds_lambda;}
  ~ViscoFcn() {}

  virtual double compute_mu(double) = 0;
  virtual double compute_lambda(double, double) = 0;
};

//------------------------------------------------------------------------------

class ConstantViscoFcn : public ViscoFcn {

public:

  ConstantViscoFcn(IoData &iod):ViscoFcn(iod) {}
  ~ConstantViscoFcn() {}

  double compute_mu(double T) { return 1.0; }
  double compute_lambda(double T, double mu) { return -twothird*mu*reynolds_lambda; }
//  double compute_lambda(double T, double mu) { return mu;}
};

//------------------------------------------------------------------------------

class SutherlandViscoFcn : public ViscoFcn {

  double alpha;
  double Ts;

public:

  SutherlandViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach; 
    Ts = iod.eqs.viscosityModel.sutherlandReferenceTemperature / iod.ref.temperature;
  }
  ~SutherlandViscoFcn() {}

  double compute_mu(double Tadim)
  {
    double T = alpha * Tadim;
    return T * sqrt(T) * (1.0 + Ts) / (T + Ts); 
  }
  
  double compute_lambda(double Tadim, double mu) { return -twothird*mu*reynolds_lambda; }
//  double compute_lambda(double T, double mu) { return mu;}

};

//------------------------------------------------------------------------------

class PrandtlViscoFcn : public ViscoFcn {

  double alpha;

public:

  PrandtlViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach;
  }
  ~PrandtlViscoFcn() {}

  double compute_mu(double Tadim) { return alpha * Tadim; }
  double compute_lambda(double Tadim, double mu) { return -twothird*mu*reynolds_lambda;}
//  double compute_lambda(double T, double mu) { return mu;}
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
