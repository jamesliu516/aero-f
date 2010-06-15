#ifndef _THERMAL_COND_FCN_H_
#define _THERMAL_COND_FCN_H_

#include <IoData.h>
#include <ViscoFcn.h>

//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	The thermal condition function is not necessarily the same for gas
//	 and liquid, thus we have added a new one.
//------------------------------------------------------------------------------

class ThermalCondFcn {

public:

  ThermalCondFcn() {}
  ~ThermalCondFcn() {}

  virtual double compute(double) = 0;

// Included (MB)
  virtual double computeDerivative(double, double, double) = 0;
  virtual void rstVar(IoData&) = 0;

};

//------------------------------------------------------------------------------

class ConstantPrandtlThermalCondFcn : public ThermalCondFcn {

  double alpha;

  ViscoFcn *viscoFcn;

public:

  ConstantPrandtlThermalCondFcn(IoData &iod, ViscoFcn *visf) 
  {
    viscoFcn = visf;
    alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.thermalCondModel.prandtl;
  }
  ~ConstantPrandtlThermalCondFcn() {}

  double compute(double Tadim) 
  { 
    return alpha * viscoFcn->compute_mu(Tadim);
  }

// Included (MB)
  double computeDerivative(double Tadim, double dTadim, double dMach)
  {
    return alpha * viscoFcn->compute_muDerivative(Tadim, dTadim, dMach);
  }
  void rstVar(IoData &iod)
  {
    viscoFcn->rstVar(iod);
  }

};

//------------------------------------------------------------------------------

/*class WaterThermalCondFcn : public ThermalCondFcn {
 
  double ...

public:

  WaterThermalCondFcn(IoData &iod, ViscoFcn *visf){
  }
  ~WaterThermalCondFcn() {}

  double compute(double Tadim) {}

};
*/

//---------------------------------------------------------------------------------

#endif
