#ifndef _FLUX_FCN_DESC_LIQUID_IN_LIQUID_H_
#define _FLUX_FCN_DESC_LIQUID_IN_LIQUID_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnLiquidInLiquidFDJacRoeEuler3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(gg, br, K1, cm, pr, new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidFDJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnLiquidInLiquidApprJacRoeEuler3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(rs, gg, br, K1, cm, pr, new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidApprJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnLiquidInLiquidExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidExactJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

//NO VAN LEER 

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnLiquidInLiquidWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidWallEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnLiquidInLiquidGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidGhidagliaEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnLiquidInLiquidInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnLiquidInLiquidOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnLiquidInLiquidInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidInternalInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnLiquidInLiquidInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnLiquidInLiquidInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnLiquidInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnLiquidInLiquidInternalOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

#endif
