#ifndef _FLUX_FCN_DESC_GAS_IN_GAS_H_
#define _FLUX_FCN_DESC_GAS_IN_GAS_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnGasInGasFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnGasInGasFDJacRoeEuler3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(gg, br, K1, cm, pr, new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasFDJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnGasInGasApprJacRoeEuler3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(rs, gg, br, K1, cm, pr, new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasApprJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnGasInGasExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasExactJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasVanLeerEuler3D : public FluxFcnVanLeerEuler3D {

public:

  FluxFcnGasInGasVanLeerEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnVanLeerEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasVanLeerEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnGasInGasWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasWallEuler3D() {}
  
  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnGasInGasGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasGhidagliaEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnGasInGasInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnGasInGasInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasInternalInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnGasInGasOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnGasInGasInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnGasInGasInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnGasInGasEuler3D(ioData), tp) {}
  ~FluxFcnGasInGasInternalOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};



#endif
