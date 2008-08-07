#ifndef _FLUX_FCN_DESC_GAS_IN_LIQUID_H_
#define _FLUX_FCN_DESC_GAS_IN_LIQUID_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnGasInLiquidFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnGasInLiquidFDJacRoeEuler3D(double gg, double br, double K1, double cm, double sr, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(gg, br, K1, cm, sr, pr, new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidFDJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnGasInLiquidApprJacRoeEuler3D(int rs, double gg, double br, double K1, double cm, double sr, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(rs, gg, br, K1, cm, sr, pr, new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidApprJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnGasInLiquidExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidExactJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

// NO VAN-LEER

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnGasInLiquidWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidWallEuler3D() {}
  
  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnGasInLiquidGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidGhidagliaEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnGasInLiquidInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnGasInLiquidInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidInternalInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnGasInLiquidOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnGasInLiquidInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnGasInLiquidInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnGasInLiquidEuler3D(ioData), tp) {}
  ~FluxFcnGasInLiquidInternalOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};



#endif
