#ifndef _FLUX_FCN_DESC_JWL_IN_GAS_H_
#define _FLUX_FCN_DESC_JWL_IN_GAS_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnJWLInGasFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnJWLInGasFDJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasFDJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnJWLInGasApprJacRoeEuler3D(int rs, double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasApprJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnJWLInGasExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasExactJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasVanLeerEuler3D : public FluxFcnVanLeerEuler3D {

public:

  FluxFcnJWLInGasVanLeerEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnVanLeerEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasVanLeerEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnJWLInGasWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasWallEuler3D() {}
  
  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnJWLInGasGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasGhidagliaEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnJWLInGasInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnJWLInGasInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasInternalInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnJWLInGasOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInGasInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnJWLInGasInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnJWLInGasEuler3D(ioData), tp) {}
  ~FluxFcnJWLInGasInternalOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};



#endif
