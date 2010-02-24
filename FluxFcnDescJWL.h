#ifndef _FLUX_FCN_DESC_JWL_H_
#define _FLUX_FCN_DESC_JWL_H_

#include <FluxFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnJWLFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnJWLFDJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, new VarFcn(ioData), tp) {}
  ~FluxFcnJWLFDJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, int = 1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnJWLApprJacRoeEuler3D(int rs, double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, new VarFcn(ioData), tp) {}
  ~FluxFcnJWLApprJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, int = 1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnJWLExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcn(ioData), tp) {}
  ~FluxFcnJWLExactJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, int = 1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnJWLWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLWallEuler3D() {}
  
  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB*)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int = 1) {}

};
//------------------------------------------------------------------------------

class FluxFcnJWLGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnJWLGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLGhidagliaEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnJWLInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int = 1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnJWLInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLInternalInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int = 1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnJWLOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int = 1) {}

};

//------------------------------------------------------------------------------

class FluxFcnJWLInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnJWLInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnJWLInternalOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int = 1) {}

};

//------------------------------------------------------------------------------

#endif

