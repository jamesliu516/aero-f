#ifndef _FLUX_FCN_DESC_WATER_COMPRESSIBLE_H_
#define _FLUX_FCN_DESC_WATER_COMPRESSIBLE_H_

#include <FluxFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleFDJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleFDJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleApprJacRoeEuler3D(int rs, double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleApprJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleExactJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

//NO VAN LEER 

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnWaterCompressibleWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleWallEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnWaterCompressibleGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleGhidagliaEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnWaterCompressibleInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------


class FluxFcnWaterCompressibleOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnWaterCompressibleOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnWaterCompressibleInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleInternalInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnWaterCompressibleInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcn(ioData), tp) {}
  ~FluxFcnWaterCompressibleInternalOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df, int fl=1) {}
  void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df, int fl=1) {}

};

//------------------------------------------------------------------------------

#endif
