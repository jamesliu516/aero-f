#ifndef _FLUX_FCN_DESC_WATER_COMPRESSIBLE_H_
#define _FLUX_FCN_DESC_WATER_COMPRESSIBLE_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleFDJacRoeEuler3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(gg, br, K1, cm, pr, new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleFDJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleApprJacRoeEuler3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(rs, gg, br, K1, cm, pr, new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleApprJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnWaterCompressibleExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleExactJacRoeEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

//NO VAN LEER 

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnWaterCompressibleWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleWallEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnWaterCompressibleGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleGhidagliaEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnWaterCompressibleInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------


class FluxFcnWaterCompressibleOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnWaterCompressibleOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnWaterCompressibleInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleInternalInflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWaterCompressibleInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnWaterCompressibleInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnWaterCompressibleEuler3D(ioData), tp) {}
  ~FluxFcnWaterCompressibleInternalOutflowEuler3D() {}

  void compute(double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

#endif
