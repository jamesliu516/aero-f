#ifndef _FLUX_FCN_DESC_TAIT_H_
#define _FLUX_FCN_DESC_TAIT_H_

#include <FluxFcnDesc.h>

class IoData;
#include "VarFcnTait.h"

//------------------------------------------------------------------------------

class FluxFcnTaitFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnTaitFDJacRoeEuler3D(double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, varFcnTait, tp) {}
  ~FluxFcnTaitFDJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnTaitApprJacRoeEuler3D(int rs, double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, varFcnTait, tp) {}
  ~FluxFcnTaitApprJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *);


};

//------------------------------------------------------------------------------

class FluxFcnTaitExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnTaitExactJacRoeEuler3D(double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, varFcnTait, tp) {}
  ~FluxFcnTaitExactJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *);


};

//------------------------------------------------------------------------------

class FluxFcnTaitWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnTaitWallEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnTaitGhidagliaEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnTaitInflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnTaitOutflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnTaitInternalInflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInternalInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);
  void computeJacobian(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnTaitInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnTaitInternalOutflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInternalOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *);
  void computeJacobian(double, double, double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

#endif
