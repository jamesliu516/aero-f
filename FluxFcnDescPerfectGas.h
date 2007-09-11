#ifndef _FLUX_FCN_DESC_PERFECT_GAS_H_
#define _FLUX_FCN_DESC_PERFECT_GAS_H_

#include <FluxFcnDesc.h>
#include <VarFcnDesc.h>

class IoData;
//------------------------------------------------------------------------------

class FluxFcnPerfectGasFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnPerfectGasFDJacRoeEuler3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(gg, br, K1, cm, pr, new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasFDJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnPerfectGasApprJacRoeEuler3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(rs, gg, br, K1, cm, pr, new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasApprJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnPerfectGasExactJacRoeEuler3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasExactJacRoeEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasVanLeerEuler3D : public FluxFcnVanLeerEuler3D {

public:

  FluxFcnPerfectGasVanLeerEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnVanLeerEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasVanLeerEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnPerfectGasWallEuler3D(IoData &ioData, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasWallEuler3D() {}
  
  void compute(double, double, double *, double, double *, double *, double *, int);

};
//------------------------------------------------------------------------------

class FluxFcnPerfectGasGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnPerfectGasGhidagliaEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasGhidagliaEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnPerfectGasInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnPerfectGasInternalInflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalInflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnPerfectGasOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnPerfectGasInternalOutflowEuler3D(IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(new VarFcnPerfectGasEuler3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalOutflowEuler3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//turbulence
//NOTE : for the coupled solver, the varFcn correspond to the FluxFcn in the sense 
//       that if varFcn = VarFcnPerfectGasSA3D, then fluxFcn = FluxFcnPerfectGasApprJacRoeSA3D
//       HOWEVER, for the segregated solver, they do not necessarily correspond since we consider
//       two different fluxes ff1 and ff2. Then there must be two varFcn (vf1 and vf2) corresponding
//       to each fluxFcn, and not one varFcn (corresponding to the physical case we are considering)
//       for both ff1 and ff2
//
//
//NOTE2:  some FluxFcn do not really need varFcn, but they were given one varFcn still (easier to implement)
//        for instance all fluxFcn of type SAturb3D and KEturb3D do not need a real varFcn. Moreover,
//        there is no corresponding varFcn for these fluxes because the varFcn always consider at least
//        the Euler variables, but never just the turbulent variables.




class FluxFcnPerfectGasFDJacRoeSA3D : public FluxFcnFDJacRoeSA3D {

 public:
  
  FluxFcnPerfectGasFDJacRoeSA3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeSA3D(gg, br, K1, cm, pr, new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasFDJacRoeSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasApprJacRoeSA3D : public FluxFcnApprJacRoeSA3D {

public:

  FluxFcnPerfectGasApprJacRoeSA3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeSA3D(rs,gg, br, K1, cm, pr, new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasApprJacRoeSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasExactJacRoeSA3D : public FluxFcnExactJacRoeSA3D {

public:

  FluxFcnPerfectGasExactJacRoeSA3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeSA3D(gg, new VarFcnPerfectGasSA3D(ioData), tp) { }
  ~FluxFcnPerfectGasExactJacRoeSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasWallSA3D : public FluxFcnWallSA3D {

public:

  FluxFcnPerfectGasWallSA3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnWallSA3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasWallSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasOutflowSA3D : public FluxFcnOutflowSA3D {

public:

  FluxFcnPerfectGasOutflowSA3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSA3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasOutflowSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalInflowSA3D : public FluxFcnInternalInflowSA3D {

public:

  FluxFcnPerfectGasInternalInflowSA3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSA3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalInflowSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalOutflowSA3D : public FluxFcnInternalOutflowSA3D {

public:

  FluxFcnPerfectGasInternalOutflowSA3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSA3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalOutflowSA3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasRoeSAturb3D : public FluxFcnRoeSAturb3D {

public:

  //FluxFcnPerfectGasRoeSAturb3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
  FluxFcnPerfectGasRoeSAturb3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnRoeSAturb3D(gg, br, K1, cm, pr, new VarFcnPerfectGasSA3D(ioData), tp) {}
    //FluxFcnRoeSAturb3D(gg, new VarFcnPerfectGasSA3D(ioData), tp) { }
  ~FluxFcnPerfectGasRoeSAturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasWallSAturb3D : public FluxFcnWallSAturb3D {

public:
  FluxFcnPerfectGasWallSAturb3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnWallSAturb3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasWallSAturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasOutflowSAturb3D : public FluxFcnOutflowSAturb3D {

public:
  FluxFcnPerfectGasOutflowSAturb3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSAturb3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasOutflowSAturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalInflowSAturb3D : public FluxFcnInternalInflowSAturb3D {

public:
  FluxFcnPerfectGasInternalInflowSAturb3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSAturb3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalInflowSAturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasInternalOutflowSAturb3D : public FluxFcnInternalOutflowSAturb3D {

public:
  FluxFcnPerfectGasInternalOutflowSAturb3D(IoData  &ioData, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSAturb3D(new VarFcnPerfectGasSA3D(ioData), tp) {}
  ~FluxFcnPerfectGasInternalOutflowSAturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasFDJacRoeKE3D : public FluxFcnFDJacRoeKE3D{

public:

  FluxFcnPerfectGasFDJacRoeKE3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeKE3D(gg, br, K1, cm, pr, new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasFDJacRoeKE3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasApprJacRoeKE3D : public FluxFcnApprJacRoeKE3D {

public:
  FluxFcnPerfectGasApprJacRoeKE3D(int rs, double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnApprJacRoeKE3D(rs, gg, br, K1, cm, pr, new VarFcnPerfectGasKE3D(ioData), tp) { }
  ~FluxFcnPerfectGasApprJacRoeKE3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasExactJacRoeKE3D : public FluxFcnExactJacRoeKE3D {

public:
  FluxFcnPerfectGasExactJacRoeKE3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeKE3D(gg, new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasExactJacRoeKE3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasWallKE3D : public FluxFcnWallKE3D {

public:

  FluxFcnPerfectGasWallKE3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnWallKE3D(new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasWallKE3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasOutflowKE3D : public FluxFcnOutflowKE3D {

public:

  FluxFcnPerfectGasOutflowKE3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKE3D(new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasOutflowKE3D() {}

  void compute(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasRoeKEturb3D : public FluxFcnRoeKEturb3D {

public:
  FluxFcnPerfectGasRoeKEturb3D(double gg, double br, double K1, double cm, int pr, IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnRoeKEturb3D(gg, br, K1, cm, pr, new VarFcnPerfectGasKE3D(ioData), tp) {}
  //FluxFcnPerfectGasRoeKEturb3D(double gg, IoData &ioData, Type tp = CONSERVATIVE) :
  //  FluxFcnRoeKEturb3D(gg, new VarFcnPerfectGasKE3D(ioData), tp) { }
  ~FluxFcnPerfectGasRoeKEturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasWallKEturb3D : public FluxFcnWallKEturb3D {

public:

  FluxFcnPerfectGasWallKEturb3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnWallKEturb3D(new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasWallKEturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnPerfectGasOutflowKEturb3D : public FluxFcnOutflowKEturb3D{

public:

  FluxFcnPerfectGasOutflowKEturb3D(IoData &ioData, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKEturb3D(new VarFcnPerfectGasKE3D(ioData), tp) {}
  ~FluxFcnPerfectGasOutflowKEturb3D() {}

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, int flag) {}
  void computeJacobian(double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

#endif

