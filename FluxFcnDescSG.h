#ifndef _FLUX_FCN_DESC_SG_H_
#define _FLUX_FCN_DESC_SG_H_

#include <FluxFcnDesc.h>

class IoData;
#include "VarFcnBase.h"
#include "VarFcnSGEuler.h"
#include "VarFcnSGSA.h"
#include "VarFcnSGKE.h"

//------------------------------------------------------------------------------

class FluxFcnSGFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnSGFDJacRoeEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnSGApprJacRoeEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnSGExactJacRoeEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGExactJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLEEuler3D : public FluxFcnFDJacHLLEEuler3D {

public:

  FluxFcnSGFDJacHLLEEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLEEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacHLLEEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLEEuler3D : public FluxFcnApprJacHLLEEuler3D {

public:

  FluxFcnSGApprJacHLLEEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLEEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacHLLEEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLCEuler3D : public FluxFcnFDJacHLLCEuler3D {

public:

  FluxFcnSGFDJacHLLCEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLCEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacHLLCEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLCEuler3D : public FluxFcnApprJacHLLCEuler3D {

public:

  FluxFcnSGApprJacHLLCEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLCEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacHLLCEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGVanLeerEuler3D : public FluxFcnVanLeerEuler3D {

public:

  FluxFcnSGVanLeerEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnVanLeerEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGVanLeerEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnSGWallEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGWallEuler3D() { vf = 0; }
  
  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB*)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};
//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnSGGhidagliaEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnSGInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnSGInternalInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInternalInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnSGOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnSGInternalOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInternalOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//turbulence
//NOTE : for the coupled solver, the varFcn correspond to the FluxFcn in the sense 
//       that if varFcn = VarFcnSGSA3D, then fluxFcn = FluxFcnSGApprJacRoeSA3D
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




class FluxFcnSGFDJacRoeSA3D : public FluxFcnFDJacRoeSA3D {

 public:
  
  FluxFcnSGFDJacRoeSA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeSA3D(ioData, gg, varFcnSGSA) {}
  ~FluxFcnSGFDJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeSA3D : public FluxFcnApprJacRoeSA3D {

public:

  FluxFcnSGApprJacRoeSA3D(int rs, double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeSA3D(ioData, rs,gg, varFcnSGSA, tp) {}
  ~FluxFcnSGApprJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeSA3D : public FluxFcnExactJacRoeSA3D {

public:

  FluxFcnSGExactJacRoeSA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeSA3D(gg, varFcnSGSA, tp) {}
  ~FluxFcnSGExactJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLESA3D : public FluxFcnFDJacHLLESA3D {

 public:

  FluxFcnSGFDJacHLLESA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnFDJacHLLESA3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGFDJacHLLESA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLESA3D : public FluxFcnApprJacHLLESA3D {

public:

  FluxFcnSGApprJacHLLESA3D(int rs, double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnApprJacHLLESA3D(ioData, rs,gg, varFcnSGSA, tp) {}
  ~FluxFcnSGApprJacHLLESA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallSA3D : public FluxFcnWallSA3D {

public:

  FluxFcnSGWallSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnWallSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGWallSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB*)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowSA3D : public FluxFcnOutflowSA3D {

public:

  FluxFcnSGOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowSA3D : public FluxFcnInternalInflowSA3D {

public:

  FluxFcnSGInternalInflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalInflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowSA3D : public FluxFcnInternalOutflowSA3D {

public:

  FluxFcnSGInternalOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGRoeSAturb3D : public FluxFcnRoeSAturb3D {

public:

  FluxFcnSGRoeSAturb3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnRoeSAturb3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGRoeSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallSAturb3D : public FluxFcnWallSAturb3D {

public:
  FluxFcnSGWallSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnWallSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGWallSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowSAturb3D : public FluxFcnOutflowSAturb3D {

public:
  FluxFcnSGOutflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowSAturb3D : public FluxFcnInternalInflowSAturb3D {

public:
  FluxFcnSGInternalInflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalInflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowSAturb3D : public FluxFcnInternalOutflowSAturb3D {

public:
  FluxFcnSGInternalOutflowSAturb3D(IoData  &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacRoeKE3D : public FluxFcnFDJacRoeKE3D{

public:

  FluxFcnSGFDJacRoeKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGFDJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeKE3D : public FluxFcnApprJacRoeKE3D {

public:
  FluxFcnSGApprJacRoeKE3D(int rs, double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnApprJacRoeKE3D(ioData, rs, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGApprJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLEKE3D : public FluxFcnFDJacHLLEKE3D{

public:

  FluxFcnSGFDJacHLLEKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLEKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGFDJacHLLEKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLEKE3D : public FluxFcnApprJacHLLEKE3D {

public:
  FluxFcnSGApprJacHLLEKE3D(int rs, double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLEKE3D(ioData, rs, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGApprJacHLLEKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeKE3D : public FluxFcnExactJacRoeKE3D {

public:
  FluxFcnSGExactJacRoeKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeKE3D(gg, varFcnSGKE, tp) {}
  ~FluxFcnSGExactJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallKE3D : public FluxFcnWallKE3D {

public:

  FluxFcnSGWallKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnWallKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGWallKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB*)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowKE3D : public FluxFcnOutflowKE3D {

public:

  FluxFcnSGOutflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGOutflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

// Included (MB)
  void computeDerivative(double, double, double *, double *, double, double, double *, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGRoeKEturb3D : public FluxFcnRoeKEturb3D {

public:
  FluxFcnSGRoeKEturb3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnRoeKEturb3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGRoeKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallKEturb3D : public FluxFcnWallKEturb3D {

public:

  FluxFcnSGWallKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnWallKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGWallKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowKEturb3D : public FluxFcnOutflowKEturb3D{

public:

  FluxFcnSGOutflowKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGOutflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

#endif

