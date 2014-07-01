#ifndef _FLUX_FCN_DESC_H_
#define _FLUX_FCN_DESC_H_

#include <FluxFcnBase.h>

#include <LowMachPrec.h>

class VarFcnBase;
class IoData;

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacRoeEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp):
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; } 

  ~FluxFcnFDJacRoeEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacRoeEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacRoeEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeEuler3D : public FluxFcnBase {

 protected:
  double gamma;

public:

  FluxFcnExactJacRoeEuler3D(double gg, VarFcnBase *vf, Type tp) : FluxFcnBase(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLEEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLEEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLEEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLEEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLCEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLCEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLCEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLCEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnVanLeerEuler3D : public FluxFcnBase {

public:

  FluxFcnVanLeerEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnVanLeerEuler3D() {}
  
protected:
  void evalFlux(double, double, double *, double, double *, double *, int);
  void evalJac(double, double, double *, double, double *, double *, int);  

// Included (MB)
  void evalDerivativeOfFlux(double, double, double, double *, double *, double, double, double *, double *, double, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWallEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnWallEuler3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnWallEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnGhidagliaEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnGhidagliaEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnModifiedGhidagliaEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnModifiedGhidagliaEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnModifiedGhidagliaEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowEuler3D : public FluxFcnBase {

public:

  FluxFcnInternalInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnPrescribedInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnPrescribedInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnPrescribedInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnOutflowEuler3D : public FluxFcnFD<5> {
  
public:

  FluxFcnOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnOutflowEuler3D() {} 
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowEuler3D : public FluxFcnBase {

public:

  FluxFcnInternalOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
//
class FluxFcnPrescribedOutflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnPrescribedOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnPrescribedOutflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
//turbulence

class FluxFcnFDJacRoeSA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;
  
 public:
  
  FluxFcnFDJacRoeSA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeSA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;
  
public:

  FluxFcnApprJacRoeSA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) : 
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeSA3D : public FluxFcnBase {

 protected:
  double gamma;

public:

  FluxFcnExactJacRoeSA3D(double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLESA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLESA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLESA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLESA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLESA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLESA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCSA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLCSA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLCSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCSA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLCSA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLCSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallSA3D : public FluxFcnFD<6> {

public:

  FluxFcnWallSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnWallSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaSA3D : public FluxFcnFD<6>{

public:

  FluxFcnGhidagliaSA3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<6>(vf,tp) {}

  ~FluxFcnGhidagliaSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSA3D : public FluxFcnBase {

public:

  FluxFcnInternalInflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowSA3D() {}
  
};

//------------------------------------------------------------------------------
//
class FluxFcnPrescribedInflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnPrescribedInflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnPrescribedInflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSA3D : public FluxFcnBase {

public:

  FluxFcnInternalOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnPrescribedOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnPrescribedOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnPrescribedOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnRoeSAturb3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnRoeSAturb3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnRoeSAturb3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaSAturb3D : public FluxFcnBase{

public:
  FluxFcnGhidagliaSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf,tp) {}

  ~FluxFcnGhidagliaSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallSAturb3D : public FluxFcnBase {

public:
  FluxFcnWallSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnWallSAturb3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnInternalInflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnPrescribedInflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnPrescribedInflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnPrescribedInflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnInternalOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnPrescribedOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnPrescribedOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnPrescribedOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacRoeKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacRoeKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeKE3D : public FluxFcnBase {

 protected:
  double gamma;

public:
  FluxFcnExactJacRoeKE3D(double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacHLLEKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLEKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacHLLEKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLEKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacHLLCKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLCKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacHLLCKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLCKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallKE3D : public FluxFcnFD<7> {

public:

  FluxFcnWallKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnWallKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaKE3D : public FluxFcnFD<7>{

public:

  FluxFcnGhidagliaKE3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<7>(vf,tp) {}

  ~FluxFcnGhidagliaKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnOutflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnOutflowKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnRoeKEturb3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnRoeKEturb3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnRoeKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallKEturb3D : public FluxFcnBase {

public:

  FluxFcnWallKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnWallKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaKEturb3D : public FluxFcnBase{

public:
  FluxFcnGhidagliaKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf,tp) {}

  ~FluxFcnGhidagliaKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowKEturb3D : public FluxFcnBase {

public:

  FluxFcnOutflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnOutflowKEturb3D() {}

};

//------------------------------------------------------------------------------

#endif
