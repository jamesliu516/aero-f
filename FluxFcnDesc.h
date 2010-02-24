#ifndef _FLUX_FCN_DESC_H_
#define _FLUX_FCN_DESC_H_

#include <FluxFcn.h>
#include <IoData.h>
#include <LowMachPrec.h>

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacRoeEuler3D(IoData &ioData, double gg, VarFcn *vf, Type tp) : 
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; } 

  ~FluxFcnFDJacRoeEuler3D() {}

 protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeEuler3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacRoeEuler3D(IoData &ioData, int rs, double gg, VarFcn *vf, Type tp) :
    FluxFcn(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacRoeEuler3D() {}
  
protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double, double *, double, double *, double *, double *);
  void computeJWL(double, double, double, double, double, double, double, double *, double, double *, double *, double *);

  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);
  void computeJacobiansBarotropicLiquid(double, double, double, double, double, double, double *, double, double *, double *, double *, double *, int);
  void computeJacobiansJWL(double, double, double, double, double, double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeEuler3D : public FluxFcn {

 protected:
  double gamma;

public:

  FluxFcnExactJacRoeEuler3D(double gg, VarFcn *vf, Type tp) : FluxFcn(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeEuler3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double *, double, double *, double *, double *, double *);
  void computeJacobiansBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLEEuler3D(IoData &ioData, double gg, VarFcn *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLEEuler3D() {}

 protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEEuler3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLEEuler3D(IoData &ioData, int rs, double gg, VarFcn *vf, Type tp) :
    FluxFcn(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLEEuler3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);
// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, 
                                     double *, double *, double *, double, double *, double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLCEuler3D(IoData &ioData, double gg, VarFcn *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLCEuler3D() {}

 protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCEuler3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLCEuler3D(IoData &ioData, int rs, double gg, VarFcn *vf, Type tp) :
    FluxFcn(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLCEuler3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, 
double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnVanLeerEuler3D : public FluxFcn {

public:

  FluxFcnVanLeerEuler3D(VarFcn *vf, Type tp) :
    FluxFcn(vf, tp) {}
  ~FluxFcnVanLeerEuler3D() {}
  
protected:
  void evalFlux(double, double, double *, double, double *, double *, int);
  void evalJac(double, double, double *, double, double *, double *, int);  

  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double *, double, double *, double *, double *, double *);

// Included (MB)
  void evalDerivativeOfFlux(double, double, double, double *, double *, double, double, double *, double *, double, double *, int);
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnWallEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnWallEuler3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnWallEuler3D() {}

protected:
  void computePerfectGas(double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double *, double, double *, double *, double *);
  void computeJWL(double *, double, double *, double *, double *);

// Included (MB*)
  void computeDerivativeOfPerfectGas(double *, double *, double, double, double *, double *, double *, double *, double *);
  void computeJacobianPerfectGas(double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnGhidagliaEuler3D(VarFcn *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnGhidagliaEuler3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);
  void computeJWL(double, double, double, double, double, double *, double, double *, double *, double *);
};

//------------------------------------------------------------------------------

class FluxFcnInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnInflowEuler3D(VarFcn *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnInflowEuler3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowEuler3D : public FluxFcn {

public:

  FluxFcnInternalInflowEuler3D(VarFcn *vf, Type tp) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalInflowEuler3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobianBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnOutflowEuler3D : public FluxFcnFD<5> {
  
public:

  FluxFcnOutflowEuler3D(VarFcn *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnOutflowEuler3D() {} 
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowEuler3D : public FluxFcn {

public:

  FluxFcnInternalOutflowEuler3D(VarFcn *vf, Type tp) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalOutflowEuler3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobianBarotropicLiquid(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------
//turbulence

class FluxFcnFDJacRoeSA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;
  
 public:
  
  FluxFcnFDJacRoeSA3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeSA3D() {}
  
protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeSA3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;
  
public:

  FluxFcnApprJacRoeSA3D(IoData &ioData, int rs, double gg, VarFcn* vf, Type tp = CONSERVATIVE) : 
    FluxFcn(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeSA3D() {}
  
protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeSA3D : public FluxFcn {

 protected:
  double gamma;

public:

  FluxFcnExactJacRoeSA3D(double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeSA3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double *, double, double *, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLESA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLESA3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLESA3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double,
 double *, double *){}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLESA3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLESA3D(IoData &ioData, int rs, double gg, VarFcn* vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLESA3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *){}

};

//------------------------------------------------------------------------------

class FluxFcnWallSA3D : public FluxFcnFD<6> {

public:

  FluxFcnWallSA3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnWallSA3D() {}
  
protected:
  void computePerfectGas(double *, double, double *, double *, double *);
  
// Included (MB*)
  void computeDerivativeOfPerfectGas(double *, double *, double, double, double *, double *, double *, double *, double *);
  void computeJacobianPerfectGas(double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnOutflowSA3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnOutflowSA3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSA3D : public FluxFcn {

public:

  FluxFcnInternalInflowSA3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalInflowSA3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSA3D : public FluxFcn {

public:

  FluxFcnInternalOutflowSA3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalOutflowSA3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnRoeSAturb3D : public FluxFcn {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnRoeSAturb3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnRoeSAturb3D() {}
  
protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnWallSAturb3D : public FluxFcn {

public:
  FluxFcnWallSAturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnWallSAturb3D() {}
  
protected:
  void computeJacobianPerfectGas(double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnOutflowSAturb3D : public FluxFcn {

public:
  FluxFcnOutflowSAturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnOutflowSAturb3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *);
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSAturb3D : public FluxFcn {

public:
  FluxFcnInternalInflowSAturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalInflowSAturb3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSAturb3D : public FluxFcn {

public:
  FluxFcnInternalOutflowSAturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnInternalOutflowSAturb3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *, int);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacRoeKE3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeKE3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeKE3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacRoeKE3D(IoData &ioData, int rs, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeKE3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);
  
// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeKE3D : public FluxFcn {

 protected:
  double gamma;

public:
  FluxFcnExactJacRoeKE3D(double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { gamma = gg; }
  ~FluxFcnExactJacRoeKE3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double *, double, double *, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double, double *, double *);

  
};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacHLLEKE3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLEKE3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double
*, double, double *, double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEKE3D : public FluxFcn {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacHLLEKE3D(IoData &ioData, int rs, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLEKE3D() {}

protected:
  void computePerfectGas(double, double, double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *, int);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double, double, double *, double *, double, double, double *, double *, double *, double
*, double, double *, double *) {}

};

//------------------------------------------------------------------------------

class FluxFcnWallKE3D : public FluxFcnFD<7> {

public:

  FluxFcnWallKE3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnWallKE3D() {}
  
protected:
  void computePerfectGas(double *, double, double *, double *, double *);

// Included (MB*)
  void computeDerivativeOfPerfectGas(double *, double *, double, double, double *, double *, double *, double *, double *);
  void computeJacobianPerfectGas(double *, double, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnOutflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnOutflowKE3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnOutflowKE3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);

// Included (MB)
  void computeDerivativeOfPerfectGas(double, double, double, double *, double *, double, double, double *, double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

class FluxFcnRoeKEturb3D : public FluxFcn {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnRoeKEturb3D(IoData &ioData, double gg, VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnRoeKEturb3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobiansPerfectGas(double, double, double, double, double *, double, double *, double *, double *, double *);
  
};

//------------------------------------------------------------------------------

class FluxFcnWallKEturb3D : public FluxFcn {

public:

  FluxFcnWallKEturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnWallKEturb3D() {}

protected:
  void computeJacobianPerfectGas(double *, double, double *, double *, double *);
  
};

//------------------------------------------------------------------------------

class FluxFcnOutflowKEturb3D : public FluxFcn {

public:

  FluxFcnOutflowKEturb3D(VarFcn *vf, Type tp = CONSERVATIVE) :
    FluxFcn(vf, tp) {}
  ~FluxFcnOutflowKEturb3D() {}

protected:
  void computePerfectGas(double, double, double *, double, double *, double *, double *);
  void computeJacobianPerfectGas(double, double, double *, double, double *, double *, double *);
  
};

//------------------------------------------------------------------------------
#endif
