#ifndef _DIST_BC_DATA_H_
#define _DIST_BC_DATA_H_

#include <DistVector.h>
#include <VarFcn.h>

class IoData;
class SubDomain;
class Communicator;
class Domain;

template<int dim> class BcData;

//------------------------------------------------------------------------------
/* 
   angles[0] in the x-z plane (flow parallel to x-axis if alpha[0]=0)
   angles[1] in the x-y plane (flow parallel to x-axis if alpha[1]=0)
*/

template<int dim>
class DistBcData {

protected:
  enum BoundaryFluid { GAS=0, TAIT=1 } boundaryFluid;

  double angles[2];

  bool gravityOn;
  double gravity;
  double depth;
  double ngravity[3];

  VarFcn *vf;

  DistSVec<double,3> Xdot;
  DistVec<double> Temp;

  DistSVec<double,dim> Ufarin;
  DistSVec<double,dim> Ufarout;
  double Uin[dim], Vin[dim];
  double Uout[dim], Vout[dim];
  double Ub[dim];

  DistSVec<double,dim> Uface;
  DistSVec<double,dim> Unode;
  DistSVec<double,dim> Uinletnode;

  int numLocSub;

  Communicator *com;
  SubDomain **subDomain;
  BcData<dim> **subBcData;
  
  map<int, RotationData *> &rotInfo;
  double tref;
  double vref;

protected:

  void finalize(VarFcn *, DistSVec<double,3> &);

public:

  DistBcData(IoData &, VarFcn *, Domain *);
  ~DistBcData();

  BcData<dim> &operator() (int i) const { return *subBcData[i]; }

  void update(DistSVec<double,3> &);
  virtual void updateFarField(DistSVec<double,3> &) {}
  virtual void computeNodeValue(DistSVec<double,3> &) {}

  DistSVec<double,3> &getVelocityVector() { return Xdot; }
  DistVec<double> &getTemperatureVector() { return Temp; }
  DistSVec<double,dim> &getInletBoundaryVector()  { return Ufarin;  }
  DistSVec<double,dim> &getOutletBoundaryVector() { return Ufarout; }

  double *getInletAngles() { return angles; }
  double *getInletConservativeState() { return Uin; }
  double *getOutletConservativeState() { return Uout; }
  double *getInterface() { return Ub; }
  double *getInletPrimitiveState() { return Vin; }

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataEuler : public DistBcData<dim> {

private:
  void setBoundaryConditionsGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsLiquid(IoData &, VarFcn *, DistSVec<double,3> &);
  void setBoundaryConditionsGasGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsLiquidLiquid(IoData &, VarFcn *, DistSVec<double,3> &);
  void setBoundaryConditionsGasLiquid(IoData &, VarFcn *, DistSVec<double,3> &);
  void setBoundaryConditionsLiquidGas(IoData &, VarFcn *, DistSVec<double,3> &);

  void updateFarField(DistSVec<double,3> &);
  void updateFarFieldGas(DistSVec<double,3> &);
  void updateFarFieldLiquid(DistSVec<double,3> &);

public:

  DistBcDataEuler(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataEuler() {}  

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataSA : public DistBcDataEuler<dim> {

  DistSVec<double,2> *tmp;
  CommPattern<double> *vec2Pat;

public:

  DistBcDataSA(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataSA();  

  void computeNodeValue(DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataKE : public DistBcDataEuler<dim> {

  DistSVec<double,3> *tmp;
  CommPattern<double> *vec3Pat;

public:

  DistBcDataKE(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataKE();

  void computeNodeValue(DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistBcData.C>
#endif

#endif
