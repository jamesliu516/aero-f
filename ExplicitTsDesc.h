#ifndef _EXPLICIT_TS_DESC_H_
#define _EXPLICIT_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>

class GeoSource;
class Domain;

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//-------------------------------------------------------------------------

template<int dim>
class ExplicitTsDesc : public TsDesc<dim> {

protected:
  LiquidModelData::Check check;
  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;
  DistSVec<double,dim> Ubc;
  bool RK4;

public:

  ExplicitTsDesc(IoData&, GeoSource&, Domain*);
  ~ExplicitTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim>&);
  void computeRKFourthOrder(DistSVec<double,dim>&);
  void computeRKSecondOrder(DistSVec<double,dim>&);

  void computeRKUpdate(DistSVec<double,dim>&, DistSVec<double,dim>&);

  virtual void computeFunctionLS(int, DistVec<double> &, DistVec<double> &,  DistVec<double> &,
                         DistSVec<double,dim> &,DistVec<double> &) { };
                                                                                                                                                         
  //-- overrides the functions implemented in ExplicitCoupledTsDesc.
  virtual void computeJacobianLS(int, DistVec<double> &, DistVec<double> &,
                         DistVec<double> &, DistSVec<double,dim> &,DistVec<double> &) { };

  virtual int solveLinearSystemLS(int, DistVec<double> &, DistVec<double> &) { };

  double reinitLS(DistVec<double> &Phi, DistSVec<double,dim> &U, int iti);
};

//-------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitTsDesc.C>
#endif

#endif
