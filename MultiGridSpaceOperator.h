#pragma once

#include <MultiGridKernel.h>

#include <MultiGridMvpMatrix.h>

#include <MultiGridDistSVec.h>

template <class Scalar,int dim>
class MultiGridSpaceOperator {

 public:

  MultiGridSpaceOperator(IoData&,Domain*, SpaceOperator<dim>*,
                         MultiGridKernel<Scalar>* pKernel,
                         SpaceOperator<dim>* = NULL);

  ~MultiGridSpaceOperator();

  void setupBcs(DistBcData<dim>* bcd);

  void computeTimeStep(int level, double cfl, 
                       MultiGridDistSVec<Scalar,dim>& V);

  void computeResidual(int level, MultiGridDistSVec<Scalar,dim>& U,
                       MultiGridDistSVec<Scalar,dim>& V,
                       MultiGridDistSVec<Scalar,dim>& res,
                       bool addDWdt = true);

  void updateStateVectors(int lvl, MultiGridDistSVec<Scalar,dim>& U) ;

  template <int neq>
  void computeJacobian(int level, MultiGridDistSVec<Scalar,dim>& U,
                       MultiGridDistSVec<Scalar,dim>& V,
                       MultiGridMvpMatrix<Scalar,neq>& mvp);

 private:

  int nLevels;

  MultiGridKernel<Scalar>* pKernel;

  MultiGridOperator<Scalar,dim>** myOperators;

  RecFcnConstant<dim> recConstant;

  FemEquationTerm* fet,*fet1; 

  FluxFcn** fluxFcn, **fluxFcn1;

  VarFcn* varFcn;
};
