/* MultiGridOperator.h

 */

#pragma once

#include <DistVector.h>
#include <SpaceOperator.h>
#include <DistMvpMatrix.h>
#include <DistTimeState.h>

template <class Scalar> class MultiGridLevel;

template<class Scalar,int dim>
class MultiGridOperator {

 public:

  MultiGridOperator(MultiGridLevel<Scalar>*, IoData& ioData, VarFcn* varFcn,
                    Domain* domain, BcFcn* = NULL);

  ~MultiGridOperator();
    
  template <class Scalar2,int neq>
  void computeJacobian(DistSVec<Scalar2,dim>& U, 
                       DistSVec<Scalar2,dim>& V,
//                       DistVec<Scalar2>& irey,
                       FluxFcn **fluxFcn, 
                       DistMvpMatrix<Scalar2,neq> &A); 

  template <class Scalar2>
  void computeResidual(DistSVec<Scalar2,dim>& V,
                       DistSVec<Scalar2,dim>& U,
  //                     DistVec<Scalar2>& irey,
                       FluxFcn** fluxFcn,
                       RecFcn* recFcn,
                       FemEquationTerm*, 
                       DistSVec<Scalar2,dim>& res,
                       bool addDWdt = true);

  template <class Scalar2>
  void applyBCsToResidual(DistSVec<Scalar2,dim>& U,
                          DistSVec<Scalar2,dim>& R);

  template <class Scalar2,int neq>
  void applyBCsToJacobian(DistSVec<Scalar2,dim>& U,
                          DistMvpMatrix<Scalar2,neq>& A);
  
  void computeGradientsLeastSquares(DistSVec<Scalar,dim>& U) ; 
 
  DistBcData<dim>& getBcData() const { return *myBcData; }
  DistTimeState<dim>& getTimeState() const { return *timeState; } 

  void computeTimeStep(double cfl, VarFcn *varFcn,
                       DistSVec<double,dim> &V);

  void updateStateVectors(DistSVec<Scalar,dim>&);

  DistSVec<Scalar,dim>& getBoundaryState() const { return *boundaryState; }

 private:

  MultiGridLevel<Scalar>* mgLevel;
    
  //SpaceOperator<dim>* mySpaceOperator;
 
  DistTimeState<dim>* timeState;
  DistBcData<dim>* myBcData;
  DistSVec<Scalar,dim>* boundaryState;

  DistSVec<Scalar,dim>* zero;
  DistSVec<Scalar,dim>* DX[3];
  DistVec<Scalar>* scalar_zero;
  DistVec<Scalar>* idti;
  DistVec<Scalar>* idtv;

  DistNodalGrad<dim,Scalar>* myNodalGrad;

  VarFcn* myVarFcn;

  BcFcn* bcFcn;

};
