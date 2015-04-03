#ifndef _EXACT_RIEMANN_SOLVER_H
#define _EXACT_RIEMANN_SOLVER_H

class IoData;
class LocalRiemann;
class VarFcn;
class SparseGridCluster;

template<class Scalar, int dim> class SVec;

#include "HigherOrderMultiFluid.h"


//------------------------------------------------------------------------------
template<int dim>
class ExactRiemannSolver{

  int numLriemann;
  LocalRiemann **lriemann;

  LocalRiemann *fsiRiemann;
	
  int iteration;
  SVec<double,dim>  &rupdate;
  Vec<double>       &weight;
  Vec<int>       &fluidIdToSet;
  SVec<double,dim-2>  &interfacialWi;
  SVec<double,dim-2>  &interfacialWj;

  int getRiemannSolverId(int i,int j) const;

  int levelSetMap[10][10];
  double levelSetSign[10][10];

  bool isHigherOrder;

  public:

  ExactRiemannSolver(IoData &, SVec<double,dim> &, Vec<double> &, 
                     SVec<double,dim-2> &, SVec<double,dim-2> &,
                     VarFcn *, SparseGridCluster *, Vec<int>& fluidIdToSet);
  ~ExactRiemannSolver();


  SVec<double,dim> &getRiemannUpdate() const { return rupdate; }
  Vec<double> &getRiemannWeight() const { return weight; }

  void storePreviousPrimitive(SVec<double,dim> &V, Vec<int> &fluidId, SVec<double,3> &X);
  int updatePhaseChange(SVec<double,dim> &V, Vec<int> &fluidId, Vec<int> &fluidIdn,
			 HigherOrderMultiFluid* );

  // for multiphase Riemann problem
  int computeRiemannSolution(double *Vi, double *Vj,
			     int IDi, int IDj, double *nphi, VarFcn *vf,
			     double *Wi, double *Wj,
			     int i, int j, int edgeNum, double dx[3],
			     int lsdim,
			     bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi, VarFcn *vf,
			      double *Wi, double *Wj,
			      int i, int j, int edgeNum, double dx[3],
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);
  
  // for structure-fluid "half-Riemann" problem
  int computeFSIRiemannSolution(double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej, int Id = 0);
  void computeFSIRiemannJacobian(double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej, double* dWdW,int Id = 0);
  void computeFSIRiemannSolution(int tag, double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej); //TODO:not needed!


  void reset(int it);
  void resetInterfacialW(int edgeNum);
	
};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExactRiemannSolver.C>
#endif

#endif
