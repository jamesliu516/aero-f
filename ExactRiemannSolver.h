#ifndef _EXACT_RIEMANN_SOLVER_H
#define _EXACT_RIEMANN_SOLVER_H

class IoData;
class LocalRiemann;
class VarFcn;
class SparseGridCluster;

template<class Scalar, int dim> class SVec;


//------------------------------------------------------------------------------
template<int dim>
class ExactRiemannSolver{

  LocalRiemann *lriemann;

  LocalRiemann *fsiRiemann;
	
  int iteration;
  SVec<double,dim>  &rupdate;
  Vec<double>       &weight;
  SVec<double,dim-2>  &interfacialWi;
  SVec<double,dim-2>  &interfacialWj;

  public:

  ExactRiemannSolver(IoData &, SVec<double,dim> &, Vec<double> &, 
                     SVec<double,dim-2> &, SVec<double,dim-2> &,
                     VarFcn *, SparseGridCluster *);
  ~ExactRiemannSolver();


  SVec<double,dim> &getRiemannUpdate() const { return rupdate; }
  Vec<double> &getRiemannWeight() const { return weight; }

  void storePreviousPrimitive(SVec<double,dim> &V, Vec<int> &fluidId, SVec<double,3> &X);
  void updatePhaseChange(SVec<double,dim> &V, Vec<int> &fluidId, Vec<int> &fluidIdn);

  // for multiphase Riemann problem
  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi, VarFcn *vf,
                              double *Wi, double *Wj,
                              int i, int j, int edgeNum, double dx[3]);

  // for structure-fluid "half-Riemann" problem
  void computeFSIRiemannSolution(double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej);
  void computeFSIRiemannSolution(int tag, double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej);


  void reset(int it);
  void resetInterfacialW(int edgeNum);
	
};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExactRiemannSolver.C>
#endif

#endif
