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

  public:

  ExactRiemannSolver(IoData &, SVec<double,dim> &, Vec<double> &, 
                     VarFcn *, SparseGridCluster *);
  ~ExactRiemannSolver();


  SVec<double,dim> &getRiemannUpdate() const { return rupdate; }
  Vec<double> &getRiemannWeight() const { return weight; }

  // for multiphase Riemann problem
  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi, VarFcn *vf,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              int i, int j, double dx[3]);

  // for structure-fluid "half-Riemann" problem
  void computeFSIRiemannSolution(double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej, int Id = 0);
  void computeFSIRiemannSolution(int tag, double *Vi, double *Vstar, double *nphi, 
                                 VarFcn *vf, double *Wstar, int nodej); //TODO:not needed!


  void reset(int it);
	
};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExactRiemannSolver.C>
#endif

#endif
