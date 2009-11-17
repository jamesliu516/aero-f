#ifndef _EXACT_RIEMANN_SOLVER_H
#define _EXACT_RIEMANN_SOLVER_H

class IoData;
class LocalRiemann;
class VarFcn;

template<class Scalar, int dim> class SVec;


//------------------------------------------------------------------------------
template<int dim>
class ExactRiemannSolver{

  LocalRiemann *lriemann;
	
  int iteration;
  SVec<double,dim>  &rupdate;
  Vec<double>       &weight;

  public:

  ExactRiemannSolver(IoData &, SVec<double,dim> &, Vec<double> &, VarFcn *);
  ~ExactRiemannSolver();


  SVec<double,dim> &getRiemannUpdate() const { return rupdate; }
  Vec<double> &getRiemannWeight() const { return weight; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi, VarFcn *vf,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              int i, int j, double dx[3]);

  void reset(int it);
	
};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExactRiemannSolver.C>
#endif

#endif
