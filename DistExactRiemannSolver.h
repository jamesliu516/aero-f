#ifndef _DIST_EXACT_RIEMANN_H_
#define _DIST_EXACT_RIEMANN_H_


class IoData;
class Domain;

template<class Scalar, int dim> class DistSVec;
template<int dim> class ExactRiemannSolver;
//------------------------------------------------------------------------------

template<int dim>
class DistExactRiemannSolver {

	// updatePhase indicates if we should use riemannupdate (cf below)
	bool updatePhase;
  bool firstpass;
  int numLocSub;

	DistSVec<double,dim> *riemannupdate;
	DistVec<double> *weight;
	// riemannupdate is used only when using GFMPAR
	// it stores the future value at each cell if
	// that cell changes phases between two time steps

  ExactRiemannSolver<dim> **subExactRiemannSolver;


public:
  DistExactRiemannSolver(IoData &iod, Domain *dom);
  ~DistExactRiemannSolver();

	bool DoUpdatePhase() { return updatePhase; }
	DistSVec<double,dim> &getRiemannUpdate() const { return *riemannupdate; }
	DistVec<double> &getRiemannWeight() const { return *weight; }

	ExactRiemannSolver<dim> &operator() (int i) 
				const { return *subExactRiemannSolver[i]; }







};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistExactRiemannSolver.C>
#endif

#endif

