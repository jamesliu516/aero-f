#ifndef _IMPLICIT_BROYDEN_TS_DESC_H_
#define _IMPLICIT_BROYDEN_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitBroydenTsDesc : public ImplicitRomTsDesc<dim> {

protected:

  int JacSkipNewton;
	int Git;	// global number of time steps 

	Vec<double> rhs;
  Vec<double> From;
  Vec<double> dFrom;
  Vec<double> Fromold;
  void computeAJ(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time
  void updateGlobalTimeSteps(int _it) {Git = _it;};	// each ROM has a different way of solving the Newton system
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	void broydenUpdate(Vec<double> &);

public:
  
  ImplicitBroydenTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitBroydenTsDesc.C>
#endif

#endif
