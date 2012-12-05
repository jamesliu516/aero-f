#ifndef _IMPLICIT_BROYDEN_TS_DESC_H_
#define _IMPLICIT_BROYDEN_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitBroydenTsDesc : public ImplicitRomTsDesc<dim> {

protected:

  int JacSkipNewton;
	int totTimeSteps;	// global number of time steps 

	Vec<double> rhs;
  Vec<double> From;
  Vec<double> dFrom;
  Vec<double> Fromold;
  void computeAJ(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time
  void updateGlobalTimeSteps(const int _it) {totTimeSteps = _it-1;};
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	void broydenUpdate(Vec<double> &);
	double *jactmp;

public:
  
  ImplicitBroydenTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitBroydenTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitBroydenTsDesc.C>
#endif

#endif
