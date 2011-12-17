#ifndef _IMPLICIT_GALERKIN_TS_DESC_H_
#define _IMPLICIT_GALERKIN_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGalerkinTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	Vec<double> rhs;
  Vec<double> From;
	void saveNewtonSystemVectors(const int totalTimeSteps)
		{this->saveNewtonSystemVectorsAction(totalTimeSteps);}
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	double *jactmp;

public:
  
  ImplicitGalerkinTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitGalerkinTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGalerkinTsDesc.C>
#endif

#endif