#ifndef _IMPLICIT_PG_TS_DESC_H_
#define _IMPLICIT_PG_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <RefVector.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitPGTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	Vec<double> rhs;
  Vec<double> From;
	void saveNewtonSystemVectors(const int totalTimeSteps) {this->saveNewtonSystemVectorsAction(totalTimeSteps);}
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	double *jactmp;

public:
  
  ImplicitPGTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitPGTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitPGTsDesc.C>
#endif

#endif
