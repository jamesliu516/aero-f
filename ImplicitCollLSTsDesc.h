#ifndef IMPLICIT_COLL_LS_TS_DESC_H_
#define IMPLICIT_COLL_LS_TS_DESC_H_

#include <ImplicitGappyTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCollLSTsDesc : public ImplicitGappyTsDesc<dim> {

protected:
  Vec<double> From;
  Vec<double> rhs;
	double *jactmp;

	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloopImplicitGalerkinTsDesc); 

public:
  
  ImplicitCollLSTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCollLSTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCollLSTsDesc.C>
#endif

#endif

