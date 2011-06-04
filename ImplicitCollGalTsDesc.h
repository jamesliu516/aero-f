#ifndef IMPLICIT_COLL_GAL_TS_DESC_H_
#define IMPLICIT_COLL_GAL_TS_DESC_H_

#include <ImplicitGappyTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCollGalTsDesc : public ImplicitGappyTsDesc<dim> {

protected:

  Vec<double> From;
  Vec<double> rhs;
	std::auto_ptr< VecSet < DistSVec<double, dim> > > podRestrict;
	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	double *jactmp;

public:
  
  ImplicitCollGalTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCollGalTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCollGalTsDesc.C>
#endif

#endif
