#ifndef IMPLICIT_COLL_GAL_TS_DESC_H_
#define IMPLICIT_COLL_GAL_TS_DESC_H_

#include <ImplicitCollLSTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCollGalTsDesc : public ImplicitCollLSTsDesc <dim> {

protected:

	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	std::auto_ptr< VecSet < DistSVec<double, dim> > > podRestrict;

public:
  
  ImplicitCollGalTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCollGalTsDesc.C>
#endif

#endif
