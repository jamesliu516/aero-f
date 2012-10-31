#ifndef IMPLICIT_COLL_LS_TS_DESC_H_
#define IMPLICIT_COLL_LS_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <RestrictionMapping.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCollLSTsDesc : public ImplicitRomTsDesc<dim> {

protected:
	int nSampleNodes;
	std::vector<int> sampleNodes;
	std::auto_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
	std::auto_ptr< DistSVec<double, dim> > ResRestrict;
  Vec<double> From;
  Vec<double> rhs;
	double *jactmp;

	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloopImplicitGalerkinTsDesc); 
	std::auto_ptr< RestrictionMapping<dim> > restrictionMapping;
	const DistInfo & getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();};
	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q);
	virtual void computeAJ(int, DistSVec<double, dim> &);
	const RestrictionMapping<dim> * restrictMapping() const { return restrictionMapping.get(); } 
	virtual bool breakloop1(const bool);
	virtual bool breakloop2(const bool);

public:
  
  ImplicitCollLSTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCollLSTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCollLSTsDesc.C>
#endif

#endif

