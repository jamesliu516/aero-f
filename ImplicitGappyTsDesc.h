#ifndef _IMPLICIT_GAPPY_TS_DESC_H_
#define _IMPLICIT_GAPPY_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <GappyOnlineTsDesc.h>
#include <DistLeastSquareSolver.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGappyTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	// Amat, Bmat, sampleNodes
	// nSampleNodes, nPodJac
	// scalapack stuff
	// nPod in ImplicitRomTsDesc

	int nSampleNodes, nPodJac;	//nPodJac specified under ROB{ NumROB2 }
	std::vector<int> sampleNodes;

  std::auto_ptr< VecSet < DistSVec < double, dim> > > Amat, Bmat;
	std::auto_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
	std::auto_ptr< DistSVec<double, dim> > ResRestrict;

	DistLeastSquareSolver leastSquaresSolver;
	std::auto_ptr< RestrictionMapping<dim> > restrictionMapping;

	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	void readSampleNodes(const char *sampleNodeFileName);
	const DistInfo & getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();};
	const RestrictionMapping<dim> * restrictMapping() const { return restrictionMapping.get(); } // TODO

	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q);
	virtual void computeAJ(int it, DistSVec<double, dim> &Q);

public:
  
  ImplicitGappyTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGappyTsDesc.C>
#endif

#endif
