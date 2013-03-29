#ifndef IMPLICIT_GAPPY_TS_DESC_H_
#define IMPLICIT_GAPPY_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGnatTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	// resMat, jacMat, sampleNodes
	// nSampleNodes, nPodJac
	// scalapack stuff
	// nPod in ImplicitRomTsDesc

	bool performRestriction;
	int nSampleNodes, nPodJac;	//nPodJac specified under ROB{ NumROB2 }
	std::vector<int> sampleNodes;
	const char *gnatPrefix;

	int numResJacMat ;	// number of matrices for A and B (1 if they use the same)
  //std::auto_ptr< VecSet < DistSVec < double, dim> > > resMat, jacMat;
  VecSet < DistSVec < double, dim> > *resMat, *jacMat;
	std::auto_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
	std::auto_ptr< DistSVec<double, dim> > ResRestrict;

	DistLeastSquareSolver leastSquaresSolver;
	std::auto_ptr< RestrictionMapping<dim> > restrictionMapping;

	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	const DistInfo & getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();};
	const RestrictionMapping<dim> * restrictMapping() const { return restrictionMapping.get(); } 

	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q);
	virtual void computeAJ(int, DistSVec<double, dim> &);
	virtual bool breakloop1(const bool);
	virtual bool breakloop2(const bool);

	double *jactmp, *column;
public:
  
  ImplicitGnatTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitGnatTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGnatTsDesc.C>
#endif

#endif
