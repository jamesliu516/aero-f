#ifndef IMPLICIT_GAPPY_TS_DESC_H_
#define IMPLICIT_GAPPY_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>
#include <NonlinearRomOnlineIII.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGnatTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	int nSampleNodes, nPodJac;	//nPodJac specified under ROB{ NumROB2 }
	std::vector<int> sampleNodes;

	int numResJacMat ;	// number of matrices for A and B (1 if they use the same)
	std::auto_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
	std::auto_ptr< DistSVec<double, dim> > ResRestrict;

	DistLeastSquareSolver leastSquaresSolver;

  void setProblemSize(DistSVec<double, dim> &);

	void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int & totalTimeSteps = 0);

	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> *R=NULL);
	virtual void computeAJ(int, DistSVec<double, dim> &);
	virtual bool breakloop1(const bool);
	virtual bool breakloop2(const bool);

  double computeGnatResidualNorm(DistSVec<double,dim> &);
  void setReferenceResidual();
  void deleteRestrictedQuantities();

	double *jactmp, *column;

  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);

public:
  
  bool checkForLastIteration(IoData &, int, double, double, DistSVec<double,dim> &);
  void monitorInitialState(int, DistSVec<double,dim> &);
  bool monitorConvergence(int, DistSVec<double,dim> &);

  ImplicitGnatTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitGnatTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGnatTsDesc.C>
#endif

#endif
