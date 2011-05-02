#ifndef _IMPLICIT_PROJ_ERROR_TS_DESC_H_
#define _IMPLICIT_PROJ_ERROR_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitProjErrorTsDesc : public ImplicitRomTsDesc<dim> {

protected:

  DistSVec<double, dim> Uold, Unew;
  DistSVec<double, dim> dUtrue, dUproj, dUerr;
  Vec<double> projDU;
	double dUerrnorm, dUrelerr;

	//jac is Phi^TPhi

	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q);
	virtual void computeAJ(int it, DistSVec<double, dim> &Q);
	virtual void postProStep(DistSVec<double,dim> &, int);	// by default, do not do post processing
	virtual void writeStateRomToDisk(int it, double cpu);	// write out projected dU
	virtual void writeErrorToDisk(int it, double cpu);	// write out projected dU

public:
  
  ImplicitProjErrorTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitProjErrorTsDesc.C>
#endif

#endif
