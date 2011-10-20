#ifndef _IMPLICIT_ROM_POSTPRO_TS_DESC_H_
#define _IMPLICIT_ROM_POSTPRO_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitRomPostproTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	FILE *readRedCoords;	// file of reduced coordinates

	void solveNewtonSystem(const int &it, double &res, bool &breakloop);
	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q);
	virtual void computeAJ(int it, DistSVec<double, dim> &Q);
  DistSVec<double, dim> Uinitial;	// solution increment at EACH NEWTON ITERATION in full coordinates
	virtual void postProStep(DistSVec<double,dim> &, int);	// by default, do not do post processing

public:
  
  ImplicitRomPostproTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomPostproTsDesc.C>
#endif

#endif
