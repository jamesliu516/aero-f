#ifndef _IMPLICIT_Gappy_TS_DESC_H_
#define _IMPLICIT_Gappy_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGappyTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	// A, B, sampleNodes
	// nSampleNodes, nJ, nR
	// scalapack stuff

	int nSampleNodes, nJ, nR;

  VecSet<DistSVec<double, dim> > A, B;

	std::vector <int> sampleNodes;

	void solveNewtonSystem(const int &it, double &res, bool &breakloop);

public:
  
  ImplicitGappyTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGappyTsDesc.C>
#endif

#endif
