#ifndef _IMPLICIT_GALERKIN_TS_DESC_H_
#define _IMPLICIT_GALERKIN_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGalerkinTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	Vec<double> rhs;
  Vec<double> From;
  void saveNewtonSystemVectors(const int _it) {this->saveNewtonSystemVectorsAction(_it);}
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);

public:
  
  ImplicitGalerkinTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGalerkinTsDesc.C>
#endif

#endif
