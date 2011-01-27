#ifndef _IMPLICIT_PG_TS_DESC_H_
#define _IMPLICIT_PG_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <ParallelRom.h>
#include <VectorSet.h>
//------------------------------------------------------------------------------

template<int dim>
class ImplicitPGTsDesc : public ImplicitRomTsDesc<dim> {

protected:

  ParallelRom<dim> parallelRom;
  Vec<double> rhs;
  Vec<double> From;
  void saveAJsol(const int);	// only implemented for PG
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);

public:
  
  ImplicitPGTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitPGTsDesc.C>
#endif

#endif
