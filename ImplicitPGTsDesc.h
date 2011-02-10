#ifndef _IMPLICIT_PG_TS_DESC_H_
#define _IMPLICIT_PG_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <ParallelRom.h>
#include <VectorSet.h>
#include <DistVector.h>
//------------------------------------------------------------------------------

template<int dim>
class ImplicitPGTsDesc : public ImplicitRomTsDesc<dim> {

private:
  typedef VecSet< DistSVec<double,dim> > SetOfVec;
  SetOfVec rhsVS;
  double **lsCoeff;

protected:

  ParallelRom<dim> parallelRom;
  Vec<double> rhs;
  Vec<double> From;
  void saveAJsol(const int);	// only implemented for PG
	void solveNewtonSystem(const int &it, double &res, bool &breakloop);

public:
  
  ImplicitPGTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitPGTsDesc();
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitPGTsDesc.C>
#endif

#endif
