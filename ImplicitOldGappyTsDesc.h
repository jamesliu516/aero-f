#ifndef _IMPLICIT_OLD_GAPPY_TS_DESC_H_
#define _IMPLICIT_OLD_GAPPY_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProdFD;


//------------------------------------------------------------------------------

template<int dim>
class ImplicitOldGappyTsDesc : public ImplicitRomTsDesc<dim> {

protected:

  int nIntNodes, dimInterpMat;
  int *globalSubSet;
  int *locNodeSet;
  FullM interpMat1;
  FullM interpMat2;
  int myNNodeInt;
  int *myLocalSubSet;
  int *myInterpNodes;

	void solveNewtonSystem(const int &it, int &res, bool &breakloop);
	void computeRestrictInfo();
	void computeAJGappy(int, DistSVec<double, dim> &, DistSVec<double, dim> &,
			Vec<double> &, VecSet<DistSVec<double, dim> > &);
	void rowPartition(int &myMinIndex, int &myMaxIndex, int nRow, int sym );

public:
  
  ImplicitOldGappyTsDesc(IoData &, GeoSource &, Domain *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitOldGappyTsDesc.C>
#endif

#endif
