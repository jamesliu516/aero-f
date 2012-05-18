#ifndef _KSP_MG_PREC_H_
#define _KSP_MG_PREC_H_

#include <DistVector.h>
#include <KspPrec.h>
#include <MultiGridLevel.h>
#include <MvpMatrix.h>
#include <MultiGridSmoothingMatrix.h>

class Domain;
class DistGeoState;
template<class Scalar,int dim> class DistSVec;

template<class Scalar, int dim, class Scalar2 = double>
class MultiGridPrec : public KspPrec<dim, Scalar2> {

  int nSmooth;

  double relaxationFactor;

  const int num_levels, agglom_size, numLocSub;
  MultiGridLevel<Scalar2> ** multiGridLevels;
  MultiGridSmoothingMatrix<Scalar2,dim> ***smoothingMatrices;

  DistSVec<Scalar2, dim> ** macroValues, ** macroValuesTmp;
  DistSVec<Scalar2, dim> ** macroR;
  DistSVec<Scalar2, dim> ** macroDX;

  DistMat<Scalar2,dim> ** macroA;

  DistGeoState& geoState;
public:

  MultiGridPrec(Domain *, DistGeoState &, int ** = 0, BCApplier* =0);
  ~MultiGridPrec();

  void setup();

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  void getData(DistMat<Scalar2,dim>& mat);

  void smooth(int level,DistSVec<Scalar2,dim>& x,
              const DistSVec<Scalar2,dim>& f);

};

#endif
