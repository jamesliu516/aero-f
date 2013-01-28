#ifndef _KSP_MG_PREC_H_
#define _KSP_MG_PREC_H_

#include <DistVector.h>
#include <KspPrec.h>

class Domain;
class DistGeoState;
template<class Scalar,int dim> class DistSVec;

template<class Scalar, int dim, class Scalar2 = double>
class MultiGridPrec : public KspPrec<dim, Scalar2> {
  const int num_levels,agglom_size,numLocSub;
  DistMacroCellSet *macroCells;
  DistSVec<Scalar2, dim> **macroValues;

  DistGeoState& geoState;
  DistInfo ** nodeDistInfo;
  DistInfo ** edgeDistInfo;
public:

  MultiGridPrec(Domain *, DistGeoState &, int ** = 0, BCApplier* =0);
  ~MultiGridPrec();

  void setup();

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);
};

#endif
