#ifndef _KSP_MG_PREC_H_
#define _KSP_MG_PREC_H_

#include <DistVector.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MultiGridLevel.h>
#include <MvpMatrix.h>
#include <MultiGridSmoothingMatrix.h>
#include <SpaceOperator.h>

class Domain;
class DistGeoState;
template<class Scalar,int dim> class DistSVec;

template<class Scalar,int dim> class MultiGridPrecMatVecProd;
template<class Scalar,int dim,class Scalar2 = double> class MultiGridPrecJacobiPrec;

template<class Scalar, int dim,class Scalar2 = double>
class MultiGridPrec : public KspPrec<dim, Scalar2>, public DistMat<Scalar2,dim> {

  int nSmooth;

  double relaxationFactor;

  const int num_levels, agglom_size, numLocSub;
  MultiGridLevel<Scalar2> ** multiGridLevels;
  MultiGridSmoothingMatrix<Scalar2,dim> ***smoothingMatrices;

  DistSVec<Scalar2, dim> ** macroValues, ** macroValuesTmp,**macroValuesOld;
  DistSVec<Scalar2, dim> ** macroR;
  DistVec<Scalar2> ** macroIrey;
  DistSVec<Scalar2, dim> ** macroDX;

  DistMat<Scalar2,dim> ** macroA;

  DistGeoState& geoState;

  KspSolver<DistSVec<Scalar2,dim>, 
            MultiGridPrecMatVecProd<Scalar2,dim> ,
            MultiGridPrecJacobiPrec<Scalar2,dim>, Communicator>* coarseSolver;

  MultiGridPrecMatVecProd<Scalar2,dim>* coarseMvp;
  MultiGridPrecJacobiPrec<Scalar2,dim>* coarsePrec;

  bool ownsFineA;

  Domain* domain;

public:

  MultiGridPrec(Domain *, DistGeoState &, PcData&,
                KspData&,bool createFineA,int ** = 0, BCApplier* =0);
  ~MultiGridPrec();

  void setup();

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  void getData(DistMat<Scalar2,dim>& mat);/*,
               DistSVec<Scalar2,dim>& V,
               SpaceOperator<dim>& spo,
               DistTimeState<dim>* timeState);*/

  void smooth(int level,DistSVec<Scalar2,dim>& x,
              const DistSVec<Scalar2,dim>& f);
  
  DistMat<Scalar2,dim> &operator= (const Scalar2 s) { return macroA[0]->operator = (s); }

  GenMat<Scalar2,dim> &operator() (int i) { return macroA[0]->operator()(i); }
};

#endif
