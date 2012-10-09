#ifndef _KSP_MG_PREC_H_
#define _KSP_MG_PREC_H_

#include <DistVector.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MultiGridLevel.h>
#include <MvpMatrix.h>
#include <MultiGridSmoothingMatrix.h>
#include <SpaceOperator.h>
#include <MultiGridKernel.h>

class Domain;
class DistGeoState;

template<class Scalar, int dim,class Scalar2 = double>
class MultiGridPrec : public KspPrec<dim, Scalar2>, public DistMat<Scalar2,dim> {

public:

  MultiGridPrec(Domain *, DistGeoState &, PcData&,
                KspData&,IoData&, VarFcn*,bool createFineA,
                DistTimeState<dim>*,int ** = 0, BCApplier* =0);
  ~MultiGridPrec();

  void initialize();

  bool isInitialized(); 

  void setOperators(SpaceOperator<dim>*);

  void setup();

  void setup(DistSVec<Scalar2,dim>&);

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  void getData(DistMat<Scalar2,dim>& mat);/*,
               DistSVec<Scalar2,dim>& V,
               SpaceOperator<dim>& spo,
               DistTimeState<dim>* timeState);*/

  DistMat<Scalar2,dim> &operator= (const Scalar2 s);// { return macroA[0]->operator = (s); }

  GenMat<Scalar2,dim> &operator() (int i);// { return macroA[0]->operator()(i); }

 private:

  PcData& pcData;

  MultiGridKernel<Scalar,dim,Scalar2>* mgKernel;
};

#endif
