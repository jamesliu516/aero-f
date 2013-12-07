#ifndef _IMPLICIT_PG_TS_DESC_H_
#define _IMPLICIT_PG_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <ParallelRom.h>
#include <RefVector.h>
#include <VectorSet.h>
#include <DistVector.h>
//------------------------------------------------------------------------------

template<int dim>
class ImplicitPGTsDesc : public ImplicitRomTsDesc<dim> {

private:
  double **lsCoeff;
  RefVec<DistSVec<double, dim> >residualRef;
  int currentProblemSize; // for local rom  

protected:

  ParallelRom<dim> *parallelRom;
  Vec<double> rhs;
  Vec<double> From;
  void saveNewtonSystemVectors(const int totalTimeSteps) {this->saveNewtonSystemVectorsAction(totalTimeSteps);}
  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int& totalTimeSteps = 0);
	int lsSolver;
	double *jactmp;
  KspPrec<dim> *pc;

  DistSVec<double, dim>* A_Uinlet;
  Vec<double>* PhiT_A_Uinlet;
  Vec<double>* PhiT_A_U;
  VecSet< Vec<double> >* PhiT_A_Phi;
  VecSet<DistSVec<double, dim> >* A_Phi;
  double regWeightConstant;
  double Kc;
  double regWeightProportional;
  double Kp;
  double regWeightIntegral;
  double Ki;
  double Ki_leak;
  double dKi;
  double dRegWeightIntegral;
  double ffError;
  double ffErrorPrev;
  double ffErrorTol;
  double minRes;

  double dt;
  double dtInit;
  double rhsNormInit;

  void setProblemSize(DistSVec<double, dim> &);

  void updateLeastSquaresWeightingVector();

  double computePGResidualNorm(DistSVec<double,dim> &);
  void setReferenceResidual();

  int controlNodeGlobalID;
  int controlNodeLocalID;
  int controlNodeSubDomain;
  int controlNodeCpuNum;

public:
  bool checkForLastIteration(IoData &, int, double, double, DistSVec<double,dim> &);
  void monitorInitialState(int, DistSVec<double,dim> &);
  bool monitorConvergence(int, DistSVec<double,dim> &);
  
  ImplicitPGTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitPGTsDesc();
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitPGTsDesc.C>
#endif

#endif
