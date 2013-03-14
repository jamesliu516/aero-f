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
  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &);
	int lsSolver;
	double *jactmp;
  KspPrec<dim> *pc;

  DistSVec<double, dim>* A_Uinit;
  Vec<double>* PhiT_A_Uinit;
  Vec<double>* PhiT_A_U;
  VecSet< Vec<double> >* PhiT_A_Phi;
  VecSet<DistSVec<double, dim> >* A_Phi;
  double regCoeff;
  double regThresh; 

  void setProblemSize(DistSVec<double, dim> &);

public:
  
  ImplicitPGTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitPGTsDesc();
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitPGTsDesc.C>
#endif

#endif
