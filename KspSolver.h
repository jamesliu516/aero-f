#ifndef _KSP_SOLVER_H_
#define _KSP_SOLVER_H_

#include <cstdio>
#include <Vector.h>
#include <DenseMatrix.h>
#include <VectorSet.h>

class KspData;
class KspConvCriterion;

#include <complex>
typedef std::complex<double> bcomp;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
class KspSolver {

protected:

  double eps;
  bool checkFinalRes;

  MatVecProdOp *mvpOp;
  PrecOp *pcOp;
  IoOp *ioOp;

  KspConvCriterion *kspConvCriterion;

  FILE *output;

public:
  int maxits;

  KspSolver() {}
  KspSolver(KspData &, MatVecProdOp *, PrecOp *, IoOp *);
  virtual ~KspSolver() { if (kspConvCriterion) delete kspConvCriterion; }

  void setup(int, int, VecType &);

  virtual int solve(VecType &, VecType &) = 0;
  virtual int solveLS(VecType &, VecType &) { return 0; };

  void printParam() { ioOp->fprintf(stderr, " solver params: %d maxits, %e eps\n", maxits, eps);  }
};

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
class RichardsonSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp> {

  VecType dx, r;

public:

  RichardsonSolver(const typename VecType::InfoType &, KspData &, 
		   MatVecProdOp *, PrecOp *, IoOp *);
  ~RichardsonSolver() {}

  int solve(VecType &, VecType &);
  int solveLS(VecType &, VecType &);

};

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
class CgSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp> {

  VecType r, Ap, y, p;

public:

  CgSolver(const typename VecType::InfoType &, KspData &,
	   MatVecProdOp *, PrecOp *, IoOp *);
  ~CgSolver() {}

  int solve(VecType &, VecType &);
  int solveLS(VecType &, VecType &);

};

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class
		ScalarT = double>
class GmresSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT> {

  int numVec;

  GenFullM<ScalarT> cs, H;
  Vec<ScalarT> g, y;

  VecSet<VecType> V;
  VecType w, r;

public:

  GmresSolver(const typename VecType::InfoType &, KspData &, 
	      MatVecProdOp *, PrecOp *, IoOp *);
  ~GmresSolver() {}

  int solve(VecType &, VecType &);
  int solveLS(VecType &, VecType &);
  int solveT(VecType &, VecType &);
  void solve(VecSet<VecType> &, VecSet<VecType> &);

  void applyPreviousRotations(int, GenFullM<ScalarT> &, GenFullM<ScalarT> &);
  void applyNewRotation(int, GenFullM<ScalarT> &, GenFullM<ScalarT> &, Vec<ScalarT> &);
  void backwardSolve(int, GenFullM<ScalarT> &, Vec<ScalarT> &, Vec<ScalarT> &);

};

//------------------------------------------------------------------------------
                                                        
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class
                ScalarT = double>
class GcrSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT> {
                                                        
                                                        
                                                        
  int numVec;
                                                        
  ScalarT alpha, beta;
                                                        
  ScalarT *ApAp;
                                                        
  ScalarT *y;                                                                                                                  
  VecSet<VecType> p, Ap;
                                                        
  VecType w, r, R, AR, temp, w0, x0;
                                                        
                                                        
                                                        
public:
                                                        
  int numCalcVec;                                                                                                                  
  GcrSolver(const typename VecType::InfoType &, KspData &,
              MatVecProdOp *, PrecOp *, IoOp *);
  ~GcrSolver() {}
                                                                                                                                                                          
  int solve(VecType &, VecType &);
  int solveMRhs(VecType &, VecType &);
                                                        
                                                        
                                                        
};
                                                                                                                                                                          
//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <KspSolver.C>
#endif

#endif
