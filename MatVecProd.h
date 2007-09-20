#ifndef _MAT_VEC_PROD_H_
#define _MAT_VEC_PROD_H_

#include <IoData.h>
#include <MvpMatrix.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <complex.h>
typedef complex<double> bcomp;

class VarFcn;
class FluxFcn;
class DistGeoState;
class MemoryPool;
class LevelSet;

template<int dim> class RecFcnConstant;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------------

template<int dim, int neq>
class MatVecProd {

public:

  MatVecProd() {}
  virtual ~MatVecProd() {}

  virtual void exportMemory(MemoryPool *mp) {}

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
			DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  virtual void evaluateLS(int, DistSVec<double,3> &, DistVec<double> &,
                          DistVec<double> &, DistVec<double> &, DistVec<double> &,
                          DistVec<double> &, DistSVec<double,dim> &,
                          DistVec<double> &) { };
  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
			                  DistSVec<double,dim> &, DistVec<double> &,
                        DistExactRiemannSolver<dim> *,  
                        DistSVec<double,dim> &)  { };

  virtual void apply(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;
  virtual void apply(DistVec<double> &, DistVec<double> &) { };

  virtual void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void applyT(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;

  virtual void applyLS(DistVec<double> &, DistVec<double> &) { };
  virtual void applyLS(DistSVec<double,neq> &, DistSVec<double,neq> &) { };
  virtual void applyLS(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) { };

// Included (MB)
  virtual void evaluateInviscid(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  virtual void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  virtual void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0) = 0;

};

//------------------------------------------------------------------------------

template<int dim, int neq>
class MatVecProdFD : public MatVecProd<dim,neq> {

  DistSVec<double,dim> Qeps;
  DistSVec<double,dim> Feps;
  DistSVec<double,neq> Qepstmp;
  DistSVec<double,neq> Fepstmp;

  SpaceOperator<dim> *spaceOp;
  RecFcnConstant<dim> *recFcnCon;
  DistSVec<double,dim> *Rn;

  DistTimeState<dim> *timeState;
  DistGeoState *geoState;
  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,neq> Q;
  DistSVec<double,neq> F;

  DistVec<double> *Phi;
  DistExactRiemannSolver<dim> *Riemann;

  Communicator *com;

// Included (MB)
  DistSVec<double,neq> Ftmp;
  IoData *iod;
  int fdOrder;
  double fdeps;

private:

  double computeEpsilon(DistSVec<double,neq> &, DistSVec<double,neq> &);

public:

// Included (MB)
  MatVecProdFD(ImplicitData &, DistTimeState<dim> *, DistGeoState *, 
	       SpaceOperator<dim> *, Domain *, IoData &, bool = false);

  ~MatVecProdFD();

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                     DistSVec<double,dim> &, DistVec<double> &,
                     DistExactRiemannSolver<dim> *, 
                     DistSVec<double,dim> &);
  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    cout << "... ERROR: ::apply function not implemented for class MatVecProdFD with complex arguments" << endl; }

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y)  {
    cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);
  
};

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
class MatVecProdH1 : public MatVecProd<dim,neq>, public DistMat<Scalar,neq> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  MvpMat<Scalar,neq> **A;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;

public:

  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *);
  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *, IoData &);
  ~MatVecProdH1();

  DistMat<Scalar,neq> &operator= (const Scalar);

  GenMat<Scalar,neq> &operator() (int i) { return *A[i]; }

  void exportMemory(MemoryPool *);

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistVec<double> &,
                DistExactRiemannSolver<dim> *, 
                DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &);

  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    cout << "... ERROR: ::apply function not implemented for class MatVecProdH1 with complex arguments" << endl; }

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" << endl; } 
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y) { 
    cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" <<
 endl; }

// Included (MB)
  void evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv, DistSVec<double,dim> &q, DistSVec<double,dim> &f) { this->com->printf(2, "*** Error: function evaluateInviscid not implemented\n");}
  void evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv, DistSVec<double,dim> &q, DistSVec<double,dim> &f) { this->com->printf(2, "*** Error: function evaluateViscous not implemented\n");}
  void applyInviscid(DistSVec<double,neq> &x, DistSVec<double,neq> &y) { this->com->printf(2, "*** Error: function applyInviscid not implemented\n");}
  void applyViscous(DistSVec<double,neq> &x, DistSVec<double,neq> &y) { this->com->printf(2, "*** Error: function applyViscous not implemented\n");}
  void rstSpaceOp(IoData &iod, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);

};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class MatVecProdH2 : public MatVecProd<dim,dim>, public DistMat<Scalar,dim> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  // format of A is the diagonal terms and then the edge contributions
  MvpMat<Scalar,dim> **A;  

  // coefficients in the linearization of the reconstructed-limited primitive states
  /* see Lesoinne et al. A linearized method for the frequency analysis of three 
     dimensional fluid/structure interaction problems in all flow regimes, Comp. Meth. Appl. Mech. Eng.
     vol. 190 (2001) pp 3121-3146 */
  DistSVec<double,dim> aij; 
  DistSVec<double,dim> aji;
  DistSVec<double,dim> bij;
  DistSVec<double,dim> bji;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;
  FluxFcn **fluxFcn;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *Q;
  DistSVec<double,dim> *F;

  // viscous flux jacobian and/or 1st order jac terms(ie face flux jac)
  MatVecProdH1<dim, Scalar ,dim> *R;

// Included (MB)
  MatVecProdFD<dim,dim> *RFD;
  DistSVec<double,dim> *vProd;
  int viscJacContrib;

public:

// Included (MB)
  MatVecProdH2(IoData &, VarFcn *, DistTimeState<dim> *, 
	       SpaceOperator<dim> *, Domain *, DistGeoState * = 0);

  ~MatVecProdH2();

  DistMat<Scalar,dim> &operator= (const Scalar);

  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

/*  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);
  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &, 
		 DistSVec<double,dim> &, DistSVec<double,dim> &); */
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluatestep1(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);

  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate2step1(int, DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &);


  void evalH(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);

  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void apply(DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &);

  void applyT(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void applyT(DistSVec<bcomp,dim> &x, DistSVec<bcomp,dim> &y);

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void applyInviscid(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void applyViscous(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);

};


//------------------------------------------------------------------------------
                                                                                                                      
template<class Scalar, int dim, int neq>
class MatVecProdLS : public MatVecProd<dim, neq>, public DistMat<Scalar,neq> {
                                                                                                                      
#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  // format of A is the diagonal terms and then the edge contributions
  MvpMat<Scalar,neq> **A;
                                                                                                                      
  // coefficients in the linearization of the reconstructed-limited primitive states
  /* see Lesoinne et al. A linearized method for the frequency analysis of three
     dimensional fluid/structure interaction problems in all flow regimes, Comp. Meth. Appl. Mech. Eng.
     vol. 190 (2001) pp 3121-3146 */

  DistSVec<double,neq> aij;
  DistSVec<double,neq> aji;
  DistSVec<double,neq> bij;
  DistSVec<double,neq> bji;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;
  DistGeoState *geoState;
	LevelSet *LS;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistVec<double> *Q;
  DistVec<double> *Qn;
  DistVec<double> *Qnm1;
  DistVec<double> *Qnm2;
  DistSVec<double,dim> *U;
  DistVec<double> *F;
  DistVec<double> Qeps;
  DistVec<double> QepsV;
  DistVec<double> QV;
  DistVec<double> Feps;

public:

  MatVecProdLS(IoData &, VarFcn *, DistTimeState<dim> *, DistGeoState *,
               SpaceOperator<dim> *, Domain *, LevelSet *);
  ~MatVecProdLS();
                                                                                                                      
  DistMat<Scalar,neq> &operator= (const Scalar);
                                                                                                                      
  GenMat<Scalar,neq> &operator() (int i) { return *A[i]; }
                                                                                                                      
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &) { };
  void evaluateLS(int, DistSVec<double,3> &, DistVec<double> &,
                DistVec<double> &, DistVec<double> &, DistVec<double> &,
                DistVec<double> &, DistSVec<double,dim> &,DistVec<double> &);
                                                                                                                      
  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &) { };
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) { };
                                                                                                                      
  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &) { };
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y) { };

  void applyLS(DistVec<double> &, DistVec<double> &);

  double computeEpsilon(DistVec<double> &, DistVec<double> &);

// Included (MB)
  void evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv, DistSVec<double,dim> &q, DistSVec<double,dim> &f) { this->com->printf(2, "*** Error: function evaluateInviscid not implemented\n");}
  void evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv, DistSVec<double,dim> &q, DistSVec<double,dim> &f) { this->com->printf(2, "*** Error: function evaluateViscous not implemented\n");}
  void applyInviscid(DistSVec<double,neq> &x, DistSVec<double,neq> &y) { this->com->printf(2, "*** Error: function applyInviscid not implemented\n");}
  void applyViscous(DistSVec<double,neq> &x, DistSVec<double,neq> &y) { this->com->printf(2, "*** Error: function applyViscous not implemented\n");}
  void rstSpaceOp(IoData &iod, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0) { this->com->printf(2, "*** Error: function rstSpaceOp not implemented\n");}

};
                                                                                                                      
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatVecProd.C>
#endif

#endif

