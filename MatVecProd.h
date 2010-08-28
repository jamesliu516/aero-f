#ifndef _MAT_VEC_PROD_H_
#define _MAT_VEC_PROD_H_

#include <IoData.h>
#include <MvpMatrix.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <complex>
typedef std::complex<double> bcomp;

class VarFcn;
class FluxFcn;
class DistGeoState;
class MemoryPool;

template<int dimLS> class LevelSet;
template<int dim> class RecFcnConstant;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim, int dimLS> class MultiPhaseSpaceOperator;
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

  virtual void apply(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;
  virtual void apply(DistVec<double> &, DistVec<double> &) { };

  virtual void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void applyT(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;

// Included (MB)
  virtual void evaluateInviscid(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &){
    std::cout<<"*** Error: function evaluateInviscid not implemented"<<std::endl;}
  virtual void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &){
    std::cout<<"*** Error: function evaluateViscous not implemented"<<std::endl;}
  virtual void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &){
    std::cout<<"*** Error: function applyInviscid not implemented"<<std::endl;}
  virtual void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &){
    std::cout<<"*** Error: function applyViscous not implemented"<<std::endl;}
  virtual void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0){
    std::cout<<"*** Error: function rstSpaceOp not implemented"<<std::endl;}

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
	       SpaceOperator<dim> *, Domain *, IoData &);

  ~MatVecProdFD();

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    std::cout << "... ERROR: ::apply function not implemented for class MatVecProdFD with complex arguments" << endl; }

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }

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
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &);

  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    std::cout << "... ERROR: ::apply function not implemented for class MatVecProdH1 with complex arguments" << endl; }

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" << endl; } 
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y) { 
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" <<
 endl; }

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);
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
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &, DistVec<int> &,
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

//----------------------------------------------------------------------------//
//                MatVecProd for Multiphase Euler equations                   //
//----------------------------------------------------------------------------//
// this class is not derived from the MatVecProd class!
// Its templates are different and additional members are different.
template<int dim, int dimLS>
class MatVecProdMultiPhase {

protected:
  DistTimeState<dim> *timeState;
  MultiPhaseSpaceOperator<dim,dimLS> *spaceOp;
  DistExactRiemannSolver<dim> *riemann;
  FluidSelector *fluidSelector;

public:

  MatVecProdMultiPhase(DistTimeState<dim> *ts, MultiPhaseSpaceOperator<dim,dimLS> *spo,
                       DistExactRiemannSolver<dim> *rsolver, FluidSelector *fs) : 
                       timeState(ts), spaceOp(spo), riemann(rsolver), fluidSelector(fs) {}
  virtual ~MatVecProdMultiPhase() { timeState=0; spaceOp = 0; riemann = 0; fluidSelector = 0; }

  virtual void exportMemory(MemoryPool *mp) {}

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
                        DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                        DistSVec<double,dim> &) = 0;

  virtual void apply(DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

};

//------------------------------------------------------------------------------
// Finite Difference method for Multiphase Euler equations
template<int dim, int dimLS>
class MatVecProdFDMultiPhase : public MatVecProdMultiPhase<dim,dimLS> {

  DistSVec<double,dim> Qeps;
  DistSVec<double,dim> Feps;
  DistSVec<double,dim> Q;
  DistSVec<double,dim> F;

  DistGeoState *geoState;
  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dimLS> *Phi;

  Communicator *com;

  double computeEpsilon(DistSVec<double,dim> &, DistSVec<double,dim> &);

public:

// Included (MB)
  MatVecProdFDMultiPhase(DistTimeState<dim> *, DistGeoState *, 
	       MultiPhaseSpaceOperator<dim,dimLS> *, DistExactRiemannSolver<dim> *,
               FluidSelector *, Domain *);

  ~MatVecProdFDMultiPhase();

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                     DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                     DistSVec<double,dim> &);
  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);

};

//------------------------------------------------------------------------------
// H1 MatVecProd for Multiphase Euler equations
template<int dim, int dimLS>
class MatVecProdH1MultiPhase : public MatVecProdMultiPhase<dim,dimLS>,
                               public DistMat<double,dim> {

  MvpMat<double,dim> **A;

public:

  MatVecProdH1MultiPhase(DistTimeState<dim> *, MultiPhaseSpaceOperator<dim,dimLS> *, 
                         DistExactRiemannSolver<dim> *, FluidSelector *, Domain *);
  ~MatVecProdH1MultiPhase();

  DistMat<double,dim> &operator= (const double);

  GenMat<double,dim> &operator() (int i) { return *A[i]; }

  void exportMemory(MemoryPool *);

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                DistSVec<double,dim> &);

  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);
};

//----------------------------------------------------------------------------//
//                MatVecProd for Level Set equations                          //
//----------------------------------------------------------------------------//
// Finite Difference MatVecProd for LevelSet equation
template<int dim, int dimLS>
class MatVecProdLS {
                                                                                                                      
  DistTimeState<dim> *timeState;
  MultiPhaseSpaceOperator<dim,dimLS> *spaceOp;
  DistGeoState *geoState;
  LevelSet<dimLS> *levelSet;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *U;
  DistVec<int> *FluidId;
  DistSVec<double,dimLS> *Q;
  DistSVec<double,dimLS> *F;
  DistSVec<double,dimLS> Qeps;
  DistSVec<double,dimLS> Feps;

  Communicator *com;

  double computeEpsilon(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

public:

  MatVecProdLS(DistTimeState<dim> *, DistGeoState *,
               MultiPhaseSpaceOperator<dim,dimLS> *, Domain *, LevelSet<dimLS> *);
  ~MatVecProdLS();
                                                                                                                      
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dimLS> &, DistSVec<double,dim> &,
                DistSVec<double,dimLS> &, DistVec<int> &);

  void apply(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);


};
                                                                                                                      
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatVecProd.C>
#endif

#endif

