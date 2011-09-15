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

  MatVecProd() : isFSI(false) {}
  virtual ~MatVecProd() {}

  virtual void exportMemory(MemoryPool *mp) {}

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
			DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  virtual void apply(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;
  virtual void apply(DistVec<double> &, DistVec<double> &) { };

  virtual void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &) { }
  
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

  // Structure to enable fluid-structure interaction computations
  struct _fsi {

    DistLevelSetStructure* LSS;
    DistVec<int>* fluidId;
    DistExactRiemannSolver<dim>* riemann;
    bool linRecAtInterface;
    DistSVec<double,3>* Nsbar;
    DistSVec<double,dim>* Wtemp;
    int Nriemann;
    DistVec<GhostPoint<dim>*>* ghostPoints;
  };

  void AttachStructure(const _fsi& f) {
    isFSI = true;
    fsi = f;
  }
 
protected:
  
  // Boolean; set to true if we are using a structure
  bool isFSI;
  _fsi fsi;
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

 
public:

private:

  double computeEpsilon(DistSVec<double,neq> &, DistSVec<double,neq> &);

public:

// Included (MB)
  MatVecProdFD(ImplicitData &, DistTimeState<dim> *, DistGeoState *, 
	       SpaceOperator<dim> *, Domain *, IoData &);
  
  ~MatVecProdFD();

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateRestrict(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &);
  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<dim> &);
  void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &);
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

  void applyViscous(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)
  {
    std::cout << " ... ERROR: MatVecProdFD::applyViscous function not implemented";
    std::cout << " for complex arguments." << std::endl;
    exit(1);
  }

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

  /// Constructor.
  /// \note UH (09/10) This constructor is only called from MatVecProdH2.
  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *);

  /// Constructor.
  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *, IoData &);

  /// Destructor.
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

  void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &);

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" << endl; } 
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y) { 
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" <<
 endl; }

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);
 
  void clearGhost();
};

//------------------------------------------------------------------------------

///
/// The MatVecProdH2 class enables matrix vector product based on the exact
/// Jacobian matrix.
/// For turbulent problems with weak coupling, the parameters dim and neq will
/// be different.
///
/// \note (09/10) The implementation for weak turbulence-model coupling is
/// not optimal because it uses a vector of local length 'dim' for the
/// inviscid part.
/// For the viscous part (laminar + turbulent), the product is using
/// an object MatVecProdH1.
///
template<int dim, class Scalar, int neq>
class MatVecProdH2 : public MatVecProd<dim,neq>, public DistMat<Scalar,dim> {

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
  MatVecProdH1<dim, Scalar, neq> *R;

  // Included (MB)
  MatVecProdFD<dim, neq> *RFD;
  DistSVec<double, neq> *vProd;

  //--------------------------------------
  /// \note (09/10) UH
  /// This nested class allows to specialize the matrix-vector product.
  /// In particular, we differentiate the cases where dim and neq are different
  /// for turbulent problems with weak coupling.
  template<int dd, int nn, class Scalar1, class Scalar2> struct Multiplier
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
      , MatVecProdH1<dd, Scalar1, nn> *R
      , MatVecProdFD<dd, nn> *RFD
      , DistSVec<Scalar2, nn> *vProd
    );
    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
    );
  };

  template<int dd, class Scalar1, class Scalar2>
  struct Multiplier<dd,dd,Scalar1,Scalar2>
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
      , MatVecProdH1<dd, Scalar1, dd> *R
      , MatVecProdFD<dd, dd> *RFD
      , DistSVec<Scalar2, dd> *vProd
    );
    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
    );
  };
  //--------------------------------------


public:

// Included (MB)
  MatVecProdH2(IoData &, VarFcn *, DistTimeState<dim> *, 
	       SpaceOperator<dim> *, Domain *, DistGeoState * = 0);

  ~MatVecProdH2();

  DistMat<Scalar,dim> &operator= (const Scalar);

  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

/*
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);
  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &, 
		 DistSVec<double,dim> &, DistSVec<double,dim> &); 
*/

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &, 
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);

  // UH (09/10)
  // The following function is never called and not implemented.
  //void evaluate(int , DistSVec<double,3> &, DistVec<double> &, DistVec<int> &,
  //              DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);

  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &);

  // UH (09/10)
  // The following functions are never called and not implemented.
  //
  //void evaluatestep1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);
  //void evaluate2step1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);

  void evalH(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);

  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &);

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y);

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  //void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &);
  //void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &);

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
                       timeState(ts), spaceOp(spo), riemann(rsolver), fluidSelector(fs), isFSI(false) {}
  virtual ~MatVecProdMultiPhase() { timeState=0; spaceOp = 0; riemann = 0; fluidSelector = 0; }

  virtual void exportMemory(MemoryPool *mp) {}

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
                        DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                        DistSVec<double,dim> &) = 0;

  virtual void apply(DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  // Structure to enable fluid-structure interaction computations
  struct _fsi {

    DistLevelSetStructure* LSS;
    DistVec<int>* fluidId;
    DistExactRiemannSolver<dim>* riemann;
    bool linRecAtInterface;
    DistSVec<double,3>* Nsbar;
    DistSVec<double,dim>* Wtemp;
    int Nriemann;
    DistVec<GhostPoint<dim>*>* ghostPoints;
  };

  void AttachStructure(const _fsi& f) {
    isFSI = true;
    fsi = f;
  }
 
protected:
  
  // Boolean; set to true if we are using a structure
  bool isFSI;
  _fsi fsi;

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

  int fdOrder;

public:

// Included (MB)
  MatVecProdFDMultiPhase(DistTimeState<dim> *, DistGeoState *, 
			 MultiPhaseSpaceOperator<dim,dimLS> *, DistExactRiemannSolver<dim> *,
			 FluidSelector *, Domain *,IoData& ioData);

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
class MatVecProdLS : public DistMat<double,dimLS> {
                                                                                                                      
  DistTimeState<dim> *timeState;
  MultiPhaseSpaceOperator<dim,dimLS> *spaceOp;
  DistGeoState *geoState;
  LevelSet<dimLS> *levelSet;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *U;
  DistVec<int> *FluidId;
  DistSVec<double,dimLS> *Q;
  DistSVec<double,dim> *V;
  DistSVec<double,dimLS> *F;
  DistSVec<double,dimLS> Qeps;
  DistSVec<double,dimLS> Feps;

  Communicator *com;

  double computeEpsilon(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

  MvpMat<double,dimLS> **A;
public:
  DistMat<double,dimLS> &operator= (const double);
  GenMat<double,dimLS> &operator() (int i) { return *A[i]; }

  MatVecProdLS(DistTimeState<dim> *, DistGeoState *,
               MultiPhaseSpaceOperator<dim,dimLS> *, Domain *, LevelSet<dimLS> *);
  ~MatVecProdLS();
                                                                                                                      
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dimLS> &, DistSVec<double,dim> &,
		DistSVec<double,dim> &,
                DistSVec<double,dimLS> &, DistVec<int> &,bool = false,DistLevelSetStructure* = NULL);

  void apply(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);


};
                                                                                                                      
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatVecProd.C>
#endif

#endif

