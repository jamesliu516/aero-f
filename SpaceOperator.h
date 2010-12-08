#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include <IoData.h>
#include <GhostPoint.h>
#include <complex>
typedef std::complex<double> bcomp;

class VarFcn;
class FluidSelector;
class BcFcn;
class RecFcn;
class FluxFcn;
class FemEquationTerm;
class VolumicForceTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class Domain;
class DistGeoState;
class Communicator;
class Timer;
class DistStructureLevelSet;

template<int dim> class DistVMSLESTerm;
template<int dim> class DistDynamicVMSTerm;
template<int dim> class DistDynamicLESTerm;
template<int dim> class DistEdgeGrad;
template<int dim> class DistExtrapolation;
template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<class Scalar, int dim> class DistSVec;
template<class Scalar, int dim> class DistMat;
template<int dim> class DistExactRiemannSolver;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim>
class SpaceOperator {

protected:

  VarFcn *varFcn;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  DistSVec<double,dim> *V;

// Included (MB)
  DistSVec<double,dim> *dU;
  DistSVec<double,dim> *dV;
  DistSVec<double,dim> *dRm;
  int jacobianAA;
  int jacobianSA;
  int reconstructionAA;
  int reconstructionSA;
  FluxFcn **fluxFcnAA;
  RecFcn *recFcnAA;
  FluxFcn **fluxFcnSA;
  RecFcn *recFcnSA;

protected:

  bool locAlloc;

  BcFcn *bcFcn;
  FluxFcn **fluxFcn;
  RecFcn *recFcn;
  //RecFcn *recFcnLS;
  DistNodalGrad<dim, double> *ngrad;
  //DistNodalGrad<dimLS, double> *ngradLS;
  DistNodalGrad<dim, bcomp> *compNodalGrad;

  DistEdgeGrad<dim> *egrad;
  DistExtrapolation<dim> *xpol;
  DistVMSLESTerm<dim> *vms;
  DistDynamicLESTerm<dim> *dles;
  DistDynamicVMSTerm<dim> *dvms;
  FemEquationTerm *fet;
  SmagorinskyLESTerm *smag;
  WaleLESTerm *wale;
  VolumicForceTerm *volForce;

  Domain *domain;

  Timer *timer;
  Communicator *com;

  bool use_modal;
  bool use_complex;
  int order;
  int failsafe;
  int rshift;

// Included (MB)
  IoData *iod;

public:

  SpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  SpaceOperator(const SpaceOperator<dim> &, bool);
  ~SpaceOperator();

  DistNodalGrad<dim, double> *getDistNodalGrad(DistSVec<double,dim> &)  { return ngrad; }
  DistNodalGrad<dim, bcomp> *getDistNodalGrad(DistSVec<bcomp,dim> &)  { return compNodalGrad; }

  int getSpaceOrder() {return order;}

  int **getNodeType() {return domain->getNodeType(); }

  BcFcn *createBcFcn(IoData &);
  FluxFcn **createFluxFcn(IoData &);
  RecFcn *createRecFcn(IoData &);
  //RecFcn *createRecFcnLS(IoData &);
  FemEquationTerm *createFemEquationTerm(IoData &);
  VolumicForceTerm *createVolumicForceTerm(IoData &);
  void setBcFcn(BcFcn *);
  void setFluxFcn(FluxFcn **);
  void setRecFcn(RecFcn *);
  void setFemEquationTerm(FemEquationTerm *);
  void fix(DistSVec<bool,2>&);
  void resetTag();

  FemEquationTerm *getFemEquationTerm() { return fet;}

// Included (MB)
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistTimeState<dim> *, bool=true);
// Included (MB)
  void computeResidual(DistExactRiemannSolver<dim> *,
                       DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistTimeState<dim> *, bool=true);
  //template<int dimLS>
  //void computeResidual(DistSVec<double,3> &, DistVec<double> &,
  //                     DistSVec<double,dim> &, DistSVec<double,dimLS> &,
  //                     DistVec<int> &, DistSVec<double,dim> &,
  //                     DistExactRiemannSolver<dim> *, int it,
  //                     DistSVec<double,dim> * = 0,
  //                     DistSVec<double,dim> * = 0);
// Kevin's FSI with FS Riemann solver
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistSVec<double,dim> &, DistLevelSetStructure *,
                       bool, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int, DistSVec<double,3> *, int it = 0,
		       DistVec<GhostPoint<dim>*> *ghostPoints=0);
// Kevin's FSI with FS Riemann solver (for thin shell problems) 
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistSVec<double,dim> &, DistLevelSetStructure *,
                       bool, DistVec<int> &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int, DistSVec<double,3> *, int it = 0);

  //template<int dimLS>
  //void computeResidualLS(DistSVec<double,3> &, DistVec<double> &,
  //                     DistSVec<double,dimLS> &, DistVec<int> &, DistSVec<double,dim> &,DistSVec<double,dimLS> &);

  void computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                                       DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                                       DistLevelSetStructure *distLSS, DistVec<int> *fluidId = 0);

  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,dim> &U,VarFcn *varFcn,DistLevelSetStructure *distLSS,DistVec<int> &tag);

  void computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X,
                           DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                           DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                           DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                           DistLevelSetStructure *distLSS, DistVec<int> *fluidId =  0);
  void updatePhaseChange(DistSVec<double,dim> &V,
                         DistSVec<double,dim> &U,
                         DistVec<double> *Weights, DistSVec<double,dim> *VWeights,
                         DistLevelSetStructure *distLSS, double* vfar, DistVec<int> *fluidId = 0);

  void computeCellAveragedStructNormal(DistSVec<double,3> &, DistLevelSetStructure *);


  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
 
  double computeRealFluidResidual(DistSVec<double, dim> &, DistSVec<double,dim> &, DistLevelSetStructure &);

  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                                     DistSVec<double,dim> &);
  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                    DistVec<int> &, DistSVec<double,dim> &);


  void computePostOpDVMS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dim> &, DistVec<double> *,
                         DistTimeState<dim> *);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		       DistTimeState<dim> *);

  void getExtrapolationValue(DistSVec<double,dim>&, DistSVec<double,dim>&, DistSVec<double,3>&);
  void applyExtrapolationToSolutionVector(DistSVec<double,dim>&, DistSVec<double,dim>&);

  template<class Scalar, int neq>
  void computeViscousJacobian(DistSVec<double,3> &, DistVec<double> &, DistMat<Scalar,neq> &);

  void applyBCsToSolutionVector(DistSVec<double,dim> &);

  void applyBCsToResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void applyBCsToJacobian(DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void applyBCsToH2Jacobian(DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar>
  void computeH1(DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistMat<Scalar,dim> &);
  template<class Scalar, int neq>
  void computeH2(DistSVec<double,3> &, DistVec<double> &,
		 DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<class Scalar1, class Scalar2>
  void applyH2(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
               DistMat<Scalar1,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2>
  void applyH2T(DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistMat<Scalar1,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar, int neq>
  void printAllMatrix(DistMat<Scalar,neq> &, int);

  void printAllVariable(DistSVec<double,3>&, DistSVec<double,dim> &, int);
  void printVariable(DistSVec<double,dim> &);


  // Included (MB)
  void rstFluxFcn(IoData &);


  // Included (MB)
  void computeDerivativeOfResidual(DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &,
                               DistSVec<double,dim> &, double, DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> * = 0);


  // Included (MB)
  void applyBCsToDerivativeOfResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);


  // Included (MB)
  void rstVarFet(IoData &ioData) 
  {
    if (fet) fet->rstVar(ioData, com);
  }


  // Included (MB)
  bool useModal() 
  {return use_modal;}


  // Included (MB)
  /// \note This function is implemented.
  /// It is called only from MatVecProdFD::evaluateInviscid and
  /// MatVecProdFD::applyInviscid.
  void computeInviscidResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> * = 0, bool=true);


  // Included (MB)
  /// \note This function is implemented.
  /// It is called only from MatVecProdFD::evaluateViscous and
  /// MatVecProdFD::applyViscous.
  void computeViscousResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> * = 0, bool=true);


  // Included (MB)
  /// \note This routine is not implemented.
  //void computeInviscidResidual
  //(
  //  DistSVec<double,3> &, DistVec<double> &,
  //  DistSVec<double,dim> &, DistVec<double> &,
  //  DistSVec<double,dim> &, bool=true
  //);


  // Included (MB)
  /// \note This routine is not implemented.
  //void computeViscousResidual
  //(
  //  DistSVec<double,3> &, DistVec<double> &,
  //  DistSVec<double,dim> &, DistVec<double> &,
  //  DistSVec<double,dim> &, bool=true
  //);


  // Included (MB)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeGradP(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);


  // Included (MB)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeDerivativeOfGradP
  (
    DistSVec<double,3> &X, DistSVec<double,3> &dX,
    DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
    DistSVec<double,dim> &U, DistSVec<double,dim> &dU
  );


  // Included (MB)
  /// \note This routine is redundant. It is never called.
  //void computeDerivativeOfGradP
  //(
  //  DistSVec<double,3> &X, DistSVec<double,3> &dX,
  //  DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
  //  DistSVec<double,dim> &U
  //);


  void computeForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X,
                        double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji);
};

//------------------------------------------------------------------------------
// derived class of SpaceOperator to run multi-phase flow problems
// the added functions allow to advance the Euler equations when
// different EOS are considered and to advance the level-set 
// advection equation(s).
template<int dim, int dimLS>
class MultiPhaseSpaceOperator : public SpaceOperator<dim> {

protected:

  RecFcn *recFcnLS;
  DistNodalGrad<dimLS, double> *ngradLS; 

  RecFcn *createRecFcnLS(IoData &);

public:

  MultiPhaseSpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  MultiPhaseSpaceOperator(const MultiPhaseSpaceOperator<dim,dimLS> &, bool);
  ~MultiPhaseSpaceOperator();


  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                       FluidSelector &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int it,
                       DistSVec<double,dim> * = 0,
                       DistSVec<double,dim> * = 0);

  void computeResidualLS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dimLS> &, DistVec<int> &, 
                         DistSVec<double,dim> &,DistSVec<double,dimLS> &);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
                       FluidSelector &, DistExactRiemannSolver<dim> *);

  template<class Scalar>
  void computeJacobianLS(DistSVec<double,3> &X,DistSVec<double,dim> &V, DistVec<double> &ctrlVol,
			 DistSVec<double,dimLS> &Phi,DistMat<Scalar,dimLS> &A,DistVec<int> &fluidId);

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SpaceOperator.C>
#endif

#endif
