#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include <IoData.h>
#include <GhostPoint.h>
#include <RestrictionMapping.h>
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

  enum DescriptorCase {
    DESCRIPTOR, HYBRID, NONDESCRIPTOR
  };
  DescriptorCase descriptorCase;



public:

  SpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  SpaceOperator(const SpaceOperator<dim> &, bool);
  ~SpaceOperator();

  DistNodalGrad<dim, double> *getDistNodalGrad(DistSVec<double,dim> &)  { return ngrad; }
  DistNodalGrad<dim, bcomp> *getDistNodalGrad(DistSVec<bcomp,dim> &)  { return compNodalGrad; }

  void updateFixes() { ngrad->updateFixes(); }

  int getSpaceOrder() {return order;}

  int **getNodeType() {return domain->getNodeType(); }

  BcFcn *createBcFcn(IoData &);
  FluxFcn **createFluxFcn(IoData &);
  RecFcn *createRecFcn(IoData &);
  //RecFcn *createRecFcnLS(IoData &);
  FemEquationTerm *createFemEquationTerm(IoData &);
  VolumicForceTerm *createVolumicForceTerm(IoData &);
  VarFcn* getVarFcn() { return varFcn; }
  void setBcFcn(BcFcn *);
  void setFluxFcn(FluxFcn **);
  void setRecFcn(RecFcn *);
  void setFemEquationTerm(FemEquationTerm *);
  void fix(DistSVec<bool,2>&);
  void resetTag();

  FemEquationTerm *getFemEquationTerm() { return fet;}
  void conservativeToPrimitive(DistSVec<double,dim> &U) {varFcn->conservativeToPrimitive(U, *V);}

// Included (MB)
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistTimeState<dim> *, bool=true);
	void computeResidualRestrict(DistSVec<double,3> &, DistVec<double> &,
			DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> *,
			RestrictionMapping<dim> &, bool=true);
// Included (MB)
  void computeResidual(DistExactRiemannSolver<dim> *,
											 DistSVec<double,3> &, DistVec<double> &,
											 DistSVec<double,dim> &, DistSVec<double,dim> &,
											 DistTimeState<dim> *, bool=true);

  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistSVec<double,dim> &, DistLevelSetStructure *,
                       bool, bool, DistVec<int> &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int, DistSVec<double,3> *, int it = 0,
                       DistVec<GhostPoint<dim>*> *ghostPoints = 0);

  void computeResidual(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, 
		  			   DistSVec<double,dim> &, DistSVec<double,dim> &, 
					   DistVec<int> &, DistVec<int> &, DistLevelSetStructure *,
                       bool, bool, DistVec<int> &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int, DistSVec<double,3> *, 
					   double, double, int it = 0, DistVec<GhostPoint<dim>*> *ghostPoints = 0);

  void updateSweptNodes(DistSVec<double,3> &X, int &phaseChangeChoice, int &phaseChangeAlg,
                        DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                        DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                        DistLevelSetStructure *distLSS, double *vfar, DistVec<int> *fluidId = 0);

  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,3> &X,DistSVec<double,dim> &U,VarFcn *varFcn,DistLevelSetStructure *distLSS,bool linFSI,DistVec<int> &tag);
  
  template <int neq>
  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,neq> &U,VarFcn *varFcn,DistLevelSetStructure *distLSS,DistVec<int> &tag) {
    fprintf(stderr,"PopulateGhostPoints<%d> not implemented!\n",neq);
    exit(-1);
  }

  void computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X,
                           DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                           DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                           DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                           DistLevelSetStructure *distLSS, DistVec<int> *fluidId =  0);

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

  template<class Scalar, int neq>
  void computeJacobian(DistExactRiemannSolver<dim> *riemann, DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		       DistTimeState<dim> *);

  template <class Scalar,int neq>
  void computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                       DistSVec<double,dim> &U,
                       DistLevelSetStructure *distLSS,
                       DistVec<int> &fluidId, 
                       DistExactRiemannSolver<dim> *riemann, 
                       int Nriemann, DistSVec<double,3> *Nsbar,
                       DistVec<GhostPoint<dim>*> *ghostPoints,
                       DistMat<Scalar,neq>& A,
                       DistTimeState<dim>*);
  
  void getExtrapolationValue(DistSVec<double,dim>&, DistSVec<double,dim>&, DistSVec<double,3>&);
  void applyExtrapolationToSolutionVector(DistSVec<double,dim>&, DistSVec<double,dim>&);

  template<class Scalar, int neq>
  void computeViscousJacobian(DistSVec<double,3> &, DistVec<double> &, DistMat<Scalar,neq> &);

  void applyBCsToSolutionVector(DistSVec<double,dim> &,DistLevelSetStructure *distLSS=0);

  void applyBCsToResidual(DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

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

  template<class Scalar>
  void applyBCsToH2Jacobian(DistSVec<double,dim> &, DistMat<Scalar,dim> &);

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


  void computeForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
                        double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
			DistVec<GhostPoint<dim>*> *ghostPoints = 0, PostFcn *postFcn = 0,DistVec<int>* fid = 0);
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

  struct {

    DistVec<int>* counts[2];
    DistSVec<double,dim>* vals[2];
    DistNodalGrad<dim,double>* cutgrads[2];    
  } higherOrderData;

public:

  MultiPhaseSpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  MultiPhaseSpaceOperator(const MultiPhaseSpaceOperator<dim,dimLS> &, bool);
  ~MultiPhaseSpaceOperator();

  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                       FluidSelector &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int it);
  void computeResidual(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, 
                       DistSVec<double,dim> &, DistLevelSetStructure *, bool, bool, DistExactRiemannSolver<dim> *, 
                       int, DistSVec<double,3> *, DistSVec<double,dimLS> &, FluidSelector &, 
                       DistSVec<double,dim> &, int, DistVec<GhostPoint<dim>*> *);

  void computeResidualLS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dimLS> &, DistVec<int> &, 
                         DistSVec<double,dim> &,DistSVec<double,dimLS> &, DistLevelSetStructure* =0, bool = true,
			 int method = 0);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
                       FluidSelector &, DistExactRiemannSolver<dim> *,DistTimeState<dim>*);

  template<class Scalar, int neq>
  void computeJacobian(DistExactRiemannSolver<dim>* riemann,
                       DistSVec<double,3>& X, DistSVec<double,dim>& U,DistVec<double>& ctrlVol,
                       DistLevelSetStructure *distLSS,
                       int Nriemann, DistSVec<double,3>* Nsbar,
                       FluidSelector &fluidSelector,
                       DistMat<Scalar,neq>& A,DistTimeState<dim>* timeState);

  template<class Scalar>
  void computeJacobianLS(DistSVec<double,3> &X,DistSVec<double,dim> &V, DistVec<double> &ctrlVol,
			 DistSVec<double,dimLS> &Phi,DistMat<Scalar,dimLS> &A,DistVec<int> &fluidId,DistLevelSetStructure* distLSS,
			 int method);

  // for phase-change update
  void extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV);
  void extrapolatePhiV2(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV);
  void updateSweptNodes(DistSVec<double,3> &X, int &phaseChangeChoice, DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                        DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                        DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                        DistLevelSetStructure *distLSS, double *vfar, bool updateWithCracking,
                        DistVec<int> *fluidId0, DistVec<int> *fluidId);
  void resetFirstLayerLevelSetFS(DistSVec<double,dimLS> &PhiV, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, 
                                 DistSVec<bool,2> &Tag);

  void findCutCells(DistSVec<double,dimLS>& phi,
		    DistVec<int>& status,
		    DistVec<int>& fluidId,
		    DistSVec<double,dim> &V,
		    DistSVec<double,3> &X);

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SpaceOperator.C>
#endif

#endif
