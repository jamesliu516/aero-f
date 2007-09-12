#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include <IoData.h>
#include <complex.h>
typedef complex<double> bcomp;

class VarFcn;
class BcFcn;
class RecFcn;
class FluxFcn;
class FemEquationTerm;
class VolumicForceTerm;
class SmagorinskyLESTerm;
class Domain;
class DistGeoState;
class Communicator;
class Timer;

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

private:

  VarFcn *varFcn;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  DistSVec<double,dim> *V;

private:

  bool locAlloc;

  BcFcn *bcFcn;
  FluxFcn **fluxFcn;
  RecFcn *recFcn;
  RecFcn *recFcnLS;
  DistNodalGrad<dim, double> *ngrad;
  DistNodalGrad<1, double> *ngradLS;
  DistNodalGrad<dim, bcomp> *compNodalGrad;

  DistEdgeGrad<dim> *egrad;
  DistExtrapolation<dim> *xpol;
  DistVMSLESTerm<dim> *vms;
  DistDynamicLESTerm<dim> *dles;
  DynamicLESTerm *dlest;
  DistDynamicVMSTerm<dim> *dvms;
  FemEquationTerm *fet;
  SmagorinskyLESTerm *smag;
  VolumicForceTerm *volForce;

  Domain *domain;

  Timer *timer;
  Communicator *com;

  bool use_modal;
  bool use_complex;
  int order; 
  int failsafe;
  int rshift;

public:

  SpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *, 
		Domain *, DistSVec<double,dim> * = 0);
  SpaceOperator(const SpaceOperator<dim> &, bool);
  ~SpaceOperator();

  DistNodalGrad<dim, double> *getDistNodalGrad(DistSVec<double,dim> &)  { return ngrad; }
  //DistNodalGrad<1, double> *getDistNodalGrad(DistVec<double> &)  { return ngradLS; }
  DistNodalGrad<dim, bcomp> *getDistNodalGrad(DistSVec<bcomp,dim> &)  { return compNodalGrad; }

  int getSpaceOrder() {return order;}

  int **getNodeType() {return domain->getNodeType(); }

  BcFcn *createBcFcn(IoData &);
  FluxFcn **createFluxFcn(IoData &);
  RecFcn *createRecFcn(IoData &);
  RecFcn *createRecFcnLS(IoData &);
  FemEquationTerm *createFemEquationTerm(IoData &);
  VolumicForceTerm *createVolumicForceTerm(IoData &);
  void setBcFcn(BcFcn *);
  void setFluxFcn(FluxFcn **);
  void setRecFcn(RecFcn *);
  void setFemEquationTerm(FemEquationTerm *);
  void fix(DistSVec<bool,2>&);
  void resetTag();

  FemEquationTerm *getFemEquationTerm() { return fet;}

  void storeGhost(DistSVec<double,dim> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistVec<double> &);

  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistTimeState<dim> *);
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistVec<double> &,
                       DistSVec<double,dim> &, 
                       DistExactRiemannSolver<dim> *, int it = 0);
  void computeResidualLS(DistSVec<double,3> &, DistVec<double> &,
                       DistVec<double> &, DistSVec<double,dim> &,DistVec<double> &);

  void storePreviousPrimitive(DistSVec<double,dim> &U, DistSVec<double,dim> &Vg,
                         DistVec<double> &Phi,
                         DistSVec<double,dim> *Vgf, DistVec<double> *weight);
  void updatePhaseChange(DistSVec<double,dim> &Vg, DistSVec<double,dim> &U,
                         DistVec<double> &Phi, DistVec<double> &Phin,
                         DistSVec<double,dim> *Vgf, DistVec<double> *weight,
                         DistExactRiemannSolver<dim> *);

  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                                     DistSVec<double,dim> &);
  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                    DistVec<double> &, DistSVec<double,dim> &);


  void computePostOpDVMS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dim> &, DistVec<double> *,
                         DistTimeState<dim> *);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		       DistTimeState<dim> *);
  
  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
                       DistVec<double> &, DistExactRiemannSolver<dim> *);

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
  template<class Scalar>
  void computeH2(DistSVec<double,3> &, DistVec<double> &,
		 DistSVec<double,dim> &, DistMat<Scalar,dim> &, 
		 DistSVec<double,dim> &, DistSVec<double,dim> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<class Scalar>
  void computeH2LS(DistSVec<double,3> &, DistVec<double> &,
                   DistVec<double> &, DistSVec<double,dim> &, 
                   DistMat<Scalar,1> &);

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

  template<class Scalar1, class Scalar2>
  void applyH2LS(DistSVec<double,3> &, DistVec<double> &ctrlVol,
                 DistVec<double> &, DistMat<Scalar1,1> &,
                 DistSVec<double,1> &, DistSVec<double,1> &,
                 DistSVec<double,1> &, DistSVec<double,1> &,
                 DistVec<Scalar2> &, DistVec<Scalar2> &);

  template<class Scalar, int neq>
  void printAllMatrix(DistMat<Scalar,neq> &, int);

  void printAllVariable(DistSVec<double,3>&, DistSVec<double,dim> &, int);
  void printVariable(DistSVec<double,dim> &);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SpaceOperator.C>
#endif

#endif
