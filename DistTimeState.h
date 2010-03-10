#ifndef _DIST_TIME_STATE_H_
#define _DIST_TIME_STATE_H_

#include <TimeData.h>
#include <DistVector.h>
#include <DistMacroCell.h>
#include <LowMachPrec.h>

class VarFcn;
class Domain;
class DistGeoState;
class FemEquationTerm;

template<int dim> class TimeState;
template<class Scalar, int dim> class DistMat;
template<int dim> class SpaceOperator;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------------

template<int dim>
class DistTimeState {

private:

  VarFcn *varFcn;
  FemEquationTerm *fet;
  DistSVec<double,dim> *V;
  DistSVec<double,dim> *Vn;
  DistSVec<double,dim> *VnBar;
  DistSVec<double,dim> *QBar;

  bool locAlloc;
  int numLocSub;

  double gam;
  double pstiff;

  TimeLowMachPrec    tprec;
  SpatialLowMachPrec sprec; //only for computation of irey
  /*bool prec;
  double beta;
  double mach;
  double cmach;
  double k1;
  double betav;
*/
  double viscousCst;

  TimeData *data;

  DistVec<double> *dt;			//actual   time stepping
  DistVec<double> *idti;		//inverse inviscid time stepping
  DistVec<double> *idtv;		//inverse viscous  time stepping
  DistVec<double> *irey;
  DistSVec<double,dim> *Un;
  DistSVec<double,dim> *Unm1;
  DistSVec<double,dim> *Unm2;
  DistSVec<double,dim> *UnBar;
  DistSVec<double,dim> *Unm1Bar;
  DistSVec<double,dim> *Unm2Bar;
  DistSVec<double,dim> *Rn;

  Domain *domain;

  TimeState<dim> **subTimeState;

// Included (MB)
  DistVec<double> *dIrey;
  DistVec<double> *dIdti;
  DistVec<double> *dIdtv;

public:

  DistTimeState(IoData &, SpaceOperator<dim> *, VarFcn *, Domain *, DistSVec<double,dim> * = 0);
  DistTimeState(const DistTimeState<dim> &, bool, IoData &);
  ~DistTimeState();

  TimeState<dim> &operator() (int i) const { return *subTimeState[i]; }

  void setResidual(DistSVec<double,dim> *rn) { if (Rn != 0) delete Rn; Rn = rn; }

  void setGlobalTimeStep (double t) { *dt = t; }

  void setup(const char *, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,dim> &, 
             IoData * = 0, DistVec<int> * = 0); //iod and fluidId needed only for multi-phase flow
  void setupUFluidIdInitialConditions(FluidModelData &fm, DistVec<int> &fluidId, SphereData &ic, int myId);
//------ to be removed (or moved) -----
  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &Ufar,
             DistSVec<double,dim> &U, IoData &iod);
  void setupUVolumesInitialConditions(IoData &iod);
  void setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X);
  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &Ufar,
             DistSVec<double,dim> &U, DistVec<int> &nodeTag, IoData &iod);
  void setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X, DistVec<int> &nodeTag);
  void setup(const char *name, DistSVec<double,dim> &Ufar, double *Ub, DistSVec<double,3> &X,
             DistSVec<double,1> &Phi, DistSVec<double,dim> &U, IoData &iod);
//-------------------------------------
  void update(DistSVec<double,dim> &);
  void update(DistSVec<double,dim> &Q, DistVec<int> &fluidId, DistVec<int> *fluidIdnm1, 
              //DistSVec<double,dim> *Vgf, DistVec<double> *Vgfweight,
              DistExactRiemannSolver<dim> *riemann);

  void writeToDisk(char *);

  double computeTimeStep(double, double*, int*, DistGeoState &, 
			 DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);
  double computeTimeStep(double, double*, int*, DistGeoState &,
                         DistVec<double> &, DistSVec<double,dim> &, DistVec<int> &);

  void computeCoefficients(double);

  void add_dAW_dt(int, DistGeoState &, DistVec<double> &, 
			  DistSVec<double,dim> &, DistSVec<double,dim> &);
  void add_dAW_dtLS(int, DistGeoState &, DistVec<double> &, 
			 DistSVec<double,1> &, DistSVec<double,1> &, DistSVec<double,1> &, 
			 DistSVec<double,1> &, DistSVec<double,1> &);

  template<class Scalar, int neq>
  void addToJacobian(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianNoPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianGasPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianLiquidPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &, Scalar);

  template<class Scalar>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  template<class Scalar>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &, Scalar);

  template<class Scalar>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &, Scalar, double);
  
  template<class Scalar>
  void addToH2LS(DistVec<double> &, DistVec<double> &, DistMat<Scalar,1> &);

  template<class Scalar>
  void addToH2Minus(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  void computeBar(bool, DistMacroCellSet *, DistGeoState &, int);
                                                                                                                          
  void get_dW_dt(bool, DistGeoState &, DistVec<double> &, DistSVec<double,dim> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                 DistMacroCellSet *, DistSVec<double,1> **, int);

  DistVec<double>* getInvReynolds(){ return irey; }
                                                                                                                          
  void multiplyByTimeStep(DistSVec<double,dim>&);
  void multiplyByTimeStep(DistSVec<double,1>&);
  void multiplyByPreconditioner(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerPerfectGas(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerLiquid(DistSVec<double,dim> &, DistSVec<double,dim>&);

  TimeData &getData() { return *data; }
  DistSVec<double,dim> &getUn() const { return *Un; }

  double getTime()  { return data->dt_n; }

// Included (MB)
  TimeData* getDataOpt() { return data; }

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  void rstVar(IoData &);

  DistVec<double>* getDerivativeOfInvReynolds(DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistTimeState.C>
#endif

#endif
