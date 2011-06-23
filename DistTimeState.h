#ifndef _DIST_TIME_STATE_H_
#define _DIST_TIME_STATE_H_

#include <TimeData.h>
#include <DistVector.h>
#include <DistMacroCell.h>
#include <LowMachPrec.h>
#include <LevelSet/LevelSetStructure.h>

struct FluidModelData;
struct InitialConditions;

class VarFcn;
class Domain;
class DistGeoState;
class FemEquationTerm;
class DistLevelSetStructure;

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

  bool isGFMPAR;
  int fvmers_3pbdf ;

private:
  void computeInitialState(InitialConditions &ic, FluidModelData &fm, double UU[dim]);
public:

  DistTimeState(IoData &, SpaceOperator<dim> *, VarFcn *, Domain *, DistSVec<double,dim> * = 0);
  DistTimeState(const DistTimeState<dim> &, bool, IoData &);
  ~DistTimeState();

  TimeState<dim> &operator() (int i) const { return *subTimeState[i]; }

  void setResidual(DistSVec<double,dim> *rn) { if (Rn != 0) delete Rn; Rn = rn; }

  void setGlobalTimeStep (double t) { *dt = t; }

  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &Ufar,
             DistSVec<double,dim> &U, IoData &iod, DistVec<int> *fluidId = 0); 
  void setupUVolumesInitialConditions(IoData &iod);
  void setupUOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X);
  void setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X);
  void setupUFluidIdInitialConditions(IoData &iod, DistVec<int> &fluidId);
  void update(DistSVec<double,dim> &,bool increasingPressure = false);
  void update(DistSVec<double,dim> &Q,  DistSVec<double,dim> &Qtilde,DistVec<int> &fluidId, DistVec<int> *fluidIdnm1, 
              DistExactRiemannSolver<dim> *riemann,class DistLevelSetStructure* = 0, bool increasingPressure = false);

  void writeToDisk(char *);

  double computeTimeStep(double, double*, int*, DistGeoState &, 
			 DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);
  double computeTimeStep(double, double*, int*, DistGeoState &,
                         DistVec<double> &, DistSVec<double,dim> &, DistVec<int> &, DistVec<double>* = NULL);

  void computeCoefficients(double);

  void add_dAW_dt(int, DistGeoState &, DistVec<double> &, 
		  DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);
  template<int dimLS>
  void add_dAW_dtLS(int, DistGeoState &, DistVec<double> &, 
			 DistSVec<double,dimLS> &, DistSVec<double,dimLS> &, DistSVec<double,dimLS> &, 
			 DistSVec<double,dimLS> &, DistSVec<double,dimLS> &,bool requireSpecialBDF = false);

  template<class Scalar, int neq>
  void addToJacobian(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianNoPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianLS(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &,bool);

  template<class Scalar, int neq>
  void addToJacobianGasPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianLiquidPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &, Scalar, double);
  
  template<class Scalar, int neq>
  void addToH2Minus(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  void computeBar(bool, DistMacroCellSet *, DistGeoState &, int);
                                                                                                                          
  void get_dW_dt(bool, DistGeoState &, DistVec<double> &, DistSVec<double,dim> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                 DistMacroCellSet *, DistSVec<double,1> **, int);

  DistVec<double>* getInvReynolds(){ return irey; }
                                                                                                                          
  void multiplyByTimeStep(DistSVec<double,dim>&);
  template<int dimLS>
  void multiplyByTimeStep(DistSVec<double,dimLS>&);
  void multiplyByPreconditioner(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerPerfectGas(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerLiquid(DistSVec<double,dim> &, DistSVec<double,dim>&);

  TimeData &getData() { return *data; }
  DistSVec<double,dim> &getUn() const { return *Un; }
  DistSVec<double,dim> &getUnm1() const { return *Unm1; }

  inline bool existsNm1() const { return data->exist_nm1; }
  inline bool useNm1() const { return data->use_nm1; }

  double getTime()  { return data->dt_n; }

  void rstVar(IoData &);

  DistVec<double>* getDerivativeOfInvReynolds(DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistTimeState.C>
#endif

#endif
