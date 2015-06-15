#ifndef _IMPLICIT_ROM_TS_DESC_H_
#define _IMPLICIT_ROM_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>
#include <NonlinearRom.h>
#include <NonlinearRomOnlineII.h>
#include <NonlinearRomOnlineIII.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;


//------------------------------------------------------------------------------

template<int dim>
class ImplicitRomTsDesc : public TsDesc<dim> {

protected:

  IoData *ioData;

  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton;
  double epsAbsIncNewton;

  MatVecProd<dim,dim> *mvp;

  DistSVec<bool,2> *tag;

  VecSet<DistSVec<double, dim> > pod;

  NonlinearRom<dim>* rom;

  int currentCluster;

  int basisUpdateFreq;  
  int tryAllFreq;

  FullM jac;
  
  int nPod;
  
  bool unsteady;
  bool systemApprox;
  bool useIncrements;
  bool tryingAllClusters;

  DistSVec<double, dim> F;	// residual
  VecSet<DistSVec<double, dim> > AJ; // Action of Jacobian (AJ) on reduced-order basis

  DistSVec<double, dim>* weightVec;     // weighting vector for least squares system
  //DistSVec<double, dim>* weightURef;	  // reference state used for weighting least squares system
  //DistSVec<double, dim>* weightFRef;	  // reference residual used for weighting least squares system

  DistVec<double>* farFieldMask;          // one for far field nodes, zero otherwise
  DistVec<double>* farFieldNeighborsMask; // one for neighbors of ff nodes, zero otherwise
  DistVec<double>* wallMask;              // one for wall nodes, zero otherwise
  DistVec<double>* wallNeighborsMask;     // one for neighbors of wall nodes, zero otherwise
  
  double interiorWeight;
  double ffWeight;
  double wallWeight;
  double bcWeightGrowthFactor;
  double levenbergMarquardtWeight;
  bool wallUp;
  bool wallDown;
  double wallWeightGrowthFactor;
  bool ffUp;
  bool ffDown;
  double ffWeightGrowthFactor;

  int allowBCWeightDecrease;
  int adjustInteriorWeight;

  double regThresh;
  double regWeight;

  // backtracking line search
  double rho;
  double c1;
  int maxItsLS; 

  std::vector<double> interpWeightsForMultiIC;

  Vec<double> dUromNewtonIt;    // set to zero before each newton iteration
  Vec<double> dUromTimeIt;      // set to zero before each time iteration
  Vec<double> dUromCurrentROB;  // set to zero after each cluster switch
  Vec<double> UromCurrentROB;   // for projection only: initialized at each cluster switch

  // dUromAccum, 

  //Vec<double> *dUnormAccum;	// accumulated contributions

  double target, res0;	// for Newton convergence

  virtual void computeAJ(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL);
  virtual void computeRedHessianSums(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time 
  virtual void computeFullResidual(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL, bool includeHomotopy = true); 
  virtual void saveNewtonSystemVectors(const int _it) {};	// only implemented for PG/Galerkin
  void saveNewtonSystemVectorsAction(const int);	// implementation for PG/Galerkin
	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &, const int &totalTimeSteps = 0) = 0;
	// each ROM has a different way of solving the Newton system
  int solveLinearSystem(int, Vec<double> &, Vec<double> &);
  virtual double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double); 
  double meritFunctionDeriv(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);
  double lineSearch(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double lineSearchBacktrack(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double zoom(double, double, double, double, double, double, double, double, double, DistSVec<double,dim>,DistSVec<double,dim>, DistSVec<double,dim>, int);
  int checkFailSafe(DistSVec<double,dim>&);
  void printBCWeightingInfo(bool);
  void resetFixesTag();
  void projectVector(VecSet<DistSVec<double, dim> >&, DistSVec<double, dim> &, Vec<double> &);
  void expandVector(Vec<double> &, DistSVec<double, dim> &);

  void loadCluster(int, bool, DistSVec<double, dim> &);

  virtual void updateLeastSquaresWeightingVector();
	//void savedUnormAccum();
	//virtual void writeStateRomToDisk(int it, double cpu);
	virtual void postProStep(DistSVec<double,dim> &, int) {};	// by default, do not do post processing
	virtual bool breakloop1(const bool);
	virtual bool breakloop2(const bool);

  virtual void setReferenceResidual() {};
  virtual void setProblemSize(DistSVec<double, dim> &) {};
  virtual void deleteRestrictedQuantities() {};

	double *projVectorTmp; // temporary vector for projectVector

  bool updateFreq;
  bool clusterSwitch;
  bool updatePerformed;

  DistSVec<double, dim>* Uinit;  // initial condition of the steady state simulation, 
                                 // stored to recalculate reference residual after 
                                 // cluster switch (new sampled mesh for GNAT) or after
                                 // changing the residual weighting 

  DistSVec<double, dim>* Uprev;  // solution at the beginning of the previous time step (needed for model II incremental bases) 

  void tryAllClusters(DistSVec<double, dim>&, const int totalTimeSteps, int*);

  void printRomResiduals(DistSVec<double, dim> &U);
  FILE *residualsFile;

protected:
  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

public:
  
  ImplicitRomTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitRomTsDesc();

  int solveNonLinearSystem(DistSVec<double, dim> &, const int _it);
  void rstVarImplicitRomTsDesc(IoData &);
  void checkLocalRomStatus(DistSVec<double, dim> &, const int);
  void setInterpWeightsForMultiIC(std::vector<double> vec) {interpWeightsForMultiIC = vec;}
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomTsDesc.C>
#endif

#endif
