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
template<int dim, int neq> class MatVecProdFD;


//------------------------------------------------------------------------------

template<int dim>
class ImplicitRomTsDesc : public TsDesc<dim> {

protected:

	IoData *ioData;

  int maxItsNewton;
  double epsNewton;

  MatVecProdFD<dim,dim> *mvpfd;

  DistSVec<bool,2> *tag;

  VecSet<DistSVec<double, dim> > pod;

  NonlinearRom<dim>* rom;

  int currentCluster;

  int basisUpdateFreq;  

  FullM jac;
  
  int nPod;

  DistSVec<double, dim> F;	// residual
  VecSet<DistSVec<double, dim> > AJ; // Action of Jacobian (AJ) on reduced-order basis

  Vec<double> dUromNewtonIt;    // set to zero before each newton iteration
  Vec<double> dUromTimeIt;      // set to zero before each time iteration
  Vec<double> dUromCurrentROB;  // set to zero after each cluster switch

  // dUromAccum, 

  //Vec<double> *dUnormAccum;	// accumulated contributions

	double target, res0;	// for Newton convergence

  virtual void computeAJ(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time
  virtual void computeFullResidual(int, DistSVec<double, dim> &, DistSVec<double, dim> *R = NULL);  

  virtual void saveNewtonSystemVectors(const int _it) {};	// only implemented for PG/Galerkin
  void saveNewtonSystemVectorsAction(const int);	// implementation for PG/Galerkin
	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &, const int &totalTimeSteps = 0) = 0;
	// each ROM has a different way of solving the Newton system
  virtual void updateGlobalTimeSteps(const int _it) {};	// broyden needs to know global time steps
  int solveLinearSystem(int, Vec<double> &, Vec<double> &);
  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double); 
  double meritFunctionDeriv(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &);
  double lineSearch(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double zoom(double, double, double, double, double, double, double, double, double, DistSVec<double,dim>,DistSVec<double,dim>, DistSVec<double,dim>, int);
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();
  void projectVector(VecSet<DistSVec<double, dim> >&, DistSVec<double, dim> &, Vec<double> &);
  void expandVector(Vec<double> &, DistSVec<double, dim> &);
  virtual void checkLocalRomStatus(DistSVec<double, dim> &, const int);
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

protected:
  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

public:
  
  ImplicitRomTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitRomTsDesc();

  int solveNonLinearSystem(DistSVec<double, dim> &, const int _it);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomTsDesc.C>
#endif

#endif
