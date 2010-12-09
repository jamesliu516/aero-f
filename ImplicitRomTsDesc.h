#ifndef _IMPLICIT_ROM_TS_DESC_H_
#define _IMPLICIT_ROM_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>

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


  int maxItsNewton;
  double epsNewton;

  MatVecProdFD<dim,dim> *mvpfd;

  DistSVec<bool,2> *tag;

  VecSet<DistSVec<double, dim> > pod;

  FullM jac;
  
  int nPod;

  DistSVec<double, dim> F;	// residual
  VecSet<DistSVec<double, dim> > AJ; // Action of Jacobian (AJ) on reduced-order basis

  Vec<double> dUrom;
  Vec<double> UromTotal;

	double target, res0;	// for Newton convergence

  virtual void computeAJ(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time
  virtual void computeFullResidual(int, DistSVec<double, dim> &);  

  virtual void saveAJsol() {};	// only implemented for PG
	virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop) = 0;
	// each ROM has a different way of solving the Newton system
  virtual void updateGlobalTimeSteps(int _it) {};	// each ROM has a different way of solving the Newton system
  int solveLinearSystem(int, Vec<double> &, Vec<double> &);
  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double); 
  double meritFunctionDeriv(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &);
  double lineSearch(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double zoom(double, double, double, double, double, double, double, double, double, DistSVec<double,dim>,DistSVec<double,dim>, DistSVec<double,dim>, int);
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();
  void projectVector(VecSet<DistSVec<double, dim> >&, DistSVec<double, dim> &, Vec<double> &);
  void expandVector(Vec<double> &, DistSVec<double, dim> &);
	virtual void writeStateRomToDisk(int it, double cpu);

public:
  
  ImplicitRomTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitRomTsDesc();

  int solveNonLinearSystem(DistSVec<double, dim> &, int _it);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomTsDesc.C>
#endif

#endif
