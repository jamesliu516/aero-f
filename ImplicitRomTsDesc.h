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
  int JacSkipNewton;

  MatVecProdFD<dim,dim> *mvpfd;

  DistSVec<bool,2> *tag;

  VecSet<DistSVec<double, dim> > pod;

  FullM jac;
  
  int nPod;

  VecSet<DistSVec<double, dim> > AJ; // Action of Jacobian (AJ) on reduced-order basis

  Vec<double> dUrom;

  Vec<double> Fromold;

  int RomSolver;

  // Gappy Pod
  int nIntNodes, dimInterpMat;
  int *globalSubSet;
  int *locNodeSet;
  FullM interpMat1;
  FullM interpMat2;
  int myNNodeInt;
  int *myLocalSubSet;
  int *myInterpNodes;

public:
  
  ImplicitRomTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitRomTsDesc();

  void computeJacobian(int, DistSVec<double, dim> &, DistSVec<double, dim> &, VecSet<DistSVec<double, dim> > &);
  void computeJacobianGappy(int, DistSVec<double, dim> &, DistSVec<double, dim> &, Vec<double> &, VecSet<DistSVec<double, dim> > &);
  void rowPartition(int &, int &, int, int sym = 0);
  void computeRestrictInfo();
  void setOperators(DistSVec<double,dim> &) {};
  int solveLinearSystem(int, Vec<double> &, Vec<double> &);
  int solveNonLinearSystem(DistSVec<double, dim> &, int _it);
  void computeFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &);  
  void recomputeFunction(DistSVec<double, dim> &, Vec<double> &);
  void broydenUpdate(Vec<double> &, Vec<double> &); 
  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double); 
  double meritFunctionDeriv(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &);
  double lineSearch(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double zoom(double, double, double, double, double, double, double, double, double, DistSVec<double,dim>,DistSVec<double,dim>, DistSVec<double,dim>, int);
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();
  void projectVector(VecSet<DistSVec<double, dim> >&, DistSVec<double, dim> &, Vec<double> &);
  void expandVector(Vec<double> &, DistSVec<double, dim> &);

  //int getMaxItsNewton() const { return maxItsNewton; }
  //double getEpsNewton() const { return epsNewton; }

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomTsDesc.C>
#endif

#endif
