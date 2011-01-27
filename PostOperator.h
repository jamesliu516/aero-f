#ifndef _POST_OPERATOR_H_
#define _POST_OPERATOR_H_

#include <PostFcn.h>
#include <VectorSet.h>
#include <map>
using std::map;


class IoData;
class VarFcn;
class SubDomain;
class Domain;
class DistGeoState;
class Communicator;
class SmagorinskyLESTerm;
class WaleLESTerm;
class DynamicLESTerm;

struct Vec3D;

template<int dim> class DistBcData;
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistVMSLESTerm;
template<int dim> class DistDynamicLESTerm;
template<int dim> class DistDynamicVMSTerm;
template<int dim> class SpaceOperator;
template<int dim> class DistTimeState;
template<int dim> class DistExactRiemannSolver;


template<int dim>
class ForceGenerator {
  public:
    virtual void getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]) = 0;
};

//------------------------------------------------------------------------------

template<int dim>
class PostOperator {

  VarFcn *varFcn;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  DistSVec<double,dim> *V;

  int numLocSub;
  SubDomain **subDomain;
  Domain *domain;
  Communicator *com;
  int numSurf;
  int numSurfHF;
  map<int,int> surfOutMap;
  map<int,int> surfOutMapHF;
  map<int,int> surfComputeMap; //AS far as I can figure out this map is never used
  ForceGenerator<dim> *forceGen;
  
 // Coefficients to Compute nodal force transfer
  double nodalForceWeights[2];

// Included (MB)
  DistSVec<double,dim> *dV;

private:

  double threshold;
  double refLengthSq;
  double pressInfty;
  DistSVec<double,2>* tmp2;
  SmagorinskyLESTerm *smag;
  WaleLESTerm *wale;  
  DistVMSLESTerm<dim> *vms;
  DistDynamicLESTerm<dim> *dles;
  DistDynamicVMSTerm<dim> *dvms;
  SpaceOperator<dim> *spaceOp;
  CommPattern<double>* vec2Pat;
  PostFcn *postFcn;
  DistVec<double> *mutOmu;
  DistVec<double> *Cs;
  DistVec<double> *CsDvms;
  DistVec<double> *CsDles;

public:

  PostOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *, 
	       Domain *, DistSVec<double,dim> * = 0);
  ~PostOperator();

  void computeNodalForce(DistSVec<double,3> &, DistSVec<double,dim> &, 
			 DistVec<double> &, DistSVec<double,3> &,
			 DistVec<int> * = 0);

  void computeNodalHeatPower(DistSVec<double,3> &, DistSVec<double,dim> &, 
			     DistVec<double> &);
  void computeNodalHeatFluxRelatedValues(DistSVec<double,3> &, DistSVec<double,dim> &,
                                               DistVec<double> &, bool includeKappa);
  void computeForceAndMoment(Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                             DistVec<int> *,
			     Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0, 
                             VecSet< DistSVec<double,3> > *mX = 0, Vec<double> *genCF = 0);
  void computeForceAndMoment(DistExactRiemannSolver<dim>&, 
                             Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                             DistVec<int> *,
                             Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0,
                             VecSet< DistSVec<double,3> > *mX = 0, Vec<double> *genCF = 0);

  void computeHeatFluxes(DistSVec<double,3> &,
                                          DistSVec<double,dim> &, double*);

  double computeInterfaceWork(DistSVec<double,3>&, DistSVec<double,dim>&, DistVec<double>&);

  void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
			     DistSVec<double,dim> &, DistVec<double> &, 
                             DistVec<double> &, DistTimeState<dim> *);
  void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistVec<double> &,
                             DistVec<double> &, DistTimeState<dim> *,
                             DistVec<int> &);
   void computeCP(DistSVec<double,3>& X, DistSVec<double,dim>& U, Vec3D &cp);
  //void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
  //                           DistSVec<double,dim> &, DistVec<double> &,
  //                           DistSVec<double,1> &, DistVec<int> &);

  void computeVectorQuantity(PostFcn::VectorType, DistSVec<double,3> &,
			     DistSVec<double,dim> &, DistSVec<double,3> &);
  void computeVectorQuantity(PostFcn::VectorType, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistSVec<double,3> &, DistVec<int> &);
  void computeForceDerivs(DistSVec<double,3> &, DistSVec<double,dim> &,
                          DistSVec<double,dim> &,Vec<double> &,VecSet< DistSVec<double, 3> > &);

  void computeForceCoefficients(Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                                Vec3D &, Vec3D &, Vec3D &, Vec3D &,  
                                VecSet< DistSVec<double,3> > *mX = 0, DistVec<double> *genCF = 0);
  int getNumSurf() { return numSurf; }
  map<int, int> &getSurfMap() { return surfOutMap; }

  int getNumSurfHF() { return numSurfHF; }
  map<int, int> &getSurfMapHF() { return surfOutMapHF; }


// Included (MB)
  void computeDerivativeOfScalarQuantity(PostFcn::ScalarDerivativeType, double [3], DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *);

  void computeDerivativeOfVectorQuantity(PostFcn::VectorDerivativeType, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &);

  void computeDerivativeOfForceAndMoment(Vec3D &, DistSVec<double,3> &, DistSVec<double,3> &,
                                                                           DistSVec<double,dim> &, DistSVec<double,dim> &, double [3],
                                                                           Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0);

  void computeDerivativeOfNodalForce(DistSVec<double,3> &, DistSVec<double,3> &,
                                                                DistSVec<double,dim> &, DistSVec<double,dim> &,
                                                                DistVec<double> &, double [3], DistSVec<double,3> &);
								
  void rstVar(IoData &iod) {pressInfty = iod.aero.pressure;}								

  void rstVarPostFcn(IoData &ioData) {postFcn->rstVar(ioData, com);}								

  void computeDerivativeOfNodalHeatPower(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double [3], DistVec<double> &);

  void checkVec(DistSVec<double,3> &);

  void setForceGenerator(ForceGenerator<dim> *fg) { forceGen = fg; }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <PostOperator.C>
#endif

#endif

