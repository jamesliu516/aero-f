#ifndef _FACE_H_
#define _FACE_H_

#include <PostFcn.h>
#include <MapFace.h>
#include <BlockAlloc.h>

/*
  The introduction of multiple element types in a code is traditionally 
  handled by creating a base class from which all other element types are 
  derived.
  
  In the base class, all of the required functions are declared as virtual 
  and then all the sub-classes can declare and implement their own variations
  of such functions.
  
  The problem we faced when wanting to generalize the code to any given set 
  of element (face) types is that the code heavily relies on templated 
  functions. C++ does not allow to have virtual templated functions. 
  Consequently we had to come up with an alternative approach.

  The approach we have used relies on the idea that a templated class can have 
  virtual functions. So if we have a templated function, that templated 
  function can construct a templated wrapper object that will handle the call 
  to the templated function on the actual element. The next question is how to 
  build that function when the type of the class is not specifically known at 
  compile time.  The way we do that is to call a virtual function to which we 
  pass a pointer to a Helper object. In fact that object is templated but each 
  of the templated classes is a subclass of a non-templated class that has a 
  virtual function to return the Wrapper we will need.

  To add an element of type E

  1) declare E as a subclass of A and give it all the required templated 
  functions.
  2) in E, add a function "void *getWrapper(int, char *memory)" That 
  function always has the exact same form. See in the example code.

  3) In GenHelper add the purely virtual function asClassE.
  4) in Helper: add the actual asClassE function. It is always the same 
  form and short

  To add a function:

  1) Add the function to the base class and in it create the helper, call 
  getWrapper on itself and call the function on the wrapper
  2) add the function to all subclasses.
*/

class VarFcn;
class FluxFcn;
class FemEquationTerm;
class EdgeSet;
class ElemSet;
class GeoState;
class BinFileHandler;
class TimeLowMachPrec;

struct Vec3D;

template<int dim> class BcData;
template<class Scalar> class Vec;
template<class Scalar, int dim> class GenMat;
template<class Scalar, int dim> class SVec;
template<class VecType> class VecSet;
template<class VecType, class VT2> class SubVecSet;

class FaceTria;

#define NOT_CORRECTED(msg) {						        \
    static int first_time = 1;                                                  \
    if (type()!=Face::TRIA && first_time) {				        \
      fprintf(stderr, "---------------------------------------------------\n"); \
      fprintf(stderr, "  WARNING: Using non triangular faces\n");               \
      fprintf(stderr, "  %s:%d: Function %s not fully corrected yet\n\n",	\
	      __FILE__, __LINE__, __FUNCTION__);				\
      fprintf(stderr, "  %s\n", msg);                                           \
      fprintf(stderr, "---------------------------------------------------\n"); \
      first_time = 0;                                                           \
    }                                                                           \
  }


//-------------- GENERAL HELPERS -----------------------------------------------
class GenFaceHelper_dim {
public:
  virtual  void *forClassTria(FaceTria *, int size, char *memorySpace) = 0;
};

class GenFaceHelper_Scalar_dim_neq {
public:
  virtual  void *forClassTria(FaceTria *, int size, char *memorySpace) = 0;
};


//-------------- GENERAL WRAPPERS ----------------------------------------------
template<int dim>
class GenFaceWrapper_dim {
public:
  virtual void computeNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
				 double *, SVec<double,dim> &, double, SVec<double,3> &, double* gradP[3]) = 0;
  virtual void computeNodalHeatPower(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
				     double*, SVec<double,dim>&, Vec<double>&) = 0;
  virtual void computeForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
				     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, 
				     Vec3D &, Vec3D &,  double* gradP[3], int, 
                                     SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF) = 0;
  virtual double computeInterfaceWork(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
				      double, double*, SVec<double,dim>&, double) = 0;
  virtual void computeScalarQuantity(PostFcn::ScalarType, ElemSet &, PostFcn *, SVec<double,3> &, 
				     Vec<double> &, double *, SVec<double,dim> &, SVec<double,2> &) = 0;
  virtual void computeGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, double *, SVec<double,dim> &, SVec<double,dim> &) = 0;
  virtual void computeForceDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
				  SVec<double,dim> &,SVec<double,dim> &, 
				  Vec<double> &, SVec<double,3> **) = 0;
  virtual void computeForceCoefficients(PostFcn *, Vec3D &, ElemSet &, SVec<double,3> &, 
					SVec<double,dim> &, Vec<double> &, 
					SVec<double, dim> &,  double, Vec3D &, Vec3D &, 
					Vec3D &, Vec3D &, double* gradP[3], VecSet< SVec<double,3> > *mX = 0,
                                        Vec<double> *genCF = 0) = 0;
  virtual void computeFDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
			      SVec<double,dim> &, Vec3D (*)) = 0;


// Included (MB)
  virtual void computeDerivativeOfNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                             Vec<double> &, double *, double *,
                                             SVec<double,dim> &, SVec<double,dim> &,
                                             double, double [3], SVec<double,3> &, 
                                             double* gradP[3], double* dGradP[3]) = 0;

  virtual void computeDerivativeOfNodalHeatPower(ElemSet&, PostFcn*, SVec<double,3>&, 
                                                 SVec<double,3>&, Vec<double>&, 
			                         double*, double*, SVec<double,dim>&, 
                                                 SVec<double,dim>&, double [3], Vec<double>&) = 0;

  virtual void computeDerivativeOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                 Vec<double> &, double *, double *,
                                                 SVec<double,dim> &, SVec<double,dim> &, double [3], 
                                                 Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &, 
                                                 double* gradP[3], double* dGradP[3], int = 0) = 0;

  virtual void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, SVec<double,3> &,
                                               Vec<double> &, double *, double *, SVec<double,dim> &, 
                                               SVec<double,dim> &, double, SVec<double,dim> &) = 0;
  
  virtual void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
                                            Vec<double> &, double *, double *, SVec<double,dim> &) = 0;

};

template<class Scalar, int dim, int neq>
class GenFaceWrapper_Scalar_dim_neq {
public:
  virtual void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
					   Vec<double> &, Vec<double> &, double *,
					   SVec<double,dim> &, GenMat<Scalar,neq> &) = 0;
};


//-------------- REAL WRAPPERS -------------------------------------------------
template<class Target, int dim>
class  FaceWrapper_dim : public GenFaceWrapper_dim<dim> {
  Target *t;

public:
  FaceWrapper_dim(Target *tt) : t(tt) { };

  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    t->computeNodalForce(elems, postFcn, X, d2wall, Vwall, V,
			 pin, F, gradP);
  }

  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    t->computeNodalHeatPower(elems, postFcn, X, 
			     d2wall, Vwall, V, P);
  }

  void computeForceAndMoment(ElemSet &elems,
			     PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			      double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    t->computeForceAndMoment(elems, postFcn, X,  d2wall, Vwall, V, 
			     x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }
  
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    return t->computeInterfaceWork(elems, postFcn, X, d2wall, ndot, Vwall, V, pin);
  }

  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    t->computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall, V, Q);
  }

  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R) {
    t->computeGalerkinTerm(elems, fet, X, d2wall, Vwall, V, R);
  }
  
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    t->computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, localMX);
  }

  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3], VecSet< SVec<double,3> > *mX = 0,
                                        Vec<double> *genCF = 0) {
    t->computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty, 
				CFi, CMi, CFv, CMv, gradP, mX, genCF);
  }

  void computeFDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    t->computeFDerivs(elems, varFcn, X, Vgl, F);
  }

// Included (MB)
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3], double* dGradP[3]) {
    t->computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, pin, dS, dF, gradP, dGradP);
  }

  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, 
                                         SVec<double,3>& dX, Vec<double>& d2wall, double* Vwall, 
                                         double* dVwall, SVec<double,dim>& V, SVec<double,dim>& dV, 
                                         double dS[3], Vec<double>& dP) {
    t->computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, dP);
  }

  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                         SVec<double,3> &X, SVec<double,3> &dX,
                                         Vec<double> &d2wall, double *Vwall, double *dVwall,
                                         SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                         Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
                                         double* gradP[3], double* dGradP[3], int hydro) {
    t->computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
  }

  void computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, 
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    t->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }

  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                    double *Vwall, double *dVwall, SVec<double,dim> &V) {
    t->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall, dVwall, V);
  }

};

template<class Target, class Scalar, int dim, int neq>
class  FaceWrapper_Scalar_dim_neq : public 
GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> {

  Target *t;

public:
  FaceWrapper_Scalar_dim_neq(Target *tt) : t(tt) { };

  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    t->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall, V, A);
  };
};


//-------------- REAL HELPERS --------------------------------------------------
template<int dim>
class FaceHelper_dim : public GenFaceHelper_dim {
public:

  void *forClassTria(FaceTria *tface, int size, char *memorySpace) {
    if(size < sizeof(FaceWrapper_dim<FaceTria, dim>) ) {
      fprintf(stderr, "Error: programming error in FaceHelper");
      exit(1);
    }
    return new (memorySpace) FaceWrapper_dim<FaceTria, dim>(tface);
  }

};


template<class Scalar, int dim, int neq>
class FaceHelper_Scalar_dim_neq : public GenFaceHelper_Scalar_dim_neq {
public:

  void *forClassTria(FaceTria *tface, int size, char *memorySpace) {
    if(size < sizeof(FaceWrapper_Scalar_dim_neq<FaceTria, Scalar, dim, neq>) ) {
      fprintf(stderr, "Error: programming error in FaceHelper");
      exit(1);
    }
    return new (memorySpace) FaceWrapper_Scalar_dim_neq<FaceTria, Scalar, dim, neq>(tface);
  }
  
};


//--------------- BASE FACE CLASS ----------------------------------------------
class Face {
public:
  enum Type {TRIA=4};

  static const int MaxNumNd = 4;

  virtual int nodeNum(int i) const = 0;

protected:
  virtual void *getWrapper_dim(GenFaceHelper_dim *, 
			       int size, char *memorySpace) = 0;
  virtual void *getWrapper_Scalar_dim_neq(GenFaceHelper_Scalar_dim_neq *, 
					  int size, char *memorySpace) = 0;
  
  int code;
  int elemNum;
  int surface_id;
  int normNum;

  virtual int* nodeNum() = 0;
  virtual int& nodeNum(int i) = 0;
  virtual int& edgeNum(int i) = 0;
  virtual int  edgeEnd(int i, int k) = 0;  

public:

  // Number of nodes
  virtual int numNodes() = 0;

  // Number of normals to be stored
  virtual int numNorms() = 0;

  // Get element type
  virtual const Type type() = 0;

  int operator[](int i) { return nodeNum(i); }  
  operator int *() { return nodeNum(); }
  operator MaxFace() { return MaxFace(numNodes(), nodeNum()); }

  /* WARNING : IS THIS THE RIGHT DEFINITION, WHEN NUMNODES()!=F.NUMNODES() ??? */
  /*           removed const ... function is actualy never used */
  bool operator<(Face &f)  {
    if(numNodes() != f.numNodes()) return numNodes() < f.numNodes();
    for (int i=0; i<numNodes(); i++) {
      if (nodeNum(i) < f.nodeNum(i)) return true;
      if (nodeNum(i) > f.nodeNum(i)) return false;
    }
    
    return false;
  }

  int* nodes(int* d = NULL) {  
    if (d)
      for (int i; i<numNodes(); i++) d[i] = nodeNum(i); 
    return nodeNum();
  }
  
  int getCode()  { return code; }
  int getSurfaceID() { return surface_id; }
  int getElementNumber() const { return elemNum; }

  void setup(int, int *, int, int surface_id = 0);
  void setType(int *);
  void setType(int t) { code = t; }
  void setNodeType(int*, int*);
  void setNodeFaceType(int*);
  void setElementNumber(int elemNum, int rotDir);
  void tagNodesOnBoundaries(Vec<bool> &);
  void tagEdgesOnBoundaries(Vec<bool> &);
  void reorder();
  void numberEdges(EdgeSet &);
  void computeEdgeNormals(SVec<double,3>&, int*, SVec<double,6>&);

  // WARNING: THIS IS A FUNCTION FOR TFACES ONLY
  static double computeVolume(Vec3D &xa_n, Vec3D &xb_n, Vec3D &xc_n, 
			      Vec3D &xa_np1, Vec3D &xb_np1, Vec3D &xc_np1);
  
  // Compute total normal and return in Vec3D  
  virtual void computeNormal(SVec<double,3> &, Vec3D &) = 0;

  // Compute subface normals and save then in Vec<Vec3D>
  virtual void computeNormal(SVec<double,3> &, Vec<Vec3D> &) = 0;
  virtual void computeNormalConfig(SVec<double,3> &, SVec<double,3> &,
                                   Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeNormalGCL1(SVec<double,3> &, SVec<double,3> &, 
				 SVec<double,3> &, Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeNormalEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
				   Vec<Vec3D> &, Vec<double> &) = 0;

  // Get face total normal from Vec<Vec3D>
  virtual Vec3D getNormal(Vec<Vec3D> &) = 0;

  // Get subface i normal from Vec<Vec3D>
  virtual Vec3D getNormal(Vec<Vec3D> &, int) = 0;

  // Get face total normal velocity from Vec<double>
  virtual double getNormalVel(Vec<double> &) = 0;

  // Get subface i normal velocity from Vec<double>
  virtual double getNormalVel(Vec<double> &, int) = 0;

  template<class NodeMap>
  void renumberNodes(NodeMap &nodemap);

  template<int dim>
  void assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U);

  template<int dim>
  void assignFreeStreamValues(double *Uin, double *Uout, double *U);
  
  template<int dim>
  void computeFaceBcValue(SVec<double,dim> &Unode, double *Uface);

  template<int dim1, int dim2>
  void computeNodeBcValue(SVec<double,3> &X, double *Uface, SVec<double,dim2> &Unode);
    
  // "Virtual template" functions (implemented through wrapper and helper functions)

  template<int dim>
  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeNodalForce(elems, postFcn, X, d2wall, Vwall, V,
			       pin, F, gradP);
  }

  template<int dim>
  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeNodalHeatPower(elems, postFcn, X, 
			     d2wall, Vwall, V, P);
  }

  template<int dim>
  void computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			     double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceAndMoment(elems, postFcn, X,  d2wall, Vwall, V, 
			     x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }
  
  template<int dim>
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    return wrapper->computeInterfaceWork(elems, postFcn, X, d2wall, ndot, Vwall, V, pin);
  }

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall, V, Q);
  }

  template<int dim>
  void computeTimeStep(VarFcn *varFcn, Vec<Vec3D> &normal, Vec<double> &normalVel,
		       SVec<double,dim> &V, Vec<double> &dt, 
		       TimeLowMachPrec &tprec, Vec<double> &Phi);

  template<int dim>
  void computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, 
		       Vec<Vec3D> &normal, Vec<double> &normalVel,
		       SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, 
		       Vec<double> &idtv, TimeLowMachPrec &tprec);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
			       Vec<double> &normalVel, SVec<double,dim> &V, 
			       double *Ub, SVec<double,dim> &fluxes);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
			       Vec<double> &normalVel, SVec<double,dim> &V,
			       double *Ub, Vec<double> &Phi, 
			       SVec<double,dim> &fluxes);

  template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				 Vec<double> &normalVel, SVec<double,dim> &V,
				 SVec<double,1> &Phi, Vec<double> &PhiF);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
				       Vec<double> &normalVel, SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				       Vec<double> &normalVel, SVec<double,dim> &V,
				       double *Ub, GenMat<Scalar,neq> &A, int* nodeType);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
				       Vec<double> &normalVel, SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A, Vec<double> &Phi);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				       Vec<double> &normalVel, SVec<double,dim> &V,
				       double *Ub, GenMat<Scalar,neq> &A, Vec<double> &Phi,
                                       int* nodeType);

  template<int dim, class Scalar>
  void computeJacobianFiniteVolumeTermLS(Vec<Vec3D> &normal,
					 Vec<double> &normalVel, SVec<double,dim> &V,
					 GenMat<Scalar,1> &A);

  template<int dim>
  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeGalerkinTerm(elems, fet, X, d2wall, Vwall, V, R);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    FaceHelper_Scalar_dim_neq<Scalar,dim,neq> h;
    char xx[64];
    GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> *wrapper=
      (GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> *)getWrapper_Scalar_dim_neq(&h, 64, xx);
    wrapper->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall, V, A);
  }

  template<int dim>
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, localMX);
  }

  template<int dim>
  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty, 
				CFi, CMi, CFv, CMv, gradP);
  }

  template<int dim>
  void computeFDerivs(ElemSet &elems,
		      VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeFDerivs(elems, varFcn, X, Vgl, F);
  }
 
// Included (MB)
  virtual void computeNormalAndDerivative(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D&) = 0;

  virtual void computeDerivativeOfNormal(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &, double &, double &) = 0;

  // Get face total normal derivative from Vec<Vec3D>
  virtual Vec3D getdNormal(Vec<Vec3D> &) = 0;

  // Get subface i normal derivative from Vec<Vec3D>
  virtual Vec3D getdNormal(Vec<Vec3D> &, int) = 0;

  // Get face total normal velocity derivative from Vec<double>
  virtual double getdNormalVel(Vec<double> &) = 0;

  // Get subface i normal velocity derivative from Vec<double>
  virtual double getdNormalVel(Vec<double> &, int) = 0;

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<Vec3D> &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
				      SVec<double,dim> &V, double *Ub,
				      double *dUb, SVec<double,dim> &dFluxes);

  template<int dim1, int dim2>
  void computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, double *Uface, double *dUface, SVec<double,dim2> &dUnode);

  template<int dim>
  void computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, double *dUfaceSA, SVec<double,dim> &dUnodeSA);

  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D>  &normals, Vec<Vec3D>  &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
			   SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			   Vec<double> &dIdti, Vec<double> &dIdtv, double dMach, 
                           TimeLowMachPrec &tprec);

  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3],  double* dGradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, pin, dS, dF, gradP, dGradP);
  }

  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, SVec<double,3>& dX, 
				 Vec<double>& d2wall, double* Vwall, double* dVwall, SVec<double,dim>& V, SVec<double,dim>& dV, double dS[3], Vec<double>& dP) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, dP);
  }

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                             SVec<double,3> &X, SVec<double,3> &dX,
                                             Vec<double> &d2wall, double *Vwall, double *dVwall,
                                             SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                             Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
				             double* gradP[3], double* dGradP[3], int hydro) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
  }

  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }
  
  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                        double *Vwall, double *dVwall, SVec<double,dim> &V) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall, dVwall, V);
  }

};


//--------------- FACE CLASS ---------------------------------------------------
class FaceDummy :  public Face {

public:

  // Ensure that when "Virtual template" are not defined in derived classes
  // an error is thrown to avoid infinite loop in wrapper function
  
  template<int dim>
  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			     double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }
  
  /* WARNING : THIS FUNCTION IS RETURNING A DOUBLE ? ... IS THIS A PROBLEM ? */
  template<int dim>
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3], VecSet< SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeFDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }

// Included (MB)
  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3], double* dGradP[3]) {
    fprintf(stderr, "Error: undifined function (computeDerivativeOfNodalForce) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, SVec<double,3>& dX, 
                                         Vec<double>& d2wall, double* Vwall, double* dVwall, SVec<double,dim>& V, 
                                         SVec<double,dim>& dV, double dS[3], Vec<double>& dP) {

    fprintf(stderr, "Error: undifined function (computeDerivativeOfNodalHeatPower) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                             SVec<double,3> &X, SVec<double,3> &dX,
                                             Vec<double> &d2wall, double *Vwall, double *dVwall,
                                             SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                             Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
				             double* gradP[3], double* dGradP[3], int hydro) {

    fprintf(stderr, "Error: undifined function (computeDerivativeOfForceAndMoment) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, 
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    fprintf(stderr, "Error: undifined function (computeDerivativeOfGalerkinTerm) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                    double *Vwall, double *dVwall, SVec<double,dim> &V) {
    fprintf(stderr, "Error: undifined function (computeBCsJacobianWallValues) for this face type\n"); exit(1);
  }

};


//------------------------------------------------------------------------------
class FaceSet {

  int numFaces;
  int numFaceNorms;

  Face **faces;
  BlockAlloc memFaces;
  
public:

  /* Need to define constructor and destructor! */
  FaceSet(int);
  ~FaceSet();

  Face &operator[](int i) { return *faces[i]; }

  void addFace(int i, Face *face) { faces[i] = face; }

  int size() const { return numFaces; }

  // Number of face normals to be stored
  int sizeNorms() const { return numFaceNorms; }

  int read(BinFileHandler &, int, int (*)[2], int *);

  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<double> &);
  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
		       SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, TimeLowMachPrec &);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, Vec<double> &, 
                               SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn **, BcData<dim> &, GeoState &, 
				 SVec<double,dim> &, SVec<double,1> &, Vec<double> &);

  // DEBUG /* Not implemented */
  template<int dim>
  void computeInviscidFluxes(FluxFcn **, BcData<dim> &, GeoState &,
                             SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &, int*);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       Vec<double> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       Vec<double> &, int*);

  template<int dim, class Scalar>
  void computeJacobianFiniteVolumeTermLS(GeoState &, 
					 SVec<double,dim> &, GenMat<Scalar,1> &);

  template<int dim>
  void computeGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, GeoState &, 
			   SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, 
				   GeoState &, SVec<double,3> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

// Included (MB)
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                           SVec<double,dim> &, SVec<double,dim> &);
  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, GeoState &,
                                       SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, 
                                       SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, BcData<dim> &, 
				   GeoState &, SVec<double,3> &, SVec<double,dim> &);
  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
			      SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, 
			      Vec<double> &, Vec<double> &, double, 
                              TimeLowMachPrec &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Face.C>
#endif

#include <FaceTria.h>

#endif
