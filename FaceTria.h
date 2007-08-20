#ifndef _FACETRIA_H_
#define _FACETRIA_H_

#include <Face.h>

//------------------------------------------------------------------------------

class FaceTria : public FaceDummy {

  static const double third;
  static const int edgeEndT[Face::MaxNumNd][2];

  int nodeNumT[3];
  int edgeNumT[3];

  // This is for triangles only (6 Vec3D's = 3 positions at t_n and t_n+1... 
  // make general function (also need to modify function calls)
  template<int dim>
  void computeForce(ElemSet &,
		    PostFcn *, SVec<double,3> &, Vec<double> &, double *, 
		    SVec<double,dim> &, double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
		    double *nodalForceWeights, int = 0);

protected:
  void *getWrapper_dim(GenFaceHelper_dim *h, 
		       int size, char *memorySpace) {
    return h->forClassTria(this,size,memorySpace);
  }
  void *getWrapper_Scalar_dim_neq(GenFaceHelper_Scalar_dim_neq *h, 
				  int size, char *memorySpace) {
    return h->forClassTria(this,size,memorySpace);
  }

  int* nodeNum() { return nodeNumT; }
  int& nodeNum(int i) { return nodeNumT[i]; }
  int& edgeNum(int i) { return edgeNumT[i]; }
  int  edgeEnd(int i, int k) { return edgeEndT[i][k]; }

public:

  // Number of nodes
  int numNodes() { return 3; }

  // Get element type
  const Type type() { return Face::TRIA; }

  // Number of normals to be stored for a triangle 
  // (the normals for all nodes equal, thus only 1)
  int numNorms() { return 1; }

  // Compute total normal and return in Vec3D
  void computeNormal(SVec<double,3> &, Vec3D &);

  // Compute subface normals and save then in Vec<Vec3D>
  void computeNormal(SVec<double,3> &, Vec<Vec3D> &);
  void computeNormalGCL1(SVec<double,3> &, SVec<double,3> &, 
			 SVec<double,3> &, Vec<Vec3D> &, Vec<double> &);
  void computeNormalEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
			   Vec<Vec3D> &, Vec<double> &);

  // Get face total normal from Vec<Vec3D>
  Vec3D  getNormal(Vec<Vec3D> &);
  double getNormalVel(Vec<double> &);

  // Get subface i normal from Vec<Vec3D>
  Vec3D  getNormal(Vec<Vec3D> &, int);
  double getNormalVel(Vec<double> &, int);

  template<int dim>
  void computeNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, 
			 Vec<double> &, double *, SVec<double,dim> &, double, 
			 SVec<double,3> &, double *nodalForceWeights);

  template<int dim>
  void computeNodalHeatPower(ElemSet &,PostFcn*, SVec<double,3>&, Vec<double>&, 
			     double*, SVec<double,dim>&, Vec<double>&);

  template<int dim>
  void computeForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
			     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, 
			     Vec3D &, Vec3D &, double *nodalForceWeights, int = 0);

  template<int dim>
  double computeInterfaceWork(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
			      double, double*, SVec<double,dim>&, double);

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType, ElemSet &, PostFcn *, SVec<double,3> &, 
			     Vec<double> &, double *, SVec<double,dim> &, SVec<double,2> &);

  template<int dim>
  void computeGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
			   Vec<double> &, double *, SVec<double,dim> &, SVec<double,dim> &);
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, double *,
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeForceDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
                          SVec<double,dim> &,SVec<double,dim> &, Vec<double> &, SVec<double,3> **);

  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, ElemSet &, SVec<double,3> &, 
				SVec<double,dim> &,  Vec<double> &, SVec<double, dim> &, 
				double, Vec3D &, Vec3D &, Vec3D &,
				Vec3D &, double *nodalForceWeights);

  template<int dim>
  void computeFDerivs(ElemSet &, VarFcn *, SVec<double,3> &, SVec<double,dim> &, Vec3D (*));

// Included (MB)
  // Get face total normal derivative from Vec<Vec3D>
  Vec3D  getdNormal(Vec<Vec3D> &);
  double getdNormalVel(Vec<double> &);

  // Get subface i normal derivative from Vec<Vec3D>
  Vec3D  getdNormal(Vec<Vec3D> &, int);
  double getdNormalVel(Vec<double> &, int);

  template<int dim>
  void computeDerivativeOfForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                       Vec<double> &, double *, double *, SVec<double,dim> &,
                                                       SVec<double,dim> &, double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double *, int = 0);
  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                                Vec<double> &, double *, double *,
                                                                SVec<double,dim> &, SVec<double,dim> &,
                                                                double, double [3], SVec<double,3> &, double *);
  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet&, PostFcn*, SVec<double,3>&, SVec<double,3>&, Vec<double>&, 
			     double*, double*, SVec<double,dim>&, SVec<double,dim>&, double [3], Vec<double>&);

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                                           Vec<double> &, double *, double *,
                                                                           SVec<double,dim> &, SVec<double,dim> &, double [3],
                                                                           Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double *, int = 0);
  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, SVec<double,3> &,
			   Vec<double> &, double *, double *, SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);
  
  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, double *, double *,
				   SVec<double,dim> &);
				   
  void computeNormalAndDerivative(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &);

  void computeDerivativeOfNormal(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &, double &, double &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <FaceTria.C>
#endif

#endif
