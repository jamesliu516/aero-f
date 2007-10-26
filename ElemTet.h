#ifndef _ELEMTET_H_
#define _ELEMTET_H_

#include <Elem.h>

//------------------------------------------------------------------------------

class ElemTet : public ElemDummy {

  static const double third;
  static const double fourth;
  static const double sixth;
  static const int edgeEndTet[6][2];
  static const int edgeFaceTet[6][2];
  static const int faceDefTet[4][3];

  int nodeNumTet[4];
  int edgeNumTet[6];

protected:
  void *getWrapper_dim(GenElemHelper_dim *h, 
		       int size, char *memorySpace) {
    return h->forClassTet(this,size,memorySpace);
  }
  void *getWrapper_Scalar_dim_neq(GenElemHelper_Scalar_dim_neq *h, 
				  int size, char *memorySpace) {
    return h->forClassTet(this,size,memorySpace);
  }

public:

  ElemTet() { volume_id = 0; }
  ~ElemTet() {}

  int* nodeNum() { return nodeNumTet; }
  int& nodeNum(int i) { return nodeNumTet[i]; }
  int& edgeNum(int i) { return edgeNumTet[i]; }
  const int  edgeEnd(int i, int k) { return edgeEndTet[i][k]; }
  const int  edgeFace(int i, int k) { return edgeFaceTet[i][k]; }
  const int  faceDef(int i, int k) { return faceDefTet[i][k]; }
  const int  faceNnd(int i) { return 3; }
  const Type type() { return Elem::TET; }

  // Number of nodes
  int numNodes() { return 4; }

  // Number of edges
  int numEdges() { return 6; }

  // Number of faces
  int numFaces() { return 4; }

  //--------------functions in TetCore.C

  double computeLongestEdge(SVec<double,3> &);
  double computeVolume(SVec<double,3> &);
  double computeControlVolumes(SVec<double,3> &, Vec<double> &);
  void printInvalidElement(int, double, int, int *, int *, SVec<double,3> &, SVec<double,3> &);
  void computeEdgeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      Vec<Vec3D> &, Vec<double> &);
  void computeEdgeNormalsEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
				Vec<Vec3D> &, Vec<double> &);
  void computeWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, 
			      SVec<double,3> &, SVec<double,3> &);
  void computeEdgeWeightsGalerkin(SVec<double,3> &, SVec<double,9> &);
  double computeGradientP1Function(SVec<double,3> &, double [4][3], double * = NULL);
  double computeGradientP1Function(Vec3D &A, Vec3D &B, Vec3D &C, Vec3D &D, 
                                   double nGrad[4][3]);
  void computeStiffAndForce(double *, double *, 
			    SVec<double, 3> &, SVec<double,3> &, double volStiff = 0.0);
  void computeStiffAndForceLIN(double *, SVec<double,3> &, SVec<double,3> &);
  void computeStiffBallVertex(double *, SVec<double, 3> &);
  void computeStiffTorsionSpring(double *, SVec<double, 3> &);

  void computeLij(double [3][3], double [3], double [6], double [5]);

  void computeBij(double [3][3], double [6], double , double [3][3], double , double [5]);

  void computeLi(double [3], double , double , double [3], double [5]);

  void computeZi(double [3], double , double , double [3], double [3], double [5], double , double); 

  void computePij(double [3][3], double [3][3]);

  void computeTemp(double *[4], double [4], double);

  void computeTempGradient(double [4][3], double [4], double [3]);

  double computeNormSij(double [3][3]);

  void computeVelocity(double *[4], double [4][3],double [3], double [4]);

  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);

  
  //-----functions in Tet.C

  template<class NodeMap>
  void renumberNodes(NodeMap &);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &, SVec<double,dim> &,
                       SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                             Vec<double> *, Vec<double> &);
                                                                                                                          
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
				 SVec<double,dim> &R);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
		          SVec<double,dim> &R);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, Vec<double> &, 
			     SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, Vec<double> &, Vec<double> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
			       SVec<double,3> &, Vec<double> &, double *, 
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
				       SVec<double,3> &, Vec<double> &, Vec<double> &, 
				       double *, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeP1Avg(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, 
		    SVec<double,3> &, SVec<double,dim> &, double, double);

  template<int dim>
  void computeCsValues(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, SVec<double,2> &, Vec<double> &,
		       SVec<double,3> &, double, double);


// Included (MB)
  double computeDerivativeOfVolume(SVec<double,3> &, SVec<double,3> &);

  double computeDerivativeOfControlVolumes(SVec<double,3> &, SVec<double,3> &, Vec<double> &);

  void computeDerivativeOfEdgeNormals(SVec<double,3> &, SVec<double,3> &, Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &);

  void computeDerivativeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);

  double computeDerivativeOfGradientP1Function(SVec<double,3> &, SVec<double,3> &, double [4][3]);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, SVec<double,3> &, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, Vec3D &,
			       SVec<double,3> &, SVec<double,3> &, Vec<double> &, double *, double *,
			       SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

// Level Set Reinitialization
  template<int dim>
  void computePsiResidual(SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                          SVec<double,dim> &ddx,SVec<double,dim> &ddy,
                          SVec<double,dim> &ddz, Vec<int> &Tag,
                          Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes,
			  int typeTracking);

  template<int dim>
  void computeDistanceCloseNodes(Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 Vec<double> &Phi,SVec<double,dim> &Psi);

  template<int dim>
  void recomputeDistanceCloseNodes(Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 Vec<double> &Phi,SVec<double,dim> &Psi);

  template<int dim>
  void computeDistanceLevelNodes(Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,dim> &Psi, Vec<double> &Phi);



private:

  //--------------functions in ElemTetCore.C

//Level Set Reinitialization
  void computePsiResidualSubTet(double psi[4], double phi[4],
                                Vec3D A, Vec3D B, Vec3D C, Vec3D D,
                                double locdphi[4], double locw[4],
                                double locbeta[4], bool debug);
  double findRootPolynomialNewtonRaphson(double f1, double f2, double fp1, double fp2);
  int findRootPolynomialLaguerre(double f1, double f2, double fp1, double fp2, double &root);
  bool computeDistancePlusPhiToOppFace(double phi[3], Vec3D Y0,
                                       Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToEdges(double phi[3], Vec3D Y0,
                                     Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToVertices(double phi[3], Vec3D Y0,
                                        Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToEdge(double phi0, double phi1,
                                    Vec3D Y0, Vec3D Y1, double &mini, bool show = false);
  int computeDistanceToAll(double phi[3],Vec3D Y0,Vec3D Y1,Vec3D Y2, double &psi);

  //--------------functions in ElemTet.C

//Level Set Reinitialization
  template<int dim>
  int findLSIntersectionPoint(Vec<double> &Phi, SVec<double,dim> &ddx,
                              SVec<double,dim> &ddy, SVec<double,dim> &ddz,
 			      SVec<double,3> &X,
                              int reorder[4], Vec3D P[4], int typeTracking);
  template<int dim>
  void findLSIntersectionPointLinear(Vec<double> &Phi, SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
				 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);
  template<int dim>
  void findLSIntersectionPointGradient(Vec<double> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
				 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);
  template<int dim>
  int findLSIntersectionPointHermite(Vec<double> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);

  template<int dim>
  void computePsiResidual0(SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                           Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes, bool debug);

  template<int dim>
  void computePsiResidual1(int reorder[4], Vec3D P[4],
                           SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                           Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes, bool debug);

  template<int dim>
  void computePsiResidual2(int reorder[4], Vec3D P[4],
                           SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                           Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes, bool debug);

  template<int dim>
  void computeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                  Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);
  template<int dim>
  void recomputeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);

  template<int dim>
  double computeDistancePlusPhi(int i, SVec<double,3> &X, SVec<double,dim> &Psi);



};

#ifdef TEMPLATE_FIX
#include <ElemTet.C>
#endif

#endif
