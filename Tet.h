#ifndef _TET_H_
#define _TET_H_

//#include <Face.h>
#include <MacroCell.h>
#include <Vector.h>
#include <Vector3D.h>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

class FemEquationTerm;
class MacroCellSet;
class VMSLESTerm;
class DynamicVMSTerm;
class SmagorinskyLESTerm;
class DynamicLESTerm;
class EdgeSet;
class BinFileHandler;
class GeoState;
class FaceSet;
class Face;


template<class Scalar, int dim> class GenMat;

#ifdef OLD_STL
typedef map<Face, int, less<Face> > MapFaces;
#else
typedef map<Face, int> MapFaces;
#endif

#include <complex.h>
typedef complex<double> bcomp;

//------------------------------------------------------------------------------

class Tet {

  static const double third;
  static const double fourth;
  static const double sixth;


  int nodeNum[4];
  int edgeNum[6];
  int volume_id;

public:

  static const int edgeEnd[6][2];
  static const int edgeFace[6][2];
  static const int faceDef[4][3];

public:

  Tet() {volume_id=-1;}
  ~Tet() {}

  int &operator[](int i) { return nodeNum[i]; }
  operator int *() { return nodeNum; }

  void setVolumeID(int i) {volume_id = i;}
  int getVolumeID() { return volume_id; }

//--------------functions in TetCore.C

  double computeLongestEdge(SVec<double,3> &);
  double computeVolume(SVec<double,3> &);
  double computeGeometricVolume(Vec3D , Vec3D , Vec3D , Vec3D );
  double computeControlVolumes(SVec<double,3> &, Vec<double> &);
  void printInvalidElement(int, double, int, int *, int *, SVec<double,3> &, SVec<double,3> &);
  void numberEdges(EdgeSet &);
  void renumberEdges(Vec<int> &);
  int countNodesOnBoundaries(Vec<bool> &);
  int setFaceToElementConnectivity(int i, Vec<bool> &, MapFaces &, FaceSet &);
  void computeEdgeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      Vec<Vec3D> &, Vec<double> &);
  void computeEdgeNormalsEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
				Vec<Vec3D> &, Vec<double> &);
  void computeWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, 
			      SVec<double,3> &, SVec<double,3> &);
  void computeEdgeWeightsGalerkin(SVec<double,3> &, SVec<double,9> &);
  double computeGradientP1Function(SVec<double,3> &, double [4][3]);
  double computeGradientP1Function(Vec3D &A, Vec3D &B, Vec3D &C, Vec3D &D, double [4][3]);
  void computeStiffAndForce(double [4][3], double [12][12], 
			    SVec<double, 3> &, SVec<double,3> &, double volStiff = 0.0);
  void computeStiffAndForceLIN(double [4][3], double [12][12],
			       SVec<double,3> &, SVec<double,3> &);
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

//-----functions in Tet.h

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
  void computePsiResidual(SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                          SVec<double,dim> &ddx,SVec<double,dim> &ddy,SVec<double,dim> &ddz,
                          Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes,
			  int typeTracking);

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
  void computePsiResidualFM(SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
                          SVec<double,dim> &ddx,SVec<double,dim> &ddy,SVec<double,dim> &ddz,
                          Vec<double> &w,Vec<double> &beta, SVec<double,dim> &PsiRes,
			  int typeTracking);
  template<int dim>
  void computeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                  Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);
  template<int dim>
  void recomputeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                  Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);
  template<int dim>
  void computeFM(SVec<double,dim> &Psi);

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
  template<int dim>
  double computeDistancePlusPhi(int i, SVec<double,3> &X, SVec<double,dim> &Psi);

};
//------------------------------------------------------------------------------

class TetSet {

  int numTets;
  Tet *tets;
  
public:

  TetSet(int);
  ~TetSet();

  Tet &operator[](int i) const { return tets[i]; }

  int read(BinFileHandler&, int, int (*)[2], int *);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &, 
			   SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                             Vec<double> *, Vec<double> &);

  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &,
                       SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &,
				 SVec<double,dim> &);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, Vec<double> &, 
			     SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, Vec<double> &, Vec<double> &);

  template<int dim>
  void computeTestFilterAvgs(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &,
		  SVec<double,3> &, SVec<double,dim> &, double, double);

  template<int dim>
  void computeCsValues(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &,
		  SVec<double,2> &, Vec<double> &, SVec<double,3> &, double, double);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &,
				   Vec<double> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computePsiResidual(SVec<double,3> &X,Vec<double> &Phi,SVec<double,dim> &Psi,
			  SVec<double,dim> &ddx, SVec<double,dim> &ddy,
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


  int size() const { return numTets; }
  
};

//------------------------------------------------------------------------------

inline
double Tet::computeGradientP1Function(SVec<double,3> &nodes, double nGrad[4][3])
{

  double jac[3][3];

  //Jacobian
  // J_ij = dx_i/dxi_j
  jac[0][0] = nodes[ nodeNum[1] ][0] - nodes[ nodeNum[0] ][0];
  jac[0][1] = nodes[ nodeNum[2] ][0] - nodes[ nodeNum[0] ][0];
  jac[0][2] = nodes[ nodeNum[3] ][0] - nodes[ nodeNum[0] ][0];
  jac[1][0] = nodes[ nodeNum[1] ][1] - nodes[ nodeNum[0] ][1];
  jac[1][1] = nodes[ nodeNum[2] ][1] - nodes[ nodeNum[0] ][1];
  jac[1][2] = nodes[ nodeNum[3] ][1] - nodes[ nodeNum[0] ][1];
  jac[2][0] = nodes[ nodeNum[1] ][2] - nodes[ nodeNum[0] ][2];
  jac[2][1] = nodes[ nodeNum[2] ][2] - nodes[ nodeNum[0] ][2];
  jac[2][2] = nodes[ nodeNum[3] ][2] - nodes[ nodeNum[0] ][2];

  // compute determinant of jac
  double dOmega = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) +
                  jac[1][0] * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) +
                  jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);

  // compute inverse matrix of jac
  // Maple code used
  double t17 = -1.0/dOmega;

  //compute shape function gradients
  nGrad[1][0] =  (-jac[1][1] * jac[2][2] + jac[1][2] * jac[2][1] ) * t17;
  nGrad[1][1] =  ( jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1] ) * t17;
  nGrad[1][2] = -( jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1] ) * t17;

  nGrad[2][0] = -(-jac[1][0] * jac[2][2] + jac[1][2] * jac[2][0] ) * t17;
  nGrad[2][1] = -( jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0] ) * t17;
  nGrad[2][2] =  ( jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0] ) * t17;

  nGrad[3][0] = -( jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0] ) * t17;
  nGrad[3][1] =  ( jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0] ) * t17;
  nGrad[3][2] = -( jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0] ) * t17;

  // Shape function gradients dN_i/dx_i = dN/dxi * transpose(jInv)
  // Note: 1st index = shape function #
  // 2nd index = direction (0=x, 1=y, 2=z)

  nGrad[0][0] = -( nGrad[1][0] + nGrad[2][0] + nGrad[3][0] );
  nGrad[0][1] = -( nGrad[1][1] + nGrad[2][1] + nGrad[3][1] );
  nGrad[0][2] = -( nGrad[1][2] + nGrad[2][2] + nGrad[3][2] );

  return sixth * dOmega;

}

//------------------------------------------------------------------------------

inline
double Tet::computeGradientP1Function(Vec3D &A, Vec3D &B, Vec3D &C, Vec3D &D, 
                                      double nGrad[4][3])
{

  //fprintf(stdout, "A = %e %e %e\n", A[0],A[1],A[2]);
  //fprintf(stdout, "B = %e %e %e\n", B[0],B[1],B[2]);
  //fprintf(stdout, "C = %e %e %e\n", C[0],C[1],C[2]);
  //fprintf(stdout, "D = %e %e %e\n", D[0],D[1],D[2]);
  
  double jac[3][3];

  //Jacobian
  // J_ij = dx_i/dxi_j
  double v = B[0];
  jac[0][0] = B[0] - A[0];
  jac[0][1] = C[0] - A[0];
  jac[0][2] = D[0] - A[0];
  jac[1][0] = B[1] - A[1];
  jac[1][1] = C[1] - A[1];
  jac[1][2] = D[1] - A[1];
  jac[2][0] = B[2] - A[2];
  jac[2][1] = C[2] - A[2];
  jac[2][2] = D[2] - A[2];

  // compute determinant of jac
  double dOmega = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) +
                  jac[1][0] * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) +
                  jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);

  // compute inverse matrix of jac
  // Maple code used
  //fprintf(stdout, "dOmega = %e\n", dOmega);
  double t17 = -1.0/dOmega;

  //compute shape function gradients
  nGrad[1][0] =  (-jac[1][1] * jac[2][2] + jac[1][2] * jac[2][1] ) * t17;
  nGrad[1][1] =  ( jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1] ) * t17;
  nGrad[1][2] = -( jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1] ) * t17;

  nGrad[2][0] = -(-jac[1][0] * jac[2][2] + jac[1][2] * jac[2][0] ) * t17;
  nGrad[2][1] = -( jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0] ) * t17;
  nGrad[2][2] =  ( jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0] ) * t17;

  nGrad[3][0] = -( jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0] ) * t17;
  nGrad[3][1] =  ( jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0] ) * t17;
  nGrad[3][2] = -( jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0] ) * t17;

  // Shape function gradients dN_i/dx_i = dN/dxi * transpose(jInv)
  // Note: 1st index = shape function #
  // 2nd index = direction (0=x, 1=y, 2=z)

  nGrad[0][0] = -( nGrad[1][0] + nGrad[2][0] + nGrad[3][0] );
  nGrad[0][1] = -( nGrad[1][1] + nGrad[2][1] + nGrad[3][1] );
  nGrad[0][2] = -( nGrad[1][2] + nGrad[2][2] + nGrad[3][2] );

  //fprintf(stdout, "dOmega = %e\n", sixth*dOmega);
  return sixth*dOmega;




}

//------------------------------------------------------------------------------
/*
inline
double Tet::computeGradientP1Function(SVec<double,3> &X, double dp1dxj[4][3])
{

  double x1 = X[ nodeNum[0] ][0];
  double y1 = X[ nodeNum[0] ][1];
  double z1 = X[ nodeNum[0] ][2];

  double x2 = X[ nodeNum[1] ][0];
  double y2 = X[ nodeNum[1] ][1];
  double z2 = X[ nodeNum[1] ][2];

  double x3 = X[ nodeNum[2] ][0];
  double y3 = X[ nodeNum[2] ][1];
  double z3 = X[ nodeNum[2] ][2];

  double x4 = X[ nodeNum[3] ][0];
  double y4 = X[ nodeNum[3] ][1];
  double z4 = X[ nodeNum[3] ][2];

  double z12 = z2 - z1;
  double z13 = z3 - z1;
  double z14 = z4 - z1;
  double z23 = z3 - z2;
  double z24 = z4 - z2;
  double z34 = z4 - z3;

  dp1dxj[0][0] = y2*z34 - y3*z24 + y4*z23;
  dp1dxj[1][0] = -y1*z34 + y3*z14 - y4*z13;
  dp1dxj[2][0] = y1*z24 - y2*z14 + y4*z12;
  dp1dxj[3][0] =-y1*z23 + y2*z13 - y3*z12;

  dp1dxj[0][1] =-x2*z34 + x3*z24 - x4*z23;
  dp1dxj[1][1] = x1*z34 - x3*z14 + x4*z13;
  dp1dxj[2][1] =-x1*z24 + x2*z14 - x4*z12;
  dp1dxj[3][1] = x1*z23 - x2*z13 + x3*z12;

  double y12 = y2 - y1;
  double y13 = y3 - y1;
  double y14 = y4 - y1;
  double y23 = y3 - y2;
  double y24 = y4 - y2;
  double y34 = y4 - y3;

  dp1dxj[0][2] = x2*y34 - x3*y24 + x4*y23;
  dp1dxj[1][2] =-x1*y34 + x3*y14 - x4*y13;
  dp1dxj[2][2] = x1*y24 - x2*y14 + x4*y12;
  dp1dxj[3][2] =-x1*y23 + x2*y13 - x3*y12;

  double vol6 = x1*dp1dxj[0][0] + x2*dp1dxj[1][0] + x3*dp1dxj[2][0] + x4*dp1dxj[3][0];

  double invvol6 = 1.0 / vol6;

  dp1dxj[0][0] *= invvol6;
  dp1dxj[0][1] *= invvol6;
  dp1dxj[0][2] *= invvol6;

  dp1dxj[1][0] *= invvol6;
  dp1dxj[1][1] *= invvol6;
  dp1dxj[1][2] *= invvol6;

  dp1dxj[2][0] *= invvol6;
  dp1dxj[2][1] *= invvol6;
  dp1dxj[2][2] *= invvol6;

  dp1dxj[3][0] *= invvol6;
  dp1dxj[3][1] *= invvol6;
  dp1dxj[3][2] *= invvol6;

  return sixth * vol6;

}
*/
//------------------------------------------------------------------------------

inline
void Tet::computeLij(double Lij[3][3], double r_u[3],
                     double r_u_u[6], double vc[5])
{
                                                                                                                                          
  if (vc[0] < 0.0000001) vc[0] = 0.0000001;
                                                                                                                                          
  Lij[0][0] = r_u_u[0] - (1.0/vc[0])*(r_u[0]*r_u[0]);
  Lij[0][1] = r_u_u[1] - (1.0/vc[0])*(r_u[0]*r_u[1]);
  Lij[0][2] = r_u_u[2] - (1.0/vc[0])*(r_u[0]*r_u[2]);
  Lij[1][1] = r_u_u[3] - (1.0/vc[0])*(r_u[1]*r_u[1]);
  Lij[1][2] = r_u_u[4] - (1.0/vc[0])*(r_u[1]*r_u[2]);
  Lij[2][2] = r_u_u[5] - (1.0/vc[0])*(r_u[2]*r_u[2]);
  Lij[2][0] = Lij[0][2];
  Lij[2][1] = Lij[1][2];
  Lij[1][0] = Lij[0][1];
  Lij[0][0] = Lij[0][0] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);
  Lij[1][1] = Lij[1][1] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);
  Lij[2][2] = Lij[2][2] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------
                                                                                                                                          
inline
void Tet::computeBij(double Bij[3][3], double r_s_p[6], double sqrt2S2,
                          double Pij[3][3], double sq_rat_delta, double vc[5])
                                                                                                                                          
{
                                                                                                                                          
   Bij[0][0] = r_s_p[0] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][0];
   Bij[0][1] = r_s_p[3] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][1];
   Bij[0][2] = r_s_p[4] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][2];
   Bij[1][1] = r_s_p[1] - sq_rat_delta*vc[0]*sqrt2S2*Pij[1][1];
   Bij[1][2] = r_s_p[5] - sq_rat_delta*vc[0]*sqrt2S2*Pij[1][2];
   Bij[2][2] = r_s_p[2] - sq_rat_delta*vc[0]*sqrt2S2*Pij[2][2];
   Bij[1][0] = Bij[0][1];
   Bij[2][0] = Bij[0][2];
   Bij[2][1] = Bij[1][2];
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------
                                                                                                                                          
inline
void Tet::computeLi(double Li[3], double r_e, double r_e_plus_p,
                   double r_u[3], double vc[5])
{
                                                                                                                                          
  if (vc[0] < 0.0000001) vc[0] = 0.0000001;
  Li[0] = (r_e+vc[4])*(r_u[0]/vc[0]) - (r_e_plus_p)*vc[1];
  Li[1] = (r_e+vc[4])*(r_u[1]/vc[0]) - (r_e_plus_p)*vc[2];
  Li[2] = (r_e+vc[4])*(r_u[2]/vc[0]) - (r_e_plus_p)*vc[3];
                                                                                                                                          
}
                                                                                                                                          
//------------------------------------------------------------------------------

inline
void Tet::computeZi(double Zi[3], double ratdelta, double sqrt2S2, double dtdxj[3],
                    double r_s_dtdxj[3], double vc[5], double gamma, double R)
{
                                                                                                                                          
  double oogamma1 = 1.0/(gamma - 1.0);
  double Cp = R*gamma*oogamma1; // specific heat at constant pressure
                                                                                                                                          
  Zi[0] = Cp * (ratdelta*vc[0]*sqrt2S2*dtdxj[0] - r_s_dtdxj[0]);
  Zi[1] = Cp * (ratdelta*vc[0]*sqrt2S2*dtdxj[1] - r_s_dtdxj[1]);
  Zi[2] = Cp * (ratdelta*vc[0]*sqrt2S2*dtdxj[2] - r_s_dtdxj[2]);
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------
                                                                                                                                          
inline
void Tet::computePij(double dudxj[3][3], double Pij[3][3])
                                                                                                                                          
{
                                                                                                                                          
  Pij[0][0] = (2.0/3.0) * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
  Pij[1][1] = (2.0/3.0) * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
  Pij[2][2] = (2.0/3.0) * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
  Pij[0][1] = dudxj[1][0] + dudxj[0][1];
  Pij[0][2] = dudxj[2][0] + dudxj[0][2];
  Pij[1][2] = dudxj[2][1] + dudxj[1][2];
  Pij[1][0] = Pij[0][1];
  Pij[2][0] = Pij[0][2];
  Pij[2][1] = Pij[1][2];
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------
                                                                                                                                          
inline
void Tet::computeTemp(double *V[4], double t[4], double R)
{
                                                                                                                                          
  double ooR = 1.0/R;
  t[0] = (V[0][4] / V[0][0]) * ooR;
  t[1] = (V[1][4] / V[1][0]) * ooR;
  t[2] = (V[2][4] / V[2][0]) * ooR;
  t[3] = (V[3][4] / V[3][0]) * ooR;
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------

inline
void Tet::computeTempGradient(double dp1dxj[4][3],
                                 double t[4], double dtdxj[3])
                                                                                                                                          
{
                                                                                                                                          
 dtdxj[0] = (dp1dxj[0][0]*t[0] +
             dp1dxj[1][0]*t[1] +
             dp1dxj[2][0]*t[2] +
             dp1dxj[3][0]*t[3]);
                                                                                                                                          
 dtdxj[1] = (dp1dxj[0][1]*t[0]+
             dp1dxj[1][1]*t[1]+
             dp1dxj[2][1]*t[2]+
             dp1dxj[3][1]*t[3]);
                                                                                                                                          
 dtdxj[2] = (dp1dxj[0][2]*t[0]+
             dp1dxj[1][2]*t[1]+
             dp1dxj[2][2]*t[2]+
             dp1dxj[3][2]*t[3]);
}
                                                                                                                                          
//-----------------------------------------------------------------------
                                                                                                                                          
inline
double Tet::computeNormSij(double duidxj[3][3])
{
                                                                                                                                          
  double S[3][3];
                                                                                                                                          
  S[0][0] = duidxj[0][0];
  S[1][1] = duidxj[1][1];
  S[2][2] = duidxj[2][2];
                                                                                                                                          
  S[0][1] = 0.5 * (duidxj[0][1] + duidxj[1][0]);
  S[0][2] = 0.5 * (duidxj[0][2] + duidxj[2][0]);
  S[1][2] = 0.5 * (duidxj[1][2] + duidxj[2][1]);
                                                                                                                                          
  S[1][0] = S[0][1];
  S[2][0] = S[0][2];
  S[2][1] = S[1][2];
                                                                                                                                          
  double S2 = (S[0][0]*S[0][0] + S[0][1]*S[0][1] + S[0][2]*S[0][2] +
               S[1][0]*S[1][0] + S[1][1]*S[1][1] + S[1][2]*S[1][2] +
               S[2][0]*S[2][0] + S[2][1]*S[2][1] + S[2][2]*S[2][2]);
                                                                                                                                          
  return sqrt(2.0 * S2);
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------

inline
void Tet::computeVelocity(double *V[4], double u[4][3],
                         double ucg[3], double ntet[4])
{
                                                                                                                                          
   u[0][0] = V[0][1]/ntet[0];
   u[0][1] = V[0][2]/ntet[0];
   u[0][2] = V[0][3]/ntet[0];
                                                                                                                                          
   u[1][0] = V[1][1]/ntet[1];
   u[1][1] = V[1][2]/ntet[1];
   u[1][2] = V[1][3]/ntet[1];
                                                                                                                                          
   u[2][0] = V[2][1]/ntet[2];
   u[2][1] = V[2][2]/ntet[2];
   u[2][2] = V[2][3]/ntet[2];
                                                                                                                                          
   u[3][0] = V[3][1]/ntet[3];
   u[3][1] = V[3][2]/ntet[3];
   u[3][2] = V[3][3]/ntet[3];
                                                                                                                                          
   ucg[0] = (1.0/4.0) * (u[0][0] + u[1][0] + u[2][0] + u[3][0]);
   ucg[1] = (1.0/4.0) * (u[0][1] + u[1][1] + u[2][1] + u[3][1]);
   ucg[2] = (1.0/4.0) * (u[0][2] + u[1][2] + u[2][2] + u[3][2]);
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------

inline
void Tet::computeVelocityGradient(double dp1dxj[4][3],
                                      double u[4][3],
                                    double dudxj[3][3])
                                                                                                                                          
{
                                                                                                                                          
  dudxj[0][0] = dp1dxj[0][0]*u[0][0] + dp1dxj[1][0]*u[1][0] +
                dp1dxj[2][0]*u[2][0] + dp1dxj[3][0]*u[3][0];
                                                                                                                                          
  dudxj[0][1] = dp1dxj[0][1]*u[0][0] + dp1dxj[1][1]*u[1][0] +
                dp1dxj[2][1]*u[2][0] + dp1dxj[3][1]*u[3][0];
                                                                                                                                          
  dudxj[0][2] = dp1dxj[0][2]*u[0][0] + dp1dxj[1][2]*u[1][0] +
                dp1dxj[2][2]*u[2][0] + dp1dxj[3][2]*u[3][0];
                                                                                                                                          
  dudxj[1][0] = dp1dxj[0][0]*u[0][1] + dp1dxj[1][0]*u[1][1] +
            dp1dxj[2][0]*u[2][1] + dp1dxj[3][0]*u[3][1];
                                                                                                                                          
  dudxj[1][1] = dp1dxj[0][1]*u[0][1] + dp1dxj[1][1]*u[1][1] +
          dp1dxj[2][1]*u[2][1] + dp1dxj[3][1]*u[3][1];
                                                                                                                                          
  dudxj[1][2] = dp1dxj[0][2]*u[0][1] + dp1dxj[1][2]*u[1][1] +
                dp1dxj[2][2]*u[2][1] + dp1dxj[3][2]*u[3][1];
                                                                                                                                          
  dudxj[2][0] = dp1dxj[0][0]*u[0][2] + dp1dxj[1][0]*u[1][2] +
              dp1dxj[2][0]*u[2][2] + dp1dxj[3][0]*u[3][2];
                                                                                                                                          
  dudxj[2][1] = dp1dxj[0][1]*u[0][2] + dp1dxj[1][1]*u[1][2] +
            dp1dxj[2][1]*u[2][2] + dp1dxj[3][1]*u[3][2];
                                                                                                                                          
  dudxj[2][2] = dp1dxj[0][2]*u[0][2] + dp1dxj[1][2]*u[1][2] +
          dp1dxj[2][2]*u[2][2] + dp1dxj[3][2]*u[3][2];
                                                                                                                                          
}
                                                                                                                                          
//-----------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <Tet.C>
#endif

#endif
