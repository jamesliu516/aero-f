#ifndef _FACE_H_
#define _FACE_H_

#include <PostFcn.h>

class VarFcn;
class FluxFcn;
class FemEquationTerm;
class EdgeSet;
class TetSet;
class GeoState;
class BinFileHandler;

struct Vec3D;

template<int dim> class BcData;
template<class Scalar> class Vec;
template<class Scalar, int dim> class GenMat;
template<class Scalar, int dim> class SVec;

//------------------------------------------------------------------------------

class Face {

  static const double third;
  static const int edgeEnd[3][2];

  int code;
  int elemNum;
  int nodeNum[3];
  int edgeNum[3];
  int surface_id;

public:

  Face() {}
  ~Face() {}

  int operator[](int i) { return nodeNum[i]; }

  operator int *() { return nodeNum; }
  int numNodes() { return(3); }
  int nodes(int* d) { d[0] = nodeNum[0]; d[1] = nodeNum[1]; d[2] = nodeNum[2]; return(3); }

  bool operator<(const Face &f) const
  {
    return nodeNum[0] < f.nodeNum[0] ||
      ((f.nodeNum[0] == nodeNum[0]) && 
       nodeNum[1] < f.nodeNum[1]) ||
      ((f.nodeNum[0] == nodeNum[0]) && 
       (f.nodeNum[1] == nodeNum[1]) && 
       nodeNum[2] < f.nodeNum[2]);
  }

  int getCode()  { return code; }
  int getSurfaceID() { return surface_id; }
  void setup(int, int *, int surface_id = 0);
  void setType(int *);
  void setNodeType(int*, int*);
  void setNodeFaceType(int*);
  void setElementNumber(int, int*);
  void tagNodesOnBoundaries(Vec<bool> &);
  void tagEdgesOnBoundaries(Vec<bool> &);
  void reorder();
  void numberEdges(EdgeSet &);

  int getElementNumber() const { return elemNum; }

  static double computeVolume(Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &);
  void computeNormal(SVec<double,3> &, Vec3D &);
  void computeEdgeNormals(SVec<double,3>&, int*, SVec<double,6>&);
  void computeNormalGCL1(SVec<double,3> &, SVec<double,3> &, 
			 SVec<double,3> &, Vec3D &, double &);
  void computeNormalEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
			   Vec3D &, double &);

  template<class NodeMap>
  void renumberNodes(NodeMap &);

  template<int dim>
  void assignFreeStreamValues2(SVec<double,dim> &, SVec<double,dim> &, double *);
  template<int dim>
  void assignFreeStreamValues(double *, double *, double *);

  template<int dim>
  void computeFaceBcValue(SVec<double,dim> &, double *);

  template<int dim1, int dim2>
  void computeNodeBcValue(SVec<double,3> &, double *, SVec<double,dim2> &);

  template<int dim>
  void computeForce(TetSet &, PostFcn *, SVec<double,3> &, Vec<double> &, double *, 
		SVec<double,dim> &, double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
		double *nodalForceWeights, int = 0);

  template<int dim>
  void computeNodalForce(TetSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
		double *, SVec<double,dim> &, double, SVec<double,3> &, double *nodalForceWeights);

  template<int dim>
  void computeNodalHeatPower(TetSet&, PostFcn*, SVec<double,3>&, Vec<double>&, 
			     double*, SVec<double,dim>&, Vec<double>&);

  template<int dim>
  void computeForceAndMoment(TetSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
			     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, 
			     Vec3D &, Vec3D &, double *nodalForceWeights, int = 0);

  template<int dim>
  double computeInterfaceWork(TetSet&, PostFcn*, SVec<double,3>&, Vec<double>&, 
			      double, double*, SVec<double,dim>&, double);

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType, TetSet&, PostFcn *, SVec<double,3> &, 
			     Vec<double> &, double *, SVec<double,dim> &, SVec<double,2> &);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, Vec3D &, double, 
			SVec<double,3> &, SVec<double,dim> &,
                       Vec<double> &, Vec<double> &, double, double, double);
  template<int dim>
  void computeTimeStep(VarFcn *, Vec3D &, double, SVec<double,dim> &,
                       Vec<double> &, double, double, double,Vec<double>&);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, Vec3D &, double, SVec<double,dim> &, 
			       double *, SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, Vec3D &, double, SVec<double,dim> &,
                               double *, Vec<double> &, SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn **, Vec3D &, double, SVec<double,dim> &, 
			       SVec<double,1> &, Vec<double> &);

  template<int dim>
  void computeFiniteVolumeBarTerm(FluxFcn **, Vec3D &, double, SVec<double,dim> &, 
			       double *, SVec<double,dim> &, SVec<double,1> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, Vec3D &, double, SVec<double,dim> &, 
				       double *, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, Vec3D &, double, SVec<double,dim> &,
                                       double *, GenMat<Scalar,neq> &, int*);

	template<int dim, class Scalar, int neq>
	void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal,
                                       double normalVel, SVec<double,dim> &V,
                                       double *Ub, GenMat<Scalar,neq> &A, Vec<double> &Phi);

	template<int dim, class Scalar, int neq>
	void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal,
			                                 double normalVel, SVec<double,dim> &V,
			                                 double *Ub, GenMat<Scalar,neq> &A,
			                                 Vec<double> &Phi, int* nodeType);

  template<int dim, class Scalar>
  void computeJacobianFiniteVolumeTermLS(Vec3D &, double, SVec<double,dim> &, 
				         GenMat<Scalar,1> &);

  template<int dim>
  void computeGalerkinTerm(TetSet &, FemEquationTerm *, SVec<double,3> &, 
			   Vec<double> &, double *, SVec<double,dim> &, SVec<double,dim> &);
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(TetSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, double *,
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeForceDerivs(TetSet &, VarFcn *, SVec<double,3> &, 
                          SVec<double,dim> &,SVec<double,dim> &, Vec<double> &, SVec<double,3> **);

  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, TetSet &, SVec<double,3> &, SVec<double,dim> &,  
                                Vec<double> &, SVec<double, dim> &, double, Vec3D &, Vec3D &, Vec3D &,
				 Vec3D &, double *nodalForceWeights);

  template<int dim>
  void computeFDerivs(TetSet &, VarFcn *, SVec<double,3> &, SVec<double,dim> &, Vec3D (*));

};

//------------------------------------------------------------------------------

class FaceSet {

  int numFaces;
  Face *faces;
	
public:

  FaceSet(int);
  ~FaceSet() { if (faces) delete[] faces; }

  Face &operator[](int i) const { return faces[i]; }

  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       double, double, double,Vec<double> &);
  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
			SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, double, double, double);

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

  template<int dim>
  void computeFiniteVolumeBarTerm(FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, SVec<double,dim> &, SVec<double,1> &);

  // DEBUG
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
	void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                       GeoState &geoState, SVec<double,dim> &V,
												               GenMat<Scalar,neq> &A, Vec<double> &Phi);


	template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
			                                 GeoState &geoState, SVec<double,dim> &V,
			                                 GenMat<Scalar,neq> &A, Vec<double> &Phi,
			                                 int* nodeType);

  template<int dim, class Scalar>
  void computeJacobianFiniteVolumeTermLS(GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,1> &);

  template<int dim>
  void computeGalerkinTerm(TetSet &, FemEquationTerm *, BcData<dim> &, GeoState &, 
			   SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(TetSet &, FemEquationTerm *, BcData<dim> &, 
				   GeoState &, SVec<double,3> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

  int read(BinFileHandler &, int, int (*)[2], int *);

  int size() const { return numFaces; }

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Face.C>
#endif

#endif
