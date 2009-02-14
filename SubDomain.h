#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include <IoData.h>
#include <PostFcn.h>
#include <Node.h>
#include <Edge.h>
#include <Face.h>
#include <Elem.h>
#include <InletNode.h>
#include <DiagMatrix.h>
#include <DistInfo.h>
#include <BCond.h>
#include <TriangulatedSurface.h>
#include <DenseMatrix.h>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

#include <complex>
typedef std::complex<double> bcomp;

class VarFcn;
class PostFcn;
class BcFcn;
class RecFcn;
class FluxFcn;
class FemEquationTerm;
class MacroCellSet;
class VMSLESTerm;
class DynamicVMSTerm;
class DynamicLESTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class SubDTopo;
class TimeData;
class GeoSource;
class GeoState;
class Connectivity;
class Communicator;
class SubMat;
class SubDiagPrec;
class MemoryPool;
class BinFileHandler;
class BCApplier;
class BCData;
class MatchNodeSet; // HB
class LevelSet;
class VolumicForceTerm;
class TriangulatedSurface;
class TimeLowMachPrec;
class EulerStructGhostFluid;

struct V6NodeData;
struct Vec3D;
struct compStruct;
struct ExtrapolationNodeData;

template<int dim, class Scalar> class NodalGrad;
template<int dim> class EdgeGrad;
template<int dim> class Extrapolation;
template<int dim> class BcData;
template<int dim> class ExactRiemannSolver;
template<class Scalar> class Vec;
template<class Scalar> class CommPattern;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class MvpMat;
template<class Scalar, int dim> class SparseMat;
template<class Scalar, int dim> class GenMat;

//------------------------------------------------------------------------------

struct EdgeDef {

  int glLeft, glRight, edgeNum, sign;

  bool operator<(const EdgeDef &e) const
  {
    return ( (glLeft < e.glLeft) || (glLeft == e.glLeft && glRight < e.glRight) );
  }

  void order()
  {
    if (glLeft < glRight) sign = 1;
    else { int tmp = glRight; glRight = glLeft; glLeft = tmp; sign = -1; }
  }

};

//------------------------------------------------------------------------------


class SubDomain {
  int testEdge;

  int locSubNum;
  int clusSubNum;
  int globSubNum;
  int numClusNodes;
  int numNodeRanges;

  char suffix[100];

  NodeSet &nodes;
  EdgeSet  edges;
  FaceSet &faces;
  ElemSet &elems;

  InletNodeSet inletNodes;

  int *locToGlobNodeMap;
  int *locToGlobFaceMap;
  int *locToGlobElemMap;

  int numNeighb;
  int *neighb;
  int *sndChannel;
  int *rcvChannel;
  Connectivity *sharedNodes;
  Connectivity *sharedInletNodes;
  Connectivity** nodesToMCNodes;
  int *numSharedEdges;
  EdgeDef **sharedEdges;
  int (*nodeRanges)[3];
  int *nodeType;
  int *nodeFaceType;

  map<int, int> bcMap;
  int *numBcNodes;

  BCondSet *mmsBCs;
  int *rotOwn;

  int **nodeFlag;
  Connectivity *NodeToNode;
  Connectivity *NodeToElem;
  Connectivity *ElemToElem;

  int **totalNeiData;
  double *gradP[3];

// Included (MB*)
  int numOffDiagEntries;
  double *dGradP[3];

//  TriangulatedSurface *triaSurf;

public:

  SubDomain(int, int, int, int, char *, NodeSet *, FaceSet *, ElemSet *,
	    int, int *, Connectivity *, int *, int *, int *, int, int (*)[3]);
  ~SubDomain();

  // topology
  int *getNodeMap()  { return locToGlobNodeMap; }
  int getGlobSubNum()  { return globSubNum; }
  int getLocSubNum()  { return locSubNum; }
  int numberEdges();

  Connectivity *createElemBasedConnectivity();
  Connectivity *createNodeToElementConnectivity();
  Connectivity *createElementToElementConnectivity();
  Connectivity *createEdgeBasedConnectivity();
  Connectivity *createNodeToMacroCellNodeConnectivity(MacroCellSet *);
  Connectivity *agglomerate(Connectivity &, int, bool *);
  void createSharedInletNodeConnectivity(int);

  compStruct *createRenumbering(Connectivity *, int, int);

  int numNodes() { return(nodes.size()); }
  int numFaces() { return(faces.size()); }
  int numElems() { return(elems.size()); }
  int numEdges() { return(edges.size()); }

  int* getElemNodeNum(int i) {return(elems[i].nodeNum()); }


  // geometry

  void setFaceType(int *);
  void setNodeType(int*, CommPattern<int> &);
  void setNodeFaceType(CommPattern<int> &);
  int* completeNodeType(int*, CommPattern<int> &);
  int* completeNodeFaceType(CommPattern<int> &);
  int setFaceToElementConnectivity();
  void getElementStatistics(int &, int &, int &, int &);
  int computeControlVolumes(int, double, SVec<double,3> &, Vec<double> &);
  void computeFaceNormals(SVec<double,3> &, Vec<Vec3D> &);
  void computeFaceEdgeNormals(SVec<double,3>&, SVec<double,6>&);
  void computeEdgeDihedralAngle(double, SVec<double,6>&, Vec<double>&);
  void propagateInfoAlongEdges(Vec<double>&);
  void computeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			  Vec<Vec3D> &, Vec<double> &, Vec<Vec3D> &, Vec<double> &);
  void computeNormalsConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                            Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel,
                            Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel);
  void computeNormalsEZGCL1(double, SVec<double,3>&, SVec<double,3>&, Vec<Vec3D>&,
			    Vec<double>&, Vec<Vec3D>&, Vec<double>&);
  void computeSmoothedSensor(SVec<double,3>&, Vec<double>&, SVec<double,3>&);
  void computeWeightsLeastSquaresEdgePart(SVec<double,3> &, SVec<double,6> &);
  void computeWeightsLeastSquaresNodePart(SVec<double,6> &);
  void computeWeightsLeastSquaresEdgePart(SVec<double,3> &, Vec<double> &,
					  SVec<int,1> &, SVec<double,6> &);
  void computeWeightsLeastSquaresNodePart(SVec<int,1> &, SVec<double,6> &);
  void computeWeightsGalerkin(SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);
  void computeEdgeWeightsGalerkin(SVec<double,3> &, Vec<double> &, SVec<double,9> &);
#define EDGE_LENGTH
#ifdef EDGE_LENGTH
  bool findTetrahedron(int, int, Vec<int>&, int**, SVec<double,3>&, V6NodeData&, bool=true, double* refLength=0);
#else
  bool findTetrahedron(int, int, Vec<int>&, int**, SVec<double,3>&, V6NodeData&, bool=true);
#endif
  bool findNormalTet1(Vec3D, Vec3D, Vec3D, Vec3D, Vec3D, int, int, ExtrapolationNodeData &, bool);
  void findNormalTet2(int*, Vec3D, Vec3D, int, SVec<double,3> , ExtrapolationNodeData &);
  bool findNormalTetrahedron(int , Vec3D , int *, int *, int ,
                             SVec<double,3> , ExtrapolationNodeData *, bool =true);
  void findEdgeTetrahedra(SVec<double,3>&, V6NodeData (*&)[2]);
  void findNormalTetrahedra(SVec<double,3>& , Vec<Vec3D>& , ExtrapolationNodeData (*&)[2]);
  void checkNormalTetrahedra(ExtrapolationNodeData (*&)[2]);
  void setInletNodes(IoData& );
  void setInletNodes2(IoData& );
  int findOppositeTet(int* , int );
  bool findNodeInTet(int, int);
  void checkInletNodes();
  void sumInletNormals(Vec<Vec3D>&,     Vec<Vec3D>&, Vec<int>&);
  void numDivideNormals(Vec<Vec3D>&, Vec<int>&);
  void getReferenceMeshPosition(SVec<double,3> &);
  void computeDisplacement(SVec<double,3> &, SVec<double,3> &);

  void assimilateCells(int, int, int*, int**, bool *, int, int *, bool *);
  MacroCellSet** findAgglomerateMesh(int, int, bool *, double);
  void createMacroCellConnectivities(MacroCellSet **, int);

  void setBCond(BCondSet *subBC) { mmsBCs = subBC; }
  void applySmoothing(Vec<double> &, Vec<double> &);
  void applySmoothing(Vec<double> &, SVec<double,2> &);
  void computeTetsConnectedToNode(Vec<int> &Ni);
  void computeLocalAvg(SVec<double,3> &, Vec<double> &, Vec<double> &);
  void computeLocalAvg(SVec<double,3> &, SVec<double,2> &, SVec<double,2> &);
  void computeFilterWidth(SVec<double,3> &, Vec<double> &);
  void finalizeTags(SVec<int,2> &);
  void setPhiForFluid1(Vec<double> &);
  void setPhiWithDistanceToGeometry(SVec<double,3> &X, double x,
                                    double y, double z, double r,
                                    double invertGasLiquid,
                                    Vec<double> &Phi);
  void setPhiByGeometricOverwriting(SVec<double,3> &X, double x,
                                    double y, double z, double r,
                                    double invertGasLiquid,
                                    Vec<double> &Phi);
  void setPhiForShockTube(SVec<double,3> &X,
                          double radius, Vec<double> &Phi);
  void setPhiForBubble(SVec<double,3> &X, double x, double y,
                       double z, double radius, double invertGasLiquid,
                       Vec<double> &Phi);
  void setupPhiVolumesInitialConditions(const int volid, Vec<double> &Phi);
  void setupPhiMultiFluidInitialConditionsSphere(SphereData &ic,
                                 SVec<double,3> &X, Vec<double> &Phi);
  void outputCsDynamicLES(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, Vec<double> &);

  void setupPhiMultiFluidInitialConditionsPlane(PlaneData &ip,
                                 SVec<double,3> &X, Vec<double> &Phi);
  template<int dim>
  void setupUVolumesInitialConditions(const int volid, FluidModelData &fm,
             VolumeInitialConditions &ic, SVec<double,dim> &U);
  template<int dim>
  void setupUMultiFluidInitialConditionsSphere(FluidModelData &fm,
             SphereData &ic, SVec<double,3> &X, SVec<double,dim> &U);
  template<int dim>
  void setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
             PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U);

  void computeLij(double [3][3], double [3], double [6], double [5]);
  void computeBij(double [3][3], double [6], double, double [3][3], double, double [5]);
  void computeZi(double [3], double, double, double [3], double [3], double [5], double);
  void computeLi(double [3], double, double, double [3], double [5]);


  // moving mesh

  void getNdAeroLists(int &, int *&, int &, int *&, int &, int *&, MatchNodeSet* matchNodes=0);

  template<class MatScalar, class PrecScalar>
  void computeStiffAndForce(DefoMeshMotionData::Element, SVec<double,3>&, SVec<double,3>&, GenMat<MatScalar,3>&,
                            GenMat<PrecScalar,3>*, double volStiff, int* ndType);

  // spatial discretization
  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
		       Vec<double> &, Vec<double> &,
                       TimeLowMachPrec &);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
		       Vec<double> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<double> &);


  template<int dim, class Scalar>
  void computeGradientsLeastSquares(SVec<double,3> &, SVec<double,6> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(SVec<double,3> &, Vec<double> &,
                                    SVec<double,6> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &);


  template<int dim, class Scalar>
  void computeGradientsGalerkin(Vec<double> &, SVec<double,3> &, SVec<double,3> &,
				SVec<double,3> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsGalerkinT(Vec<double> &, SVec<double,3> &,
                SVec<double,3> &, SVec<double,3> &, SVec<Scalar,dim> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim>
  void computeMinMaxStencilValues(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
//  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, SVec<double,3> &, Vec<double> &,
  void computeMultiDimLimiter(RecLimiter *, SVec<double,3> &, Vec<double> &,

			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &);

  template<int dim>
  void computePressureSensor(SVec<double,3>&, SVec<double,dim>&, SVec<double,dim>&,
			     SVec<double,dim>&, SVec<double,dim>&, SVec<double,3>&);

  template<int dim>
  int computeFiniteVolumeTerm(Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&,
                              EdgeGrad<dim>*, SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, Vec<double> &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              NodalGrad<1>&,
                              SVec<double,dim>&, int, SVec<double,dim> *,
                              SVec<double,dim> *,
                              SVec<int,2>&, int, int);
  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, LevelSetStructure &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              SVec<double,dim>&, int, SVec<int,2>&, int, int);

 template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, BcData<dim>&, GeoState&,
                               SVec<double,3>&, SVec<double,dim>&,
                               NodalGrad<dim>&, NodalGrad<1>&, EdgeGrad<dim>*, SVec<double,1>&,
                               Vec<double>&);
  template<int dim>
  int computeFiniteVolumeBar_Step1(Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&, SVec<double,3>& ,
                                    SVec<double,dim>&, NodalGrad<dim> &, EdgeGrad<dim>* , SVec<double,dim>&,
                                    SVec<int,2> &, int, int);

  template<int dim>
  void computeFiniteVolumeBar_Step2(MacroCellSet **,SVec<double,1> &, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void computeVolumeChangeTerm(Vec<double> &ctrlVol, GeoState &geoState,
                               SVec<double,dim> &U, SVec<double,dim> &R);
  void computeVolumeChangeTerm(Vec<double> &ctrlVol, GeoState &geoState,
                               Vec<double> &Phi, Vec<double> &dPhi);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       Vec<double> &, SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       CommPattern<double> *);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                                       FluxFcn **, BcData<dim> &, GeoState &,
                                       NodalGrad<dim> &, NodalGrad<1> &,
                                       SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       Vec<double> &, CommPattern<double> *);
  template<int dim>
  void recomputeRHS(VarFcn*, SVec<double,dim>& ,SVec<double,dim>& , Extrapolation<dim>*,
                                        BcData<dim>&, GeoState&, SVec<double,3> &);
  template<int dim>
  void recomputeRHS(VarFcn*, SVec<double,dim>& ,Vec<double> &, SVec<double,dim>&,
                    Extrapolation<dim>*, BcData<dim>&, GeoState&, SVec<double,3> &);

  template<int dim>
  void recomputeResidual(SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeRealFluidResidual(SVec<double, dim> &, SVec<double,dim> &, Vec<double> &);


  template<class Scalar,int dim>
  void checkRHS(Scalar (*)[dim]);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
			   SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeVolumicForceTerm(VolumicForceTerm *, Vec<double> &,
                               SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &,
				 SVec<double,dim> &);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeMutOMuSmag(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuVMS(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuDynamicVMS(DynamicVMSTerm *, SVec<double,dim> &, SVec<double,3> &,
                               SVec<double,dim> &, Vec<double> &, Vec<double> &);

  template<int dim>
  void computeMutOMuWale(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuDynamicLES(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &,
                               SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeTestFilterAvgs(SVec<double,dim> &,  SVec<double,16> &, SVec<double,6> &, Vec<double> &,
                             SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double);

  template<int dim>
  void computeCsValues(SVec<double,dim> &,  SVec<double,16> &, SVec<double,6> &, Vec<double> &,
                       SVec<double,8> &, SVec<double,2> &, Vec<int> &, SVec<double,3> &, double, double);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, SVec<double,dim> &,
                             SVec<double,dim> &);

  template<int dim>
  void computeVMSLES_Step1(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeVMSLES_Step2(SVec<double,1> &, MacroCellSet *, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void computeGalerkinBar_Step1(FemEquationTerm *, BcData<dim> &, GeoState &, SVec<double,3> &,
                                SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeGalerkinBar_Step2(MacroCellSet **, SVec<double,1> &, SVec<double,dim> &,
                                SVec<double,dim> &, int);

  template<int dim>
  void computeMBarAndM_Step1(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &,
                                      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeMBarAndM_Step2(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &,
                             SVec<double,dim> &, SVec<double,dim> &, int, int);

  template<int dim>
  void computeCsDeltaSq(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                                   SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                                   Vec<double> &, Vec<double> &, int);

  template<int dim>
  void computeDynamicVMSTerm_Step1(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                                   SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                                   Vec<double> *, Vec<double> &);
  template<int dim>
  void computeDynamicVMSTerm_Step2(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void computedWBar_dt(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
				   SVec<double,3> &, Vec<double> &, SVec<double,dim> &,
				   GenMat<Scalar,neq> &) ;

  template<int dim, class Scalar, int neq>
  void computeJacobianVolumicForceTerm(VolumicForceTerm *, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &) ;

  template<class Scalar, int neq>
  void finishJacobianGalerkinTerm(Vec<double> &, GenMat<Scalar,neq> &) ;

  template<int dim>
  void getExtrapolationValue(Extrapolation<dim>*,SVec<double,dim> &, SVec<double,dim> &, VarFcn*,
                                         BcData<dim>&, GeoState&, SVec<double,3>&);

  template<int dim>
  void applyExtrapolationToSolutionVector(Extrapolation<dim>*, SVec<double,dim> &,
					  SVec<double,dim> &);
  template<int dim>
  void applyBCsToSolutionVector(BcFcn *, BcData<dim> &, SVec<double,dim> &);

  template<int dim>
  void applyBCsToResidual(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToH2Jacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int dim>
  SparseMat<Scalar,dim> *createMaskJacobian(int *, MemoryPool *);

  template<class Scalar, int dim>
  MvpMat<Scalar,dim> *createMaskMatVecProd(bool flag = false);

  template<class Scalar, int dim>
  DiagMat<Scalar,dim> *createMaskDiagonal(typename DiagMat<Scalar,dim>::Type, int *);

  template<class Scalar, int dim>
  SparseMat<Scalar,dim> *createMaskILU(int, int, int *);

  template<class Scalar, int dim>
  void computeH1(FluxFcn **, BcData<dim> &, GeoState &, Vec<double> &,
                 SVec<double,dim> &, GenMat<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeH2(FluxFcn **, RecFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
		 SVec<double,dim> &, NodalGrad<dim> &, GenMat<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeH2LS(GeoState &, SVec<double,3> &, SVec<double,dim> &,
                          NodalGrad<dim> &, GenMat<Scalar,1> &);

  template<class Scalar, int dim>
  void precomputeRec(RecFcn *, SVec<double,3> &, SVec<double,dim> &,
		     NodalGrad<dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
		     SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<class Scalar, int dim>
  void computeMatVecProdH1(bool *, GenMat<Scalar,dim> &, SVec<double,dim> &,
			   SVec<double,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(RecFcn *, SVec<double,3> &, Vec<double> &,
			   GenMat<Scalar1,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			   SVec<double,dim> &, SVec<double,dim> &, SVec<Scalar2,dim> &,
			   NodalGrad<dim, Scalar2> &, SVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2>
  void computeMatVecProdH2LS(RecFcn *, SVec<double,3> &, Vec<double> &,
			   GenMat<Scalar1,1> &, SVec<double,1> &, SVec<double,1> &,
			   SVec<double,1> &, SVec<double,1> &, Vec<Scalar2> &,
			   Vec<Scalar2> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2T(RecFcn *, SVec<double,3> &, Vec<double> &,
                GenMat<Scalar1,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                SVec<double,dim> &, SVec<double,dim> &, SVec<Scalar2,dim> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2Tb(RecFcn *, SVec<double,3> &, Vec<double> &,
                GenMat<Scalar1,dim> &, NodalGrad<dim, Scalar2> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);


  // I/O

  template<class Scalar, int dim>
  double readTagFromFile(const char *, int, int *, int *);

  template<class Scalar, int dim>
  void openFileForWriting(const char *, int);

  template<class Scalar, int dim>
  void writeTagToFile(const char *, int, double);

  template<class Scalar, int dim>
  void readVectorFromFile(const char *, int, int, SVec<Scalar,dim> &, Scalar *);

  template<class Scalar, int dim>
  void writeVectorToFile(const char *, int, SVec<Scalar,dim> &, Scalar *);

  template<int dim>
  void assignFreeStreamValues2(SVec<double,dim> &, SVec<double,dim> &,
			       SVec<double,dim> &, SVec<double,dim> &);
  template<int dim>
  void assignFreeStreamValues(double *, double *, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void setNodeBcValue(double*, SVec<double,dim>&);

  template<int dim>
  void computeFaceBcValue(SVec<double,dim> &, SVec<double,dim> &);

  template<int dim1, int dim2>
  void computeNodeBcValue(SVec<double,3> &, SVec<double,dim1> &, SVec<double,dim2> &);

  template<int dim>
  void computeNodalForce(PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
		SVec<double,dim> &, Vec<double> &, SVec<double,3> &);

  template<int dim>
  void computeNodalHeatPower(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&,
			     SVec<double,dim>&, Vec<double>&);

  template<int dim>
  void computeForceAndMoment(map<int,int> &surfIndexMap, PostFcn *,
                             BcData<dim> &, GeoState &, SVec<double,3> &,
			     SVec<double,dim> &, Vec3D &, Vec3D *, Vec3D *,
			     Vec3D *, Vec3D *, int ,
                             SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX = 0,
                             Vec<double> *genCF = 0);

  template<int dim>
  double computeInterfaceWork(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&,
			      SVec<double,dim>&, Vec<double>&);
  template<int dim>
  void computeFaceScalarQuantity(PostFcn::ScalarType, PostFcn *, BcData<dim> &, GeoState &,
				 SVec<double,3> &, SVec<double,dim> &, SVec<double,2> &);

  template<int dim>
  void computeNodeScalarQuantity(PostFcn::ScalarType, PostFcn *,
				 SVec<double,dim> &, SVec<double,3> &, Vec<double> &);
  template<int dim>
  void computeXP(PostFcn *, SVec<double,dim> &V, SVec<double,3> &X, Vec<double> &XP, int);

  template<int dim>
  void computeNodeScalarQuantity(PostFcn::ScalarType, PostFcn *,
                                 SVec<double,dim> &, SVec<double,3> &,
				 Vec<double> &, Vec<double> &);

  template<int dim>
  void computeForceDerivs(VarFcn *, SVec<double,3> &, SVec<double,dim> &,
                          SVec<double,dim> &, Vec<double> &, SVec<double, 3> **);

  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, GeoState &, BcData<dim> &, SVec<double,3> &,
                                SVec<double,dim> &, double, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
                                VecSet< SVec<double,3> > *mX = 0 , Vec<double> *genCF = 0);

  // communication

  void markLenNodes(DistInfo &distInfo) { distInfo.setLen(locSubNum, nodes.size()); }
  void markLenEdges(DistInfo &distInfo) { distInfo.setLen(locSubNum, edges.size()); }
  void markLenFaces(DistInfo &distInfo) { distInfo.setLen(locSubNum, faces.size()); }
  void markLenFaceNorms(DistInfo &distInfo) { distInfo.setLen(locSubNum, faces.sizeNorms()); }
  void markLenInletNodes(DistInfo &distInfo) {distInfo.setLen(locSubNum, inletNodes.size()); }
  void markLenNull(DistInfo &distInfo) {distInfo.setLen(locSubNum, 0); }
  void makeMasterFlag(DistInfo &);
  void setChannelNums(SubDTopo &);
  void identifyEdges(CommPattern<int> &);
  void sndNormals(CommPattern<double> &, Vec3D *, double *);
  void rcvNormals(CommPattern<double> &, Vec3D *, double *);
  void sndEdgeInfo(CommPattern<int> &);
  void rcvEdgeInfo(CommPattern<int> &);

  template<class Scalar>
  void setComLenNodes(int dim, CommPattern<Scalar> &);

  template<class Scalar>
  void setComLenInletNodes(int dim, CommPattern<Scalar> &);

  template<class Scalar>
  void setComLenEdges(int dim, CommPattern<Scalar> &);

  template<class Scalar, int dim>
  void sndData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void sndInletData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvInletData(CommPattern<Scalar> &, Scalar (*)[dim], bool = false);

  template<class Scalar, int dim>
  void sndInletRhsData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvInletRhsData(CommPattern<Scalar> &, Scalar(*)[dim]);

  template<class Scalar, int dim>
  void minRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void maxRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar1, class Scalar2, int dim1, int dim2>
  void TagPsiExchangeData(CommPattern<Scalar1> &splevel, Scalar1 (*level)[dim1],
                          CommPattern<Scalar2> &sppsi, Scalar2 (*psi)[dim2]);

  template<class Scalar, int dim>
  void sndDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void sndDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void checkDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void sndEdgeData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvEdgeData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void sndOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  // test

  void testNormals(Vec<Vec3D> &, Vec<double> &, Vec<Vec3D> &, Vec<double> &);

  template<int dim>
  int checkSolution(VarFcn *, SVec<double,dim> &);

  template<int dim>
  int checkSolution(VarFcn *, Vec<double> &, SVec<double,dim> &, Vec<double> &, Vec<double> &);

  template<int dim>
  void restrictionOnPhi(SVec<double,dim> &initial, Vec<double> &Phi,
                        SVec<double,dim> &restriction, int sign);

  template<int dim>
  void checkFailSafe(VarFcn*, SVec<double,dim>&, SVec<bool,2>&, Vec<double> * = 0);

  template<int dim, int neq>
  int clipSolution(TsData::Clipping, BcsWallData::Integration, VarFcn*,
		   double*, bool*, SVec<double,dim>&, int*, int*, double*);

  template<int dim>
  void checkGradientsSetUp(SVec<double,3> &, SVec<double,dim> &);

  template<int dim>
  void checkGradientsWrite(SVec<double,3> &, NodalGrad<dim> &);

  template<int dim>
  void checkMatVecProd(SVec<double,dim> &, const char *);

#ifdef CXFS
  void printInfo(FILE *);
#endif

  template<int dim>
  void zeroInternalVals(SVec<double, dim> &);

  int *getRotSurfaceOwnership(CommPattern<int> &, map<int,SurfaceData *> &surfaceMap);
  void completeRotateSurfaceOwnership(CommPattern<int> &);

  int *getSlipSurfOwnership(CommPattern<int> &, map<int,SurfaceData *> &surfaceMap);

  void createSlipSurfProjection(int *surfaceOwn, CommPattern<int> &, BCApplier *,
                                SurfaceData *slipSurfaces[]);

  int* getMeshMotionDofType(map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP, MatchNodeSet* matchNodes=0 );

  void completeMeshMotionDofType(int* DofType, CommPattern<int> &ntP);

  template<int dim>
  void zeroMeshMotionBCDofs(SVec<double,dim> &x, int* DofType);

  int *getRotOwn() { return rotOwn; }
  NodeSet &getNodes() { return nodes; }

  void TagInterfaceNodes(Vec<int> &Tag, Vec<double> &Phi, int level);
  void FinishReinitialization(Vec<int> &Tag, SVec<double,1> &Psi, int level);

	template<int dim>
  void storePrimitive(SVec<double,dim> &Vg, SVec<double,dim> &Vgf,
                      Vec<double> &weight, Vec<double> &Phi);

  template<int dim>
  void storeGhost(SVec<double,dim> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computePsiResidual(SVec<double,3> &X, NodalGrad<dim> &grad,
                          Vec<double> &Phi, SVec<double,dim> &Psi,
                          Vec<int> &Tag,
                          Vec<double> &w, Vec<double> &beta,
                          SVec<double,dim> &PsiRes, int typeTracking);
  template<int dim>
  void computePsiResidual2(Vec<int> &Tag, Vec<double> &w, Vec<double> &beta,
                           SVec<double,dim> &PsiRes);
  template<int dim>
  void computePsiResidual3(double bmax, Vec<int> &Tag, Vec<double> &w, Vec<double> &beta,
                           SVec<double,dim> &PsiRes,bool localdt);
  template<int dim>
  void copyCloseNodes(int level, Vec<int> &Tag,Vec<double> &Phi,SVec<double,dim> &Psi);
  template<int dim>
  void computeDistanceCloseNodes(Vec<int> &Tag, SVec<double,3> &X,
                                 NodalGrad<dim> &grad,
                                 Vec<double> &Phi,SVec<double,dim> &Psi);
  template<int dim>
  void recomputeDistanceCloseNodes(Vec<int> &Tag, SVec<double,3> &X,
                                 NodalGrad<dim> &grad,
                                 Vec<double> &Phi,SVec<double,dim> &Psi);
  template<int dim>
  double computeDistanceLevelNodes(Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,dim> &Psi, Vec<double> &Phi);

  template<int dim>
  void checkNodePhaseChange(SVec<double,dim> &X);
  template<int dim>
  void getSignedDistance(SVec<double,dim> &Psi, Vec<double> &Phi);
  template<int dim>
  void checkWeights(Vec<double> &Phi, Vec<double> &Phin,
                    SVec<double,dim> &Update, Vec<double> &Weight);

  template<int dim>
  void checkExtrapolationValue(SVec<double,dim>&,  VarFcn*,
                               BcData<dim>&, GeoState&);

  template<int dim>
  void printVariable(SVec<double,dim>&, VarFcn *vf);

  template<int dim>
  void printInletVariable(SVec<double,dim>&);

  template<int dim>
  void printAllVariable(Vec<int> &, SVec<double,dim>&, int , int);

  void printPhi(SVec<double,3> &, Vec<double>&, int );

  template<class Scalar, int neq>
  void printAllMatrix(GenMat<Scalar,neq> &, int );

  template<int dim>
  void hardyInterpolationLogMap(SVec<double, dim> ***, SVec<double, dim> **, int, int, int, FullM &, FullM &);

  template<int dim>
   void padeReconstruction(SVec<double, dim> **, SVec<double, dim> **, int *, double *, double, int, int, int, int );
  void buildPadeMatrix(bcomp *, int *, int, double *, bcomp *, int, int, int );

  void buildPadeRhs(bcomp *, int *, int, bcomp *, int, int, int );
  void solveLinearSystem(bcomp *, bcomp *, int );

  void padeSolution(bcomp *, int *, double , int , int , int , double , int , bcomp *, int* , int , double , double );

  void buildDeltaFreq(double *, int , double *, int *);

  void extractElementsRelativeToAComponentAndAMode(double* , bcomp* , int , int , int , int , int , int , double );

  template<int dim>
  void extractElementsRelativeToANode(SVec<double, dim> **, double *, int , int );

  template<int dim>
  void extractElementsRelativeToANodeAndAVector(SVec<double, dim> ***, double *, int , int , int , int );

  template<int dim>
  void snapshotsConstruction(SVec<double, dim> **, bcomp* , int , int , int, int , int , double );

  void multiPointsDeltaFreq(int, double*, int , double* , int *);

  int multiPointsFreq(int , int , double *, int , int );

  void multiPade(bcomp *, int *, double *, bcomp *, bcomp *, int , int , int , double , double , bcomp *, double *);
// Included (MB)
  int computeDerivativeOfControlVolumes(int, double, SVec<double,3> &, SVec<double,3> &, Vec<double> &);

  void computeDerivativeOfNormals(SVec<double,3> &, SVec<double,3> &, Vec<Vec3D> &,
                        Vec<Vec3D> &, Vec<double> &, Vec<double> &, Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &);

  void computeDerivativeOfWeightsLeastSquaresEdgePart(SVec<double,3> &, SVec<double,3> &, SVec<double,6> &, SVec<double,6> &);

  void computeDerivativeOfWeightsLeastSquaresNodePart(SVec<double,6> &, SVec<double,6> &);

  void computeDerivativeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares(SVec<double,3> &, SVec<double,3> &,
                    SVec<double,6> &, SVec<double,6> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsGalerkin(Vec<double> &, Vec<double> &, SVec<double,3> &, SVec<double,3> &,
                SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfMinMaxStencilValues(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfMultiDimLimiter(RecLimiter *, SVec<double,3> &, SVec<double,3> &, Vec<double> &, Vec<double> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                  SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(Vec<double> &, Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
			       SVec<double,3>&, SVec<double,3>&, SVec<double,dim>&, SVec<double,dim>&,
			       NodalGrad<dim>&, EdgeGrad<dim>*, double, SVec<double,dim>&);

  template<int dim1, int dim2>
  void computeDerivativeOfNodeBcValue(SVec<double,3> &, SVec<double,3> &, SVec<double,dim1> &, SVec<double,dim1> &, SVec<double,dim2> &);

  template<int dim>
  void computeDerivativeOfNodalForce(PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &, SVec<double,3> &,
                                                                 SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, double [3],
                                                                 SVec<double,3> &);

  template<int dim>
  void computeDerivativeOfForceAndMoment(map<int,int> &surfIndexMap, PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &, SVec<double,3> &,
                                                                           SVec<double,dim> &, SVec<double,dim> &, double [3],
                                                                           Vec3D &, Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
			   SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void applyBCsToDerivativeOfResidual(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfNodeScalarQuantity(PostFcn::ScalarDerivativeType, PostFcn *, double [3], SVec<double,dim> &, SVec<double,dim> &, SVec<double,3> &, SVec<double,3> &, Vec<double> &);

  template<int dim, class Scalar>
  void applyBCsToH2Jacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobianWallValues(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeBCsJacobianWallValues(FemEquationTerm *, BcData<dim> &, GeoState &,
				   SVec<double,3> &, SVec<double,dim> &) ;

  template<int dim>
  void computeNodeBCsWallValues(SVec<double,3> &, SVec<double,1> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar>
  void applyBCsToProduct(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfNodalHeatPower(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&, SVec<double,3>&,
			     SVec<double,dim>&, SVec<double,dim>&, double [3], Vec<double>&);

  template<int dim>
  void computeDerivativeOfVolumicForceTerm(VolumicForceTerm *, Vec<double> &, Vec<double> &,
                               SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
                                SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
                                Vec<double> &, Vec<double> &, double, TimeLowMachPrec &);

  void checkVec(SVec<double,3> &);

  template<int dim>
  int fixSolution(VarFcn *, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void getGradP(NodalGrad<dim>&);

  template<int dim>
  void getDerivativeOfGradP(NodalGrad<dim>&);

//  void getTriangulatedSurfaceFromFace( SVec<double,3> &);

  void getTriangulatedSurfaceFromFace( TriangulatedSurface* );

  void getTriangulatedSurfaceFromFace( SVec<double,3> &,  TriangulatedSurface* );

//  void printTriangulatedSurface();

  void getGhostNodes(double*, int*, int&);

  bool isINodeinITet(Vec3D, int, SVec<double,3>&);

  void  localCoord(Vec3D, int, SVec<double, 3>&, Vec3D&);

  int* getNeiElemOfNode(int, int, int&);

  void getNodeCoords(int, SVec<double,3> &, double&, double&, double&);

  void computeCharacteristicEdgeLength(SVec<double,3>&, double&, double&, double&, int&, const double, const double, const double, const double, const double, const double);

  double specifyBandwidth(Vec<double> &);

  double scalarNormalExtrap(double*, Vec3D, Vec3D, int, SVec<double,3> &, bool);

  bool insideOutside(double*, const double, const double, const double, const double,
                     const double, const double);

  double getMeshInBoundingBox(SVec<double,3> &, const double, const double,
                              const double, const double,
                              const double, const double,
                              int*, int&, int*,
                              int&, int (*)[4]);

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SubDomain.C>
#endif

#endif
