#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <IoData.h>
#include <DistInfo.h>
#include <Timer.h>
#include <VectorSet.h>
#include <Vector.h>
#include <complex.h>
typedef complex<double> bcomp;

class VarFcn;
class BcFcn;
class RecFcn;
class FemEquationTerm;
class VolumicForceTerm;
class FluxFcn;
class TimeData;
class SubDTopo;
class SubDomain;
class GeoSource;
class DistGeoState;
class DistMacroCellSet;
class DynamicVMSTerm;
class VMSLESTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class DynamicLESTerm;
class ViscoFcn;
class PostFcn;
class LevelSet;

class BCApplier; //HB
class MatchNodeSet;

struct Vec3D;

template<int dim> class RecFcnLtdMultiDim;
template<int dim> class DistEdgeGrad;
template<int dim> class DistExtrapolation;
template<int dim> class DistBcData;
template<class T> class CommPattern;
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<class Scalar, int dim> class DistMat;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif


//------------------------------------------------------------------------------

class Domain {

  int numLocSub;

  SubDomain **subDomain;
  SubDTopo *subTopo;

  int **nodeType;

  int **nodeFaceType;

  DistInfo *nodeDistInfo;
  DistInfo *edgeDistInfo;
  DistInfo *faceDistInfo;
  DistInfo *faceNormDistInfo;
  DistInfo *inletNodeDistInfo;

  CommPattern<double> *vecPat;
  CommPattern<bcomp> *compVecPat;
  CommPattern<double> *vec3DPat;
  CommPattern<double> *volPat;
  CommPattern<int> *levelPat;

  CommPattern<double> *weightPat;
  CommPattern<double> *edgePat;
  CommPattern<double> *momPat;
  CommPattern<double> *csPat;
  CommPattern<int> *fsPat;

  CommPattern<double> *inletVec3DPat;
  CommPattern<int> *inletCountPat;
  CommPattern<double> *inletRhsPat;

  Communicator *com;
  Communicator *strCom;
  Communicator *heatCom;
  Communicator *globCom;

  Timer *timer;
  Timer *strTimer;
  Timer *heatTimer;

  DistVec<double> *Delta;
  DistVec<double> *CsDelSq;
  DistVec<double> *PrT;
  DistVec<double> *WCsDelSq;
  DistVec<double> *WPrT;

  DistSVec<int,2> *tag;
  DistSVec<int,2> *tagBar;
                                                                                                                          
  BCApplier* meshMotionBCs; //HB

// Included (MB)
  CommPattern<double> *weightDerivativePat;

public:

  Domain();
  ~Domain();

  int getNumLocSub() const { return numLocSub; }
  SubDomain **getSubDomain() const { return subDomain; }
  SubDTopo *getSubTopo() const { return subTopo; }
  int **getNodeType() const { return nodeType; }
  int **getNodeFaceType() const { return nodeFaceType; }
  int **getNodeTypeExtrapolation() const{
    if (inletRhsPat) return nodeType;
    int** empty = 0;
    return empty;
  }
  BCApplier* getMeshMotionBCs() const { return meshMotionBCs; } //HB
  CommPattern<double> *getVecPat() const { return vecPat; }
  CommPattern<bcomp> *getCompVecPat() const { return compVecPat; }
  CommPattern<double> *getVec3DPat() const { return vec3DPat; }
  CommPattern<double> *getVolPat() const { return volPat; }
  CommPattern<double> *getWeightPat() const { return weightPat; }
  CommPattern<double> *getMomPat() const { return momPat; }
  CommPattern<double> *getCsPat() const { return csPat; }
 

  template<int dim>
  CommPattern<double> *getCommPat(DistSVec<double,dim> &vec) { return vecPat; }
  template<int dim>
  CommPattern<bcomp> *getCommPat(DistSVec<bcomp,dim> &vec) { return compVecPat; }
  CommPattern<double> *getCommPat(DistVec<double> &vec) { return vecPat; }

  Communicator *getCommunicator() const { return com; }
  Communicator *getStrCommunicator() { return strCom; }
  Communicator *getHeatCommunicator() { return heatCom; }
  Timer *getTimer() const { return timer; }
  Timer *getStrTimer() const { return strTimer; }
  Timer *getHeatTimer() const { return heatTimer; }

  DistInfo &getNodeDistInfo() const { return *nodeDistInfo; }
  DistInfo &getEdgeDistInfo() const { return *edgeDistInfo; }
  DistInfo &getFaceDistInfo() const { return *faceDistInfo; }
  DistInfo &getFaceNormDistInfo() const { return *faceNormDistInfo; }
  DistInfo &getInletNodeDistInfo() const { return *inletNodeDistInfo; }

  void getGeometry(GeoSource &, IoData&);
  void createRhsPat(int, IoData&);
  void createVecPat(int, IoData * = 0);
  void numberEdges();
  void setNodeType(IoData &);
  void setInletNodes(IoData &);
  void makeRotationOwnership(IoData &);
  void setFaceToElementConnectivity();
  void printElementStatistics();
  template<int dim>  
  void storeGhost(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &);
  int computeControlVolumes(double, DistSVec<double,3> &, DistVec<double> &);
  void computeFaceNormals(DistSVec<double,3> &, DistVec<Vec3D> &);
  void computeNormalsGCL1(DistSVec<double,3> &, DistSVec<double,3> &, 
			  DistSVec<double,3> &, DistVec<Vec3D> &, DistVec<double> &,
			  DistVec<Vec3D> &, DistVec<double> &);
  void computeNormalsGCL2(TimeData &, DistVec<Vec3D> &, DistVec<Vec3D> &, 
			  DistVec<double> &, DistVec<double> &, DistVec<Vec3D> &, 
			  DistVec<Vec3D> &, DistVec<double> &, DistVec<double> &);
  void computeNormalsEZGCL1(double, DistSVec<double,3>&, DistSVec<double,3>&, DistVec<Vec3D>&, 
			    DistVec<double>&, DistVec<Vec3D>&, DistVec<double>&);
  void computeNormalsEZGCL2(TimeData&, DistVec<double>&, DistVec<double>&, 
			    DistVec<double>&, DistVec<double>&);
  void computeNormalsEZGCL3(TimeData&, DistVec<double>&, DistVec<double>&, DistVec<double>&, 
			    DistVec<double>&, DistVec<double>&, DistVec<double>&);
  void computeInletNormals(DistVec<Vec3D>&, DistVec<Vec3D>&, DistVec<int> &);
  void computeVelocities(ImplicitData::Velocities, TimeData &, DistSVec<double,3> &, 
			 DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &, 
			 DistSVec<double,3> &);
  void computeWeightsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &);
  void computeWeightsLeastSquares(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,6> &);
  void computeWeightsGalerkin(DistSVec<double,3> &, DistSVec<double,3> &, 
			      DistSVec<double,3> &, DistSVec<double,3> &);
  void getReferenceMeshPosition(DistSVec<double,3> &);
  void getNdAeroLists(int *&, int **&, int *&, int **&, int *&, int **&, MatchNodeSet** matchNodes=0);

  void computeDelRatios(DistMacroCellSet *, DistVec<double> &, int, double *, double *, double *, double *, int *);
  void applySmoothing(DistVec<double> &, DistVec<double> &);
                                                                                                                          
  template<class MatScalar, class PrecScalar>
  void computeStiffAndForce(DefoMeshMotionData::Element, DistSVec<double,3>&, 
			    DistSVec<double,3>&, DistMat<MatScalar,3>&, 
			    DistMat<PrecScalar,3>*, double volStiff, int** ndType);

  template<int dim>  
  void computeTimeStep(double, double, FemEquationTerm *, VarFcn *, DistGeoState &, 
		       DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &, 
		       DistVec<double> &, DistVec<double> &, double, double, double, double);

  template<int dim>
  void computeTimeStep(double, double, FemEquationTerm *, VarFcn *, DistGeoState &, DistVec<double> &,
                       DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &,
		       DistVec<double> &, double, double, double, 
		       DistVec<double> &);


  template<int dim, class Scalar>
  void computeGradientsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &, 
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, 
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(DistSVec<double,3> &, DistVec<double> &,
                                    DistSVec<double,6> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsGalerkin(DistVec<double> &, DistSVec<double,3> &, 
				DistSVec<double,3> &, DistSVec<double,3> &, 
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsGalerkinT(DistVec<double> &, DistSVec<double,3> &,
                DistSVec<double,3> &, DistSVec<double,3> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim>
  void computePointWiseSourceTerm(DistGeoState &geoState, DistVec<double> &ctrlVol, 
				  DistNodalGrad<dim> &dVdxj, DistSVec<double,dim> &VV, 
				  DistSVec<double,dim> &RR);
  template<int dim>
  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &, 
			      DistVec<double> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, 
			      DistSVec<double,dim> &, DistSVec<double,dim> &, 
			      DistSVec<double,dim> &, DistSVec<double,dim> &);
  template<int dim>
  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &, 
			      DistVec<double> &, DistSVec<bcomp,dim> &,
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &, 
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &, 
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &)  {
       cout << "computeMultiDimLimiter not implemented for complex data" << endl; }

  template<int dim>
  void computePressureSensor(double, DistSVec<double,3>&,
                             DistSVec<double,dim>&, DistSVec<double,dim>&,
                             DistSVec<double,dim>&, DistSVec<double,dim>&,
                             DistSVec<double,3>&, DistVec<double>&);
  template<int dim>
  void computePressureSensor(double, DistSVec<double,3>&,
                             DistSVec<bcomp,dim>&, DistSVec<bcomp,dim>&,
                             DistSVec<bcomp,dim>&, DistSVec<bcomp,dim>&,
                             DistSVec<bcomp,3>&, DistVec<bcomp>&) { cout << "computePressureSensor not implemented for complex operations" <<endl; }

  template<int dim>
  double reinitLS(DistSVec<double,3> &X, DistVec<double> &Phi, DistSVec<double,dim> &U, int iti);
                                                                                                              
  void solveLSequation(DistSVec<double,3> &X, DistVec<double> &Phi, DistSVec<double,6> &ddx,DistSVec<double,6> &ddy, DistSVec<double,6> &ddz, DistVec<double> &PhiF);

  void solveLS(DistVec<double> &b, DistVec<double> &dPhi, double);
 
  template<int dim>
  void computeFiniteVolumeTerm(DistVec<double> &, DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&, 
			       DistSVec<double,3>&, DistSVec<double,dim>&, 
			       DistNodalGrad<dim>&, DistEdgeGrad<dim>*, 
			       DistSVec<double,dim>&, int, int);

  template<int dim>
  void computeFiniteVolumeTerm(DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&,
                               DistSVec<double,3>&, DistSVec<double,dim>&,
                               DistVec<double> &,
                               DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
                               DistSVec<double,dim>&, int, int);

  template<int dim>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, DistBcData<dim>&, DistGeoState&,
                               DistSVec<double,3>&, DistSVec<double,dim>&,
                               DistNodalGrad<dim>&, DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
                               DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeBarTerm(DistVec<double> &, DistVec<double> &, FluxFcn**, 
                                  RecFcn*, DistBcData<dim>&, DistGeoState&, DistSVec<double,3>&,
                                  DistMacroCellSet *, DistSVec<double,dim> &, DistSVec<double,1> &,
                                  DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
                                  DistSVec<double,dim>&, int, int, int, int);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, DistBcData<dim> &, DistGeoState &, 
				       DistVec<double> &,
				       DistVec<double> &, DistSVec<double,dim> &, 
				       DistMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, DistBcData<dim> &, DistGeoState &,
                                       DistVec<double> &, DistSVec<double,dim> &,
                                       DistMat<Scalar,neq> &, DistVec<double> &);
  template<int dim>
  void recomputeRHS(VarFcn*, DistSVec<double,dim> &, DistSVec<double,dim> &, DistExtrapolation<dim>*,
                   DistBcData<dim>&, DistGeoState&, DistSVec<double,3> &);
  template<int dim>
  void recomputeRHS(VarFcn*, DistSVec<double,dim> &, DistVec<double> &,
                   DistSVec<double,dim> &, DistExtrapolation<dim>*,
                   DistBcData<dim>&, DistGeoState&, DistSVec<double,3> &);
  template<int dim>
  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, DistBcData<dim> &, 
			   DistGeoState &, DistSVec<double,3> &, 
			   DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeVolumicForceTerm(VolumicForceTerm *, DistVec<double> &,
                               DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, DistSVec<double,3> &,
				 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeMutOMuSmag(SmagorinskyLESTerm *, DistVec<double> &, DistSVec<double,3> &,
                               DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void computeMutOMuWale(WaleLESTerm *, DistVec<double> &, DistSVec<double,3> &,
                         DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void computeGalerkinBarTerm(bool, FemEquationTerm *, DistBcData<dim> &, DistGeoState &, DistSVec<double,3> &,
                              DistMacroCellSet *, DistSVec<double,dim> &, DistSVec<double,1> &,
                              DistSVec<double,dim> &, int, int);
                                                                                                                          
  template<int dim>
  void computeVBar(DistMacroCellSet *, bool, DistGeoState &, DistSVec<double,dim> &,
                   DistSVec<double,dim> &, int, int);
                                                                                                                          
  template<int dim>
  void computeBarTerm(DistMacroCellSet *, bool, DistVec<double> &, DistSVec<double,dim> **,
                      DistSVec<double,1> **, DistSVec<double,3> &, DistSVec<double,dim> &, int, int);
                                                                                                                          
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, DistMacroCellSet *, bool, DistVec<double> &,
                             DistSVec<double,dim> **, DistSVec<double,1> **, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &,
                             int, int, int, DistVec<double> * = 0);
                                                                                                                          
  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, DistMacroCellSet *, bool,
                           DistVec<double> &, DistSVec<double,dim> &,
                           DistSVec<double,1> &, DistSVec<double,3> &,
                           DistSVec<double,dim> &, DistSVec<double,dim> &, int);

  template<int dim>
  void computeTestFilterValues(DistSVec<double,dim> &,
                  DistSVec<double,16> &, DistSVec<double,6> &,
                  DistSVec<double,2> &, DistVec<double> &, DistSVec<double,3> &, DistSVec<double,dim> &, double, double);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, DistSVec<double,2> &, DistVec<double> &,
                             DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  void computeDynamicLESTerm(DynamicLESTerm *, DistSVec<double,2> &,
                             DistSVec<double,3> &, DistVec<double> &, DistVec<double> &); // Note Function Overloading

  template<int dim>
  void assemble_dWdt(DistSVec<double, dim> &, DistSVec<double, dim> &);
                                                                                                                          
  template<int dim>
  void computedWBar_dt(DistSVec<double, dim> &, DistSVec<double, dim> &, DistMacroCellSet *, DistSVec<double,1> **, int);
                                                                                                                          
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, DistBcData<dim> &, 
				   DistGeoState &, DistSVec<double,3> &, 
				   DistVec<double> &, DistSVec<double,dim> &,
				   DistMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianVolumicForceTerm(VolumicForceTerm *, DistVec<double> &,
                                       DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void finishJacobianGalerkinTerm(DistVec<double> &, DistMat<Scalar,neq> &);

  template<int dim>
  void getExtrapolationValue(DistExtrapolation<dim> *, DistSVec<double,dim>&, DistSVec<double,dim>&,
                             VarFcn*, DistBcData<dim>&, DistGeoState&, DistSVec<double,3>&);

  template<int dim>
  void applyExtrapolationToSolutionVector(DistExtrapolation<dim> *, DistSVec<double,dim>&,
                                          DistSVec<double,dim>&);

  template<int dim>
  void applyBCsToSolutionVector(BcFcn *, DistBcData<dim> &, DistSVec<double,dim> &);

  template<int dim>
  void applyBCsToResidual(BcFcn *, DistBcData<dim> &, 
			  DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobian(BcFcn *, DistBcData<dim> &, 
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToH2Jacobian(BcFcn *, DistBcData<dim> &, 
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int dim>
  void computeH1(FluxFcn **, DistBcData<dim> &, DistGeoState &,
                 DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeH2(FluxFcn **, RecFcn *, DistBcData<dim> &, DistGeoState &, 
		 DistSVec<double,3> &, DistSVec<double,dim> &, DistNodalGrad<dim, double> &, 
		 DistMat<Scalar,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim, class Scalar>
  void computeH2LS( DistGeoState &,
                       DistSVec<double,3> &, DistSVec<double,dim> &,
                       DistNodalGrad<dim, double> &, DistMat<Scalar,1> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(RecFcn *, DistSVec<double,3> &, DistVec<double> &,
                           DistMat<Scalar1,dim> &, DistSVec<double,dim> &,
                           DistSVec<double,dim> &, DistSVec<double,dim> &,
                           DistSVec<double,dim> &, DistSVec<Scalar2,dim> &,
                           DistNodalGrad<dim, Scalar2> &, DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2>
  void computeMatVecProdH2LS(RecFcn *, DistSVec<double,3> &,
             DistVec<double> &, DistMat<Scalar1,1> &,
             DistSVec<double,1> &, DistSVec<double,1> &,
             DistSVec<double,1> &, DistSVec<double,1> &,
             DistVec<Scalar2> &, 
             DistVec<Scalar2> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2T(RecFcn *, DistSVec<double,3> &,
                DistVec<double> &, DistMat<Scalar1,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2Tb(RecFcn *, DistSVec<double,3> &,
                DistVec<double> &, DistMat<Scalar1,dim> &,
                DistNodalGrad<dim, Scalar2> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<int dim, class Scalar>
  void assemble(CommPattern<Scalar> *, DistSVec<Scalar,dim> &);
                                                                                                                          
  template<class Scalar>
  void assemble(CommPattern<Scalar> *, DistVec<double> &);
                                                                                                                          
  template<class Scalar, int dim>
  bool readVectorFromFile(const char *, int, double *, DistSVec<Scalar,dim> &, Scalar* = 0);

  template<class Scalar, int dim>
  void writeVectorToFile(const char *, int, double, DistSVec<Scalar,dim> &, Scalar* = 0);

  template<int dim>
  void computeForceDerivs(VarFcn *, DistSVec<double,3> &, DistSVec<double,dim> &, 
                          DistSVec<double ,dim> &, Vec<double> &, VecSet< DistSVec<double,3> > &);
/*
  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, DistSVec<double,3> &, 
                                DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, Vec3D &);
*/

  void testNormals(DistVec<Vec3D> &, DistVec<double> &, DistVec<Vec3D> &, DistVec<double> &);

  template<int dim>
  int checkSolution(VarFcn *, DistSVec<double,dim> &);

  template<int dim>
  int checkSolution(VarFcn *, DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void checkFailSafe(VarFcn*, DistSVec<double,dim>&, DistSVec<bool,2>&);

  template<int dim, int neq>
  int clipSolution(TsData::Clipping, BcsWallData::Integration, 
		   VarFcn*, double*, DistSVec<double,dim>&);

  template<int dim>
  void checkGradients(DistSVec<double,3> &, DistVec<double> &, 
		      DistSVec<double,dim> &, DistNodalGrad<dim> &);

  template<int dim>
  void checkMatVecProd(DistSVec<double,dim> &, const char *);

  template<int dim>
  void zeroInternalVals(DistSVec<double, dim> &);

  template<int dim>
  void printVariable(DistSVec<double,dim>&, VarFcn *);
                                                                                             
  template<int dim>
  void printInletVariable(DistSVec<double,dim>&);
  template<int dim>
  void printAllVariable(DistSVec<double,3> &, DistSVec<double,dim>&, int );

  void printPhi(DistSVec<double,3> &, DistVec<double> &, int);

  template<int dim>
  void checkExtrapolationValue(DistSVec<double,dim>&, VarFcn*,
                               DistBcData<dim>&, DistGeoState&);
  template<class Scalar, int neq>
  void printAllMatrix( DistMat<Scalar,neq> &, int );

  // XML Debugging routine
  void checkNormalSums(IoData &);

  template<int dim>
  void padeReconstruction(VecSet<DistSVec<double, dim> >&, VecSet<DistSVec<double, dim> >&, int*, double*, double, int, int, int, int );

// Included (MB)
  int computeDerivativeOfControlVolumes(double, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &);

  void computeDerivativeOfNormals(DistSVec<double,3> &, DistSVec<double,3> &, DistVec<Vec3D> &, DistVec<Vec3D> &,
                         DistVec<double> &, DistVec<double> &, DistVec<Vec3D> &, DistVec<Vec3D> &, DistVec<double> &, DistVec<double> &);

  void computeDerivativeOfWeightsLeastSquares(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,6> &);

  void computeDerivativeOfWeightsGalerkin(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistSVec<double,3> &, DistSVec<double,3> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares(DistSVec<double,3> &, DistSVec<double,3> &,
                    DistSVec<double,6> &, DistSVec<double,6> &,
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsGalerkin(DistVec<double> &, DistVec<double> &,
                DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                  DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(DistVec<double> &, DistVec<double> &, DistVec<double> &, DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&,
			       DistSVec<double,3>&, DistSVec<double,3>&, DistSVec<double,dim>&, DistSVec<double,dim>&,
			       DistNodalGrad<dim>&, DistEdgeGrad<dim>*, double,
			       DistSVec<double,dim>&);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, DistSVec<double,dim> &, double, DistSVec<double,dim> &);

  template<int dim>
  void applyBCsToDerivativeOfResidual(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfSmagorinskyLESTerm(SmagorinskyLESTerm *, DistSVec<double,3> &,
				 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeOnlyGalerkinTerm(FemEquationTerm *, DistBcData<dim> &, 
			   DistGeoState &, DistSVec<double,3> &, 
			   DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim, class Scalar>
  void applyBCsToH2Jacobian(BcFcn *, DistBcData<dim> &, 
			  DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  template<int dim>
  void computeBCsJacobianWallValues(FemEquationTerm *, DistBcData<dim> &, 
				   DistGeoState &, DistSVec<double,3> &, 
				   DistSVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobianWallValues(BcFcn *, DistBcData<dim> &, 
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<int dim, class Scalar2>
  void applyBCsToProduct(BcFcn *, DistBcData<dim> &, DistSVec<double,dim> &, DistSVec<Scalar2,dim> &);

  template<int dim>
  void computeDerivativeOfVolumicForceTerm(VolumicForceTerm *, DistVec<double> &, DistVec<double> &,
                               DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>  
  void computeDerivativeOfInvReynolds(FemEquationTerm *, VarFcn *, DistGeoState &, 
			     DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, 
			     DistVec<double> &, DistVec<double> &, DistVec<double> &, DistVec<double> &, DistVec<double> &, double, double, double, double, double, double);

  template<int dim>
  void fixSolution(VarFcn *, DistSVec<double,dim> &, DistSVec<double,dim> &);

 };

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Domain.C>
#endif

#endif
