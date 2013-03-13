#ifndef _NONLINEAR_ROM_H_
#define _NONLINEAR_ROM_H_

#include <RestrictionMapping.h>
#include <vector>


template <int dim>
class NonlinearRom {

// This class handles all input and output for the nonlinear ROM database.
// Currently, the ROM databse is implemented as a simple file system.  In
// the future, if it is necessary to implement a higher performance
// database, it would only be necessary to change the methods in this class.
// Note that every other nonlinear ROM class is derived from this base class.

  protected:
  Domain& domain;
  Communicator* com; 
  IoData* ioData;

  // IO directory information
  const char* databasePrefix;
  const char* databaseName;
  const char* clusterName;
  const char* sensitivityClusterName;

  // State snapshot clusters
  char* stateSnapsName;
  char* mapName;  
  char* indexName;
  char* connName;
  char* centersName;
  char* nearestName;

  // State bases
  char* stateBasisPrefix;
  char* stateBasisName;
  char* stateSingValsName;
  char* updateInfoName;
  char* stateFastDistCalcInfoName;
  char* projErrorName;
  char* refStateName;

  // Krylov snaps
  char* krylovSnapsName;

  // Krylov bases
  char* krylovBasisPrefix;
  char* krylovBasisName;
  char* krylovSingValsName;
  char* krylovFastDistCalcInfoName;

  // Sensitivities
  char* sensitivitySnapsName;

  // Sensitivity Basis
  char* sensitivityBasisPrefix;
  char* sensitivityBasisName;
  char* sensitivitySingValsName;

  // Residual snaps
  char* residualSnapsName;

  // Residual Bases
  char* residualBasisPrefix;
  char* residualBasisName;
  char* residualSingValsName;

  // Action-of-Jacobian snaps
  char* jacActionSnapsName;

  // Action-of-Jacobian Bases
  char* jacActionBasisPrefix;
  char* jacActionBasisName;
  char* jacActionSingValsName;

  // GNAT quantities
  char* sampledNodesName;
  char* sampledNodesFullCoordsName;
  char* sampledStateBasisName;
  char* sampledResidualBasisName;
  char* sampledJacActionBasisName;
  char* sampledMeshName;
  char* sampledSolutionName;
  char* sampledWallDistName;
  char* gappyJacActionName;
  char* gappyResidualName;

  // Surface quantities
  char* surfaceStateBasisName;
  char* surfaceSolutionName;
  char* surfaceWallDistName;
  char* surfaceMeshName;

  // ROM database data
  VecSet< DistSVec<double, dim> >* snap; // snap(nTotSnaps, domain.getNodeDistInfo())
  VecSet< DistSVec<double, dim> >* clusterCenters; // average of all snapshots in a cluster 
  VecSet< DistSVec<double, dim> >* nearestSnapsToCenters; // closest snapshot to each cluster center 
  int* snapsInCluster; // number of snaps in each cluster
  int* clusterIndex; // stores original cluster association for each snapshot (before any overlap is introduced)
  int** clusterSnapshotMap;  // one vector per cluster, lists snapshots to include in the cluster (including overlapping snapshots) 
  int** clusterNeighbors;  // one vector per cluster, lists neighboring clusters
  int* clusterNeighborsCount; // stores number of neighbors for each cluster
  int* clusterNewtonCount;  // counts number of residuals/jacActions stored in each cluster
  int* clusterKrylovCount;  // counts number of krylov vectors stored in each cluster
  DistSVec<double, dim>* snapRefState; 

  // thin svd update quantities
  double* columnSumsV;
  double* sVals;
  DistSVec<double, dim>* Uinit; 

  // private database IO functions
  void createDirectories();
  void outputClusterBasis(int, int, char*);
  void outputClusterSnapshots(char*);
  void outputClusterReferenceState(int, DistSVec<double, dim> &);  // automatically stores reference state used to form snapshots for each cluster.
  void readClusterReferenceState(int);  // read the reference state that was automatically stored for each cluster.
  void readClusterSnapshots(int, bool, char*, int first = 0, int last = 0);
  void readNearestSnapsToCenters();
  void readReferenceState();  // read a reference state specified by the user

  // for local GNAT preprocessing

  // for local GNAT online simulations
  int nSampleNodes;
  std::vector<int> sampleNodes;
  int numResJacMat;
  VecSet<DistSVec<double, dim> >* resMat;
  VecSet<DistSVec<double, dim> >* jacMat;
  std::auto_ptr< RestrictionMapping<dim> > restrictionMapping;
  void readClusterSampleNodes(int);
  void readClusterGappyMatrix(int, char*); 

  // ASCII output files 
  FILE* clustUsageFile;
  FILE* reducedCoordsFile;

  int nState, nKrylov, nSens;

  public:

  NonlinearRom(Communicator *, IoData &, Domain &);
  ~NonlinearRom();

  int nClusters;
  VecSet< DistSVec<double, dim> >* basis;

  // basic functions
  double distance(DistSVec<double, dim> &, DistSVec<double, dim> &);
  void closestCenterFull(DistSVec<double, dim> &, int* index1 = NULL, int* index2 = NULL, double* dist1 = NULL, double* dist2 = NULL);

  // public database IO functions
	void determineFileName(const char*, const char*, const char*, char*&);
  void determinePath(char*, int, char*&); // passing -1 as the cluster number indicates that the file resides in the top-level database directory
  void readClusterBasis(int, char*);
  void readClusterCenters();
  void writeClusteredBinaryVectors(int, DistSVec<double,dim> *, DistSVec<double,dim> *, DistSVec<double,dim> *);
  void initializeClusteredOutputs(); 

  // for online ROMs (both with and without hyper-reduction)
  virtual void readDistanceCalcInfo() {};
  virtual void closestCenter(DistSVec<double, dim> &, int*) {}; 
  virtual void updateBasis(int, DistSVec<double, dim> &) {};
  virtual void appendNonStateDataToBasis(int, char*) {};
  virtual void readClusterOnlineQuantities(int) {};
  void writeReducedCoords(const int, bool, bool, int, Vec<double>); 

  // for online ROMs with hyper-reduction
  int getNumResJacMat() {return numResJacMat;}
  VecSet<DistSVec<double,dim> >* getResMat() {return resMat;}
  VecSet<DistSVec<double,dim> >* getJacMat() {if (numResJacMat==2) { return jacMat; } else { return resMat;} }
  const DistInfo& getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();}
  RestrictionMapping<dim>* restrictMapping() { return restrictionMapping.get(); } 
  // no const in last function prototype

  //virtual int getNumResJacMat() { return 0; }
  //virtual VecSet<DistSVec<double,dim> >* getResMat() { return NULL; }
  //virtual VecSet<DistSVec<double,dim> >* getJacMat() { return NULL; }
  //const virtual DistInfo& getRestrictedDistInfo () { return domain.getNodeDistInfo(); }
  //virtual RestrictionMapping<dim>* restrictMapping() { return NULL; }

 // virtual void appendVectorToBasis(DistSVec<double, dim>*, int numVec = 0) {};

};

#include "NonlinearRom.C"
#endif
