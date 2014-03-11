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
  Timer *timer;

  // IO directory information
  const char* databasePrefix;
  const char* databaseName;
  const char* clusterName;
  const char* sensitivityClusterName;

  NonlinearRomFilesData* romFiles;

  // State snapshot clusters
  char* stateSnapsName;
  char* mapName;  
  char* indexName;
  char* connName;
  char* centersName;
  char* nearestName;
  char* centerNormsName;
  char* distanceMatrixName;

  // State bases
  char* stateBasisPrefix;
  char* stateBasisName;
  char* stateSingValsName;
  char* simpleUpdateInfoName;
  char* exactUpdateInfoPrefix;      // only user-specified value for exact update file names
  char* basisBasisProductsName;     // for exact updates (exactUpdateInfoPrefix.exactUpdates_F)
  char* basisUrefProductsName;      // for exact updates (exactUpdateInfoPrefix.exactUpdates_e)
  char* basisUicProductsName;       // for exact updates (exactUpdateInfoPrefix.exactUpdates_d)
  char* urefUicProductsName;        // for exact updates (exactUpdateInfoPrefix.exactUpdates_c)
  char* urefUrefProductsName;       // for exact updates (exactUpdateInfoPrefix.exactUpdates_g)   
  char* urefComponentwiseSumsName;  // for exact updates with uniform IC (exactUpdateInfoPrefix.exactUpdates_UrefComponentwiseSums)
  char* basisComponentwiseSumsName; // for exact updates with uniform IC (exactUpdateInfoPrefix.exactUpdates_StateBasisComponentwiseSums)
  //char* approxUpdateInfoName; // treating this as a GNAT quantity
  char* stateDistanceComparisonInfoName;
  char* centerComponentWiseSumsName; // for fast distance calcs when using either no updates or exact updates ()
  char* stateDistanceComparisonInfoExactUpdatesName;
  char* projErrorName;
  char* refStateName;

  // Krylov snaps
  char* krylovSnapsName;

  // Krylov bases
  char* krylovBasisPrefix;
  char* krylovBasisName;
  char* krylovSingValsName;
  char* krylovDistanceComparisonInfoName;

  // Sensitivities
  char* sensitivitySnapsName;

  // Sensitivity Basis
  char* sensitivityBasisPrefix;
  char* sensitivityBasisName;
  char* sensitivitySingValsName;
  char* sensitivityDistanceComparisonInfoName;

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
  char* sampledCentersName;
  char* sampledStateBasisName;
  char* sampledKrylovBasisName;
  char* sampledSensitivityBasisName;
  char* sampledResidualBasisName;
  char* sampledJacActionBasisName;
  char* sampledMeshName;
  char* sampledSolutionName;
  char* sampledRefStateName;
  char* sampledWallDistName;
  char* gappyJacActionName;
  char* gappyResidualName;
  char* approxMetricLowRankName;
  char* approxMetricLowRankFullCoordsName;
  // Surface quantities
  //char* surfaceStateBasisName;
  //char* surfaceSolutionName;
  //char* surfaceWallDistName;
  //char* surfaceMeshName;

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
  // 1: common to all update methods (simple updates, exact updates, approx updates)
  double rTol;
  std::vector<double>* columnSumsV;
  std::vector<double>* sVals;
  DistSVec<double, dim>* Uref; 
  // 2: unique to exact updates
  double uicNorm;
  std::vector<std::vector<std::vector<std::vector<double> > > > basisBasisProducts;  // [iCluster][pCluster][:][:]
  std::vector<std::vector<std::vector<double> > > basisUrefProducts;  // [Cluster_Basis][Cluster_Uref][:]
  std::vector<std::vector<double> > basisUicProducts;  // [iCluster][1:nPod] only precomputed if Uic specified
  std::vector<double> urefUicProducts; // [iCluster] only precomputed if Uic specified
  std::vector<std::vector<double> > urefUrefProducts; //[iCluster][jCluster] symmetric (lower triangular)
  std::vector<std::vector<double> > urefComponentwiseSums; //[iCluster][1:dim]
  std::vector<std::vector<std::vector<double> > > basisComponentwiseSums;  // [iCluster][iVec][1:dim]
  std::vector<double> exactUpdatesAlpha;  // [jVec]
  std::vector<std::vector<double> > exactUpdatesBeta; //[iCluster][jVec]  
  std::vector<std::vector<std::vector<double> > > exactUpdatesN;  // [iCluster][iVec][jVec]
  double exactUpdatesAlphaSwitch;  // scalar
  std::vector<double> exactUpdatesBetaSwitch; //[iCluster]
  std::vector<std::vector<double> > exactUpdatesNSwitch;  // [iCluster][iVec]
  DistSVec<double, dim>* Uic;
  // 3: unique to approximate updates
  VecSet<DistSVec<double, dim> >* lowRankFactor; // low rank factor for approx metric

  // fast distance calculation quantities
  // 1: common to all GNAT update methods (no updates, exact updates, approx updates)
  bool specifiedIC;
  SVec<double, dim>* uniformIC;  // value of uniform initial condition at node 0 (should be representative)
  std::vector<std::vector<double> > centerNorms;
  std::vector<std::vector<std::vector<std::vector<double> > > > stateBasisCentersProduct;  //[iCluster][mCenter][pCenter][:]
  std::vector<std::vector<std::vector<std::vector<double> > > > krylovBasisCentersProduct; //[iCluster][mCenter][pCenter][:]
  std::vector<std::vector<std::vector<double> > > sensitivityBasisCentersProduct;          //[mCenter][pCenter][:]
  std::vector<std::vector<double> > distanceComparisons;  // this is "z_(m,p)" from Amsallem et al., INJME 2012, but with p<m
  void checkUniformInitialCondition(DistSVec<double, dim> &);
  void checkForSpecifiedInitialCondition();
  // 2: unique to exact updates
  std::vector<std::vector<double> > initialConditionCentersProduct; 
  std::vector<std::vector<std::vector<double> > > refStateCentersProduct;
  // 3: unique to approximate updates
  double ***hForFastDistComp;
  double ***cForFastDistComp;
  
  // non-database IO function
  int readSnapshotFiles(char *, bool);
  std::vector<int> stateSnapsFromFile;   // stateSnapsFromFile[iFile] = number of snapshots taken from file iFile
  std::vector<std::vector<double> > stateSnapshotTags; // stateSnapshotInfo[iFile][iSnap] = tag associated with snapshot
                                                       // iSnap from file iFile

  // database IO functions
  void createDirectories();
  void outputClusteredSnapshots(char*);
  void readClusteredSnapshots(int, bool, char*, int first = 0, int last = 0);
  void outputClusteredBasis(int, int, char*);  // readClusteredBasis is public
  void outputClusteredReferenceState(int, DistSVec<double, dim> &);  // automatically stores snapshot reference state
  void readClusteredReferenceState(int, char*);  // read the reference state that was automatically stored for each cluster.
  void readNearestSnapsToCenters();
  void readReferenceState();  // read a reference state specified by the user
  void readClusteredSampleNodes(int iCluster, bool deleteExistingRestrictionMapping = true);
  void readClusteredGappyMatrix(int, char*);

  void outputCenterNorms(std::vector<std::vector<double> > &);
  void readCenterNorms();
  void outputClusteredInfoASCII(int, char*, std::vector<double>* vec1 = NULL, 
                                std::vector<std::vector<double> >* vec2 = NULL,
                                std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                                std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void readClusteredInfoASCII(int, char*, std::vector<double>* vec1 = NULL,
                              std::vector<std::vector<double> >* vec2 = NULL,
                              std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                              std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void writeMultiVecASCII(char*, std::vector<double>* vec1 = NULL, 
                          std::vector<std::vector<double> >* vec2 = NULL,
                          std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                          std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void readMultiVecASCII(char*, std::vector<double>* vec1 = NULL,
                         std::vector<std::vector<double> >* vec2 = NULL,
                         std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                         std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);


  // for local GNAT preprocessing

  // for local GNAT online simulations
  int nSampleNodes;
  std::vector<int> sampleNodes;
  int numResJacMat;
  VecSet<DistSVec<double, dim> >* resMat;
  VecSet<DistSVec<double, dim> >* jacMat;
  RestrictionMapping<dim>* restrictionMapping;
  // for storing all online quantities in memory (optional)
  bool storedAllOnlineQuantities;
  std::vector<int>** allSampleNodes;
  VecSet<DistSVec<double, dim> >** allResMat;
  VecSet<DistSVec<double, dim> >** allJacMat;
  VecSet<DistSVec<double, dim> >** allStateBases;
  VecSet<DistSVec<double, dim> >** allKrylovBases;
  VecSet<DistSVec<double, dim> >* sensitivityBasis;
  std::vector<double>** allStateSVals; 
  std::vector<double>** allKrylovSVals;
  std::vector<double>* sensitivitySVals;
  VecSet<DistSVec<double, dim> >* allRefStates;
  std::vector<double>** allColumnSumsV; 
  RestrictionMapping<dim>** allRestrictionMappings;

  // ASCII output files 
  FILE* clustUsageFile;
  FILE* reducedCoordsFile;

  int nState, nKrylov, nSens;

  public:

  NonlinearRom(Communicator *, IoData &, Domain &);
  ~NonlinearRom();

  int nClusters;
  int nFullMeshNodes;
  int nLowRankFactors;
  VecSet< DistSVec<double, dim> >* basis;


  // online selection of closest cluster center (calls either closestCenterFull or closestCenterFast)
  void closestCenter(DistSVec<double, dim> &, int* index1=NULL);

  // calculate closest center to current state using full vectors
  void closestCenterFull(DistSVec<double, dim> &, int* index1=NULL, int* index2=NULL, double* dist1=NULL, double* dist2=NULL);
  void distancesToCentersFull(DistSVec<double, dim> &, std::vector<double> &, int* closest=NULL);
  double distanceFull(DistSVec<double, dim> &, DistSVec<double, dim> &);

  // calculate closest center to current state without using full vectors (approach depends on ROB update method)
  void closestCenterFast(int* index1=NULL);
  void initializeDistanceComparisons(DistSVec<double, dim> &);
  void resetDistanceComparisonQuantitiesApproxUpdates();
  void incrementDistanceComparisons(Vec<double> &, int);  // calls one of the following three functions
  void incrementDistanceComparisonsForNoUpdates(Vec<double> &, int);
  void incrementDistanceComparisonsForExactUpdates(Vec<double> &, int);
  void incrementDistanceComparisonsForApproxUpdates(Vec<double> &, int);

  // public database IO functions
	void determineFileName(const char*, const char*, const char*, char*&);
	void determinePrefixName(const char*, const char*, char*&);
  void determinePath(char*, int, char*&); // top-level database directory is cluster "-1", sensitivity basis is cluster "-2"
  void readClusteredBasis(int, char*, bool relProjError = false);
  void readClusteredColumnSumsV(int, char*);
  void readClusteredUpdateInfo(int, char*);
  void readNonClusteredUpdateInfo(char*);
  void readExactUpdateInfo();
  void readClusterCenters(char*);
  void readAllClusteredOnlineQuantities();
  void readApproxMetricLowRankFactor(char *);
  void readDistanceComparisonInfo(char*); 
  void writeClusteredBinaryVectors(int, DistSVec<double,dim> *, DistSVec<double,dim> *, DistSVec<double,dim> *);
  void initializeClusteredOutputs(); 

  // for online ROMs (both with and without hyper-reduction)
  virtual bool updateBasis(int, DistSVec<double, dim> &, Vec<double>* coords = NULL) {return false;};
  virtual void appendNonStateDataToBasis(int, char*, bool relProjError = false) {};
  virtual void readClusteredOnlineQuantities(int) {};
  void writeReducedCoords(const int, bool, bool, int, Vec<double>); 
  void initializeFastExactUpdatesQuantities(DistSVec<double, dim> &);

  // for online ROMs with hyper-reduction
  void determineNumResJacMat(); 
  void deleteRestrictedQuantities();
  int getNumResJacMat() {return numResJacMat;}
  VecSet<DistSVec<double,dim> >* getResMat() {return resMat;}
  VecSet<DistSVec<double,dim> >* getJacMat() {if (numResJacMat==2) { return jacMat; } else { return resMat;} }
  const DistInfo& getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();}
  RestrictionMapping<dim>* restrictMapping() { return restrictionMapping; } 

  virtual void appendVectorToBasis(DistSVec<double, dim>&, int numVec = 0) {};

};

#include "NonlinearRom.C"
#endif
