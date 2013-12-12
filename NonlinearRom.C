#include <TsInput.h>
#include <cmath>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

using std::stable_sort;

extern "C"      {
   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}


template<int dim> 
NonlinearRom<dim>::NonlinearRom(Communicator *_com, IoData &_ioData, Domain &_domain)  : 
com(_com), ioData(&_ioData), domain(_domain)
{ 
  //NOTE: ioData->example, com->example, domain.example

  nClusters = ioData->romDatabase.nClusters;  // overwritten later if there are actually fewer clusters
  nFullMeshNodes = 0;  // read from centerNorms file if necessary

  // Directory information
  databasePrefix = ioData->romDatabase.directories.prefix;
  databaseName = ioData->romDatabase.directories.databaseName;
  clusterName = ioData->romDatabase.directories.clusterName;
  sensitivityClusterName  = ioData->romDatabase.directories.sensitivityClusterName;

  romFiles = &(ioData->romDatabase.files); 

  // State snapshot clusters
  determineFileName(romFiles->stateSnapsName, "snaps", romFiles->statePrefix, stateSnapsName);
  determineFileName(romFiles->mapName, "map", romFiles->statePrefix, mapName);
  determineFileName(romFiles->indexName, "index", romFiles->statePrefix, indexName);
  determineFileName(romFiles->connName, "conn", romFiles->statePrefix, connName);
  determineFileName(romFiles->centersName, "centers", romFiles->statePrefix, centersName);
  determineFileName(romFiles->nearestName, "nearest", romFiles->statePrefix, nearestName);
  determineFileName(romFiles->centerNormsName, "centerNorms", romFiles->statePrefix, centerNormsName);
  determineFileName(romFiles->distanceMatrixName, "distanceMatrix", romFiles->statePrefix, distanceMatrixName);

  // State bases
  determinePrefixName(romFiles->stateBasisPrefix, romFiles->statePrefix, stateBasisPrefix);
  determineFileName(romFiles->stateBasisName, "rob", stateBasisPrefix, stateBasisName);
  determineFileName(romFiles->stateSingValsName, "svals", stateBasisPrefix, stateSingValsName);
  determineFileName(romFiles->projErrorName, "proj", stateBasisPrefix, projErrorName);
  determineFileName(romFiles->refStateName, "refState", stateBasisPrefix, refStateName); 

  // Update info for state bases (this is a bit tricky)
  determineFileName(romFiles->simpleUpdateInfoName, "simpleUpdates", stateBasisPrefix, simpleUpdateInfoName);
  determineFileName(romFiles->exactUpdateInfoName, "exactUpdates", stateBasisPrefix, exactUpdateInfoName);
  determineFileName(romFiles->approxUpdateInfoName, "approxUpdates", stateBasisPrefix, approxUpdateInfoName);
  // TODO: define file names for each individual exact update precomputed quantity (_a, _b, _c, _d, _e, _F, _g)
  determineFileName(romFiles->stateDistanceComparisonInfoName, "dist", stateBasisPrefix, stateDistanceComparisonInfoName);
  determineFileName(romFiles->stateDistanceComparisonInfoExactUpdatesName, "distExactUpdates", stateBasisPrefix, stateDistanceComparisonInfoExactUpdatesName);
 


  // Krylov snaps
  determineFileName(romFiles->krylovSnapsName, "snaps", romFiles->krylovPrefix, krylovSnapsName);

  // Krylov bases
  determinePrefixName(romFiles->krylovBasisPrefix, romFiles->krylovPrefix, krylovBasisPrefix);
  determineFileName(romFiles->krylovBasisName, "rob", krylovBasisPrefix, krylovBasisName);
  determineFileName(romFiles->krylovSingValsName, "svals", krylovBasisPrefix, krylovSingValsName);
  determineFileName(romFiles->krylovDistanceComparisonInfoName, "dist", krylovBasisPrefix, krylovDistanceComparisonInfoName);

  // Sensitivity snaps
  determineFileName(romFiles->sensitivitySnapsName, "snaps", romFiles->sensitivityPrefix, sensitivitySnapsName);

  // Sensitivity basis
  determinePrefixName(romFiles->sensitivityBasisPrefix, romFiles->sensitivityPrefix, sensitivityBasisPrefix);
  determineFileName(romFiles->sensitivityBasisName, "rob", sensitivityBasisPrefix, sensitivityBasisName);
  determineFileName(romFiles->sensitivitySingValsName, "svals", sensitivityBasisPrefix, sensitivitySingValsName);
  determineFileName(romFiles->sensitivityDistanceComparisonInfoName, "dist", sensitivityBasisPrefix, sensitivityDistanceComparisonInfoName);

  // Residual snaps
  determineFileName(romFiles->residualSnapsName, "snaps", romFiles->residualPrefix, residualSnapsName);

  // Residual bases
  determinePrefixName(romFiles->residualBasisPrefix, romFiles->residualPrefix, residualBasisPrefix);
  determineFileName(romFiles->residualBasisName, "rob", residualBasisPrefix, residualBasisName);
  determineFileName(romFiles->residualSingValsName, "svals", residualBasisPrefix, residualSingValsName);

  // Action-of-Jacobian snaps
  determineFileName(romFiles->jacActionSnapsName, "snaps", romFiles->jacActionPrefix, jacActionSnapsName);

  // Action-of-Jacobian bases
  determinePrefixName(romFiles->jacActionBasisPrefix, romFiles->jacActionPrefix, jacActionBasisPrefix);
  determineFileName(romFiles->jacActionBasisName, "rob", jacActionBasisPrefix, jacActionBasisName);
  determineFileName(romFiles->jacActionSingValsName, "svals", jacActionBasisPrefix, jacActionSingValsName);

  // GNAT quantities
  determineFileName(romFiles->sampledNodesName, "sampledNodes", romFiles->gnatPrefix, sampledNodesName);
  determineFileName(romFiles->sampledNodesFullCoordsName, "sampledNodesFullCoords", romFiles->gnatPrefix, sampledNodesFullCoordsName);
  determineFileName(romFiles->sampledStateBasisName, "sampledStateROB", romFiles->gnatPrefix, sampledStateBasisName);
  determineFileName(romFiles->sampledSensitivityBasisName, "sampledSensitivityROB", romFiles->gnatPrefix, sampledSensitivityBasisName);
  determineFileName(romFiles->sampledKrylovBasisName, "sampledKrylovROB", romFiles->gnatPrefix, sampledKrylovBasisName);
  determineFileName(romFiles->sampledResidualBasisName, "sampledResROB", romFiles->gnatPrefix, sampledResidualBasisName);
  determineFileName(romFiles->sampledJacActionBasisName, "sampledJacROB", romFiles->gnatPrefix, sampledJacActionBasisName);
  determineFileName(romFiles->sampledMeshName, "top", romFiles->gnatPrefix, sampledMeshName);
  determineFileName(romFiles->sampledSolutionName, "sol", romFiles->gnatPrefix, sampledSolutionName);
  determineFileName(romFiles->sampledRefStateName, "sampledRefState", romFiles->gnatPrefix, sampledRefStateName);
  determineFileName(romFiles->sampledWallDistName, "dwall", romFiles->gnatPrefix, sampledWallDistName);
  determineFileName(romFiles->gappyJacActionName, "gappyJac", romFiles->gnatPrefix, gappyJacActionName);
  determineFileName(romFiles->gappyResidualName, "gappyRes", romFiles->gnatPrefix, gappyResidualName);
  determineFileName(romFiles->approxMetricLowRankName, "approxMetric", romFiles->gnatPrefix, approxMetricLowRankName);
  // Surface quantities
  //determineFileName(romFiles->gappyResidualName, "sampledStateROB", romFiles->surfacePrefix, surfaceStateBasisName);
  //determineFileName(romFiles->gappyResidualName, "sol", romFiles->surfacePrefix, surfaceSolutionName);
  //determineFileName(romFiles->gappyResidualName, "dwall", romFiles->surfacePrefix, surfaceWallDistName);
  //determineFileName(romFiles->gappyResidualName, "top", romFiles->surfacePrefix, surfaceMeshName);


  basis = NULL;
  snap = NULL; // snap(nTotSnaps, domain.getNodeDistInfo())
  clusterCenters = NULL; // average of all snapshots in a cluster 
  nearestSnapsToCenters = NULL; // closest snapshot to each cluster center 
  snapsInCluster = NULL; // number of snaps in each cluster
  clusterIndex = NULL; // stores original cluster association for each snapshot (before any overlap is introduced)
  clusterSnapshotMap = NULL;  // one vector per cluster, lists snapshots to include in the cluster (including overlapping snapshots) 
  clusterNeighbors = NULL;  // one vector per cluster, lists neighboring clusters
  clusterNeighborsCount = NULL; // stores number of neighbors for each cluster
  snapRefState = NULL;
  columnSumsV = NULL;
  sVals = NULL;
  Uref = NULL;
  clusterNewtonCount = NULL;
  clusterKrylovCount = NULL;
 
  nSampleNodes = 0;
  sampleNodes.clear();
  numResJacMat = 0;
  resMat = NULL;
  jacMat = NULL;
  restrictionMapping = NULL;

  clustUsageFile = NULL;
  reducedCoordsFile = NULL;

  // initialize ASCII outputs
  if (!(strcmp(ioData->output.rom.clusterUsage,"")==0) && nClusters>0) {
    char *fullClustUsageName = new char[strlen(ioData->output.rom.prefix) + 1 + strlen(ioData->output.rom.clusterUsage) + 1];
    sprintf(fullClustUsageName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.clusterUsage);
    if (this->com->cpuNum() == 0)  clustUsageFile = fopen(fullClustUsageName, "wt");
    delete [] fullClustUsageName;
  }

  if (strcmp(ioData->output.rom.reducedCoords,"")==0) {
    if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT)
      this->com->fprintf(stderr, "\n*** Warning: Reduced coordinates output file not specified\n\n");
  } else {
    char *fullReducedCoordsName = new char[strlen(ioData->output.rom.prefix) + 1 + strlen(ioData->output.rom.reducedCoords) + 1];
    sprintf(fullReducedCoordsName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.reducedCoords);
    if (this->com->cpuNum() == 0)  reducedCoordsFile = fopen(fullReducedCoordsName, "wt");
    delete [] fullReducedCoordsName;
  }

  nState = 0;
  nKrylov = 0;
  nSens = 0;

  storedAllOnlineQuantities = false;
  allSampleNodes = NULL;
  allResMat = NULL;
  allStateBases = NULL;
  allKrylovBases = NULL;
  sensitivityBasis = NULL;
  allStateSVals = NULL;
  allKrylovSVals = NULL;
  sensitivitySVals = NULL;
  allRefStates = NULL;
  allColumnSumsV = NULL; 
  allRestrictionMappings = NULL;

  // for fast distance calculations
  specifiedIC = true;
  uniformIC = NULL;

}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRom<dim>::~NonlinearRom() 
{
  delete [] stateSnapsName;
  delete [] mapName;  
  delete [] indexName;
  delete [] connName;
  delete [] centersName;
  delete [] nearestName;
  delete [] stateBasisName;
  delete [] stateSingValsName;
  delete [] simpleUpdateInfoName;
  delete [] exactUpdateInfoName;
  delete [] approxUpdateInfoName;
  delete [] stateDistanceComparisonInfoName;
  delete [] stateDistanceComparisonInfoExactUpdatesName;
  delete [] refStateName;
  delete [] projErrorName;
  delete [] krylovSnapsName;
  delete [] krylovBasisName;
  delete [] krylovSingValsName;
  delete [] krylovDistanceComparisonInfoName;
  delete [] residualSnapsName;
  delete [] jacActionSnapsName;
  delete [] residualBasisName;
  delete [] residualSingValsName;
  delete [] jacActionBasisName;
  delete [] jacActionSingValsName;
  delete [] sampledNodesName;
  delete [] sampledStateBasisName;
  delete [] sampledKrylovBasisName;
  delete [] sampledSensitivityBasisName;
  delete [] sampledResidualBasisName;
  delete [] sampledJacActionBasisName;
  delete [] sampledMeshName;
  delete [] sampledSolutionName;
  delete [] sampledRefStateName;
  delete [] sampledWallDistName;
  delete [] gappyJacActionName;
  delete [] gappyResidualName;
  delete [] approxMetricLowRankName;
  //delete [] surfaceStateBasisName;
  //delete [] surfaceSolutionName;
  //delete [] surfaceWallDistName;
  //delete [] surfaceMeshName;

  delete [] stateBasisPrefix;
  delete [] krylovBasisPrefix;
  delete [] sensitivityBasisPrefix;
  delete [] residualBasisPrefix;
  delete [] jacActionBasisPrefix;

  if (basis) delete basis;
  if (snap) delete snap; 
  if (clusterCenters) delete clusterCenters;
  if (nearestSnapsToCenters) delete nearestSnapsToCenters;
  if (snapsInCluster) delete [] snapsInCluster;
  if (clusterIndex) delete [] clusterIndex;
  if (clusterNeighborsCount) delete [] clusterNeighborsCount;
  if (clusterNewtonCount) delete [] clusterNewtonCount;
  if (clusterKrylovCount) delete [] clusterKrylovCount;
  if (snapRefState) delete snapRefState;
  if (columnSumsV) delete columnSumsV; 
  if (sVals) delete sVals;
  if (Uref) delete Uref;
  
  //TODO
  //clusterSnapshotMap
  //clusterNeighbors

  if (resMat) delete resMat;
  if (jacMat) delete jacMat;

  if (this->com->cpuNum() == 0) {
    if (clustUsageFile) fclose(clustUsageFile);
    if (reducedCoordsFile) fclose(reducedCoordsFile);
  }

  if (storedAllOnlineQuantities) {
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      if (allSampleNodes) delete allSampleNodes[iCluster];
      if (allResMat) delete allResMat[iCluster];
      if (allJacMat) delete allJacMat[iCluster];
      if (allStateBases) delete allStateBases[iCluster];
      if (allStateSVals) delete allStateSVals[iCluster];
      if (allKrylovBases) delete allKrylovBases[iCluster];
      if (allKrylovSVals) delete allKrylovSVals[iCluster];
      if (allColumnSumsV) delete allColumnSumsV[iCluster];
      if (allRestrictionMappings) delete allRestrictionMappings[iCluster];
    }
    if (allSampleNodes) delete [] allSampleNodes;
    if (allResMat) delete [] allResMat;
    if (allJacMat) delete [] allJacMat;
    if (allStateBases) delete [] allStateBases;
    if (allStateSVals) delete [] allStateSVals;
    if (allKrylovBases) delete [] allKrylovBases;
    if (allKrylovSVals) delete [] allKrylovSVals;
    if (allColumnSumsV) delete [] allColumnSumsV;
    if (allRestrictionMappings) delete [] allRestrictionMappings; 
    if (sensitivityBasis) delete sensitivityBasis;
    if (sensitivitySVals) delete sensitivitySVals;
    if (allRefStates) delete allRefStates;
  } else {
    if (restrictionMapping) delete restrictionMapping;
  }

  if (uniformIC) delete uniformIC;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenter(DistSVec<double, dim> &vec, int *index1) {

  // for determining the closest cluster center during an online ROM simulation

  if (ioData->romOnline.distanceComparisons) {
    closestCenterFast(index1);
  } else {
    closestCenterFull(vec, index1);
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenterFull(DistSVec<double, dim> &vec, int *index1, int *index2, double *dist1, double *dist2) {
// Computes the closest cluster center to vec using full scale vectors.

  double tmp;
  
  int i1 = 0; // index of closest cluster
  int i2 = 0; // index of second closest cluster
  double d1 = 0;  // distance to closest cluster center
  double d2 = 0;  // distance to second closest cluster center

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    tmp = distanceFull( vec, (*clusterCenters)[iCluster]);
    if ((tmp<d1) || (iCluster==0)) {
      d2 = d1;
      i2 = i1;
      d1 = tmp;
      i1 = iCluster;
    } else if ((tmp<d2) || (iCluster==1)) {
      d2 = tmp;
      i2 = iCluster;
    }
  }
 
  if (index1) *index1 = i1;
  if (index2) *index2 = i2;
  if (dist1) *dist1 = d1;
  if (dist2) *dist2 = d2;

}

//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRom<dim>::distanceFull(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {
// distance calculation using full-size state vectors

  DistSVec<double, dim> diff(domain.getNodeDistInfo());
  diff = (U1 - U2);
  double dist = diff.norm();

  return dist;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenterFast(int *index1) {
// returns the index of the closest cluster center (using precomputed distance quantities)

  int closestCluster = 0;

  for (int iCluster=1; iCluster<nClusters; ++iCluster) {
    closestCluster = (distanceComparisons[iCluster][closestCluster] < 0) ? iCluster : closestCluster;
  }

  *index1 = closestCluster;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::checkUniformInitialCondition(DistSVec<double, dim> &ic) {

    if (uniformIC) return;
    
    // double check that the IC is indeed uniform
    double minIC[dim], maxIC[dim];
    ic.min(minIC);
    ic.max(maxIC);

    for (int iDim=0; iDim<dim; ++iDim) {
      if (minIC[iDim]!=maxIC[iDim]) { // should be exact, no tolerance needed
       com->fprintf(stderr, " *** Error: a uniform initial condition was assumed during preprocessing");
       exit(-1);
      }
    }
 
    uniformIC = new SVec<double, dim>( ic(0) ); // value of initial condition at node 0 (should be representative)

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeDistanceComparisons(DistSVec<double, dim> &ic) {

  std::vector<double> centerMinusICNorms;
  centerMinusICNorms.resize(nClusters);

  if (specifiedIC) { // preprocessing was performed for a specified IC 

    // centerNorms is nClusters-by-1 with format:  || Center_0 - IC ||^2;  || Center_1 - IC ||^2; ...
    
    for (int mCenter=0; mCenter<nClusters; ++mCenter) {
      centerMinusICNorms[mCenter] = centerNorms[mCenter][0]; 
    }

  } else { // preprocesing was performed assuming a uniform initial condition for the online ROM

    // centerNorms is nClusters-by-(dim+1) with format:
    //    ||Center_0||^2, sum( Center_0(iDim=0) ), ... , sum( Center_0(iDim=dim-1) ); ||Center_1||^2, ...

    checkUniformInitialCondition(ic); // checks if IC is uniform, sets variable uniformIC

    for (int mCenter=0; mCenter<nClusters; ++mCenter) {
      centerMinusICNorms[mCenter] = centerNorms[mCenter][0];
      for (int iDim=0; iDim<dim; ++iDim) { 
        centerMinusICNorms[mCenter] -= 2.0 * (*uniformIC)[0][iDim] * centerNorms[mCenter][iDim+1];
        centerMinusICNorms[mCenter] += ((double) nFullMeshNodes) * (*uniformIC)[0][iDim] * (*uniformIC)[0][iDim];
      }
    }

    if (ioData->romOnline.basisUpdates == NonlinearRomOnlineData::UPDATES_FAST_EXACT) {
      initialConditionCentersProduct.resize(nClusters);
      for (int mCenter=1; mCenter<nClusters; ++mCenter) {
        initialConditionCentersProduct[mCenter].resize(mCenter);
        for (int pCenter=0; pCenter<mCenter; ++pCenter) { 
          initialConditionCentersProduct[mCenter][pCenter] = 0.0;
          for (int iDim=0; iDim<dim; ++iDim) {
            initialConditionCentersProduct[mCenter][pCenter] += 2.0 * (*uniformIC)[0][iDim] * 
                                                                (centerNorms[pCenter][iDim+1] - centerNorms[mCenter][iDim+1]);
          }
        }
      }
    }
  }

  // the diagonal elements are all zero, but element [0][0] is inserted to make the indexing more intuitive ([mCenter][pCenter])
  distanceComparisons.clear();
  distanceComparisons.resize(nClusters);  

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    distanceComparisons[mCenter].reserve(mCenter);
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      distanceComparisons[mCenter].push_back(centerMinusICNorms[mCenter] - centerMinusICNorms[pCenter]);
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisons(Vec<double> &dUTimeIt, int currentCluster) {
 
    switch (ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        incrementDistanceComparisonsForNoUpdates(dUTimeIt, currentCluster);
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        com->fprintf(stderr, "*** Error: fast distance comparisons are incompatible with simple ROB updates (use Exact)\n");
        exit(-1);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        incrementDistanceComparisonsForExactUpdates(dUTimeIt, currentCluster);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        incrementDistanceComparisonsForApproxUpdates(dUTimeIt, currentCluster);
        break;
      default:
        this->com->fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
    }
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForNoUpdates(Vec<double> &dUTimeIt, int currentCluster) {
 
  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      for (int iState=0; iState<nState; ++iState) {
        distanceComparisons[mCenter][pCenter] += stateBasisCentersProduct[currentCluster][mCenter][pCenter][iState] 
                                                 * dUTimeIt[iState];
      }
      for (int iKrylov=0; iKrylov<nKrylov; ++iKrylov) {
        // account for Krylov bases
      }
      for (int iSens=0; iSens<nSens; ++iSens) {
        // account for sensitivity bases
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForExactUpdates(Vec<double> &dUTimeIt, int currentCluster) {

  return;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForApproxUpdates(Vec<double> &dUTimeIt, int currentCluster) {

  return;
}

//----------------------------------------------------------------------------------

template<int dim>
int NonlinearRom<dim>::readSnapshotFiles(char* snapType, bool preprocess) {

  // Check for snapshot command file
  char *vecFile;
  bool typeIsState = false;

  if (strcmp(snapType, "state")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.stateSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.stateSnapFile);
    typeIsState = true;
  } else if (strcmp(snapType,"sensitivity")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.sensitivitySnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.sensitivitySnapFile);
  } else if (strcmp(snapType,"approxMetricState")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.approxMetricStateSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.approxMetricStateSnapFile);
//  } else if (strcmp(snapType,"residual")==0) {
//    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.residualSnapFile) + 1];
//    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.residualSnapFile);
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", snapType);
    exit (-1);
  }

  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    this->com->fprintf(stderr, "*** Error: No snapshots FILES in %s\n", vecFile);
    exit (-1);
  }

  delete [] vecFile;
  vecFile = NULL;

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  this->com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

  if (typeIsState) stateSnapsFromFile.resize(nData, 0);

  char** snapFile = new char*[nData];
  for (int iData = 0; iData < nData; ++iData)
    snapFile[iData] = new char[500];
  char snapFile1[500];
  int* numSnaps = new int[nData];
  int* startSnaps = new int[nData];
  int* endSnaps = new int[nData];
  int* sampleFreq = new int[nData];
  double* snapWeight = new double[nData];
  int nSnap, iStart, iEnd, iFreq;
  double weight;

  if (typeIsState && ioData->romOffline.rob.state.snapshots.incrementalSnaps)
    this->com->fprintf(stderr, "*** Warning: Incremental snapshots is not supported for multiple bases (yet) \n");

  if (typeIsState && (ioData->romOffline.rob.state.dataCompression.energyOnly == DataCompressionData::ENERGY_ONLY_TRUE))
    this->com->fprintf(stderr, "*** Warning: EnergyOnly is not supported for multiple bases\n");

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %lf", snapFile1,&iStart,&iEnd,&iFreq,&weight);
    if (iStart < 1) iStart = 1;
    if (iEnd < 0) iEnd = 0;
    if (iFreq < 1) iFreq = 1;
    //numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart - 1;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    this->com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);
  }

  // compute the total number of snapshots
  int nTotSnaps = 0;
  int dummyStep = 0;
  double dummyTag = 0.0;
  for (int iData = 0; iData < nData; ++iData) {
    bool status = this->domain.template readTagFromFile<double, dim>(snapFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));

    if (!status) {
      this->com->fprintf(stderr, "*** Error: could not read snapshots from %s \n", snapFile[iData]);
      exit(-1);
    }

    if ((endSnaps[iData]==0) || (endSnaps[iData]>numSnaps[iData]))
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
        if (typeIsState) ++stateSnapsFromFile[iData]; 
      }
    }
  }

  bool incrementalSnaps = false;
  bool subtractRefSol = false;
  if (preprocess) {
    if (ioData->romOffline.rob.relativeProjectionError.projectIncrementalSnaps) {
      incrementalSnaps = true;
      nTotSnaps -= nData;
    } else if (ioData->romOffline.rob.relativeProjectionError.subtractRefSol) {
      subtractRefSol = true;
      if (!(ioData->input.stateSnapRefSolution)) {
        this->com->fprintf(stderr, "*** Error: Reference solution not found \n");
        exit (-1);
      }
      this->readReferenceState();
    }
  }

  if (snap) delete snap;
  snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufOld = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufNew = new DistSVec<double, dim>(this->domain.getNodeDistInfo());

  if (typeIsState) {
    stateSnapshotTags.resize(nData);
    for (int iFile=0;iFile<nData;++iFile)
      stateSnapshotTags[iFile].resize(stateSnapsFromFile[iFile], -1.0);
  }

  *snapBufOld = 0.0;
  *snapBufNew = 0.0;

  double* tags = new double [nTotSnaps];
  double tagOld;
  double tagNew;

  int numCurrentSnapshots = 0;
  bool status;

  if (subtractRefSol) {
    *snapBufOld = *(this->snapRefState);
    delete (this->snapRefState);
    (this->snapRefState) = NULL;
  }

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        if ((iSnap == startSnaps[iData]) && incrementalSnaps) {
          status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tagOld, *snapBufOld);
        } else {
          status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tagNew, *snapBufNew);
          tags[numCurrentSnapshots] = tagNew;
          (*snap)[numCurrentSnapshots] = *snapBufNew - *snapBufOld;  //snapBufOld = 0 if not using incremental snaps
          if (incrementalSnaps) *snapBufOld = *snapBufNew;
          if (snapWeight[iData]) (*snap)[numCurrentSnapshots] *= snapWeight[iData]; //CBM--check
          ++numCurrentSnapshots;
        }
      }
    }
  }

  if (typeIsState) {
    int snapCount=0;
    for (int iFile=0;iFile<nData;++iFile) {
      for (int iSnap=0;iSnap<stateSnapsFromFile[iFile];++iSnap) {
        stateSnapshotTags[iFile][iSnap] = tags[snapCount];
        ++snapCount;
      }
    }
  }

  delete snapBufOld;
  snapBufOld = NULL;
  delete snapBufNew;
  snapBufNew = NULL;

  for (int iData=0; iData < nData; ++iData) {
    delete [] snapFile[iData];
  }
  delete [] snapFile;
  snapFile = NULL;
  delete [] numSnaps;
  numSnaps = NULL;
  delete [] startSnaps;
  startSnaps = NULL;
  delete [] endSnaps;
  endSnaps = NULL;
  delete [] sampleFreq;
  sampleFreq = NULL;
  delete [] snapWeight;
  snapWeight = NULL;
  delete [] tags;
  tags = NULL;

  return nData;

}
//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredSnapshots(char* snapType)  { 
 
  int nTotSnaps = snap->numVectors();

  if (strcmp(snapType, "state")==0) { 

    createDirectories();
  
    if (com->cpuNum() == 0) {
  
      // output cluster index as ASCII file (original clusters before overlap)
      FILE *clustIndexFile;
      char *clustIndexPath = 0; 
      determinePath(indexName, -1, clustIndexPath); 
      com->fprintf(stdout, "\nWriting cluster index to disk\n");
  
      clustIndexFile = fopen(clustIndexPath, "wt");
      if (clustIndexFile) {
         com->fprintf(clustIndexFile,"Snapshot# OriginalCluster\n");
      } else {  
         com->fprintf(stdout,"***Error: Cannot open %s\n",clustIndexPath);
         exit(-1);
      }

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        com->fprintf(clustIndexFile,"%d %d\n", iSnap, clusterIndex[iSnap]);
      }
      fclose (clustIndexFile);
      delete [] clustIndexPath;
      clustIndexPath = NULL;
 
      // output cluster snapshot map as ASCII file (clusters after overlap)
      if (clusterSnapshotMap) {
        FILE *clustMapFile;
        char *clustMapPath = 0;
        determinePath(mapName, -1, clustMapPath);
        com->fprintf(stdout, "\nWriting cluster-snapshot map to disk\n");
  
        clustMapFile = fopen(clustMapPath, "wt");
        com->fprintf(clustMapFile,"Snapshot# Cluster\n");
  
        for (int iCluster=0; iCluster<nClusters; ++iCluster) {
          for (int iSnap=0; iSnap<snapsInCluster[iCluster]; ++iSnap) {
            com->fprintf(clustMapFile,"%d %d\n", clusterSnapshotMap[iCluster][iSnap], iCluster);
          }
        }
        fclose (clustMapFile);
        delete [] clustMapPath;
        clustMapPath = NULL;
      }

      // output cluster connectivity as ASCII file
      if (clusterNeighbors) {
        FILE *clustConnFile;
        char *clustConnPath = 0;
        determinePath(connName, -1, clustConnPath);
        com->fprintf(stdout, "\nWriting cluster connectivity to disk\n");
  
        clustConnFile = fopen(clustConnPath, "wt");
        com->fprintf(clustConnFile,"Cluster#  ...Neighbors...\n");
        for (int iCluster=0; iCluster<nClusters; ++iCluster) {
          com->fprintf(clustConnFile,"%d", iCluster);
          for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {
            com->fprintf(clustConnFile," %d", clusterNeighbors[iCluster][iNeighbor]);
          }
          com->fprintf(clustConnFile,"\n");
        }
        fclose (clustConnFile);
        delete [] clustConnPath;
        clustConnPath = NULL;
      }
    }

    com->barrier();
  
    // output cluster centers
    char *clustCentersPath = 0;
    determinePath(centersName, -1, clustCentersPath);
    com->fprintf(stdout, "\nWriting cluster centers to disk\n");
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->barrier();
      domain.writeVectorToFile(clustCentersPath, iCluster, double(snapsInCluster[iCluster]), (*clusterCenters)[iCluster]);
    }
    delete [] clustCentersPath;
    clustCentersPath = NULL;  

    // output nearest snap to each cluster
    char *nearestSnapsPath = 0;
    determinePath(nearestName, -1, nearestSnapsPath);
    com->fprintf(stdout, "\nWriting nearest snapshot to each center to disk\n");
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->barrier();
      domain.writeVectorToFile(nearestSnapsPath, iCluster, double(snapsInCluster[iCluster]), (*nearestSnapsToCenters)[iCluster]);
    }
    delete [] nearestSnapsPath;
    nearestSnapsPath = NULL;  

    // output clustered snapshots
    for (int iCluster=0; iCluster<nClusters; iCluster++) {
  
      char *snapshotsPath = 0;
      determinePath(stateSnapsName, iCluster, snapshotsPath);
      com->fprintf(stdout, "\nWriting %d snapshots to cluster %d\n", snapsInCluster[iCluster], iCluster);
  
      int numWritten = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        for (int jSnap=0; jSnap<snapsInCluster[iCluster]; ++jSnap) {
          if (iSnap == clusterSnapshotMap[iCluster][jSnap]) {
            domain.writeVectorToFile(snapshotsPath, numWritten, double(numWritten), (*snap)[iSnap] );
            ++numWritten;
          }
        }
      } 
      delete [] snapshotsPath;
      snapshotsPath = NULL;
    }
  
    com->fprintf(stdout, "\nFreeing memory for parallel SVD; read in snapshots as needed\n");
  
    delete snap;
    snap = NULL;
    delete clusterCenters;
    clusterCenters = NULL;
    delete nearestSnapsToCenters;
    nearestSnapsToCenters = NULL;
    delete [] clusterIndex;
    clusterIndex = NULL;
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) delete [] clusterSnapshotMap[iCluster];
    delete [] clusterSnapshotMap;
    clusterSnapshotMap = NULL;
    delete [] snapsInCluster;
    snapsInCluster = NULL;
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) delete [] clusterNeighbors[iCluster];
    delete [] clusterNeighbors;
    clusterNeighbors = NULL;
    delete [] clusterNeighborsCount;
    clusterNeighborsCount = NULL;

  } else if (strcmp(snapType, "sensitivity")==0) {

    // output sensitivities

    char *sensitivityPath = 0;
    determinePath(sensitivitySnapsName, -2, sensitivityPath);
    com->fprintf(stdout, "\nWriting %d snapshots to sensitivity cluster\n", nTotSnaps);

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      domain.writeVectorToFile(sensitivityPath, iSnap, double(iSnap), (*snap)[iSnap] );
    } 
    delete [] sensitivityPath;
    sensitivityPath = NULL;  

    com->fprintf(stdout, "\nFreeing memory for parallel SVD; read in snapshots as needed\n");

    delete snap;
    snap = NULL;
  }
}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredSnapshots(int iCluster, bool preprocess, char *basisType, int first, int last) {

  int nTotSnaps;
  int normalizeSnaps;

  // reconstruct cluster name
  char *snapshotsPath = 0;
  if (strcmp(basisType,"state")==0) {
    if (strcmp(stateSnapsName,"")==0) return;
    determinePath(stateSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : snapsInCluster[iCluster];
    com->fprintf(stdout, "\nReading %d snapshots from cluster %d \n", nTotSnaps, iCluster);
    normalizeSnaps = ioData->romOffline.rob.state.snapshots.normalizeSnaps;
    if ((ioData->romOffline.rob.state.snapshots.subtractRefState) && 
        !(ioData->romOffline.rob.state.snapshots.subtractCenters || ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters)){
      if (!(ioData->input.stateSnapRefSolution)) {
        com->fprintf(stderr, "*** Error: reference solution file not specified\n");
        exit (-1);
      }
      readReferenceState();
    }
  } else if (strcmp(basisType,"residual")==0) {
    if (strcmp(this->residualSnapsName,"")==0) return;
    determinePath(residualSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading residual snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.residual.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"jacAction")==0) {
    if (strcmp(this->jacActionSnapsName,"")==0) return;
    determinePath(jacActionSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading action-of-Jacobian snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.jacAction.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"krylov")==0) {
    if (strcmp(this->krylovSnapsName,"")==0) return;
    determinePath(krylovSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading krylov snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.krylov.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"sensitivity")==0) {
    if (strcmp(this->sensitivitySnapsName,"")==0) return;
    determinePath(sensitivitySnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading sensitivity from snapshots\n");
    normalizeSnaps = ioData->romOffline.rob.sensitivity.snapshots.normalizeSnaps;
  } else {
    exit (-1);
  }

  if (snap) delete snap;
  snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());

  // read in Snapshot Vectors
  double tmp;
  bool status = true; 
  int snapCount = 0;
  int snapIndex = first; 
  DistSVec<double, dim>* tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {
    status = domain.readVectorFromFile(snapshotsPath, snapIndex, &tmp, *tmpSnap);
    if (!status) {
      if (snapCount<nTotSnaps) {
        nTotSnaps = snapCount;
        VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (*snapNew)[iSnap] = (*snap)[iSnap];
        }
        delete snap;
        snap = snapNew;
        snapNew = NULL;
      }
      break;
    }
    if (snapCount==nTotSnaps) {
      if (last>=0) break;
      ++nTotSnaps;
      VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
      for (int iSnap=0; iSnap<(nTotSnaps-1); ++iSnap) {
        (*snapNew)[iSnap] = (*snap)[iSnap];
      }
      delete snap;
      snap = snapNew;
      snapNew = NULL;
    }
    (*snap)[snapCount] = *tmpSnap;
    ++snapCount;
    ++snapIndex;
  }

  delete tmpSnap;
  tmpSnap = NULL;
  delete [] snapshotsPath;
  snapshotsPath = NULL;

  if (preprocess) {
    if (strcmp(basisType,"state")==0) { 
      if (ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters) {
        com->fprintf(stdout, " ... subtracting nearest snapshot to center from snapshots \n");
        if (ioData->romOffline.rob.state.snapshots.subtractCenters)
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring subtractCenters command \n");
        if (ioData->romOffline.rob.state.snapshots.subtractRefState) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
        if (ioData->romOffline.rob.state.snapshots.incrementalSnaps) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
        --nTotSnaps; //
        VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
        tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());
        outputClusteredReferenceState(iCluster, (*nearestSnapsToCenters)[iCluster]);
        com->barrier();
        snapCount = 0;
        for (int iSnap=0; iSnap<=nTotSnaps; ++iSnap) {
          *tmpSnap = (*snap)[iSnap] - (*nearestSnapsToCenters)[iCluster];
          if ((*tmpSnap).norm() > 1e-6) { // one will be exactly zero; don't include this one
            (*snapNew)[snapCount] = *tmpSnap;
            ++snapCount;
          }
        }
        delete tmpSnap;
        tmpSnap = NULL;
        delete snap;
        snap = snapNew;
        snapNew = NULL;
      } else if (ioData->romOffline.rob.state.snapshots.subtractCenters) {
        com->fprintf(stdout, " ... subtracting cluster center from snapshots \n");
        if (ioData->romOffline.rob.state.snapshots.subtractRefState) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
        if (ioData->romOffline.rob.state.snapshots.incrementalSnaps) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
        outputClusteredReferenceState(iCluster,(*clusterCenters)[iCluster]);
        com->barrier();
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (*snap)[iSnap] -= (*clusterCenters)[iCluster];
        }
      } else if (ioData->romOffline.rob.state.snapshots.subtractRefState) {
        com->fprintf(stdout, " ... subtracting reference solution from snapshots \n");
        if (ioData->romOffline.rob.state.snapshots.incrementalSnaps) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (*snap)[iSnap] -= *snapRefState;
        }
        outputClusteredReferenceState(iCluster,(*snapRefState));
        com->barrier();
        delete snapRefState;
        snapRefState = NULL;
      } else {
        if (ioData->romOffline.rob.state.snapshots.incrementalSnaps)
          com->fprintf(stderr, "*** Warning: incremental snapshots not supported for local ROM construction\n");
        com->fprintf(stdout, " ... using raw states as snapshots\n");
        DistSVec<double, dim> zeroVec(domain.getNodeDistInfo(), 0);
        outputClusteredReferenceState(iCluster, zeroVec);
        com->barrier();
      }
    }

    if (normalizeSnaps) {
      com->fprintf(stdout, " ... normalizing snapshots \n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        double tmpNorm = ((*snap)[iSnap]).norm();
        if (tmpNorm != 0) {
          double normalize = 1.0 / tmpNorm;
          (*snap)[iSnap] *= normalize;
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredReferenceState(int iCluster, DistSVec<double, dim>& ref) {
  // Automatically stores the reference snapshot for cluster iCluster.
  // This functionality is in place to reduce the user's workload.

  char *refStatePath = 0;
  determinePath(refStateName, iCluster, refStatePath);

  com->fprintf(stdout, " ... storing reference state to disk for later use\n", iCluster);
  domain.writeVectorToFile(refStatePath, 0, double(iCluster), ref);

  delete [] refStatePath;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredReferenceState(int iCluster, char* refType) {
  // This function reads in the automatically stored reference snapshot for a cluster.
  // By storing these reference snapshots and reading them automatically it reduces the user's workload.

  if (Uref) {
    delete Uref;
    Uref = NULL;
  }

  if (storedAllOnlineQuantities) {
    com->fprintf(stdout, " ... loading snapshot reference state for cluster %d\n", iCluster);
    Uref = new DistSVec<double, dim>(domain.getNodeDistInfo());
    *Uref = (*allRefStates)[iCluster];
    return;
  }

  char *refStatePath = 0;
  if (strcmp(refType,"state")==0) {
      determinePath(refStateName, iCluster, refStatePath);
  } else if (strcmp(refType,"sampledState")==0) {
      determinePath(sampledRefStateName, iCluster, refStatePath);
  } else {
      exit (-1);
  }

  double tmp;
  bool status;
 
  com->fprintf(stdout, "Reading reference snapshot %s\n", refStatePath);
  Uref = new DistSVec<double, dim>(domain.getNodeDistInfo());
  status = domain.readVectorFromFile(refStatePath, 0, &tmp, *Uref);

  delete [] refStatePath;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredCenters() {

  // never called by GNAT so there's no need to support reading sampled cluster centers

  char *clustCentersPath = 0;
  determinePath(centersName, -1, clustCentersPath);

  if (nClusters <= 0) {
    com->fprintf(stderr, "\n*** Error: invalid value for NumClusters (%d)\n", nClusters);
    exit(-1);
  }

  // read centers

  if (clusterCenters) delete clusterCenters;  
  clusterCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  if (snapsInCluster) delete [] snapsInCluster;
  snapsInCluster = new int[nClusters];

  com->fprintf(stdout, "\nReading cluster centers\n");
  bool status;
  int expectedClusters = nClusters;
  nClusters = 0;
  double tmp;
  DistSVec<double, dim>* tmpVec = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {

    status = domain.readVectorFromFile(clustCentersPath, nClusters, &tmp, *tmpVec);
    if (!status) break;
 
    ++nClusters;

    if (nClusters > expectedClusters) {
      com->fprintf(stderr, "\n*** Error: found more clusters than expected (NumClusters was specified as %d)\n", nClusters);
      exit(-1);
    }   

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
    int* snapsInClusterNew = new int[nClusters];

    for (int iCluster=0; iCluster<(nClusters-1); ++iCluster) {
      (*clusterCentersNew)[iCluster] = (*clusterCenters)[iCluster];
      snapsInClusterNew[iCluster] = snapsInCluster[iCluster];
    }

    (*clusterCentersNew)[nClusters-1] = *tmpVec;
    snapsInClusterNew[nClusters-1] = int(tmp);

    delete clusterCenters;
    delete [] snapsInCluster;

    clusterCenters = clusterCentersNew;
    snapsInCluster = snapsInClusterNew;

    clusterCentersNew = NULL;
    snapsInClusterNew = NULL;
  }

  com->fprintf(stdout, "\n ... found %d clusters in this database\n", nClusters);
  com->barrier();
  delete tmpVec;
  tmpVec = NULL;
  delete [] clustCentersPath;
  clustCentersPath = NULL;

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readNearestSnapsToCenters() {
  // NOTE: this function must only be run after readClusteredCenters (this is because the
  // cluster centers file defines snapsInCluster and nClusters)

  char *nearestSnapsPath = 0;
  determinePath(nearestName, -1, nearestSnapsPath);

  if (nearestSnapsToCenters) delete nearestSnapsToCenters; 
  nearestSnapsToCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  bool status;

  com->fprintf(stdout, "\nReading closest snapshot to each center\n");

  // read in Snapshot Vectors
  double tmp;
  for (int iCluster = 0; iCluster<nClusters; ++iCluster) {
    status = domain.readVectorFromFile(nearestSnapsPath, iCluster, &tmp, (*nearestSnapsToCenters)[iCluster]);
    if (!status) {
       com->fprintf(stderr,"*** Error: readNearestSnapsToCenters attempted to read %d vecs from a file with only %d.\n", nClusters, iCluster-1);
       exit(-1);
    }
  }

  delete [] nearestSnapsPath;
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredBasis(int iCluster, char* basisType, bool relProjError) {

  if (storedAllOnlineQuantities) {
    if (basis) delete basis;
    if (sVals) delete sVals;
    if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      com->fprintf(stdout, " ... loading state ROB for cluster %d\n", iCluster);
      nState = allStateBases[iCluster]->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nState, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nState; ++iVec)
        (*basis)[iVec] = (*(allStateBases[iCluster]))[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *(allStateSVals[iCluster]);
    } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0) ) {
      com->fprintf(stdout, " ... loading sensitivity ROB\n");
      nSens = sensitivityBasis->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nSens, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nSens; ++iVec)
        (*basis)[iVec] = (*sensitivityBasis)[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *sensitivitySVals;
    } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      com->fprintf(stdout, " ... loading Krylov ROB for cluster %d\n", iCluster);
      nKrylov = allKrylovBases[iCluster]->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nKrylov, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nKrylov; ++iVec)
        (*basis)[iVec] = (*(allKrylovBases[iCluster]))[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *(allKrylovSVals[iCluster]);
    }
    return;
  }

  int maxDimension;
  int minDimension;
  double energyTol;

  if (ioData->problem.type[ProblemData::NLROMOFFLINE]) {
    if (relProjError) {
      if (strcmp(basisType,"state")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"residual")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"jacAction")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"sensitivity")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.sensitivity.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.sensitivity.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.sensitivity.energy;
      } else if (strcmp(basisType,"krylov")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.krylov.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.krylov.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.krylov.energy;
      }
    } else { 
      if (strcmp(basisType,"state")==0) {
        maxDimension = ioData->romOffline.gnat.maxDimensionState;
        minDimension = ioData->romOffline.gnat.minDimensionState;
        energyTol = ioData->romOffline.gnat.energyState;
      } else if (strcmp(basisType,"sensitivity")==0) {
        maxDimension = ioData->romOffline.gnat.maxDimensionSensitivity;
        minDimension = ioData->romOffline.gnat.minDimensionSensitivity;
        energyTol = ioData->romOffline.gnat.energySensitivity;
      } else if (strcmp(basisType,"krylov")==0) {
        maxDimension = ioData->romOffline.gnat.maxDimensionKrylov;
        minDimension = ioData->romOffline.gnat.minDimensionKrylov;
        energyTol = ioData->romOffline.gnat.energyKrylov;
      } else if (strcmp(basisType,"residual")==0) {
        maxDimension = ioData->romOffline.gnat.maxDimensionResidual;
        minDimension = ioData->romOffline.gnat.minDimensionResidual;
        energyTol = ioData->romOffline.gnat.energyResidual;
      } else if (strcmp(basisType,"jacAction")==0) {
        if (ioData->romOffline.gnat.maxDimensionJacAction <= 0) {
          com->fprintf(stdout, "*** Warning: JacAction greedy parameters not specified; using Residual parameters\n");
          maxDimension = ioData->romOffline.gnat.maxDimensionResidual;
          minDimension = ioData->romOffline.gnat.minDimensionResidual;
          energyTol = ioData->romOffline.gnat.energyResidual;
        } else {
          maxDimension = ioData->romOffline.gnat.maxDimensionJacAction;
          minDimension = ioData->romOffline.gnat.minDimensionJacAction;
          energyTol = ioData->romOffline.gnat.energyJacAction;
        }
      }
    }
  } else { // ONLINE
    if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      maxDimension = ioData->romOnline.maxDimension;
      minDimension = ioData->romOnline.minDimension;
      energyTol = ioData->romOnline.energy;
    } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0) ) {
      maxDimension = ioData->romOnline.sensitivity.maxDimension;
      minDimension = ioData->romOnline.sensitivity.minDimension;
      energyTol = ioData->romOnline.sensitivity.energy;
    } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      maxDimension = ioData->romOnline.krylov.maxDimension;
      minDimension = ioData->romOnline.krylov.minDimension;
      energyTol = ioData->romOnline.krylov.energy;
    }
  }

  char* singValsPath = 0;
  char* basisPath = 0;
  
  if (strcmp(basisType,"state")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      determinePath(stateBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sampledState")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      determinePath(sampledStateBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"residual")==0) { 
      determinePath(residualSingValsName, iCluster, singValsPath);
      determinePath(residualBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"jacAction")==0) {  
      determinePath(jacActionSingValsName, iCluster, singValsPath);
      determinePath(jacActionBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"krylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(krylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sampledKrylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(sampledKrylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sensitivity")==0) {
      iCluster = -2;
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sensitivityBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sampledSensitivity")==0) {
      iCluster = -2;
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sampledSensitivityBasisName, iCluster, basisPath);
  } else {
      exit (-1);
  }

  int step = 0;
  double tag = 0.0;
  int numSteps = 0;
  int status = domain.readTagFromFile<double, dim>(basisPath, step, &tag, &numSteps);

  if (!status)  {
    com->fprintf(stderr, "*** Error: unable to open file %s\n", basisPath);
    exit (-1);
  }

  if (maxDimension<=0) maxDimension = numSteps;

  FILE *singValFile = fopen(singValsPath, "r");

  if (!singValFile)  {
    com->fprintf(stderr, "*** Error: unable to open file %s\n", singValsPath);
    exit (-1);
  }

  delete [] singValsPath;
  singValsPath = NULL;

  int _n;
  int vecNumber;
  double percentEnergy;
  double tmpSVal;

  if (sVals) delete sVals;
  sVals = new std::vector<double>;

  for (int iVec=0; iVec<maxDimension; ++iVec) {
    _n = fscanf(singValFile,"%d %le %le", &vecNumber, &tmpSVal, &percentEnergy);
    if (_n == 3) {
      sVals->push_back(tmpSVal);
      if ((percentEnergy>=energyTol)&&((iVec+1)>=minDimension)) break;
    } else if (feof(singValFile)) {
      break;
    } else {
      com->fprintf(stderr, "*** Error: fscanf interrupted by non-EOF error\n");
      exit(-1);
    }
  }

  int basisSize = sVals->size();
  fclose(singValFile);

  // read in basis vectors

  com->fprintf(stdout, "\nReading basis %d\n", iCluster);

  if (basis) delete basis;
  basis = new VecSet< DistSVec<double, dim> >(basisSize, domain.getNodeDistInfo());

  int numBasisVecs = 0;

  for (int iVec = 0; iVec<basisSize; ++iVec) {
    status = domain.readVectorFromFile(basisPath, iVec, &tmpSVal, (*basis)[iVec]);
    if (!status) break;
    ++numBasisVecs; 
  }

  delete [] basisPath;
  basisPath = NULL;

  if (numBasisVecs != basisSize) {
    VecSet< DistSVec<double, dim> >* basisNew =  new VecSet< DistSVec<double, dim> >(numBasisVecs, domain.getNodeDistInfo());

    for (int iVec=0; iVec<numBasisVecs; ++iVec) {
      (*basisNew)[iVec]=(*basis)[iVec];
    }

    delete basis;
    basis = basisNew;
    basisNew = NULL;    
  }    

  // for ASCII output files (clusterUsage and reducedCoords)
  if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      nState = numBasisVecs;
  } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      nKrylov = numBasisVecs;
  } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0)) {
      nSens = numBasisVecs;
  }

  com->fprintf(stdout, "\n");

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredUpdateInfo(int iCluster, char* basisType) {

  // always need Uref and ColumnSumsV regardless of update algorithm
  readClusteredReferenceState(iCluster, basisType);
  readClusteredColumnSumsV(iCluster, basisType);

  if (false) readClusteredExactUpdateInfo(iCluster, basisType);
  if (false) readClusteredApproxUpdateInfo(iCluster, basisType);

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredExactUpdateInfo(int iCluster, char* basisType) {



}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredApproxUpdateInfo(int iCluster, char* basisType) {


}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredColumnSumsV(int iCluster, char* basisType) {

  if (columnSumsV) {
    delete columnSumsV;
    columnSumsV = NULL;
  }

  if ((strcmp(basisType, "state") == 0) || (strcmp(basisType, "sampledState") == 0) ) { 

    if (storedAllOnlineQuantities) {
      com->fprintf(stdout, " ... loading columnSumsV for cluster %d\n", iCluster);
      columnSumsV = new std::vector<double>;
      *columnSumsV = *(allColumnSumsV[iCluster]);
      return;
    }

    char *basisUpdatePath = NULL;
    determinePath(simpleUpdateInfoName, iCluster, basisUpdatePath);

    com->fprintf(stdout, "\nReading update info for basis %d \n", iCluster);

    FILE *basisUpdateFile = fopen(basisUpdatePath, "r");

    if (!basisUpdateFile)  {
      com->fprintf(stderr, "*** Error: unable to open file %s\n", basisUpdatePath);
      exit (-1);
    }

    columnSumsV = new std::vector<double>;
    double tmp;

    while (true) {
      int _n = fscanf(basisUpdateFile, "%le", &tmp);
      if (_n == 1) {
        columnSumsV->push_back(tmp);
      } else if (feof(basisUpdateFile)) {
        break;
      } else {
        com->fprintf(stderr, "*** Error: fscanf interrupted by non-EOF error\n");
        exit(-1);
      }
    }

    fclose(basisUpdateFile);
    delete [] basisUpdatePath;
    basisUpdatePath = NULL;

  } else {
    com->fprintf(stderr, "*** Error: unexpected ROB type (%s)\n", basisType);
    exit(-1);
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredBasis(int iCluster, int nTotSnaps, char* basisType) {

  int podSize = basis->numVectors();

  char* singValsPath = 0;
  char* basisPath = 0;

  if (strcmp(basisType,"state")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      determinePath(stateBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"residual")==0) {
      determinePath(residualSingValsName, iCluster, singValsPath);
      determinePath(residualBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"jacAction")==0) {
      determinePath(jacActionSingValsName, iCluster, singValsPath);
      determinePath(jacActionBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"krylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(krylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sensitivity")==0) {
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sensitivityBasisName, iCluster, basisPath);
  } else {
      exit (-1);
  }

  // output vectors
  com->fprintf(stdout, "\nWriting basis %d to disk\n", iCluster);
  for (int iVec=0; iVec<podSize; ++iVec) {
      domain.writeVectorToFile(basisPath, iVec, double(iVec), (*basis)[iVec] );
  }

  delete basis;
  basis = NULL;
  delete [] basisPath;
  basisPath = NULL;

  // write singular values to ASCII file
  com->fprintf(stdout, "\nWriting singular values to disk\n");

  double sValSqrdSum = 0;
  for (int iSnap = 0; iSnap<sVals->size(); ++iSnap) 
    sValSqrdSum += ((*sVals)[iSnap] * (*sVals)[iSnap]);
  double invSValSqrdSum = 1/sValSqrdSum;

  double sValSqrdPartialSum = 0;

  if (com->cpuNum() == 0) {
    FILE *singValFile = fopen(singValsPath, "wt");
    //com->fprintf(singValFile,"Vector# SVal SValSum\n");

    for (int iVec=0; iVec<sVals->size(); ++iVec) {
      sValSqrdPartialSum += ((*sVals)[iVec]*(*sVals)[iVec]*invSValSqrdSum);
      com->fprintf(singValFile,"%d %1.12e %1.12e\n", iVec, (*sVals)[iVec], sValSqrdPartialSum);
     }

    fclose (singValFile);
  }

 
  delete [] singValsPath;
  singValsPath = NULL;
  delete sVals;
  sVals = NULL;

  if (strcmp(basisType, "state") == 0) { 
    // write thinSVD update info to ASCII file
    com->fprintf(stdout, "\nWriting sums of right singular vectors to disk (for basis updates)\n");

    char *basisUpdatePath = 0;
    determinePath(simpleUpdateInfoName, iCluster, basisUpdatePath);

    if (com->cpuNum() == 0) {
      FILE *basisUpdateFile = fopen(basisUpdatePath, "wt");

      for (int iVec=0; iVec<(columnSumsV->size()); ++iVec) {
        com->fprintf(basisUpdateFile,"%le\n", (*columnSumsV)[iVec]);
      }

      fclose (basisUpdateFile);
    }

    delete [] basisUpdatePath;
    basisUpdatePath = NULL;
    delete columnSumsV;
    columnSumsV = NULL;
  }


}

//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determineFileName(const char* fileNameInput, const char* fileNameExtension, const char* prefix, char*& fileName) { 

  if (strcmp(fileNameInput,"") == 0) {
    if (strcmp(prefix,"") == 0) {
      fileName = new char[1];
      fileName[0] = 0;
    } else {
      fileName = new char[strlen(prefix) + 1 + strlen(fileNameExtension) + 1];
      sprintf(fileName, "%s.%s", prefix, fileNameExtension);
    }
  } 
  else {
    fileName = new char[strlen(fileNameInput) + 1];
    sprintf(fileName, "%s", fileNameInput); 
  }  
} 

//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determinePrefixName(const char* prefixInput, const char* prefixDefault, char*& prefix) { 

  if (strcmp(prefixInput,"") == 0) {
    if (strcmp(prefixDefault,"") == 0) {
      prefix = new char[1];
      prefix[0] = 0;
    } else {
      prefix = new char[strlen(prefixDefault) + 1];
      sprintf(prefix, "%s", prefixDefault);
    }
  } 
  else {
    prefix = new char[strlen(prefixInput) + 1];
    sprintf(prefix, "%s", prefixInput); 
  }  
} 


//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determinePath(char* fileName, int iCluster, char*& path) { 

  if (iCluster == -1) { // top level of ROM database
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s", databasePrefix, databaseName, fileName);
  } else if (iCluster == -2) { // sensitivity cluster
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(sensitivityClusterName) + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s/%s", databasePrefix, databaseName, sensitivityClusterName, fileName);
  } else { // path to appropriate cluster directory
    int addedDigits = 1;
    if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(clusterName) + addedDigits + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s%d/%s", databasePrefix, databaseName, clusterName, iCluster, fileName);
  }  
} 
//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::createDirectories() { 

  if (com->cpuNum() == 0) {
    char *fullDatabaseName = new char[strlen(databasePrefix) + 1 + strlen(databaseName)];
    sprintf(fullDatabaseName, "%s%s", databasePrefix, databaseName);
    int status;
    status = mkdir(fullDatabaseName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (int iCluster=0; iCluster<nClusters; iCluster++) {
      int addedDigits = 1;
      if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
      char *fullClusterName = new char[strlen(fullDatabaseName) + 1 + strlen(clusterName) + addedDigits + 1];
      sprintf(fullClusterName, "%s/%s%d", fullDatabaseName, clusterName, iCluster);
      int status;
      status = mkdir(fullClusterName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      delete [] fullClusterName;
      fullClusterName = NULL;
    }

    if (strcmp(sensitivityClusterName,"")!=0) {
      char *fullSensitivitiesClusterName = new char[strlen(fullDatabaseName) + 1 + strlen(sensitivityClusterName)  + 1];
      sprintf(fullSensitivitiesClusterName, "%s/%s", fullDatabaseName, sensitivityClusterName);
      int status;
      status = mkdir(fullSensitivitiesClusterName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      delete [] fullSensitivitiesClusterName;
      fullSensitivitiesClusterName = NULL;
    }

    delete [] fullDatabaseName;
    fullDatabaseName = NULL;
  }
  com->barrier();
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readReferenceState() {

  double tmp;
  bool status;
  char* fullRefName;
  const char* refSolName = ioData->input.stateSnapRefSolution;
  fullRefName = new char[strlen(ioData->input.prefix) + strlen(refSolName) + 1];
  sprintf(fullRefName, "%s%s", ioData->input.prefix, refSolName);
  com->fprintf(stdout, "\nReading reference solution for snapshots from %s\n", fullRefName);
  if (snapRefState) delete snapRefState;
  snapRefState = new DistSVec<double, dim>(domain.getNodeDistInfo());
  status = domain.readVectorFromFile(fullRefName, 0, &tmp, *snapRefState);
  delete [] fullRefName;
  fullRefName = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeClusteredOutputs()
{ // This function is called prior to clustering any non-state snapshots (by Model II (online) simulations
  // or, if using snapshot collection method 0, during GNAT preprocessing). It allows the code
  // to accumulate residual and jacaction snapshots over restarts / multiple simulations

  if (strcmp(residualSnapsName,"")!=0) {
    clusterNewtonCount = new int[nClusters];
 
    if (ioData->output.rom.overwriteModelIISnaps == ROMOutputData::OVERWRITE_OFF) {
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
        double tag = 0.0;
        int numSteps = 0;
        int step = 0;
        char *snapshotsPath = 0;
        determinePath(residualSnapsName, iCluster, snapshotsPath); 
        bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);  // if file DNE, returns false, tag=0, and numSteps=0
        clusterNewtonCount[iCluster] = numSteps;
        delete [] snapshotsPath;
        snapshotsPath = NULL;
        if (strcmp(jacActionSnapsName,"")) { // if outputting jac-action snaps, check that num_residuals == num_jac_actions
          determinePath(jacActionSnapsName, iCluster, snapshotsPath);
          tag = 0;
          numSteps = 0;
          bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);
          delete [] snapshotsPath;
          snapshotsPath = NULL;
          if (clusterNewtonCount[iCluster] != numSteps) {
            com->fprintf(stderr, "*** Error: %d residual snapshots found in cluster %d, %d jac-action snapshots found (should match)",
              clusterNewtonCount[iCluster], iCluster, numSteps);
            exit(-1);
          }
        }
      }
    } else {
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) clusterNewtonCount[iCluster] = 0;
    }
  } else if (strcmp(jacActionSnapsName,"")!=0) {
    com->fprintf(stderr, "*** Error: Aero-F assumes that if jacAction snapshots are being collected, then residual snapshots are also being collected");
    exit(-1);
  }

  if (strcmp(krylovSnapsName,"")!=0) {
    clusterKrylovCount = new int[nClusters];
    for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
      double tag = 0.0;
      int numSteps = 0;
      int step = 0;
      char *snapshotsPath = 0;
      determinePath(krylovSnapsName, iCluster, snapshotsPath);
      bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);  // if file DNE, returns false, tag=0, and numSteps=0
      clusterKrylovCount[iCluster] = numSteps;
      delete [] snapshotsPath;
      snapshotsPath = NULL;
    }  
  } 
}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeClusteredBinaryVectors(int iCluster, DistSVec<double,dim> *U1 = NULL, DistSVec<double,dim> *U2 = NULL, 
                                                     DistSVec<double,dim> *U3 = NULL)
{ 
    
  if (strcmp(residualSnapsName,"") && U1)  { // for both pg residuals and fom residuals
    char *residualSnapsPath = 0;
    determinePath(residualSnapsName, iCluster, residualSnapsPath);
    domain.writeVectorToFile(residualSnapsPath, clusterNewtonCount[iCluster], 0.0, *U1);
    delete [] residualSnapsPath;
    residualSnapsPath = NULL;

    if (strcmp(jacActionSnapsName,"") && U2)  {
      char *jacActionSnapsPath = 0;
      determinePath(jacActionSnapsName, iCluster, jacActionSnapsPath);
      domain.writeVectorToFile(jacActionSnapsPath, clusterNewtonCount[iCluster], 0.0, *U2);
      delete [] jacActionSnapsPath;
      jacActionSnapsPath = NULL;
    }
    
    ++(clusterNewtonCount[iCluster]);
  }

  if (strcmp(krylovSnapsName,"") && U3) {
    char *krylovSnapsPath = 0;
    determinePath(krylovSnapsName, iCluster, krylovSnapsPath);
    domain.writeVectorToFile(krylovSnapsPath, clusterKrylovCount[iCluster], 0.0, *U3);
    delete [] krylovSnapsPath;
    krylovSnapsPath = NULL;

    ++(clusterKrylovCount[iCluster]);
  }

  // storing the entire JPhi is no longer supported since it's infeasible in practice

}


//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determineNumResJacMat() { 

  // tests for jacMat; sets numResJacMat
  double tag = 0.0;
  int numSteps = 0;
  int tmp = 0;
  char *jacMatPath = 0;
  determinePath(gappyJacActionName, tmp, jacMatPath);  // check cluster 0
  numResJacMat = (domain.template readTagFromFile<double, dim>(jacMatPath, tmp, &tag, &numSteps)) ? 2 : 1;
  delete [] jacMatPath;
  jacMatPath = NULL;

} 


//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredSampleNodes(int iCluster, bool deleteExistingRestrictionMapping) {

  if (storedAllOnlineQuantities) {
    com->fprintf(stdout, " ... loading sampled nodes for cluster %d\n", iCluster);
    sampleNodes.clear();
    sampleNodes = *(allSampleNodes[iCluster]);
    nSampleNodes = sampleNodes.size();
    restrictionMapping = allRestrictionMappings[iCluster];  //setting pointer
    restrictionMapping->recomputeConnectedTopology();       //need to reinitialize the domain for this restriction mapping
    return;
  }

  char *sampleNodesPath = 0;
  determinePath(sampledNodesName, iCluster, sampleNodesPath);
  FILE *sampleNodeFile = fopen(sampleNodesPath, "r");

  if (!sampleNodeFile)  {
     com->fprintf(stderr, "*** Error: unable to open file %s\n", sampleNodesPath);
     exit (-1);
  }

  int _n;
  _n = fscanf(sampleNodeFile, "%d", &nSampleNodes);  // first entry is the number of sample nodes

  sampleNodes.clear();
  sampleNodes.reserve(nSampleNodes);  // know it will be nSampleNodes long (efficiency)

  int index, currentSampleNode;
  for (int i = 0; i < nSampleNodes; ++i){
    _n = fscanf(sampleNodeFile, "%d", &index);
    _n = fscanf(sampleNodeFile, "%d", &currentSampleNode);
    sampleNodes.push_back(currentSampleNode-1); // reads in the sample node plus one (TODO - change this after code validation)
    if (_n != 1) {
      com->fprintf(stderr, "*** Error: unexpected file format encountered while reading %s\n", sampleNodesPath);
      exit (-1);
    }
  }

  delete [] sampleNodesPath;

  fclose(sampleNodeFile);

  if (restrictionMapping && deleteExistingRestrictionMapping) delete restrictionMapping;
  restrictionMapping = new RestrictionMapping<dim>(&domain, sampleNodes.begin(), sampleNodes.end());

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::deleteRestrictedQuantities() {

  if (resMat) {
    delete resMat;
    resMat = NULL;
  }

  if (jacMat) {
    delete jacMat;
    jacMat = NULL;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredGappyMatrix(int iCluster, char* matrixType) {

  if (storedAllOnlineQuantities) {
    if (strcmp(matrixType,"resMatrix")==0) { 
      com->fprintf(stdout, " ... loading gappy residual matrix for cluster %d\n", iCluster);
      if (resMat) delete resMat;
      int numVecs = (allResMat[iCluster])->numVectors();
      resMat = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());
      for (int iVec=0; iVec<numVecs; ++iVec) 
        (*resMat)[iVec] = (*(allResMat[iCluster]))[iVec];
    } else if (strcmp(matrixType,"jacMatrix")==0) {
      com->fprintf(stdout, " ... loading gappy jacAction matrix for cluster %d\n", iCluster);
      if (jacMat) delete jacMat;       
      int numVecs = (allJacMat[iCluster])->numVectors();
      jacMat = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());
      for (int iVec=0; iVec<numVecs; ++iVec) 
        (*jacMat)[iVec] = (*(allJacMat[iCluster]))[iVec]; 
    }
    return;
  }

  char* matrixPath = 0;

  if (strcmp(matrixType,"resMatrix")==0) {
      if (resMat) {
        delete resMat;
        resMat = NULL;
      }
      determinePath(gappyResidualName, iCluster, matrixPath);
      com->fprintf(stdout, "\nReading gappy residual matrix for cluster %d\n", iCluster);
  } else if (strcmp(matrixType,"jacMatrix")==0) {
      if (jacMat) {
        delete jacMat;
        jacMat = NULL;
      }
      determinePath(gappyJacActionName, iCluster, matrixPath);
      com->fprintf(stdout, "\nReading gappy jacAction matrix for cluster %d\n", iCluster);
  } else {
      exit (-1);
  }

  double tag = 0.0;
  int numVecs = 0;
  int step = 0;
  bool status = domain.readTagFromFile<double, dim>(matrixPath, step, &tag, &numVecs);  // if file DNE, returns false, tag=0, and numSteps=0

  if (!status) {
    com->fprintf(stderr, "\nCould not open file %s\n", matrixPath);
    exit(-1);
  }

  VecSet<DistSVec<double, dim> > gappyMatrix(numVecs, domain.getNodeDistInfo());

  for (int iVec=0; iVec<numVecs; ++iVec) {
    status = domain.readVectorFromFile(matrixPath, iVec, &tag, gappyMatrix[iVec]);
  }

  delete [] matrixPath;

  // restrictionMapping is set during readClusteredSampleNodes
  VecSet<DistSVec<double, dim> >* restrictedGappyMatrix = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());

  for (int i = 0; i < numVecs; ++i) {
    restrictionMapping->restriction(gappyMatrix[i],(*restrictedGappyMatrix)[i]);
  }

  if (strcmp(matrixType,"resMatrix")==0) {
      resMat = restrictedGappyMatrix;
  } else if (strcmp(matrixType,"jacMatrix")==0) {
      jacMat = restrictedGappyMatrix;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readAllOnlineQuantities() {

// stores all online quantities at the beginning of the simulation
// (as opposed to reading information at each cluster switch)

// note: fast distance calculation info is handled separately

// initial allocation

  allStateBases = new VecSet< DistSVec<double, dim> >*[nClusters];
  allStateSVals = new std::vector<double>*[nClusters];

  if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
    allKrylovBases = new VecSet< DistSVec<double, dim> >*[nClusters];
    allKrylovSVals = new std::vector<double>*[nClusters];
  }

  if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
    allRefStates = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
    allColumnSumsV = new std::vector<double>*[nClusters];
    if (false) { 
      // additional update info
    }
  }

  if (ioData->romOnline.systemApproximation==NonlinearRomOnlineData::GNAT){
    allSampleNodes = new std::vector<int>*[nClusters];
    allRestrictionMappings = new RestrictionMapping<dim>*[nClusters];
    allResMat = new VecSet< DistSVec<double, dim> >*[nClusters];
    if (numResJacMat==2) allJacMat = new VecSet< DistSVec<double, dim> >*[nClusters];
  }

// read and store online info for each cluster

  switch (ioData->romOnline.systemApproximation) {
    case (NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE):
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {

        // read state ROB and sVals
        readClusteredBasis(iCluster, "state");
        allStateBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*(allStateBases[iCluster]))[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        allStateSVals[iCluster] = new vector<double>;
        *(allStateSVals[iCluster]) = *sVals;
        delete sVals;
        sVals = NULL;

        // read update info
        if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
          readClusteredUpdateInfo(iCluster, "state");
          (*allRefStates)[iCluster] = *Uref;
          delete Uref;
          Uref = NULL;
          allColumnSumsV[iCluster] = new vector<double>;
          *(allColumnSumsV[iCluster]) = *columnSumsV;
          delete columnSumsV;
          columnSumsV = NULL;
          if (false) { 
            // read additional update information
          }
        }

        // read Krylov ROB and sVals
        if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
          readClusteredBasis(iCluster, "krylov");
          allKrylovBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
          for (int iVec=0; iVec<basis->numVectors(); ++iVec)
            (*(allKrylovBases[iCluster]))[iVec] = (*basis)[iVec];
          delete basis;
          basis = NULL;
          allKrylovSVals[iCluster] = new vector<double>;
          *(allKrylovSVals[iCluster]) = *sVals;
          delete sVals;
          sVals = NULL;
        }

      }

      // read sensitivity ROB and sVals
      if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
        readClusteredBasis(-2, "sensitivity");
        sensitivityBasis = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*sensitivityBasis)[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        sensitivitySVals = new vector<double>;
        *sensitivitySVals = *sVals;
        delete sVals;
        sVals = NULL;
      }
            
      break;
    case (NonlinearRomOnlineData::GNAT):
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        // read sample nodes
        readClusteredSampleNodes(iCluster, false); // resets restriction map
        allSampleNodes[iCluster] = new vector<int>;
        *(allSampleNodes[iCluster]) = sampleNodes;
        allRestrictionMappings[iCluster] = restrictionMapping; // sets pointer to dynamically allocated memory
         
        // read gappy POD matrix for residual
        readClusteredGappyMatrix(iCluster, "resMatrix");
        allResMat[iCluster] = new VecSet< DistSVec<double, dim> >(resMat->numVectors(), restrictionMapping->restrictedDistInfo());
        for (int iVec=0; iVec<resMat->numVectors(); ++iVec)
          (*(allResMat[iCluster]))[iVec] = (*resMat)[iVec];
        delete resMat;
        resMat = NULL;

        // read gappy POD matrix for jacobian
        if (numResJacMat==2) {
          readClusteredGappyMatrix(iCluster, "jacMatrix");
          allJacMat[iCluster] = new VecSet< DistSVec<double, dim> >(jacMat->numVectors(), restrictionMapping->restrictedDistInfo());
          for (int iVec=0; iVec<jacMat->numVectors(); ++iVec)
            (*(allJacMat[iCluster]))[iVec] = (*jacMat)[iVec];
          delete jacMat;
          jacMat = NULL;
        }

        // read sampled state ROB and sVals
        readClusteredBasis(iCluster, "sampledState");
        allStateBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*(allStateBases[iCluster]))[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        allStateSVals[iCluster] = new vector<double>;
        *(allStateSVals[iCluster]) = *sVals;
        delete sVals;
        sVals = NULL;

        // read ROB update information
        if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
          readClusteredUpdateInfo(iCluster, "sampledState");
          (*allRefStates)[iCluster] = *Uref;
          delete Uref;
          Uref = NULL;
          allColumnSumsV[iCluster] = new vector<double>;
          *(allColumnSumsV[iCluster]) = *columnSumsV;
          delete columnSumsV;
          columnSumsV = NULL;
          if (false) {
            // read additional update information
          }
        }

        // read sampled Krylov ROB
        if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) { 
          readClusteredBasis(iCluster, "sampledKrylov");
          allKrylovBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
          for (int iVec=0; iVec<basis->numVectors(); ++iVec)
            (*(allKrylovBases[iCluster]))[iVec] = (*basis)[iVec];
          delete basis;
          basis = NULL;
          allKrylovSVals[iCluster] = new vector<double>;
          *(allKrylovSVals[iCluster]) = *sVals;
          delete sVals;
          sVals = NULL;
        }
 
      }

      // read sampled sensitivity ROB
      if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
        readClusteredBasis(-2, "sampledSensitivity");
        sensitivityBasis = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*sensitivityBasis)[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        sensitivitySVals = new vector<double>;
        *sensitivitySVals = *sVals;
        delete sVals;
        sVals = NULL;
      }

      break;
    default:
      com->fprintf(stderr, "*** Error:  Unexpected system approximation type\n");
      exit(-1);
  }

  storedAllOnlineQuantities = true;

}


//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeReducedCoords(const int totalTimeSteps, bool clusterSwitch, bool updateFreq, int iCluster, Vec<double> dUromTimeIt) {

  int nPod = basis->numVectors();

  if (this->com->cpuNum() == 0) {
    if (clustUsageFile)  // for plotting only
      this->com->fprintf(clustUsageFile,"%d %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
    if (reducedCoordsFile) {
      if (clusterSwitch) {
        if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
          this->com->fprintf(reducedCoordsFile,"%d switch update %d %d %d %d %d\n",
            totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        } else {
          this->com->fprintf(reducedCoordsFile,"%d switch noUpdate %d %d %d %d %d\n",
            totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        }
      } else if (updateFreq && (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF)) {
        this->com->fprintf(reducedCoordsFile,"%d noSwitch update %d %d %d %d %d\n",
          totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      } else { 
        this->com->fprintf(reducedCoordsFile,"%d noSwitch noUpdate %d %d %d %d %d\n",
          totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      }
      for (int iPod=0; iPod<nPod; ++iPod) {
        this->com->fprintf(reducedCoordsFile, "%23.15e\n", dUromTimeIt[iPod]);
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputCenterNorms(std::vector<std::vector<double> > &vec) {

  char *infoPath = 0;
  determinePath(centerNormsName, -1, infoPath);   

  int nFluidMeshNodes = domain.getNumGlobNode();    

  if (com->cpuNum() == 0) {
    // format:  numFullMeshNodes numClusters vecSize; element1_1; element1_2; ...
    FILE *outputFile = fopen(infoPath, "wt");

    com->fprintf(outputFile, "%d %d %d\n", nFluidMeshNodes, nClusters, vec[0].size());

    for (int iCenter=0; iCenter<nClusters; ++iCenter) {
      for (int iVec=0; iVec<vec[iCenter].size(); ++iVec) {
        com->fprintf(outputFile,"%23.15e\n", vec[iCenter][iVec]);
      }
    }
 
    fclose (outputFile);
  }

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readCenterNorms() {

  centerNorms.clear();

  char *infoPath = 0;
  determinePath(centerNormsName, -1, infoPath);   

  // format:  numFullMeshNodes numClusters vecSize; element1_1; element1_2; ...
  FILE *inputFile = fopen(infoPath, "r");
  int _n;

  int expectedNClusters;
  int vecSize;

  _n = fscanf(inputFile, "%d %d %d\n", &nFullMeshNodes, &expectedNClusters, &vecSize);

  assert(expectedNClusters == nClusters);

  if (vecSize>1) {
    specifiedIC=false; 
    assert(vecSize==(dim+1));
  }

  centerNorms.resize(nClusters);
  for (int iCenter=0; iCenter<nClusters; ++iCenter)
    centerNorms[iCenter].reserve(vecSize);

  double tmpVal;

  for (int iCenter=0; iCenter<nClusters; ++iCenter) {
    for (int iVec=0; iVec<vecSize; ++iVec) {
      _n = fscanf(inputFile,"%le", &tmpVal);
      if (_n == 1) {
        centerNorms[iCenter].push_back(tmpVal);
      } else if (feof(inputFile)) {
        break;
      } else {
        com->fprintf(stderr, "*** Error: fscanf interrupted by non-EOF error\n");
        exit(-1);
      }
    }
  }

  fclose (inputFile);

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredInfoASCII(int iCluster, char* type, std::vector<double>* vec1,
                                             std::vector<std::vector<double> >* vec2,
                                             std::vector<std::vector<std::vector<double> > >* vec3,
                                             std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // interface for outputting small ASCII quantities

  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    com->fprintf(stderr, "*** Error: only one precomputed distance comparison quantity can be output at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    com->fprintf(stderr, "*** Error: no precomputed distance comparison quantity specified\n");
    exit(-1); 
  }

  char *infoPath = NULL;
  if (strcmp(type, "refState") == 0) { // 2*(U_center_p - U_center_m)^T U_ref
    determinePath(stateDistanceComparisonInfoExactUpdatesName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "initialCondition") == 0) { // 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "state") == 0) { // 2*(U_center_p - U_center_m)^T V_state_k
    determinePath(stateDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "krylov") == 0) { // 2*(U_center_p - U_center_m)^T V_krylov_k
    determinePath(krylovDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "sensitivity") == 0) { // 2*(U_center_p - U_center_m)^T V_sens
    determinePath(sensitivityDistanceComparisonInfoName, -2, infoPath);
    assert(vec3);
  } else if (strcmp(type, "distanceMatrix") == 0) { // A_ij = ||U_i - U_j||_2 
    determinePath(distanceMatrixName, -1, infoPath);
    assert(vec2);
  /*} else if (strcmp(type, "basisBasisProducts") == 0) {               
    determinePath(basisBasisProductsName, -1, infoPath);
    assert(vec4);
  } else if (strcmp(type, "basisUrefProducts") == 0) {               
    determinePath(basisUrefProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "normUic") == 0) {               
    determinePath(normUicName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "normUref") == 0) {               
    determinePath(normUrefName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "productUrefUic") == 0) {               
    determinePath(productUrefUicName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "basisUicProducts") == 0) {               
    determinePath(basisUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "urefComponentwiseSums") == 0) {               
    determinePath(urefComponentwiseSumsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisComponentwiseSums") == 0) {               
    determinePath(basisComponentwiseSumsName, -1, infoPath);
    assert(vec3);*/
  } else {
    exit(-1);
  }
  
  writeMultiVecASCII(infoPath, vec1, vec2, vec3, vec4);

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredInfoASCII(int iCluster, char* type, std::vector<double>* vec1,
                                             std::vector<std::vector<double> >* vec2,
                                             std::vector<std::vector<std::vector<double> > >* vec3,
                                             std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {

  // interface for reading small ASCII quantities
                                    
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    com->fprintf(stderr, "*** Error: only one precomputed distance comparison quantity can be output at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    com->fprintf(stderr, "*** Error: no precomputed distance comparison quantity specified\n");
    exit(-1); 
  }

  char *infoPath = NULL;

  if (strcmp(type, "refState") == 0) {// 2*(U_center_p - U_center_m)^T U_ref
    determinePath(stateDistanceComparisonInfoExactUpdatesName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "initialCondition") == 0) {// 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "state") == 0) {// 2*(U_center_p - U_center_m)^T V_state_k
    determinePath(stateDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "krylov") == 0) {// 2*(U_center_p - U_center_m)^T V_krylov_k
    determinePath(krylovDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "sensitivity") == 0) {// 2*(U_center_p - U_center_m)^T V_sens
    determinePath(sensitivityDistanceComparisonInfoName, -2, infoPath);
    assert(vec3);
  /*} else if (strcmp(type, "basisBasisProducts") == 0) {
    determinePath(basisBasisProductsName, -1, infoPath);
    assert(vec4);
  } else if (strcmp(type, "basisUrefProducts") == 0) {
    determinePath(basisUrefProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "normUic") == 0) {
    determinePath(normUicName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "normUref") == 0) {               
    determinePath(normUrefName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "productUrefUic") == 0) {               
    determinePath(productUrefUicName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "basisUicProducts") == 0) {               
    determinePath(basisUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "urefComponentwiseSums") == 0) {               
    determinePath(urefComponentwiseSumsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisComponentwiseSums") == 0) {               
    determinePath(basisComponentwiseSumsName, -1, infoPath);
    assert(vec3);*/
  } else {
    exit(-1);
  }

  readMultiVecASCII(infoPath, vec1, vec2, vec3, vec4);

  delete [] infoPath;
  infoPath = NULL;

} 

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeMultiVecASCII(char* path, std::vector<double>* vec1,
                                           std::vector<std::vector<double> >* vec2,
                                           std::vector<std::vector<std::vector<double> > >* vec3,
                                           std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // general IO function for multi-vector quantities with arbitrary dimensions
 
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    com->fprintf(stderr, "*** Error: only one multivec can be output at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    com->fprintf(stderr, "*** Error: no multivec specified\n");
    exit(-1); 
  }

  if (com->cpuNum() == 0) {

    int dim1, dim2, dim3, dim4, multiVecType;  

    FILE *outputFile = fopen(path, "wt");

    if (vec1) {
      dim1 = vec1->size();
      multiVecType = 1;
    } else if (vec2) {
      dim1 = vec2->size();
      multiVecType = 2;
    } else if (vec3) {
      dim1 = vec3->size();
      multiVecType = 3;
    } else {
      dim1 = vec4->size();
      multiVecType = 4;
    }

    com->fprintf(outputFile, "MultiVecType: %d\n", multiVecType);
    com->fprintf(outputFile, "Dimension#1: %d\n", dim1);
    for (int i=0; i<dim1; ++i){
      if (vec1) {
        com->fprintf(outputFile,"%23.15e\n", (*vec1)[i]);
      } else {
        if (vec2) dim2 = (*vec2)[i].size();
        if (vec3) dim2 = (*vec3)[i].size();
        if (vec4) dim2 = (*vec4)[i].size();
        com->fprintf(outputFile, "Dimension#2: %d\n", dim2);
        for (int j=0; j<dim2; ++j) {
          if (vec2) {
            com->fprintf(outputFile,"%23.15e\n", (*vec2)[i][j]);
          } else {
            if (vec3) dim3 = (*vec3)[i][j].size();
            if (vec4) dim3 = (*vec4)[i][j].size();
            com->fprintf(outputFile, "Dimension#3: %d\n", dim3);
            for (int k=0; k<dim3; ++k) {
              if (vec3) {
                com->fprintf(outputFile,"%23.15e\n", (*vec3)[i][j][k]);
              } else {
                dim4 = (*vec4)[i][j][k].size();
                com->fprintf(outputFile, "Dimension#4: %d\n", dim4);
                for (int l=0; l<dim4; ++l) {
                  com->fprintf(outputFile,"%23.15e\n", (*vec4)[i][j][k][l]);
                }
              }
            }
          }
        }
      }
    }
  
    fclose (outputFile);

  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readMultiVecASCII(char* path, std::vector<double>* vec1,
                                          std::vector<std::vector<double> >* vec2,
                                          std::vector<std::vector<std::vector<double> > >* vec3,
                                          std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // general IO function for multi-vector quantities with arbitrary dimensions
                                    
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    com->fprintf(stderr, "*** Error: only one multivec can be read at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    com->fprintf(stderr, "*** Error: no multivec specified\n");
    exit(-1); 
  }

  FILE *inputFile = fopen(path, "r");

  int expectedMultiVecType;

  if (vec1) { 
    vec1->clear();
    expectedMultiVecType = 1;
  } else if (vec2) {
    vec2->clear();
    expectedMultiVecType = 2;
  } else if (vec3) {
    vec3->clear();
    expectedMultiVecType = 3;
  } else {
    vec4->clear();
    expectedMultiVecType = 4;
  }

  int _n, dim1, dim2, dim3, dim4, multiVecType;
  double tmpVal;

  fscanf(inputFile, "MultiVecType: %d\n", &multiVecType);
  assert(expectedMultiVecType == multiVecType);
  fscanf(inputFile, "Dimension#1: %d\n", &dim1);

  if (vec1) vec1->reserve(dim1);
  if (vec2) vec2->resize(dim1);
  if (vec3) vec3->resize(dim1);
  if (vec4) vec4->resize(dim1);

  for (int i=0; i<dim1; ++i){
    if (vec1) {
      _n = fscanf(inputFile,"%le\n", &tmpVal);
      vec1->push_back(tmpVal);
    } else {
      _n = fscanf(inputFile, "Dimension#2: %d\n", &dim2);
      if (vec2) (*vec2)[i].reserve(dim2);
      if (vec3) (*vec3)[i].resize(dim2);
      if (vec4) (*vec4)[i].resize(dim2);
      for (int j=0; j<dim2; ++j) {
        if (vec2) {
          _n = fscanf(inputFile,"%le\n", &tmpVal);
          (*vec2)[i].push_back(tmpVal);
        } else {
          _n = fscanf(inputFile, "Dimension#3: %d\n", &dim3);
          if (vec3) (*vec3)[i][j].reserve(dim3);
          if (vec4) (*vec4)[i][j].resize(dim3);
          for (int k=0; k<dim3; ++k) {
            if (vec3) {
              _n = fscanf(inputFile,"%le\n", &tmpVal);
              (*vec3)[i][j].push_back(tmpVal);
            } else {
              _n = fscanf(inputFile, "Dimension#4: %d\n", &dim4);
              (*vec4)[i][j][k].reserve(dim4);
              for (int l=0; l<dim4; ++l) {
                com->fprintf(inputFile,"%23.15e\n", (*vec4)[i][j][k][l]);
                _n = fscanf(inputFile,"%le\n", &tmpVal);
                (*vec4)[i][j][k].push_back(tmpVal);
              }
            }
          }
        }
      }
    }
  }

  fclose (inputFile);

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readDistanceComparisonInfo(char* updateType) {

  readCenterNorms();

  if ((strcmp(updateType, "noUpdates") == 0) || (strcmp(updateType, "exactUpdates") == 0)){

    stateBasisCentersProduct.resize(nClusters);

    if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
      krylovBasisCentersProduct.resize(nClusters); 

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      readClusteredInfoASCII(iCluster, "state", NULL, NULL, &stateBasisCentersProduct[iCluster]);
    
      if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
        readClusteredInfoASCII(iCluster, "krylov", NULL, NULL, &krylovBasisCentersProduct[iCluster]); 
    }

    if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
      readClusteredInfoASCII(-2, "sensitivity", NULL, NULL, &sensitivityBasisCentersProduct);      

    if (strcmp(updateType, "exactUpdates") == 0) {
      if (specifiedIC) {
        readClusteredInfoASCII(-1, "initialConditon", NULL, &initialConditionCentersProduct);
      } else {
        // this will be constructed during initializeDistanceComparisons 
      }

      refStateCentersProduct.resize(nClusters);
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        readClusteredInfoASCII(iCluster, "refState", NULL, &refStateCentersProduct[iCluster]);
      }
    }
  } else if (strcmp(updateType, "approxUpdates") == 0) {   


  } else {
    com->fprintf(stderr, "*** Error: unexpected update method\n");
    exit(-1);
  }

}

//------------------------------------------------------------------------------

/*
  if (staterom) {
    if (it0 != 0)
      fpStateRom = backupAsciiFile(staterom);
    if (it0 == 0 || fpStateRom == 0) {
      fpStateRom = fopen(staterom, "w");
      if (!fpStateRom) {
  fprintf(stderr, "*** Error: could not open \'%s\'\n", staterom);
  exit(1);
      }
      fprintf(fpStateRom, "# TimeIteration ElapsedTime StatePodCoords \n");
    }
    fflush(fpStateRom);
  }


  if (iod.output.rom.staterom[0] != 0) {
    staterom = new char[sprom + strlen(iod.output.rom.staterom)];
    sprintf(staterom, "%s%s", iod.output.rom.prefix, iod.output.rom.staterom);
  }
  else
    staterom = 0;

  if (fpStateRom) fclose(fpStateRom);
*/

