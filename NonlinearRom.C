#include <TsInput.h>
#include <math.h>
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

  // Directory information
  databasePrefix = ioData->romDatabase.directories.prefix;
  databaseName = ioData->romDatabase.directories.databaseName;
  clusterName = ioData->romDatabase.directories.clusterName;
  sensitivityClusterName  = ioData->romDatabase.directories.sensitivityClusterName;

  // State snapshot clusters
  determineFileName(ioData->romDatabase.files.stateSnapsName, "snaps", ioData->romDatabase.files.statePrefix, stateSnapsName);
  determineFileName(ioData->romDatabase.files.mapName, "map", ioData->romDatabase.files.statePrefix, mapName);
  determineFileName(ioData->romDatabase.files.indexName, "index", ioData->romDatabase.files.statePrefix, indexName);
  determineFileName(ioData->romDatabase.files.connName, "conn", ioData->romDatabase.files.statePrefix, connName);
  determineFileName(ioData->romDatabase.files.centersName, "centers", ioData->romDatabase.files.statePrefix, centersName);
  determineFileName(ioData->romDatabase.files.nearestName, "nearest", ioData->romDatabase.files.statePrefix, nearestName);

  // State bases
  determineFileName(ioData->romDatabase.files.stateBasisPrefix, "", ioData->romDatabase.files.statePrefix, stateBasisPrefix);
  determineFileName(ioData->romDatabase.files.stateBasisName, "rob", stateBasisPrefix, stateBasisName);
  determineFileName(ioData->romDatabase.files.stateSingValsName, "svals", stateBasisPrefix, stateSingValsName);
  determineFileName(ioData->romDatabase.files.updateInfoName, "update", stateBasisPrefix, updateInfoName);
  determineFileName(ioData->romDatabase.files.stateFastDistCalcInfoName, "dist", stateBasisPrefix, stateFastDistCalcInfoName);
  determineFileName(ioData->romDatabase.files.projErrorName, "proj", stateBasisPrefix, projErrorName);
  determineFileName(ioData->romDatabase.files.refStateName, "refState", stateBasisPrefix, refStateName); 

  // Krylov snaps
  determineFileName(ioData->romDatabase.files.krylovSnapsName, "snaps", ioData->romDatabase.files.krylovPrefix, krylovSnapsName);

  // Krylov bases
  determineFileName(ioData->romDatabase.files.krylovBasisPrefix, "", ioData->romDatabase.files.krylovPrefix, krylovBasisPrefix);
  determineFileName(ioData->romDatabase.files.krylovBasisName, "rob", krylovBasisPrefix, krylovBasisName);
  determineFileName(ioData->romDatabase.files.krylovSingValsName, "svals", krylovBasisPrefix, krylovSingValsName);
  determineFileName(ioData->romDatabase.files.krylovFastDistCalcInfoName, "dist", krylovBasisPrefix, krylovFastDistCalcInfoName);

  // Sensitivity snaps
  determineFileName(ioData->romDatabase.files.sensitivitySnapsName, "snaps", ioData->romDatabase.files.sensitivityPrefix, sensitivitySnapsName);

  // Sensitivity basis
  determineFileName(ioData->romDatabase.files.sensitivityBasisPrefix, "", ioData->romDatabase.files.sensitivityPrefix, sensitivityBasisPrefix);
  determineFileName(ioData->romDatabase.files.sensitivityBasisName, "rob", sensitivityBasisPrefix, sensitivityBasisName);
  determineFileName(ioData->romDatabase.files.sensitivitySingValsName, "svals", sensitivityBasisPrefix, sensitivitySingValsName);

  // Residual snaps
  determineFileName(ioData->romDatabase.files.residualSnapsName, "snaps", ioData->romDatabase.files.residualPrefix, residualSnapsName);

  // Residual bases
  determineFileName(ioData->romDatabase.files.residualBasisPrefix, "", ioData->romDatabase.files.residualPrefix, residualBasisPrefix);
  determineFileName(ioData->romDatabase.files.residualBasisName, "rob", residualBasisPrefix, residualBasisName);
  determineFileName(ioData->romDatabase.files.residualSingValsName, "svals", residualBasisPrefix, residualSingValsName);

  // Action-of-Jacobian snaps
  determineFileName(ioData->romDatabase.files.jacActionSnapsName, "snaps", ioData->romDatabase.files.jacActionPrefix, jacActionSnapsName);

  // Action-of-Jacobian bases
  determineFileName(ioData->romDatabase.files.jacActionBasisPrefix, "", ioData->romDatabase.files.jacActionPrefix, jacActionBasisPrefix);
  determineFileName(ioData->romDatabase.files.jacActionBasisName, "rob", jacActionBasisPrefix, jacActionBasisName);
  determineFileName(ioData->romDatabase.files.jacActionSingValsName, "svals", jacActionBasisPrefix, jacActionSingValsName);

  // GNAT quantities
  determineFileName(ioData->romDatabase.files.sampledNodesName, "sampledNodes", ioData->romDatabase.files.gnatPrefix, sampledNodesName);
  determineFileName(ioData->romDatabase.files.sampledNodesFullCoordsName, "sampledNodesFullCoords", ioData->romDatabase.files.gnatPrefix, sampledNodesFullCoordsName);
  determineFileName(ioData->romDatabase.files.sampledStateBasisName, "sampledStateROB", ioData->romDatabase.files.gnatPrefix, sampledStateBasisName);
  determineFileName(ioData->romDatabase.files.sampledResidualBasisName, "sampledResROB", ioData->romDatabase.files.gnatPrefix, sampledResidualBasisName);
  determineFileName(ioData->romDatabase.files.sampledJacActionBasisName, "sampledJacROB", ioData->romDatabase.files.gnatPrefix, sampledJacActionBasisName);
  determineFileName(ioData->romDatabase.files.sampledMeshName, "top", ioData->romDatabase.files.gnatPrefix, sampledMeshName);
  determineFileName(ioData->romDatabase.files.sampledSolutionName, "sol", ioData->romDatabase.files.gnatPrefix, sampledSolutionName);
  determineFileName(ioData->romDatabase.files.sampledWallDistName, "dwall", ioData->romDatabase.files.gnatPrefix, sampledWallDistName);
  determineFileName(ioData->romDatabase.files.gappyJacActionName, "gappyJac", ioData->romDatabase.files.gnatPrefix, gappyJacActionName);
  determineFileName(ioData->romDatabase.files.gappyResidualName, "gappyRes", ioData->romDatabase.files.gnatPrefix, gappyResidualName);

  // Surface quantities
  determineFileName(ioData->romDatabase.files.gappyResidualName, "sampledStateROB", ioData->romDatabase.files.surfacePrefix, surfaceStateBasisName);
  determineFileName(ioData->romDatabase.files.gappyResidualName, "sol", ioData->romDatabase.files.surfacePrefix, surfaceSolutionName);
  determineFileName(ioData->romDatabase.files.gappyResidualName, "dwall", ioData->romDatabase.files.surfacePrefix, surfaceWallDistName);
  determineFileName(ioData->romDatabase.files.gappyResidualName, "top", ioData->romDatabase.files.surfacePrefix, surfaceMeshName);


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
  Uinit = NULL;
  clusterNewtonCount = NULL;

  nSampleNodes = 0;
  sampleNodes.clear();
  numResJacMat = 0;
  resMat = NULL;
  jacMat = NULL;

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
  delete [] updateInfoName;
  delete [] stateFastDistCalcInfoName;
  delete [] refStateName;
  delete [] projErrorName;
  delete [] krylovSnapsName;
  delete [] krylovBasisName;
  delete [] krylovSingValsName;
  delete [] krylovFastDistCalcInfoName;
  delete [] residualSnapsName;
  delete [] jacActionSnapsName;
  delete [] residualBasisName;
  delete [] residualSingValsName;
  delete [] jacActionBasisName;
  delete [] jacActionSingValsName;
  delete [] sampledNodesName;
  delete [] sampledStateBasisName;
  delete [] sampledResidualBasisName;
  delete [] sampledJacActionBasisName;
  delete [] sampledMeshName;
  delete [] sampledSolutionName;
  delete [] sampledWallDistName;
  delete [] gappyJacActionName;
  delete [] gappyResidualName;
  delete [] surfaceStateBasisName;
  delete [] surfaceSolutionName;
  delete [] surfaceWallDistName;
  delete [] surfaceMeshName;

  delete [] stateBasisPrefix;
  delete [] krylovBasisPrefix;
  delete [] sensitivityBasisPrefix;
  delete [] residualBasisPrefix;
  delete [] jacActionBasisPrefix;

  if (snap) delete snap; 
  if (clusterCenters) delete clusterCenters;
  if (nearestSnapsToCenters) delete nearestSnapsToCenters;
  if (snapsInCluster) delete [] snapsInCluster;
  if (clusterIndex) delete [] clusterIndex;
  if (clusterNeighborsCount) delete [] clusterNeighborsCount;
  if (clusterNewtonCount) delete [] clusterNewtonCount;
  if (clusterKrylovCount) delete [] clusterKrylovCount;
  if (snapRefState) delete snapRefState;
  if (columnSumsV) delete [] columnSumsV; 
  if (sVals) delete [] sVals;
  if (Uinit) delete Uinit;
  
  //TODO
  //clusterSnapshotMap
  //clusterNeighbors

  if (resMat) delete resMat;
  if (jacMat) delete jacMat;
  restrictionMapping.reset();

  if (this->com->cpuNum() == 0) {
    if (clustUsageFile) fclose(clustUsageFile);
    if (reducedCoordsFile) fclose(reducedCoordsFile);
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
    tmp = distance( vec, (*clusterCenters)[iCluster]);
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
double NonlinearRom<dim>::distance(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {

// TODO: Implement several different distance functions
  DistSVec<double, dim> diff(domain.getNodeDistInfo());
  diff = (U1 - U2);
  double dist = diff.norm();

  return dist;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusterSnapshots(char* snapType)  { 
 
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
      com->fprintf(clustIndexFile,"Snapshot# OriginalCluster\n");
  
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
void NonlinearRom<dim>::readClusterSnapshots(int iCluster, bool preprocess, char *basisType, int first, int last) {

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

  snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());

  // read in Snapshot Vectors
  double tmp;
  bool moreVecs = true; 
  int snapCount = 0;
  int snapIndex = first; 
  DistSVec<double, dim>* tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {
    moreVecs = domain.readVectorFromFile(snapshotsPath, snapIndex, &tmp, *tmpSnap);
    if (!moreVecs) {
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
        outputClusterReferenceState(iCluster, (*nearestSnapsToCenters)[iCluster]);
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
        outputClusterReferenceState(iCluster,(*clusterCenters)[iCluster]);
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
        outputClusterReferenceState(iCluster,(*snapRefState));
        com->barrier();
        delete snapRefState;
        snapRefState = NULL;
      } else {
        if (ioData->romOffline.rob.state.snapshots.incrementalSnaps)
          com->fprintf(stderr, "*** Warning: incremental snapshots not supported for local ROM construction\n");
        com->fprintf(stdout, " ... using raw states as snapshots\n");
        DistSVec<double, dim> zeroVec(domain.getNodeDistInfo(), 0);
        outputClusterReferenceState(iCluster, zeroVec);
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
void NonlinearRom<dim>::outputClusterReferenceState(int iCluster, DistSVec<double, dim>& ref) {
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
void NonlinearRom<dim>::readClusterReferenceState(int iCluster) {
  // This function reads in the automatically stored reference snapshot for a cluster.
  // By storing these reference snapshots and reading them automatically it reduces the user's workload.

  char *refStatePath = 0;
  determinePath(refStateName, iCluster, refStatePath);

  double tmp;
  bool moreVecs;
 
  com->fprintf(stdout, "Reading reference snapshot %s\n", refStatePath);
  Uinit = new DistSVec<double, dim>(domain.getNodeDistInfo());
  moreVecs = domain.readVectorFromFile(refStatePath, 0, &tmp, *Uinit);

  delete [] refStatePath;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusterCenters() {

  char *clustCentersPath = 0;
  determinePath(centersName, -1, clustCentersPath);

  if (nClusters <= 0) {
    com->fprintf(stderr, "\n*** Error: invalid value for NumClusters (%d)\n", nClusters);
    exit(-1);
  }

  // read centers
  clusterCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
  snapsInCluster = new int[nClusters];

  com->fprintf(stdout, "\nReading cluster centers\n");
  bool moreVecs;
  int expectedClusters = nClusters;
  nClusters = 0;
  double tmp;
  DistSVec<double, dim>* tmpVec = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {

    moreVecs = domain.readVectorFromFile(clustCentersPath, nClusters, &tmp, *tmpVec);
    if (!moreVecs) break;
 
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
  // NOTE: this function must only be run after readClusterCenters (this is because the
  // cluster centers file defines snapsInCluster and nClusters)

  char *nearestSnapsPath = 0;
  determinePath(nearestName, -1, nearestSnapsPath);
 
  nearestSnapsToCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  bool moreVecs;

  com->fprintf(stdout, "\nReading closest snapshot to each center\n");

  // read in Snapshot Vectors
  double tmp;
  for (int iCluster = 0; iCluster<nClusters; ++iCluster) {
    moreVecs = domain.readVectorFromFile(nearestSnapsPath, iCluster, &tmp, (*nearestSnapsToCenters)[iCluster]);
    if (!moreVecs) {
       com->fprintf(stderr,"*** Error: readNearestSnapsToCenters attempted to read %d vecs from a file with only %d.\n", nClusters, iCluster-1);
       exit(-1);
    }
  }

  delete [] nearestSnapsPath;
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusterBasis(int iCluster, char* basisType) {

  if (basis) {
    delete basis;
    basis = NULL;
  }

  int maxDimension;
  int minDimension;
  double energyTol;

  if (ioData->problem.type[ProblemData::NLROMOFFLINE]) {
    // GNAT preprocessing TODO - add krylov and sensitivities
    if ((ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_) ||
        (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_1_) ||
        (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_2_)) {
      if (strcmp(basisType,"state")==0) {
        if (ioData->romOffline.gnat.maxDimensionState <= 0) {
          com->fprintf(stderr, "*** Error: Maximum dimension for state basis  must be specified (and greater than zero)\n");
          exit (-1);
        } else {
          maxDimension = ioData->romOffline.gnat.maxDimensionState;
          minDimension = ioData->romOffline.gnat.minDimensionState;
          energyTol = ioData->romOffline.gnat.energyState;
        }
      } else if (strcmp(basisType,"residual")==0) {
        if (ioData->romOffline.gnat.maxDimensionResidual <= 0) {
          com->fprintf(stderr, "*** Error: Maximum dimension for residual basis must be specified (and greater than zero)\n");
          exit (-1);
        } else {
          maxDimension = ioData->romOffline.gnat.maxDimensionResidual;
          minDimension = ioData->romOffline.gnat.minDimensionResidual;
          energyTol = ioData->romOffline.gnat.energyResidual;
        }
      } else if (strcmp(basisType,"jacAction")==0) {
        if (ioData->romOffline.gnat.maxDimensionJacAction <= 0) {
          if (ioData->romOffline.gnat.maxDimensionResidual <= 0) {  // defaults to the value of maxDimensionResidual
            com->fprintf(stderr, "*** Error: Maximum dimension for residual must be specified (and greater than zero)\n");
            exit (-1);
          } else {
            com->fprintf(stdout, "*** Warning: JacAction greedy parameters not specified. Using Residual paraemters for both bases\n");
            maxDimension = ioData->romOffline.gnat.maxDimensionResidual;
            minDimension = ioData->romOffline.gnat.minDimensionResidual;
            energyTol = ioData->romOffline.gnat.energyResidual;
          }
        } else {
          maxDimension = ioData->romOffline.gnat.maxDimensionJacAction;
          minDimension = ioData->romOffline.gnat.minDimensionJacAction;
          energyTol = ioData->romOffline.gnat.energyJacAction;
        }
      }
    } else { /* //TODO  relProjError
      if (ioData->romOffline.rob.dataCompression.maxBasisSize <= 0) {
        com->fprintf(stderr, "*** Error: Maximum ROM dimension must be specified (and greater than zero)\n");
        exit (-1);
      } else {
        maxDimension = ioData->romOffline.rob.dataCompression.maxBasisSize;
        minDimension = ioData->romOffline.rob.dataCompression.minBasisSize;
        energyTol = ioData->romOffline.rob.dataCompression.maxEnergyRetained;
      }*/
    }
  } else { // TODO ONLINE
    if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      if (ioData->romOnline.maxDimension <= 0) {
        com->fprintf(stderr, "*** Error: Maximum ROM dimension must be specified (and greater than zero)\n");
        exit (-1);
      } else {
        maxDimension = ioData->romOnline.maxDimension;
        minDimension = ioData->romOnline.minDimension;
        energyTol = ioData->romOnline.energy;
      }
    } else if (strcmp(basisType,"sensitivity")==0) {
      if (ioData->romOnline.sensitivity.maxDimension <= 0) {
        com->fprintf(stderr, "*** Error: Maximum ROM dimension must be specified (and greater than zero)\n");
        exit (-1);
      } else {
        maxDimension = ioData->romOnline.sensitivity.maxDimension;
        minDimension = ioData->romOnline.sensitivity.minDimension;
        energyTol = ioData->romOnline.sensitivity.energy;
      }
    } else if (strcmp(basisType,"krylov")==0) {
      if (ioData->romOnline.krylov.maxDimension <= 0) {
        com->fprintf(stderr, "*** Error: Maximum ROM dimension must be specified (and greater than zero)\n");
        exit (-1);
      } else {
        maxDimension = ioData->romOnline.krylov.maxDimension;
        minDimension = ioData->romOnline.krylov.minDimension;
        energyTol = ioData->romOnline.krylov.energy;
      }
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
  } else if (strcmp(basisType,"sensitivity")==0) {
      iCluster = -2;
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sensitivityBasisName, iCluster, basisPath);
  } else {
      exit (-1);
  }

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
  sVals = new double[maxDimension];

  int basisSize = maxDimension;
  int moreSVals = 0;
  for (int iVec=0; iVec<maxDimension; ++iVec) {
   moreSVals = fscanf(singValFile,"%d %le %le", &vecNumber, &(sVals[iVec]), &percentEnergy);
   if ((moreSVals!=0) && (percentEnergy>=energyTol)&&((iVec+1)>=minDimension)) {
     basisSize = (iVec+1);
     com->fprintf(stdout, " ... using %d vectors from basis %d (energy tolerance = %e)\n", basisSize, iCluster, energyTol);
     double* sValsNew = new double[basisSize];
     for (int iSVal = 0; iSVal<basisSize; ++iSVal)
       sValsNew[iSVal] = sVals[iSVal];
     delete [] sVals;
     sVals = sValsNew;
     sValsNew = NULL;
     break;
   }
  }

  fclose(singValFile);

  // read in basis vectors

  com->fprintf(stdout, "\nReading basis %d\n", iCluster);

  basis = new VecSet< DistSVec<double, dim> >(basisSize, domain.getNodeDistInfo());
  bool moreVecs;
  int numBasisVecs = 0;
  double tmpSVal;

  for (int iVec = 0; iVec<basisSize; ++iVec) {
    moreVecs = domain.readVectorFromFile(basisPath, iVec, &tmpSVal, (*basis)[iVec]);
    if (!moreVecs) break;
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

  // also read in the basis update quantities
  if (strcmp(basisType, "state") == 0) { 

    char *basisUpdatePath = NULL;
    determinePath(updateInfoName, iCluster, basisUpdatePath);

    com->fprintf(stdout, "\nReading update info for basis %d \n", iCluster);

    FILE *basisUpdateFile = fopen(basisUpdatePath, "r");

    if (!basisUpdateFile)  {
      com->fprintf(stderr, "*** Error: unable to open file %s\n", basisUpdatePath);
      exit (-1);
    }

    columnSumsV = new double[snapsInCluster[iCluster]];
    for (int iVec = 0; iVec < snapsInCluster[iCluster]; ++iVec){
      _n = fscanf(basisUpdateFile, "%le", &(columnSumsV[iVec]));
    }

    fclose(basisUpdateFile);

    delete [] basisUpdatePath;
    basisUpdatePath = NULL;
  }

  com->barrier();

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusterBasis(int iCluster, int nTotSnaps, char* basisType) {

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
  for (int iSnap = 0; iSnap<nTotSnaps; ++iSnap) sValSqrdSum += (sVals[iSnap] * sVals[iSnap]);
  double invSValSqrdSum = 1/sValSqrdSum;

  double sValSqrdPartialSum = 0;

  if (com->cpuNum() == 0) {
    FILE *singValFile = fopen(singValsPath, "wt");
    //com->fprintf(singValFile,"Vector# SVal SValSum\n");

    for (int iVec=0; iVec<nTotSnaps; ++iVec) {
      sValSqrdPartialSum += (sVals[iVec]*sVals[iVec]*invSValSqrdSum);
      com->fprintf(singValFile,"%d %1.12e %1.12e\n", iVec, sVals[iVec], sValSqrdPartialSum);
     }

    fclose (singValFile);
  }

  com->barrier();
 
  delete [] singValsPath;
  singValsPath = NULL;
  delete [] sVals;
  sVals = NULL;

  if (strcmp(basisType, "state") == 0) { 
    // write thinSVD update info to ASCII file
    com->fprintf(stdout, "\nWriting thin SVD update info to disk\n");

    char *basisUpdatePath = 0;
    determinePath(updateInfoName, iCluster, basisUpdatePath);

    if (com->cpuNum() == 0) {
      FILE *basisUpdateFile = fopen(basisUpdatePath, "wt");

      for (int iVec=0; iVec<snapsInCluster[iCluster]; ++iVec) {
        com->fprintf(basisUpdateFile,"%le\n", columnSumsV[iVec]);
      }

      fclose (basisUpdateFile);
    }

    delete [] basisUpdatePath;
    basisUpdatePath = NULL;
    delete [] columnSumsV;
    columnSumsV = NULL;
  }

  com->barrier();

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
  bool endOfFile;
  char* fullRefName;
  const char* refSolName = ioData->input.stateSnapRefSolution;
  fullRefName = new char[strlen(ioData->input.prefix) + strlen(refSolName) + 1];
  sprintf(fullRefName, "%s%s", ioData->input.prefix, refSolName);
  com->fprintf(stdout, "\nReading reference solution for snapshots from %s\n", fullRefName);
  snapRefState = new DistSVec<double, dim>(domain.getNodeDistInfo());
  endOfFile = domain.readVectorFromFile(fullRefName, 0, &tmp, *snapRefState);
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

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusterSampleNodes(int iCluster) {

  char *sampleNodesPath = 0;
  determinePath(sampledNodesName, iCluster, sampleNodesPath);
  FILE *sampleNodeFile = fopen(sampleNodesPath, "r");
  delete [] sampleNodesPath;

  int count;
  count = fscanf(sampleNodeFile, "%d", &nSampleNodes);  // first entry is the number of sample nodes

  sampleNodes.reserve(nSampleNodes);  // know it will be nSampleNodes long (efficiency)

  int index, currentSampleNode;
  for (int i = 0; i < nSampleNodes; ++i){
    count = fscanf(sampleNodeFile, "%d", &index);
    count = fscanf(sampleNodeFile, "%d", &currentSampleNode);
    sampleNodes.push_back(currentSampleNode-1); // reads in the sample node plus one (TODO - change this after code validation)
  }

  fclose(sampleNodeFile);

  restrictionMapping.reset(new RestrictionMapping<dim>(&domain, sampleNodes.begin(), sampleNodes.end()));

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusterGappyMatrix(int iCluster, char* matrixType) {

  char* matrixPath = 0;

  if (strcmp(matrixType,"resMatrix")==0) {
      if (resMat) {
        delete resMat;
        resMat = NULL;
      }
      determinePath(gappyResidualName, iCluster, matrixPath);
  } else if (strcmp(matrixType,"jacMatrix")==0) {
      if (jacMat) {
        delete jacMat;
        jacMat = NULL;
      }
      determinePath(gappyJacActionName, iCluster, matrixPath);
  } else {
      exit (-1);
  }

  double tag = 0.0;
  int numVecs = 0;
  int step = 0;
  bool status = domain.readTagFromFile<double, dim>(matrixPath, step, &tag, &numVecs);  // if file DNE, returns false, tag=0, and numSteps=0

  VecSet<DistSVec<double, dim> > gappyMatrix(numVecs, domain.getNodeDistInfo());

  for (int iVec=0; iVec<numVecs; ++iVec) {
    status = domain.readVectorFromFile(matrixPath, iVec, &tag, gappyMatrix[iVec]);
  }

  delete [] matrixPath;

  // restrictionMapping is set during readClusterSampleNodes
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
void NonlinearRom<dim>::writeReducedCoords(const int totalTimeSteps, bool clusterSwitch, bool updateFreq, int iCluster, Vec<double> dUrom) {

  int nPod = basis->numVectors();

  if (this->com->cpuNum() == 0) {
    if (clustUsageFile)  // for plotting only
      this->com->fprintf(clustUsageFile,"%d %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
    if (reducedCoordsFile) {
      if (clusterSwitch) {
        if (this->ioData->romOnline.basisUpdates) {
          this->com->fprintf(reducedCoordsFile,"%d switch update %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        } else {
          this->com->fprintf(reducedCoordsFile,"%d switch noUpdate %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        }
      } else if (updateFreq && this->ioData->romOnline.basisUpdates) {
        this->com->fprintf(reducedCoordsFile,"%d noSwitch update %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      } else { 
        this->com->fprintf(reducedCoordsFile,"%d noSwitch noUpdate %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      }
      for (int iPod=0; iPod<nPod; ++iPod) {
        this->com->fprintf(reducedCoordsFile, "%23.15e\n", dUrom[iPod]);
      }
    }
  }

  //this->com->barrier();

}



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

