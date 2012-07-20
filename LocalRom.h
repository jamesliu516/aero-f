#ifndef _LOCAL_ROM_H_
#define _LOCAL_ROM_H_

template <int dim>
class LocalRom {

  private:
  Domain &domain;
  Communicator *com; 
  IoData *ioData;

  // read from ioData->snapshots.clustering
  int nClusters;
  int minClusterSize;
  int kMeansRandSeed;  
  double percentOverlap;
  double kMeansTol;

  // read from ioData->snapshots.dataCompression
  int maxBasisSize;
  double singValTolerance;
  double maxPercentEnergyRetained;

  // file and directory names
  const char* prefix;
  const char* indexName;
  const char* mapName;
  const char* connName;
  const char* centersName;
  const char* nearestName;
  const char* snapsName;
  const char* basisName;
  const char* clusterName;
  const char* databaseName;
  const char* singValName;
  const char* projErrorName;
  const char* basisUpdateName;

  // intermediate quantities
  int nTotSnaps;
  VecSet< DistSVec<double, dim> >* snap; // snap(nTotSnaps, domain.getNodeDistInfo())
  VecSet< DistSVec<double, dim> >* clusterCenters; // average of all snapshots in a cluster 
  VecSet< DistSVec<double, dim> >* nearestSnapsToCenters; // closest snapshot to each cluster center 
  int* snapsInCluster; // number of snaps in each cluster
  int* clusterIndex; // stores original cluster association for each snapshot (before any overlap is introduced)
  int** clusterSnapshotMap;  // one vector per cluster, lists snapshots to include in the cluster (including overlapping snapshots) 
  int** clusterNeighbors;  // one vector per cluster, lists neighboring clusters
  int* clusterNeighborsCount; // stores number of neighbors for each cluster
  VecSet<Vec<double> >* projErrorLog;
  DistSVec<double, dim>* snapRefSol; 

  // thin svd update quantities
  double* columnSumsV;
  double* sVals;
  DistSVec<double, dim>* Uinit; 

  // private functions
  void readSnapshotFile(bool);
  double distance(DistSVec<double, dim> &, DistSVec<double, dim> &);
  double calcResidual(VecSet< DistSVec<double, dim> > &, VecSet< DistSVec<double, dim> > &);
  void localPod();
  void kmeans();
  void outputClusterSnapshots();
  void outputClusterBasis(int);
  void readClusterSnapshots(int, bool);
  void readNearestSnapsToCenters();
  void localRelProjError();
  void writeProjErrorToDisk();
  void readReferenceSolution();
  void outputClusterReferenceSnapshot(int, DistSVec<double, dim> &);
  void readClusterReferenceSnapshot(int);
 
 
  public:

  LocalRom(Communicator *, IoData &, Domain &);
  ~LocalRom();

  VecSet< DistSVec<double, dim> >* basis;

  void createAndAnalyzeDatabase();
  void closestCenter(DistSVec<double, dim> &, int &, int &, double &, double &);
  void readClusterCenters();
  void readClusterBasis(int);
  void updateBasis(int, DistSVec<double, dim> &);

};
#include "LocalRom.C"
#endif
