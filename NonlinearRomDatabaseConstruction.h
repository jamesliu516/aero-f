#ifndef _NONLINEAR_ROM_DATABASE_CONSTRUCTION_H_
#define _NONLINEAR_ROM_DATABASE_CONSTRUCTION_H_

#include <NonlinearRomOnlineII.h>

template <int dim>
class NonlinearRomDatabaseConstruction : public NonlinearRomOnlineII<dim> {

  protected:

  RelativeProjectionErrorData* ioDataProjError;

  GeoSource &geoSource;

  VecSet<Vec<double> >* projErrorLog;
//  DistSVec<double, dim>* snapRefSol; 

  // for snapshot collection method 0 (requires offline clustering of FOM residuals and Krylov vectors)
  int nSnapshotFiles;         // number of snapshot files; should be the same for all FOM snapshots (state, residual, etc.)
  int* stateSnapsFromFile;    // stateSnapsFromFile[iFile] = number of snapshots taken from file iFile;
  double** stateSnapshotTags; // stateSnapshotInfo[iFile][iSnap] = tag associated with snapshot iSnap from file iFile
  bool*** stateSnapshotClustersAfterOverlap;  // stateSnapshotClustersAfterOverlap[iFile][iSnap][iCluster] = true or false
                                              // depending on whether snap iSnap from file iFile is a member of cluster iCluster

  // private functions
  double calcResidual(VecSet< DistSVec<double, dim> > &, VecSet< DistSVec<double, dim> > &);
  void localPod(char *);
  void kmeans();
  void localRelProjError(); 

  // fast distance calculation preprocessing
  bool arbitraryUniformIC;
  DistSVec<double, dim>* initialCondition;
  void readInitialCondition();

  void preprocessForDistanceComparisons();
  void productOfBasisAndCenterDifferences(int, char*);
  void productOfVectorAndCenterDifferences(int, char*);

  // IO functions that are independent of database structure
  void writeProjErrorToDisk();
  void readSnapshotFile(char *, bool);
  void placeNonStateSnapshotsInClusters(char *);

  public:

  NonlinearRomDatabaseConstruction(Communicator *, IoData &, Domain &, GeoSource &);
  ~NonlinearRomDatabaseConstruction();

  void constructDatabase();

};

#include "NonlinearRomDatabaseConstruction.C"
#endif
