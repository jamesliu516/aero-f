#ifndef _NONLINEAR_ROM_DATABASE_CONSTRUCTION_H_
#define _NONLINEAR_ROM_DATABASE_CONSTRUCTION_H_

#include <NonlinearRomOnlineII.h>

template <int dim>
class NonlinearRomDatabaseConstruction : public NonlinearRomOnlineII<dim> {

  protected:

  ROBConstructionData* robConstruction;
  RelativeProjectionErrorData* projError;

  GeoSource &geoSource;

  VecSet<Vec<double> >* projErrorLog;
//  DistSVec<double, dim>* snapRefSol; 

  // for snapshot collection method 0 (requires offline clustering of FOM residuals and Krylov vectors)
  int nSnapshotFiles;  // number of snapshot files; should be the same for all FOM snapshots (state, residual, etc.)
  std::vector<std::vector<std::vector<bool> > > stateSnapshotClustersAfterOverlap;  // [iFile][iSnap][iCluster] = true or false
                                              // depending on whether snap iSnap from file iFile is a member of cluster iCluster

  // private functions
  double calcResidual(VecSet< DistSVec<double, dim> > &, VecSet< DistSVec<double, dim> > &);
  void localPod(const char*);
  void SVD(VecSet< DistSVec<double, dim> >*&, VecSet< DistSVec<double, dim> > &, double *S, FullM &, int, int, bool computeV=true);
  void scalapackSVD(VecSet< DistSVec<double, dim> >*&, VecSet< DistSVec<double, dim> > &, double *S, FullM &, int, bool computeV=true);
  void probabilisticSVD(VecSet< DistSVec<double, dim> >*&, VecSet< DistSVec<double, dim> > &, double *S, FullM &V, int, bool computeV=true);
  void initializeClusterCenters();
  void kmeans();
  void kmeansWithBounds();
  void computeClassicalMultiDimensionalScaling();
  void localRelProjError(); 
  void localRelProjErrorSweep();

  // fast distance calculation preprocessing
  bool arbitraryUniformIC;
  DistSVec<double, dim>* initialCondition;
  void readInitialCondition();
  void preprocessForDistanceComparisons();
  void productOfBasisAndCenterDifferences(int, const char*);
  void productOfVectorAndCenterDifferences(int, const char*);

  // gnat basis updates preprocessing
  void preprocessForExactBasisUpdates();

  // IO functions that are independent of database structure
  void writeProjErrorToDisk();
  void writeProjErrorSweepToDisk(std::vector<std::vector<int> >, std::vector<std::vector<double> >);
  void placeNonStateSnapshotsInClusters(const char*);

  public:

  NonlinearRomDatabaseConstruction(Communicator *, IoData &, Domain &, GeoSource &);
  ~NonlinearRomDatabaseConstruction();

  void constructDatabase();

};

#include "NonlinearRomDatabaseConstruction.C"
#endif
