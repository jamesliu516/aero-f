#ifndef __REINITIALIZE_DISTANCE_TO_WALL_H__
#define __REINITIALIZE_DISTANCE_TO_WALL_H__

#include <DistVector.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector.h>

class Domain;
class SubDomain;
class DistGeoState;
class DistLevelSetStructure;
class GeoState;
class LevelSetStructure;

template <int dimLS>
class ReinitializeDistanceToWall
{
  IoData &iod;
  Domain &dom;
  // DistVec<bool> done;           // could delete?
  DistSVec<double, 1> d2wall;
  DistVec<int> sortedNodes;
  int *nSortedNodes, *firstCheckedNode;

  DistVec<int> tag;
  DistSVec<double, dimLS> dummyPhi;

  int predictorActive;
  int *nPredictors;
  double *tPredictors;
  double ***d2wPredictors;
  DistSVec<double, 1> d2wallnm1, d2wallnm2;
  double tnm1, tnm2;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, double t);

private:
  // sjg, 04/2017: wall distance predictor
  int UpdatePredictorsCheckTol(double t);
  void ReinitializePredictors(int update, double t, DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

  void InitializeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, int *nPredLoc, double t);
  void GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  int PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, int iterativeLevel, int *nPredLoc);
  void GetLevelsFromInterfaceAndMarchForward2(int max_level, DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

  void ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState);
  void PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

// private:
  // void InitializeWallFunction(SubDomain &subD, LevelSetStructure &LSS, SVec<double, 3> &X, SVec<double, 1> &d2w, Vec<int> &tag);
}; //Vec<bool> &done,

#endif
