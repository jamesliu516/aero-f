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
  DistSVec<double, 1> d2wall;
  DistVec<int> tag;
  DistVec<int> sortedNodes;
  int *nSortedNodes, *nActiveNodes, *firstCheckedNode;
  // DistSVec<double, dimLS> dummyPhi;

  // predictors
  double *predictorTime;
  DistSVec<double, 1> d2wnm1, d2wnm2;
  int countreinits, nPredTot;

  double tolmin, tolmax;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, double t);

private:
  // sjg, 04/2017: wall distance predictor
  int UpdatePredictorsCheckTol(double t, DistLevelSetStructure &LSS, DistGeoState &distGeoState);
  void ReinitializePredictors(int update, double t, DistLevelSetStructure *LSS, DistSVec<double, 3> &X);

  int PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, int iterativeLevel);
  void IterativeMethodUpdate(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);
  // void InitializeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, double t);
  // void GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

  void ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState);
  void PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
};

#endif
