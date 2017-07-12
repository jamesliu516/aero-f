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

template <int dimLS, int dim>
class ReinitializeDistanceToWall
{
  IoData &iod;
  Domain &dom;
  DistSVec<double, 1> d2wall;

  // default wall distance
  DistVec<int> tag;
  DistVec<int> sortedNodes;
  int *nActiveNodes, *firstCheckedNode; //,*nSortedNodes
  // DistSVec<double, dimLS> dummyPhi;

  // FEM wall distance
  int **activeElemList, **elemTag, **knownNodes;
  int *nSortedElems, *nSortedNodes, *firstCheckedElem, *nUnsortedNodes;
  DistVec<int> isSharedNode, nodeTag, unsortedNodes, unsortedTag;

  // predictors
  double predictorTime[3];
  DistSVec<double, 1> d2wnm1, d2wnm2;
  int countreinits, nPredTot;

  double tolmin, tolmax;
  // double cw1;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, DistSVec<double, dim> &V, double t);

private:
  // sjg, 04/2017: wall distance predictor
  int UpdatePredictorsCheckTol(DistLevelSetStructure &LSS, DistGeoState &distGeoState, DistSVec<double, dim> &V, double t);
  void ReinitializePredictors(int update, double t, DistLevelSetStructure *LSS, DistSVec<double, 3> &X);

  // FEM wall distance
  void PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);

  // default wall distance
  void PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, int iterativeLevel);
  void IterativeMethodUpdate(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);
  // void InitializeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, double t);
  // void GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

  void ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState);
  void PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
};

#endif
