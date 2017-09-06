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
template<int dim> class SpaceOperator;

template <int dimLS, int dim>
class ReinitializeDistanceToWall
{
  IoData &iod;
  Domain &dom;
  DistSVec<double, 1> d2wall;

  DistVec<int> tag, sortedNodes, isSharedNode; // acts as an active nodes list with computed distances
  int *firstCheckedNode, *nSortedNodes;
  // DistSVec<double, dimLS> dummyPhi; // depreciated iterative wall distance

  // predictors
  SpaceOperator<dim> &spaceOp;
  DistVec<int> predictorTag;
  double predictorTime[3];
  DistSVec<double, 1> d2wnm1, d2wnm2;
  DistVec<double> SAsensitiv;
  double SAsensitivScale;
  int countReinits, nPredTot;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain, SpaceOperator<dim> &spaceOp);
  ~ReinitializeDistanceToWall();

  int ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, const double t);

private:
  // wall distance algorithms
  void PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);
  void IterativeMethodUpdate(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);

  // wall distance predictor
  int UpdatePredictorsCheckTol(DistLevelSetStructure &LSS, DistGeoState &distGeoState, const double t);
  void ReinitializePredictors(DistGeoState &distGeoState, DistLevelSetStructure *LSS, DistSVec<double, 3> &X);

  // depreciated methods for old iterative implementation
  // void InitializeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, double t);
  // void GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);

  // miscellaneous debugging tools
  void ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void ComputePercentChange(DistLevelSetStructure &LSS, DistGeoState &distGeoState);
  void PrintIntersectedValues(DistLevelSetStructure *LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
};

#endif
