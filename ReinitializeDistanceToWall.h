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
  DistVec<int> tag, sortedNodes, isSharedNode;
  int *firstCheckedNode, *nSortedNodes;

  // for distance error estimation
  SpaceOperator<dim> &spaceOp;
  DistVec<int> predictorTag;
  double predictorTime[3];
  DistSVec<double, 1> d2wnm1, d2wnm2;
  DistVec<double> SAsensitiv;
  double SAsensitivScale;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain, SpaceOperator<dim> &spaceOp);
  ~ReinitializeDistanceToWall();

  int ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, const double t);

private:
  // wall distance algorithm methods
  void PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);
  void IterativeMethodUpdate(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);

  // wall distance error estimator methods
  int UpdatePredictorsCheckTol(DistLevelSetStructure &LSS, DistGeoState &distGeoState, const double t);
  void ReinitializePredictors(DistGeoState &distGeoState, DistLevelSetStructure *LSS, DistSVec<double, 3> &X);

  // debugging methods
  void ComputeExactErrors(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
  void PrescribedValues(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState);
};

#endif
