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

  // default wall distance
  DistVec<int> tag, sortedNodes, isSharedNode; // acts as an active nodes list with computed distances
  int *firstCheckedNode, *nSortedNodes;
  // DistSVec<double, dimLS> dummyPhi; // depreciated iterative wall distance

  // FEM wall distance
  int **activeElemList, **elemTag, **knownNodes;
  int *nSortedElems, *firstCheckedElem;
  DistVec<int> nodeTag;

  // predictors
  SpaceOperator<dim> &spaceOp;
  DistVec<int> predictorTag;
  double predictorTime[3];
  DistSVec<double, 1> d2wnm1, d2wnm2;
  int countReinits, nPredTot;
  DistVec<double> SAsensi;

  double tolmin, tolmax; // to be read from parser!

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain &domain, SpaceOperator<dim> &spaceOp);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure &LSS, DistSVec<double, 3> &X, DistGeoState &distGeoState, const double t); //DistSVec<double, dim> &V, double t);

private:
  // wall distance predictor
  int UpdatePredictorsCheckTol(DistLevelSetStructure &LSS, DistGeoState &distGeoState, const double t); // DistSVec<double, dim> &V,
  void ReinitializePredictors(int update, const double t, DistLevelSetStructure *LSS, DistSVec<double, 3> &X);

  // FEM wall distance
  void PseudoFastMarchingMethodFEM(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);

  // default wall distance
  void PseudoFastMarchingMethod(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);
  void IterativeMethodUpdate(DistLevelSetStructure &LSS, DistSVec<double, 3> &X);

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
