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

template<int dimLS>
class ReinitializeDistanceToWall
{
  Domain& dom;
  DistVec<bool> done;
  DistSVec<double,1> d2wall;

  DistVec<int> tag;
  DistSVec<double,dimLS> dummyPhi;

public:
  ReinitializeDistanceToWall(Domain& domain);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);

private:
  void InitializeWallFunction(SubDomain& subD,LevelSetStructure& LSS,Vec<bool>& done,SVec<double,3>& X,SVec<double,1>& d2w,Vec<int>& tag,Vec<ClosestPoint>& closestPoint);
};

#endif
