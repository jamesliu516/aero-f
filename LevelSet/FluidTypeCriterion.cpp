/*
 * LevelSetCriterion.cpp
 *
 *  Created on: Feb 16, 2009
 *      Author: Michel Lesoinne
 */

#include <LevelSet/FluidTypeCriterion.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector.h>
#include <DistVector.h>

class FluidTypeFromIntersect : public FluidTypeCriterion {
    const LevelSetStructure &intersector;
  public:
    FluidTypeFromIntersect(const LevelSetStructure &p) : intersector(p) {}
    bool isSameFluid(int i, int j) const {
      return !intersector.edgeIntersectsStructure(0, i, j);
    }
    int myPhase(int i) const {
      return intersector.fluidModel(0.0, i);
    }
};

bool
FluidTypeFromLevelSet::isSameFluid(int i, int j) const
{ return phi[i]*phi[j] > 0; }

int
FluidTypeFromLevelSet::myPhase(int i) const 
{ return (phi[i]>0) ? 1 : -1; }

DistFluidTypeFromLevelSet::DistFluidTypeFromLevelSet(const DistVec<double> &phi) {
  ft = new FluidTypeFromLevelSet*[phi.numLocSub()];
  for(int i = 0; i < phi.numLocSub(); ++i)
    ft[i] = new FluidTypeFromLevelSet(phi(i));
}

FluidTypeCriterion &
DistFluidTypeFromLevelSet::operator ()(int i) const {
  return *ft[i];
}

FluidTypeCriterion &
DistFluidTypeFromIntersect::operator()(int i) const {
  std::map<int, FluidTypeCriterion *>::iterator it = ftcMap->find(i);
  if(it != ftcMap->end())
    return *it->second;
  FluidTypeFromIntersect *ftfi = new FluidTypeFromIntersect(dLS(i));
  (*ftcMap)[i] = ftfi;
  return *ftfi;
}
