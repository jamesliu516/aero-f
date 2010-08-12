#include <LevelSet/LevelSetStructure.h>
#include <DistVector.h>
#include <Vector.h>

DistVec<double> &DistLevelSetStructure::getPhi() {
  for (int iSub=0; iSub<numLocSub; iSub++)
    (*this)(iSub).computePhi((*pseudoPhi)(iSub));
  return *pseudoPhi;
}

void LevelSetStructure::computePhi(Vec<double> &phi) {
      for(int i = 0; i < phi.size(); ++i)
        phi[i] = isActive(0,i) ? 1 : -1;
    }
