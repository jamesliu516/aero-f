#include <LevelSet/LevelSetStructure.h>
#include <LevelSet/IntersectionFactory.h>
#include <Communicator.h>
#include <DistVector.h>
#include <Vector.h>

void DistLevelSetStructure::clearTotalForce()
{
  for (int iSub=0; iSub<numLocSub; iSub++)
    (this)[iSub].totalForce[0] = (this)[iSub].totalForce[1] = (this)[iSub].totalForce[2] = 0.0;
  totalForce[0] = totalForce[1] = totalForce[2] = 0.0;
}

//-----------------------------------------------------------------------------------

Vec3D DistLevelSetStructure::getTotalForce()
{
  for (int iSub=0; iSub<numLocSub; iSub++) {
    totalForce[0] += (this)[iSub].totalForce[0];
    totalForce[1] += (this)[iSub].totalForce[1];
    totalForce[2] += (this)[iSub].totalForce[2];
  }
  Communicator *com = IntersectionFactory::getCommunicator();
  com->globalSum(3, totalForce);

 // com->fprintf(stderr,"Total Force on Structure Surface = [%f, %f, %f]\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
 // com->fprintf(forceFile, "%lf %lf %lf\n", totalForce[0]*pref, totalForce[1]*pref, totalForce[2]*pref);
  return totalForce;
}

DistVec<double> &DistLevelSetStructure::getPhi() {
  for (int iSub=0; iSub<numLocSub; iSub++)
    (*this)(iSub).computePhi((*pseudoPhi)(iSub));
  return *pseudoPhi;
}

void LevelSetStructure::computePhi(Vec<double> &phi) {
      for(int i = 0; i < phi.size(); ++i)
        phi[i] = isActive(0,i) ? 1 : -1;
    }
