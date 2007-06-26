#include <LevelSet.h>

#include <DistVector.h>
#include <Vector.h>


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <alloca.h>








template<int dim>
void LevelSet::conservativeToPrimitive(DistSVec<double,dim> &U)
{
  int numLocSub = Phi.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*phi) = Phi.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      phi[i] /= u[i][0];
  }
                                                                                                                
}
                                                                                                                
//-------------------------------------------------------------------------
template<int dim>
void LevelSet::primitiveToConservative(DistSVec<double,dim> &U)
{
  int numLocSub = Phi.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*u)[dim] = U.subData(iSub);
    double (*phi) = Phi.subData(iSub);
    for (int i=0; i<U.subSize(iSub); i++)
      phi[i] *= u[i][0];
  }
                                                                                                                
}



