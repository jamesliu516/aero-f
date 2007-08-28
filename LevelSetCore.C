//
//       Class LevelSet   created by Francois Courty (C.A.S. Boulder)
//       March 2004
//
///////////////////////////////////////////////////////////////////////
#include <LevelSet.h>



#include <SubDomain.h>

#include <Vector3D.h>
#include <Vector.h>
#include <DistVector.h>
#include <Communicator.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <alloca.h>





LevelSet::LevelSet(DistSVec<double,3> & X, IoData &iod, Domain *domain) : Phi(domain->getNodeDistInfo()),  PhiRet(domain->getNodeDistInfo()), Phi1(domain->getNodeDistInfo()), Phi2(domain->getNodeDistInfo())

{ 
 int numLocSub = Phi.numLocSub();
 com       = domain->getCommunicator();
 double cenx, ceny, cenz;

 nspheres = 0;
 spheres[nspheres][0] = iod.mf.icd.s1.cen_x;
 spheres[nspheres][1] = iod.mf.icd.s1.cen_y;
 spheres[nspheres][2] = iod.mf.icd.s1.cen_z;
 spheres[nspheres][3] = iod.mf.icd.s1.r;
 ++nspheres;

 double dist[10];
 int j = 0;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*x)[3] = X.subData(iSub);
    double (*phi) = Phi.subData(iSub);
    for (int i=0; i<X.subSize(iSub); i++){
      //for bubble
      dist[j] = - spheres[j][3] + sqrt((x[i][0]-spheres[j][0])*(x[i][0]-spheres[j][0])+
                                         (x[i][1]-spheres[j][1])*(x[i][1]-spheres[j][1])+
                                         (x[i][2]-spheres[j][2])*(x[i][2]-spheres[j][2]));
      //for shock tube (Charbel's regular quasi-1D mesh with x from 0 to 1)
      //cf DistTimeState.C
      //dist[j]  = x[i][0]  -0.5  +0.001;

      phi[i]  = dist[j];
//      com->printf(2, "phi init = %f\n", phi[i]);
    }
 }

 Phi1 = Phi;
 Phi2 = Phi;
 com->printf(2,"Norm, min, max of Phi = %6.4f %6.4f %6.4f\n",Phi.norm(),Phi.min(), Phi.max());
 com->printf(2,"Center of bubble = %6.4f %6.4f %6.4f\n",iod.mf.icd.s1.cen_x, iod.mf.icd.s1.cen_y, iod.mf.icd.s1.cen_z);
}
//-------------------------------------------------------------------------
LevelSet::~LevelSet()
{
}
//---------------------------------------------------------------------------------------------------------
DistVec<double>& LevelSet::getPhi()
{
    return Phi;
}



