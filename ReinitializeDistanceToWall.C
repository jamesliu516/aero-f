#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>

//------------------------------------------------------------------------------

template<int dimLS>
ReinitializeDistanceToWall<dimLS>::ReinitializeDistanceToWall(Domain& domain)
  : dom(domain),done(domain.getNodeDistInfo()),d2wall(domain.getNodeDistInfo()),tag(domain.getNodeDistInfo()),dummyPhi(domain.getNodeDistInfo())
{}

//------------------------------------------------------------------------------

template<int dimLS>
ReinitializeDistanceToWall<dimLS>::~ReinitializeDistanceToWall()
{}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::ComputeWallFunction(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{
    if(true){ // Just use prescribed values
      double mind=1e10,maxd=-1e10;
#pragma omp parallel for
      for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
        for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i){
          d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
          mind=min(mind,d2wall(iSub)[i][0]);
          maxd=max(maxd,d2wall(iSub)[i][0]);
        }
      }
      dom.getCommunicator()->globalMin(1,&mind);
      dom.getCommunicator()->globalMax(1,&maxd);
      dom.getCommunicator()->fprintf(stderr,"Min: %e\t\tMax: %e\n",mind,maxd);
      return;
    }

    done=false;
    tag=0;
    DistVec<ClosestPoint>& closestPoint(LSS.getClosestPoints());

    // Fill with initial guess
#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
#if 1
      d2wall=1e10;
#else
      // Just temporary - To be removed
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
            d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
#endif
      // ------------------------------
      InitializeWallFunction(*dom.getSubDomain()[iSub],LSS(iSub),done(iSub),X(iSub),d2wall(iSub),tag(iSub),closestPoint(iSub));
      dom.getSubDomain()[iSub]->sndData(*dom.getVolPat(),d2wall(iSub).data());
    }

    dom.getVolPat()->exchange();

    int min_level=0,level=1,max_level=0;
    while(min_level<=0){ // Tag every level
      dom.TagInterfaceNodes(1,tag,dummyPhi,level);

      min_level=1;
      for(int iSub = 0; iSub < dom.getNumLocSub(); ++iSub){// TODO(jontg): Parallelize better!
        for(int i = 0; i < done(iSub).len ; ++i){min_level=min(min_level,tag(iSub)[i]);
            max_level=max(max_level,tag(iSub)[i]);}
			}
      dom.getCommunicator()->globalMin(1,&min_level);
      ++level;
    }
    dom.getCommunicator()->globalMax(1,&max_level);
    dom.getCommunicator()->fprintf(stderr,"There are %d levels\n",level);

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      dom.getSubDomain()[iSub]->minRcvData(*dom.getVolPat(), d2wall(iSub).data());
    }

// Just temporary - To be removed -------------
#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i) {
        if (tag(iSub)[i]>1 && tag(iSub)[i]<max_level) {d2wall(iSub)[i][0] = 1e10;}
      }
    }
// --------------------------------------------

    // Propagate information outwards
    MultiFluidData::CopyCloseNodes copy=MultiFluidData::FALSE;
    for(int ilvl=2;ilvl<max_level;++ilvl){
    // for(int ilvl=2;ilvl<1;++ilvl){
      double res = 1.0, resn = 1.0, resnm1 = 1.0;
      int it = 0;
      while(res>1e-4 || it<=1){
        resnm1 = resn;
        dom.computeDistanceLevelNodes(1,tag,ilvl,X,d2wall,resn,dummyPhi,copy);
        dom.getCommunicator()->globalMax(1,&resn);
        it++;
        res = fabs((resn-resnm1)/(resn+resnm1));
        // dom.getCommunicator()->fprintf(stderr,"Level %d: Iter %d Res %f\n",ilvl,it,res);
      }
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
        for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
            distGeoState(iSub).getDistanceToWall()[i]=d2wall(iSub)[i][0];}
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::InitializeWallFunction(SubDomain& subD,LevelSetStructure& LSS,Vec<bool>& done,SVec<double,3>& X,SVec<double,1>& d2w,Vec<int>& tag,Vec<ClosestPoint>& closestPoint)
{
#if 1
    for(int i=0;i<closestPoint.size();++i)
        if(closestPoint[i].nearInterface() || closestPoint[i].dist >= 0){
            done[i]=true;tag[i]=1;
            d2w[i][0]=closestPoint[i].dist;}
#else
//    d2w=1e10; // A sufficiently large distance
    int (*ptrEdge)[2]=subD.getEdges().getPtr();
    // Just temporary - To be removed ---------
    for(int l=0;l<subD.getEdges().size();++l){
        if(LSS.edgeIntersectsStructure(0,l)){
            int i=ptrEdge[l][0],j=ptrEdge[l][1];
            done[i]=true;tag[i]=1;
            LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0,i,j);
            d2w[i][0]=LSS.isPointOnSurface(X[i],lsRes.trNodes[0],lsRes.trNodes[1],lsRes.trNodes[2]);

            done[j]=true;tag[j]=1;
            lsRes = LSS.getLevelSetDataAtEdgeCenter(0,j,i);
            d2w[j][0]=LSS.isPointOnSurface(X[j],lsRes.trNodes[0],lsRes.trNodes[1],lsRes.trNodes[2]);
        }
    }
    // ----------------------------------------
    for(int l=0;l<subD.getEdges().size();++l){
        if(LSS.edgeIntersectsStructure(0,l)){
            int i=ptrEdge[l][0],j=ptrEdge[l][1];
            done[i]=true;tag[i]=1;
            LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0,i,j);
            d2w[i][0]=min(d2w[i][0],LSS.isPointOnSurface(X[i],lsRes.trNodes[0],lsRes.trNodes[1],lsRes.trNodes[2]));

            done[j]=true;tag[j]=1;
            lsRes = LSS.getLevelSetDataAtEdgeCenter(0,j,i);
            d2w[j][0]=min(d2w[j][0],LSS.isPointOnSurface(X[j],lsRes.trNodes[0],lsRes.trNodes[1],lsRes.trNodes[2]));
        }
    }

    // MPI OR on done, MIN on d2wall
#endif
}

//------------------------------------------------------------------------------
template class ReinitializeDistanceToWall<1>;
template class ReinitializeDistanceToWall<2>;
