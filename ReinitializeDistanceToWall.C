#include <ReinitializeDistanceToWall.h>
#include <LevelSet/LevelSetStructure.h>
#include <Domain.h>
#include <DistVector.h>

//------------------------------------------------------------------------------

template<int dimLS>
ReinitializeDistanceToWall<dimLS>::ReinitializeDistanceToWall(IoData &ioData, Domain& domain)
  : iod(ioData),dom(domain),done(domain.getNodeDistInfo()),d2wall(domain.getNodeDistInfo()),tag(domain.getNodeDistInfo()),dummyPhi(domain.getNodeDistInfo()),sortedNodes(domain.getNodeDistInfo())
{
  int nSub         = dom.getNumLocSub();
  nSortedNodes     = new int[nSub];
  firstCheckedNode = new int[nSub];
}

//------------------------------------------------------------------------------

template<int dimLS>
ReinitializeDistanceToWall<dimLS>::~ReinitializeDistanceToWall()
{
  delete[] nSortedNodes;
  delete[] firstCheckedNode;
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::ComputeWallFunction(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{

  if (iod.eqs.tc.tm.d2wall.type ==  WallDistanceMethodData::ITERATIVE) {
    DistanceToClosestPointOnMovingStructure(LSS,X,distGeoState);
    GetLevelsFromInterfaceAndMarchForward(LSS,X,distGeoState);
  }
  else if (iod.eqs.tc.tm.d2wall.type ==  WallDistanceMethodData::NONITERATIVE) {
    PseudoFastMarchingMethod(LSS,X,distGeoState);
  }
  else {
  // PrescribedValues(LSS,X,distGeoState);
    fprintf(stderr," *** Error ***, Unknown wall distance method\n");
    exit(1);
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
      for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
          distGeoState(iSub).getDistanceToWall()[i]=d2wall(iSub)[i][0];
  }
  
#ifdef ERROR_CHECK
  computeExactErrors(LSS,X,distGeoState);
#endif
  return;
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::DistanceToClosestPointOnMovingStructure(
	DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{
  done=false;
  tag=0;
  DistVec<ClosestPoint> *closestPoint = LSS.getClosestPointsPointer();

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {

// Fill with initial guess
#if 1
    d2wall=1e10;
#else
    for (int i = 0; i < distGeoState(iSub).getDistanceToWall().size(); ++i)
      d2wall(iSub)[i][0] = distGeoState(iSub).getDistanceToWall()[i];
#endif

    InitializeWallFunction(*dom.getSubDomain()[iSub],LSS(iSub),done(iSub),X(iSub),d2wall(iSub),tag(iSub),closestPoint?&((*closestPoint)(iSub)):NULL);
    dom.getSubDomain()[iSub]->sndData(*dom.getVolPat(),d2wall(iSub).data());
  }

  dom.getVolPat()->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < dom.getNumLocSub(); ++iSub) {
    dom.getSubDomain()[iSub]->minRcvData(*dom.getVolPat(), d2wall(iSub).data());
  }
  return;
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::InitializeWallFunction(SubDomain& subD,LevelSetStructure& LSS,Vec<bool>& done,SVec<double,3>& X,SVec<double,1>& d2w,Vec<int>& tag,Vec<ClosestPoint> *closestPoint)
{
  int (*ptrEdge)[2]=subD.getEdges().getPtr();
  for(int l=0;l<subD.getEdges().size();++l){
    if(LSS.edgeIntersectsStructure(0,l)){
      int i=ptrEdge[l][0],j=ptrEdge[l][1];
      done[i]=true;tag[i]=1;
      done[j]=true;tag[j]=1;
      d2w[i][0] = closestPoint?(*closestPoint)[i].dist:0.0;
      d2w[j][0] = closestPoint?(*closestPoint)[j].dist:0.0;
    }
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::GetLevelsFromInterfaceAndMarchForward(
	DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{

  int min_level=0,level=1,max_level=1;
  dummyPhi = 1.0;
  while(min_level<=0){ // Tag every level
    dom.TagInterfaceNodes(0,tag,dummyPhi,level,&LSS);
    min_level=1;
    for(int iSub = 0; iSub < dom.getNumLocSub(); ++iSub){// TODO(jontg): Parallelize better!
      for(int i = 0; i < done(iSub).len ; ++i) {
	  min_level=min(min_level,tag(iSub)[i]);
          max_level=max(max_level,tag(iSub)[i]);
      }
    }
    dom.getCommunicator()->globalMin(1,&min_level);
    ++level;
  }
  dom.getCommunicator()->globalMax(1,&max_level);

// Propagate information outwards
//
  MultiFluidData::CopyCloseNodes copy=MultiFluidData::FALSE;
  bool printwarning = false;
  double maxres = -FLT_MAX; 
  int maxreslvl = 1;
  for(int ilvl=2;ilvl<=max_level;++ilvl){
    double res = 1.0, resn = 1.0, resnm1 = 1.0;
    int it = 0;
    while(res>iod.eqs.tc.tm.d2wall.eps && it<iod.eqs.tc.tm.d2wall.maxIts){
      resnm1 = resn;
      dom.computeDistanceLevelNodes(1,tag,ilvl,X,d2wall,resn,dummyPhi,copy);
      dom.getCommunicator()->globalMax(1,&resn);
      it++;
      res = fabs((resn-resnm1)/(resn+resnm1));
      // dom.getCommunicator()->fprintf(stderr,"Level %d: Iter %d Res %f\n",ilvl,it,res);
    }
    if (res>iod.eqs.tc.tm.d2wall.eps) {
      printwarning = true;
      if (res > maxres) {
	maxres = res;
	maxreslvl = ilvl;
      }
    }
  }
  if (printwarning) dom.getCommunicator()->fprintf(stderr, "*** Warning: Distance to wall computation (Max residual: %e at level: %d, target: %e, total levels: %d)\n", maxres,maxreslvl,iod.eqs.tc.tm.d2wall.eps,max_level);

  return;
}

//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::PseudoFastMarchingMethod(
	DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{
  // The following is an adaptation of the Fast Marching Method to Embedded Turbulent computation. 
  // Adam 2012.09
  sortedNodes  =-1;
  int nSub     = dom.getNumLocSub();
  d2wall       = 1.0e10;
  tag          = -1;
  int isDone   = 0;

  // First, get an accurate distance to the wall for the first layer
  /*
  DistVec<ClosestPoint>& closestPoint(LSS.getClosestPoints());
  #pragma omp parallel for
  for(int iSub=0;iSub<nSub;iSub++) {
    nSortedNodes[iSub]     = 0;
    firstCheckedNode[iSub] = 0;
    for(int i=0;i<closestPoint(iSub).size();++i) {
      if(!LSS(iSub).isActive(0.0,i))
      {
	d2wall(iSub)[i][0] = 0.0;
	tag(iSub)[i]       = 0;
	sortedNodes(iSub)[nSortedNodes[iSub]] = i;
	nSortedNodes[iSub]++;
      }
    }
    dom.getSubDomain()[iSub]->sndData(*dom.getLevelPat(),reinterpret_cast<int (*)[1]>(tag.subData(iSub)));
    dom.getSubDomain()[iSub]->sndData(*dom.getVolPat(),d2wall(iSub).data());
  }
  #pragma omp parallel for
  for(int iSub=0;iSub<nSub;iSub++) {
    dom.getSubDomain()[iSub]->maxRcvDataAndCountUpdates(*dom.getLevelPat(),reinterpret_cast<int (*)[1]>(tag.subData(iSub)),nSortedNodes[iSub],sortedNodes(iSub));
    dom.getSubDomain()[iSub]->minRcvData(*dom.getVolPat(), d2wall(iSub).data());
    firstCheckedNode[iSub] = nSortedNodes[iSub];
  }
 

  #pragma omp parallel for
  for(int iSub=0;iSub<nSub;iSub++) {
    for(int i=0;i<closestPoint(iSub).size();++i) {
      if(tag(iSub)[i]<0) {
        if(closestPoint(iSub)[i].nearInterface() && closestPoint(iSub)[i].dist >= 0)
        {
  	  d2wall(iSub)[i][0] = closestPoint(iSub)[i].dist;
	  tag(iSub)[i]       = 1;
	  sortedNodes(iSub)[nSortedNodes[iSub]] = i;
	  nSortedNodes++;
	}
      }
    }
    dom.getSubDomain()[iSub]->sndData(*dom.getLevelPat(),reinterpret_cast<int (*)[1]>(tag.subData(iSub)));
    dom.getSubDomain()[iSub]->sndData(*dom.getVolPat(),d2wall(iSub).data());
  }
  #pragma omp parallel for
  for(int iSub=0;iSub<nSub;iSub++) {
    dom.getSubDomain()[iSub]->maxRcvDataAndCountUpdates(*dom.getLevelPat(),reinterpret_cast<int (*)[1]>(tag.subData(iSub)),nSortedNodes[iSub],sortedNodes(iSub));
    dom.getSubDomain()[iSub]->minRcvData(*dom.getVolPat(), d2wall(iSub).data());
  }
  */
  DistVec<ClosestPoint> *closestPoint = LSS.getClosestPointsPointer();
  int level   =  0; // Level 0 (inActive nodes) and 1 (Embedded surface neighbors) 
  while(isDone == 0){ // Tag every level
    dom.pseudoFastMarchingMethod<1>(tag,X,d2wall,level,sortedNodes,nSortedNodes,firstCheckedNode,&LSS,closestPoint);
    // I don't think it is a good idea to OMP parallelize this loop. nSub should be small, though!
    isDone = 1;
    for(int iSub = 0; iSub < nSub; ++iSub){
      if(nSortedNodes[iSub] != tag(iSub).len) {isDone=0; break;}
    }
    dom.getCommunicator()->globalMin(1,&isDone);
    ++level;
  }
  dom.getCommunicator()->globalMax(1,&level);
//  dom.getCommunicator()->fprintf(stderr,"There are %d levels\n",--level);
  return;
}
//------------------------------------------------------------------------------
//
template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::computeExactErrors(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{
  double localError = 0.0;
  int nSub = dom.getNumLocSub();
  double **errors;
  int      nDofs[nSub];
  bool   **masterFlag;
  errors = new double*[nSub];
  masterFlag = new bool*[nSub];

#pragma omp parallel for 
  for (int iSub=0;iSub<nSub;++iSub) {
    errors[iSub]     = new double[3];
    errors[iSub][0] = 0.0; 
    errors[iSub][1] = 0.0;
    errors[iSub][2] = 0.0;
    masterFlag[iSub] = d2wall.info().getMasterFlag(iSub);
    nDofs[iSub]      = 0; 
    for (int i=0; i<distGeoState(iSub).getDistanceToWall().size();++i) {
      if(LSS(iSub).isActive(0.0,i) && masterFlag[iSub][i]) {
	nDofs[iSub]++;
        localError = fabs(
		     d2wall(iSub)[i][0] 
		   - sqrt(X(iSub)[i][0]*X(iSub)[i][0]+
			  X(iSub)[i][1]*X(iSub)[i][1]+
			  X(iSub)[i][2]*X(iSub)[i][2]) 
		   + 1.0);
	errors[iSub][0] += localError;
	errors[iSub][1] += localError*localError;
        errors[iSub][2] = errors[iSub][2]<localError?localError:errors[iSub][2];
      }
    }
  }

  for(int iSub=1;iSub<nSub;++iSub) {
    errors[0][0] += errors[iSub][0]; 
    errors[0][1] += errors[iSub][1]; 
    errors[0][2]  = errors[iSub][2]>errors[0][2]?errors[iSub][2]:errors[0][2];
    nDofs [0]    += nDofs[iSub];
  }

  dom.getCommunicator()->globalSum(1,nDofs);
  dom.getCommunicator()->globalSum(2,errors[0]);
  dom.getCommunicator()->globalMax(1,errors[0]+2);

  errors[0][0] /= nDofs[0];
  errors[0][1] /= nDofs[0];
  errors[0][1] = sqrt(errors[0][1]);

  dom.getCommunicator()->fprintf(stderr,"Distance to the wall computation Error for the Embedded Cylinder\n d2wall Error: %12.8e %12.8e %12.8e\n",errors[0][0],errors[0][1],errors[0][2]);
  for(int iSub=0;iSub<nSub;iSub++) delete[] errors[iSub];
  delete[] errors;
  return; 
}
//------------------------------------------------------------------------------

template<int dimLS>
void ReinitializeDistanceToWall<dimLS>::PrescribedValues(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState)
{
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

//------------------------------------------------------------------------------
template class ReinitializeDistanceToWall<1>;
template class ReinitializeDistanceToWall<2>;
