#include "MultiGridLevelSetStructure.h"

MultiGridLevelSetStructure::
MultiGridLevelSetStructure(DistMultiGridLevelSetStructure& lss,
			   SubDomain& sub,
			   Vec<int>& status,Vec<double>& distance,Vec<bool>& is_swept,Vec<bool>& is_active,Vec<bool>& is_occluded,Vec<bool>& edge_intersects,
			   LevelSetStructure* parent, int mySub,MultiGridLevel<double>* myLevel)
  : LevelSetStructure(status, distance, is_swept, is_active, is_occluded, edge_intersects),
    distLSS(lss), parent(parent),
    subD(sub), mySub(mySub), myLevel(myLevel), edges(*myLevel->getEdges()[mySub]),
    status(status), distance(distance), is_swept(is_swept), is_active(is_active),
    is_occluded(is_occluded), edge_intersects(edge_intersects)
    {}

void MultiGridLevelSetStructure::
recompute() {

  status = -1;

  Vec<int>& stat = parent->getStatus();

  int N = stat.size();
  
  Vec<int>& nodeMapping =  myLevel->getNodeMapping()(mySub);
  for (int i = 0; i < N; ++i) {

    status[nodeMapping[i] ] = std::max<int>(status[nodeMapping[i] ] ,
					    stat[i]);
  }
}

void MultiGridLevelSetStructure::
computeEdgeCrossing() {

  int N = edges.size();
  
  int (*ptr)[2] = edges.getPtr();
  for (int i = 0; i < N; ++i) {

    edge_intersects[i] = (status[ptr[i][0]] != status[ptr[i][1]]);
    
  }

  
}

LevelSetResult MultiGridLevelSetStructure::
getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j) {
  if (!edge_intersects[l]) {
    int (*ptr)[2] = edges.getPtr();
    int i=i_less_j ? ptr[l][0] : ptr[l][1],
        j=i_less_j ? ptr[l][1] : ptr[l][0];
    //fprintf(stderr,"%02d There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
    //               globIndex,locToGlobNodeMap[i]+1, status[i],is_occluded[i], locToGlobNodeMap[j]+1, status[j],is_occluded[j],l);
    fprintf(stderr,"There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
                   i, status[i],is_occluded[i], j, status[j],is_occluded[j],l);
    exit(-1);
  }
  
  LevelSetResult lsRes;

  lsRes.alpha = 0.5;
  lsRes.xi[0] = 1.0/3.0;
  lsRes.xi[1] = 1.0/3.0;
  lsRes.xi[2] = 1.0/3.0;
  
  lsRes.trNodes[0] = -1;
  lsRes.trNodes[1] = -1;
  lsRes.trNodes[2] = -1;

  lsRes.normVel = 0;

  lsRes.gradPhi = 0; 

  return lsRes;
}

int MultiGridLevelSetStructure::numOfFluids() {
  
  return distLSS.numOfFluids();
}


int
DistMultiGridLevelSetStructure::
recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    subLSS[iSub]->recompute();
  }

  myLevel->assembleMax(*status);
  
  int nOfF = numOfFluids();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    int N = (*is_active)(iSub).size();
    for (int i = 0; i < N; ++i) 
      (*is_active)(iSub)[i] = (*status)(iSub)[i];
  }

  *is_occluded = false;
  
  
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    subLSS[iSub]->computeEdgeCrossing();
  }

  return 0;
}

DistMultiGridLevelSetStructure::
DistMultiGridLevelSetStructure(IoData &iod, Communicator *comm,
			       DistLevelSetStructure* parent,
			       MultiGridLevel<double>* level) : DistLevelSetStructure(),
								myLevel(level),
								parent(parent),
								com(comm) {

  

}

void DistMultiGridLevelSetStructure::
initialize(Domain * d, DistSVec<double,3> &X, DistSVec<double,3> &Xn, 
	   IoData &iod, DistVec<int> *point_based_id, 
	   DistVec<int>* oldStatus) {

  
  domain = d;
  numLocSub = d->getNumLocSub();
  subLSS = new MultiGridLevelSetStructure*[numLocSub];

  //closest = new DistVec<ClosestPoint>(myLevel->getNodeDistInfo()); //needed only for multi-phase cracking.
  status = new DistVec<int>(myLevel->getNodeDistInfo());
  distance = new DistVec<double>(myLevel->getNodeDistInfo());
  is_swept = new DistVec<bool>(myLevel->getNodeDistInfo());
  is_active = new DistVec<bool>(myLevel->getNodeDistInfo());
  is_occluded = new DistVec<bool>(myLevel->getNodeDistInfo());
  edge_intersects = new DistVec<bool>(myLevel->getEdgeDistInfo());


#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i)
    subLSS[i] = new MultiGridLevelSetStructure(*this,*domain->getSubDomain()[i],
					       (*status)(i), (*distance)(i),
					       (*is_active)(i),
					       (*is_swept)(i),
					       (*is_occluded)(i),
					       (*edge_intersects)(i),
					       &(*parent)(i), i, 
					       myLevel);


  *distance=0.0;

}

LevelSetStructure &
DistMultiGridLevelSetStructure::
operator()(int subNum) const {
  return *subLSS[subNum];
}
