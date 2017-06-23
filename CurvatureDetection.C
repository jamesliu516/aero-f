#include <CurvatureDetection.h>
#include <Domain.h>
#include <SubDomain.h>
#include <Vector3D.h>
#include <DistVector.h>

//------------------------------------------------------------------------------

CurvatureDetection::CurvatureDetection(Domain* domain)
{

  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();
//  tag = new DistVec<double>(domain->getEdgeDistInfo());
//  vec1 = domain->getVolPat();
  normals = new DistSVec<double,6>(domain->getEdgeDistInfo());
  vec6 = new CommPattern<double>(domain->getSubTopo(),domain->getCommunicator(),
				 CommPattern<double>::CopyOnSend);
#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenEdges(6, *vec6);
  vec6->finalize();

  originAndTag = new DistSVec<double,5>(domain->getNodeDistInfo());
  vec5 = new CommPattern<double>(domain->getSubTopo(),domain->getCommunicator(),
				 CommPattern<double>::CopyOnSend);
#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenNodes(5, *vec5);
  vec5->finalize();
}

//------------------------------------------------------------------------------

CurvatureDetection::~CurvatureDetection()
{

//  if (tag) delete tag;
  if (normals) delete normals;
  if (vec6) delete vec6;
  if (originAndTag) delete originAndTag;
  if (vec5) delete vec5;

}

//------------------------------------------------------------------------------

void CurvatureDetection::compute(double threshold, int nLayers, double maxDist,
                                 DistSVec<double,3>& X, DistVec<bool>& t)
{

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeFaceEdgeNormals(X(iSub), (*normals)(iSub));
    subDomain[iSub]->sndEdgeData(*vec6, normals->subData(iSub));
  }
  vec6->exchange();
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->addRcvEdgeData(*vec6, normals->subData(iSub));
    subDomain[iSub]->computeEdgeDihedralAngle(threshold, X(iSub),  (*normals)(iSub), (*originAndTag)(iSub));
    subDomain[iSub]->sndData(*vec5, reinterpret_cast<double (*)[5]>(originAndTag->subData(iSub)));
  }

  for (int iL=0;iL<nLayers;iL++) {
    vec5->exchange();
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->otRcvData(*vec5,
          reinterpret_cast<double (*)[5]>(originAndTag->subData(iSub)));
      subDomain[iSub]->propagateInfoAlongEdges(maxDist, X(iSub),
                                               (*originAndTag)(iSub));
    }
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->sndData(*vec5,
          reinterpret_cast<double (*)[5]>(originAndTag->subData(iSub)));
    }
  }
  vec5->exchange();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->otRcvData(*vec5, reinterpret_cast<double (*)[5]>(originAndTag->subData(iSub)));
    double (*_tag)[5] = originAndTag->subData(iSub);
    bool* _t = t.subData(iSub);
    for (int i=0; i<t.subSize(iSub); ++i) {
      if (_tag[i][3] > 0.0)
	_t[i] = true;
    }
  }
}

//------------------------------------------------------------------------------
