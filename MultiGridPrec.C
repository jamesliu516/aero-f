#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, int **nodeType, BCApplier *bcs)
    : num_levels(5), agglom_size(8), numLocSub(dom->getNumLocSub()), multiGridLevels(new MultiGridLevel<Scalar2>*[num_levels+1]),
    macroValues(new DistSVec<Scalar2,dim>*[num_levels]), geoState(distGeoState)
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    dom->getSubDomain()[iSub]->getEdges().updateLength(distGeoState.getXn()(iSub));
    dom->getSubDomain()[iSub]->makeMasterFlag(dom->getNodeDistInfo());
  }

  multiGridLevels[0] = new MultiGridLevel<Scalar2>(dom->getNodeDistInfo(), dom->getEdgeDistInfo());
  multiGridLevels[0]->copyRefinedState(dom->getNodeDistInfo(), dom->getEdgeDistInfo(), geoState.getXn(), *dom);

  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1] = new MultiGridLevel<Scalar2>(multiGridLevels[level]->getNodeDistInfo(), multiGridLevels[level]->getEdgeDistInfo());
    multiGridLevels[level+1]->agglomerate(multiGridLevels[level]->getNodeDistInfo(),
                                          multiGridLevels[level]->getIdPat(),
                                          multiGridLevels[level]->getSharedNodes(),
                                          multiGridLevels[level]->getConnectivity(),
                                          multiGridLevels[level]->getEdges(),
                                          *dom);
  }

  macroValues = new DistSVec<Scalar2,dim>*[num_levels+1];
  for(int level = 0; level <= num_levels; ++level)
    macroValues[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
{
#pragma omp parallel for
  for (int level = 0; level <= num_levels; ++level) {
    delete macroValues[level];
    delete multiGridLevels[level];
  }
  delete []macroValues;
  delete []multiGridLevels;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup()
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    multiGridLevels[0]->setCtrlVolumes(iSub, geoState(iSub).getCtrlVol_n());

  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1]->computeRestrictedQuantities(multiGridLevels[level]->getVol(), multiGridLevels[level]->getX());
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  *macroValues[0] = x;
  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroValues[level], *macroValues[level+1]);
  }
  for(int level = num_levels; level > 0; --level) {
    multiGridLevels[level]->Prolong(*multiGridLevels[level-1], *macroValues[level], *macroValues[level], *macroValues[level-1]);
  }
  Px = *macroValues[0];
}

//------------------------------------------------------------------------------

template class MultiGridPrec<float, 1, double>;
template class MultiGridPrec<float, 2, double>;
template class MultiGridPrec<float, 5, double>;
template class MultiGridPrec<float, 6, double>;
template class MultiGridPrec<float, 7, double>;
