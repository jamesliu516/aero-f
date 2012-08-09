#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <DistMacroCell.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, int **nodeType, BCApplier *bcs)
    : num_levels(8),agglom_size(8),numLocSub(dom->getNumLocSub()), macroCells(0), geoState(distGeoState),
    nodeDistInfo(0), edgeDistInfo(0)
{
  bool ** masterFlag = new bool *[numLocSub];

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; iSub++)
    masterFlag[iSub] = dom->getNodeDistInfo().getMasterFlag(iSub);

  macroCells = new DistMacroCellSet(dom, 1.0, masterFlag, agglom_size, num_levels);
  macroValues = new DistSVec<Scalar2,dim>*[num_levels];
  nodeDistInfo = new DistInfo*[num_levels];

  const DistInfo& nDistInfo(dom->getNodeDistInfo());

  for(int level = 0; level < num_levels; ++level) {
    nodeDistInfo[level] = new DistInfo(nDistInfo.numLocThreads, nDistInfo.numLocSub,
                                       nDistInfo.numGlobSub,    nDistInfo.locSubToGlobSub, nDistInfo.com);
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub)
      nodeDistInfo[level]->setLen(iSub, macroCells->obtainMacroCell(iSub,level)->size());
    nodeDistInfo[level]->finalize(true);

    macroValues[level] = new DistSVec<Scalar2,dim>(*nodeDistInfo[level]);
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; iSub++) {
    delete nodeDistInfo[iSub];
    delete edgeDistInfo[iSub];
    delete macroValues[iSub];
  }

  delete []nodeDistInfo;
  delete []edgeDistInfo;
  delete []macroValues;
  delete macroCells;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup()
{
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  static bool doInitialTasks = true;

  macroCells->computeVBar(doInitialTasks, geoState, x, *macroValues[0], 0, 1);
  for(int level = 1; level < num_levels; ++level) {
  }
  doInitialTasks = false;

  // Not really even a preconditioner at all!
  Px = x;
}

//------------------------------------------------------------------------------

template class MultiGridPrec<float, 1, double>;
template class MultiGridPrec<float, 2, double>;
template class MultiGridPrec<float, 5, double>;
template class MultiGridPrec<float, 6, double>;
template class MultiGridPrec<float, 7, double>;
