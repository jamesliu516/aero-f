#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <MatVecProd.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, int **nodeType, BCApplier *bcs)
    : num_levels(5), agglom_size(8), numLocSub(dom->getNumLocSub()), multiGridLevels(new MultiGridLevel<Scalar2>*[num_levels+1]),
    macroValues(new DistSVec<Scalar2,dim>*[num_levels]), geoState(distGeoState)
{

  smoothingMatrices = new MultiGridSmoothingMatrix<Scalar2,dim>**[num_levels]; 

  nSmooth = 1;

  relaxationFactor = 1.0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    dom->getSubDomain()[iSub]->getEdges().updateLength(distGeoState.getXn()(iSub));
    dom->getSubDomain()[iSub]->makeMasterFlag(dom->getNodeDistInfo());
  }

  multiGridLevels[0] = new MultiGridLevel<Scalar2>(*dom, dom->getNodeDistInfo(), dom->getEdgeDistInfo());
  multiGridLevels[0]->copyRefinedState(dom->getNodeDistInfo(), dom->getEdgeDistInfo(), geoState, *dom);

  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1] = new MultiGridLevel<Scalar2>(*dom, multiGridLevels[level]->getNodeDistInfo(), multiGridLevels[level]->getEdgeDistInfo());
    multiGridLevels[level+1]->agglomerate(multiGridLevels[level]->getNodeDistInfo(),
                                          multiGridLevels[level]->getIdPat(),
                                          multiGridLevels[level]->getDistGeoState(),
                                          multiGridLevels[level]->getSharedNodes(),
                                          multiGridLevels[level]->getConnectivity(),
                                          multiGridLevels[level]->getEdges(),
                                          *dom,dim);
  }

  macroA = new DistMat<Scalar2,dim>*[num_levels+1];

  macroValues = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroR = new DistSVec<Scalar2,dim>*[num_levels+1];
  macroDX = new DistSVec<Scalar2,dim>*[num_levels+1];
  for(int level = 0; level <= num_levels; ++level) {
    macroValues[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroR[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
    macroDX[level] = new DistSVec<Scalar2,dim>(multiGridLevels[level]->getNodeDistInfo());
  }
  
  for(int level = 0; level < num_levels; ++level) {

    smoothingMatrices[level] = new MultiGridSmoothingMatrix<Scalar2,dim>*[numLocSub];
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub)
      smoothingMatrices[level][iSub] = new
        MultiGridSmoothingMatrix<Scalar2,dim>(MultiGridSmoothingMatrix<Scalar2,dim>::LineJacobi,iSub,
                                              multiGridLevels[level]->getNodeDistInfo().subSize(iSub),
                                              multiGridLevels[level]->getEdgeDistInfo().subSize(iSub),
                                              0,
                                              multiGridLevels[level+1],
                                              multiGridLevels[level]);
  }                      
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
{
#pragma omp parallel for
  for (int level = 0; level <= num_levels; ++level) {
    delete macroValues[level];
    delete macroR[level];
    delete macroDX[level];
    delete multiGridLevels[level];
  }
  delete []macroA;
  delete []macroValues;
  delete []macroR;
  delete []macroDX;
  delete []multiGridLevels;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup()
{
  for(int level = 0; level < num_levels; ++level) {
    multiGridLevels[level+1]->computeRestrictedQuantities(multiGridLevels[level]->getDistGeoState());
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  //*macroValues[0] = x;
  *macroValues[0] = 0.0;
  smooth(0, *macroValues[0], x);
  for(int level = 0; level < num_levels; ++level) {
  //  std::cout << "level = " << level << std::endl;
    multiGridLevels[level+1]->Restrict(*multiGridLevels[level], *macroValues[level], *macroValues[level+1]);
  }
  for(int level = num_levels; level > 0; --level) {
    multiGridLevels[level]->Prolong(*multiGridLevels[level-1], *macroValues[level], *macroValues[level], *macroValues[level-1]);
  }
  smooth(0, *macroValues[0], x);

  Px = *macroValues[0];
}

//------------------------------------------------------------------------------
template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::smooth(int level,DistSVec<Scalar2,dim>& x,
                                               DistSVec<Scalar2,dim>& f) {

  double rcurr;
  for (int i = 0; i < nSmooth; ++i) {

//    std::cout << "Smooth step " << i << std::endl;
    MatVecProdH1<dim,Scalar2,dim>* mvp = dynamic_cast< MatVecProdH1<dim,Scalar2,dim>*>(macroA[level]);

    if (mvp == NULL) {

      std::cout << "Error: multigrid preconditioner must be used with an H1 product" << std::endl;
      exit(1);
    }
 
    // Compute the residual
    mvp->apply(x, *macroR[level]);
    (*macroR[level]) *= -1.0;
    (*macroR[level]) += f;

    double rnew = (*macroR[level]).norm();
    // bailout?
    if (i > 0 && rnew > rcurr && x.norm() > 0.0) {
      x -= relaxationFactor*(*macroDX[level]);
      break;
    } else {
      rcurr = rnew;
    }

//    std::cout << "i = " << i << " r = " << (*macroR[level]).norm() << std::endl;
   // std::cout << "dot = " << (*macroR[level])*f / (sqrt((*macroR[level]).norm() * f.norm())) << std::endl;

#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      
      smoothingMatrices[0][iSub]->smooth((*macroR[level])(iSub),(*macroDX[level])(iSub)); 
    }

    multiGridLevels[level]->assemble(*macroDX[level]);
  //  std::cout << "DX = " << (*macroDX[level]).norm() << std::endl; 

//    *macroDX[level] = *macroR[level];
    x += relaxationFactor*(*macroDX[level]);
  }
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::getData(DistMat<Scalar2,dim>& mat) {

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    smoothingMatrices[0][iSub]->getData(mat(iSub));
  }

  macroA[0] = &mat;
  
}


template class MultiGridPrec<float, 1, double>;
template class MultiGridPrec<float, 2, double>;
template class MultiGridPrec<float, 5, double>;
template class MultiGridPrec<float, 6, double>;
template class MultiGridPrec<float, 7, double>;

