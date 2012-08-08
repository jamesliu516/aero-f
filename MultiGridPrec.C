#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <MatVecProd.h>
#include <DistMvpMatrix.h>

#include <IoData.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, PcData& pcData,KspData& coarseSolverData, IoData& ioData, VarFcn* varFcn,bool createFineA, 
DistTimeState<dim>* ts,int **nodeType, BCApplier *bcs) : pcData(pcData), DistMat<Scalar2,dim>(dom)
{

  mgKernel = new MultiGridKernel<Scalar,dim,Scalar2>(dom, distGeoState, coarseSolverData, ioData, varFcn,createFineA,pcData.num_multigrid_levels, ts,nodeType, bcs);
  mgKernel->setParameters(pcData.num_multigrid_smooth1, 
                          pcData.num_multigrid_smooth2,
                          pcData.num_fine_sweeps,
                          pcData.mg_smooth_relax, 
                          pcData.mg_output);

}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::initialize() {

  mgKernel->initialize();
  
  if (pcData.mg_type == PcData::MGGEOMETRIC)
    mgKernel->setGeometric();
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
{
  delete mgKernel;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setOperators(SpaceOperator<dim>* spo) {

  mgKernel->setOperators(spo);
}

//------------------------------------------------------------------------------
template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup() {

  std::cout << "Error: wrong setup called.  For multigrid you should call "
               "MultiGridPrec<Scalar,dim,Scalar2>::setup(DistSVec<Scalar2,dim>& V)" << std::endl;
  exit(1);
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup(DistSVec<Scalar2,dim>& U)
{
  if (pcData.mg_type == PcData::MGGEOMETRIC)
    mgKernel->setupGeometric(U);
  else
    mgKernel->setupAlgebraic();
}

template<class Scalar, int dim, class Scalar2>
bool MultiGridPrec<Scalar,dim,Scalar2>::isInitialized() {

  return mgKernel->isInitialized();
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  mgKernel->cycle(x, Px);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::getData(DistMat<Scalar2,dim>& mat)/* ,
                                                DistSVec<Scalar2,dim>& V,
                                                SpaceOperator<dim>& spo,
                                                DistTimeState<dim>* timeState)*/
{
//  VarFcn* varFcn = spo.getVarFcn();

  mgKernel->getData(mat);
  
/*  *macroValues[0] = V;
  if (timeState)
    *macroIrey[0] = *timeState->getInvReynolds();
  else
    *macroIrey[0] = 0.0;
*/
  //std::cout << "Vnorm 0 = " << V*V <<std::endl;
  //smooth(0, *macroValues[0], V);
/*#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < macroValues[0]->subSize(iSub); ++i)
      std::cout << i << " " << (*macroValues[0])(iSub)[i][0] << " " << 
        (*timeState)(iSub).getDt()[i] << std::endl;
  }
*/
}

template<class Scalar, int dim, class Scalar2>
DistMat<Scalar2,dim>& MultiGridPrec<Scalar,dim,Scalar2>::
operator= (const Scalar2 s) {

  return (mgKernel->getFineMatrix() = s);
}

template<class Scalar, int dim, class Scalar2>
GenMat<Scalar2,dim>& MultiGridPrec<Scalar,dim,Scalar2>::
operator() (int i) {

  return mgKernel->getFineMatrix()(i);
}

//------------------------------------------------------------------------------
#define INSTANTIATION_HELPER(T,dim) \
    template class MultiGridPrec<T, dim, double>;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
INSTANTIATION_HELPER(float,7);

INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
