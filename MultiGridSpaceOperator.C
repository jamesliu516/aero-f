#include <MultiGridSpaceOperator.h>

template <class Scalar,int dim>
MultiGridSpaceOperator<Scalar,dim>::
MultiGridSpaceOperator(IoData& ioData,Domain* domain, 
                       SpaceOperator<dim>* spo,MultiGridKernel<Scalar>* pKernel,
                       SpaceOperator<dim>* spo1) : pKernel(pKernel) {

  nLevels = pKernel->numLevels();

  myOperators = new MultiGridOperator<Scalar,dim>*[nLevels];
  BcFcn* bcFcn1 = spo->getBcFcn();
  if (spo1)
    bcFcn1 = spo1->getBcFcn();
  for (int i = 1; i < nLevels; ++i) {

    myOperators[i] = 
      new MultiGridOperator<Scalar,dim>(pKernel->getLevel(i),
                                        ioData,spo->getVarFcn(),
                                        domain, spo->getBcFcn(),
                                        bcFcn1);
                                        
  }

  fluxFcn = spo->getFluxFcn();

  fet = spo->getFemEquationTerm();

  if (spo1) {
    fluxFcn1 = spo1->getFluxFcn();
    fet1 = spo1->getFemEquationTerm();
  }
  else {
    fluxFcn1 = fluxFcn;
    fet1 = fet;
  }

  varFcn = spo->getVarFcn();
}

template <class Scalar,int dim>
MultiGridSpaceOperator<Scalar,dim>::
~MultiGridSpaceOperator() {

  for (int i = 1; i < nLevels; ++i) {
    delete myOperators[i];
  }

  delete [] myOperators;
}

template <class Scalar,int dim>
void MultiGridSpaceOperator<Scalar,dim>::setupBcs(DistBcData<dim>* bcd) {

  for(int level = 1; level < nLevels; ++level) {                                                  
                                                                                                      
    if (level == 1)
      pKernel->getLevel(level)->setupBcs(*bcd, myOperators[level]->getBcData(),myOperators[level]->getBoundaryState());
    else
      pKernel->getLevel(level)->setupBcs(myOperators[level-1]->getBcData(), myOperators[level]->getBcData(),myOperators[level]->getBoundaryState());
  }
}

template <class Scalar,int dim>
void MultiGridSpaceOperator<Scalar,dim>::
computeTimeStep(int level, double cfl, 
                MultiGridDistSVec<Scalar,dim>& V) {

  myOperators[level]->computeTimeStep(cfl, varFcn, V(level));
}


template <class Scalar,int dim> 
void MultiGridSpaceOperator<Scalar,dim>::
computeResidual(int level, MultiGridDistSVec<Scalar,dim>& U,
                MultiGridDistSVec<Scalar,dim>& V,
                MultiGridDistSVec<Scalar,dim>& res,
                bool addDWdt) {

  myOperators[level]->computeResidual(V(level), U(level),
                                      fluxFcn, &recConstant,
                                      fet,res(level), addDWdt);
  myOperators[level]->applyBCsToResidual(U(level),res(level));
}

template <class Scalar,int dim> 
void MultiGridSpaceOperator<Scalar,dim>::
updateStateVectors(int level, MultiGridDistSVec<Scalar,dim>& U) {

  myOperators[level]->updateStateVectors(U(level));
}

template <class Scalar,int dim> 
template <int neq>
void MultiGridSpaceOperator<Scalar,dim>::
computeJacobian(int lvl, MultiGridDistSVec<Scalar,dim>& U,
                MultiGridDistSVec<Scalar,dim>& V,
                MultiGridMvpMatrix<Scalar,neq>& mvp) {

  myOperators[lvl]->computeJacobian(U(lvl), V(lvl),
                                      fluxFcn1, fet1,mvp(lvl));

  myOperators[lvl]->applyBCsToJacobian(U(lvl),mvp(lvl));
}

template class MultiGridSpaceOperator<double,1>;
template class MultiGridSpaceOperator<double,2>;
template class MultiGridSpaceOperator<double,5>;
template class MultiGridSpaceOperator<double,6>;
template class MultiGridSpaceOperator<double,7>;

#define INST_HELPER(dim,neq) \
template void MultiGridSpaceOperator<double,dim>:: \
computeJacobian(int level, MultiGridDistSVec<double,dim>&,\
MultiGridDistSVec<double,dim>&, MultiGridMvpMatrix<double,neq>&);

INST_HELPER(5,5);
INST_HELPER(6,6);
INST_HELPER(7,7);

INST_HELPER(6,5);
INST_HELPER(7,5);
