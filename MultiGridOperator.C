/* MultiGridOperator.C

 */

#include <MultiGridOperator.h>
#include <MultiGridLevel.h>

template <class Scalar,int dim>
MultiGridOperator<Scalar,dim>::MultiGridOperator(MultiGridLevel<Scalar>* mg_lvl,
                                                 IoData& ioData, VarFcn* varFcn, 
                                                 Domain* domain) {

  mgLevel = mg_lvl;
  myBcData = new DistBcData<dim>(ioData, varFcn, domain, &mg_lvl->getNodeDistInfo(),
                                 &mg_lvl->getInletNodeDistInfo(), &mg_lvl->getFaceDistInfo());
  timeState = new DistTimeState<dim>(ioData, NULL, varFcn, domain, mg_lvl->getNodeDistInfo(),
                                     NULL);

  timeState->createSubStates();

  zero = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  scalar_zero = new DistVec<Scalar>(mg_lvl->getNodeDistInfo());

  idti = new DistVec<Scalar>(mg_lvl->getNodeDistInfo()); 
  idtv = new DistVec<Scalar>(mg_lvl->getNodeDistInfo()); 
}

template<class Scalar,int dim>
template <class Scalar2,int neq>
void MultiGridOperator<Scalar,dim>::computeJacobian(DistSVec<Scalar2,dim>& U, DistSVec<Scalar2,dim>& V,
    //                                         DistVec<Scalar2>& irey,
                                             FluxFcn **fluxFcn,
                                             //FemEquationTerm* fet,
                                             DistMvpMatrix<Scalar2,neq>& matrices) {

  DistVec<Scalar2>& irey = *scalar_zero; 
#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    matrices(iSub) = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    mgLevel->getEdges()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, (mgLevel->getGeoState())(iSub),irey(iSub),
                                                 mgLevel->getXn()(iSub),
                                                 mgLevel->getCtrlVol()(iSub), 
                                                 V(iSub), matrices(iSub));
    mgLevel->getFaces()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn,(*myBcData)(iSub), (mgLevel->getGeoState())(iSub),
                                                 V(iSub), matrices(iSub));
 
    Vec<double>& ctrlVol = mgLevel->getCtrlVol()(iSub); 
    for (int i=0; i<ctrlVol.size(); ++i) {
      Scalar voli = 1.0 / ctrlVol[i];
      Scalar2 *Aii = matrices(iSub).getElem_ii(i);
      for (int k=0; k<neq*neq; ++k)
        Aii[k] *= voli;
    }

  }

  mgLevel->assemble(matrices);

  if (timeState) {

    timeState->addToJacobian(mgLevel->getCtrlVol(), matrices, U);
  } 
}

template<class Scalar,int dim>
template <class Scalar2>
void MultiGridOperator<Scalar,dim>::computeResidual(DistSVec<Scalar2,dim>& V,
                                                DistSVec<Scalar2,dim>& U,
  //                                           DistVec<Scalar2>& irey,
                                             FluxFcn** fluxFcn,
                                             RecFcn* recFcn,
                                             DistSVec<Scalar2,dim>& res,
                                             bool addDWdt) {


  ElemSet dummy;
  res = 0.0;
  DistVec<Scalar2>& irey = *scalar_zero; 
#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    NodalGrad<dim,Scalar2> ngrad((*zero)(iSub), (*zero)(iSub), (*zero)(iSub));
    mgLevel->getEdges()[iSub]->template computeFiniteVolumeTerm<dim>(NULL, irey(iSub), fluxFcn,
                                         recFcn, dummy, (mgLevel->getGeoState())(iSub),
                                         (mgLevel->getGeoState().getXn())(iSub), V(iSub), ngrad, NULL, res(iSub),
                                         (mgLevel->getFVCompTag())(iSub), 0,0); 


    mgLevel->getFaces()[iSub]->computeFiniteVolumeTerm(fluxFcn, (*myBcData)(iSub), 
                                         (mgLevel->getGeoState())(iSub), V(iSub),
                                         res(iSub));
    
  }
  
  mgLevel->assemble(res);

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    for (int i = 0; i < V.subSize(iSub); ++i) {

      for (int j = 0; j < dim; ++j)
        res(iSub)[i][j] /= (mgLevel->getCtrlVol())(iSub)[i];
    }
  }
  
  if (timeState && addDWdt) {

    timeState->add_dAW_dt(0, mgLevel->getGeoState(),
                          mgLevel->getCtrlVol(), U, res);
  } 
  
  
}                                             

template<class Scalar,int dim>
void MultiGridOperator<Scalar,dim>::updateStateVectors(DistSVec<Scalar,dim>& U) {

  timeState->update(U);
}

template<class Scalar,int dim>
void MultiGridOperator<Scalar,dim>::computeTimeStep(double cfl, VarFcn *varFcn,
                                                    DistSVec<double,dim> &V) {

  int iSub;
  DistVec<double>& dt = timeState->getDt();
  dt = 0.0;
  *idti = 0.0;
  *idtv = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < V.numLocSub(); ++iSub) {
    mgLevel->getEdges()[iSub]->computeTimeStep(NULL, varFcn,  (mgLevel->getGeoState())(iSub),
                                 (mgLevel->getGeoState().getXn())(iSub),
                                 V(iSub), (*idti)(iSub), (*idtv)(iSub), 
                                 timeState->getTimeLowMachPrec());
    mgLevel->getFaces()[iSub]->computeTimeStep(NULL, varFcn,  (mgLevel->getGeoState())(iSub),
                                 (mgLevel->getGeoState().getXn())(iSub),
                                 V(iSub), (*idti)(iSub), (*idtv)(iSub), 
                                 timeState->getTimeLowMachPrec());
  }

  mgLevel->assemble(*idti);
  *idtv = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < V.numLocSub(); ++iSub) {
    double (*idtimev) = idtv->subData(iSub);
    double (*idtimei) = idti->subData(iSub);
    double (*dtime) = dt.subData(iSub);
    //double (*ireynolds) = irey.subData(iSub);
    double (*volume) = mgLevel->getCtrlVol().subData(iSub);
    for (int i = 0; i < mgLevel->getCtrlVol().subSize(iSub); ++i) {
      //   idtimev[i] = idtimev[i] / volume[i];
      dtime[i] = cfl *volume[i]/(-1.0*idtimei[i]/* + viscous*idtimev[i]*/);
      //ireynolds[i] = -sprec.getViscousRatio()*idtimev[i] / idtimei[i];
    }
  }

  double dt_glob = dt.min();
  timeState->computeCoefficients(dt_glob);
}

#define INST_HELPER(S,D) \
template \
void MultiGridOperator<S,D>::computeJacobian(DistSVec<double,D>& U, DistSVec<double,D>& V,\
                                             FluxFcn **fluxFcn, \
                                             DistMvpMatrix<double,D>& matrices); \
template void MultiGridOperator<S,D>::computeResidual(DistSVec<double,D>& V, \
                                                DistSVec<double,D>& U, \
                                             FluxFcn** fluxFcn,\
                                             RecFcn* recFcn,\
                                             DistSVec<double,D>& res,bool);

template class MultiGridOperator<double,1>;
template class MultiGridOperator<double,2>;
template class MultiGridOperator<double,5>;
template class MultiGridOperator<double,6>;
template class MultiGridOperator<double,7>;

INST_HELPER(double,1);
INST_HELPER(double,2);
INST_HELPER(double,5);
INST_HELPER(double,6);
INST_HELPER(double,7);

