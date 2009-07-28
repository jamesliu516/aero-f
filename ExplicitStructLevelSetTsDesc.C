
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include <Domain.h>

#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
ExplicitStructLevelSetTsDesc<dim>::
ExplicitStructLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  StructLevelSetTsDesc<dim>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()),
  k3(this->getVecInfo()), k4(this->getVecInfo()),
  p1(this->getVecInfo()), p2(this->getVecInfo()),
  p3(this->getVecInfo()), p4(this->getVecInfo()),
  U0(this->getVecInfo()), Phi0(this->getVecInfo())
{
  timeType = ioData.ts.expl.type;

  //initialize mmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) {
    MeshMotionHandler *_mmh = 0;
    _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    this->mmh = _mmh;
  } else this->mmh = 0;
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitStructLevelSetTsDesc<dim>::~ExplicitStructLevelSetTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitStructLevelSetTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim>& U)
// the key player.
{solveNLSystemOneBlock(U);} //let's just do a Forward Euler for now.

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{solveNLAllFE(U);} // Forward Euler.

//------------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLAllFE(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if(this->LS) 
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  if(this->TYPE==1 && this->mmh) { //recompute intersections and update phase change.

    if (this->Weights && this->VWeights)
      if (this->phaseChangeChoice==0)
        this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                       *this->VWeights, this->distLSS);
      else if (this->phaseChangeChoice==1)
        this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij, 
                                                       *this->Wstarji, *this->Weights, *this->VWeights,
                                                       this->distLSS);

    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts); //TODO: should do this only for the unsteady case.

    if (this->Weights && this->VWeights)
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS);
  }

  computeRKUpdate(U, k1, 1);  

  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->timer->addFluidSolutionTime(t0);

  if (this->TYPE==2) {
    if(this->LS) {
      this->spaceOp->storePreviousPrimitive(U0, this->Vg, this->Phi,
                                            this->Vgf, this->Vgfweight);
      t0 = this->timer->getTime();
  
      computeRKUpdateLS(this->Phi, p1, U);
      this->Phi = this->Phi - p1;

      this->timer->addLevelSetSolutionTime(t0);

      U = U0;
      this->updateFSInterface();
      this->updateNodeTag();
      // Riemann overwrite using the value of Phi_{n+1}
      this->spaceOp->updatePhaseChange(this->Vg, U, this->Phi,
                                       this->LS->Phin, this->Vgf,
                                       this->Vgfweight, this->riemann);
    } else {
      U = U0;
      this->updateFSInterface();
      this->updateNodeTag();
      this->spaceOp->updatePhaseChange(*this->X, U, *this->Wstarij, *this->Wstarji, this->distLSS, this->nodeTag0, this->nodeTag);
    }

    checkSolution(U);

    if(this->LS) {
      this->boundaryFlux  = *this->tmpDistSVec;
      this->interfaceFlux = *this->tmpDistSVec2;
    }
  } else { //fluid-fullbody
    U = U0;
    checkSolution(U);
  }

}

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                  DistSVec<double,dim>& dU, int it)
{
  // manually cast DistVec<int> to DistVec<double>. TODO: should avoid doing this.
  if (this->TYPE==2) {
    DistVec<double> nodeTagCopy(this->getVecInfo());
#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      Vec<double> &subTagCopy = nodeTagCopy(iSub);
      Vec<int> &subTag = this->nodeTag(iSub);
      for (int iNode=0; iNode<subTag.size(); iNode++) subTagCopy[iNode] = (double)subTag[iNode];
    }

    if(this->LS) {
      *this->tmpDistSVec  = 0.0;
      *this->tmpDistSVec2 = 0.0;
    }
    this->spaceOp->applyBCsToSolutionVector(Ulocal);
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   nodeTagCopy, dU, this->riemann,it);
    this->timeState->multiplyByTimeStep(dU);

  } else {
    
    this->spaceOp->applyBCsToSolutionVector(Ulocal);
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   dU, this->riemann,it);
    this->timeState->multiplyByTimeStep(dU);
  }
}


//------------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::computeRKUpdateLS(DistVec<double> &Philocal,
                                  DistVec<double> &dPhi, DistSVec<double,dim> &U)
{
  if (this->TYPE!=2) {
    fprintf(stderr,"in ExplicitStructLevelSetTsDesc::computeRKUpdateLS, shouldn't call me! abort.\n");
    exit(-1);
  }
  this->spaceOp->computeResidualLS(*this->X, *this->A, Philocal, U, dPhi);
  // for RK2 on moving grids
  this->timeState->multiplyByTimeStep(dPhi);

}

//------------------------------------------------------------------------------





