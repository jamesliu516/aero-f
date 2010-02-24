
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

  if(ioData.ts.expl.type == ExplicitData::FORWARD_EULER) FE = true;
  else FE = false;

  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i];
   
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
{
  solveNLSystemOneBlock(U);
  return 1;
} //so far only consider one system to solve (there is no level-set)

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{
  if(!FE && this->TYPE==1) {
    //this->com->fprintf(stderr,"Using Runge-Kutta 2.\n");
    solveNLAllRK2(U);
  } else
    solveNLAllFE(U); 
} 

//------------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLAllFE(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());

  //----------------------------------------------------
  if(this->TYPE==1 && this->mmh) { 
    //recompute intersections and update phase change.
    double tw = this->timer->getTime();
    if (this->Weights && this->VWeights)
      if (this->phaseChangeChoice==0)
        this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                       *this->VWeights, this->distLSS);
      else if (this->phaseChangeChoice==1)
        this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij, 
                                                       *this->Wstarji, *this->Weights, *this->VWeights,
                                                       this->distLSS);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);


    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

    tw = this->timer->getTime();
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts); //TODO: should do this only for the unsteady case.
    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    tw = this->timer->getTime();
    if (this->Weights && this->VWeights)
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, vfar);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  }
  //----------------------------------------------------

  computeRKUpdate(U, k1, 1);  

  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U0);

  //----------------------------------------------------
  if (this->TYPE==2) { 
    // fluid-shell-fluid
    U = U0;
    this->updateFSInterface();
    this->updateNodeTag();
    this->spaceOp->updatePhaseChange(*this->X, U, *this->Wstarij, *this->Wstarji, this->distLSS, this->nodeTag0, this->nodeTag);
    checkSolution(U);
  } 
  //--------------------------------------------------

  else { //fluid-fullbody
    U = U0;
    checkSolution(U);
  }
  this->timer->addFluidSolutionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::solveNLAllRK2(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();
  DistSVec<double,dim> Ubc(this->getVecInfo());

  //recompute Intersections and updatePhaseChange
  if(this->TYPE==1 && this->mmh) {
    // 1. compute weights
    double tw = this->timer->getTime();
    if (this->Weights && this->VWeights)
      if (this->phaseChangeChoice==0)
        this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                       *this->VWeights, this->distLSS);
      else if (this->phaseChangeChoice==1)
        this->spaceOp->computeRiemannWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Wstarij,
                                                       *this->Wstarji, *this->Weights, *this->VWeights,
                                                       this->distLSS);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    // 2. get dts (mmh->update(...) does nothing but return dts)
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);
    
    // 3. recompute intersections
    tw = this->timer->getTime();
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts); 
    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    // 4. update phase change
    tw = this->timer->getTime();
    if (this->Weights && this->VWeights)
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, vfar);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  }

  computeRKUpdate(U, k1, 1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  computeRKUpdate(U0, k2, 1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);

  U = U - 1.0/2.0 * (k1 + k2);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);

  checkSolution(U);
  this->timer->addFluidSolutionTime(t0);
}

//-----------------------------------------------------------------------------

template<int dim>
void ExplicitStructLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                  DistSVec<double,dim>& dU, int it)
{
/*  if(it==1) {//set vfar
    DistSVec<double,dim> Vlocal(Ulocal);
    this->varFcn->conservativeToPrimitive(Ulocal, Vlocal); 
    for(int iDim=0; iDim<dim; iDim++)
      vfar[iDim] = Vlocal(0)[0][iDim];
    sleep(1);
    fprintf(stderr,"on Proc %d, vfar is set to %e %e %e %e %e\n", this->com->cpuNum(), vfar[0], vfar[1], vfar[2], vfar[3], vfar[4]);
    sleep(1);
  }
*/
  if (this->TYPE==2) {

    this->spaceOp->applyBCsToSolutionVector(Ulocal);
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   this->nodeTag, dU, this->riemann,it);
    this->timeState->multiplyByTimeStep(dU);

  } else {
    
    this->spaceOp->applyBCsToSolutionVector(Ulocal);
    this->spaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   this->linRecAtInterface, dU, this->riemann, this->riemannNormal, it);
    this->timeState->multiplyByTimeStep(dU);
    this->timeState->multiplyByPreconditioner(Ulocal,dU);
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





