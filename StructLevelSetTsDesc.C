#include <StructLevelSetTsDesc.h>
#include <DistExactRiemannSolver.h>
#include <LevelSet/IntersectionFactory.h>
#include <FSI/DynamicNodalTransfer.h>

#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif


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
StructLevelSetTsDesc<dim>::
StructLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), nodeTag(this->getVecInfo()), nodeTag0(this->getVecInfo()),
  Phi(this->getVecInfo()), Vg(this->getVecInfo()),
  PhiV(this->getVecInfo()), boundaryFlux(this->getVecInfo()),
  computedQty(this->getVecInfo()), interfaceFlux(this->getVecInfo()),TYPE(ioData.eqs.numPhase)
{

  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;
  forceApp = ioData.strucIntersect.forceApproach;

  this->postOp->setForceGenerator(this);
  if (TYPE==1) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-STRUCTURE SIMULATION --------\n");
  if (TYPE==2) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-SHELL-FLUID SIMULATION --------\n");

  timeStep = 0.0;
  fsiPosition = 0.0;
  fsiNormal = 0.0;
  fsiVelocity = 0.0;

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain);
  distLSS = 0;

  const char *intersectorName = ioData.strucIntersect.intersectorName;
  if(intersectorName != 0)
    distLSS = IntersectionFactory::getIntersectionObject(intersectorName, *this->domain);

  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());

  *Wstarij = 0.0;
  *Wstarji = 0.0;

  pressureRef = ioData.ref.rv.pressure;
//------ load structure mesh information ----------------------
  numStructNodes = 0;
  Fs = 0;
  numStructNodes = distLSS->getNumStructNodes();
  if (numStructNodes>0) {
    this->com->fprintf(stderr,"# of struct nodes: %d.\n", numStructNodes);
    // We allocate Fs from memory that allows fast one-sided MPI communication
    Fs = new (*this->com) double[numStructNodes][3];
  } else this->com->fprintf(stderr,"Warning: failed loading structure mesh information!\n");
//-------------------------------------------------------------

// for FF interface
//  LS = new LevelSet(ioData, this->domain);
  LS = 0; // for debug!
  Vgf = 0;
  Vgfweight = 0;
  if(ioData.mf.typePhaseChange == MultiFluidData::EXTRAPOLATION){
    Vgf = new DistSVec<double,dim>(this->getVecInfo());
    Vgfweight = new DistVec<double>(this->getVecInfo());
    *Vgf =-1.0;
    *Vgfweight =0.0;
  }

  //multiphase conservation check (only for FF interface)
  boundaryFlux  = 0.0;
  interfaceFlux = 0.0;
  computedQty   = 0.0;
  tmpDistSVec   = 0;
  tmpDistSVec2  = 0;
  for(int i=0; i<dim; i++){
    expectedTot[i] = 0.0; expectedF1[i] = 0.0; expectedF2[i] = 0.0;
    computedTot[i] = 0.0; computedF1[i] = 0.0; computedF2[i] = 0.0;
  }

  frequencyLS = ioData.mf.frequency;
  interfaceTypeFF = ioData.mf.interfaceType;
  if(this->domain->getEmbedCommunicator()) {
    Communicator *intraCom =this->domain->getEmbedCommunicator()->merge(true);
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *intraCom);
  } else
    dynNodalTransfer = 0;
}

//------------------------------------------------------------------------------

template<int dim>
StructLevelSetTsDesc<dim>::~StructLevelSetTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;

  if (Fs) delete[] Fs;

  if (LS) delete LS;
  if (Vgf) delete Vgf;
  if (Vgfweight) delete Vgfweight;
  if(tmpDistSVec)  delete tmpDistSVec;
  if(tmpDistSVec2) delete tmpDistSVec2;
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
//from TsDesc::setupTimeStepping. 
{
  this->geoState->setup2(this->timeState->getData());

  if (TYPE==1) {
    this->timeState->setup(this->input->solutions, this->bcData->getInletBoundaryVector(), *this->X, *U);
    this->distLSS->initialize(this->domain,*this->X);
  } else if (TYPE==2) {
    // load the FS interface.
    if (ioData.mf.initialConditions.nplanes != 1) {
      fprintf(stderr,"number of planes != 1! Abort...\n"); exit(-1);}
    fsiPosition[0] = ioData.mf.initialConditions.p1.cen_x;
    fsiPosition[1] = ioData.mf.initialConditions.p1.cen_y;
    fsiPosition[2] = ioData.mf.initialConditions.p1.cen_z;
    fsiNormal[0] = ioData.mf.initialConditions.p1.nx;
    fsiNormal[1] = ioData.mf.initialConditions.p1.ny;
    fsiNormal[2] = ioData.mf.initialConditions.p1.nz;
    fsiNormal /= fsiNormal.norm();
    fsiVelocity = 10.0; //TODO: should read from the input file.

    fprintf(stderr,"Position = [%e, %e, %e];  Normal = [%e, %e, %e];  Velo = %e.\n", fsiPosition[0], fsiPosition[1], fsiPosition[2], fsiNormal[0], fsiNormal[1], fsiNormal[2], fsiVelocity);

    this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), *U, nodeTag, ioData);
    nodeTag0 = nodeTag;
    this->distLSS->initialize(this->domain,*this->X);

    if(LS) LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData);
  }
}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U)
{
  // manually cast DistVec<int> to DistVec<double>. TODO: should avoid doing this.
  if (TYPE==2) {
    DistVec<double> nodeTagCopy(this->getVecInfo());
#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      Vec<double> &subTagCopy = nodeTagCopy(iSub);
      Vec<int> &subTag = nodeTag(iSub);
      for (int iNode=0; iNode<subTag.size(); iNode++) subTagCopy[iNode] = (double)subTag[iNode];
    }

    double t0 = this->timer->getTime();
    this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);
    int numSubCycles = 1;
    double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                              &numSubCycles, *this->geoState, *this->A, U,nodeTagCopy);
    if (this->problemType[ProblemData::UNSTEADY])
      this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                        dt*this->refVal->time, numSubCycles);

    this->timer->addFluidSolutionTime(t0);
    this->timer->addTimeStepTime(t0);

    timeStep = dt;
    return dt;
  } else {
    double t0 = this->timer->getTime();
    this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);
    int numSubCycles = 1;
    double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                              &numSubCycles, *this->geoState, *this->X, *this->A, U);

    if (this->problemType[ProblemData::UNSTEADY])
      this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                        dt*this->refVal->time, numSubCycles);

    this->timer->addFluidSolutionTime(t0);
    this->timer->addTimeStepTime(t0);

    return dt;
  }

}

//---------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  if(LS) {
    LS->update(Phi);
    this->timeState->update(U, LS->Phin, LS->Phinm1, LS->Phinm2,
                            Vgf, Vgfweight, riemann);
    if(frequencyLS > 0 && it%frequencyLS == 0){
      LS->conservativeToPrimitive(Phi,PhiV,U);
      LS->reinitializeLevelSet(*this->geoState,*this->X, *this->A, U, PhiV);
      LS->primitiveToConservative(PhiV,Phi,U);
    }
  } else this->timeState->update(U);
}

//-----------------------------------------------------------------------------

template<int dim>
int StructLevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{
  int ierr = 0;
  if (TYPE==2) {
    if(LS) ierr = this->domain->checkSolution(this->varFcn, *this->A, U, Phi, LS->Phin);
    else   ierr = this->domain->checkSolution(this->varFcn, U, nodeTag);
  } else ierr = this->domain->checkSolution(this->varFcn, U); //also check ghost nodes.

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::conservationErrors(DistSVec<double,dim> &U, int it)
{
  if(!tmpDistSVec)  tmpDistSVec  = new DistSVec<double,dim>(this->getVecInfo());
// computes the total mass, momentum and energy
// 1- in the whole domain regardless of which phase they belong to
// 2- in the positive phase (phi>=0.0 <=> sign= 1)
// 3- in the negative phase (phi< 0.0 <=> sign=-1)
// The computed total domain mass, momentum and energy can be compared
//     to expected values (since only what goes out and comes in need to be
//     to be accounted for).
// When considering one phase, only the mass can be compared to an expected
//     value since the flux is zero at the interface between the two fluids.
//     For the momentum and the energy, the physical flux depends on the
//     pressure at this interface. More to be said here (numerical flux among
//     other things)
  int fluid1 =  1;
  int fluid2 = -1;

  // total mass, total mass of fluid1, total mass of fluid2
  // note that there are two types of conservation errors:
  // 1 - fluxes do not cancel at the interface
  // 2 - populating of cell that changes fluid
  // Both have an effect on the total mass
  // Only the type-2 error has an effect on the total mass of fluid-i

  computedQty = U;
  computedQty *= *this->A;
  computedQty.sum(computedTot);

  this->domain->restrictionOnPhi(computedQty, Phi, *tmpDistSVec, fluid1); //tmpDistSVec reinitialized to 0.0inside routine
  tmpDistSVec->sum(computedF1);

  this->domain->restrictionOnPhi(computedQty, Phi, *tmpDistSVec, fluid2); //tmpDistSVec reinitialized to 0.0inside routine
  tmpDistSVec->sum(computedF2);

  // expected total mass, mass in fluid1, mass in fluid2
  // an expected mass is computed iteratively, that is
  // m_{n+1} = m_{n} + fluxes_at_boundaries
  if(it==0){
    for (int i=0; i<dim; i++){
       expectedTot[i] = computedTot[i];
       expectedF1[i]  = computedF1[i];
       expectedF2[i]  = computedF2[i];
    }
    return;
  }

  const double dt = this->timeState->getTime();
  double bcfluxsum[dim];

  boundaryFlux.sum(bcfluxsum);
  for (int i=0; i<dim; i++)
    expectedTot[i] += dt*bcfluxsum[i];

  this->domain->restrictionOnPhi(boundaryFlux, Phi, *tmpDistSVec, fluid1); //tmpDistSVec reinitialized to 0.0 inside routine
  tmpDistSVec->sum(bcfluxsum);
  for (int i=0; i<dim; i++)
    expectedF1[i] += dt*bcfluxsum[i];

  this->domain->restrictionOnPhi(boundaryFlux, Phi, *tmpDistSVec, fluid2); //tmpDistSVec reinitialized to 0.0 inside routine
  tmpDistSVec->sum(bcfluxsum);
  for (int i=0; i<dim; i++)
    expectedF2[i] += dt*bcfluxsum[i];

}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt, int it, double t,
                                                  DistSVec<double,dim> &U)
{ // comes from TsDesc::setupOutputToDisk
//  this->com->barrier(); this->com->fprintf(stderr,"it = %d, maxIts = %d.\n", it, this->data->maxIts);

  if (it == this->data->maxIts)
    *lastIt = true;
  else monitorInitialState(it, U);

  this->output->openAsciiFiles();
  this->timer->setSetupTime();

  // manually cast DistVec<int> to DistVec<double>. TODO: should avoid doing this.
  DistVec<double> nodeTagCopy(this->getVecInfo());
  if (TYPE==2) {
#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      Vec<double> &subTagCopy = nodeTagCopy(iSub);
      Vec<int> &subTag = nodeTag(iSub);
      for (int iNode=0; iNode<subTag.size(); iNode++) subTagCopy[iNode] = (double)subTag[iNode];
    }
  }
  // ----------------------------------------------------

  if (it == 0) {
    // First time step: compute GradP before computing forces
    this->spaceOp->computeGradP(*this->X, *this->A, U);
    if (TYPE==2) {
      this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTagCopy);
      this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTagCopy);
      this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTagCopy);
      this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTagCopy);
    } else {
      this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    }
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    if (TYPE==2) this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTagCopy);
    else this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, distLSS->getPhi());

    if(LS)  this->output->writeConservationErrors(ioData, it, t, expectedTot, expectedF1, expectedF2, computedTot, computedF1, computedF2);

    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                                             double t, double dt, DistSVec<double,dim> &U)
{
  if(LS) conservationErrors(U,it);

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;

  // manually cast DistVec<int> to DistVec<double>. TODO: should avoid doing this.
  DistVec<double> nodeTagCopy(this->getVecInfo());
  if (TYPE==2) {
#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      Vec<double> &subTagCopy = nodeTagCopy(iSub);
      Vec<int> &subTag = nodeTag(iSub);
      for (int iNode=0; iNode<subTag.size(); iNode++) subTagCopy[iNode] = (double)subTag[iNode];
    }
  }
  // ----------------------------------------------------
  if (TYPE==2) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTagCopy);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTagCopy);
    this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTagCopy);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTagCopy);
    this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTagCopy);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  } else {
    this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,
                                     distLSS->getPhi());
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }
  if(LS) {
    this->output->writeConservationErrors(ioData, it, t, expectedTot, expectedF1, expectedF2, computedTot, computedF1, computedF2);
    this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);
  } else
    this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, 0);

  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)
{ 
  double cpu = this->timer->getRunTime();
  if(LS)
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &Phi);
  else
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputPositionVectorToDisk()
{}
//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{//TODO: need to be re-written.
/*   //TODO: mmh not needed?
  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  if (_mmh)
    _mmh->resetOutputToStructure(postOp, *X, U);
*/
}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{
  // manually cast DistVec<int> to DistVec<double>. TODO: should avoid doing this.
  DistVec<double> nodeTagCopy(this->getVecInfo());
  if (TYPE==2) {
#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      Vec<double> &subTagCopy = nodeTagCopy(iSub);
      Vec<int> &subTag = nodeTag(iSub);
      for (int iNode=0; iNode<subTag.size(); iNode++) subTagCopy[iNode] = (double)subTag[iNode];
    }
  }

  if (TYPE==2) 
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, nodeTagCopy, *this->R, this->riemann, 0);
  else 
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, *this->R, this->riemann, 0);

  this->spaceOp->applyBCsToResidual(U, *this->R);
  double res = 0.0;
  res = this->spaceOp->computeRealFluidResidual(*this->R, *this->Rreal, *distLSS);
  return sqrt(res);
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));

  if (!this->problemType[ProblemData::UNSTEADY]) {
    double trhs = this->timer->getTimeSyncro();
    this->com->printf(2, "Getting residual norm\n");
    this->data->residual = computeResidualNorm(U);
    trhs = this->timer->getTimeSyncro() - trhs;
    if (it == 0)
      this->restart->residual = this->data->residual;
    if (this->data->resType == -1)
      this->com->printf(2, "Spatial residual norm = %.12e\n", this->data->residual);
    else
      this->com->printf(2, "Spatial residual norm[%d] = %.12e\n", this->data->resType, this->data->residual);
    this->com->printf(2, "Time for one residual evaluation: %f s\n", trhs);
  }

  this->com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateFSInterface()
{
  if (TYPE!=2) {
    fprintf(stderr,"In StructLevelSetTsDesc::updateFSInterface: Shouldn't call me! Abort.\n");
    exit(-1);
  }
  if (timeStep<=0.0) {fprintf(stderr,"in StructLevelSetTsDesc, timeStep <= 0.0! Abort...\n"); exit(-1);}
  double displacement = timeStep*fsiVelocity;
  fsiPosition += displacement*fsiNormal;
//  fsiPosition += 0.12; //for debug only.
  fprintf(stderr,"Interface: [%e %e %e], normal [%e %e %e].\n", fsiPosition[0], fsiPosition[1], fsiPosition[2], fsiNormal[0], fsiNormal[1], fsiNormal[2]);
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateNodeTag()
{
  if (TYPE!=2) {
    fprintf(stderr,"In StructLevelSetTsDesc::updateNodeTag: Shouldn't call me! Abort.\n");
    exit(-1);
  }
  if(LS) {fprintf(stderr,"Not ready for FF interface. Abort...\n"); exit(-1);}
  nodeTag0 = nodeTag;
  nodeTag = 0;

//  domain->updateNodeTag(*this->X, distLSS, nodeTag0, nodeTag);

  // temporarily.
#pragma omp parallel for
  for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    SVec<double,3> &subX = (*this->X)(iSub);
    Vec<int> &subNodeTag = nodeTag(iSub);
    for (int i=0; i<subX.size(); i++) {
      double scalar = fsiNormal[0]*(subX[i][0] - fsiPosition[0])
                    + fsiNormal[1]*(subX[i][1] - fsiPosition[1])
                    + fsiNormal[2]*(subX[i][2] - fsiPosition[2]);
      subNodeTag[i] = (scalar>0.0) ? -1 : 1;
    }
  }
}

//-------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::computeForceLoad()
{
  if (!Fs) {fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); return;}
  for (int i=0; i<numStructNodes; i++) Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
  this->spaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X, Fs, numStructNodes, distLSS, *Wstarij, *Wstarji);
/*
  double tempFs[numStructNodes*3];
  for (int i=0; i<numStructNodes; i++)
    for (int j=0; j<3; j++)
      tempFs[i*3+j] = Fs[i][j];
  this->com->globalSum(3*numStructNodes, tempFs);
  for (int i=0; i<numStructNodes; i++) {
    for (int j=0; j<3; j++)
      Fs[i][j] = tempFs[i*3+j];
  }

  FILE* fx = fopen("Fx.sol","w");
  FILE* fy = fopen("Fy.sol","w");
  FILE* fz = fopen("Fz.sol","w");
  for (int i=0; i<numStructNodes; i++)
    this->com->fprintf(fx, "%e\n", Fs[i][0]);
  for (int i=0; i<numStructNodes; i++)
    this->com->fprintf(fy, "%e\n", Fs[i][1]);
  for (int i=0; i<numStructNodes; i++)
    this->com->fprintf(fz, "%e\n", Fs[i][2]);
  fclose(fx);
  fclose(fy);
  fclose(fz);
*/
}

//-------------------------------------------------------------------------------

template <int dim>
void StructLevelSetTsDesc<dim>::getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]) {
  computeForceLoad();
  F[0] = F[1] = F[2] = 0.0;
  for (int i=0; i<numStructNodes; i++) {
    F[0]+=Fs[i][0]; F[1]+=Fs[i][1]; F[2]+=Fs[i][2];}

  M[0] = M[1] = M[2] = 0;
  Vec<Vec3D>& Xstruc = distLSS->getStructPosition();
  for (int i = 0; i < Xstruc.size(); ++i) {
     M[0] += Xstruc[i][1]*Fs[i][2]-Xstruc[i][2]*Fs[i][1];
     M[1] += Xstruc[i][2]*Fs[i][0]-Xstruc[i][0]*Fs[i][2];
     M[2] += Xstruc[i][0]*Fs[i][1]-Xstruc[i][1]*Fs[i][0];
  }
}

//-------------------------------------------------------------------------------

template <int dim>
void StructLevelSetTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft, DistSVec<double,dim> &U)
{
	  /*       // TODO: mmh not needed?
	  if (mmh) {
	    double work[2];
	    mmh->computeInterfaceWork(dt, postOp, geoState->getXn(), this->timeState->getUn(), *X, U, work);
	    this->restart->energy[0] += work[0];
	    restart->energy[1] += work[1];
	  }

	  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
	  if (_mmh)
	    _mmh->updateOutputToStructure(dt, dtLeft, postOp, *X, U);

	  if (hth)
	    hth->updateOutputToStructure(dt, dtLeft, postOp, *X, U);
	*/
  if(dynNodalTransfer) {
    computeForceLoad();
    // Now send the force to the structure
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->sendForce(v);
  }
}
