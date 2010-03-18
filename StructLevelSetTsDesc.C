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
  Phi(this->getVecInfo()), Vg(this->getVecInfo()), Vtemp(this->getVecInfo()),
  PhiV(this->getVecInfo()), boundaryFlux(this->getVecInfo()),
  computedQty(this->getVecInfo()), interfaceFlux(this->getVecInfo()),numFluid(ioData.eqs.numPhase)
{

  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;
  forceApp = ioData.strucIntersect.forceApproach;
  pressureChoice = ioData.strucIntersect.pressureChoice;
  phaseChangeChoice = ioData.strucIntersect.phaseChangeChoice;

  this->postOp->setForceGenerator(this);
  if (numFluid==1) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-STRUCTURE SIMULATION --------\n");
  if (numFluid==2) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-SHELL-FLUID SIMULATION --------\n");

  interpolatedNormal = (ioData.strucIntersect.normal==StructureIntersect::INTERPOLATED) ? true : false;
  linRecAtInterface  = (ioData.strucIntersect.reconstruct==StructureIntersect::LINEAR) ? true : false;

  if (ioData.strucIntersect.riemannNormal==StructureIntersect::FLUID) 
    riemannNormal = 1;
  else if (ioData.strucIntersect.riemannNormal==StructureIntersect::STRUCTURE)
    riemannNormal = 0;
  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);
  distLSS = 0;

  const char *intersectorName = ioData.strucIntersect.intersectorName;
  if(intersectorName != 0)
    distLSS = IntersectionFactory::getIntersectionObject(intersectorName, *this->domain);
  distLSS->setNumOfFluids(numFluid);

  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *Wstarij = 0.0;
  *Wstarji = 0.0;
  nodeTag0 = 0;
  nodeTag = 0;

  Weights = 0;
  VWeights = 0;

//------ load structure mesh information ----------------------
  numStructNodes = 0;
  Fs = 0;
  numStructNodes = distLSS->getNumStructNodes();
  if (numStructNodes>0) {
    this->com->fprintf(stderr,"# of struct nodes: %d.\n", numStructNodes);
    // We allocate Fs from memory that allows fast one-sided MPI communication
    Fs = new (*this->com) double[numStructNodes][3];
  } else 
    this->com->fprintf(stderr,"Warning: failed loading structure mesh information!\n");

  FsComputed = false;
//-------------------------------------------------------------

// for FF interface
  //LS = 0; // for debug!
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

//------------- For Fluid-Structure Interaction -------------------------
  if(ioData.embeddedStructure.type != EmbeddedStructureInfo::NONE) {
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *this->domain->getCommunicator(), *this->domain->getStrCommunicator(),
                                                this->domain->getTimer());

    //for updating phase change
    Weights  = new DistVec<double>(this->getVecInfo());
    VWeights = new DistSVec<double,dim>(this->getVecInfo());
  } else
    dynNodalTransfer = 0;

//----------------------------------------------------------------------

}

//------------------------------------------------------------------------------

template<int dim>
StructLevelSetTsDesc<dim>::~StructLevelSetTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;

  if (dynNodalTransfer) delete dynNodalTransfer;
  if (Fs) delete[] Fs;

  //if (LS) delete LS;
  if (Vgf) delete Vgf;
  if (Vgfweight) delete Vgfweight;
  if(tmpDistSVec)  delete tmpDistSVec;
  if(tmpDistSVec2) delete tmpDistSVec2;
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  distLSS->initialize(this->domain,*this->X, interpolatedNormal);
  if(this->numFluid>1) //initialize nodeTags for multi-phase flow
    nodeTag0 = nodeTag = distLSS->getStatus();

  this->geoState->setup2(this->timeState->getData());
  this->timeState->setup(this->input->solutions, this->bcData->getInletBoundaryVector(), 
                         *this->X, *U, &ioData, &nodeTag);
  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->mmh);
  if(_mmh) {
    double *tmax = &(this->data)->maxTime;
    _mmh->setup(tmax); //obtain maxTime from structure
  }

  //compute force
  // construct Wij, Wji from U. Then use them for force calculation.
  DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> VV(this->getVecInfo());
  *Wij = 0.0;
  *Wji = 0.0;
  this->varFcn->conservativeToPrimitive(*U,VV,&nodeTag);
  SubDomain **subD = this->domain->getSubDomain();

#pragma omp parallel for
  for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    SVec<double,dim> &subWij = (*Wij)(iSub);
    SVec<double,dim> &subWji = (*Wji)(iSub);
    SVec<double,dim> &subWstarij = (*Wstarij)(iSub);
    SVec<double,dim> &subWstarji = (*Wstarji)(iSub);
    SVec<double,dim> &subVV = VV(iSub);
    int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();
    if (subWij.size()!=(subD[iSub]->getEdges()).size()) {fprintf(stderr,"WRONG!!!\n"); exit(-1);}

    for (int l=0; l<subWij.size(); l++) {
      int i = ptr[l][0];
      int j = ptr[l][1];
      for (int k=0; k<dim; k++) {
        subWij[l][k] = subVV[i][k];
        subWji[l][k] = subVV[j][k];
      }
    }
  }

  computeForceLoad(Wij, Wji);
  delete Wij;
  delete Wji;
  FsComputed = true;
  // Now "accumulate" the force for the embedded structure
  if(dynNodalTransfer){
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(0.0, 0.0, v); //dt=dtLeft=0.0-->They are not used!
  }
}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U)
{
  if(!FsComputed) fprintf(stderr,"WARNING: FSI force not computed!\n");
  FsComputed = false; //reset FsComputed at the beginning of a fluid iteration

  double t0 = this->timer->getTime();
  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);
  int numSubCycles = 1;

  double dt;
  if(numFluid==1)
    dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->X, *this->A, U);
  else //numFLuid>1
    dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->A, U, nodeTag);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  dtf = dt;
  dtfLeft = *dtLeft + dt;

  return dt;
}

//---------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U); 
}

//-----------------------------------------------------------------------------

template<int dim>
int StructLevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{
  int ierr = 0;
  if(numFluid==1)
    ierr = this->domain->checkSolution(this->varFcn, U); //also check ghost nodes.
  else
    ierr = this->domain->checkSolution(this->varFcn, U, nodeTag);

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
{
  if (it == this->data->maxIts)
    *lastIt = true;
  else 
    monitorInitialState(it, U);

  this->output->openAsciiFiles();
  this->timer->setSetupTime();

  if (it == 0) {
    // First time step: compute GradP before computing forces
    this->spaceOp->computeGradP(*this->X, *this->A, U);
    if (numFluid==2) {
      this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
      this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
      this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
      this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    } else {
      this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
      this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    }
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    if (numFluid==2){
      DistSVec<double,1> tempPhi(this->domain->getNodeDistInfo());
      tempPhi = distLSS->getPhi();
      this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, tempPhi, nodeTag);
      //this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, distLSS->getPhi(), nodeTag);
    }
    else 
      this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                                             double t, double dt, DistSVec<double,dim> &U)
{

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;

  if (numFluid==2) {
    this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
    DistSVec<double,1> tempPhi(this->domain->getNodeDistInfo());
    tempPhi = distLSS->getPhi();
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, tempPhi, nodeTag);
    //this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, distLSS->getPhi(), nodeTag);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  } 
  else { //numFluid == 1

    this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
    this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

  // TODO
  //this->restart->writeToDisk<dim,1>(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState);
  this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, this->distLSS->getStructPosition_n());

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
  if(this->numFluid==1)
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  else
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &this->nodeTag);
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{
  if(this->numFluid==1)
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface, *this->R, this->riemann, riemannNormal, 0);
  else //numFluid>1
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface,  nodeTag, *this->R, this->riemann, riemannNormal, 0);

  this->spaceOp->applyBCsToResidual(U, *this->R);

  double res = 0.0;
  if(this->numFluid==1)
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
void StructLevelSetTsDesc<dim>::computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji)
{
  if (!Fs) {fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); return;}
  for (int i=0; i<numStructNodes; i++) 
    Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
  this->spaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X, Fs, numStructNodes, distLSS, *Wij, *Wji);
  //at this stage Fs is NOT globally assembled!
}

//-------------------------------------------------------------------------------

template <int dim>
void StructLevelSetTsDesc<dim>::getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]) 
{
  if (!FsComputed) 
    if (pressureChoice==0) //(default) use p* in force calculation.
      computeForceLoad(this->Wstarij, this->Wstarji);

    else if (pressureChoice==1) {
      // construct Wij, Wji from U. Then use them for force calculation. 
      DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
      DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
      DistSVec<double,dim> VV(this->getVecInfo());
      *Wij = 0.0;
      *Wji = 0.0;
      this->varFcn->conservativeToPrimitive(U,VV,&nodeTag);
      SubDomain **subD = this->domain->getSubDomain();

#pragma omp parallel for
      for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
        SVec<double,dim> &subWij = (*Wij)(iSub);
        SVec<double,dim> &subWji = (*Wji)(iSub);
        SVec<double,dim> &subWstarij = (*Wstarij)(iSub);
        SVec<double,dim> &subWstarji = (*Wstarji)(iSub);
        SVec<double,dim> &subVV = VV(iSub);
        int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();    
        if (subWij.size()!=(subD[iSub]->getEdges()).size()) {fprintf(stderr,"WRONG!!!\n"); exit(-1);}
  
        for (int l=0; l<subWij.size(); l++) {
          int i = ptr[l][0];
          int j = ptr[l][1];
          if (subWstarij[l][0]>1e-10 && subWstarij[l][4]>1e-10)  
            for (int k=0; k<dim; k++) subWij[l][k] = subVV[i][k];
          if (subWstarji[l][0]>1e-10 && subWstarji[l][4]>1e-10)  
            for (int k=0; k<dim; k++) subWji[l][k] = subVV[j][k];
        }
      } 

      computeForceLoad(Wij, Wji);
      delete Wij;
      delete Wji;
    }

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
  if(dynNodalTransfer) {

    if (pressureChoice==0) //(default) use p* in force calculation.
      computeForceLoad(this->Wstarij, this->Wstarji);
    
    else if (pressureChoice==1) {
      // construct Wij, Wji from U. Then use them for force calculation. 
      DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
      DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
      DistSVec<double,dim> VV(this->getVecInfo());
      *Wij = 0.0;
      *Wji = 0.0;
      this->varFcn->conservativeToPrimitive(U,VV,&nodeTag);
      SubDomain **subD = this->domain->getSubDomain();

#pragma omp parallel for
      for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
        SVec<double,dim> &subWij = (*Wij)(iSub);
        SVec<double,dim> &subWji = (*Wji)(iSub);
        SVec<double,dim> &subWstarij = (*Wstarij)(iSub);
        SVec<double,dim> &subWstarji = (*Wstarji)(iSub);
        SVec<double,dim> &subVV = VV(iSub);
        int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();
        if (subWij.size()!=(subD[iSub]->getEdges()).size()) {fprintf(stderr,"WRONG!!!\n"); exit(-1);}
  
        for (int l=0; l<subWij.size(); l++) {
          int i = ptr[l][0];
          int j = ptr[l][1];
          if (subWstarij[l][0]>1e-10)  
            for (int k=0; k<dim; k++) subWij[l][k] = subVV[i][k];
          if (subWstarji[l][0]>1e-10)  
            for (int k=0; k<dim; k++) subWji[l][k] = subVV[j][k];
        }
      }

      computeForceLoad(Wij, Wji);
      delete Wij;
      delete Wji;

    }

    FsComputed = true; //to avoid redundant computation of Fs.
    // Now "accumulate" the force for the embedded structure
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(dt, dtLeft, v);
  }
}

//-------------------------------------------------------------------------------











