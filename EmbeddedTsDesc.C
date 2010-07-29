#include <EmbeddedTsDesc.h>
#include <DistExactRiemannSolver.h>
#include <FSI/DynamicNodalTransfer.h>

#ifdef DO_EMBEDDED
#include <IntersectorFRG/PhysBAMIntersect.h>
//#include <IntersectorPhysBAM/blablabla>
#endif

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
EmbeddedTsDesc<dim>::
EmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), nodeTag(this->getVecInfo()), nodeTag0(this->getVecInfo()),
  Vtemp(this->getVecInfo()), numFluid(ioData.eqs.numPhase)
{

  simType         = (ioData.problem.type[ProblemData::UNSTEADY]) ? 1 : 0;
  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;

  this->postOp->setForceGenerator(this);
  if (numFluid==1) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-STRUCTURE SIMULATION --------\n");
  if (numFluid>=2) this->com->fprintf(stderr,"-------- EMBEDDED FLUID-SHELL-FLUID SIMULATION --------\n");

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
  forceApp           = (ioData.embed.forceAlg==EmbeddedFramework::RECONSTRUCTED_SURFACE) ? 3 : 1;
  linRecAtInterface  = (ioData.embed.reconstruct==EmbeddedFramework::LINEAR) ? true : false;
  if (ioData.embed.riemannNormal!=EmbeddedFramework::AUTO)
    riemannNormal = (int)ioData.embed.riemannNormal;
  else //auto
    riemannNormal = simType ? 1 : 0;
      
  if(orderOfAccuracy==1) //first-order everywhere...
    linRecAtInterface = false; 

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  //for phase-change update
  Weights = 0;
  VWeights = 0;

//------------- For Fluid-Structure Interaction ---------------
  if(ioData.problem.type[ProblemData::FORCED] || ioData.problem.type[ProblemData::AERO]) {
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *this->domain->getCommunicator(), *this->domain->getStrCommunicator(),
                                                this->domain->getTimer());
    //for updating phase change
    Weights  = new DistVec<double>(this->getVecInfo());
    VWeights = new DistSVec<double,dim>(this->getVecInfo());
  } else
    dynNodalTransfer = 0;
//-------------------------------------------------------------

#ifdef DO_EMBEDDED
  switch (ioData.embed.intersectorName) {
    case EmbeddedFramework::FRG :
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();
        double (*xyz)[3] = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();
//        dynNodalTransfer->getEmbeddedMesh(nNodes,xyz,nElems,abc);
        distLSS = new DistPhysBAMIntersector(ioData, this->com, nNodes, xyz, nElems, abc);
      } else
        distLSS = new DistPhysBAMIntersector(ioData, this->com);
      break;
    case EmbeddedFramework::PHYSBAMLITE :  
      //distLSS = new ....
      this->com->fprintf(stderr,"ERROR: PhysBAM-Lite Intersector hasn't been integrated into the code.\n");
      exit(-1);
      break;
    default:
      this->com->fprintf(stderr,"ERROR: No valid intersector specified! Check input file\n");
      exit(-1);
  }
#else
  this->com->fprintf(stderr,"ERROR: Embedded framework is NOT compiled! Check your makefile.\n");
  exit(-1);
#endif

  //Riemann solution stored on edges
  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *Wstarij = 0.0;
  *Wstarji = 0.0;

  //cell-averaged structure normals
  if(riemannNormal==2) {
    Nsbar        = new DistSVec<double,3>(this->domain->getNodeDistInfo());
    *Nsbar = 0.0;
  } else 
    Nsbar = 0;

  //TODO: should be merged with fluidId in TsDesc
  nodeTag0 = 0;
  nodeTag = 0;

  //only for IncreasePressure
  Prate = ioData.mf.Prate;
  Pinit = ioData.mf.Pinit;
  Pscale = ioData.ref.rv.pressure;
  tmax = (ioData.bc.inlet.pressure - Pinit)/Prate;

  globIt = -1;
  inSubCycling = false;

  //store farfield state for phase-change update for fluid-fullbody
  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i];

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

  // Adam 04/06/2010
  ghostPoints = 0;
  switch (ioData.eqs.type)
    {
    case EquationsData::EULER:
      eqsType = EmbeddedTsDesc<dim>::EULER;
      break;
    case EquationsData::NAVIER_STOKES:
      eqsType = EmbeddedTsDesc<dim>::NAVIER_STOKES;
      ghostPoints  = new DistVec<GhostPoint<dim> *>(this->getVecInfo());
      break;
    }
}

//------------------------------------------------------------------------------

template<int dim>
EmbeddedTsDesc<dim>::~EmbeddedTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;

  if (dynNodalTransfer) delete dynNodalTransfer;
  if (Fs) delete[] Fs;
  if(ghostPoints) delete ghostPoints;
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  distLSS->initialize(this->domain,*this->X, ioData);
  if(riemannNormal==2)
    this->spaceOp->computeCellAveragedStructNormal(*Nsbar, distLSS);

  if(this->numFluid>1) //initialize nodeTags for multi-phase flow
    nodeTag0 = nodeTag = distLSS->getStatus();

  this->geoState->setup2(this->timeState->getData());
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), 
                         *U, ioData, &nodeTag);
  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->mmh);
  if(_mmh) {
    double *tMax = &(this->data)->maxTime;
    _mmh->setup(tMax); //obtain maxTime from structure
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
double EmbeddedTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U)
{
  if(!FsComputed&&simType) this->com->fprintf(stderr,"WARNING: FSI force not computed!\n");
  FsComputed = false; //reset FsComputed at the beginning of a fluid iteration

  //check if it's in subcycling with iCycle>1.
  if(globIt==it)
    inSubCycling = true;
  else {
    globIt = it;
    inSubCycling = false;
  }

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
void EmbeddedTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U); 
}

//-----------------------------------------------------------------------------

template<int dim>
int EmbeddedTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
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
void EmbeddedTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt, int it, double t,
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
    if (numFluid>=2) {
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
    if (numFluid>=2){
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
void EmbeddedTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                                             double t, double dt, DistSVec<double,dim> &U)
{

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;

  if (numFluid>=2) {
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

  TsRestart *restart2 = this->restart; // Bug: compiler does not accept this->restart->writeToDisk<dim,1>(...)
                                       //      it does not seem to understand the template
  restart2->writeToDisk<dim,1>(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState);
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
void EmbeddedTsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
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
void EmbeddedTsDesc<dim>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{  
  // Ghost-Points Population
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->nullifyPointers();
      this->spaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,this->nodeTag);
    }
  if(this->numFluid==1)
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface, *this->R, this->riemann, riemannNormal, Nsbar, 0, ghostPoints);
  else //numFluid>1
    this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface,  nodeTag, *this->R, this->riemann, riemannNormal, Nsbar, 0);

  this->spaceOp->applyBCsToResidual(U, *this->R);

  double res = 0.0;
  if(this->numFluid==1)
    res = this->spaceOp->computeRealFluidResidual(*this->R, *this->Rreal, *distLSS);

  return sqrt(res);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
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
void EmbeddedTsDesc<dim>::computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji)
{
  if (!Fs) {fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); return;}
  for (int i=0; i<numStructNodes; i++) 
    Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
  this->spaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X, Fs, numStructNodes, distLSS, *Wij, *Wji);
  //at this stage Fs is NOT globally assembled!
}

//-------------------------------------------------------------------------------

template <int dim>
void EmbeddedTsDesc<dim>::getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           double F[3], double M[3]) 
{
  if (!FsComputed) 
    computeForceLoad(this->Wstarij, this->Wstarji);

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
void EmbeddedTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft, DistSVec<double,dim> &U)
{
  if(dynNodalTransfer) {
    computeForceLoad(this->Wstarij, this->Wstarji);
    FsComputed = true; //to avoid redundant computation of Fs.

    // Now "accumulate" the force for the embedded structure
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(dt, dtLeft, v);
  }
}

//-------------------------------------------------------------------------------

template<int dim>
bool EmbeddedTsDesc<dim>::IncreasePressure(double dt, double t, DistSVec<double,dim> &U)
{

  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) {// max pressure was reached, so now we solve
//    this->com->fprintf(stdout, "max pressure reached\n"); 
    return true;
  } 

  // max pressure not reached, so we do not solve. Instead, we just increase pressure and let structure react
  
  // construct Wij, Wji from U. 
  DistSVec<double,dim> VV(this->getVecInfo());
  this->varFcn->conservativeToPrimitive(U,VV,&nodeTag);
  SubDomain **subD = this->domain->getSubDomain();

#pragma omp parallel for
  for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    SVec<double,dim> &subWstarij = (*Wstarij)(iSub);
    SVec<double,dim> &subWstarji = (*Wstarji)(iSub);
    SVec<double,dim> &subVV = VV(iSub);
    int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();

    for (int l=0; l<subWstarij.size(); l++) {
      int i = ptr[l][0];
      int j = ptr[l][1];
      for (int k=0; k<dim; k++) {
        subWstarij[l][k] = subVV[i][k];
        subWstarji[l][k] = subVV[j][k];
      }
    }
  }

  if(this->mmh && !inSubCycling) {
    double tw = this->timer->getTime();
    //store previous states for phase-change update
    if(this->numFluid==1)
      this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                     *this->VWeights, this->distLSS);
    else //numFluid>1
      this->spaceOp->computeWeightsForEmbeddedStruct(*this->X, U, this->Vtemp, *this->Weights,
                                                         *this->VWeights, this->distLSS, &this->nodeTag);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);
    //recompute intersections
    tw = this->timer->getTime();
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts);
    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
    //update nodeTags (only for numFluid>1)
    if(numFluid>1) {
      nodeTag0 = this->nodeTag;
      nodeTag = this->distLSS->getStatus();
    }

    //update phase-change
    tw = this->timer->getTime();
    if(this->numFluid==1)
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, vfar);
    else //numFluid>1
      this->spaceOp->updatePhaseChange(this->Vtemp, U, this->Weights, this->VWeights, this->distLSS, vfar, &this->nodeTag);
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  } 
  
  this->com->fprintf(stdout, "about to increase pressure to %e\n", (Pinit+t*Prate)*Pscale);
  this->domain->IncreasePressure(Pinit+t*Prate, this->varFcn, U, nodeTag);

  return false;

}

//------------------------------------------------------------------------------










