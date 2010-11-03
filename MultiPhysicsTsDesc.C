#include <DistExactRiemannSolver.h>
#include <Domain.h>
#include <FluidSelector.h>
#include <FSI/DynamicNodalTransfer.h>
#include <GeoSource.h>
#include <LevelSet.h>
#include <MultiPhysicsTsDesc.h>

#ifdef DO_EMBEDDED
#include <IntersectorFRG/IntersectorFRG.h>
#include <IntersectorPhysBAM/IntersectorPhysBAM.h>
#endif

#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhysicsTsDesc<dim,dimLS>::
MultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), V0(this->getVecInfo()),
  PhiV(this->getVecInfo()), PhiWeights(this->getVecInfo()), 
  fluidSelector(ioData.eqs.numPhase, ioData, dom), //memory allocated for fluidIds
  Vtemp(this->getVecInfo()), numFluid(ioData.eqs.numPhase)
{
  simType = (ioData.problem.type[ProblemData::UNSTEADY]) ? 1 : 0;
  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;

  multiPhaseSpaceOp = new MultiPhaseSpaceOperator<dim,dimLS>(ioData, this->varFcn, this->bcData, this->geoState, 
                                                             this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, multiPhaseSpaceOp, this->varFcn, this->domain, this->V);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  setupEmbeddedFSISolver(ioData);
  setupMultiPhaseFlowSolver(ioData);

  // only for IncreasePressure
  Prate = ioData.mf.Prate;
  Pinit = ioData.mf.Pinit;
  Pscale = ioData.ref.rv.pressure;
  tmax = (ioData.bc.inlet.pressure - Pinit)/Prate;

  // miscellaneous
  globIt = -1;
  inSubCycling = false;
  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i]; //for phase-change update only
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupEmbeddedFSISolver(IoData &ioData)
{

  this->postOp->setForceGenerator(this);

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
  forceApp           = (ioData.embed.forceAlg==EmbeddedFramework::RECONSTRUCTED_SURFACE) ? 3 : 1;
  linRecAtInterface  = (ioData.embed.reconstruct==EmbeddedFramework::LINEAR) ? true : false;
  if (ioData.embed.riemannNormal!=EmbeddedFramework::AUTO)
    riemannNormal = (int)ioData.embed.riemannNormal;
  else //auto
    riemannNormal = simType ? 1 : 0;

  if(orderOfAccuracy==1) //first-order everywhere...
    linRecAtInterface = false;

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

//-------------------- Setup Intersector ----------------------
#ifdef DO_EMBEDDED
  switch (ioData.embed.intersectorName) {
    case EmbeddedFramework::FRG : // use IntersectorFRG
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
        this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();
        double (*xyz)[3] = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();
        distLSS = new DistIntersectorFRG(ioData, this->com, nNodes, xyz, nElems, abc);
      } else
        distLSS = new DistIntersectorFRG(ioData, this->com);
      break;
    case EmbeddedFramework::PHYSBAM : // use IntersectorPhysBAM
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
        this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();
        double (*xyz)[3] = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();
        distLSS = new DistIntersectorPhysBAM(ioData, this->com, nNodes, xyz, nElems, abc);
      } else
        distLSS = new DistIntersectorPhysBAM(ioData, this->com);
      break;
    default: // use exit(-1)... 
      this->com->fprintf(stderr,"ERROR: No valid intersector specified! Check input file\n");
      exit(-1);
  }
#else
  this->com->fprintf(stderr,"ERROR: Embedded framework is NOT compiled! Check your makefile.\n");
  exit(-1);
#endif
//-------------------------------------------------------------

  //Riemann solution stored on edges
  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *Wstarij = 0.0;
  *Wstarji = 0.0;

  //cell-averaged structure normals
  if(riemannNormal==2) {
    Nsbar = new DistSVec<double,3>(this->domain->getNodeDistInfo());
    *Nsbar = 0.0;
  } else
    Nsbar = 0;

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

//------ Setup ghostPoints for Viscous Flows ------------------
  ghostPoints = 0;
  switch (ioData.eqs.type)
    {
    case EquationsData::EULER:
      eqsType = MultiPhysicsTsDesc<dim,dimLS>::EULER;
      break;
    case EquationsData::NAVIER_STOKES:
      eqsType = MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES;
      ghostPoints  = new DistVec<GhostPoint<dim> *>(this->getVecInfo());
      ghostPoints->nullifyPointers();
      break;
    }
//-------------------------------------------------------------

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupMultiPhaseFlowSolver(IoData &ioData)
{
  LS = new LevelSet<dimLS>(ioData, this->domain);
  frequencyLS = ioData.mf.frequency;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhysicsTsDesc<dim,dimLS>::~MultiPhysicsTsDesc()
{
  if (riemann) delete riemann;

  // Embedded FSI
  if (distLSS) delete distLSS;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;
  if (dynNodalTransfer) delete dynNodalTransfer;
  if (Fs) delete[] Fs;
  if(ghostPoints) {
    ghostPoints->deletePointers();
    delete ghostPoints;
  }

  // Level-set
  if (LS) delete LS;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  // Setup fluid mesh geometry
  this->geoState->setup2(this->timeState->getData());
  // Initialize intersector and compute intersections
  DistVec<int> point_based_id(this->domain->getNodeDistInfo());
  distLSS->initialize(this->domain,*this->X, ioData, &point_based_id);
  if(riemannNormal==2){
    this->multiPhaseSpaceOp->computeCellAveragedStructNormal(*Nsbar, distLSS);
  }
  // Initialize fluid state vector
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(),
                         *U, ioData, &point_based_id); //populate U by i.c. or restart data.
  // Initialize level-sets 
  LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData);

  // Initialize or reinitialize (i.e. at the beginning of a restart) fluid Ids
  if(this->input->levelsets[0] == 0) // init
    fluidSelector.initializeFluidIds(distLSS->getStatus(), *this->X, ioData);
  else //restart
    fluidSelector.reinitializeFluidIds(distLSS->getStatus(), Phi);

  // Initialize the embedded FSI handler
  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->mmh);
  if(_mmh) {
    double *tMax = &(this->data)->maxTime;
    _mmh->setup(tMax); //obtain maxTime from structure
  }

  //compute fluid indueced force on the structural surface
  // construct Wij, Wji from U. Then use them for force calculation.
  DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> VV(this->getVecInfo());
  *Wij = 0.0;
  *Wji = 0.0;
  this->varFcn->conservativeToPrimitive(*U,VV,fluidSelector.fluidId);
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

template<int dim, int dimLS>
double MultiPhysicsTsDesc<dim,dimLS>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U)
{
  if(!FsComputed&&dynNodalTransfer) this->com->fprintf(stderr,"WARNING: FSI force not computed!\n");
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
  dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                          &numSubCycles, *this->geoState, *this->A, U, *(fluidSelector.fluidId));

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  dtf = dt;
  dtfLeft = *dtLeft + dt;

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;
}

//---------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);

  LS->update(Phi);
  fluidSelector.update();

  this->timeState->update(U, *(fluidSelector.fluidIdn), fluidSelector.fluidIdnm1, riemann);
                            //fluidIdn, fluidIdnm1 and riemann are used only for implicit time-integrators
  if(frequencyLS > 0 && it%frequencyLS == 0){
    this->com->printf(5, "LevelSet norm before reinitialization = %e\n", Phi.norm());
    LS->conservativeToPrimitive(Phi,PhiV,U);
    LS->reinitializeLevelSet(*this->geoState,*this->X, *this->A, U, PhiV);
    LS->primitiveToConservative(PhiV,Phi,U);
    this->com->printf(5, "LevelSet norm after reinitialization = %e\n", Phi.norm());
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int MultiPhysicsTsDesc<dim,dimLS>::checkSolution(DistSVec<double,dim> &U)
{
  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, *fluidSelector.fluidId, *fluidSelector.fluidIdn);
                             // fluidIdn is only used for screen output when an error is found
  return ierr;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                          int it, double t, DistSVec<double,dim> &U)
{
//  conservationErrors(U,it); /* not checking conservation errors for the moment...*/

  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U); // Phi?

//  this->output->setMeshMotionHandler(ioData, this->mmh); /*currently not supported by embedded framework.*/

  this->output->openAsciiFiles();
  this->timer->setSetupTime();

  if (it == 0) {
    this->multiPhaseSpaceOp->computeGradP(*this->X, *this->A, U); /*really used???*/
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId, &Phi);
//    this->output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
//    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{
  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;

  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId, &Phi);
//  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);

  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);
  this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, distLSS->getStructPosition_n());

  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)
{
  double cpu = this->timer->getRunTime();
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{/* currently only consider the embedded framework. No mesh motion. */}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  /* currently only consider the embedded framework. No mesh motion. */
/*  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, fluidSelector.fluidId);
*/
}

//------------------------------------------------------------------------------

template <int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::updateOutputToStructure(double dt, double dtLeft, DistSVec<double,dim> &U)
{
  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(),
                                    this->timeState->getUn(), *this->X, U,
                                    work, fluidSelector.fluidIdn, fluidSelector.fluidId); //FF interface
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  if(dynNodalTransfer) {
    computeForceLoad(this->Wstarij, this->Wstarji);
    FsComputed = true; //to avoid redundant computation of Fs.

    // Now "accumulate" the force for the embedded structure
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(dt, dtLeft, v);
  }
}

//-------------------------------------------------------------------------------

template<int dim, int dimLS>
double MultiPhysicsTsDesc<dim,dimLS>::computeResidualNorm(DistSVec<double,dim>& U)
{
  // Ghost-Points Population (for Navier-Stokes only)
  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES) {
    this->ghostPoints->deletePointers();
    this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,*(this->fluidSelector.fluidId));
  }
  LS->conservativeToPrimitive(Phi,PhiV,U);
  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface, this->riemann, riemannNormal, Nsbar, PhiV, fluidSelector, *this->R, 0, 0);

  this->multiPhaseSpaceOp->applyBCsToResidual(U, *this->R);

  double res = 0.0;
  if(this->numFluid==1)
    res = this->multiPhaseSpaceOp->computeRealFluidResidual(*this->R, *this->Rreal, *distLSS);

  return sqrt(res);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::monitorInitialState(int it, DistSVec<double,dim> &U)
{
  /* only used for steady simulations. Is it meaningful when a ff interface is present */
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

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji)
{
  if (!Fs) {fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); return;}
  double t0 = this->timer->getTime();
  for (int i=0; i<numStructNodes; i++)
    Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
  this->multiPhaseSpaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X, Fs, numStructNodes, distLSS, *Wij, *Wji);
  this->timer->addEmbeddedForceTime(t0);
  //at this stage Fs is NOT globally assembled!
}

//-------------------------------------------------------------------------------

template <int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::getForcesAndMoments(DistSVec<double,dim> &U, DistSVec<double,3> &X,
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
template<int dim, int dimLS>
bool MultiPhysicsTsDesc<dim,dimLS>::IncreasePressure(double dt, double t, DistSVec<double,dim> &U)
{ return true; }
//-------------------------------------------------------------------------------









