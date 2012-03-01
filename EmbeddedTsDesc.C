#include <EmbeddedTsDesc.h>
#include <DistExactRiemannSolver.h>
#include <FSI/DynamicNodalTransfer.h>
#include <FSI/CrackingSurface.h>
#include <Domain.h>

#ifdef DO_EMBEDDED
#include <IntersectorFRG/IntersectorFRG.h>
#include <IntersectorPhysBAM/IntersectorPhysBAM.h>
#endif

#include <cmath>

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
  Vtemp(this->getVecInfo()), numFluid(ioData.eqs.numPhase), Wtemp(this->getVecInfo())
{

  simType         = (ioData.problem.type[ProblemData::UNSTEADY]) ? 1 : 0;
  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;

  this->postOp->setForceGenerator(this);

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
  forceApp           = (ioData.embed.forceAlg==EmbeddedFramework::RECONSTRUCTED_SURFACE) ? 3 : 1;


  // Debug - To be deleted
  std::ifstream forceCalculationType("forceCalculationType.txt");
  if(forceCalculationType) {
    int newVersion;
    forceCalculationType>>newVersion;
    if(newVersion) forceApp++;
  }
//  this->com->fprintf(stderr,"*************************************** ForceApproach: %d *************************************\n",forceApp);


  linRecAtInterface  = (ioData.embed.reconstruct==EmbeddedFramework::LINEAR) ? true : false;
  riemannNormal = (int)ioData.embed.riemannNormal;
      
  if(orderOfAccuracy==1) //first-order everywhere...
    linRecAtInterface = false; 

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  //for phase-change update
  Weights = 0;
  VWeights = 0;

//------------- For Fluid-Structure Interaction ---------------
  withCracking = false;
  if(ioData.problem.type[ProblemData::FORCED] || ioData.problem.type[ProblemData::AERO]) {
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *this->domain->getCommunicator(), *this->domain->getStrCommunicator(),
                                                this->domain->getTimer());
    withCracking = dynNodalTransfer->cracking(); 

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
        //this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();

        if(withCracking) {
          this->com->fprintf(stderr,"ERROR: IntersectorFRG is not capable of handling cracking structures!\n");
          this->com->fprintf(stderr,"       Try IntersectorPhysBAM instead!\n");
          exit(-1);
        }

        double *xyz = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();
        distLSS = new DistIntersectorFRG(ioData, this->com, nNodes, xyz, nElems, abc);
      } else
        distLSS = new DistIntersectorFRG(ioData, this->com);
      break;

    case EmbeddedFramework::PHYSBAM : 
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
        //this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();
        double *xyz = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();

        if(withCracking) {
          this->com->fprintf(stderr,"Note: Topology change of the embedded surface will be considered...\n"); 
          distLSS = new DistIntersectorPhysBAM(ioData, this->com, nNodes, xyz, nElems, abc, dynNodalTransfer->getCrackingSurface());
        } else
          distLSS = new DistIntersectorPhysBAM(ioData, this->com, nNodes, xyz, nElems, abc);

      } else
        distLSS = new DistIntersectorPhysBAM(ioData, this->com);
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

  //copies for fail safe
  WstarijCopy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  WstarjiCopy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *WstarijCopy = 0.0;
  *WstarjiCopy = 0.0;



  if (this->timeState->useNm1()) {
    Wstarij_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    Wstarji_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    *Wstarij_nm1 = 0.0;
    *Wstarji_nm1 = 0.0;

    Wstarij_nm1Copy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    Wstarji_nm1Copy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    *Wstarij_nm1Copy = 0.0;
    *Wstarji_nm1Copy = 0.0;
  } else {
    Wstarij_nm1 = 0;
    Wstarji_nm1 = 0;

    Wstarij_nm1Copy = 0;
    Wstarji_nm1Copy = 0;
  }
  UCopy = new DistSVec<double,dim>(this->domain->getNodeDistInfo());

  //cell-averaged structure normals
  if(riemannNormal==2) {
    Nsbar        = new DistSVec<double,3>(this->domain->getNodeDistInfo());
    *Nsbar = 0.0;
  } else 
    Nsbar = 0;

  //TODO: should be merged with fluidId in TsDesc
  nodeTag0 = 0;
  nodeTag = 0;

 nodeTagCopy = new DistVec<int>(this->getVecInfo());
 *nodeTagCopy = 0;
 nodeTag0Copy = new DistVec<int>(this->getVecInfo());
 *nodeTag0Copy = 0;


//---------------------------------------------------------------------
  // for IncreasePressure
  if(ioData.implosion.type==ImplosionSetup::LINEAR)
    implosionSetupType = EmbeddedTsDesc<dim>::LINEAR;
  else if(ioData.implosion.type==ImplosionSetup::SMOOTHSTEP)
    implosionSetupType = EmbeddedTsDesc<dim>::SMOOTHSTEP;
  else {
    this->com->fprintf(stderr,"ERROR: ImplosionSetup::Type = %d. Code not supported!\n", ioData.implosion.type);
    implosionSetupType = EmbeddedTsDesc<dim>::NONE;
    exit(-1);
  }

  Prate = ioData.implosion.Prate;
  Pinit = ioData.implosion.Pinit;
  Pfinal = ioData.bc.inlet.pressure;

  Pscale = ioData.ref.rv.pressure;
  intersector_freq = ioData.implosion.intersector_freq;
  if(ioData.implosion.type==ImplosionSetup::LINEAR)
    tmax = (ioData.bc.inlet.pressure - Pinit)/Prate;
  else if(ioData.implosion.type==ImplosionSetup::SMOOTHSTEP) {
    tmax = ioData.implosion.tmax;
    if(tmax<=0.0) {
      this->com->fprintf(stderr,"ERROR: ImplosionSetup::Tmax = %e!\n", tmax);
      exit(-1);
    }
  }

  if(intersector_freq<1) {
    this->com->fprintf(stderr,"ERROR: InterfaceTrackingFrequency must be larger than 0. Currently it is %d.\n", intersector_freq);
    exit(-1);
  }

  recomputeIntersections = true;
  unifPressure[0] = unifPressure[1] = Pinit;
//---------------------------------------------------------------------

  globIt = -1;
  inSubCycling = false;

  //store farfield state for phase-change update for fluid-fullbody
  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i];

//------ load structure mesh information ----------------------
  Fs = 0;
  numStructNodes = distLSS->getNumStructNodes();
  totStructNodes = dynNodalTransfer ? dynNodalTransfer->totStNodes() : numStructNodes;

  if (numStructNodes>0) {
    this->com->fprintf(stderr,"- Embedded Structure Surface: %d (%d) nodes\n", numStructNodes, totStructNodes);
    // We allocate Fs from memory that allows fast one-sided MPI communication
    Fs = new (*this->com) double[totStructNodes][3];
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
      ghostPoints->nullifyPointers();
      break;
    }

  increasingPressure = false;

}


//------------------------------------------------------------------------------

template<int dim>
EmbeddedTsDesc<dim>::~EmbeddedTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (Wstarij_nm1) delete Wstarij_nm1;
  if (Wstarji_nm1) delete Wstarji_nm1;
  if (WstarijCopy) delete WstarijCopy;
  if (WstarjiCopy) delete WstarjiCopy;
  if (Wstarij_nm1Copy) delete Wstarij_nm1Copy;
  if (Wstarji_nm1Copy) delete Wstarji_nm1Copy;
  if (UCopy) delete UCopy;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;

  if (dynNodalTransfer) delete dynNodalTransfer;
  //PJSA if (Fs) delete[] Fs;
  if(ghostPoints) 
    {
      ghostPoints->deletePointers();
      delete ghostPoints;
    }
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  // Setup fluid mesh geometry
  this->geoState->setup2(this->timeState->getData());
  // Initialize intersector and compute intersections
  DistVec<int> point_based_id(this->domain->getNodeDistInfo());
  distLSS->initialize(this->domain,*this->X, ioData, &point_based_id);
  if(riemannNormal==2){
    this->spaceOp->computeCellAveragedStructNormal(*Nsbar, distLSS);
  }
  // Initialize fluid state vector
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(),
                         *U, ioData, &point_based_id); //populate U by i.c. or restart data.
  // Initialize fluid Ids
  nodeTag0 = nodeTag = distLSS->getStatus();
  // Initialize the embedded FSI handler
  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->mmh);
  if(_mmh) {
    double *tMax = &(this->data)->maxTime;
    _mmh->setup(tMax); //obtain maxTime from structure
  }

  // If 'IncreasePressure' is activated, re-initialize the fluid state
  if(Pinit>=0.0 && Prate>=0.0 && this->getInitialTime()<tmax) {
    increasingPressure = true;
    this->domain->IncreasePressure(currentPressure(this->getInitialTime()), this->varFcn, *U, nodeTag);
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

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
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


  // Ghost-Points Population
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers(); // Not needed cause it has already been done in the constructor.
      this->spaceOp->populateGhostPoints(this->ghostPoints,*U,this->varFcn,this->distLSS,this->nodeTag);
    }
  // Population of spaceOp->V for the force computation
  this->spaceOp->conservativeToPrimitive(*U);


  computeForceLoad(Wij, Wji);
  delete Wij;
  delete Wji;


  FsComputed = true;
  // Now "accumulate" the force for the embedded structure
  if(dynNodalTransfer){
    numStructNodes = dynNodalTransfer->numStNodes();
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(0.0, 0.0, v); //dt=dtLeft=0.0-->They are not used!
  }
}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
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

  this->com->barrier();
  double t0 = this->timer->getTime();
  int numSubCycles = 1;
  double dt=0.0;

  if(TsDesc<dim>::failSafeFlag == false){
    if(TsDesc<dim>::timeStepCalculation == TsData::CFL || it==1){
      this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);
      if(numFluid==1)
        dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                                &numSubCycles, *this->geoState, *this->X, *this->A, U);
      else {//numFLuid>1
        dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                                &numSubCycles, *this->geoState, *this->A, U, nodeTag);
      }
    }
    else  //time step size with error estimation
      dt = this->timeState->computeTimeStep(it, dtLeft, &numSubCycles);
  }
  else    //if time step is repeated
    dt = this->timeState->computeTimeStepFailSafe(dtLeft, &numSubCycles);

  if(TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION && it == 1)
    this->timeState->setDtMin(dt * TsDesc<dim>::data->getCflMinOverCfl0());


  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  dtf = dt;
  dtfLeft = *dtLeft + dt;

//  fprintf(stderr,"dt = %e, dtfLeft = %e, dtLeft = %e.\n", dt, dtfLeft, *dtLeft);
  return dt;
}

//---------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U,increasingPressure); 
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
  this->output->cleanProbesFile();

  if (it == 0) {
    // First time step: compute GradP before computing forces
    this->spaceOp->computeGradP(*this->X, *this->A, U);
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeMaterialVolumesToDisk(it, 0.0, *this->A, &nodeTag);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTag);
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

  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeMaterialVolumesToDisk(it, t, *this->A, &nodeTag);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTag);
  this->output->writeProbesToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,nodeTag, this->distLSS, ghostPoints);
  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);

  TsRestart *restart2 = this->restart; // Bug: compiler does not accept this->restart->writeToDisk<dim,1>(...)
                                       //      it does not seem to understand the template
  restart2->template writeToDisk<dim,1>(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState);
  this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, this->distLSS->getStructPosition_n()); //KW: must be after writeToDisk

  this->output->updatePrtout(t);
  this->restart->updatePrtout(t);
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
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES) {
    this->ghostPoints->deletePointers();
    this->spaceOp->populateGhostPoints(this->ghostPoints,U,this->varFcn,this->distLSS,this->nodeTag);
  }

  this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface,  nodeTag, *this->R, this->riemann, riemannNormal, Nsbar, 0, ghostPoints);

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
  double t0 = this->timer->getTime();
  if(dynNodalTransfer)
    numStructNodes = dynNodalTransfer->numStNodes();

  if(!increasingPressure || recomputeIntersections) {
    for (int i=0; i<numStructNodes; i++) 
      Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
    this->spaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X,*this->A, Fs, numStructNodes, distLSS, *Wij, *Wji, 
                                    ghostPoints, this->postOp->getPostFcn(), &nodeTag);
  } else {
    this->com->fprintf(stderr, "OK. I am cheap...\n");
    if(unifPressure[0]==0) {
      this->com->fprintf(stderr,"ERROR: Detected pressure p = %e in Implosion Setup.\n", unifPressure[0]);
      exit(-1);
    }
    for (int i=0; i<numStructNodes; i++) {
      Fs[i][0] *= unifPressure[1]/unifPressure[0]; 
      Fs[i][1] *= unifPressure[1]/unifPressure[0]; 
      Fs[i][2] *= unifPressure[1]/unifPressure[0]; 
    }
  }

  this->timer->addEmbeddedForceTime(t0);
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
  if(dynNodalTransfer)
    numStructNodes = dynNodalTransfer->numStNodes();
  for (int i=0; i<numStructNodes; i++) {
    F[0]+=Fs[i][0]; F[1]+=Fs[i][1]; F[2]+=Fs[i][2];}

  M[0] = M[1] = M[2] = 0;
  Vec<Vec3D>& Xstruc = distLSS->getStructPosition();
  for (int i = 0; i < numStructNodes; ++i) {
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
    numStructNodes = dynNodalTransfer->numStNodes();
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(dt, dtLeft, v);
  }
}

//-------------------------------------------------------------------------------

template<int dim>
bool EmbeddedTsDesc<dim>::IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U)
{

  increasingPressure = false;
  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) {// max pressure was reached, so now we solve
//    this->com->fprintf(stdout, "max pressure reached\n"); 
    return true;
  } 

  increasingPressure = true;

  // max pressure not reached, so we do not solve and we increase pressure and let structure react
  
  if(this->mmh && !inSubCycling) {
    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);
    //recompute intersections
    double tw = this->timer->getTime();
    if((it-1)%intersector_freq==0) {
      this->com->fprintf(stderr,"recomputing fluid-structure intersections.\n");
      recomputeIntersections = true;
      this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag); 
    } else
      recomputeIntersections = false;

    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    nodeTag0 = this->nodeTag;
    nodeTag = this->distLSS->getStatus();

    //store previous states for phase-change update
    tw = this->timer->getTime();
    if(recomputeIntersections)
      this->spaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, U, this->Vtemp,
            *this->Weights, *this->VWeights, *this->Wstarij, *this->Wstarji,
            this->distLSS, (double*)this->vfar, (this->numFluid == 1 ? (DistVec<int>*)0 : &this->nodeTag));
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  } 
  
  // construct Wij, Wji from U. 
  DistSVec<double,dim> VV(this->getVecInfo());
  this->varFcn->conservativeToPrimitive(U,VV,&nodeTag);
  SubDomain **subD = this->domain->getSubDomain();

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
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

  double pnow = currentPressure(t);;
  this->com->fprintf(stdout, "about to increase pressure to %e\n", pnow*Pscale);
  this->domain->IncreasePressure(pnow, this->varFcn, U, nodeTag);
  unifPressure[0] = unifPressure[1];
  unifPressure[1] = pnow;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU) {

  if (this->fixSol == 1)
    this->domain->fixSolution(this->varFcn,U,dU,&this->nodeTag);
}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::currentPressure(double t)
{
  double p;
  if(implosionSetupType==EmbeddedTsDesc<dim>::LINEAR)
    p = Pinit + t*Prate;
  else if(implosionSetupType==EmbeddedTsDesc<dim>::SMOOTHSTEP) {
    double tbar = t/tmax;
    //fprintf(stderr,"t = %e, tmax = %e, tbar = %e.\n", t, tmax, tbar);
    p = Pinit + (Pfinal - Pinit)*tbar*tbar*tbar*(10.0-15.0*tbar+6.0*tbar*tbar);
  } else {
    this->com->fprintf(stderr,"ERROR! ImplosionSetup::Type = %d NOT recognized!\n", implosionSetupType);
    exit(-1);
  }
  return p;
}

