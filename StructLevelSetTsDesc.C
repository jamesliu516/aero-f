#include <StructLevelSetTsDesc.h>
#include <DistExactRiemannSolver.h>
#include <LevelSet/IntersectionFactory.h>

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
  TsDesc<dim>(ioData, geoSource, dom)
{

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain);
  distLSS = 0;

  const char *intersectorName = ioData.strucIntersect.intersectorName;
  if(intersectorName != 0)
    distLSS = IntersectionFactory::getIntersectionObject(intersectorName, *this->domain);
}

//------------------------------------------------------------------------------

template<int dim>
StructLevelSetTsDesc<dim>::~StructLevelSetTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
//from TsDesc::setupTimeStepping. Skipped alot mesh motion handlers...
{
  this->geoState->setup2(this->timeState->getData());
  //TODO: timeState->setup different as in LevelSetTsDesc.
  this->timeState->setup(this->input->solutions, this->bcData->getInletBoundaryVector(), *this->X, *U);
  this->distLSS->initialize(this->domain,*this->X);
}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U)
// same with both TsDEsc::computeTimeStep and LevelSetTsDesc::computeTimeStep, except "printf" line.
{
//  this->com->barrier(); this->com->fprintf(stderr,"in SLSTsDesc::computeTimeStep.\n");
  double t0 = this->timer->getTime();
//  this->com->barrier(); this->com->fprintf(stderr,"in SLSTsDesc::computeTimeStep (1).\n");
  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);
//  this->com->barrier(); this->com->fprintf(stderr,"in SLSTsDesc::computeTimeStep (2).\n");
  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft,
                            &numSubCycles, *this->geoState, *this->X, *this->A, U);
  //TODO: create new computeTimeStep in TimeState. Set dt in ghost region to be a large constant.
//  this->com->barrier(); this->com->fprintf(stderr,"dt = %f.\n", dt);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;
}

//---------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{ // same as in TsDesc::updateStateVectors.
//  this->com->barrier(); this->com->fprintf(stderr,"in SLSTsDesc::updateStateVectors.\n");
  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U);
}

//------------------------------------------------------------------------------

template<int dim>
int StructLevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
// same as in LevelSetTsDesc::checkSolution. TODO: no phi. how to check?
{
  int ierr = 0;
  ierr = this->domain->checkSolution(this->varFcn, U); //also check ghost nodes.
  return ierr;
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

  if (it == 0) {
    // First time step: compute GradP before computing forces
    this->spaceOp->computeGradP(*this->X, *this->A, U);

    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, distLSS->getPhi());
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

//  this->com->barrier(); this->com->fprintf(stderr,"DONE.\n");
}

//------------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                                             double t, double dt, DistSVec<double,dim> &U)
{ // comes from TsDesc::outputToDisk
  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;

  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,
                                   distLSS->getPhi());
  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
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
{ //TODO:simply copied from TsDesc::outputForces. need to be modified to work.
  // add a new writeForcesToDisk in TsOutput. then go into PostOperator and add a new computeForceAndMoment.
  double cpu = this->timer->getRunTime();
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

//-----------------------------------------------------------------------------

template<int dim>
void StructLevelSetTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft,
                                          DistSVec<double,dim> &U)
{ //TODO: need to be re-written.
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
}

//------------------------------------------------------------------------------

template<int dim>
double StructLevelSetTsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{ //from TsDesc::computeResidualNorm. only compute residual for Phi>0 TODO: shouldn't do this for shell.
  this->spaceOp->computeResidual(*this->X, *this->A, U, distLSS, *this->R, this->riemann, 0);
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





























