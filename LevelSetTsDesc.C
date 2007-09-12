#include <LevelSetTsDesc.h>


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
LevelSetTsDesc<dim>::
LevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), Vg(this->getVecInfo()),
	PhiV(this->getVecInfo())
{

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  LS = new LevelSet(ioData, this->domain);

	Vgf = 0;
  Vgfweight = 0;
	if(ioData.mf.typePhaseChange == MultiFluidData::EXTRAPOLATION){
    Vgf = new DistSVec<double,dim>(this->getVecInfo());
    Vgfweight = new DistVec<double>(this->getVecInfo());
    *Vgf =-1.0;
    *Vgfweight =0.0;
  }

	frequencyLS = ioData.mf.frequency;

}

//------------------------------------------------------------------------------

template<int dim>
LevelSetTsDesc<dim>::~LevelSetTsDesc()
{

  if (LS) delete LS;
  if (Vgf) delete Vgf;

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{

  this->geoState->setup2(this->timeState->getData());

  // initalize solution
	if(ioData.mf.problem == MultiFluidData::SHOCKTUBE)
    this->timeState->setup(this->input->solutions, this->bcData->getOutletBoundaryVector(), this->bcData->getInterface(), *this->X, *U, ioData);
	else if(ioData.mf.problem == MultiFluidData::BUBBLE)
    this->timeState->setup(this->input->solutions, this->bcData->getInletBoundaryVector(), this->bcData->getInterface(), *this->X, *U, ioData);

  LS->setup(this->input->levelsets, *this->X, Phi, *U, ioData);
	this->com->printf(2, "*** Warning: setup of shocktube problem uses xb and yb for initialization ***\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U);
  else if (this->hth)
    this->hth->setup(&this->restart->frequency, &this->data->maxTime);

	*this->Xs = *this->X;

  this->timer->setSetupTime();
}

//------------------------------------------------------------------------------

template<int dim>
double LevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
				    DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);

  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft, &numSubCycles, *this->geoState, *this->A, U, Phi);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*this->refVal->time, numSubCycles);

  this->timer->addFluidSolutionTime(t0);

  return dt;

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  this->geoState->update(*this->X, *this->A);
  LS->update(Phi);
  this->timeState->update(U, LS->Phin, LS->Phinm1, LS->Phinm2);
  if(frequencyLS > 0 && it%frequencyLS == 0){
    LS->conservativeToPrimitive(Phi,PhiV,U);
    LS->reinitializeLevelSet(*this->geoState,*this->X, *this->A, U, PhiV);
    LS->primitiveToConservative(PhiV,Phi,U);
  }

}

//------------------------------------------------------------------------------

template<int dim>
int LevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, U, Phi);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                                                    int it, double t, DistSVec<double,dim> &U)
{
 
  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U);

  this->output->setMeshMotionHandler(dynamic_cast<RigidMeshMotionHandler *>(this->mmh));
  this->output->openAsciiFiles();

  if (it == 0) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, Phi);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void LevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
                                                                                                      
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, Phi);
  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS);

  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }

}
