#include <ExplicitLevelSetTsDesc.h>

#include <math.h>
                                                                                                        
#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

#include <GeoSource.h>
#include <DistTimeState.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <NewtonSolver.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <PostOperator.h>
#include <DistVector.h>
//#include <Timer.h>


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
ExplicitLevelSetTsDesc<dim>::
ExplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  ExplicitTsDesc<dim>(ioData, geoSource, dom), Vg(dom->getNodeDistInfo())
{

  Communicator *com = dom->getCommunicator();
  DistSVec<double,3> X0(this->domain->getNodeDistInfo());
  this->domain->getReferenceMeshPosition(X0);
  LS = new LevelSet(X0, ioData, this->domain);
  Vgf = new DistSVec<double,dim>(this->domain->getNodeDistInfo());
  *Vgf =100.0;
  ImplicitData &implicitData = ioData.ts.implicit;
  this->ns = new NewtonSolver<ExplicitTsDesc<dim> >(this);
  mvpLS  = new MatVecProdLS<MatScalar,dim,1>(ioData, this->varFcn, this->timeState, this->geoState, this->spaceOp, this->domain);
  kspLS = createKrylovSolverLS(this->getVecInfo(), implicitData.newton.ksp.lsi, mvpLS, this->pc, this->com);
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitLevelSetTsDesc<dim>::~ExplicitLevelSetTsDesc()
{

  if (LS) delete LS;
  if (mvpLS) delete mvpLS;
  if (kspLS) delete kspLS;

}

//------------------------------------------------------------------------------
template<int dim>
template<int neq>
KspSolver<DistVec<double>, MatVecProd<dim,1>, KspPrec<neq,double>, Communicator, double> *
ExplicitLevelSetTsDesc<dim>::createKrylovSolverLS(const DistInfo &info, KspData &kspdata,
                                        MatVecProd<dim,1> *_mvpLS, KspPrec<neq,double> *_pc,
                                        Communicator *_com)
{
                                                                                                                                                         
  KspSolver<DistVec<double>, MatVecProd<dim,1>,
    KspPrec<neq>, Communicator> *_ksp = 0;
                                                                                                                                                         
  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);
  else if (kspdata.type == KspData::GCR)
    _ksp = new GcrSolver<DistVec<double>, MatVecProd<dim,1>,
                           KspPrec<neq>, Communicator>(info, kspdata, _mvpLS, _pc, _com);                                                                                                                                          
  return _ksp;
                                                                                                                                                         
}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
         
  this->geoState->setup2(this->timeState->getData());

  // initalize solution
  this->timeState->setup(this->input->solutions, this->bcData->getInletBoundaryVector(), this->bcData->getInterface(), *this->X, *U, ioData);

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U);
  else if (this->hth)
    this->hth->setup(&this->restart->frequency, &this->data->maxTime);

  this->timer->setSetupTime();
 
}

//------------------------------------------------------------------------------

template<int dim>
double ExplicitLevelSetTsDesc<dim>::computeTimeStep(int it, double *dtLeft, DistSVec<double,dim> &U)
{

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual);

  int numSubCycles = 1;
  double dt = this->timeState->computeTimeStep(this->data->cfl, dtLeft, &numSubCycles, *this->geoState, *this->A, U, LS->Phi);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*this->refVal->time, numSubCycles);

  return dt;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U)
{

  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U, LS->Phi, LS->Phi1, LS->Phi2);

}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitLevelSetTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, U, LS->Phi);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt,
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
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, LS->Phi);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

  LS->primitiveToConservative(U);
  LS->Phi1  = LS->Phi;
  LS->Phi2  = LS->Phi1; 
}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{
                                                                                                      
  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;
                                                                                                      
  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
                                                                                                      
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  LS->conservativeToPrimitive(U);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, LS->Phi);
  LS->primitiveToConservative(U);
  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState);
                                                                                                      
  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }
                                                                                                      
}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(this->postOp, *this->X, U);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitLevelSetTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft,
					  DistSVec<double,dim> &U)
{

  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(), this->timeState->getUn(), *this->X, U, work);
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U);

  if (this->hth)
    this->hth->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U);

}

//------------------------------------------------------------------------------
template<int dim>
int ExplicitLevelSetTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U)  {


  double t0 = this->timer->getTime(); 
    
// 1st step of Runge Kutta
  computeRKUpdate(U, this->k1, this->Ubc);
  this->spaceOp->getExtrapolationValue(U, this->Ubc, *this->X);
                                                                                                                                                         
  this->U0 = U - 0.5 * this->k1;
  this->spaceOp->applyExtrapolationToSolutionVector(this->U0, this->Ubc);
//2nd step
  computeRKUpdate(this->U0, this->k2, this->Ubc);
  this->spaceOp->getExtrapolationValue(this->U0, this->Ubc, *this->X);
                                                                                                                                                         
  this->U0 = U - 0.5 * this->k2;
  this->spaceOp->applyExtrapolationToSolutionVector(this->U0, this->Ubc);
                                                                                                                                                         
//3rd step
  computeRKUpdate(this->U0, this->k3, this->Ubc);
  this->spaceOp->getExtrapolationValue(this->U0, this->Ubc, *this->X);
                                                                                                                                                         
  this->U0 = U - this->k3;
  this->spaceOp->applyExtrapolationToSolutionVector(this->U0, this->Ubc);
                                                                                                                                                         
//4th step
  computeRKUpdate(this->U0, this->k4, this->Ubc);
  this->spaceOp->getExtrapolationValue(this->U0, this->Ubc, *this->X);

  U -= 1.0/6.0 * (this->k1 + 2.0 * (this->k2 + this->k3) + this->k4);
  this->spaceOp->applyExtrapolationToSolutionVector(U, this->Ubc);

  if (this->varFcn->getType() == VarFcn::GASINLIQUID) this->spaceOp->storeGhost(U, LS->Phi, *Vgf);

  this->spaceOp->applyBCsToSolutionVector(U);


  if (this->check == LiquidModelData::YES){
    int ierr = checkSolution(U);
    if (ierr > 0) exit(1);
  }else{
    int ierr2 = checkSolution(U);
    if (ierr2 > 0) this->com->fprintf(stderr, "pressure or density out of range\n");
  }

// Store the primitive variables
  this->varFcn->conservativeToPrimitive(U, Vg, &(LS->Phi));
// solve the Level Set equation
  int itsLS = this->ns->solveLS(LS->Phi, LS->Phi1, LS->Phi2, U);
//  for (int i=0; i<LS->Phi.size(); i++){
//      LS->Phi[i]  = LS->Phi[i]  +0.001;
//  }
  LS->Phi2  = LS->Phi1;
  LS->Phi1  = LS->Phi;

  if (this->varFcn->getType() == VarFcn::GASINLIQUID) {
    LS->conservativeToPrimitive(U);
// Restore the conservative state from the primitive state
    this->varFcn->primitiveToConservative(Vg, U, &(LS->Phi), &(LS->Phi2), Vgf);
    LS->primitiveToConservative(U);
    LS->Phi1  = LS->Phi;
  }
  else {
    this->varFcn->primitiveToConservative(Vg, U, &(LS->Phi), &(LS->Phi2), &Vg);
  }

  this->timer->addLevelSetSolutionTime(t0);

  return 1;

}

//------------------------------------------------------------------------------
                                                                                                                                                         
template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& U,
				DistSVec<double,dim>& dU, DistSVec<double,dim>& Ubc)
{
  this->spaceOp->applyBCsToSolutionVector(U);
  this->spaceOp->computeResidual(*this->X, *this->A, U, LS->Phi, dU);
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeFunctionLS(int it, DistVec<double> &Phi, 
						    DistVec<double> &Phi1, 
						    DistVec<double> &Phi2, 
                                                    DistSVec<double,dim> &U,
						    DistVec<double> &PhiF)
{
  this->spaceOp->computeResidualLS(*this->X, *this->A, Phi, U, PhiF);

  this->timeState->add_dAW_dtLS(it, *this->geoState, *this->A, Phi, Phi1, Phi2, PhiF);

}
//------------------------------------------------------------------------------
template<int dim>
void ExplicitLevelSetTsDesc<dim>::computeJacobianLS(int it, DistVec<double> &Phi,
                                                           DistVec<double> &Phi1,
                                                           DistVec<double> &Phi2,
                                                           DistSVec<double,dim> &U,
                                                           DistVec<double> &PhiF)
{

  mvpLS->evaluateLS(it, *this->X, *this->A, Phi, Phi1, Phi2, U, PhiF);

}
//------------------------------------------------------------------------------

template<int dim>
int ExplicitLevelSetTsDesc<dim>::solveLinearSystemLS(int it, DistVec<double> &b,
                                                  DistVec<double> &dQ)
{
 
  double t0 = this->timer->getTime();

  dQ = 0.0;
 
  kspLS->setup(it, 2, b);
 
  int lits = kspLS->solveLS(b, dQ);
  
  return lits;

}

