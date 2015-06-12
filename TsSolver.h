#ifndef _TS_SOLVER_H_
#define _TS_SOLVER_H_

# include<IoData.h>
#include <ErrorHandler.h>

class IoData;

//------------------------------------------------------------------------------
/** Class which handles the algorithmic organization of the solution for all problems */
template<class ProblemDescriptor>
class TsSolver {

  ProblemDescriptor *probDesc;

  int resolve(typename ProblemDescriptor::SolVecType &,
              IoData &);

  int solveWithMultipleICs(typename ProblemDescriptor::SolVecType &,
              IoData &);
public:

  TsSolver(ProblemDescriptor *);
  ~TsSolver() {}

  int solve(IoData &);

  int fsoSolve(IoData &);
  int fsisoSolve(IoData &);
// Included (MB)
  int fsaSolve(IoData &);


};

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
TsSolver<ProblemDescriptor>::TsSolver(ProblemDescriptor *prbd)
{

  probDesc = prbd;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::solve(IoData &ioData)
{
  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());

  int status;
  if (ioData.problem.solveWithMultipleICs) { 
    // solve the system starting from multiple initial conditions
    status = solveWithMultipleICs(U,ioData);
  } else { 
    // standard solve
    probDesc->setupTimeStepping(&U, ioData);  // initialize solutions and geometry
    status = resolve(U, ioData);
  }
  return status;
  
}

//------------------------------------------------------------------------------
template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::solveWithMultipleICs(typename ProblemDescriptor::SolVecType &U, IoData &iod)
{
  // Solve the system starting from multiple initial conditions.
  // This is intended to help with training ROM simulations

  FILE *inFP = fopen(iod.input.multiSolutions,"r");
  if (!inFP)  {
    //probDesc->com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", iod.input.multiSolutions);
    exit (-1);
  }
  int nVecs, _n;
  _n = fscanf(inFP, "%d",&nVecs);

  char solnFile[500];
  double tmp;
  int status;

  int initialIt = probDesc->getInitialIteration();
  double initialTime = probDesc->getInitialTime();

  for (int iVec=0; iVec < nVecs; ++iVec) {
    _n = fscanf(inFP, "%s", solnFile);
    probDesc->readICFromDisk(solnFile, iVec, nVecs, U);
    probDesc->setRestartIterationAndTime(initialIt, initialTime);
    probDesc->setupTimeStepping(&U, iod);
    status = resolve(U, iod);
  }
  fclose(inFP);

  return status;  
}

//------------------------------------------------------------------------------

// Included (MB)
template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::fsaSolve(IoData &ioData)
{

  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());


  //
  // Check that an input file for the solution is specified
  //
  if (ioData.input.solutions[0] == 0)
  {
    probDesc->fsaPrintTextOnScreen("\n !!! SensitivityAnalysis requires an input solution !!!\n\n");
    exit(1);
  }


  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);

  probDesc->fsaPrintTextOnScreen("**********************************\n");
  probDesc->fsaPrintTextOnScreen("*** Fluid Sensitivity Analysis ***\n");
  probDesc->fsaPrintTextOnScreen("**********************************\n");
  
  probDesc->fsaHandler(ioData, U);

  return 0;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::fsoSolve(IoData &ioData)
{

  probDesc->printf(0,"******************************************\n");
  probDesc->printf(0,"*** Fluid Shape Optimization Interface ***\n");
  probDesc->printf(0,"******************************************\n");

  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());

  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);
  probDesc->fsoInitialize(ioData, U);
  resolve(U, ioData);

  probDesc->fsoHandler(ioData, U);

  return 0;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::fsisoSolve(IoData &ioData)
{

  probDesc->printf(0,"*******************************************************\n");
  probDesc->printf(0,"*** Steady Aeroelastic Shape Optimization Interface ***\n");
  probDesc->printf(0,"*******************************************************\n");

  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());
  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);
  probDesc->fsoInitialize(ioData, U);
  resolve(U, ioData);

  ioData.sa.fsiFlag = true;
  probDesc->fsoHandler(ioData, U);
  probDesc->printf(0," ***** fsisoSolve is done ********\n");

  return 0;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::resolve(typename ProblemDescriptor::SolVecType &U,
                                         IoData &ioData)
{
  int verboseFlag = ioData.problem.verbose;
  bool lastIt = false;

  typename ProblemDescriptor::SolVecType *dU = NULL;
  typename ProblemDescriptor::SolVecType *dUPrev = NULL;
  typename ProblemDescriptor::SolVecType *dUPrevPrev = NULL;
  double angle = -2.0;

  if (ioData.ts.cfl.strategy == CFLData::DIRECTION || ioData.ts.cfl.strategy == CFLData::HYBRID){
    dU = new typename ProblemDescriptor::SolVecType(probDesc->getVecInfo());
    dUPrev = new typename ProblemDescriptor::SolVecType(probDesc->getVecInfo());
    dUPrevPrev = new typename ProblemDescriptor::SolVecType(probDesc->getVecInfo());
    (*dU) = 0.0;
    (*dUPrev) = 0.0;
    (*dUPrevPrev) = 0.0;
  }

  typename ProblemDescriptor::SolVecType *UPrev = new typename ProblemDescriptor::SolVecType(probDesc->getVecInfo());
  (*UPrev) = 0.0;

  // dts is structural time step
  double dt, dts;
  int it = probDesc->getInitialIteration();
  double t = probDesc->getInitialTime();

  // setup solution output files
  probDesc->setupOutputToDisk(ioData, &lastIt, it, t, U);

  /** for embedded method: send force (if it>0) and receive disp (from Struct). */
  dts = probDesc->computePositionVector(&lastIt, it, t, U); // [F] receive displacement from structure ...

  // For an embedded viscous simulation with turbulence model, compute the distance to the wall
  probDesc->computeDistanceToWall(ioData);

  if (lastIt)
    probDesc->outputPositionVectorToDisk(U);

  while (!lastIt) {
    probDesc->resetOutputToStructure(U);
    int stat = 0;
    int itSc = 0;
    int itNl = 0;
    int itNlLS = 0;

    // initialize remaining time in fluid subcycling
    double dtLeft = dts;
    it++;
    
    *(probDesc->getTimeIt()) = it;
 
    bool solveOrNot = true;
    
    bool repeat;
    do { // Subcycling
      (*UPrev) = U;

      repeat = false;
      double dtLeftPrev = dtLeft;
      stat = 0;
      itSc++;
      probDesc->setCurrentTime(t,U);
 
      if(probDesc->structureSubcycling() || //in this case computeTimeStep is called in computePositionVector
         (it>1 && probDesc->willNotSolve(dtLeft,t)) ) {//in this case AERO-F should never subcycle
        probDesc->setFluidSubcycling(false);
        dt = dtLeft;
        dtLeft = 0.0;
      }
      else{
        if (dU && dUPrev && dUPrev->norm() != 0) angle = ((*dU) * (*dUPrev))/(dU->norm()*dUPrev->norm());
        else angle = -2.0;
        dt = probDesc->computeTimeStep(it, &dtLeft, U, angle);
      }
      
      t += dt;

      probDesc->setCurrentTimeStep(dt);

      // update coefficients for enforcing the Farfield BC.
      probDesc->updateFarfieldCoeffs(dt);
      // estimate mesh position in subcycle
      probDesc->interpolatePositionVector(dt, dtLeft);
      // compute control volumes and velocities
      probDesc->computeMeshMetrics();
      // Fluid Solution
      solveOrNot = probDesc->IncreasePressure(it,dt,t,U);
      if (ioData.problem.solvefluid == ProblemData::OFF) {
        solveOrNot = false;
      }
      if(solveOrNot){
        if (dU && dUPrev){
          *dUPrevPrev = *dUPrev;
          *dUPrev = *dU;
          *dU = -1.0*U;
        }
        if(probDesc->getErrorHandler()) probDesc->getErrorHandler()->clearError(ErrorHandler::ALL);
        probDesc->checkLocalRomStatus(U, it);
        stat = probDesc->solveNonLinearSystem(U, it);
        if(probDesc->getErrorHandler()) probDesc->getErrorHandler()->reduceError();

        if(probDesc->getTsParams()) probDesc->getTsParams()->resolveErrors();

        //probDesc->getErrorHandler()->printError();

        // stat = -10 signals that the time iteration must be redone!
        // Regardless of what happens elsewhere, this should work.
        if (probDesc->getErrorHandler()->globalErrors[ErrorHandler::REDO_TIMESTEP]){ // || stat == -10 // must redo iteration with a different CFL number, undo everything we have done so far 
          probDesc->getErrorHandler()->globalErrors[ErrorHandler::REDO_TIMESTEP]=0;
          probDesc->printf(1,"Repeating time-step.\n");
	  probDesc->printf(1, "WARNING: Alex modified this... hopefully to fix a bug.\n");
          //probDesc->setFailSafe(true);
          U = (*UPrev); // Reset U to its previous state
          repeat = true;
          // Reset directions for direction strategy
          if (dU && dUPrev){
            *dU = *dUPrev;
            *dUPrev = *dUPrevPrev;
          }
          // undo time step and subcycling
          t -= dt;
          itSc--;
          dtLeft = dtLeftPrev;
          continue;
        }

        if (dU && dUPrev) *dU += U;
        if(stat>0){
          itNl += stat;
          // compute the current aerodynamic force
          probDesc->updateOutputToStructure(dt, dtLeft, U);
          probDesc->updateStateVectors(U, it);
        }
        else{
          if(itSc > 200){ 
            probDesc->printf(1, "Fail safe failed! \n",itSc);
            exit(-1);
          }
          probDesc->printf(1, "stat: %i \n",stat);
          probDesc->printf(1, "itSc:  %i \n",itSc);
          t -= dt;
          probDesc->setFailSafe(true);
        }
      } else {
        probDesc->updateOutputToStructure(dt, dtLeft, U);
        probDesc->updateStateVectors(U, it);
      }

    } while (repeat || dtLeft != 0.0 || stat<0);

// Modified (MB)
    lastIt = probDesc->checkForLastIteration(ioData, it, t, dt, U);

    probDesc->outputForces(ioData, &lastIt, it, itSc, itNl, t, dt, U);
    dts = probDesc->computePositionVector(&lastIt, it, t, U); // [F] send force to structure

  // For an embedded viscous simulation with turbulence model and moving object, compute the distance to the wall
    if ( (ioData.problem.framework == ProblemData::EMBEDDED) || 
         (ioData.problem.framework == ProblemData::EMBEDDEDALE) )
      if (ioData.problem.alltype == ProblemData::_UNSTEADY_AEROELASTIC_ ||
          ioData.problem.alltype == ProblemData::_ACC_UNSTEADY_AEROELASTIC_ ||
          ioData.problem.alltype == ProblemData::_FORCED_) {
        if (!lastIt) probDesc->computeDistanceToWall(ioData);
      }

    probDesc->outputToDisk(ioData, &lastIt, it, itSc, itNl, t, dt, U);

  }
  return 0;

}

//------------------------------------------------------------------------------

#endif
