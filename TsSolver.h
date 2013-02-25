#ifndef _TS_SOLVER_H_
#define _TS_SOLVER_H_

# include<IoData.h>

class IoData;

//------------------------------------------------------------------------------
/** Class which handles the algorithmic organization of the solution for all problems */
template<class ProblemDescriptor>
class TsSolver {

  ProblemDescriptor *probDesc;

  int resolve(typename ProblemDescriptor::SolVecType &,
              IoData &);

public:

  TsSolver(ProblemDescriptor *);
  ~TsSolver() {}

  int solve(IoData &);

	int fsoSolve(IoData &);
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
  
  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);
  int status = resolve(U, ioData);
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

  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());

  //
  // Check that an input file for the solution is specified
  //
/*  if (ioData.input.solutions[0] == 0)
  {
    probDesc->fsoPrintTextOnScreen("\n !!! Shape Optimization Handler requires an input solution !!!\n\n");
    exit(1);
  }
*/

  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);

	probDesc->moveMesh(

  probDesc->fsoPrintTextOnScreen("******************************************\n");
  probDesc->fsoPrintTextOnScreen("*** Fluid Shape Optimization Interface ***\n");
  probDesc->fsoPrintTextOnScreen("******************************************\n");

  probDesc->fsoHandler(ioData, U);

  return 0;

}


//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int TsSolver<ProblemDescriptor>::resolve(typename ProblemDescriptor::SolVecType &U,
                                         IoData &ioData)
{
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

  // dts is structural time step
  double dt, dts;
  int it = probDesc->getInitialIteration();
  double t = probDesc->getInitialTime();

  // For an embedded viscous simulation with turbulence model, compute the distance to the wall
  probDesc->computeDistanceToWall(ioData);

  // setup solution output files
  probDesc->setupOutputToDisk(ioData, &lastIt, it, t, U);

  /** for embedded method: send force (if it>0) and receive disp (from Struct). */
  dts = probDesc->computePositionVector(&lastIt, it, t, U);

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
     
    bool solveOrNot = true;
    
    // For an embedded viscous simulation with turbulence model and moving object, compute the distance to the wall
    if (ioData.problem.alltype == ProblemData::_UNSTEADY_AEROELASTIC_ ||
        ioData.problem.alltype == ProblemData::_ACC_UNSTEADY_AEROELASTIC_ ||
        ioData.problem.alltype == ProblemData::_FORCED_) {
      probDesc->computeDistanceToWall(ioData);
    }
    bool repeat;
    do { // Subcycling
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
        stat = probDesc->solveNonLinearSystem(U, it);
        if (stat == -10){ // must redo iteration with a different CFL number, undo everything we have done so far
          probDesc->printf(1,"Found unphysical solution. Re-calculating CFL number and repeating iteration.\n");
          repeat = true;
          if (dU && dUPrev){
            *dU = *dUPrev;
            *dUPrev = *dUPrevPrev;
          }
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
    dts = probDesc->computePositionVector(&lastIt, it, t, U);
    probDesc->outputToDisk(ioData, &lastIt, it, itSc, itNl, t, dt, U);

  }
  return 0;

}

//------------------------------------------------------------------------------

#endif
