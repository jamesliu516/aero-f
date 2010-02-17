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
int TsSolver<ProblemDescriptor>::resolve(typename ProblemDescriptor::SolVecType &U,
                                         IoData &ioData)
{

  bool lastIt = false;

  // dts is structural time step
  double dt, dts;

  int it = probDesc->getInitialIteration();
  double t = probDesc->getInitialTime();
  // setup solution output files

  probDesc->setupOutputToDisk(ioData, &lastIt, it, t, U);

  /** for embedded method: send force (if it>0) and send receive disp (from Struct). */
  dts = probDesc->computePositionVector(&lastIt, it, t);

  if (lastIt)
    probDesc->outputPositionVectorToDisk(U);

  while (!lastIt) {
    probDesc->resetOutputToStructure(U);
    int itSc = 0;
    int itNl = 0;
    int itNlLS = 0;

    // initialize remaining time in fluid subcycling
    double dtLeft = dts;
    it++;

    do {
      itSc++;
      dt = probDesc->computeTimeStep(it, &dtLeft, U);
      t += dt;

      // estimate mesh position in subcycle
      probDesc->interpolatePositionVector(dt, dtLeft);
      // compute control volumes and velocities
      probDesc->computeMeshMetrics();
      // Fluid Solution
      bool solveOrNot = probDesc->IncreasePressure(dt,t,U);
      if(solveOrNot){
        //fprintf(stdout, "No increase in pressure\n");
        itNl += probDesc->solveNonLinearSystem(U);
      }
      // compute the current aerodynamic force
      probDesc->updateOutputToStructure(dt, dtLeft, U);
      probDesc->updateStateVectors(U, it);
    } while (dtLeft != 0.0);


// Modified (MB)
    lastIt = probDesc->checkForLastIteration(ioData, it, t, dt, U);

    probDesc->outputForces(ioData, &lastIt, it, itSc, itNl, t, dt, U);
    dts = probDesc->computePositionVector(&lastIt, it, t);
    probDesc->outputToDisk(ioData, &lastIt, it, itSc, itNl, t, dt, U);
  }
  return 0;

}

//------------------------------------------------------------------------------

#endif
