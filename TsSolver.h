#ifndef _TS_SOLVER_H_
#define _TS_SOLVER_H_

# include<IoData.h>

class IoData;

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
class TsSolver {

  ProblemDescriptor *probDesc;

  int resolve(typename ProblemDescriptor::SolVecType &, 
              IoData &);

public:

  TsSolver(ProblemDescriptor *);
  ~TsSolver() {}

  int solve(IoData &);

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

  dts = probDesc->computePositionVector(&lastIt, it, t);

  if (lastIt)
    probDesc->outputPositionVectorToDisk();

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

      // compute fluid subcyling time step
      //fprintf(stdout, "TsSolver.h 1\n");
      dt = probDesc->computeTimeStep(it, &dtLeft, U);
      t += dt;

      // estimate mesh position in subcycle
      //fprintf(stdout, "TsSolver.h 2\n");
      probDesc->interpolatePositionVector(dt, dtLeft);

      // compute control volumes and velocities
      //fprintf(stdout, "TsSolver.h 3\n");
      probDesc->computeMeshMetrics();

      //fprintf(stdout, "TsSolver.h 4\n");
      // Fluid Solution
      itNl += probDesc->solveNonLinearSystem(U);

      // compute the current aerodynamic force
//      fprintf(stdout, "TsSolver.h 5\n");
      probDesc->updateOutputToStructure(dt, dtLeft, U);

//      fprintf(stdout, "TsSolver.h 6\n");
      probDesc->updateStateVectors(U, it);

    } while (dtLeft != 0.0);

    lastIt = probDesc->checkForLastIteration(it, t, dt, U);

    probDesc->outputForces(ioData, &lastIt, it, itSc, itNl, t, dt, U);
    dts = probDesc->computePositionVector(&lastIt, it, t);

    probDesc->outputToDisk(ioData, &lastIt, it, itSc, itNl, t, dt, U);

  }

  return 0;

}

//------------------------------------------------------------------------------

#endif
