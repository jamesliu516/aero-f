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
  double residual0 = probDesc->computeResidualNorm(U);
  fprintf(stderr, "right after 'setupTimeStepping', data->residual = %lf.\n", residual0);
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

  FILE *forceFile = fopen("/home/icmewang/Simulations/Prolate3D_force/ghost/results/ghostforce","w");
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
//      if (ioData.input.ghostsolid.runit == 1)  
//        probDesc->updateGhostFluid(U);
      // compute fluid subcyling time step
      dt = probDesc->computeTimeStep(it, &dtLeft, U);
      t += dt;

      // estimate mesh position in subcycle
      probDesc->interpolatePositionVector(dt, dtLeft);
      // compute control volumes and velocities
      probDesc->computeMeshMetrics();
      // Fluid Solution
      itNl += probDesc->solveNonLinearSystem(U);
      // compute the current aerodynamic force
      probDesc->updateOutputToStructure(dt, dtLeft, U);
      probDesc->updateStateVectors(U, it);
      if (ioData.input.ghostsolid.runit == 1) {
        Vec3D totalForce;
        probDesc->updateGhostFluid(U,totalForce,dt);
        fprintf(stderr,"%f, %f, %f.\n", totalForce[0],totalForce[1],totalForce[2]);
        fprintf(forceFile, "%d  %f  %f  %f  %f\n", it, t, totalForce[0], totalForce[1], totalForce[2]);
      }
    } while (dtLeft != 0.0);

// Modified (MB)
    lastIt = probDesc->checkForLastIteration(ioData, it, t, dt, U);

    probDesc->outputForces(ioData, &lastIt, it, itSc, itNl, t, dt, U);
    dts = probDesc->computePositionVector(&lastIt, it, t);

//    if (ioData.input.ghostsolid.runit == 1)
//      probDesc->updateGhostFluid(U);
 
    probDesc->outputToDisk(ioData, &lastIt, it, itSc, itNl, t, dt, U);

  }
  fclose(forceFile);
  return 0;

}

//------------------------------------------------------------------------------

#endif
