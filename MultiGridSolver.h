#pragma once

#include<IoData.h>

class IoData;

//------------------------------------------------------------------------------
/** Class which handles the algorithmic organization of the solution for all problems */
template<class ProblemDescriptor>
class MultiGridSolver {

  ProblemDescriptor *probDesc;

  int resolve(typename ProblemDescriptor::SolVecType &,
              IoData &);

  int num_presmooth,num_postsmooth;
  int num_finesweeps;

 
  typename ProblemDescriptor::MultiGridKernelType* multiGridKernel; 

  typedef typename ProblemDescriptor::MultiGridKernelType::MultiGridSmoother
    MultiGridSmoother;

  typename ProblemDescriptor::MultiGridKernelType* pOwner; 
 
  class MySmoother : public MultiGridSmoother {

      ProblemDescriptor* probDesc;
      FluxFcn** fluxFcn;
      VarFcn* varFcn;

      double R0;
      bool hasR0;

    public:

     MySmoother(typename ProblemDescriptor::MultiGridKernelType* owner,
                ProblemDescriptor* probDesc) :
       MultiGridSmoother(owner), probDesc(probDesc) {

       fluxFcn = probDesc->getSpaceOperator()->getFluxFcn();
       varFcn = probDesc->getSpaceOperator()->getVarFcn();
       hasR0 = false;
     }

     ~MySmoother() {

     }

     void smooth(int level, typename ProblemDescriptor::SolVecType& x,
                 const typename ProblemDescriptor::SolVecType& f,int steps) {

       typename ProblemDescriptor::SolVecType& R = this->pOwner->getResidual(level);
       typename ProblemDescriptor::SolVecType& dx = this->pOwner->getDX(level);
       double dummy = 0.0;
       for (int i = 0; i < steps; ++i) {

         if (level == 0) {

           probDesc->computeTimeStep(0,&dummy, x);
           std::cout << "CFL = " << probDesc->data->cfl << std::endl;
           probDesc->computeFunction(0, x, R);
           probDesc->computeJacobian(0, x, R);
           probDesc->setOperators(x);
         } else {
 
           this->pOwner->getLevel(level)->RestrictFaceVector(*this->pOwner->getLevel(level-1),
                                                   (level == 1 ? probDesc->getSpaceOperator()->getDistBcData()->getFaceStateVector() : 
                                                     this->pOwner->getOperator(level-1)->getBcData().getFaceStateVector()),
                                                   this->pOwner->getOperator(level)->getBcData().getFaceStateVector());
           this->varFcn->conservativeToPrimitive(x, this->pOwner->getTemporaryV(level));
           this->pOwner->getOperator(level)->updateStateVectors(x);
           this->pOwner->getOperator(level)->computeTimeStep(probDesc->data->cfl*0.75,
                                                             varFcn,this->pOwner->getTemporaryV(level));
           this->pOwner->getOperator(level)->computeResidual(this->pOwner->getTemporaryV(level),
                                                     x, fluxFcn,
                                                     probDesc->getConstantRecFcn(),
                                                     R);
           this->pOwner->getOperator(level)->computeJacobian(x,
                                                       this->pOwner->getTemporaryV(level),
                                                       fluxFcn, this->pOwner->getA(level));
           this->pOwner->setupPreconditioner(level);
         }
         R = f-R;
         if (!hasR0) {

           R0 = sqrt(R*R);
           hasR0 = true;
         }
         std::cout << "Residual = " << sqrt(R*R)/R0 << std::endl;

         dx = 0.0;
         if (level == 0) {
 
           probDesc->solveLinearSystem(0, R,dx);
         } else {
           char fn[32];
           sprintf(fn,"R10%d",i);
           this->pOwner->getLevel(level)->writeXpostFile(fn,R, 0);
           this->pOwner->kspSolve(level, R, dx);
           this->pOwner->getLevel(level)->writeXpostFile("U10p",x, 0);
         }

         //if (level == 0)
           x += dx;
         if (level == 0) {
           probDesc->updateStateVectors(x, 0);
           probDesc->monitorConvergence(0, x);
           R = f-probDesc->getCurrentResidual();
         } else {
           char fn[32];
           sprintf(fn,"Udx%d",i);
           this->pOwner->getLevel(level)->writeXpostFile(fn,dx, 0);
         }
       }
         if (level == 0) {
           probDesc->monitorConvergence(0, x);
           R = f-probDesc->getCurrentResidual();
         }
     }

     void applyOperator(int level, typename ProblemDescriptor::SolVecType& f,
                        typename ProblemDescriptor::SolVecType& x) {

       if (level == 0) {

         probDesc->computeFunction(0, f, x);
       } else {

         this->varFcn->conservativeToPrimitive(f, this->pOwner->getTemporaryV(level));
         this->pOwner->getOperator(level)->computeResidual(this->pOwner->getTemporaryV(level),
                                                     f, fluxFcn,
                                                     probDesc->getConstantRecFcn(),
                                                     x, false);
       }
     }
  };
  
  MySmoother* mySmoother;

public:

  MultiGridSolver(ProblemDescriptor *);
  ~MultiGridSolver() {}

  int solve(IoData &);

// Included (MB)
  int fsaSolve(IoData &);

};

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
MultiGridSolver<ProblemDescriptor>::MultiGridSolver(ProblemDescriptor *prbd)
{

  probDesc = prbd;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int MultiGridSolver<ProblemDescriptor>::solve(IoData &ioData)
{
  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());

  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);
  mySmoother = new MySmoother(probDesc->getMultiGridKernel(),probDesc);
  
  int status = resolve(U, ioData);
  return status;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int MultiGridSolver<ProblemDescriptor>::resolve(typename ProblemDescriptor::SolVecType &U,
                                                IoData &ioData)
{
  bool lastIt = false;

  // dts is structural time step
  double dt, dts;
  int it = probDesc->getInitialIteration();
  double t = probDesc->getInitialTime();
  // setup solution output files
  probDesc->setupOutputToDisk(ioData, &lastIt, it, t, U);

  /** for embedded method: send force (if it>0) and receive disp (from Struct). */
  dts = probDesc->computePositionVector(&lastIt, it, t, U);

  if (lastIt)
    probDesc->outputPositionVectorToDisk(U);
    
  // compute control volumes and velocities
  probDesc->computeMeshMetrics();
    
  typename ProblemDescriptor::SolVecType zero(probDesc->getVecInfo());
  zero = 0.0;

  multiGridKernel = probDesc->getMultiGridKernel();
  multiGridKernel->setSmoother(mySmoother);
  multiGridKernel->setParameters(2,0,2,1.0,1);
  multiGridKernel->initialize();
  multiGridKernel->setGeometric();

  while (!lastIt) {
    probDesc->resetOutputToStructure(U);
    it++;
    
    bool solveOrNot = true;
    probDesc->setCurrentTime(t,U);

    multiGridKernel->cycle(zero, U);

    // compute the current aerodynamic force
    probDesc->updateOutputToStructure(0.0, 0.0, U);

// Modified (MB)
    lastIt = probDesc->checkForLastIteration(ioData, it, t, 0.00, U);

    probDesc->outputForces(ioData, &lastIt, it, 1, 1, t, 0.0, U);
    dts = probDesc->computePositionVector(&lastIt, it, t, U);
    probDesc->outputToDisk(ioData, &lastIt, it, 1, 1, t, 0.0, U);
  }
  return 0;

}

//------------------------------------------------------------------------------
