#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
using std::sort;


#include <Domain.h>
#include <Modal.h>
#include <IoData.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <MatVecProd.h>
#include <Timer.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <PostOperator.h>

#ifdef DO_MODAL
#include <arpack++/include/ardsmat.h>
//#include <arpack++/include/ardnsmat.h>
#include <arpack++/include/ardssym.h>
#endif

#include <LinkF77.h>
#include <string.h>
extern "C"      {

   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);

   void F77NAME(thinsvd)(int &, int &, int &, int &, int&, int&, int&, int&, int&, double *, int &,
                        int &, int &, int &, int &, int &, int &, double *U, double *S, double *V,
                        int &, double *, int &);
  void F77NAME(lworksize)(int &, int &, int &, int &, int&, int&, int&, int&, int&, int &);

}

//----------------------------------------------------------------------------------

template <int dim>
ModalSolver<dim>::ModalSolver(Communicator *_com, IoData &_ioData, Domain &dom) : 
           domain(dom), mX(0, dom.getNodeDistInfo() ), Xref(dom.getNodeDistInfo()), 
           Uref(dom.getNodeDistInfo()), DX(0, dom.getNodeDistInfo()), 
           DE(0, dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()), controlVolComp(dom.getNodeDistInfo())  {

  com = _com;
  double f = 0;
  double pi = 3.14159265358979;
  ioData = &_ioData; 
  const char *modeFile = ioData->linearizedData.strModesFile;

  DistSVec<double, dim> tmpVec(domain.getNodeDistInfo());
  DistSVec<double, 3> Xtmp(domain.getNodeDistInfo());

  if (strcmp(modeFile, "") != 0)  {
    com->fprintf(stderr, " ... Reading Modefile %s\n", modeFile);
    modeFile = ioData->linearizedData.strModesFile;
    domain.readVectorFromFile(modeFile, 0, &f, Xtmp);
    nStrMode = int(f);
    if (ioData->linearizedData.numStrModes > nStrMode)  {
      com->fprintf(stderr, " *** WARNING: Setting number of structural modes to number in file: %d\n",
                 nStrMode);
    }
    else
      nStrMode = ioData->linearizedData.numStrModes;
      
  }
  else  {
    com->fprintf(stderr, " ... Running Fluid Alone without structural mode file\n");
    nStrMode = 0;
  }

  K = new double[nStrMode];

  // We read the modal deformations
  mX.resize(nStrMode);
  com->fprintf(stderr, " ... Number of modes in file: %d\n", nStrMode);
  
  for(int iMode = 0; iMode < nStrMode; ++iMode) {
    domain.readVectorFromFile(modeFile, iMode+1, &f, mX[iMode]);
    com->fprintf(stderr, "Mode %d f = %e disp^2 = %e\n", iMode+1, f, mX[iMode]*mX[iMode]);
    K[iMode] = f*f*pi*pi*4.0;
  }

  tInput = new TsInput(*ioData);

  tOutput = 0;
  tRestart = 0;
  bcData = 0;
  geoState = 0;
  tState = 0;
  spaceOp = 0;
  postOp = 0;
  HOp = 0;
  HOp2 = 0;
  HOp2step1 = 0;
  HOpstep2 = 0;
  HOpstep3 = 0;
  pc = 0;
  ksp = 0;
  ksp2 = 0;
  ksp3 = 0;
  kspComp = 0;
  kspCompGcr = 0;

}

//-------------------------------------------------------------------------------
template <int dim> 
void ModalSolver<dim>::solve()  {

  // set up Timers
  Timer *modalTimer = domain.getTimer();
  double t0;
  
  if (ioData->problem.alltype == ProblemData::_INTERPOLATION_)
    interpolatePOD();
  else  {
    preProcess();
    modalTimer->setSetupTime();
    if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_) {
      t0 = modalTimer->getTime();
      constructPOD();
      modalTimer->addPodConstrTime(t0);
    }
    else {
      t0 = modalTimer->getTime();
      solveInTimeDomain();
      //modalTimer->addRomSolTime(t0);
      modalTimer->addFluidSolutionTime(t0);
    }
  }

  modalTimer->setRunTime();
  modalTimer->print(modalTimer);
}

//-------------------------------------------------------------------------------
template <int dim> 
void ModalSolver<dim>::solveInTimeDomain()  {

  Timer *modalTimer = domain.getTimer();
  
  double sdt = ioData->linearizedData.stepsize;

  // Initial Conditions
  double *delY = new double [nStrMode];;
  double *delU = new double [nStrMode];;

  //Initial Conditions
  int i;
  for (i=0; i<nStrMode; ++i) {
    delY[i] = 0.0;
    delU[i] = 0.0;
  }

  int nSteps = ioData->ts.maxIts;

  if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_ || 
      ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_ ||
      ioData->linearizedData.type == LinearizedData::FORCED)  {
    VecSet<DistSVec<double, dim> > snaps(0, domain.getNodeDistInfo());

    int modeNum = ioData->linearizedData.modeNumber - 1;

    if (ioData->linearizedData.type == LinearizedData::FORCED)
      com->fprintf(stderr, " ... Running Forced Oscillations w/no initializaitons\n");
    else if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_)
      com->fprintf(stderr, " ... Running Unsteady Linearized w/no structural initializaitons\n");
    else if (ioData->linearizedData.initCond == LinearizedData::DISPLACEMENT)  {
      com->fprintf(stderr, " ... Initializing with Displacement Mode: %d\n", ioData->linearizedData.modeNumber);
      delU[modeNum] = ioData->linearizedData.amplification;
    }
    else if (ioData->linearizedData.initCond == LinearizedData::VELOCITY)  {
      com->fprintf(stderr, " ... Initializing with Velocity Mode: %d\n", ioData->linearizedData.modeNumber);
      delY[modeNum] = ioData->linearizedData.amplification;
    }
    else  {
      com->fprintf(stderr, "*** Error: Bad Initial Condition: %d\n", ioData->linearizedData.initCond);
      exit (-1);
    }

    int dummySnap = 0;
    timeIntegrate(snaps, nSteps, 0, delU, delY, dummySnap, sdt);

  }

  if (ioData->problem.alltype == ProblemData::_ROM_ || 
      ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {

    int nPodVecs = ioData->linearizedData.numPOD;
    VecSet<DistSVec<double, dim> > podVecs(nPodVecs, domain.getNodeDistInfo());
    readPodVecs(podVecs, nPodVecs);

    // print out ROM if requested
    if (ioData->output.transient.romFile[0] != 0)  {
      
      if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {
        VecSet<Vec<double> > outputRom(nPodVecs, nStrMode);
        formOutputRom(outputRom, podVecs, nPodVecs);

        evalAeroSys(outputRom, podVecs, nPodVecs);
      }
      else
        evalFluidSys(podVecs, nPodVecs);
    }


    // timeIntegrate ROM
    int modeNum = ioData->linearizedData.modeNumber - 1;
    if (ioData->problem.alltype == ProblemData::_ROM_)
      com->fprintf(stderr, " ... Running Unsteady Linearized ROM w/no structural initializations\n");
    else if (ioData->linearizedData.initCond == LinearizedData::DISPLACEMENT)  {
      com->fprintf(stderr, " ... Initializing with Displacement Mode: %d\n", ioData->linearizedData.modeNumber);
      delU[modeNum] = ioData->linearizedData.amplification;
    }
    else if (ioData->linearizedData.initCond == LinearizedData::VELOCITY)  {
      com->fprintf(stderr, " ... Initializing with Velocity Mode: %d\n", ioData->linearizedData.modeNumber);
      delY[modeNum] = ioData->linearizedData.amplification;
    }
    else  {
      com->fprintf(stderr, "*** Error: Bad Initial Condition: %d\n", ioData->linearizedData.initCond);
      exit (-1);
    }

    // Form ROM operators
    VecSet<Vec<double> > ecVecs(nStrMode, nPodVecs);
    VecSet<Vec<double> > gVecs(nStrMode, nPodVecs);

    VecSet<Vec<double> > romOperator0(nPodVecs, nPodVecs);
    double *romOperator = new double[nPodVecs*nPodVecs];
    double *romOperator1 = new double[nPodVecs*nPodVecs];
    double *romOperator2 = new double[nPodVecs*nPodVecs];
    
    double t0 = modalTimer->getTime();
    constructROM2(romOperator, romOperator0, romOperator1, romOperator2, ecVecs, gVecs, podVecs, nPodVecs);
    modalTimer->addRomConstrTime(t0);

    modalTimer->getTime(); 
    timeIntegrateROM(romOperator, romOperator0, romOperator1, romOperator2, ecVecs, gVecs, podVecs, nSteps, nPodVecs, delU, delY, sdt);
    modalTimer->addRomTimeIntegTime(t0);

  }
}

//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::createPODInTime()  {

  if (nStrMode == 0) return;

  Timer *modalTimer = domain.getTimer();
  double t0 = modalTimer->getTime();
#ifdef DO_MODAL
  double sdt = ioData->linearizedData.stepsize;
  int nSteps = ioData->ts.maxIts;
  int numSnapsTot = 0;

  int iSnap = 0;

  VecSet<DistSVec<double, dim> > snaps(numSnapsTot, domain.getNodeDistInfo());

  // set up Snapshot Matrix
  int nModes = nStrMode;
  int numSnapsPerImpulse = nSteps;
  int snapFreq = nSteps / numSnapsPerImpulse;
  //char *suffix = "mode";
  numSnapsTot = 2*nModes*numSnapsPerImpulse;
  snaps.resize(numSnapsTot);
  com->fprintf(stderr, " ... Allocating for %d total snapshots\n", numSnapsTot);

  // Time Integration Loop for EC impulses
  int i;

  double *delU = new double[nModes];
  double *delY = new double[nModes];

  for (i = 0; i < nModes; i++)  {
    com->fprintf(stderr, " ... Impulsing mode %d\n", i);
    //delY[i] = sqrt(K[i]) * ioData->linearizedData.amplification;
    delY[i] = ioData->linearizedData.amplification;
    timeIntegrate(snaps, nSteps, snapFreq, delU, delY, iSnap, sdt);
    com->fprintf(stderr, " ... Completed Snapshot %d\n", iSnap);
    delY[i] = 0.0;

    delU[i] = ioData->linearizedData.amplification;
    timeIntegrate(snaps, nSteps, snapFreq, delU, delY, iSnap, sdt);
    com->fprintf(stderr, " ... Completed Snapshot %d\n", iSnap);
    delU[i] = 0.0;
  }
  t0 = modalTimer->getTime();

  int nPOD = ioData->linearizedData.numPOD;
  int nconv = 0;

  com->fprintf(stderr, " ... Forming Correlation Matrix\n");
  // allocate for upper half of sym. eigprob
  double *rVals = new double[numSnapsTot*(numSnapsTot+1)/2];
  for (i = 0; i < numSnapsTot; i++)
    for (int j = 0; j <= i; j++)
      rVals[(i+1)*i/2 + j] = snaps[j] * snaps[i];

  double tolerance = ioData->linearizedData.tolerance;

  ARdsSymMatrix<double> pod(numSnapsTot, rVals, 'U');

  com->fprintf(stderr, " ... Factoring Correlation Matrix\n");

  pod.FactorA();
  int ncv = 2*nPOD;
  if (ncv > numSnapsTot)  ncv = numSnapsTot-1;
  ARluSymStdEig<double> podEigProb(nPOD, pod, "LM", ncv, tolerance, 300*nPOD);

  com->fprintf(stderr, " ... Solving EigenProblem\n");
  nconv = podEigProb.FindEigenvectors();

  com->fprintf(stderr, " ... Got %d converged eigenvectors out of %d snaps\n", nconv, numSnapsTot);
  outputPODVectors(podEigProb, snaps, nPOD, numSnapsTot);

#else
  com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
  exit(-1);
#endif

  modalTimer->addMeshSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::timeIntegrate(VecSet<DistSVec<double, dim> > &snaps, 
                       int nSteps, int snapFreq, double *delU, double *delY, 
                       int &iSnap, double sdt, char *snapFile)  {

  Timer *modalTimer = domain.getTimer();
  double t0;

  DistSVec<double,dim> FF(domain.getNodeDistInfo());
  VarFcn *varFcn = new VarFcn(*ioData); 
  spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);

  // basic initializations
  DistSVec<double,3> deltmp(domain.getNodeDistInfo());
  DistSVec<double,dim> delW(domain.getNodeDistInfo());
  DistSVec<double,dim> delWtmp(domain.getNodeDistInfo());
  DistSVec<double,dim> rhs(domain.getNodeDistInfo());
  DistSVec<double,dim> rhsA(domain.getNodeDistInfo());
  DistSVec<double,dim> rhsB(domain.getNodeDistInfo());
  DistSVec<double,dim> rhsC(domain.getNodeDistInfo());
  DistSVec<double,dim> rhsD(domain.getNodeDistInfo());
  DistSVec<double,dim> delWnm1(domain.getNodeDistInfo());
  DistSVec<double,dim> delWint(domain.getNodeDistInfo());

  // Uncomment for forced oscillations
  double pi = 3.14159265358979;
  double freq = 2.0*pi*ioData->linearizedData.frequency;
  double dMax = ioData->linearizedData.amplification;

  int i;
  deltmp = 0.0;

  double sdt0 = sdt*dt0/dt;

  double r = dt/(2*dt0);
  double *prevU = new double[nStrMode];
  double *prevY = new double[nStrMode];
  for (i = 0; i < nStrMode; ++i)  {
    deltmp += (delU[i]+sdt0/2*delY[i])*mX[i];
    prevU[i] = delU[i];
    prevY[i] = delY[i];
  }
  
  deltmp += Xref;

  // Time Loop
  ksp->printParam();
  ksp2->printParam();
  ksp3->printParam(); 
 
  int printFreq = nSteps / 20;
  if (nSteps < 20)
    printFreq = 1;
  
  com->fprintf(stderr, " ... Doing Time Integration with refVals: F= %e, rho = %e, P=%e, M=%e\n", ioData->ref.rv.force, ioData->ref.density, ioData->ref.pressure, ioData->ref.mach); 

  // Init delW
  char *nlSolFile = 0;
  if (ioData->input.perturbed[0] == 0)
    nlSolFile = tInput->solutions;
  else  {
    int sp = strlen(ioData->input.prefix) + 1;
    nlSolFile = new char[sp + strlen(ioData->input.perturbed)];
    sprintf(nlSolFile, "%s%s", ioData->input.prefix, ioData->input.perturbed);
  }
  domain.readVectorFromFile(nlSolFile, 0, 0, delW);
  com->fprintf(stderr, " ... Read Perturbed solution: W = %e\n", delW.norm());

  tOutput->writeForcesToDisk(0, 0, 0, 0, 0.0, com->cpuNum(), tRestart->energy, deltmp, delW);
  tOutput->writeBinaryVectorsToDisk(false, 0, 0, deltmp, controlVol, delW, tState);

  delW -= Uref;
  delWnm1 = delW;
  delWint = delW;

  // Compute Reference Modal Force
  DistSVec<double,3> refNodalForce(domain.getNodeDistInfo());
  DistVec<double> Pin(domain.getFaceDistInfo());
  Vec<double> refModalF(nStrMode);
  com->fprintf(stderr, " ... Computing Ref Nodal Force w/pressure: %e\n", ioData->aero.pressure);
  Pin = ioData->aero.pressure;
  postOp->computeNodalForce(Xref, Uref, Pin, refNodalForce);
  
  CommPattern<double> *vpat = domain.getVecPat();
  domain.assemble(vpat, refNodalForce);

  for (i= 0; i < nStrMode; i++)
    refModalF[i] = mX[i]*refNodalForce;

  int cntr = 0;
  int cntp1;

  // every step corresponds to solving for 2x each member of the BDF scheme
  FF *= (-2.0); 
  for (int cnt = 0; cnt < nSteps; ++cnt) {

    cntp1 = cnt+1;

    t0 = modalTimer->getTime();  
    if (cnt == 0) {

      HOp2step1->apply(delWint,rhs);

      rhsA = 0.0;
      rhsB = 0.0;
      rhsC = 0.0;
      rhsD = 0.0;
      for (i = 0; i < nStrMode; ++i)  {
        rhsA += DX[i]*delU[i];
        rhsB += DX[i]*delY[i];
        rhsC += DE[i]*delY[i];

      }

      rhsA *= 2.0;
      rhsB *= sdt0;
      rhsC *= 2.0;

      rhs -= rhsA;
      rhs -= rhsB;
      rhs -= rhsC;

      rhs += FF;
      delWint = rhs;
      delWint /= controlVol;
      delWint *= dt0/2.0;

      delW = delWint;

    }
    else if (cnt == 1) {
 
      rhs = delW*(2.0+2.0/r) - delWnm1*(2.0*r/(1.0+r));

      rhs *= controlVol;

      rhs *= 1.0/dt0;

      rhsA = 0.0;
      rhsB = 0.0;
      rhsC = 0.0;
      rhsD = 0.0;
      for (i = 0; i < nStrMode; ++i)  {
         rhsA += DX[i]*delU[i];
         rhsB += DX[i]*delY[i];
         rhsC += DE[i]*delY[i];
         rhsD += DE[i]*prevY[i];
      }

      rhsA *= 2.0;
      rhsB *= sdt;
      rhsC *= 2.0*(1.0+r);
      rhsD *= 2.0*r;

      rhs -= rhsA;
      rhs -= rhsB;
      rhs -= rhsC;
      rhs += rhsD;

      delWnm1 = delW;

      rhs += FF;

      ksp2->solve(rhs, delW);
      modalTimer->addKspTime(t0);      

    }
    else if (cnt == 2) {
      rhs = delW*3.0 - delWnm1*(1.0/3.0);

      rhs *= controlVol;

      rhs *= 1.0/dt;

      rhsA = 0.0;
      rhsB = 0.0;
      rhsC = 0.0;
      rhsD = 0.0;
      for (i = 0; i < nStrMode; ++i)  {
         rhsA += DX[i]*delU[i];
         rhsB += DX[i]*delY[i];
         rhsC += DE[i]*delY[i];
         rhsD += DE[i]*prevY[i];

      }

      rhsA *= 2.0;
      rhsB *= sdt;
      rhsC *= 3.0;

      rhs -= rhsA;
      rhs -= rhsB;
      rhs -= rhsC;
      rhs += rhsD;

      delWnm1 = delW;

      rhs += FF;

      ksp3->solve(rhs, delW);
      modalTimer->addKspTime(t0);
    }

    else {
      rhs = (delW*4.0) - delWnm1;

      rhs *= controlVol;

      rhs *= 1.0/dt;

      rhsA = 0.0;  rhsB = 0.0; rhsC = 0.0; rhsD = 0.0;
      for (i = 0; i < nStrMode; ++i)  {
        rhsA += DX[i]*delU[i];
        rhsB += DX[i]*delY[i];
        rhsC += DE[i]*delY[i];
        rhsD += DE[i]*prevY[i];
      }
      rhsA *= 2.0;
      rhsB *= sdt;
      rhsC *= 3.0;
      rhs -= rhsA;
      rhs -= rhsB;
      rhs -= rhsC;
      rhs += rhsD;
      delWnm1 = delW;

      rhs += FF;
      ksp->solve(rhs, delW);
      modalTimer->addKspTime(t0);
    }

    // for forced oscillations          
    if (ioData->linearizedData.type == LinearizedData::FORCED)  {
      if (cnt <= 1) {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*(cnt)*sdt0)*dMax;   
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*(cnt)*sdt0)*dMax;
      } 
      else {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*((cnt-1)*sdt+sdt0))*dMax;
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*((cnt-1)*sdt+sdt0))*dMax;  
      }
    }

    if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {
    
      for (int kk = 0; kk < nStrMode; kk++)  {
        prevU[kk] = delU[kk];
        prevY[kk] = delY[kk];
      }
      t0 = modalTimer->getTime();
     
      if (cnt == 0) 
        computeModalDispStep1(sdt0, deltmp, delW, delU, delY, refModalF);
      else
        computeModalDisp(sdt, deltmp, delW, delU, delY, refModalF);
  
      modalTimer->addStructUpdTime(t0);

    }
  
    // compute updated position
    deltmp = 0.0;
    for (i = 0; i < nStrMode; ++i) {
        deltmp += (delU[i]+sdt/2*delY[i])*mX[i];
    }
    deltmp += Xref;

    if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
      if (cnt % snapFreq == 0 || cnt == 0) {
        snaps[iSnap] = delW;
        if (snapFile) domain.writeVectorToFile(snapFile, iSnap, 0, snaps[iSnap]);
        iSnap++;
      }
    }

    delWtmp = Uref + delW;
    int ierr = domain.checkSolution(varFcn, delWtmp);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, cnt*dt+dt0, com->cpuNum(), tRestart->energy, deltmp, delWtmp);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, cnt*dt+dt0, deltmp, controlVol, delWtmp, tState);


    if (ierr)  {
      com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
      cntr++;
      if (cntr > 4)  {
        com->barrier();
        exit(-1);
      }
    }

    if (cnt % 20 == 1)
      com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, (cnt-1)*sdt+sdt0);   
 
    if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_ || 
        ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
      for (i = 0; i < nStrMode; ++i) {
        delU[i] = 0.0;
        delY[i] = 0.0;
      }
    }
  }
}

//----------------------------------------------------------------------------------

template<int dim>
void
ModalSolver<dim>::timeIntegrateROM(double *romOp, VecSet<Vec<double> > &romOp0, double *romOp1, double *romOp2, VecSet<Vec<double> > &ecMat, VecSet<Vec<double> > &gMat, VecSet<DistSVec<double, dim> > &podVecs, int nSteps, int nPodVecs, double *delU, double *delY, double sdt)  {

#ifdef DO_MODAL
  VarFcn *varFcn = new VarFcn(*ioData); 

  // basic initializations
  DistSVec<double,3> deltmp(domain.getNodeDistInfo());
  DistSVec<double, dim> delWFull(domain.getNodeDistInfo());
  Vec<double> delWRom(nPodVecs);
  Vec<double> delWRomTemp(nPodVecs);
  Vec<double> delWnm1Rom(nPodVecs);
  Vec<double> rhs(nPodVecs);
  Vec<double> rhsTemp(nPodVecs);
  Vec<double> modalF(nStrMode);
  modalF = 0.0;

  delWRom = 0.0;
  deltmp = 0.0;

  double sdt0 = sdt*dt0/dt;

  double r = dt/(2.0*dt0);

  double *prevU = new double[nStrMode];
  double *prevY = new double[nStrMode];

  int i;
  for (i = 0; i < nStrMode; ++i)
    deltmp += (delU[i]+0.5*sdt0*delY[i])*mX[i];

  deltmp += Xref;

  // Init delW
  char *nlSolFile = 0;
  if (ioData->input.perturbed[0] == 0)
    nlSolFile = tInput->solutions;
  else  {
    int sp = strlen(ioData->input.prefix) + 1;
    nlSolFile = new char[sp + strlen(ioData->input.perturbed)+1];
    sprintf(nlSolFile, "%s%s", ioData->input.prefix, ioData->input.perturbed);
  }
  domain.readVectorFromFile(nlSolFile, 0, 0, delWFull);
  com->fprintf(stderr, " ... Read Perturbed solution: W = %e, Uref = %e\n", delWFull.norm(), Uref.norm());

  delWFull -= Uref;

  // construct initial delWRom
  for (i = 0; i < nPodVecs; i++)
    delWRom[i] = podVecs[i] * delWFull;

  delWnm1Rom = delWRom;

  // Compute Reference Modal Force
  DistSVec<double,3> refNodalForce(domain.getNodeDistInfo());
  DistVec<double> Pin(domain.getFaceDistInfo());
  Vec<double> refModalF(nStrMode);
  com->fprintf(stderr, " ... Computing Ref Nodal Force w/pressure: %e\n", ioData->aero.pressure);
  Pin = ioData->aero.pressure;
  postOp->computeNodalForce(Xref, Uref, Pin, refNodalForce);

  CommPattern<double> *vpat = domain.getVecPat();
  domain.assemble(vpat,refNodalForce);

  // Form ROM operators
  VecSet<Vec<double> > ecOpMat(nStrMode, nPodVecs);
  VecSet<Vec<double> > gOpMat(nStrMode, nPodVecs);
  VecSet<Vec<double> > ecOpMat1(nStrMode, nPodVecs);
  VecSet<Vec<double> > gOpMat1(nStrMode, nPodVecs);
  VecSet<Vec<double> > ecOpMat2(nStrMode, nPodVecs);
  VecSet<Vec<double> > gOpMat2(nStrMode, nPodVecs);

  VecSet<Vec<double> > romOpMinus(nPodVecs, nPodVecs);
  VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);
  VecSet<Vec<double> > romOperator1(nPodVecs, nPodVecs);
  VecSet<Vec<double> > romOperator2(nPodVecs, nPodVecs);

  int iVec;
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    romOpMinus[iVec] = 0.0;
    romOpMinus[iVec][iVec] = 1.0;
  } 


  ARdsNonSymMatrix<double, double> romOpPlus(nPodVecs, romOp);
  romOpPlus.FactorA();
  ARdsNonSymMatrix<double, double> romOpPlus1(nPodVecs, romOp1);
  romOpPlus1.FactorA();
  ARdsNonSymMatrix<double, double> romOpPlus2(nPodVecs, romOp2);
  romOpPlus2.FactorA();

  com->fprintf(stderr, " ... Factored Matrix\n");

  com->fprintf(stderr, " ... Forming Reduced ROM Operator\n");

  for (iVec = 0; iVec < nPodVecs; iVec++) {
    romOpPlus.MultInvv(romOpMinus[iVec].data(), romOperator[iVec].data());
    romOpPlus1.MultInvv(romOpMinus[iVec].data(), romOperator1[iVec].data());
    romOpPlus2.MultInvv(romOpMinus[iVec].data(), romOperator2[iVec].data());
  }

  for (iVec = 0; iVec < nStrMode; iVec++)  {
    romOpPlus.MultInvv(ecMat[iVec].data(), ecOpMat[iVec].data());
    romOpPlus.MultInvv(gMat[iVec].data(), gOpMat[iVec].data());
    romOpPlus1.MultInvv(ecMat[iVec].data(), ecOpMat1[iVec].data());
    romOpPlus1.MultInvv(gMat[iVec].data(), gOpMat1[iVec].data());
    romOpPlus2.MultInvv(ecMat[iVec].data(), ecOpMat2[iVec].data());
    romOpPlus2.MultInvv(gMat[iVec].data(), gOpMat2[iVec].data());
  }
  
  for (i = 0; i < nStrMode; i++)
    refModalF[i] = mX[i]*refNodalForce;

  com->fprintf(stderr, " ... RefModalF = %e\n", refModalF.norm());
  com->fprintf(stderr, " ... Initial delWRom: %e\n", delWRom.norm());

  int cntp1;
  delWFull = 0.0;
  for (i = 0; i < nPodVecs; ++i)
    delWFull += podVecs[i] * delWRom[i];
  delWFull += Uref;
  tOutput->writeForcesToDisk(0, 0, 0, 0, 0.0, com->cpuNum(), tRestart->energy, deltmp, delWFull);
  tOutput->writeBinaryVectorsToDisk(false, 0, 0, deltmp, controlVol, delWFull, tState);

  for (int cnt = 0; cnt < nSteps; ++cnt) {

    cntp1 = cnt+1;

    if (cnt == 0) {

      rhs = 0.0;
      rhsTemp = 0.0;

      delWRomTemp = delWRom;

      for (i = 0; i < nPodVecs; ++i)
        rhs += delWRomTemp[i]*romOp0[i];

      for (i = 0; i < nStrMode; ++i)
        rhsTemp += - delU[i]*gMat[i] - 0.5*sdt0*delY[i]*gMat[i] - delY[i]*ecMat[i];
    
      // compute new delWRom
      delWnm1Rom = delWRom;
      delWRom = dt0 *(rhs + rhsTemp);

      // project soltn into full space
      delWFull = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWFull += podVecs[i] * delWRom[i];

    }

    else if (cnt == 1){
     
      rhs = 0.0;
      rhsTemp = 0.0;

      delWRomTemp = (2.0+2.0/r)*delWRom;
      delWRomTemp -= 2.0*r/(1.0+r)*delWnm1Rom;

      delWRomTemp *= 1.0/dt0;
      for (i = 0; i < nPodVecs; ++i)
        rhs += delWRomTemp[i]*romOperator1[i];

      for (i = 0; i < nStrMode; ++i) 
        rhsTemp += - 2.0*delU[i]*gOpMat1[i] - sdt*delY[i]*gOpMat1[i] - 2.0*(1.0+r)*delY[i]*ecOpMat1[i] + 2.0*r*prevY[i]*ecOpMat1[i];

      // compute new delWRom
      delWnm1Rom = delWRom;
      delWRom = rhs + rhsTemp;

      // project soltn into full space
      delWFull = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWFull += podVecs[i] * delWRom[i];
    }

    else if (cnt == 2){

      rhs = 0.0;
      rhsTemp = 0.0;

      delWRomTemp = 3.0*delWRom;
      delWRomTemp -= 1.0/3.0*delWnm1Rom;

      delWRomTemp *= 1.0/dt;
      for (i = 0; i < nPodVecs; ++i)
        rhs += delWRomTemp[i]*romOperator2[i];

      for (i = 0; i < nStrMode; ++i)
        rhsTemp += - 2.0*delU[i]*gOpMat2[i] - sdt*delY[i]*gOpMat2[i] - 3.0*delY[i]*ecOpMat2[i] + prevY[i]*ecOpMat2[i];

      // compute new delWRom
      delWnm1Rom = delWRom;
      delWRom = rhs + rhsTemp;

      // project soltn into full space
      delWFull = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWFull += podVecs[i] * delWRom[i];

    }
        
    else {

      rhs = 0.0;
      rhsTemp = 0.0;

      delWRomTemp = 4.0*delWRom;
      delWRomTemp -= delWnm1Rom;

      delWRomTemp *= 1.0/dt; 
      for (i = 0; i < nPodVecs; ++i)  
        rhs += delWRomTemp[i]*romOperator[i];

      for (i = 0; i < nStrMode; ++i) 
        rhsTemp += - 2.0*delU[i]*gOpMat[i] - sdt*delY[i]*gOpMat[i] - 3.0*delY[i]*ecOpMat[i] + prevY[i]*ecOpMat[i];

      // compute new delWRom
      delWnm1Rom = delWRom;
      delWRom = rhs + rhsTemp;
      // project soltn into full space
      delWFull = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWFull += podVecs[i] * delWRom[i];

    }
  
    // Output the reduced order vector
    if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {
      for (int kk = 0; kk < nStrMode; kk++)  {
        prevU[kk] = delU[kk];
        prevY[kk] = delY[kk];
      }
  
      if (cnt == 0)
        computeModalDispStep1(sdt0,deltmp, delWFull, delU, delY, refModalF);
      else
        computeModalDisp(sdt, deltmp, delWFull, delU, delY, refModalF);
    }

    // compute Cl, Cm
    deltmp = 0.0;
    for (i = 0; i < nStrMode; ++i)
      deltmp += (delU[i]+sdt/2*delY[i])*mX[i];

    deltmp += Xref;
    delWFull += Uref;

    int ierr = domain.checkSolution(varFcn, delWFull);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, cnt*dt+dt0, com->cpuNum(), tRestart->energy, deltmp, delWFull);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, cnt*dt+dt0, deltmp, controlVol, delWFull, tState);

    if (ierr)  {
      com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
      exit(-1);
    }

    if (cnt % 20 == 0)
      com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);

  }
#else
  com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
  exit(-1);

#endif
}

//--------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::preProcess()  {

  // setup solvers
  VarFcn *varFcn = new VarFcn(*ioData);
  geoState = new DistGeoState(*ioData, &domain);
  geoState->setup1(tInput->positions, &Xref, &controlVol);
  bcData = new DistBcDataEuler<dim>(*ioData, varFcn, &domain, Xref);

  spaceOp = new SpaceOperator<dim>(*ioData, varFcn, bcData, geoState, &domain);
  postOp = new PostOperator<dim> (*ioData, varFcn, bcData, geoState, &domain);

  // Temporal operator contains Uref for us
  tState = new DistTimeState<dim>(*ioData, spaceOp, varFcn, &domain);
  tState->setup(tInput->solutions,  Xref, bcData->getInletBoundaryVector(), Uref, *ioData);

  RefVal *refVal = new RefVal(ioData->ref.rv);
  tOutput = new TsOutput<dim>(*ioData, refVal, &domain, postOp);
  tRestart = new TsRestart(*ioData, refVal);

  HOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  HOp2 = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  HOp2step1 = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  HOpstep2 = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  HOpstep3 = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);

  // Stuff for solver
  PcData pcData = ioData->ts.implicit.newton.ksp.ns.pc;
  KspData &kspData = ioData->ts.implicit.newton.ksp.ns;

  if (pcData.type == PcData::IDENTITY)
    pc = new IdentityPrec<dim, double>();
  else if (pcData.type == PcData::JACOBI)
    pc = new JacobiPrec<PreScalar,dim>(DiagMat<PreScalar,dim>::DENSE, &domain);
  else if (pcData.type == PcData::AS || pcData.type == PcData::RAS ||
           pcData.type == PcData::ASH || pcData.type == PcData::AAS)
    pc = new IluPrec<PreScalar,dim>(pcData, &domain);


  if (kspData.type == KspData::GMRES) {
    ksp = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com);
    ksp2 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com);
    ksp3 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com);
    }
  else if (kspData.type == KspData::GCR) {
    ksp = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com);
    ksp2 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com);
    ksp3 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com);
  }
                                                        
  if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
    HOpC = new MatVecProdH2<5,bcomp,5>(*ioData,  varFcn, tState, spaceOp, &domain);

                                                        
    pcComplex = new IluPrec<bcomp ,dim, bcomp>(pcData, &domain);
   if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
                                                        
                                                                                                                 
     kspCompGcr = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim, 5>,
                    KspPrec<dim, bcomp>, Communicator, bcomp>
                    ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com);
    }
    else {
                                                        
      if (kspData.type == KspData::GMRES)
                                                        
        kspComp = new GmresSolver<DistSVec<bcomp,dim>, MatVecProd<dim, 5>,
                      KspPrec<dim, bcomp>, Communicator, bcomp>
                     ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com);
      else if (kspData.type == KspData::GCR)
                                                        
        kspComp = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim, 5>,
                      KspPrec<dim, bcomp>, Communicator, bcomp>
                     ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com);
                                                        
                                                        
    }
  }


  geoState->setup2(tState->getData());

  //Setup Time state
  double dummyTime = 0.0;
  int dummySubCycles = 1;
  tState->computeTimeStep(1.0, &dummyTime, &dummySubCycles, *geoState, Xref, controlVol, Uref);
  dt = ioData->linearizedData.stepsize/ioData->ref.rv.time;    //ts.timestep/  (ref.length/velocity)
  
  //time step for the first iteration (Forward Euler) has to be much smaller : here sdt0 = sdt^2
  
  dt0 = (ioData->linearizedData.stepsize)*(ioData->linearizedData.stepsize)/ioData->ref.rv.time; 
  com->fprintf(stderr, " ... Running Modal Solver w/%d struct modes, fluid step: %f, RefTime = %f (%f/%f)\n", nStrMode, dt, ioData->ref.rv.time, ioData->ref.rv.length, ioData->ref.rv.velocity);

  geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);
  geoState->update(Xref, controlVol);
  geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);

  // F is not used
  DistSVec<double,dim> FF(domain.getNodeDistInfo());
  spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);
  FF /= controlVol;

  com->fprintf(stderr, " ... Norm FF: %e\n", FF.norm());
  if (FF.norm() > 1e-4)
    com->fprintf(stderr, " *** WARNING: Check Steady State Convergence for this Mach Number ***\n");

  // Compute Derivatives***********************************
  //***Variable Initializations***
  DistSVec<double,3> Xnp1(domain.getNodeDistInfo());
  DistVec<double> CV1(domain.getNodeDistInfo());
  DistSVec<double,dim> F1(domain.getNodeDistInfo());
  DistVec<double> CV2(domain.getNodeDistInfo());
  DistSVec<double,dim> F2(domain.getNodeDistInfo());
  VecSet< DistVec<double> > DA(nStrMode, domain.getNodeDistInfo());
  DistSVec<double,dim> E(domain.getNodeDistInfo());
  VecSet< DistSVec<double,dim> > DV(nStrMode, domain.getNodeDistInfo());

  DX.resize(nStrMode);
  DE.resize(nStrMode);

  double eps = ioData->linearizedData.eps;
  double alpha = dt/ioData->linearizedData.eps2;
  double eps2 = eps*alpha;
  com->fprintf(stderr, "Alpha = %f\n", alpha);

  // Loop over the modes
  
  for(int ic=0;ic<nStrMode;++ic) {

    //First DFDX & DADX***********************************
    //***Computing F(X1,V=0)
    Xnp1 = Xref -  eps*mX[ic]; //1.0*Xref;
    
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    geoState->update(Xnp1, CV1);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    
    spaceOp->computeResidual(Xnp1, CV1, Uref, F1, tState);
    spaceOp->applyBCsToResidual(Uref, F1);

    //***Computing F(X2,V=0)
    Xnp1 = Xref + eps*(mX[ic]);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV2);
    geoState->update(Xnp1, CV2);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV2);

    spaceOp->computeResidual(Xnp1, CV2, Uref, F2, tState);
    spaceOp->applyBCsToResidual(Uref, F2);

    DX[ic] = 1.0/(2.0*eps)*(F2 - F1);
    DA[ic] = 1.0/(2.0*eps)*(CV2 - CV1);
    com->fprintf(stderr, " ... Norm DX, Mode %i: %e\n",ic, DX[ic].norm());
    com->fprintf(stderr, " ... Norm DA, Mode %i: %e\n",ic, DA[ic].norm());

    // Now Compute w*dA/dX***********************************
    int numLocSub = Uref.numLocSub();

    #pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = DA[ic].subData(iSub);
      double (*r)[dim] = Uref.subData(iSub);
      double (*ans)[dim] = E.subData(iSub);

      for (int i=0; i<DA[ic].subSize(iSub); ++i) {
        for (int j=0; j<dim; ++j)
          ans[i][j] = cv[i]*r[i][j];
      }
    }
    
    // need to mult. by time to maintain a dimensional vel. (conversion from adim time to dim time)
    DE[ic] = ioData->ref.rv.time*E;
    
    // Then DFDXdot*****************************************
    // Computing F(X0,V1) 
        
    Xnp1 = Xref - eps2*(mX[ic]);
    
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    geoState->update(Xnp1, CV1);
    
    Xnp1 = Xref + eps2*(mX[ic]);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

    //spaceOp->computeResidual(Xnp1, CV1, Uref, F1, tState);
    spaceOp->computeResidual(Xref, controlVol, Uref, F1, tState);
    spaceOp->applyBCsToResidual(Uref, F1);
 
    // Computing F(X0, -V1)
    Xnp1 = Xref + eps2*(mX[ic]);
  
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    geoState->update(Xnp1, CV1);

    Xnp1 = Xref - eps2*(mX[ic]);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

    //spaceOp->computeResidual(Xnp1, CV1, Uref, F2, tState);
    spaceOp->computeResidual(Xref, controlVol, Uref, F2, tState);
    spaceOp->applyBCsToResidual(Uref, F2);


    // no dt here because the spatial computation will divide the
    // given displacement change by dt to define its velocity.
    // Thus with Xnp1 = X0+dt*eps*Xm and Xn = X0-dt*eps*Xm: V_n+1/2 = 2*eps*Xm
    
    DV[ic] = (1.0/(2*eps))*(F1-F2);
    DV[ic] *= ioData->ref.rv.time;

    com->fprintf(stderr, " ... Norm DV, Mode %i: %e\n",ic, DV[ic].norm());

    // create E+C
    DE[ic] += DV[ic];
    com->fprintf(stderr, " ... Norm DE, Mode %i: %e\n\n",ic, DE[ic].norm());
  }

  // Now setup H matrices

  double delt0;
    delt0 = dt0;

  double r = dt/(2*delt0);
  // ***This is (c1*A+c2*dt*H)
  HOp->evaluate(3, Xref, controlVol, Uref, FF,0.0);

  //***This is (c1*A-c2*dt*H)
  HOp2->evaluate2(0, Xref, controlVol, Uref, FF);
  
  tState->setGlobalTimeStep(delt0);
  HOp2step1->evaluate(5, Xref, controlVol, Uref, FF,0.0);

  tState->setGlobalTimeStep(delt0);
  HOpstep2->evaluate(6, Xref, controlVol, Uref, FF,r);

  tState->setGlobalTimeStep(dt);
  HOpstep3->evaluate(7, Xref, controlVol, Uref, FF,0.0);

  //set up preconditioner, GMRES Solver
  DistMat<PreScalar,dim> *_pc = dynamic_cast<DistMat<PreScalar,dim> *>(pc);

  if (_pc) {
    spaceOp->computeH1(Xref, controlVol, Uref, *_pc);
    tState->addToH1(controlVol, *_pc);
    spaceOp->applyBCsToJacobian(Uref, *_pc);
  }

  pc->setup();
  ksp->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);
  ksp2->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);
  ksp3->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);

  tOutput->openAsciiFiles();
}

//-------------------------------------------------------------------------------
/*
template<int dim>
void ModalSolver<dim>::constructROM(VecSet<Vec<double> > &romOperator,
                       VecSet<Vec<double> > &ecVecs, VecSet<Vec<double> > &gVecs,
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  com->fprintf(stderr, " ... Forming Reduced Integrator \n");
  int iVec;

  // form H ROM op
  DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
  DistSVec<double, dim> tmpRomVec(domain.getNodeDistInfo());
  ksp->printParam();
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    HOp2->apply(podVecs[iVec], tmpVec);
    ksp->solve(tmpVec, tmpRomVec);
    for (int jVec = 0; jVec < nPodVecs; jVec++)
      romOperator[iVec][jVec] = podVecs[jVec] * tmpRomVec;
  }
  com->fprintf(stderr, " ... Computed Rom Operator\n");

  // form coupling matrix ROMs
  DistSVec<double,dim> tmpECvec(domain.getNodeDistInfo());
  DistSVec<double,dim> tmpGvec(domain.getNodeDistInfo());
  for (iVec = 0; iVec < nStrMode; iVec++)  {

    tmpECvec = 0.0; 
    tmpGvec = 0.0;
    // form (A+.5*dt*H)^-1 B
    ksp->solve(DE[iVec], tmpECvec);
    ksp->solve(DX[iVec], tmpGvec);

    for (int jVec = 0; jVec < nPodVecs; jVec++)  {
      ecVecs[iVec][jVec] = podVecs[jVec] * tmpECvec;
      gVecs[iVec][jVec] = podVecs[jVec] * tmpGvec;
    }
  }
  com->fprintf(stderr, " ... Computed Coupling Rom Vectors\n");

}
*/
//-------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::constructROM2(double *romOpPlusVals, VecSet<Vec<double> > &romOperator0, double *romOpPlusVals1, double *romOpPlusVals2,
                       VecSet<Vec<double> > &ecVecs, VecSet<Vec<double> > &gVecs,
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {
#ifdef DO_MODAL
  com->fprintf(stderr, " ... Forming Reduced System Directly \n");

  // form fluid ROM
  DistSVec<double,dim> FF(domain.getNodeDistInfo());

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcn(*ioData); 
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  onlyHOp->evalH(0, Xref, controlVol, Uref);

  double r = dt/(2.0*dt0);

  double invDt0 = 1.0/dt0;
  double invDt3 = 3.0/dt;
  double invDtp = 8.0/(3.0*dt);
  double invDtr = (4.0*r+2.0)/(r*(r+1.0)*dt0); 
  int iVec, jVec;

  // form H ROM op
  DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
  

  double romVal;
  com->barrier();
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    onlyHOp->apply(podVecs[iVec], tmpVec);
    tmpVec /= controlVol;
    for (jVec = 0; jVec < nPodVecs; jVec++)  {
      romVal = podVecs[jVec] * tmpVec;
      romOpPlusVals[iVec*nPodVecs+jVec] = 2.0*romVal;
      romOpPlusVals1[iVec*nPodVecs+jVec] = 2.0*romVal;
      romOpPlusVals2[iVec*nPodVecs+jVec] = 2.0*romVal;
      romOperator0[iVec][jVec] = -romVal;
    }

    romOpPlusVals[iVec*nPodVecs+iVec] += invDt3;
    romOpPlusVals1[iVec*nPodVecs+iVec] += invDtr;
    romOpPlusVals2[iVec*nPodVecs+iVec] += invDtp;
    romOperator0[iVec][iVec] += invDt0;

    //spaceOp has to be reset because it has been modified by the apply function
    DistSVec<double,dim> FF(domain.getNodeDistInfo());
    spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);  
  }
  
  // form coupling matrix ROMs
  DistSVec<double,dim> tmpECvec(domain.getNodeDistInfo());
  DistSVec<double,dim> tmpGvec(domain.getNodeDistInfo());

  Vec<double> tmpECrom(nPodVecs);
  Vec<double> tmpGrom(nPodVecs);

  for (iVec = 0; iVec < nStrMode; iVec++)  {

    tmpECvec = DE[iVec]; 
    tmpGvec = DX[iVec];
    
    tmpECvec /= controlVol;
    tmpGvec /= controlVol;
 
    for (jVec = 0; jVec < nPodVecs; jVec++)  {
      tmpECrom[jVec] = podVecs[jVec] * tmpECvec;
      tmpGrom[jVec] = podVecs[jVec] * tmpGvec;
    }

    ecVecs[iVec] =  tmpECrom;
    gVecs[iVec] = tmpGrom;

  }
  com->fprintf(stderr, " ... Computed Coupling Rom Vectors\n");

#else
  com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
  exit(-1);

#endif

}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::computeModalDisp(double sdt, DistSVec<double, 3> &xPos, DistSVec<double, dim> &delW, double *delU, double *delY, Vec<double> &refModalF)  {

  Vec<double> modalF(nStrMode);
  postOp->computeForceDerivs(xPos, Uref, delW, modalF, mX);
  modalF += refModalF;
  modalF *= ioData->ref.rv.force;
  // form struct rhs and solve u dof
  double sdt2 = sdt*sdt;
  double prevU, prevY;

  double *srhs = new double[nStrMode];

  for (int i = 0; i < nStrMode; i++)  {

    srhs[i] = (1.0 - 0.25*sdt2*K[i])*delU[i] + sdt*delY[i] + 0.5*sdt2*modalF[i];
    prevU = delU[i];
    delU[i] = srhs[i] / (1+0.25*sdt2*K[i]);
    prevY = delY[i];
    delY[i] = 2.0 * (delU[i] - prevU)/sdt - prevY;
  }

}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::computeModalDispStep1(double sdt, DistSVec<double, 3> &xPos, DistSVec<double, dim> &delW, double *delU, double *delY, Vec<double> &refModalF)  {

  Vec<double> modalF(nStrMode);
  postOp->computeForceDerivs(xPos, Uref, delW, modalF, mX);
  modalF += refModalF;
  modalF *= ioData->ref.rv.force;

  // form struct rhs and solve u dof
  double sdt2 = sdt*sdt;
  double prevU, prevY;

  double *srhs = new double[nStrMode];

  for (int i = 0; i < nStrMode; i++)  {

    srhs[i] = delU[i] + sdt*delY[i] + 0.5*sdt2*modalF[i];
    prevU = delU[i];
    delU[i] = srhs[i] / (1+0.5*sdt2*K[i]);
    prevY = delY[i];
    delY[i] = 2.0 * (delU[i] - prevU)/sdt - prevY;
  }
}


//-----------------------------------------------------------------------------

template <int dim>
void ModalSolver<dim>::constructPOD()  {

  Timer *modalTimer = domain.getTimer();
  double t0;

  // Initial Conditions
  VecSet<DistSVec<bcomp,dim> > prevW(nStrMode, domain.getNodeDistInfo());

  int i;
  for (i = 0; i < nStrMode; ++i)
    prevW[i] = 0.0;

  int nSteps = ioData->ts.maxIts;

  double refLength = ioData->linearizedData.refLength;

  if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {

    // adjust (E+C) matrices for freq. domain
    double invT = 1.0/ioData->ref.rv.time;
    double refConst = ioData->ref.rv.velocity/refLength;
    for (int iMode = 0; iMode < nStrMode; iMode++)
      DE[iMode] *= invT;

    com->fprintf(stderr, " ... Adjusting E+C with vel/length: %e (%e)\n", invT, refConst);

    DistSVec<double,dim> FF(domain.getNodeDistInfo());
    DistMat<bcomp,dim> *_pcC = dynamic_cast<DistMat<bcomp,dim> *>(pcComplex);
    int iSnap = 0;
    
    int nSteps = ioData->ts.maxIts;
    int nSnaps = (2*nSteps+1)*nStrMode;
    double deltaFreq = ioData->linearizedData.freqStep;
    int nSnapsCoarse = 0;

    int *stepParam = new int[4];                                                               
    stepParam[0] = nSteps;
    stepParam[2] = nSnaps;
    nPadeDeriv = -1;

    com->fprintf(stderr, " ... Running POD from reduced freq. 0.0 using %d steps with refLength %f\n", nSteps, refLength);

    VecSet<DistSVec<double, dim> > snaps(nSnaps, domain.getNodeDistInfo());

    int nStepsCoarse = 0;
    double *coarseFreq = new double[11];
    int L, M, nPoints;

    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)  {
      L = ioData->linearizedData.pade.degNum;
      M = ioData->linearizedData.pade.degDen;
      nPoints = ioData->linearizedData.pade.nPoints;
      if (ioData->linearizedData.pade.freq[0] >= 0.0) coarseFreq[nStepsCoarse++] = ioData->linearizedData.pade.freq[0];
      for (i=1; i< ioData->linearizedData.pade.num; i++) 
        if (ioData->linearizedData.pade.freq[i] > 0.0) coarseFreq[nStepsCoarse++] = ioData->linearizedData.pade.freq[i];
      nPadeDeriv = int(ceil((L+M+1.0)/((double)nPoints))-1.0);
      nSnapsCoarse = nStepsCoarse*2*(nPadeDeriv+1)*nStrMode;
      if (coarseFreq[0] == 0.0)
        nSnapsCoarse -= nStrMode; // No imaginary part for the snapshot at k=0.0
      nSteps = nStepsCoarse-1;
      com->fprintf(stderr, " ... Computing the POD and their derivatives at %d coarse freq.\n", nStepsCoarse);                                                
    }
    VecSet<DistSVec<bcomp,dim> > prevDerivW(nStrMode*(nPadeDeriv+1), domain.getNodeDistInfo());
    if (nPadeDeriv+1 > 0) {
      for (i = 0; i < nStrMode*(nPadeDeriv+1); ++i)
        prevDerivW[i] = 0.0;
    }
    stepParam[1] = nStepsCoarse;
    stepParam[3] = nSnapsCoarse;
    VecSet<DistSVec<double, dim> > snapsCoarse(nSnapsCoarse, domain.getNodeDistInfo());
    double kFreq;
    for (i = 0; i <= nSteps; i++)  {

      if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)
        kFreq = coarseFreq[i];
      else
        kFreq = i*deltaFreq;
      double hzFreq = kFreq*invT/6.283185307; 
      bcomp shift(0.0, kFreq);
      HOpC->evaluate(0, Xref, controlVol, Uref, FF, shift);
      com->fprintf(stderr, " ... Shifting by  k = %f(%f Hz) i\n", kFreq, hzFreq);
      if (_pcC) {
        spaceOp->computeH1(Xref, controlVol, Uref, *_pcC);
        tState->addToH1(controlVol, *_pcC, shift);
        spaceOp->applyBCsToJacobian(Uref, *_pcC);
      }
      pcComplex->setup();

      //t0 = modalTimer->getTime();
      if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
        freqIntegrateMultipleRhs(snapsCoarse, kFreq, iSnap, prevDerivW);
      }
      
      else
        freqIntegrate(snaps, kFreq, iSnap, prevW);
      //modalTimer->addSnapsLinSolvTime(t0); 

    }

    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)  {
      com->fprintf(stderr, " ... Computing the Pade Reconstruction\n");
      nSteps = ioData->ts.maxIts;
      t0 = modalTimer->getTime();
      domain.padeReconstruction(snapsCoarse, snaps, stepParam, coarseFreq, deltaFreq, nStrMode, L, M, nPoints);
      modalTimer->addPadeReconstrTime(t0);
    }                                                          
    delete[] coarseFreq;
    delete[] stepParam;
    makeFreqPOD(snaps, nSnaps);
  }
  else
    createPODInTime();

}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::freqIntegrate(VecSet<DistSVec<double, dim> >&snaps,
                       double kFreq, int &iSnap,
                       VecSet<DistSVec<bcomp,dim> > &prevW)  {
  
  Timer *modalTimer = domain.getTimer();
  double t0;
                                                      
  // loop over modes
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp zeroComp(0.0,0.0);
  bcomp kImag(0.0, kFreq);
                                                        
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  //DistSVec<double, dim> delWReal(domain.getNodeDistInfo());
  Vec3D x0(0.0, 0.0, 0.0);
                                                        
  Vec<double> modalF(nStrMode);

  rhs = oneReal*DX[0] + kImag*DE[0];


  t0 = modalTimer->getTime();
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->setup(1, 40, rhs);
    kspCompGcr->printParam();
    kspCompGcr->numCalcVec = 0;
  }
  else {
    kspComp->setup(1, 40, rhs);
    kspComp->printParam();
  }
  VecSet<Vec<bcomp> > aeroOp(nStrMode, nStrMode);
                                                        
  int iMode;
  for (iMode = 0; iMode < nStrMode; iMode++)  {
                                                        
    // form rhs
    rhs = -oneReal*DX[iMode] - kImag*DE[iMode];
    HOpC->applyT(rhs, delW);                                                      
    // solve
    com->fprintf(stderr, " ... Solving for mode %d, w/rhs norm %e\n", iMode+1, rhs.norm());
                                                        
    delW = prevW[iMode];
    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
      kspCompGcr->solve(rhs, delW);
    }
    else {
      kspComp->solve(rhs, delW);
    }
    if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
      snaps[iSnap++].getReal(delW);
                                                        
      if (kFreq != 0.0)
        snaps[iSnap++].getImag(delW);
    }
                                                        
    prevW[iMode] = delW;
                                                        
  }
  modalTimer->addSnapsLinSolvTime(t0);
}

//-----------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::freqIntegrateMultipleRhs(VecSet<DistSVec<double, dim> >&snapsCoarse,
                       double kFreq, int &iSnap,
                       VecSet<DistSVec<bcomp,dim> > &prevDerivW)  {
                                                               
        
  Timer *modalTimer = domain.getTimer();
  double t0;
                                                       
  // loop over modes
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp kImag(0.0, kFreq);
                                                               
  int nRhsCalc = 0;
  int nRecycleVect = 10;
                                                               
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  //DistSVec<double, dim> delWReal(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> prevWtemp(domain.getNodeDistInfo());
  Vec3D x0(0.0, 0.0, 0.0);
                                                               
                                                               
  Vec<double> modalF(nStrMode);
  rhs = oneReal*DX[0] + kImag*DE[0];

  t0 = modalTimer->getTime();
  kspCompGcr->setup(1, 40, rhs);
  kspCompGcr->printParam();
  VecSet<Vec<bcomp> > aeroOp(nStrMode, nStrMode);
                                                               
                                                               
  kspCompGcr->numCalcVec = 0;
  int iMode;
  int iPadeDeriv;
  for (iMode = 0; iMode < nStrMode; iMode++)  {
    // form rhs
    kspCompGcr->numCalcVec = 0;
    rhs = -oneReal*DX[iMode] - kImag*DE[iMode];
                                                               
                                                               
    // solve
    com->fprintf(stderr, " ... Solving for mode %d, w/rhs norm %e\n", iMode+1, rhs.norm());
    
    delW = prevDerivW[iMode*(nPadeDeriv +1)];
    kspCompGcr->solveMRhs(rhs, delW);
    nRhsCalc += 1;
    if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
      snapsCoarse[iSnap++].getReal(delW);
                                                               
                                                               
      if (kFreq != 0.0)
        snapsCoarse[iSnap++].getImag(delW);
    }
                                                               
                                                               
    prevDerivW[iMode*(nPadeDeriv +1)] = delW;
                                                               
    if (nPadeDeriv > 0 ) {
                                                               
      controlVolComp = oneImag*controlVol;
      prevWtemp = prevDerivW[iMode*(nPadeDeriv +1)];
      prevWtemp *= controlVolComp;
                                                               
      rhs = -oneImag*DE[iMode] - prevWtemp;
                                                               
      // solve for the derivative
      com->fprintf(stderr, " ... Solving for the first derivative for mode %d\n", iMode+1);
      delW = prevDerivW[iMode*(nPadeDeriv +1)+1];
      kspCompGcr->solveMRhs(rhs, delW);
      nRhsCalc += 1;
      if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
        snapsCoarse[iSnap++].getReal(delW);
                                                               
                                                               
        snapsCoarse[iSnap++].getImag(delW);
      }
                                                               
      prevDerivW[iMode*(nPadeDeriv +1)+1] = delW;
                                                               
      if (nPadeDeriv > 1) {
        for (iPadeDeriv = 2; iPadeDeriv < nPadeDeriv+1; iPadeDeriv++) {
                                                               
          bcomp iPadeDerivComp(0.0,iPadeDeriv);
                                                               
          prevWtemp = prevDerivW[iMode*(nPadeDeriv +1)+(iPadeDeriv-1)];
          prevWtemp *= controlVolComp;
                                                               
                                                               
          rhs = -iPadeDerivComp*prevWtemp;
                                                               
          // solve for the derivative
          if (iPadeDeriv == 2)
            com->fprintf(stderr, " ... Solving for the second derivative for mode %d\n", iMode+1);
          else if (iPadeDeriv == 3)
            com->fprintf(stderr, " ... Solving for the third derivative for mode %d\n", iMode+1);
          else
            com->fprintf(stderr, " ... Solving for the %dth derivative for mode %d\n", iPadeDeriv, iMode+1);
          delW = prevDerivW[iMode*(nPadeDeriv +1)+iPadeDeriv];
          kspCompGcr->solveMRhs(rhs, delW);
          nRhsCalc += 1;
                                                               
          if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
            snapsCoarse[iSnap++].getReal(delW);
                                                                                                                               
            snapsCoarse[iSnap++].getImag(delW);
          }
                                                               
          prevDerivW[iMode*(nPadeDeriv +1)+iPadeDeriv] = delW;
        }
      }
    }
                                                               
  }
  modalTimer->addSnapsLinSolvTime(t0);     
                                                         
}
                                                               
                                                               
//-----------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::makeFreqPOD(VecSet<DistSVec<double, dim> > &snaps, int nSnaps){
  
  Timer *modalTimer = domain.getTimer();

  int nPOD = ioData->linearizedData.numPOD;
  int nconv = 0;

  com->fprintf(stderr, " ... Forming %d x %d Correlation Matrix\n", nSnaps, nSnaps);

  VecSet<DistSVec<double, dim> > Utrue(nSnaps, domain.getNodeDistInfo());
  Vec<double> singVals(nSnaps);
  FullM Vtrue(nSnaps);

  double t0 = modalTimer->getTime();
  parallelSVD(snaps, Utrue, singVals.data(), Vtrue, nSnaps);
  modalTimer->addEigSolvTime(t0);
#ifdef DO_MODAL
#ifdef DO_MODAL
  outputPODVectors(Utrue, singVals, nPOD);
#endif

/*#ifdef DO_MODAL
  // allocate for upper half of sym. eigprob
  double *eigVals = new double[nPOD];
  double *eigVecs = new double[nPOD*nSnaps];
  double *rVals = new double[nSnaps*(nSnaps+1)/2];
  for (int i = 0; i < nSnaps; i++)
    for (int j = 0; j <= i; j++)
      rVals[(i+1)*i/2 + j] = snaps[j] * snaps[i];

  double tolerance = ioData->linearizedData.tolerance;

  ARdsSymMatrix<double> pod(nSnaps, rVals, 'U');
  com->fprintf(stderr, " ... Factoring Correlation Matrix\n");
  pod.FactorA();
  ARluSymStdEig<double> podEigProb(nPOD, pod, "LM", nSnaps-1, tolerance, 300*nPOD);
  modalTimer->addCorrelMatrixTime(t0);

  t0 = modalTimer->getTime();
  nconv = podEigProb.FindEigenvectors();
  modalTimer->addEigSolvTime(t0);

  com->fprintf(stderr, " ... Got %d converged eigenvectors out of %d snaps\n", nconv, nSnaps);
  outputPODVectors(podEigProb, snaps, nPOD, nSnaps);
*/
#else
  com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH SCALAPACK \n");
  exit(-1);

#endif
}

//---------------------------------------------------------------------------------------
template<int dim>
void
ModalSolver<dim>::outputModalDisp(double *delU, double *delY, double sdt, int cnt,
                                  int nStrMode, FILE *fp) {

  for (int j = 0; j < nStrMode; j++)
    com->fprintf(fp, "%e  ", delU[j]);
  com->fprintf(fp, "\n ");

}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::setTransfer(int* locDomSize,int* cpuNodes,int &nCpuBl,int blockFactor) {

  int numLocSub = domain.getNumLocSub();

  // loop over all domain nodes and compute the cpu to send to
  // -1 indicates that it's a slave node and no sending
  DistVec<int> cpuDestination(domain.getNodeDistInfo());
  cpuDestination = -1;

  SubDomain **subDomain = domain.getSubDomain();

  int nTotCpus = com->size();
  int thisCPU = com->cpuNum();

  // count non-overlapping subdomain sizes
  for (int iCpu = 0; iCpu < nTotCpus; iCpu++)
    locDomSize[iCpu] = 0;

  #pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    bool *locMasterFlag = cpuDestination.getMasterFlag(iSub);
    for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); iNode++)  {
      if (locMasterFlag[iNode])
        locDomSize[thisCPU]++;
    }
  }

  int maxDomainSize = locDomSize[thisCPU];
  com->globalMax(1, &maxDomainSize);
  //add globally
  com->globalSum(nTotCpus, locDomSize);

  // compute the number of blocks needed bu cpu
  int numTotalNodes = 0;
  for (int iCpu=0; iCpu < nTotCpus; ++iCpu)
    numTotalNodes += locDomSize[iCpu];
  //com->fprintf(stderr, "There are %d nodes in total \n",numTotalNodes);

  //total number of blocks required
  int numBlocks = int(ceil(double(numTotalNodes)/double(blockFactor)));
  //size of the last block
  int nLastBlock = blockFactor - (numBlocks*blockFactor - numTotalNodes);
  //number of blocks per Cpu
  nCpuBl = int(ceil(double(numBlocks)/double(nTotCpus)));
  int *cpuBlocks = new int[nTotCpus];
  for (int iCpu=0; iCpu < nTotCpus; ++iCpu)
    cpuBlocks[iCpu] = nCpuBl;
  for (int i=1; i<=nTotCpus*nCpuBl-numBlocks; ++i)
    cpuBlocks[nTotCpus-i] -= 1;
  int lastCpu;
  if (nLastBlock>0){
    lastCpu = nTotCpus - (nTotCpus*nCpuBl-numBlocks) -1;
  }
  else
    lastCpu = nTotCpus - (nTotCpus*nCpuBl-numBlocks);

  // compute the number of nodes needed by cpu
  for (int iCpu = 0; iCpu < nTotCpus; ++iCpu)
    cpuNodes[iCpu] = cpuBlocks[iCpu]*blockFactor;
  if (nLastBlock < blockFactor)
    cpuNodes[lastCpu] += (nLastBlock-blockFactor);

  delete[] cpuBlocks;
}
//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::parallelSVD(VecSet< DistSVec<double, dim> > &snaps, VecSet<DistSVec<double, dim> > &Utrue, double *S, FullM &Vtrue, int nSnaps) {

#ifdef DO_SCALAPACK
  int numLocSub = domain.getNumLocSub();

  // specify the block size (in terms of nodes)
  // blockfactor > nSnaps => the right singular vector matrix V is computed by cpu 0
  int blockFactor = int(floor(2*nSnaps/dim));

  SubDomain **subDomain = domain.getSubDomain();

  int nTotCpus = com->size();
  int thisCPU = com->cpuNum();

  int *cpuNodes = new int[nTotCpus];
  int *locDomSize = new int[nTotCpus];
  int *locSendReceive = new int[nTotCpus];
  //set the transfer parameters
  int nCpuBl=0;
  setTransfer(locDomSize,cpuNodes,nCpuBl,blockFactor);
  //create the matrix of snapshots with the right distribution
  double *subMat = new double[blockFactor*dim*nCpuBl*nSnaps];
  for (int i=0; i < blockFactor*dim*nCpuBl*nSnaps; ++i)
    subMat[i] =0.0;


  //transfer the extra nodes where needed and fill subMat
  transferData(snaps,subMat,locDomSize,cpuNodes,locSendReceive,nSnaps,nCpuBl,blockFactor);

  com->barrier();
  // Allocate for svd call
  // allocate for eigenvectors and eigenvalues
  int locLLD = dim * cpuNodes[thisCPU];
  int maxLLD = locLLD;
  com->globalMax(1, &maxLLD);
  int maxLocC = nSnaps;
  int locLLD_V = 1;
  if (thisCPU == 0)
    locLLD_V = nSnaps;

  double *U = new double[locLLD*maxLocC];
  double *V = new double[locLLD_V*maxLocC];
  for (int iSnaps=0; iSnaps<nSnaps; ++iSnaps)
    S[iSnaps] = 0.0;
  int nprow = nTotCpus;
  int npcol = 1;
  int globNumRows = dim * snaps[0].nonOverlapSize();
  // call svd
  int lwork = 0;
  int myrow, mycol, ictxt, info;
  int rowIndex = 1;
  int colIndex = 1;
  blockFactor *= dim;


  F77NAME(lworksize)(ictxt, nprow, npcol, myrow, mycol, globNumRows, nSnaps,
                   blockFactor, blockFactor, lwork);

  com->barrier();
  double *work = new double[lwork + 2];
  work[0] = -99.0;
  
  com->barrier();
  F77NAME(thinsvd)(ictxt, nprow, npcol, myrow, mycol, globNumRows, nSnaps,
                   blockFactor, blockFactor, subMat, maxLLD, locLLD, nSnaps,
                   locLLD_V, thisCPU, rowIndex, colIndex, U, S, V, lwork, work,
                   info);

  com->barrier();

  for (int i = 0; i < nSnaps; i++){
    for (int j = 0; j < nSnaps; j++) {
      Vtrue[i][j] = 0.0;
      if (thisCPU==0)
        Vtrue[i][j] = V[nSnaps*i+j];
    }
  }
  delete[] V;

  com->globalSum(nSnaps*nSnaps, Vtrue.data());

  //Permute back matrix U
  com->barrier();

  for (int i = 0; i < nSnaps; ++i)
    Utrue[i] = 0.0;
  transferDataBack(U, Utrue , locDomSize, cpuNodes, locSendReceive, nSnaps);

  com->barrier();
  delete[] locDomSize;
  delete[] cpuNodes;
  delete[] locSendReceive;
  delete[] work;
  delete[] U;

#else
  com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH SCALAPACK and DO_SCALAPACK Flag\n");
  exit(-1);
#endif
}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::transferData(VecSet< DistSVec<double, dim> > &snaps, double* subMat, int *locDomSize, int *cpuNodes, int *locSendReceive, int nSnaps, int nCpuBl, int blockFactor ) {


  int numLocSub = domain.getNumLocSub();
  DistInfo &nodeDistInfo = domain.getNodeDistInfo();
  int nTotCpus = com->size();
  int thisCPU = com->cpuNum();

  int *locDomSizeCpy = new int[nTotCpus];
  for (int iCpu=0; iCpu<nTotCpus; ++iCpu)
    locDomSizeCpy[iCpu] = locDomSize[iCpu];

  int subMatRowsNum = cpuNodes[thisCPU]*dim;
  double *recData;
  double *buffer;

  // send/receive for each cpu (if positive, it has received, if negative, it has sent)
  for (int iCpu = 0; iCpu < nTotCpus; iCpu++)
    locSendReceive[iCpu] = 0;

  // loop over all domain nodes and compute the cpu to send to
  // -1 indicates that it's a slave node and no sending
  DistVec<int> cpuDestination(domain.getNodeDistInfo());
  cpuDestination = -1;

  SubDomain **subDomain = domain.getSubDomain();
  int *totalSentNodes = new int[nTotCpus];
  int nSendNodes,buffLen,index;
  int nRecNodes = 0;
  for (int iCpu = 0; iCpu < nTotCpus; ++iCpu){
    totalSentNodes[iCpu] = 0;
    if (locDomSizeCpy[iCpu] - cpuNodes[iCpu] > 0){
      for (int jCpu = 0; jCpu < nTotCpus; ++jCpu){
        if (locDomSizeCpy[iCpu] - cpuNodes[iCpu] > 0 && locDomSizeCpy[jCpu] - cpuNodes[jCpu] < 0) {
          nSendNodes = min(locDomSizeCpy[iCpu] - cpuNodes[iCpu],-(locDomSizeCpy[jCpu] - cpuNodes[jCpu]));
          //send Here from iCpu to jCpu
          if (thisCPU==iCpu)
            locSendReceive[jCpu] = -nSendNodes;
          else if(thisCPU==jCpu)
            locSendReceive[iCpu] = nSendNodes;

          buffLen = nSendNodes*dim*nSnaps;
          // create array for received data
          if (thisCPU==jCpu)
            recData = new double[buffLen];

          if (thisCPU==iCpu) {
            // create buffer
            buffer = new double[buffLen];
            //fill buffer
            for (int iSnap = 0; iSnap < nSnaps; ++iSnap)  {
              index = 0;
              for (int iSub = 0; iSub < numLocSub; ++iSub) {
                // create buffer
                double (*locSnap)[dim] = snaps[iSnap].subData(iSub);
                bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
                for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
                  if (locMasterFlag[iNode]) {
                    if (index == nSendNodes+totalSentNodes[iCpu])  break;
                    index++;
                    if (index-1 >= totalSentNodes[iCpu])  {
                      for (int j = 0; j < dim; ++j)
                        buffer[(index-totalSentNodes[iCpu]-1)*dim+j+iSnap*nSendNodes*dim] = locSnap[iNode][j];
                    }
                  }
                }
              }
            }
            totalSentNodes[iCpu] = index;

            // send data in buffer
            //fprintf(stderr, "*** CPU #%d sending to CPU #%d: %d nodes = %d entries\n", thisCPU, jCpu, nSendNodes, buffLen);
            com->sendTo(jCpu, thisCPU*nTotCpus+jCpu, buffer, buffLen);
            com->waitForAllReq();
          }

          // receive data and populate submatrix
          if (thisCPU==jCpu) {
            com->recFrom(iCpu, iCpu*nTotCpus+thisCPU, recData, buffLen);
            for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
              for (int k = nRecNodes; k < nSendNodes+nRecNodes; ++k)  {
                for (int j = 0; j < dim; ++j)
                  subMat[iSnap*subMatRowsNum + k*dim+j] = recData[(k-nRecNodes)*dim+j+iSnap*nSendNodes*dim];
              }
            }
            nRecNodes += nSendNodes;
          }
          locDomSizeCpy[jCpu] += nSendNodes;
          locDomSizeCpy[iCpu] -= nSendNodes;
          if (thisCPU==iCpu)
            delete[] buffer;
          else if (thisCPU==jCpu)
            delete[] recData;
        }
      }
      com->barrier();
    }
  }

  //for each CPU fill the rest of subMat with its own data
  for (int iCpu = 0; iCpu < nTotCpus; ++iCpu) {
    if (thisCPU == iCpu) {
      int k = nRecNodes;
      for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
        k = nRecNodes;
        index = 0;
        for (int iSub = 0; iSub < numLocSub; ++iSub) {
          double (*locSnap)[dim] = snaps[iSnap].subData(iSub);
          bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
          for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
            if (locMasterFlag[iNode]) {
              //if not already sent
              if (index >= totalSentNodes[iCpu])  {
                for (int j = 0; j < dim; ++j)
                  subMat[iSnap*subMatRowsNum + k*dim+j] = locSnap[iNode][j];
                if (iSnap==nSnaps-1)
                  locDomSizeCpy[iCpu]--;
                k++;
              }
              index++;
            }
          }
        }
      }
    }
  }
  delete[] totalSentNodes;
}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::transferDataBack(double *U, VecSet< DistSVec<double, dim> > &Utrue , int *locDomSize, int *cpuNodes, int *locSendReceive, int nSnaps) {

  int numLocSub = domain.getNumLocSub();
  DistInfo &nodeDistInfo = domain.getNodeDistInfo();
  int nTotCpus = com->size();
  int thisCPU = com->cpuNum();

  double *recData;
  double *buffer;

  int locLLD = dim * cpuNodes[thisCPU];
  // loop over all domain nodes and compute the cpu to send to
  // -1 indicates that it's a slave node and no sending
  DistVec<int> cpuDestination(domain.getNodeDistInfo());
  cpuDestination = -1;

  SubDomain **subDomain = domain.getSubDomain();
  int *totalSentNodes = new int[nTotCpus];
  int nSendNodes,buffLen,index;
  int transf = 0;
  int nRecNodes = 0;
  com->barrier();
  for (int iCpu=0; iCpu < nTotCpus; ++iCpu) {
    totalSentNodes[iCpu] = 0;
    for (int jCpu=0; jCpu < nTotCpus; ++jCpu) {
      nSendNodes = 0;
      buffLen = 0;
      transf = 0;
      if (thisCPU==iCpu) {
        if (locSendReceive[jCpu] > 0) {
          nSendNodes = locSendReceive[jCpu];
          buffLen = nSendNodes*dim*nSnaps;
          transf = 1;
        }
      }
      else if (thisCPU==jCpu) {
        if (locSendReceive[iCpu] < 0) {
          nSendNodes = -locSendReceive[iCpu];
          buffLen = nSendNodes*dim*nSnaps;
          transf = 1;
        }
      }
      com->barrier();
      com->globalSum(1, &transf);
      if (transf > 1) {
        // create array for received data
        if (thisCPU==jCpu)
          recData = new double[buffLen];
        com->barrier();
        //create buffer
        if (thisCPU==iCpu) {
          buffer = new double[buffLen];
           //fill buffer
          for (int iSnap = 0; iSnap < nSnaps; ++iSnap)  {
            index = 0;
            for (int iSub = 0; iSub < numLocSub; ++iSub) {
              // create buffer
              bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
              for (int iNode = 0; iNode < cpuNodes[iCpu]; ++iNode)  {
                if (index == nSendNodes+totalSentNodes[iCpu])  break;
                index++;
                if (index-1 >= totalSentNodes[iCpu])  {
                  for (int j = 0; j < dim; ++j)
                    buffer[(index-totalSentNodes[iCpu]-1)*dim+j+iSnap*nSendNodes*dim] = U[locLLD*iSnap+dim*iNode+j];
                }
              }
            }
          }
          totalSentNodes[iCpu] = index;
          // send data in buffer
          //fprintf(stderr, "*** CPU #%d sending back to CPU #%d: %d nodes = %d entries\n", thisCPU, jCpu, nSendNodes, buffLen);
          com->sendTo(jCpu, 10*thisCPU*nTotCpus+jCpu, buffer, buffLen);
        }
        com->barrier();
        // receive data and populate submatrix
        if (thisCPU==jCpu) {
          com->recFrom(iCpu, 10*iCpu*nTotCpus+thisCPU, recData, buffLen);
          for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
            index = 0;
            int k=0;
            for (int iSub = 0; iSub < numLocSub; ++iSub) {
              double (*locU)[dim] = Utrue[iSnap].subData(iSub);
              bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
              for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
                if (locMasterFlag[iNode]) {
                  if (k>=nRecNodes && k <nSendNodes+nRecNodes) {
                    for (int j = 0; j < dim; ++j)
                      locU[iNode][j] = recData[(k-nRecNodes)*dim+j+iSnap*nSendNodes*dim];
                  }
                  k++;
                }
              }
            }
          }
          nRecNodes += nSendNodes;
          delete[] recData;

        }
        com->barrier();
        if (thisCPU==iCpu)
          delete[] buffer;
      }
      com->barrier();
    }
  }
  //Completing the rest of the matrices
  for (int iCpu = 0; iCpu < nTotCpus; ++iCpu) {
    if (thisCPU == iCpu) {
      int k;
      for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
        k = totalSentNodes[iCpu];
        index = 0;
        for (int iSub = 0; iSub < numLocSub; ++iSub) {
          double (*locU)[dim] = Utrue[iSnap].subData(iSub);
          bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
          for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
            if (locMasterFlag[iNode]) {
              //if not already sent
              if (index >= nRecNodes)  {
                for (int j = 0; j < dim; ++j)
                  locU[iNode][j] = U[iSnap*locLLD+k*dim+j];
                k++;
              }
              index++;
            }
          }
        }
      }
    }
  }

  //complete the slave nodes
  CommPattern<double> *vpat = domain.getVecPat();
  for (int iSnap=0; iSnap < nSnaps; ++iSnap)
    domain.assemble(vpat, Utrue[iSnap]);
}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::interpolatePOD()  {


  Timer *modalTimer = domain.getTimer();
  com->fprintf(stderr, " ... Interpolating POD on a tangent space to the Grassmann manifold \n");

  char *vecFile = tInput->podFile;
  if (!vecFile)
    vecFile = "podFiles.in";
  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
    exit (-1);
  }
  int nData;
  fscanf(inFP, "%d",&nData);

  char **podFile = new char *[nData];

  for (int iData = 0; iData < nData; ++iData){
    podFile[iData] = new char[500];
    //char *podFile1 = new char[500];
    fscanf(inFP, "%s", podFile[iData]);
    //podFile[iData] = podFile1;
    com->fprintf(stderr, " ... Reading POD from %s \n", podFile[iData]);
  }


  if (ioData->output.transient.podFile[0] == 0)  {
    com->fprintf(stderr, "*** ERROR: POD Basis File not specified\n");
    exit (-1);
  }
  int sp = strlen(ioData->output.transient.prefix);
  char *outFile = new char[sp + strlen(ioData->output.transient.podFile)+1];
  sprintf(outFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);

  int numPod = ioData->linearizedData.numPOD;

  double *mach = new double[nData];
  double *angle =new double[nData];
  double newMach;
  double newAngle;
  for (int iData = 0; iData < nData; ++iData)
    fscanf(inFP, "%lf", mach+iData);
  fscanf(inFP, "%lf", &newMach);
  for (int iData = 0; iData < nData; ++iData)
    fscanf(inFP, "%lf", angle+iData);
  fscanf(inFP, "%lf", &newAngle);

  com->fprintf(stderr, " ... Interpolating new POD basis at Mach = %f and angle of attack = %f from :\n",newMach, newAngle);
  for (int iData = 0; iData < nData; ++iData)
    com->fprintf(stderr,"      -> Mach = %f and angle of attack = %f\n", mach[iData], angle[iData]);


  VecSet< DistSVec<double, dim> > **pod = new VecSet< DistSVec<double, dim> >*[nData];
  for (int iData=0; iData < nData; ++iData){
    pod[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
    // read in Pod Vectors
    double *eig = new double[numPod];

    domain.readVectorFromFile(podFile[iData], 0, &eig[0], (*pod[iData])[0] );
    if (numPod > eig[0])  {
      com->fprintf(stderr, "*** Warning: Resetting number of interpolated POD vectors from %d to %d\n", numPod, (int) eig[0]);
      numPod = (int) eig[0];
    }

    for (int iPod = 0; iPod < numPod; ++iPod)
      domain.readVectorFromFile(podFile[iData], iPod+1, &eig[iPod], (*pod[iData])[iPod]);
      delete [] eig;
  }

  double t0 = modalTimer->getTime();

  //build the reduced coordinates
  double *reducedMach = new double[nData];
  double *reducedAngle = new double[nData];
  if (newAngle!=0.0) {
    for (int iData=0; iData < nData; ++iData) {
      reducedMach[iData] = mach[iData]/newMach-1;
      reducedAngle[iData] = angle[iData]/newAngle-1;
    }
  }
  else {
    for (int iData=0; iData < nData; ++iData) {
      reducedMach[iData] = mach[iData]-newMach;
      reducedAngle[iData] = angle[iData]-newAngle;
    }
  }
  //find the central point
  int iDataMin = 0;
  double radMinSq = reducedMach[0]*reducedMach[0] + reducedAngle[0]*reducedAngle[0];
  double radSq = 0.0;
  for (int iData=1; iData < nData; ++iData) {
     radSq = reducedMach[iData]*reducedMach[iData] + reducedAngle[iData]*reducedAngle[iData];
     if (radSq < radMinSq){
       iDataMin = iData;
       radMinSq = radSq;
     }
   }
  //store the reference pod basis
  VecSet< DistSVec<double, dim> > podRef(*(pod[iDataMin]));

  //Logarithmic mappings
  VecSet< DistSVec<double, dim> > **projMap = new VecSet< DistSVec<double, dim> >*[nData];
  FullM matVals(numPod);

  //compute the matrices (Phi0'*Phi)^(-1)
  for (int iData=0; iData < nData; ++iData){
    if (iData!=iDataMin) {
      for (int j = 0; j < numPod; ++j) {
        for (int k = 0; k < numPod; ++k)
          matVals[j][ k] = podRef[j] * ((*pod[iData])[k]);
      }
      com->barrier();
      matVals.invert();

      //compute the projection mapping = Phi*(Phi0'*Phi)^(-1)-Phi0
      projMap[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
      com->barrier();
      for (int iPod = 0; iPod < numPod; ++iPod) {
        (*projMap[iData])[iPod] = 0.0;
        for (int jPod = 0; jPod < numPod; ++jPod) {
          (*projMap[iData])[iPod] += (*pod[iData])[jPod] * matVals[jPod][iPod];
        }
      }
      com->barrier();
      for (int iPod = 0; iPod < numPod; ++iPod)
        (*projMap[iData])[iPod] = (*projMap[iData])[iPod] - podRef[iPod];
      com->barrier();
   }
  }

  com->barrier();
  for (int iData=0; iData<nData;++iData){
    delete pod[iData];
    com->barrier();
  }

  com->barrier();
  delete [] pod;
  //SVD decomposition of projMap
  VecSet<DistSVec<double, dim> > U(numPod, domain.getNodeDistInfo());
  double *Sigma = new double[numPod];
  FullM V(numPod);
  double *Theta = new double[numPod];
  VecSet<DistSVec<double, dim> > U2(numPod, domain.getNodeDistInfo());
  double *Sigma2 = new double[numPod];
  FullM V2(numPod);
  //Logarithmic mapping
  VecSet< DistSVec<double, dim> > **logMap = new VecSet< DistSVec<double, dim> >*[nData];
  for (int iData=0; iData < nData; ++iData) {
    if (iData!=iDataMin) {
      //compute SVD
      parallelSVD(*projMap[iData],U,Sigma,V,numPod);//call SVD here
      com->barrier();
      delete projMap[iData];

      //compute tan^(-1) Sigma
      for (int iPod = 0; iPod < numPod; ++iPod)
        Theta[iPod] = atan(Sigma[iPod]);
      //compute the logarithmic mapping

      //build logMap
      logMap[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
      for (int iPod = 0; iPod < numPod; ++iPod) {
        (*logMap[iData])[iPod] = 0.0;
        for (int jPod = 0; jPod < numPod; ++jPod){
          (*logMap[iData])[iPod] += ((Theta[jPod]*V[iPod][jPod])*U[jPod]);  // logMap = U tan^(-1) Sig V^T
        }
      }
   }
  }
  delete [] Theta;
  delete [] Sigma;

  VecSet< DistSVec<double, dim> > logMapInterp(numPod, domain.getNodeDistInfo());
  // Interpolation (Hardy's multiquadrics method is here used)
  // see "Two Dimensional Spline Interpolation Algorithms", H. Spath, A.K. Peters Ltd  
  double q = 0.25;
  double Rsq;

  //find Rsq
  double minReducedM = reducedMach[0];
  double maxReducedM = reducedMach[0];
  double minReducedA = reducedAngle[0];
  double maxReducedA = reducedAngle[0];
  for (int iData=1; iData < nData; ++iData) {
    if (reducedMach[iData] < minReducedM)
      minReducedM = reducedMach[iData];
    if (reducedMach[iData] > maxReducedM)
      maxReducedM = reducedMach[iData];
    if (reducedAngle[iData] < minReducedA)
      minReducedA = reducedAngle[iData];
    if (reducedAngle[iData] > maxReducedA)
      maxReducedA = reducedAngle[iData];
  }
  double maxDelM = maxReducedM - minReducedM;
  double maxDelA = maxReducedA - minReducedA;
  double maxConst = maxDelM > maxDelA ? maxDelM : maxDelA;
  Rsq = 0.1*maxConst;

  FullM B(nData), b(nData,1);


  double rsq;
  for (int j=0; j < nData; ++j) {
    for (int k=0; k < nData; ++k) {
      rsq = pow(reducedMach[j]-reducedMach[k],2)+pow(reducedAngle[j]-reducedAngle[k],2);
      B[j][k] = pow(rsq+Rsq,q);
    }
    rsq = pow(reducedMach[j],2)+pow(reducedAngle[j],2);
    b[j][0] = pow(rsq+Rsq,q);
  }
  B.invert();

  //compute the interpolated POD basis
  domain.hardyInterpolationLogMap(logMap,logMapInterp,nData,numPod,iDataMin,B,b);
  //Exponential mapping
  VecSet<DistSVec<double, dim> > UInterp(numPod, domain.getNodeDistInfo());
  double *SigInterp = new double[numPod];
  FullM VInterp(numPod);

  parallelSVD(logMapInterp,UInterp,SigInterp,VInterp,numPod);//call SVD here

  //compute podRefVInterp
  VecSet<DistSVec<double, dim> > podRefVInterp(numPod, domain.getNodeDistInfo());
  for (int iPod = 0; iPod < numPod; ++iPod) {
    podRefVInterp[iPod] = 0.0;
    for (int jPod = 0; jPod < numPod; ++jPod)
      podRefVInterp[iPod] += VInterp[jPod][iPod]*podRef[jPod];
 }

  VecSet<DistSVec<double, dim> >newPOD(numPod,domain.getNodeDistInfo());

  for (int iPod = 0; iPod < numPod; ++iPod) {
    newPOD[iPod] = cos(SigInterp[iPod])*podRefVInterp[iPod]+ sin(SigInterp[iPod])*UInterp[iPod];

  //Do a Gramm-Schmidt to finalize the POD reconstruction
  for (int j = 0; j < iPod; ++j)
    newPOD[iPod] -= (newPOD[iPod] * newPOD[j])*newPOD[j];
    newPOD[iPod] *= 1.0/newPOD[iPod].norm();

  }

  com->fprintf(stderr, " ... Writing Interpolated POD to %s\n", outFile);
  // output interpolated pod vectors
  domain.writeVectorToFile(outFile, 0, numPod, newPOD[0]);
  double newEig=-1.0;

  for (int jPod = 0; jPod < numPod; ++jPod)  {
    com->fprintf(stderr, " ... Writing pod vec %d\n", jPod+1);
    domain.writeVectorToFile(outFile, jPod+1, newEig, newPOD[jPod]);
  }
  modalTimer->addMeshSolutionTime(t0);

  delete[] SigInterp;
  for (int iData=0; iData< nData; ++iData)
    delete [] podFile[iData];
  delete[] podFile;
  delete[] outFile;

  modalTimer->setRunTime();


}
//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::evalFluidSys(VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcn(*ioData); 
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  onlyHOp->evalH(0, Xref, controlVol, Uref);

  VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);

  // form H ROM op
  DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());

  double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;

  double timeDimConst =  ioData->ref.mach * sqrt(gamma) / ioData->ref.length;
  com->fprintf(stderr, " ... Using M = %f, Gamma = %f, rho = %e, L = %f for dim constant: %16.10e\n",  
       ioData->ref.mach, gamma, ioData->ref.density, ioData->ref.length, timeDimConst);
 
  int iVec, jVec;

  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    onlyHOp->apply(podVecs[iVec], tmpVec);
    tmpVec /= controlVol;
    for (jVec = 0; jVec < nPodVecs; jVec++)
      romOperator[iVec][jVec] = (podVecs[jVec] * tmpVec) * timeDimConst;
  }

  com->fprintf(stderr, " ... Computed Rom Operator\n");

  com->fprintf(stderr, " ... Created Fluid System\n");
  int sp = strlen(ioData->output.transient.prefix);
  char *romFile = new char[sp + strlen(ioData->output.transient.romFile)+1];
  sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
  FILE *romFP = fopen(romFile, "w");
  com->barrier();

  com->fprintf(romFP, "%d 0\n", nPodVecs);
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    for (jVec = 0; jVec < nPodVecs; jVec++)
      com->fprintf(romFP, "%.16e ", romOperator[jVec][iVec]);
    com->fprintf(romFP, "\n");
  }

  delete[] romFile;

  //spaceOp has to be reset because it has been modified by the apply function
  DistSVec<double,dim> FF(domain.getNodeDistInfo());
  spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);
  
}

//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::evalAeroSys(VecSet<Vec<double> > &outRom, 
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcn(*ioData); 
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  onlyHOp->evalH(0, Xref, controlVol, Uref);

  VecSet<Vec<double> > ecVecs(nStrMode, nPodVecs);
  VecSet<Vec<double> > gVecs(nStrMode, nPodVecs);
  VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);

  // form H ROM op
  DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
  DistSVec<double,dim> tmpVec2(domain.getNodeDistInfo());

  double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;

  double timeDimConst = -ioData->ref.mach * sqrt(gamma) / ioData->ref.length;
  com->fprintf(stderr, " ... Using M = %f, Gamma = %f, rho = %e, L = %f for dim constant: %16.10e\n",  ioData->ref.mach, gamma, ioData->ref.density, ioData->ref.length, timeDimConst);
  
  int iVec, jVec;

  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    onlyHOp->apply(podVecs[iVec], tmpVec);
    tmpVec /= controlVol;
    for (jVec = 0; jVec < nPodVecs; jVec++)
      romOperator[iVec][jVec] = (podVecs[jVec] * tmpVec) * timeDimConst;
  }

  com->fprintf(stderr, " ... Computed Rom Operator\n");

  // form coupling matrix ROMs
  // form Structural Component: K
  double *structSys = new double[nStrMode];

  double invDt = -1.0/ioData->ref.rv.time;

  for (iVec = 0; iVec < nStrMode; iVec++)  {
    //tmpVec = DE[iVec];
    tmpVec = DE[iVec] * invDt; 
    tmpVec2 = DX[iVec] ; 
    tmpVec /= controlVol;
    tmpVec2 /= controlVol;
 
    for (jVec = 0; jVec < nPodVecs; jVec++)  {
      ecVecs[iVec][jVec] = podVecs[jVec] * tmpVec;
      gVecs[iVec][jVec] = podVecs[jVec] * tmpVec2 * timeDimConst;
    }

    structSys[iVec] = -K[iVec];
  }

  com->fprintf(stderr, " ... Computed Coupling Input Rom Vectors\n");

  int sysSize = 2*nStrMode + nPodVecs;
  double *sysVals = new double[sysSize*sysSize];
  for (int j = 0; j < sysSize*sysSize; j++)
    sysVals[j] = 0.0;
  com->fprintf(stderr, " ... Forming Aeroelastic System of size: %d\n", sysSize);

  // populate 1,1 & 2,1 block
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    for (jVec = 0; jVec < nPodVecs; jVec++)
      sysVals[iVec*sysSize+jVec] = romOperator[iVec][jVec];
    for (jVec = 0; jVec < nStrMode; jVec++)
      sysVals[iVec*sysSize + nPodVecs + jVec] = outRom[iVec][jVec];  
  }

  // populate 1,2 and 3,2 block

  for (iVec = 0; iVec < nStrMode; iVec++)  {
    for (jVec = 0; jVec < nPodVecs; jVec++)
      sysVals[sysSize*nPodVecs + iVec*sysSize + jVec] = ecVecs[iVec][jVec];
    for (jVec = 0; jVec < nStrMode; jVec++) {
      if (jVec == iVec)
        sysVals[sysSize*nPodVecs + iVec*sysSize + nPodVecs+nStrMode + jVec] = 1.0;
      else
        sysVals[sysSize*nPodVecs + iVec*sysSize + nPodVecs+nStrMode + jVec] = 0.0;
    }
  }

  // populate 1,3 and 2,3 block
  for (iVec = 0; iVec < nStrMode; iVec++)  {
    for (jVec = 0; jVec < nPodVecs; jVec++)
      sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + jVec] = gVecs[iVec][jVec];
    for (jVec = 0; jVec < nStrMode; jVec++) {
      if (jVec == iVec)
        sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + nPodVecs + jVec] = structSys[iVec];
      else
        sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + nPodVecs + jVec] = 0.0;
    }
  }

  com->fprintf(stderr, " ... Created Aeroelastic System\n");
  int sp = strlen(ioData->output.transient.prefix);
  char *romFile = new char[sp + strlen(ioData->output.transient.romFile)+1];
  sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
  FILE *romFP = fopen(romFile, "w");
  com->barrier();

  com->fprintf(romFP, "%d %d\n", nPodVecs, nStrMode);
  for (iVec = 0; iVec < sysSize; iVec++)  {
    for (jVec = 0; jVec < sysSize; jVec++){
      com->fprintf(romFP, "%.16e ", sysVals[jVec*sysSize+iVec]);
   }       
   com->fprintf(romFP, "\n");
  }

  delete[] romFile;

  //spaceOp has to be reset because it has been modified by the apply function
  DistSVec<double,dim> FF(domain.getNodeDistInfo());
  spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);

}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::formOutputRom(VecSet<Vec<double> > &outputRom,
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  Vec<double> modalF(nStrMode);
  modalF = 0.0;

  double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;
  double machSquare = ioData->ref.mach * ioData->ref.mach;
  double lengthSquare = ioData->ref.length * ioData->ref.length;

  for (int iVec = 0; iVec < nPodVecs; iVec++)  {
    postOp->computeForceDerivs(Xref, Uref, podVecs[iVec], modalF, mX);
    outputRom[iVec] = machSquare * gamma * lengthSquare * modalF;
  }
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void ModalSolver<dim>::readPodVecs(VecSet<DistSVec<Scalar, dim> > &podVecs,
                       int &nPod)  {

  // read in POD Vectors
  char *vecFile = tInput->podFile;

  double eigValue; 
  int nPodVecs;

  // read number of vecs
  DistSVec<Scalar,dim> tmpVec(domain.getNodeDistInfo());
  domain.readVectorFromFile(vecFile, 0, &eigValue, tmpVec);

  nPodVecs = (int) eigValue;
  com->fprintf(stderr, " ... There are %d total podVecs \n", nPodVecs);

  if (nPod > nPodVecs)  {
    com->fprintf(stderr, " ... There are only %d POD Vectors \n", nPodVecs);
    nPod = nPodVecs;
  }
  else
    nPodVecs = nPod;

  com->fprintf(stderr, " ... Reading %d POD Vectors from file %s\n", nPodVecs, vecFile);

  int iVec;
  double firstEig;
  domain.readVectorFromFile(vecFile, 1, &firstEig, podVecs[0]);
  for (iVec = 1; iVec < nPodVecs; iVec++)
    domain.readVectorFromFile(vecFile, iVec+1, &eigValue, podVecs[iVec]);

  com->fprintf(stderr, " ... Eigenvalue Ratio: (%e/%e) = %e\n", eigValue, firstEig, eigValue/firstEig);

}

//------------------------------------------------------------------------------
#ifdef DO_MODAL

template<int dim>
void ModalSolver<dim>::outputPODVectors(VecSet<DistSVec<double, dim> > &podVecs, 
                       Vec<double> &sVals, int nPOD)  {

   Timer *modalTimer = domain.getTimer();

  // write header
  int sp = strlen(ioData->output.transient.prefix) + 1;

  char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
  sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);

  com->fprintf(stderr, " ... Writing %d (%f)podVecs to File\n", nPOD, (double) nPOD);

  domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVecs[0]);

  double t0;
  int jj, kk;
  DistSVec<double, dim> tmpVec(domain.getNodeDistInfo());

  for (jj = 0; jj < nPOD; jj++)  {
    t0 = modalTimer->getTime();

    // do gram-schmidt on podVecs
    tmpVec = 0.0;
    for (int jVec = 0; jVec < jj; jVec++)  {
      //podVecs[jj] -= podVecs[jVec] * (podVecs[jj]*podVecs[jVec]);
      tmpVec += podVecs[jVec] * (podVecs[jj]*podVecs[jVec]);
    }
    podVecs[jj] -= tmpVec;

    double norm = sqrt(podVecs[jj]*podVecs[jj]);
    podVecs[jj] *= (1.0 / norm);
    modalTimer->addGramSchmidtTime(t0);
    domain.writeVectorToFile(podFileName, jj+1, sVals[jj], podVecs[jj]);
    com->fprintf(stderr, "%d %e, residual = %e\n", jj, sVals[jj], tmpVec.norm());
  }
}

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::outputPODVectors(ARluSymStdEig<double> &podEigProb,
                       VecSet<DistSVec<double, dim> > &snaps, int nPOD, int numSnapsTot)  {
 
   Timer *modalTimer = domain.getTimer();

  // output POD Vectors
  VecSet<DistSVec<double, dim> > podVecs(nPOD, domain.getNodeDistInfo());

  // write header
  int sp = strlen(ioData->output.transient.prefix) + 1;

  char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
  sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);

  podVecs[0] = 0;
  com->fprintf(stderr, " ... Writing %d (%f)podVecs to File\n", nPOD, (double) nPOD);

  domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVecs[0]);

  double t0;
  double *rawEigVec = new double[numSnapsTot];
  int jj, kk;
  for (jj = 0; jj < nPOD; jj++)  {
    t0 = modalTimer->getTime();
    for (kk = 0; kk < numSnapsTot; kk++)
      rawEigVec[kk] = podEigProb.Eigenvector(nPOD-jj-1, kk);

    Vec<double> pVec(numSnapsTot, rawEigVec);
    podVecs[jj] = 0.0;
    for (kk = 0; kk < numSnapsTot; kk++)
      podVecs[jj] += snaps[kk] * pVec[kk];

    double eig = podEigProb.Eigenvalue(nPOD-jj-1);
    podVecs[jj] *= 1.0/sqrt(eig);

    // do gram-schmidt on podVecs
    for (int jVec = 0; jVec < jj; jVec++)
      podVecs[jj] -= podVecs[jVec] * (podVecs[jj]*podVecs[jVec]);

    double norm = sqrt(podVecs[jj]*podVecs[jj]);
    podVecs[jj] *= (1.0 / norm);
    modalTimer->addGramSchmidtTime(t0);
    domain.writeVectorToFile(podFileName, jj+1, eig, podVecs[jj]);
    com->fprintf(stderr, "%d %e\n", jj, eig);

  }
}
#endif


