#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
using std::sort;


#include <Domain.h>
#include <Modal.h>
#include <IoData.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <MatVecProd.h>
#include <Timer.h>
#include <VarFcnDesc.h>
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
  pc = 0;
  ksp = 0;
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
    if (ioData->output.transient.romFile != "")  {

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
      com->fprintf(stderr, " ... Running Unsteady Linearized ROM w/no structural initializaitons\n");
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

    VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);
    double t0 = modalTimer->getTime();
    if (ioData->linearizedData.type == LinearizedData::ROM)
      constructROM(romOperator, ecVecs, gVecs, podVecs, nPodVecs);
    else if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_ ||
             ioData->problem.alltype == ProblemData::_ROM_)
      constructROM2(romOperator, ecVecs, gVecs, podVecs, nPodVecs);
    modalTimer->addRomConstrTime(t0);

    modalTimer->getTime(); 
    timeIntegrateROM(romOperator, ecVecs, gVecs, podVecs, nSteps, nPodVecs, delU, delY, sdt);
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
  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 
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

  // Uncomment for forced oscillations
  double pi = 3.14159265358979;
  double freq = 2.0*pi*ioData->linearizedData.frequency;
  double dMax = ioData->linearizedData.amplification;

  int i;
  deltmp = 0.0;
  double *prevU = new double[nStrMode];
  double *prevY = new double[nStrMode];
  for (i = 0; i < nStrMode; ++i)  {
    deltmp += (delU[i]+sdt/2*delY[i])*mX[i];
    prevU[i] = delU[i];
    prevY[i] = delY[i];
  }
  
  deltmp += Xref;

  // Time Loop
  ksp->printParam();
  
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

  // Compute Reference Modal Force
  DistSVec<double,3> refNodalForce(domain.getNodeDistInfo());
  DistVec<double> Pin(domain.getFaceDistInfo());
  Vec<double> refModalF(nStrMode);
  com->fprintf(stderr, " ... Computing Ref Nodal Force w/pressure: %e\n", ioData->aero.pressure);
  Pin = ioData->aero.pressure;
  postOp->computeNodalForce(Xref, Uref, Pin, refNodalForce);
  
  CommPattern<double> *vpat = domain.getVecPat();
  domain.assemble(vpat, refNodalForce);

  for (i = 0; i < nStrMode; i++)
    refModalF[i] = mX[i]*refNodalForce;

  int cntr = 0;
  int cntp1;
  for (int cnt = 0; cnt < nSteps; ++cnt) {

    cntp1 = cnt+1;

    t0 = modalTimer->getTime();    

    rhs = (delW*4.0) - delWnm1;

    rhs *= controlVol;

    rhs *= 1.0/dt;

    rhsA = 0.0;  rhsB = 0.0; rhsC = 0.0; rhsD = 0.0;
    // We drop the developed term of G*v_n bc it creates instabilities, so this alg. is formerly at least 1st
    // order accurate
    for (i = 0; i < nStrMode; ++i)  {
      rhsA += DX[i]*delU[i];
      //rhsB += DX[i]*delY[i];
      rhsC += DE[i]*delY[i];
      rhsD += DE[i]*prevY[i];
    }
    rhsA *= 2.0;
    //rhsB *= dt;
    rhsC *= 3.0;

    rhs -= rhsA;
    //rhs -= rhsB;
    rhs -= rhsC;
    rhs += rhsD;

    delWnm1 = delW;

    rhs += FF;
    ksp->solve(rhs, delW);
    modalTimer->addKspTime(t0);

    // for forced oscillations
    if (ioData->linearizedData.type == LinearizedData::FORCED)  {
      delU[ioData->linearizedData.modeNumber-1] = sin(freq*(cnt)*sdt)*dMax;
      delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*(cnt)*sdt)*dMax;
    }

    if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {
    
      for (int kk = 0; kk < nStrMode; kk++)  {
        prevU[kk] = delU[kk];
        prevY[kk] = delY[kk];
      }
      t0 = modalTimer->getTime();
      computeModalDisp(sdt, deltmp, delW, delU, delY, refModalF);
      modalTimer->addStructUpdTime(t0);

    }
  
    // compute updated position
    deltmp = 0.0;
    for (i = 0; i < nStrMode; ++i)
      deltmp += (delU[i]+sdt/2*delY[i])*mX[i];
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

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, cntp1*dt, com->cpuNum(), tRestart->energy, deltmp, delWtmp);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, cntp1*dt, deltmp, controlVol, delWtmp, tState);

    if (ierr)  {
      com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
      cntr++;
      if (cntr > 4)  {
        com->barrier();
        exit(-1);
      }
    }

    if (cnt % 20 == 0)
      com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);

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
ModalSolver<dim>::timeIntegrateROM(VecSet<Vec<double> > &romOp,
                  VecSet<Vec<double> > &ecMat, VecSet<Vec<double> > &gMat,
                  VecSet<DistSVec<double, dim> > &podVecs, int nSteps, int nPodVecs,
                  double *delU, double *delY, double sdt)  {

  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 

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

  double *prevU = new double[nStrMode];
  double *prevY = new double[nStrMode];

  int i;
  for (i = 0; i < nStrMode; ++i)
    deltmp += (delU[i]+0.5*sdt*delY[i])*mX[i];

  deltmp += Xref;
  // Init delW
  char *nlSolFile = 0;
  if (ioData->input.perturbed[0] == 0)
    nlSolFile = tInput->solutions;
  else  {
    int sp = strlen(ioData->input.prefix) + 1;
    nlSolFile = new char[sp + strlen(ioData->input.perturbed)];
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

  for (int cnt = 0; cnt < nSteps; ++cnt) {

    cntp1 = cnt+1;
    rhs = 0.0;
    rhsTemp = 0.0;


    delWRomTemp = 4.0*delWRom;
    delWRomTemp -= delWnm1Rom;

    delWRomTemp *= 1.0/dt; 
    for (i = 0; i < nPodVecs; ++i)  
      rhs += delWRomTemp[i]*romOp[i];

    // not exactly the right terms... see ::timeIntegrate for full comment
    for (i = 0; i < nStrMode; ++i)
      rhsTemp += prevY[i]*ecMat[i] - 3.0*delY[i]*ecMat[i] - 2.0*delU[i]*gMat[i];
      //rhsTemp += prevU[i]*gMat[i] - 3.0*delY[i]*ecMat[i] - 3.0*delU[i]*gMat[i] + prevY[i]*ecMat[i];

    // compute new delWRom
    delWnm1Rom = delWRom;
    delWRom = rhs + rhsTemp;

    // project soltn into full space
    delWFull = 0.0;
    for (i = 0; i < nPodVecs; ++i)
      delWFull += podVecs[i] * delWRom[i];

    for (int kk = 0; kk < nStrMode; kk++)  {
      prevU[kk] = delU[kk];
      prevY[kk] = delY[kk];
    }
    computeModalDisp(sdt, deltmp, delWFull, delU, delY, refModalF);

    // compute Cl, Cm
    deltmp = 0.0;
    for (i = 0; i < nStrMode; ++i)
      deltmp += (delU[i]+sdt/2*delY[i])*mX[i];

    deltmp += Xref;
    delWFull += Uref;

    int ierr = domain.checkSolution(varFcn, delWFull);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, cntp1*dt, com->cpuNum(), tRestart->energy, deltmp, delWFull);

    if (ierr)  {
      com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
      exit(-1);
    }

    if (cnt % 20 == 0)
      com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);

  }
}

//--------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::preProcess()  {

  // setup solvers
  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 
  geoState = new DistGeoState(*ioData, &domain);
  geoState->setup1(tInput->positions, &Xref, &controlVol);
  bcData = new DistBcDataEuler<dim>(*ioData, varFcn, &domain, Xref);

  spaceOp = new SpaceOperator<dim>(*ioData, varFcn, bcData, geoState, &domain);
  postOp = new PostOperator<dim> (*ioData, varFcn, bcData, geoState, &domain);

  // Temporal operator contains Uref for us
  tState = new DistTimeState<dim>(*ioData, spaceOp, varFcn, &domain);
  tState->setup(tInput->solutions,  bcData->getInletBoundaryVector(), Xref, Uref);

  RefVal *refVal = new RefVal(ioData->ref.rv);
  tOutput = new TsOutput<dim>(*ioData, refVal, &domain, postOp);
  tRestart = new TsRestart(*ioData, refVal);

  HOp = new MatVecProdH2<double, 5>(*ioData,  varFcn, tState, spaceOp, &domain);
  HOp2 = new MatVecProdH2<double,5>(*ioData,  varFcn, tState, spaceOp, &domain);

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

  /*ksp = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com);

  if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
    HOpC = new MatVecProdH2<bcomp,5>(*ioData,  varFcn, tState, spaceOp, &domain);

    pcComplex = new IluPrec<bcomp ,dim, bcomp>(pcData, &domain);

    kspComp = new GmresSolver<DistSVec<bcomp,dim>, MatVecProd<dim, 5>,
                  KspPrec<dim, bcomp>, Communicator, bcomp>
                  ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com);

  }
  */

  if (kspData.type == KspData::GMRES)
    ksp = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com);
  else if (kspData.type == KspData::GCR)
    ksp = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, 5>,
        KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com);
                                                        
  if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
    HOpC = new MatVecProdH2<bcomp,5>(*ioData,  varFcn, tState, spaceOp, &domain);
                                                        
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
  double sdt = ioData->linearizedData.eps2;
  dt = ioData->linearizedData.stepsize/ioData->ref.rv.time;
  com->fprintf(stderr, " ... Running Modal Solver w/%d struct modes, coupling timestep: %f and fluid step: %f, RefTime = %f (%f/%f)\n", nStrMode, sdt, dt, ioData->ref.rv.time, ioData->ref.rv.length, ioData->ref.rv.velocity);

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
    // Computing F(X0,V1) sdt is dimensional
    double bigeps = 0.5*eps/sdt;
    //double bigeps = 0.5*eps/ioData->ref.rv.time;
        
    Xnp1 = Xref - sdt*bigeps*(mX[ic]);
    
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    geoState->update(Xnp1, CV1);
    
    Xnp1 = Xref + sdt*bigeps*(mX[ic]);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

    spaceOp->computeResidual(Xnp1, CV1, Uref, F1, tState);
    tState->add_dAW_dt(1,*geoState,CV1,Uref,F1);
    spaceOp->applyBCsToResidual(Uref, F1);
 
    // Computing F(X0, -V1)
    Xnp1 = Xref + sdt*bigeps*(mX[ic]);
  
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
    geoState->update(Xnp1, CV1);

    Xnp1 = Xref - sdt*bigeps*(mX[ic]);
    geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

    spaceOp->computeResidual(Xnp1, CV1, Uref, F2, tState);
    tState->add_dAW_dt(1,*geoState,CV1,Uref,F2);
    spaceOp->applyBCsToResidual(Uref, F2);


    // no dt here because the spatial computation will divide the
    // given displacement change by dt to define its velocity.
    // Thus with Xnp1 = X0+dt*eps*Xm and Xn = X0-dt*eps*Xm: V_n+1/2 = 2*eps*Xm
    
    DV[ic] = (1.0/(4*bigeps))*(F1-F2);

    com->fprintf(stderr, " ... Norm DV, Mode %i: %e\n",ic, DV[ic].norm());

    // create E+C
    DE[ic] += DV[ic];
    com->fprintf(stderr, " ... Norm DE, Mode %i: %e\n\n",ic, DE[ic].norm());
  }

  // Now setup H matrices

  // ***This is (A+dt/2*H)
  HOp->evaluate(0, Xref, controlVol, Uref, FF);

  //***This is (A-dt/2*H)
  HOp2->evaluate2(0, Xref, controlVol, Uref, FF);

  //set up preconditioner, GMRES Solver
  DistMat<PreScalar,dim> *_pc = dynamic_cast<DistMat<PreScalar,dim> *>(pc);

  if (_pc) {
    spaceOp->computeH1(Xref, controlVol, Uref, *_pc);
    tState->addToH1(controlVol, *_pc);
    spaceOp->applyBCsToJacobian(Uref, *_pc);
  }

  pc->setup();
  ksp->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);

  tOutput->openAsciiFiles();
}

//-------------------------------------------------------------------------------

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

//-------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::constructROM2(VecSet<Vec<double> > &romOperator,
                       VecSet<Vec<double> > &ecVecs, VecSet<Vec<double> > &gVecs,
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {
#ifdef DO_MODAL
  com->fprintf(stderr, " ... Forming Reduced System Directly \n");

  // form fluid ROM
  DistSVec<double,dim> FF(domain.getNodeDistInfo());

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 
  MatVecProdH2<double, dim> *onlyHOp = new MatVecProdH2<double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  onlyHOp->evalH(0, Xref, controlVol, Uref);

  double invDt = 1.0/dt;
  double invDt3 = 3.0/dt;
  int iVec, jVec;

  // form H ROM op
  DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
  
  double *romOpPlusVals = new double[nPodVecs*nPodVecs];
  VecSet<Vec<double> > romOpMinus(nPodVecs, nPodVecs);
  VecSet<Vec<double> > romOpPlusInv(nPodVecs, nPodVecs);

  double romVal;
  com->barrier();
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    onlyHOp->apply(podVecs[iVec], tmpVec);
    //com->fprintf(stderr," iVec = %d norm = %e\n",iVec,tmpVec.norm());
    tmpVec /= controlVol;
    for (int jVec = 0; jVec < nPodVecs; jVec++)  {
      romVal = 2.0*(podVecs[jVec] * tmpVec);
      romOpPlusVals[iVec*nPodVecs+jVec] = romVal;
      //romOpMinus[iVec][jVec] = -romVal;
    }

    romOpPlusVals[iVec*nPodVecs+iVec] += invDt3;
    romOpMinus[iVec] = 0.0;
    romOpMinus[iVec][iVec] = 1.0;

    //spaceOp has to be reset because it has been modified by the apply function
    DistSVec<double,dim> FF(domain.getNodeDistInfo());
    spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);  
  }
  ARdsNonSymMatrix<double, double> romOpPlus(nPodVecs, romOpPlusVals);
  romOpPlus.FactorA();

  com->fprintf(stderr, " ... Factored Matrix\n");

  Vec<double> res(nPodVecs);

  com->fprintf(stderr, " ... Forming Reduced ROM Operator\n");

  for (iVec = 0; iVec < nPodVecs; iVec++)
    romOpPlus.MultInvv(romOpMinus[iVec].data(), romOperator[iVec].data());

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

    romOpPlus.MultInvv(tmpECrom.data(), ecVecs[iVec].data());
    romOpPlus.MultInvv(tmpGrom.data(), gVecs[iVec].data());

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
  bcomp kImag(0.0, kFreq);
                                                        
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  DistSVec<double, dim> delWReal(domain.getNodeDistInfo());
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
  DistSVec<double, dim> delWReal(domain.getNodeDistInfo());
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

#ifdef DO_MODAL
  int nPOD = ioData->linearizedData.numPOD;
  int nconv = 0;

  com->fprintf(stderr, " ... Forming %d x %d Correlation Matrix\n", nSnaps, nSnaps);
  double t0 = modalTimer->getTime();

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

#else
  com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
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
void
ModalSolver<dim>::interpolatePOD()  {

  Timer *modalTimer = domain.getTimer();
  int nData = 2;
  com->fprintf(stderr, " ... Interpolating POD\n");
  char **podFile = new char *[nData];

  char *vecFile = tInput->podFile;
  if (!vecFile)  vecFile = "podFiles.in";
  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
    exit (-1);
  }

  char podFile1[500], podFile2[500]; 
  fscanf(inFP, "%s %s", podFile1, podFile2);
  podFile[0] = podFile1;
  podFile[1] = podFile2;
  com->fprintf(stderr, " ... Reading POD from %s and %s\n", podFile[0], podFile[1]);

  if (ioData->output.transient.podFile[0] == 0)  {
    com->fprintf(stderr, "*** ERROR: POD Basis File not specified\n");
    exit (-1);
  }
  int sp = strlen(ioData->output.transient.prefix);
  char *outFile = new char[sp + strlen(ioData->output.transient.podFile)];
  sprintf(outFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);

  int numPod = ioData->linearizedData.numPOD;

  double *mach = new double[nData];
  double newMach = ioData->ref.mach;
  fscanf(inFP, "%lf  %lf  %lf", mach+0, mach+1, &newMach);

  com->fprintf(stderr, " ... Interpolating new Mach %f between %f and %f\n", newMach, mach[0], mach[1]);

  VecSet< DistSVec<double, dim> > pod1(numPod, domain.getNodeDistInfo());
  VecSet< DistSVec<double, dim> > pod2(numPod, domain.getNodeDistInfo());

  // read in Pod Vectors

// Comment/Uncomment here
  int iPod;
  int j, k;
  double *eig1 = new double[numPod];
  double *eig2 = new double[numPod];

  domain.readVectorFromFile(podFile[0], 0, &eig1[0], pod1[0]);
  if (numPod > eig1[0])  {
    com->fprintf(stderr, "*** Warning: Resetting number of interpolated POD vectors from %d to %d\n", numPod, (int) eig1[0]);
    numPod = (int) eig1[0];
  }
  domain.readVectorFromFile(podFile[1], 0, &eig2[0], pod2[0]);
  if (numPod > eig2[0])  {
    com->fprintf(stderr, "*** Warning: Resetting number of interpolated POD vectors from %d to %d\n", numPod, (int) eig2[0]);
    numPod = (int) eig2[0];
  }

  for (iPod = 0; iPod < numPod; iPod++)  {
    domain.readVectorFromFile(podFile[0], iPod+1, &eig1[iPod], pod1[iPod]);
    domain.readVectorFromFile(podFile[1], iPod+1, &eig2[iPod], pod2[iPod]);
  }

  double t0 = modalTimer->getTime();
  // compute principle angles and vectors between pod subspaces
  double *matVals = new double[numPod*numPod];
  for (j = 0; j < numPod; j++)
    for (k = 0; k < numPod; k++)
      matVals[j*numPod + k] = pod1[k] * pod2[j];

  double *sigma = new double[numPod];
  double *error = new double[numPod];
  double *work = new double[numPod];
  int info;
  double *zVec = new double[numPod*numPod]; // right singular vectors
  double *yVec = new double[numPod*numPod]; // left singular vectors

  F77NAME(dsvdc)(matVals, numPod, numPod, numPod, sigma, error,
                 yVec, numPod, zVec, numPod, work, 11, info);
  com->fprintf(stderr, " ... Computed SVD with code %d for mach %f\n", info, newMach);

  // compute principle vectors and interpolate on angles
  double machCoef = (newMach-mach[0])/(mach[1]-mach[0]);
  DistSVec<double, dim> U(domain.getNodeDistInfo());  // left principle vector
  Vec<double> Utmp(numPod);  
  DistSVec<double, dim> V(domain.getNodeDistInfo());  // right principle vector
  Vec<double> Vtmp(numPod);  
  double pAngle;  // interpolated principle angle

  VecSet<DistSVec<double, dim> > pod(numPod, domain.getNodeDistInfo());
  DistSVec<double, dim> orthoVec(domain.getNodeDistInfo());

  com->barrier();
  for (j = 0; j < numPod; j++)  {

    U = 0.0;  V = 0.0;
    pAngle = acos(sigma[j]) * machCoef;
    com->fprintf(stderr, " ... Interpolating vec %d to angle %e (%e)\n", j, pAngle, acos(sigma[j])); 
     
    for (k = 0; k < numPod; k++)  {
      Utmp[k] = yVec[j*numPod+k]; 
      Vtmp[k] = zVec[j*numPod+k]; 
    }
    for (k = 0; k < numPod; k++)  {
      U += pod1[k] * Utmp[k];
      V += pod2[k] * Vtmp[k];
    }

    orthoVec = V - (U*V)*U;
    orthoVec *= 1.0/orthoVec.norm();
    pod[j] = U*cos(pAngle) + orthoVec*sin(pAngle);
  }

  com->fprintf(stderr, " ... Writing Interpolated POD to %s\n", outFile);
  // output interpolated pod vectors
  DistSVec<double, dim> newPOD(domain.getNodeDistInfo());
  domain.writeVectorToFile(outFile, 0, numPod, newPOD); 
  double newEig;
  for (j = 0; j < numPod; j++)  {
    newPOD = 0.0;
    for (k = 0; k < numPod; k++)
      newPOD += pod[k] * zVec[k*numPod+j];

    newEig = machCoef*(eig2[j]-eig1[j]);
    com->fprintf(stderr, " ... Writing pod vec %d with newEig %f(%f-%f)\n", j+1, newEig, eig2[j],eig1[j]);
    domain.writeVectorToFile(outFile, j+1, newEig, newPOD); 
  }
  modalTimer->addMeshSolutionTime(t0);

// Comment/Uncomment here
// THis block of code is for testing purposed ONLY.  It does interpolation using Lagranges method
// directly on  the POD basis vectors and not the subspace angles.
/*

  VecSet<DistSVec<double, dim> > newPOD(numPod, domain.getNodeDistInfo());
  VecSet< DistSVec<double, dim> > pod(nData, domain.getNodeDistInfo());
  DistSVec<double, dim>  testPod(domain.getNodeDistInfo());
  DistSVec<double, dim>  switchPod(domain.getNodeDistInfo());

  double *eig = new double [nData];
  double newEig;
  domain.writeVectorToFile(outFile, 0, numPod, newPOD[0]); 
  int iPod, j, k;

  for (iPod = 0; iPod < numPod; iPod++)  {

    for (j = 0; j < nData; j++)
      domain.readVectorFromFile(podFile[j], iPod+1, eig+j, pod[j]);

    for (j = 1; j < nData; j++)  {
      testPod = pod[j] - pod[j-1];
      if (testPod.norm() > 1.0)  {
        double refNorm = testPod.norm();
        pod[j] *= -1.0;
        testPod = pod[j] - pod[j-1];
        if (testPod.norm() > refNorm)
          pod[j] *= -1.0;
      }
    }

    newPOD[iPod] = 0.0;
    newEig = 0.0;

    double dM = 1.0;
    for (j = 0; j < nData; j++)  {

      for (k = 0; k < nData; k++)  {
        if (k == j) continue;
        dM *= (newMach - mach[k]) / (mach[j] - mach[k]);
      }
      newPOD[iPod] += pod[j]*dM;
      newEig += eig[j]*dM;
    }

    for (j = 0; j < iPod; j++)
      newPOD[iPod] -= (newPOD[iPod] * newPOD[j])*newPOD[j];
       
    newPOD[iPod] *= 1.0/newPOD[iPod].norm();

    domain.writeVectorToFile(outFile, iPod+1, newEig, newPOD[iPod]); 
  }
*/
  //modalTimer->setRunTime();
  //modalTimer->print(modalTimer);

}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::evalFluidSys(VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 
  MatVecProdH2<double, dim> *onlyHOp = new MatVecProdH2<double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
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
  char *romFile = new char[sp + strlen(ioData->output.transient.romFile)];
  sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
  FILE *romFP = fopen(romFile, "w");
  com->barrier();

  com->fprintf(romFP, "%d 0\n", nPodVecs);
  for (iVec = 0; iVec < nPodVecs; iVec++)  {
    for (jVec = 0; jVec < nPodVecs; jVec++)
      com->fprintf(romFP, "%e ", romOperator[iVec][jVec]);
    com->fprintf(romFP, "\n ");
  }

  //spaceOp has to be reset because it has been modified by the apply function
  DistSVec<double,dim> FF(domain.getNodeDistInfo());
  spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);
  
}

//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::evalAeroSys(VecSet<Vec<double> > &outRom, 
                       VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

  // Allocate ROM operators
  VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData); 
  MatVecProdH2<double, dim> *onlyHOp = new MatVecProdH2<double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
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
    //tmpVec = DE[iVec] * invDt; #BUG#
    tmpVec = DE[iVec];
    tmpVec2 = DX[iVec];
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
      sysVals[sysSize*nPodVecs + iVec*sysSize + jVec] = -ecVecs[iVec][jVec];
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
  char *romFile = new char[sp + strlen(ioData->output.transient.romFile)];
  sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
  FILE *romFP = fopen(romFile, "w");
  com->barrier();

  com->fprintf(romFP, "%d %d\n", nPodVecs, nStrMode);
  for (iVec = 0; iVec < sysSize; iVec++)  {
    for (jVec = 0; jVec < sysSize; jVec++)
      com->fprintf(romFP, "%e ", sysVals[jVec*sysSize+iVec]);
    com->fprintf(romFP, "\n ");
  }

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
  //const char *podFileName = ioData->output.transient.podFile;

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


