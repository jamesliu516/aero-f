#include <cstdio>
#include <cmath>
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
//#include <VarFcnDesc.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <PostOperator.h>
#include <ParallelRom.h>

#ifdef DO_MODAL
#include <arpack++/include/ardsmat.h>
//#include <arpack++/include/ardnsmat.h>
#include <arpack++/include/ardssym.h>
#endif

#include <cstring>
extern "C"      {

	void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
				               double *, double *, int &, double *, int &,
											 double *, const int &, int &);

  void F77NAME(thinsvd)(int &, int &, int &, int &, int&, int&, int&, int&, int&, double *, int &,
						           int &, int &, int &, int &, int &, int &, double *U, double *S, double *V,
											int &, double *, int &, int &);
	void F77NAME(lworksize)(int &, int &, int &, int &, int&, int&, int&, int&, int&, int &);

}

//----------------------------------------------------------------------------------

template <int dim>
ModalSolver<dim>::ModalSolver(Communicator *_com, IoData &_ioData, Domain &dom) : 
          domain(dom), mX(0, dom.getNodeDistInfo() ), Xref(dom.getNodeDistInfo()), 
          Uref(dom.getNodeDistInfo()), DX(0, dom.getNodeDistInfo()), 
          DE(0, dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()), 
					controlVolComp(dom.getNodeDistInfo())  {

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
 totalEnergy = 0.0;	// for POD basis construction

}

//-------------------------------------------------------------------------------
template <int dim> 
void ModalSolver<dim>::solve()  {

 // set up Timers
 Timer *modalTimer = domain.getTimer();
 double t0;

 const char *snapsFile = tInput->snapFile;
 modalTimer->setSetupTime();//CBM

 if (ioData->problem.alltype == ProblemData::_INTERPOLATION_)
   interpolatePOD();
 else if (ioData->problem.alltype == ProblemData::_ROB_CONSTRUCTION_ && snapsFile){
   t0 = modalTimer->getTime(); //CBM--check
   buildGlobalPOD();
   modalTimer->addPodConstrTime(t0);
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_){
	 geoState = new DistGeoState(*ioData, &domain);
	 geoState->setup1(tInput->positions, &Xref, &controlVol);
	 GnatPreprocessing<dim> gappy(com,*ioData,domain,geoState);
	 gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_1_){
	 geoState = new DistGeoState(*ioData, &domain);
	 geoState->setup1(tInput->positions, &Xref, &controlVol);
	 GnatPreprocessingStep1<dim> gappy(com,*ioData,domain,geoState);
	 gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_2_){
	 geoState = new DistGeoState(*ioData, &domain);
	 geoState->setup1(tInput->positions, &Xref, &controlVol);
	 GnatPreprocessingStep2<dim> gappy(com,*ioData,domain,geoState);
	 gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_SURFACE_MESH_CONSTRUCTION_){
	 geoState = new DistGeoState(*ioData, &domain);
	 geoState->setup1(tInput->positions, &Xref, &controlVol);
	 SurfMeshGen<dim> surfMeshBuilder(com,*ioData,domain,geoState);
	 surfMeshBuilder.buildReducedModel();	
 }
 else if (ioData->problem.alltype == ProblemData::_SAMPLE_MESH_SHAPE_CHANGE_){
	 //KTC: CHANGE!!!
	 geoState = new DistGeoState(*ioData, &domain);
	 geoState->setup1(tInput->positions, &Xref, &controlVol);
	 ReducedMeshShapeChanger<dim> reducedMeshShapeChanger(com,*ioData,domain,geoState);
	 reducedMeshShapeChanger.buildReducedModel();	
 }
 else  {
   preProcess();

   //modalTimer->setSetupTime();
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

/*
 // set up Timers
 Timer *modalTimer = domain.getTimer();
 double t0;

 if (ioData->problem.alltype == ProblemData::_INTERPOLATION_)
   interpolatePOD();
 else  {
   preProcess();
   modalTimer->setSetupTime();
   if (ioData->problem.alltype == ProblemData::_ROB_CONSTRUCTION_) {
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
*/
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
     com->fprintf(stderr, " ... Running Forced Oscillations w/no initializations\n");
   else if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_)
     com->fprintf(stderr, " ... Running Unsteady Linearized w/no structural initializations\n");
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

    delete [] romOperator;
    delete [] romOperator1;
    delete [] romOperator2;

 }
  delete [] delY;
  delete [] delU;
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

 //  set up Snapshot Matrix
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
 double *prevU = new double[nStrMode];
 double *prevY = new double[nStrMode];
  double *prevA = new double[nStrMode];

 for (i = 0; i < nStrMode; ++i)  {
    deltmp += (delU[i]+sdt*delY[i])*mX[i];
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

  //Compute initial acceleration 
  Vec<double> modalF(nStrMode);
  postOp->computeForceDerivs(Xref, Uref, delW, modalF, mX);
  modalF += refModalF;
  modalF *= ioData->ref.rv.force;

  for (i = 0; i < nStrMode; i++)
    prevA[i] = modalF[i] - K[i]*prevU[i];

 int cntr = 0;
 int cntp1;

 // every step corresponds to solving for 2x each member of the BDF scheme //
 
 FF *= (-2.0); 

 for (int cnt = 0; cnt < nSteps; ++cnt) {

   cntp1 = cnt+1;

   t0 = modalTimer->getTime();  
   if (cnt == 0) {

     rhsA = 0.0;
     rhsB = 0.0;
     rhsC = 0.0;
     rhsD = 0.0;
      rhs  = 0.0;
     for (i = 0; i < nStrMode; ++i)  {
       rhsA += DX[i]*delU[i];
       rhsB += DX[i]*delY[i];
       rhsC += DE[i]*delY[i];
        rhsD += DE[i]*prevA[i];
     }

     //Predictor step
      HOp2step1->apply(delWint,rhs);
      rhs+= rhsA + rhsC - 0.5*FF;
      rhs*= 0.5*dt;
      rhs /= controlVol;
      delW = delWint - rhs;
     
      rhs = 0.0;
      //Corrector step 
      HOp2step1->apply(delW,rhs);
      rhs+= rhsA + 0.5*sdt*rhsB + rhsC + 0.5*sdt*rhsD - 0.5*FF;
      rhs*= dt;
      rhs /= controlVol;
      delW = delWint - rhs;
   }
   else if (cnt == 1) {
      
      rhs = delW*(6.0/dt) - delWnm1*(2.0/(3.0*dt));
     rhs *= controlVol;

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

     ksp2->solve(rhs, delW);
     modalTimer->addKspTime(t0);
   }
   else if (cnt == 2) {
  
      rhs = delW*(6.0/dt) - delWnm1*(8.0/(3.0*dt));
     rhs *= controlVol;

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
     
      rhs = delW*(4.0/dt) - delWnm1*(1.0/dt);
     rhs *= controlVol;

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

     ksp->solve(rhs, delW);
     modalTimer->addKspTime(t0);
   }

   // for forced oscillations          
   if (ioData->linearizedData.type == LinearizedData::FORCED)  {
     if (cnt <= 1) {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*(cnt)*sdt)*dMax;   
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*(cnt)*sdt)*dMax;
     } 
     else {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*((cnt)*sdt))*dMax;
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*((cnt)*sdt))*dMax;  
     }
   }

   if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {

     for (int kk = 0; kk < nStrMode; kk++)  {
       prevU[kk] = delU[kk];
       prevY[kk] = delY[kk];
     }
     t0 = modalTimer->getTime();


      computeModalDisp(sdt, deltmp, delW, delU, delY, refModalF, cnt);
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

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, (cnt+1)*dt, com->cpuNum(), tRestart->energy, deltmp, delWtmp);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, (cnt+1)*dt, deltmp, controlVol, delWtmp, tState);

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

  for (i = 0; i < nStrMode; ++i) {
    delU[i] = 0.5*(delU[i] + prevU[i]);
    delY[i] = 0.5*(delY[i] + prevY[i]);
  }
  if (!(ioData->input.perturbed[0] == 0))
    delete [] nlSolFile;
  delete [] prevU;
  delete [] prevY;
  delete [] prevA;
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
  Vec<double> prevWRom(nPodVecs);
  Vec<double> pprevWRom(nPodVecs);
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
   nlSolFile = new char[sp + strlen(ioData->input.perturbed)+1];
   sprintf(nlSolFile, "%s%s", ioData->input.prefix, ioData->input.perturbed);
 }
 domain.readVectorFromFile(nlSolFile, 0, 0, delWFull);
 com->fprintf(stderr, " ... Read Perturbed solution: W = %e, Uref = %e\n", delWFull.norm(), Uref.norm());
 delWFull -= Uref;

 // construct initial delWRom
 for (i = 0; i < nPodVecs; i++)
   delWRom[i] = podVecs[i] * delWFull;

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

  pprevWRom = delWRom; 
  prevWRom = delWRom;

 
  //Precompute P*PHI
  VecSet<Vec<double> > PtimesPhi(nPodVecs, nStrMode);

  for (int iVec = 0; iVec < nPodVecs; iVec++)  {
    postOp->computeForceDerivs(Xref, Uref, podVecs[iVec], modalF, mX);
    PtimesPhi[iVec] = ioData->ref.rv.force*modalF;
  }

  //Time integration loop
  for (int cnt = 0; cnt < nSteps+1; ++cnt) {

    cntp1 = cnt+1;

    if (cnt == 0){
      delWRomTemp = prevWRom;
      for (i = 0; i < nPodVecs; ++i)
        delWRomTemp += 0.5*dt*delWRom[i]*romOp0[i];

      for (i = 0; i < nStrMode; ++i){
        delWRomTemp -= 0.5*dt*(delU[i]*gMat[i] + delY[i]*ecMat[i]);
        delWRom -= dt*((delU[i] + 0.5*sdt*delY[i])*gMat[i] + delY[i]*ecMat[i]);
      }

     for (i = 0; i < nPodVecs; ++i)
        delWRom += dt*delWRomTemp[i]*romOp0[i];
   }
  
    else if (cnt == 1) {
     delWRom = 0.0;

     for (i = 0; i < nPodVecs; ++i)
      delWRom += (1/dt)*(3*prevWRom[i] - (1/3)*pprevWRom[i])*romOperator[i];

      for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat[i] );
    } 
 
    else if (cnt == 2) {
      delWRom = 0.0;
     for (i = 0; i < nPodVecs; ++i)
        delWRom += (1/dt)*(3*prevWRom[i] - (4/3)*pprevWRom[i])*romOperator1[i];
 
     for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat1[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat1[i] );
   }

   else {
      delWRom = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWRom += (1/dt)*(2*prevWRom[i] - 0.5*pprevWRom[i])*romOperator2[i];
 
      for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat2[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat2[i] );
    }

    pprevWRom = prevWRom;
    prevWRom = delWRom;

    // project soltn into full space
    delWFull = 0.0;
    for (i = 0; i < nPodVecs; ++i)
      delWFull += podVecs[i] * delWRom[i];


   if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_){
     for (int kk=0; kk < nStrMode; kk++){
       prevU[kk] = delU[kk];
       prevY[kk] = delY[kk];
     }
   }

    // Output the reduced order vector
    computeModalDisp(sdt, delWRom, delU, delY, refModalF, PtimesPhi, nPodVecs, cnt);
   
   // compute Cl, Cm
   deltmp = 0.0;
   for (i = 0; i < nStrMode; ++i)
     deltmp += (delU[i]+sdt/2*delY[i])*mX[i];

   deltmp += Xref;
   delWFull += Uref;

   int ierr = domain.checkSolution(varFcn, delWFull);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, (cnt+1)*dt, com->cpuNum(), tRestart->energy, deltmp, delWFull);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, (cnt+1)*dt, deltmp, controlVol, delWFull, tState);

   if (ierr)  {
     com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
     exit(-1);
   }

   if (cnt % 20 == 0)
     com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);
  }

  for (i = 0; i < nStrMode; ++i) {
    delU[i] = 0.5*(prevU[i] + delU[i]);
    delY[i] = 0.5*(prevY[i] + delY[i]);
  } 
  delete [] prevU;
  delete [] prevY;
  if (!(ioData->input.perturbed[0] == 0))
    delete [] nlSolFile;

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

 HOp = new MatVecProdH2<dim,double, dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOp2 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOp2step1 = new MatVecProdH2<dim,double, dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOpstep2 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOpstep3 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 

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
   ksp = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com); 
   ksp2 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com); 
   ksp3 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>, 
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com); 
   }
 else if (kspData.type == KspData::GCR) {
   ksp = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com); 
   ksp2 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com); 
   ksp3 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com); 
 }

 if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
   HOpC = new MatVecProdH2<dim,bcomp,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 

   pcComplex = new IluPrec<bcomp ,dim, bcomp>(pcData, &domain);
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {


    kspCompGcr = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim,dim>,
                   KspPrec<dim, bcomp>, Communicator, bcomp>
                   ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com); 
   }
   else {

     if (kspData.type == KspData::GMRES)

       kspComp = new GmresSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>,
                     KspPrec<dim, bcomp>, Communicator, bcomp>
                    ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com); 
     else if (kspData.type == KspData::GCR)

       kspComp = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>,
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

  tState->setGlobalTimeStep(dt);
  HOp2step1->evalH(0, Xref, controlVol, Uref); 

  tState->setGlobalTimeStep(dt);
  HOpstep2->evaluate(8, Xref, controlVol, Uref, FF,0.0);

 tState->setGlobalTimeStep(dt);
  HOpstep3->evaluate(9, Xref, controlVol, Uref, FF,0.0);

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
  //MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); // PJSA 11/04/2010
 onlyHOp->evalH(0, Xref, controlVol, Uref);

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
      romOpPlusVals[iVec*nPodVecs+jVec] = romVal;
      romOpPlusVals1[iVec*nPodVecs+jVec] = romVal;
      romOpPlusVals2[iVec*nPodVecs+jVec] = romVal;
     romOperator0[iVec][jVec] = -romVal;
   }

    romOpPlusVals[iVec*nPodVecs+iVec] += 8/(3*dt);
    romOpPlusVals1[iVec*nPodVecs+iVec] += 5/(3*dt);
    romOpPlusVals2[iVec*nPodVecs+iVec] += 3/(2*dt);

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
void ModalSolver<dim>::computeModalDisp(double sdt, Vec<double> &delWRom, double *delU, double *delY, Vec<double> &refModalF, VecSet<Vec<double> > &PtimesPhi, int nPodVecs, int timeIt) {
//For ROM
Vec<double> modalF(nStrMode);
modalF = ioData->ref.rv.force*refModalF;
for (int i = 0; i < nPodVecs; i++)
  modalF += delWRom[i]*PtimesPhi[i];


  updateModalValues(sdt, delU, delY, modalF, timeIt);
}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::computeModalDisp(double sdt, DistSVec<double, 3> &xPos, DistSVec<double, dim> &delW, double *delU, double *delY, Vec<double> &refModalF, int timeIt) {
//For FOM

 Vec<double> modalF(nStrMode);
 postOp->computeForceDerivs(xPos, Uref, delW, modalF, mX);
 modalF += refModalF;
 modalF *= ioData->ref.rv.force;

  updateModalValues(sdt, delU, delY, modalF, timeIt);
}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::updateModalValues(double sdt, double *delU, double *delY, Vec<double> &modalF, int timeIt){
//For both ROM and FOM
  double prevU, prevY;
 double *srhs = new double[nStrMode];
double sdt2 = sdt*sdt;

 if (timeIt == 0){
    for (int i = 0; i < nStrMode; i++) {
       prevU = delU[i];
       prevY = delY[i];
       srhs[i] =  modalF[i] - K[i]*prevU + (1/sdt - 0.5*sdt*K[i])*prevY;
       delY[i] = srhs[i]/(1/sdt + 0.5*sdt*K[i]);
       delU[i] = 0.5*sdt*(prevY + delY[i]) + prevU;
    }
  }else {
    for (int i = 0; i < nStrMode; i++)  {
       prevU = delU[i];
       prevY = delY[i];
       srhs[i] =  modalF[i] - K[i]*prevU + (1/sdt - 0.25*sdt*K[i])*prevY;
       delY[i] = srhs[i]/(1/sdt + 0.25*sdt*K[i]);
       delU[i] = 0.5*sdt*(prevY + delY[i]) + prevU;
     }
  }
} 
//---------------------------------------------------------------------------------

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
void ModalSolver<dim>::makeFreqPOD(VecSet<DistSVec<double, dim> > &snaps, int nSnaps, int nPOD = 0, bool outputToDisk = true){

 Timer *modalTimer = domain.getTimer();

 if (nPOD == 0)
	 nPOD = (ioData->linearizedData.numPOD) ? ioData->linearizedData.numPOD : ioData->rom.dimension; // CBM--check   
 if (nPOD == 0)
	 nPOD = nSnaps;

 int nconv = 0;
 if (nPOD > nSnaps) {
	 com->fprintf(stderr, " ... WARNING: nPOD > nSnaps. Changing nPOD from %d to %d\n", nPOD, nSnaps);
	 nPOD = nSnaps;
 }

 com->fprintf(stderr, " ... Computing %d POD vectors from %d snapshots\n", nPOD, nSnaps);

 if (podMethod == 0) {	// svd
#ifdef DO_SCALAPACK
	 VecSet<DistSVec<double, dim> > Utrue(nSnaps, domain.getNodeDistInfo());
	 Vec<double> singVals(nSnaps);
	 FullM VtrueDummy(1);	// do not need VtrueDummy

	 double t0 = modalTimer->getTime();
	 ParallelRom<dim> parallelRom(domain,com);
	 parallelRom.parallelSVD(snaps, Utrue, singVals.data(), VtrueDummy, nSnaps, false);
	 modalTimer->addEigSolvTime(t0);
	 if (outputToDisk)
		 outputPODVectors(Utrue, singVals, nPOD);
	 else {	// overwrite snaps with Utrue and return
		 for (int i = 0; i < nPOD; ++i) {
			 snaps[i] = Utrue[i]*singVals[i];
		 }
	 }
#else
 com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH SCALAPACK \n");
 exit(-1);

#endif
	}
 else if (podMethod == 1) {	// eig
#ifdef DO_MODAL

  // allocate for upper half of sym. eigprob
  double *rVals = new double[nSnaps*(nSnaps+1)/2];
  for (int i = 0; i < nSnaps; i++){
    com->fprintf(stderr," ... processing snap %d\n",i);//CBM
    for (int j = 0; j <= i; j++)
      rVals[(i+1)*i/2 + j] = snaps[j] * snaps[i];
  }

  double tolerance = ioData->snapshots.dataCompression.tolerance;

	com->barrier();
  ARdsSymMatrix<double> pod(nSnaps, rVals, 'U');
  com->fprintf(stderr, " ... Factoring Correlation Matrix\n");

	double t0 = modalTimer->getTime();
	int iSnap;
	pod.FactorA();
  ARluSymStdEig<double> podEigProb(nPOD, pod, "LM", nSnaps-1, tolerance, 300*nPOD);
  modalTimer->addCorrelMatrixTime(t0);

  t0 = modalTimer->getTime();
  nconv = podEigProb.FindEigenvectors();
  modalTimer->addEigSolvTime(t0);

	com->fprintf(stderr, " ... Got %d converged eigenvectors out of %d snaps\n", nconv, nSnaps);
	outputPODVectors(podEigProb, snaps, nPOD, nSnaps);
	delete [] rVals;

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);

#endif
 }
}

//---------------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::buildGlobalPOD() {

	char *vecFile = tInput->snapFile;
	if (!vecFile) strcpy(vecFile,"snapshotFiles.in");
	FILE *inFP = fopen(vecFile, "r");
	if (!inFP)  {
		com->fprintf(stderr, "*** Warning: No snapshots FILES in %s\n", vecFile);
		exit (-1);
	}

	int nData, _n;
	_n = fscanf(inFP, "%d",&nData);
	com->fprintf(stderr, "Building a global POD basis out of %d solution files \n",nData);

	/* // open snapshot reference solution if it exists
		 char *snapRefSolFile = tInput->snapRefSolutionFile;
		 DistSVec<double,dim> snapRefSol(domain.getNodeDistInfo());
	//if (inRSFP)  {
	com->fprintf(stderr, "Reading reference solution for snapshots in %s\n", snapRefSolFile);
	domain.readVectorFromFile(snapRefSolFile, 0, 0, snapRefSol); 
	}
	else{
	com->fprintf(stderr, "*** Warning: No snapshots FILES in %s\n", snapRefSolFile);
	exit (-1); 
	}
	*/
	char **snapFile = new char *[nData];
	for (int iData=0; iData < nData; ++iData)
		snapFile[iData] = new char[500];
	char snapFile1[500];
	int *numSnaps = new int[nData];
	int *numSkip = new int[nData];
	int sampleFreq = ioData->snapshots.sampleFreq;
	if (sampleFreq > 1)
		com->fprintf(stderr, " ... skipping %d snapshots at a time (taking every %dn snapshots for n=1,...)\n", sampleFreq-1,sampleFreq);
	double *snapWeight = new double[nData];
	int nSnap, nSkip;
	double weight;

	for (int iData = 0; iData < nData; ++iData){
		_n = fscanf(inFP, "%s %d %d %lf", snapFile1,&nSnap,&nSkip,&weight);
		strcpy(snapFile[iData],snapFile1);
		numSnaps[iData] = nSnap;
		numSkip[iData] = nSkip;
		snapWeight[iData] = weight;
		com->fprintf(stderr, " ... Reading snapshots from %s \n", snapFile[iData]);
	}

	if (ioData->output.transient.podFile[0] == 0)  {
		com->fprintf(stderr, "*** ERROR: POD Basis File not specified\n");
		exit (-1);
	}

	// open snapshot reference solution if it exists
	char *snapRefSolFile = tInput->snapRefSolutionFile;
	char **refSnapFile = new char *[nData];
	for (int iData=0; iData < nData; ++iData)
		refSnapFile[iData] = new char[500]; 
	char refSnapFile1[500];
	FILE *inRSFP = fopen(snapRefSolFile, "r");

	if (inRSFP){  
		com->fprintf(stderr, "Reading reference solution for snapshots in %s\n", snapRefSolFile);

		for (int iData = 0; iData < nData; ++iData){
			_n = fscanf(inRSFP, "%s", refSnapFile1);
			strcpy(refSnapFile[iData],refSnapFile1);
			com->fprintf(stderr, " ... Reading reference solution for snapshots from %s \n", refSnapFile[iData]);
		}
	}

	DistSVec<double,dim> snapRefSol(domain.getNodeDistInfo());
	//compute the total number of snapshots
	int nTotSnaps = 0;
	for (int iData = 0; iData < nData; ++iData) {
		for (int iSnap = 0; iSnap < numSnaps[iData]; ++iSnap) {
			if (iSnap>=numSkip[iData] && iSnap % sampleFreq == 0 ) {
				++nTotSnaps;
			}
		}
	}

	int incrementalSnapshots = ioData->snapshots.incrementalSnaps;	// change in solution should be used as a snapshot
	bool subtractInitialCondition = false;
	if (ioData->snapshots.subtractIC==SnapshotsData::SUBTRACT_IC_TRUE)
		subtractInitialCondition = true;
	if (incrementalSnapshots || subtractInitialCondition) {	// do not use last snapshot if incremental approach
		nTotSnaps -= nData;
		//for (int iData = 0 ; iData < nData; ++iData)
		//	--numSnaps[iData];
		if (incrementalSnapshots) com->fprintf(stderr, " ... using incremental snapshots (increment in state at each time step) ...\n");
	}


	if (ioData->snapshots.dataCompression.podMethod == DataCompressionData::Eig) {
		podMethod = 1;	// eigenvalue decomposition
		com->fprintf(stderr, " ... using eigenvalue decomposition to compute POD basis ...\n");
	}
	else {
		podMethod = 0;	// SVD
		com->fprintf(stderr, " ... using SVD to compute POD basis ...\n");
	}

	int maxVecStorage = ioData->snapshots.dataCompression.maxVecStorage;
	if (maxVecStorage > 0 && podMethod == 1){	// eigenvalue problem
		com->fprintf(stderr, "*** Warning: Not using limited memory POD (not implemented with eigenvalue decomposition) \n");
		maxVecStorage = 0;
	}

	bool limitedMemorySVD = (maxVecStorage == 0 ) ? false : true;
	if (limitedMemorySVD && maxVecStorage > nTotSnaps) {	// do not need to use limited memory algorithm
		limitedMemorySVD = false;
	}

	bool energyOnly = false;
	if (ioData->snapshots.dataCompression.energyOnly == DataCompressionData::ENERGY_ONLY_TRUE)
		energyOnly = true;
	bool computeSVD = true;
	if (energyOnly)
		computeSVD = false;


	int nSVD, iSVD, nStoredSnaps, nHandledSnaps;	//nSVD: # of sections for svd; nKeep: # vectors to keep from each section
	int *nKeep;
	nStoredSnaps = nTotSnaps;

	int nSnapsFirstSVD = 0;
	if (limitedMemorySVD) {
		nSVD = (int) ceil((double) nTotSnaps / (double) maxVecStorage);
		nKeep = new int [nSVD];
		for (int i = 0; i < nSVD; ++i)
			nKeep[i] = 0;
		nSnapsFirstSVD = nTotSnaps - ((nSVD -1)* maxVecStorage);
		iSVD = 0;	// start with the second SVD (want first to be smallest window)
		nHandledSnaps = 0;
		while (nHandledSnaps < maxVecStorage) {
			if ( iSVD == 0  && nKeep[0] == nSnapsFirstSVD )
				iSVD = 1;
			++nKeep[iSVD];
			++iSVD;
			++nHandledSnaps;
			if ( iSVD == nSVD )
				iSVD = 0;
		}

		nStoredSnaps = maxVecStorage;
		iSVD = 0;
	}

	nHandledSnaps = 0;

	VecSet< DistSVec<double, dim> > snap(nStoredSnaps, domain.getNodeDistInfo());
	VecSet< DistSVec<double, dim> > fullSnaps(maxVecStorage, domain.getNodeDistInfo());	// for limited memory
	DistSVec<double, dim> snapBuf(domain.getNodeDistInfo());
	double *eig = new double [nStoredSnaps];
	double eigBuf;

	int numCurrentSnapshots = 0;
	bool endOfFile;
	int normalizeSnapshots = ioData->snapshots.normalizeSnaps;
	if (normalizeSnapshots) com->fprintf(stderr, " ... snapshots normalized ...\n");
	int iStoredSnaps = 0;	// how many snapshots have been stored so far

	if (subtractInitialCondition && inRSFP)
		com->fprintf(stderr, " ... WARNING: not using reference solution file. Subtracting initial condition instead ...\n");

	for (int iData=0; iData < nData; ++iData){
		// read in Snapshot Vectors
		if (inRSFP)
			domain.readVectorFromFile(refSnapFile[iData], 0, 0, snapRefSol);
		for (int iSnap = 0; iSnap < numSnaps[iData]; ++iSnap) {
			if (!(incrementalSnapshots == 1 && sampleFreq == 1)) // always read unless incremental snapshots and not skipping at all
				endOfFile = domain.readVectorFromFile(snapFile[iData], iSnap, &eigBuf, snapBuf);
			if (subtractInitialCondition && iSnap == 0)
				snapRefSol = snapBuf;
			if (iSnap >= numSkip[iData] && iSnap % sampleFreq == 0  && !(subtractInitialCondition && iSnap == 0)) {	// use snapshot
				snap[numCurrentSnapshots] = snapBuf;
				eig[numCurrentSnapshots] = eigBuf;
				if ((inRSFP || subtractInitialCondition) && incrementalSnapshots == 0) snap[numCurrentSnapshots] -= snapRefSol;	// do not need to subtract snapRefSol if using incremental snapshots
				if (incrementalSnapshots == 1) {
					endOfFile = domain.readVectorFromFile(snapFile[iData], iSnap, &eigBuf, snapBuf);
					snap[numCurrentSnapshots] -= snapBuf;
				}
				if (normalizeSnapshots > 0) normalizeSnap(snap[numCurrentSnapshots], iSnap, numSnaps[iData]);
				if (snapWeight[iData]) snap[numCurrentSnapshots] *= snapWeight[iData]; //CBM--check
				totalEnergy += pow(snap[numCurrentSnapshots].norm(),2);	// total energy is square of frobenius norm
				++numCurrentSnapshots;
				++nHandledSnaps;
				if (limitedMemorySVD && (numCurrentSnapshots == maxVecStorage || nHandledSnaps == nTotSnaps || (iSVD == 0 && numCurrentSnapshots == nSnapsFirstSVD ))) {
					// compute SVD of snaps and put into fullSnaps matrix;
					if (computeSVD) makeFreqPOD(snap, numCurrentSnapshots, nKeep[iSVD], false);
					for (int iFullSnap = 0; iFullSnap < nKeep[iSVD] ; ++iFullSnap) {
						fullSnaps[iStoredSnaps] = snap[iFullSnap];
						++iStoredSnaps;
					}
					numCurrentSnapshots = 0;
					++iSVD;
				}
			}
		}
	}

	com->fprintf(stderr, " ... Total snapshot energy is %e \n", totalEnergy );

	if (limitedMemorySVD && computeSVD) {
		makeFreqPOD(fullSnaps, iStoredSnaps);
		delete [] nKeep;
	}
	else if (computeSVD)
		makeFreqPOD(snap, nStoredSnaps);

	delete [] numSnaps;
	delete [] snapWeight;
	delete [] eig;

}

//---------------------------------------------------------------------------------------


// template<int dim>
// void ModalSolver<dim>::projectFullSoltn() {		// this function is no longer supported
// 
//  // read in POD Basis
//  int nPodVecs = ioData->rom.dimension;
//  VecSet<DistSVec<double, dim> > podVecs(nPodVecs, domain.getNodeDistInfo());
//  readPodVecs(podVecs, nPodVecs);
// 
//  // setup solvers
//  VarFcn *varFcn = new VarFcn(*ioData); 
//  //VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData);
//  geoState = new DistGeoState(*ioData, &domain);
//  geoState->setup1(tInput->positions, &Xref, &controlVol);
//  bcData = new DistBcDataSA<dim>(*ioData, varFcn, &domain, Xref);
// 
//  spaceOp = new SpaceOperator<dim>(*ioData, varFcn, bcData, geoState, &domain);
//  postOp = new PostOperator<dim> (*ioData, varFcn, bcData, geoState, &domain);
// 
//  // Temporal operator contains Uref for us
//  tState = new DistTimeState<dim>(*ioData, spaceOp, varFcn, &domain);
//  tState->setup(tInput->solutions,  bcData->getInletBoundaryVector(), Xref, Uref);
//  geoState->setup2(tState->getData());
// 
//  RefVal *refVal = new RefVal(ioData->ref.rv);
//  tOutput = new TsOutput<dim>(*ioData, refVal, &domain, postOp);
//  tRestart = new TsRestart(*ioData, refVal);
// 
//  tOutput->openAsciiFiles();
// 
//  // initialization
//  double* Wr = new double[nPodVecs];
//  DistSVec<double, dim> WFull(domain.getNodeDistInfo());
// 
// 
//  // read in State Vectors
//  char *snapFile = 0;
//  int sp = strlen(ioData->input.prefix) + 1;
//  snapFile = new char[sp + strlen(ioData->input.stateVecFile)];
//  sprintf(snapFile, "%s%s", ioData->input.prefix, ioData->input.stateVecFile);
// 
//  VecSet< DistSVec<double, dim> > snap(1, domain.getNodeDistInfo());
//  //double *tag = new double[1];
//  double tag = 0.0;
// 
//  // project and output forces
//  int iSnap = 0;
//  int j,ierr;
//  bool endOfFile = true;
//  while(endOfFile) {
//    endOfFile = domain.readVectorFromFile(snapFile, iSnap, &tag, snap[0]);
//    if(endOfFile) { //CBM -- BAD fix
//      for (j = 0; j < nPodVecs; ++j) {
//        Wr[j] = podVecs[j] * snap[0];
//      }
//      WFull = 0.0;
//      for (j = 0; j < nPodVecs; ++j) {
//        WFull += Wr[j] * podVecs[j];
//      }
// 
//      ierr = domain.checkSolution(varFcn, WFull);
//      if (ierr)
//        com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
// 
//      tOutput->writeForcesToDisk(0, iSnap, 0, 0, tag, com->cpuNum(), tRestart->energy, Xref, WFull);
//      tOutput->writeLiftsToDisk((*ioData), 0, iSnap, 0, 0, tag, com->cpuNum(), tRestart->energy, Xref, WFull);
//      iSnap++;
//    }
//  }
// }
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
void ModalSolver<dim>::interpolatePOD()  {


	Timer *modalTimer = domain.getTimer();
	com->fprintf(stderr, " ... Interpolating POD on a tangent space to the Grassmann manifold \n");

	char *vecFile = tInput->podFile;
	if (!vecFile)
	{
		string str = "podFiles.in";
		vecFile    =  new char [str.size()+1];
		strcpy (vecFile, str.c_str());
		// vecFile now contains a c-string copy of str. Adam 2010.08.23 : g++4.3 was complaining
	}

	FILE *inFP = fopen(vecFile, "r");
	if (!inFP)  {
		com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
		exit (-1);
	}
	int nData, _n;
	_n = fscanf(inFP, "%d",&nData);

	char **podFile = new char *[nData];

	for (int iData = 0; iData < nData; ++iData){
		podFile[iData] = new char[500];
		//char *podFile1 = new char[500];
		_n = fscanf(inFP, "%s", podFile[iData]);
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
		_n = fscanf(inFP, "%lf", mach+iData);
	_n = fscanf(inFP, "%lf", &newMach);
	for (int iData = 0; iData < nData; ++iData)
		_n = fscanf(inFP, "%lf", angle+iData);
	_n = fscanf(inFP, "%lf", &newAngle);

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
				for (int k = 0; k < numPod; ++k) {
					matVals[j][ k] = podRef[j] * ((*pod[iData])[k]);
				}
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
			ParallelRom<dim> parallelRom(domain,com); 
			parallelRom.parallelSVD(*projMap[iData],U,Sigma,V,numPod);//call SVD here 

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
	ParallelRom<dim> parallelRom(domain,com);
	parallelRom.parallelSVD(logMapInterp,UInterp,SigInterp,VInterp,numPod);//call SVD here

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


	com->fprintf(stderr,"End of InterpolatePOD\n'");
	modalTimer->setRunTime();


}
//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::evalFluidSys(VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

 // Allocate ROM operators
 VarFcn *varFcn = new VarFcn(*ioData); 
 MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain);
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
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain);
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
   tmpVec = DE[iVec] * invDt;
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
   for (jVec = 0; jVec < sysSize; jVec++) {
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

template<int dim>
void ModalSolver<dim>::outputPODVectors(VecSet<DistSVec<double, dim> > &podVecs, 
                      Vec<double> &sVals, int nPOD)  {

  Timer *modalTimer = domain.getTimer();

 // write header
 int sp = strlen(ioData->output.transient.prefix) + 1;
 const char *sValExtension = ".singularVals";

 char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
 char *sValsFileName = new char[sp + strlen(ioData->output.transient.podFile)+strlen(sValExtension)];
 sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);
 if (com->cpuNum() == 0)
	 sprintf(sValsFileName, "%s%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile, sValExtension);
 FILE *sValsFile;
 if (com->cpuNum() == 0)	// only open with cpu 0
	 sValsFile = fopen(sValsFileName, "wt");

 com->fprintf(sValsFile,"%d\n", nPOD);
 com->fprintf(stderr, " ... Writing %d (%f)POD vectors to File\n", nPOD, (double) nPOD);

 com->fprintf(sValsFile,"Singular values\n");
 domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVecs[0]);

 for (int jj = 0; jj < nPOD; ++jj)
	 com->fprintf(sValsFile,"%e ", sVals[jj]);
 com->fprintf(sValsFile,"\n");
 computeRelativeEnergy(sValsFile, sVals, nPOD);

 //const int waitTime = 60;

 for (int jj = 0; jj < nPOD; ++jj) {
   com->fprintf(stderr, "%d %e\n", jj, sVals[jj]);
	 //com->fprintf(stderr, " ... waiting %d seconds ...\n", waitTime);
	 com->barrier();
	 //wait(waitTime);	// avoid file system crash due to excessive I/O
   domain.writeVectorToFile(podFileName, jj+1, sVals[jj]*sVals[jj], podVecs[jj]);
 }

}

#ifdef DO_MODAL
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::outputPODVectors(ARluSymStdEig<double> &podEigProb,
                      VecSet<DistSVec<double, dim> > &snaps, int nPOD, int numSnapsTot)  {

  Timer *modalTimer = domain.getTimer();

 DistSVec<double, dim> podVec(domain.getNodeDistInfo());
 podVec = 0.0;

 // write header
 int sp = strlen(ioData->output.transient.prefix) + 1;
 const char *sValExtension = ".singularVals";

 char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
 char *sValsFileName = new char[sp + strlen(ioData->output.transient.podFile)+strlen(sValExtension)];
 sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);
 if (com->cpuNum() == 0)
	 sprintf(sValsFileName, "%s%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile, sValExtension);
 FILE *sValsFile;
 if (com->cpuNum() == 0)	// only open with cpu 0
	 sValsFile = fopen(sValsFileName, "wt");
 delete [] sValsFileName;

 com->fprintf(sValsFile,"%d\n", nPOD);
 com->fprintf(stderr, " ... Writing %d (%f) POD vectors to File\n", nPOD, (double) nPOD);
 //const char *podFileName = ioData->output.transient.podFile;

 com->fprintf(sValsFile,"Singular values\n");
 domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVec);	// dummy vector

 double t0;
 double *rawEigVec = new double[numSnapsTot];
 int jj, kk;
 Vec<double> sVals(nPOD);
 //const int waitTime = 60;
 for (jj = 0; jj < nPOD; jj++)  {
   t0 = modalTimer->getTime();
   for (kk = 0; kk < numSnapsTot; kk++)
     rawEigVec[kk] = podEigProb.Eigenvector(nPOD-jj-1, kk);

   Vec<double> pVec(numSnapsTot, rawEigVec);
   podVec = 0.0;
   for (kk = 0; kk < numSnapsTot; kk++)
     podVec += snaps[kk] * pVec[kk];

   double eig = podEigProb.Eigenvalue(nPOD-jj-1);
   if (eig < 0.0) {
     com->fprintf(stderr,"Negative eigenvalue: %e\n",eig);
     break;
   }
   podVec *= 1.0/sqrt(eig);

   // do gram-schmidt on podVec
   //for (int jVec = 0; jVec < jj; jVec++)
     //podVec[jj] -= podVec[jVec] * (podVec[jj]*podVec[jVec]);

   //norm = podVec.norm();
   //podVec *= (1.0 / norm);
   modalTimer->addGramSchmidtTime(t0);
	 //com->barrier();
	 //com->fprintf(stderr, " ... waiting %d seconds ...\n", waitTime);
	 //wait(waitTime);	// avoid file system crash
	 domain.writeVectorToFile(podFileName, jj+1, eig, podVec);
	 com->fprintf(stderr, "%d %e\n", jj, sqrt(eig));
	 sVals[jj] = sqrt(eig);
 }
 delete [] rawEigVec;
 delete [] podFileName;

 for (int jj = 0; jj < nPOD; ++jj) {
	 com->fprintf(sValsFile,"%e ", sVals[jj]);
 }
 com->fprintf(sValsFile,"\n");
 computeRelativeEnergy(sValsFile, sVals, nPOD);

}
#endif

template<int dim>
void ModalSolver<dim>::computeRelativeEnergy(FILE *sValsFile, const Vec<double> &sVals, const int nPod){

	// TODO: optimize!

	com->fprintf(sValsFile,"Relative energy: s(i)^2/sum(s(1:end).^2)\n");
	std::vector<double> relEnergy;
	int nSnap = sVals.size();
	/*	// totalEnergy now computed when snapshots are read
	if (totalEnergy == 0.0) {
		for (int i = 0; i < nSnap; ++i)
			totalEnergy += pow(sVals[i],2);
	}
	*/

	for (int i = 0; i < nSnap; ++i) {
		double currentRelEnergy = pow(sVals[i],2)/totalEnergy;
		com->fprintf(sValsFile,"%d %e\n", i+1, currentRelEnergy);
		relEnergy.push_back(currentRelEnergy);
	}
	com->fprintf(sValsFile,"Cumulative energy: sum(s(1:k).^2)/sum(s(1:end).^2)\n");
	double cumulativeEnergy = 0.0;
	double criteria [10] = {0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999};
	int energyIndex [10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int handledCriteria  [10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int critCounter = 0;
	for (int i = 0; i < nSnap; ++i) {
		cumulativeEnergy+=relEnergy[i];
		com->fprintf(sValsFile,"%d %e\n", i+1, cumulativeEnergy);

		int critCounterTmp = 0;
		for (int j = critCounter; j < 10;++j) {
			if (cumulativeEnergy >= criteria[j] && handledCriteria[j] == 0) {
				energyIndex[j] = i;
				handledCriteria[j] = 1;
				++critCounterTmp;
			}
		}
		critCounter +=critCounterTmp;
	}
	com->fprintf(sValsFile,"Cumulative energy indices\n");

	for (int i = 0; i < 10; ++i) {
		com->fprintf(sValsFile,"%e: %d\n", criteria[i], energyIndex[i]+1);
	}
}

template<int dim>
void ModalSolver<dim>::normalizeSnap(DistSVec<double, dim> &snap, const int iSnap, const int nSnaps){

	// PURPOSE: normalize snaphots

	double scalingFactor, magnitude;

	const double smallestScaling = 0.1;	// sets smallest scaling of any snaphot
	const double tau = pow((double)nSnaps-1,2)/log(smallestScaling);	

	//com->fprintf(stderr, " ... Distance weights:\n");
	magnitude = snap.norm();
	if (magnitude == 0.0)
		magnitude = 1.0;
	scalingFactor = 1.0/magnitude;
	if (ioData->snapshots.snapshotWeights == SnapshotsData::RBF) {
		double distanceWeight = exp(pow((double)iSnap,2)/tau);
		scalingFactor *= distanceWeight;
		com->fprintf(stderr, "%d: %e\n", iSnap, distanceWeight);
	}
	snap *= scalingFactor;
}


template<int dim>
void ModalSolver<dim>::wait(const int seconds )
{
	clock_t endwait;
	endwait = clock () + seconds * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}

