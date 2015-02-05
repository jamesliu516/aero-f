#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <NewtonSolver.h>

#include <TsOutput.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::ImplicitRomTsDesc(IoData &_ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(_ioData, geoSource, dom), pod(0, dom->getNodeDistInfo()), F(dom->getNodeDistInfo()), AJ(0, dom->getNodeDistInfo()),
  rom(0) {

	ioData = &_ioData;
  tag = 0;

  maxItsNewton = ioData->ts.implicit.newton.maxIts;
  epsNewton = ioData->ts.implicit.newton.eps;  
  epsAbsResNewton = ioData->ts.implicit.newton.epsAbsRes;  
  epsAbsIncNewton = ioData->ts.implicit.newton.epsAbsInc;  

  this->timeState = new DistTimeState<dim>(*ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;

  mvpfd = new MatVecProdFD<dim,dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, *ioData);

  switch (ioData->romOnline.systemApproximation) {
    case (NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE):
      rom = new NonlinearRomOnlineII<dim>(dom->getCommunicator(), _ioData, *dom);
      systemApprox=false;
      break;
    case (NonlinearRomOnlineData::GNAT):
      rom = new NonlinearRomOnlineIII<dim>(dom->getCommunicator(), _ioData, *dom);
      systemApprox=true;
      break;
    default:
      this->com->fprintf(stderr, "*** Error:  Unexpected system approximation type\n");
      exit (-1);
  }

  currentCluster = -1;

  // nPod = 0 ?
  nPod = ioData->romOnline.maxDimension;

  // necessary? 
  pod.resize(nPod);
  AJ.resize(nPod);
  dUromNewtonIt.resize(nPod);
  dUromTimeIt.resize(nPod);
  dUromCurrentROB.resize(nPod);
  UromCurrentROB.resize(0);
  dUromNewtonIt = 0.0;
  dUromTimeIt = 0.0;
  dUromCurrentROB = 0.0;
 
  //dUnormAccum = new Vec<double> [2];
  //for (int i = 0 ; i < 2; ++i) { 
  //  dUnormAccum[i].resize(nPod);	// before any time step, it is zero
  //  dUnormAccum[i] = 0.0;	// before any time step, it is zero
  //}

  MemoryPool mp;
  this->mmh = this->createMeshMotionHandler(*ioData, geoSource, &mp);
 
  basisUpdateFreq = ioData->romOnline.basisUpdateFreq;
  tryAllFreq = ioData->romOnline.tryAllFreq;

  updateFreq = false;
  clusterSwitch = false;
  updatePerformed = false;

  if (ioData->romOnline.weightedLeastSquares!=NonlinearRomOnlineData::WEIGHTED_LS_FALSE) {
    weightVec = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  } else {
    weightVec = NULL;
  }

/*  if (ioData->romOnline.weightedLeastSquares==NonlinearRomOnlineData::WEIGHTED_LS_STATE) {
    weightURef = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  } else {
    weightURef = NULL;
  }

  if (ioData->romOnline.weightedLeastSquares==NonlinearRomOnlineData::WEIGHTED_LS_RESIDUAL) {
    weightFRef = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  } else {
    weightFRef = NULL;
  }*/

  interiorWeight = 1.0;
  ffWeight = this->ioData->romOnline.ffWeight;
  wallWeight = this->ioData->romOnline.wallWeight;
  bcWeightGrowthFactor = this->ioData->romOnline.bcWeightGrowthFactor;
  levenbergMarquardtWeight = this->ioData->romOnline.levenbergMarquardtWeight;
  wallUp = false;
  wallDown = false;
  wallWeightGrowthFactor = bcWeightGrowthFactor;
  ffUp = false;
  ffDown = false;
  ffWeightGrowthFactor = bcWeightGrowthFactor;

 if (ioData->romOnline.weightedLeastSquares==NonlinearRomOnlineData::WEIGHTED_LS_BOCOS) {
    farFieldMask = new DistVec<double>(this->domain->getNodeDistInfo());
    farFieldNeighborsMask = new DistVec<double>(this->domain->getNodeDistInfo());
    this->domain->setFarFieldMask(*farFieldMask, *farFieldNeighborsMask);
    wallMask = new DistVec<double>(this->domain->getNodeDistInfo());
    wallNeighborsMask = new DistVec<double>(this->domain->getNodeDistInfo());
    this->domain->setWallMask(*wallMask, *wallNeighborsMask);
    updateLeastSquaresWeightingVector();
  } else {
    farFieldMask = NULL;
    farFieldNeighborsMask = NULL;
    wallMask = NULL;
    wallNeighborsMask = NULL;
  }

  allowBCWeightDecrease = this->ioData->romOnline.allowBCWeightDecrease;
  adjustInteriorWeight = this->ioData->romOnline.adjustInteriorWeight;

  regThresh = this->ioData->romOnline.regThresh;
  regWeight = 0.0;

  Uinit = NULL;

  useIncrements = (ioData->romDatabase.avgIncrementalStates==NonlinearRomFileSystemData::AVG_INCREMENTAL_STATES_TRUE) ? true : false;

  if ((!systemApprox) && useIncrements) {
    Uprev = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  } else {
    Uprev = NULL;
  }

  unsteady = this->problemType[ProblemData::UNSTEADY];

  tryingAllClusters = false;

  rho = 0.5;
  c1 = 0.25;
  maxItsLS = 10;

  const char *residualsFileName = ioData->output.rom.residualsForCoordRange;
  if (strcmp(residualsFileName, "") != 0)  {
    char *residualsPath = new char[strlen(ioData->output.rom.prefix) + strlen(residualsFileName)+1];
    if (this->com->cpuNum() == 0) sprintf(residualsPath, "%s%s",ioData->output.rom.prefix, residualsFileName);
    if (this->com->cpuNum() == 0) residualsFile = fopen(residualsPath, "wt");
    delete [] residualsPath;
  } else { 
    residualsFile = NULL;
  }


}

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::~ImplicitRomTsDesc()
{

	// output net dUnormAccum
/*	for (int iPod = 0; iPod < nPod; ++iPod)
		dUnormAccum[1][iPod] = sqrt(dUnormAccum[1][iPod]);	// to complete 2-norm definition

	const char *dUnormAccumFile = ioData->output.rom.dUnormAccum;
	if (strcmp(dUnormAccumFile, "") != 0)  {
		char *outdUnormAccumFile = new char[strlen(ioData->output.rom.prefix) +
			strlen(dUnormAccumFile)+1];
		if (this->com->cpuNum() ==0) sprintf(outdUnormAccumFile, "%s%s",
				ioData->output.rom.prefix, dUnormAccumFile);
		FILE *writingFile;
		if (this->com->cpuNum() ==0) writingFile = fopen(outdUnormAccumFile, "wt");
		this->com->fprintf(writingFile, "PodCoefficient NetContribution(1norm) NetContribution(2norm)\n");
		for (int i = 0; i < nPod; ++i) { 
			this->com->fprintf(writingFile, "%d %e %e \n", i + 1, dUnormAccum[0][i],dUnormAccum[1][i]);
		}
		delete [] outdUnormAccumFile;
	}*/

	//delete [] dUnormAccum;

  if (Uinit && residualsFile) printRomResiduals(*Uinit);
  if (this->com->cpuNum()==0 && residualsFile) fclose(residualsFile);

  if (tag) delete tag;
  if (weightVec) delete weightVec;
  //if (weightURef) delete weightURef;
  //if (weightFRef) delete weightFRef;
  if (farFieldMask) delete farFieldMask;
  if (farFieldNeighborsMask) delete farFieldNeighborsMask;
  if (wallMask) delete wallMask;
  if (wallNeighborsMask) delete wallNeighborsMask;
  if (Uinit) delete Uinit;
  if (Uprev) delete Uprev;

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitRomTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
{

  KspPrec<neq> *_pc = 0;

  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<Scalar,neq>(DiagMat<Scalar,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS ||
     pcdata.type == PcData::RAS ||
     pcdata.type == PcData::ASH ||
     pcdata.type == PcData::AAS)
    _pc = new IluPrec<Scalar,neq>(pcdata, dom);
//  else if (pcdata.type == PcData::MG)
//    _pc = new MultiGridPrec<Scalar,neq>(dom, *this->geoState);

  return _pc;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::printRomResiduals(DistSVec<double, dim> &U)  {
  // calculates the residual corresponding for the initial condition of the simulation plus an
  // increment in the first basis vector, for a range of increments defined by min (minimum generalized coord value),
  // max (maximum generalized coord value), and res (number of values between min and max).  

  double min = ioData->romOnline.residualsCoordMin; 
  double max = ioData->romOnline.residualsCoordMax;
  double res = ioData->romOnline.residualsCoordRes;

  if ((max-min)<=0.0) return;

  DistSVec<double, dim>* Ueval = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  DistSVec<double, dim>* Reval = new DistSVec<double, dim>(this->domain->getNodeDistInfo());

  double delta = (max - min) / res;

  for (int i=0; i<=res; ++i) {
   double coef = min + i*delta;
   *Ueval = U + (*rom->basis)[0]*coef;
   int dummyIt = 0; // assume spatial only
   computeFullResidual(dummyIt, *Ueval, false, Reval);
   double unweightedResNormSquared = pow(Reval->norm(),2);
   computeFullResidual(dummyIt, *Ueval, true, Reval);
   double weightedResNormSquared = pow(Reval->norm(),2);
   if (this->com->cpuNum() ==0)
     this->com->fprintf(residualsFile, "%1.15e %1.15e %1.15e\n", coef, unweightedResNormSquared, weightedResNormSquared); 
  }

  delete Ueval;
  delete Reval;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::tryAllClusters(DistSVec<double, dim> &U, const int totalTimeSteps, int* bestCluster) {

  this->com->fprintf(stdout, " ... trying to reduce residual with all clusters ...\n");

  double bestResReduction = -1.0;
  
  DistSVec<double, dim> Ubackup(this->domain->getNodeDistInfo());
  Ubackup = U;

  int dummyInt = 0;
  bool dummyTrue = true;
  bool dummyFalse = false;

  tryingAllClusters=true;

  for (int iCluster=0; iCluster<(rom->nClusters); ++iCluster) {

    loadCluster(iCluster, dummyTrue, U);

    computeFullResidual(dummyInt, U, true);
    double resInit = (this->F)*(this->F);

    solveNonLinearSystem(U, totalTimeSteps);

    computeFullResidual(dummyInt, U, true);
    double resFinal = (this->F)*(this->F);

    if (resFinal==0.0) {
      *bestCluster = iCluster;
      return;
    }

    double resReduction = resFinal / resInit;

    this->com->fprintf(stdout, " ... reduced residual from %e to %e (ratio of %e) using cluster %d\n",
                              resInit, resFinal, resReduction, iCluster);
  
    if ((resReduction<bestResReduction) || (bestResReduction<0.0)) {
      *bestCluster = iCluster;
      bestResReduction = resReduction;
    }

    U = Ubackup;

  }

  tryingAllClusters=false;
  currentCluster = -1;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::checkLocalRomStatus(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  // checks whether the local ROM needs to be modified

  if ((rom->nClusters > 1) || (basisUpdateFreq>0) || (currentCluster == -1)) {

    if (currentCluster == -1) { // first iteration
      if (ioData->romOnline.distanceComparisons)
        rom->initializeDistanceComparisons(U);
      if (ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_FAST_EXACT) {
        rom->initializeFastExactUpdatesQuantities(U);
      } else if (ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
        rom->initializeProjectionQuantities(U);
      }
    }

    int closestCluster;

    if (rom->nClusters > 1) {
      if (useIncrements) {
        bool tryAllNow = ((tryAllFreq > 0) && (totalTimeSteps%tryAllFreq == 0)) ? true : false;
        if ((currentCluster==-1) || (tryAllNow) || (dUromTimeIt.norm()<ioData->romOnline.incrementCoordsTol)) { 
          tryAllClusters(U, totalTimeSteps, &closestCluster);
        } else {
          if (systemApprox) {
            rom->closestCenter(U, &closestCluster); // U isn't actually used (and Uprev is NULL)
          } else {
            DistSVec<double, dim> Uincrement(this->domain->getNodeDistInfo());
            Uincrement = U - *Uprev; 
            rom->closestCenter(Uincrement, &closestCluster);
          }
        }
      } else {
        rom->closestCenter(U, &closestCluster);
      }
      this->com->fprintf(stdout, " ... using cluster %d\n", closestCluster);
    } else {
      closestCluster = 0;
    }

    updateFreq = ((basisUpdateFreq > 0) && (totalTimeSteps%basisUpdateFreq == 0)) ? true : false;
    clusterSwitch = (currentCluster != closestCluster) ? true : false;
    updatePerformed = false;

    if (updateFreq || clusterSwitch) loadCluster(closestCluster, clusterSwitch, U);

  }
  dUromTimeIt = 0.0;

}

//-----------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::loadCluster(int closestCluster, bool clusterSwitch, DistSVec<double, dim> &U) {

  if (clusterSwitch) {
    deleteRestrictedQuantities(); // only defined for GNAT
    int prevCluster = currentCluster;
    currentCluster = closestCluster;
    rom->readClusteredOnlineQuantities(currentCluster);  // read state basis, update info, and (if applicable) gappy matrices
    if (this->ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
      if (this->ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_OFF) {
        rom->projectSwitchStateOntoAffineSubspace(currentCluster, prevCluster, U, UromCurrentROB);
      } else {
        this->com->fprintf(stderr, "*** Warning: Updates and projection were both specified; not performing projection\n");
      }
    }
  }

  if (this->ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
    updatePerformed = rom->updateBasis(currentCluster, U, &dUromCurrentROB);
    if (this->ioData->romOnline.bufferEnergy!=0.0) rom->truncateBufferedBasis();
    if (this->ioData->romOnline.krylov.include) rom->appendNonStateDataToBasis(currentCluster,"krylov");
    if (this->ioData->romOnline.sensitivity.include) rom->appendNonStateDataToBasis(currentCluster,"sensitivity");
  }
  
  nPod = rom->basis->numVectors();
  pod.resize(nPod);
  for (int iVec=0; iVec<nPod; ++iVec) {
    pod[iVec] = (*(rom->basis))[iVec];
  }
  AJ.resize(nPod);
  dUromNewtonIt.resize(nPod);
  dUromTimeIt.resize(nPod);
  dUromCurrentROB.resize(nPod);
  dUromCurrentROB = 0.0;
  setProblemSize(U);  // defined in derived classes
  // TODO also set new reference residual if the weighting changes
  if (clusterSwitch && !unsteady) setReferenceResidual(); // for steady gnat (reference residual is restricted to currently active nodes)

}

//-----------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::printBCWeightingInfo(bool updateWeights) {

  int numFFNodes = 0;
  double resNormSquaredFFNodes = 0.0;

  int numFFNeighborNodes = 0;
  double resNormSquaredFFNeighborNodes = 0.0;

  int numWallNodes = 0;
  double resNormSquaredWallNodes = 0.0;

  int numWallNeighborNodes = 0;
  double resNormSquaredWallNeighborNodes = 0.0;

  int numInteriorNodes = 0;
  double resNormSquaredInteriorNodes = 0.0;

  // weight residual
  int numLocSub = this->F.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    bool *masterFlag = this->domain->getNodeDistInfo().getMasterFlag(iSub);
    double *ffMask = farFieldMask->subData(iSub); // vector with nonzero entries at farfield nodes
    double *ffNeighborsMask =  farFieldNeighborsMask->subData(iSub); // vector with nonzero entries at neighbors of farfield nodes
    double *wMask = wallMask->subData(iSub);   // vector with nonzero entries at wall nodes
    double *wNeighborsMask =  wallNeighborsMask->subData(iSub); // vector with nonzero entries at neighbors of wall nodes
    double (*unweightedResidual)[dim] = this->F.subData(iSub);
    for (int i=0; i<F.subSize(iSub); ++i) {
      if (!masterFlag[i]) continue;
      if (ffMask[i]>0) {
        ++numFFNodes;
        for (int j=0; j<dim; ++j)
          resNormSquaredFFNodes += pow(unweightedResidual[i][j],2);
      } else if (ffNeighborsMask[i]>0) {
        ++numFFNeighborNodes;
        for (int j=0; j<dim; ++j)
          resNormSquaredFFNeighborNodes += pow(unweightedResidual[i][j],2);
      } else if (wMask[i]>0) {
        ++numWallNodes;
        for (int j=0; j<dim; ++j)
          resNormSquaredWallNodes += pow(unweightedResidual[i][j],2);
      } else if (wNeighborsMask[i]>0) {
        ++numWallNeighborNodes;
        for (int j=0; j<dim; ++j)
          resNormSquaredWallNeighborNodes += pow(unweightedResidual[i][j],2);
      } else {
        ++numInteriorNodes;
        for (int j=0; j<dim; ++j)
          resNormSquaredInteriorNodes += pow(unweightedResidual[i][j],2);
      }
    }
  }

  this->com->globalSum(1, &numFFNodes);
  this->com->globalSum(1, &numFFNeighborNodes);
  this->com->globalSum(1, &numWallNodes);
  this->com->globalSum(1, &numWallNeighborNodes);
  this->com->globalSum(1, &numInteriorNodes);
  this->com->globalSum(1, &resNormSquaredFFNodes);
  this->com->globalSum(1, &resNormSquaredFFNeighborNodes);
  this->com->globalSum(1, &resNormSquaredWallNodes);
  this->com->globalSum(1, &resNormSquaredWallNeighborNodes);
  this->com->globalSum(1, &resNormSquaredInteriorNodes);

  if (numFFNodes==0) numFFNodes=1;
  if (numFFNeighborNodes==0) numFFNeighborNodes=1;
  if (numWallNodes==0) numWallNodes=1;
  if (numWallNeighborNodes==0) numWallNeighborNodes=1; 

  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "(||R_interior||)^2   = %e\n", resNormSquaredInteriorNodes);
  this->com->fprintf(stdout, "epsilon_interior = %e\n", interiorWeight);
  this->com->fprintf(stdout, "(||W_interior*R_interior||)^2 = %e\n", resNormSquaredInteriorNodes*interiorWeight);
  this->com->fprintf(stdout, "(||R_interior||)^2 / numInteriorNodes = %e\n", resNormSquaredInteriorNodes/((double) numInteriorNodes));
  this->com->fprintf(stdout, "(||W_interior*R_interior||)^2 / numInteriorNodes = %e\n",  resNormSquaredInteriorNodes*interiorWeight/ ((double) numInteriorNodes));
  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "(||R_ff||)^2 = %e\n", resNormSquaredFFNodes);
  this->com->fprintf(stdout, "epsilon_ff = %e\n", ffWeight);
  this->com->fprintf(stdout, "(||W_ff*R_ff||)^2 = %e\n", resNormSquaredFFNodes*ffWeight);
  this->com->fprintf(stdout, "(||R_ff||)^2 / numFarFieldNodes = %e\n", resNormSquaredFFNodes / ((double) numFFNodes));
  this->com->fprintf(stdout, "(||W_ff*R_ff||)^2 / numFarFieldNodes = %e\n",  resNormSquaredFFNodes*ffWeight/ ((double) numFFNodes));
  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "(||R_ff_neighbors||)^2 = %e\n", resNormSquaredFFNeighborNodes);
  this->com->fprintf(stdout, "(||R_ff_neighbors||)^2 / numFarFieldNeighborNodes = %e\n", resNormSquaredFFNeighborNodes / ((double) numFFNeighborNodes));
  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "(||R_wall||)^2 = %e\n", resNormSquaredWallNodes);
  this->com->fprintf(stdout, "epsilon_wall = %e\n", wallWeight);
  this->com->fprintf(stdout, "(||W_wall*R_wall||)^2 = %e\n", resNormSquaredWallNodes*wallWeight);
  this->com->fprintf(stdout, "(||R_wall||)^2 / numWallNodes = %e\n", resNormSquaredWallNodes / ((double) numWallNodes));
  this->com->fprintf(stdout, "(||W_wall*R_wall||)^2 / numWallNodes = %e\n",  resNormSquaredWallNodes*wallWeight/ ((double) numWallNodes));
  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "(||R_wall_neighbors||)^2 = %e\n", resNormSquaredWallNeighborNodes);
  this->com->fprintf(stdout, "(||R_wall_neighbors||)^2 / numWallNeighborNodes = %e\n", resNormSquaredWallNeighborNodes / ((double) numWallNeighborNodes));
  this->com->fprintf(stdout, "-----------------------------------------------\n");
  this->com->fprintf(stdout, "Ninterior=%d, Nff=%d, Nwall=%d, Ntot=%d\n", numInteriorNodes, numFFNodes, numWallNodes, numInteriorNodes+numFFNodes+numWallNodes);
  this->com->fprintf(stdout, "TotWeight=%e\n", interiorWeight*(double)numInteriorNodes+(double)numWallNodes*wallWeight+(double)numFFNodes*ffWeight);

  if (updateWeights) {

    double relErrFF = (resNormSquaredFFNodes/((double)numFFNodes) - resNormSquaredFFNeighborNodes/((double)numFFNeighborNodes))/
                      (resNormSquaredFFNeighborNodes/((double)numFFNeighborNodes));

    double relErrWall = (resNormSquaredWallNodes/((double)numWallNodes) - resNormSquaredWallNeighborNodes/((double)numWallNeighborNodes))/
                        (resNormSquaredWallNeighborNodes/((double)numWallNeighborNodes));

    relErrFF = abs(relErrFF);
    relErrWall = abs(relErrWall);

    double totRelError = relErrFF+relErrWall;

    double ffRatio = 10.0;
    double wallRatio = 10.0;

    ffWeightGrowthFactor = 1.0 + (bcWeightGrowthFactor-1.0)*relErrFF/totRelError;
    wallWeightGrowthFactor = 1.0 + (bcWeightGrowthFactor-1.0)*relErrWall/totRelError;

    if ((resNormSquaredFFNodes/((double)numFFNodes))/(resNormSquaredFFNeighborNodes/((double)numFFNeighborNodes))>ffRatio && interiorWeight>0.0 
         && resNormSquaredFFNodes/((double)numFFNodes)>1e-12) {
      if (ffDown) ffWeightGrowthFactor = (ffWeightGrowthFactor-1.0)/2 + 1.0;
      ffUp = true;
      ffDown = false;
      ffWeight *= ffWeightGrowthFactor;
    } else if ((resNormSquaredFFNodes/((double)numFFNodes))/(resNormSquaredFFNeighborNodes/((double)numFFNeighborNodes))<=ffRatio && ffWeight>1.0
               && allowBCWeightDecrease) {
      if (ffUp) ffWeightGrowthFactor = (ffWeightGrowthFactor-1.0)/2 + 1.0;
      ffUp = false;
      ffDown = true;
      ffWeight /= ffWeightGrowthFactor;
    }

    if ((resNormSquaredWallNodes/((double)numWallNodes))/(resNormSquaredWallNeighborNodes/((double)numWallNeighborNodes))>wallRatio && interiorWeight>0.0
        && resNormSquaredWallNodes/((double)numWallNodes)>1e-12) {
      if (wallDown) wallWeightGrowthFactor = (wallWeightGrowthFactor-1.0)/2 + 1.0;
      wallUp = true;
      wallDown = false;
      wallWeight *= wallWeightGrowthFactor; 
    } else if ((resNormSquaredWallNodes/((double)numWallNodes))/(resNormSquaredWallNeighborNodes/((double)numWallNeighborNodes))<=wallRatio && wallWeight>1.0 
                && allowBCWeightDecrease) {
      if (wallUp) wallWeightGrowthFactor = (wallWeightGrowthFactor-1.0)/2 + 1.0;
      wallUp = false;
      wallDown = true;
      wallWeight /= wallWeightGrowthFactor;
    }

    ffWeight = (ffWeight<1.0) ? 1.0 : ffWeight;
    wallWeight = (wallWeight<1.0) ? 1.0 : wallWeight;
    if (adjustInteriorWeight) {
      interiorWeight = ((double)numInteriorNodes+(double)numWallNodes*(1.0-wallWeight)+(double)numFFNodes*(1.0-ffWeight))/((double)numInteriorNodes);
      interiorWeight = (interiorWeight<0) ? 0.0 : interiorWeight;
    }
  }

  /*  double relErrFF = (resNormSquaredFFNodes/((double)numFFNodes) - resNormSquaredInteriorNodes/((double)numInteriorNodes))/
                      (resNormSquaredInteriorNodes/((double)numInteriorNodes));

    double relErrWall = (resNormSquaredWallNodes/((double)numWallNodes)  - resNormSquaredInteriorNodes/((double)numInteriorNodes))/
                      (resNormSquaredInteriorNodes/((double)numInteriorNodes));

    relErrFF = abs(relErrFF);
    relErrWall = abs(relErrWall);

    double totRelError = relErrFF+relErrWall;

    ffWeightGrowthFactor = 1.0 + (bcWeightGrowthFactor-1.0)*relErrFF/totRelError;
    wallWeightGrowthFactor = 1.0 + (bcWeightGrowthFactor-1.0)*relErrWall/totRelError;
 
    double ffRatio = 1.0;
    double wallRatio = 10.0;

    if ((resNormSquaredFFNodes/((double)numFFNodes))/(resNormSquaredInteriorNodes/((double)numInteriorNodes))>ffRatio && interiorWeight>0.0
         && resNormSquaredFFNodes/((double)numFFNodes)>1e-12) {
      if (ffDown) ffWeightGrowthFactor = (ffWeightGrowthFactor-1.0)/2 + 1.0;
      ffUp = true;
      ffDown = false;
      ffWeight *= ffWeightGrowthFactor;
    } else if ((resNormSquaredFFNodes/((double)numFFNodes))/(resNormSquaredInteriorNodes/((double)numInteriorNodes))<=ffRatio && ffWeight>1.0
               && allowBCWeightDecrease) {
      if (ffUp) ffWeightGrowthFactor = (ffWeightGrowthFactor-1.0)/2 + 1.0;
      ffUp = false;
      ffDown = true;
      ffWeight /= ffWeightGrowthFactor;
    }

    if ((resNormSquaredWallNodes/((double)numWallNodes))/(resNormSquaredInteriorNodes/((double)numInteriorNodes))>wallRatio && interiorWeight>0.0
        && resNormSquaredWallNodes/((double)numWallNodes)>1e-12) {
      if (wallDown) wallWeightGrowthFactor = (wallWeightGrowthFactor-1.0)/2 + 1.0;
      wallUp = true;
      wallDown = false;
      wallWeight *= wallWeightGrowthFactor;
    } else if ((resNormSquaredWallNodes/((double)numWallNodes))/(resNormSquaredInteriorNodes/((double)numInteriorNodes))<=wallRatio && wallWeight>1.0
                && allowBCWeightDecrease) {
      if (wallUp) wallWeightGrowthFactor = (wallWeightGrowthFactor-1.0)/2 + 1.0;
      wallUp = false;
      wallDown = true;
      wallWeight /= wallWeightGrowthFactor;
    }

    ffWeight = (ffWeight<1.0) ? 1.0 : ffWeight;
    wallWeight = (wallWeight<1.0) ? 1.0 : wallWeight;
    if (adjustInteriorWeight) {
      interiorWeight = 1.0 - (ffWeight-1.0)*ffRatio*(double)numFFNodes/(double)numInteriorNodes - (wallWeight-1.0)*wallRatio*(double)numWallNodes/(double)numInteriorNodes;
      interiorWeight = (interiorWeight<0.0) ? 0.0 : interiorWeight;
    }
  }*/
}

//------------------------------------------------------------------------------
template<int dim>
int ImplicitRomTsDesc<dim>::solveNonLinearSystem(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  if (Uprev) *Uprev = U;

  // initializations 

  double t0 = this->timer->getTime();

  int it = 0;
  int fsIt = 0;

  DistSVec<double, dim> dUfull(this->domain->getNodeDistInfo());	// solution increment at EACH NEWTON ITERATION in full coordinates
  dUfull = 0.0;	// initial zero increment

  double res;
  bool breakloop = false;
  bool breakloopNow = false;

	// line search variables
  double alpha;
  bool convergeFlag=0;

  postProStep(U,totalTimeSteps);

  updateLeastSquaresWeightingVector(); //only updated at the start of Newton

  //if (totalTimeSteps==1) printRomResiduals(U);

  for (it = 0; it < maxItsNewton; it++)  {

   // TODO remove all of the U / R weighting machinery
   // if (it==0) {
   //   if (ioData->romOnline.weightedLeastSquares==NonlinearRomOnlineData::WEIGHTED_LS_STATE)
   //     *weightURef = U;
   //
   //   if (ioData->romOnline.weightedLeastSquares==NonlinearRomOnlineData::WEIGHTED_LS_RESIDUAL)
   //     *weightFRef = this->F;
   // }

    double tRes = this->timer->getTime();
    computeFullResidual(it, U, false);
    this->timer->addResidualTime(tRes);

    if (this->ioData->romOnline.onlineResiduals.include 
        && (it%(this->ioData->romOnline.onlineResiduals.newtonFreq)==0)
        && (totalTimeSteps%(this->ioData->romOnline.onlineResiduals.timeFreq)==0)) {
      //TODO: spatial only!
      // we want the actual flux, so inv(control volumes) * spatial residual

      if (it==0) {
        rom->appendVectorToBasis(this->F, nPod);
      } else {
        rom->appendVectorToBasis(this->F);
      }
      setProblemSize(U);
    }

    computeAJ(it, U, true);	// skipped some times for Broyden

    if (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE) {
      tRes = this->timer->getTime();
      computeFullResidual(it, U, true);
      this->timer->addResidualTime(tRes);
    }

    solveNewtonSystem(it, res, breakloop, U, totalTimeSteps);	// 1) check if residual small enough, 2) solve 
      // INPUTS: AJ, F
      // OUTPUTS: dUromNewtonIt, res, breakloop

    breakloopNow = breakloop1(breakloop);
    if (breakloopNow) break;

    if (this->ioData->romOnline.lineSearch==NonlinearRomOnlineData::LINE_SEARCH_WOLF) { 
      // do line search (linesearch exits with alpha=0 and convergenceFlag if convergence criteria is satisfied)
      alpha = lineSearch(U,dUromNewtonIt,it,AJ,epsNewton, convergeFlag);
      if (it > 0 && convergeFlag == 1) break;
      dUromNewtonIt *= alpha;
    } else if (this->ioData->romOnline.lineSearch==NonlinearRomOnlineData::LINE_SEARCH_BACKTRACKING) {
      expandVector(dUromNewtonIt, dUfull);
      for (int itLS=0; itLS<maxItsLS; ++itLS) {
        if (itLS>0){
          alpha *= rho;
          if (itLS==1) {
            dUfull *= (rho-1);
            dUromNewtonIt *= (rho-1);
          } else {
            dUfull *= rho;
            dUromNewtonIt *= rho;
          }
        } else {
          alpha = 1.0;
        }

        // increment or backtrack from previous trial 
        U += dUfull;          
        double restrial = 1.0e20;
        if (!this->checkSolution(U)) {
          computeFullResidual(it, U, true);
          restrial = (this->F)*(this->F);
        }

        if (restrial>=0.0) {
          if (sqrt(restrial) < sqrt(1-2.0*alpha*c1)*res)// || dQ.norm() <= epsAbsInc)
            break;
        }
        if (itLS == maxItsLS-1 && maxItsLS != 1)
          this->com->printf(1, "*** Warning: Line Search reached %d its ***\n", maxItsLS);
      }
    }

    double tSol = this->timer->getTime();    
    dUromTimeIt += dUromNewtonIt; // solution increment in reduced coordinates (initialized to zero in checkLocalRomStatus)
    if (this->ioData->romOnline.lineSearch!=NonlinearRomOnlineData::LINE_SEARCH_BACKTRACKING) {
      expandVector(dUromNewtonIt, dUfull); // solution increment in full coordinates
      U += dUfull;
    }
    this->timer->addSolutionIncrementTime(tSol);

    saveNewtonSystemVectors(totalTimeSteps);	// only implemeted for PG rom

    // verify that the solution is physical
    if (this->checkSolution(U)) {
      if (checkFailSafe(U) && fsIt < 5) {
        this->com->fprintf(stderr, "*** Warning: Not yet implemented\n");
        //fprintf(stderr,"*** Warning: Newton solver redoing iteration %d\n", it+1);
        //Q = rhs;
        //--it;
        //++fsIt;
      }
      else{
        this->com->fprintf(stderr, "*** Exiting\n");
        exit(-1);
      }
    }
    breakloopNow = breakloop2(breakloop);
    if (breakloopNow) break;
  }	// end Newton loop

  if (this->ioData->romOnline.weightedLeastSquares == NonlinearRomOnlineData::WEIGHTED_LS_BOCOS) {
    computeFullResidual(it, U, false);
    printBCWeightingInfo(true);
  }
	
  //savedUnormAccum();
  if (fsIt > 0 && checkFailSafe(U) == 1)
  resetFixesTag();

  if (it == maxItsNewton && maxItsNewton!=1 && maxItsNewton!=0) {
    this->com->fprintf(stderr, "*** Warning: ROM Newton solver reached %d its", maxItsNewton);
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  if (it==0) {
    this->com->fprintf(stderr, "*** Warning: ROM converged on first iteration");
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  this->timer->addFluidSolutionTime(t0);

  this->com->fprintf(stdout, "\n");

  // If trying all clusters, exit before writing reduced coordinates
  if (tryingAllClusters) return (maxItsNewton == 0 || it==0) ? 1 : it;

  // output POD coordinates
  dUromCurrentROB += dUromTimeIt;
  if (UromCurrentROB.size() == dUromTimeIt.size()) UromCurrentROB += dUromTimeIt;

  rom->writeReducedCoords(totalTimeSteps, clusterSwitch, updatePerformed, currentCluster, dUromTimeIt); 

  if (ioData->romOnline.distanceComparisons)
    rom->advanceDistanceComparisons(currentCluster, dUromTimeIt, UromCurrentROB);

  return (maxItsNewton == 0 || it==0) ? 1 : it;

}
                                                                                                           
//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::solveLinearSystem(int it , Vec<double> &rhs, Vec<double> &sol)
{

  double *x = rhs.data();
  FullM myjac(jac);
  myjac.factor();
  myjac.reSolve(x);
  sol = rhs;  

  return 0;

}

//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitRomTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy)
{
  if (R==NULL) R = &F;  // make R an alias for F

  this->spaceOp->computeResidual(*this->X, *this->A, Q, *R, this->timeState);

  if (includeHomotopy)
      this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, *R);

  this->spaceOp->applyBCsToResidual(Q, *R);  // wall BCs only

  if (applyWeighting && (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE)) {
    // weight residual
    double weightExp = this->ioData->romOnline.weightingExponent;
    double weightNorm = weightVec->norm();
    weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;
    int numLocSub = R->numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double (*weight)[dim] = weightVec->subData(iSub);
      double (*r)[dim] = R->subData(iSub);
      for (int i=0; i<weightVec->subSize(iSub); ++i) {
        for (int j=0; j<dim; ++j) {
          //weight[i][j] = pow(abs(weight[i][j])/weightNorm, weightExp);
          r[i][j] = r[i][j] * weight[i][j];
        }
      }
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &R, double stepLength)  {
	// merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  computeFullResidual(it,newQ,true,&R);

  double merit = 0.0;
  merit += R.norm();	// merit function = 1/2 * (norm of full-order residual)^2
  merit *= merit;
//  merit *= 0.5;

  if (this->ioData->romOnline.lsSolver == 2) {
    DistSVec<double, dim>* A_Uerr = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *A_Uerr = 0.0;
    int numLocSub = A_Uerr->numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = this->A->subData(iSub); // vector of control volumes
      double (*auerr)[dim] = A_Uerr->subData(iSub);
      double (*u)[dim] = newQ.subData(iSub);
      for (int i=0; i<this->A->subSize(iSub); ++i) {
        if (cv[i]>regThresh) {
          for (int j=0; j<dim; ++j)
            auerr[i][j] = cv[i] * pow(u[i][j] - (this->bcData->getInletConservativeState())[j],2);
        }
      }
    }

    DistSVec<double, dim>* ones = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *ones = 1.0;

    double regTerm = (*A_Uerr)*(*ones);
    regTerm *= regWeight;
    merit += regTerm;

    delete ones;
    delete A_Uerr;
  }

  return merit;

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunctionDeriv(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &p, DistSVec<double, dim> &R, double currentMerit)  { 

  double eps = mvpfd->computeEpsilon(Q,p);
  DistSVec<double, dim> newR(this->domain->getNodeDistInfo());

  double newMerit = meritFunction(it, Q, p, newR, eps);

  double meritDeriv = (newMerit - currentMerit) / eps;

  return meritDeriv;

/*
  double eps = mvpfd->computeEpsilon(Q,p);
  DistSVec<double, dim> newF(this->domain->getNodeDistInfo());

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + eps*p;
  computeFullResidual(it,newQ,&newF);

  newF -= F;  // overwrite new Flux with finite difference
  newF *= (1.0/eps);

  double meritDeriv = 0.0;
  meritDeriv = newF*F; // Take inner product
    // merit function = norm of full-order residual
  return meritDeriv;
*/
}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::lineSearch(DistSVec<double, dim> &Q, Vec<double> &prom, int it, VecSet<DistSVec<double, dim> > &leftProj, double eps, bool &convergeFlag){
  // This function (along with "zoom"), finds a steplength satisfying the Strong Wolfe conditions in two steps:
  // 1) Find a valid bracket containing a point that satisfies these conditions
  // 2) Perform interpolation/bisection within the bracket to locate this valid point
  // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

  // Parameters for Wolfe conditions

  double c1 = 0.00001;  // for Armijo (sufficient decrease) condition
  double c2 = 0.0001; // relative gradient condition.
  //NOTE: We require 0<c1<c2<1 
  //NOTE: the smaller the value of c2, the closer the steplength will be to a local minimizer

  // Parameters for bracketing step (step 1)

  int maxIter = 100;  // max iterations for finding a bracket
  DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo()); // new state vector initialization
  double beta = 1.2;  // amount to increment alpha each time
  double alpha_max = 100000000.0;    // bound on step size alpha  (this was about 50)
  int maxedFlag = 0;  // flag to determine whether or not alpha_max has been reached
  int count = 0;  // counter for number of iterations
  double alpha = 1.0; // initial alpha is just the Newton step
  double merit, meritZero, meritOld;
  double meritDeriv, meritDerivZero,targetMeritDerivZero, meritDerivOld;

  DistSVec<double,dim> dQ(this->domain->getNodeDistInfo());
  DistSVec<double,dim> R(this->domain->getNodeDistInfo());
  Vec<double> From(nPod);
  expandVector(prom, dQ);

  //----- Obtain an initial bracket -----//

  meritZero = meritFunction(it, Q, dQ, R, 0.0);
  meritOld = meritZero;
  meritDerivZero = meritFunctionDeriv(it, Q, dQ, R, meritZero);

  // Reverse the search direction if it is not a descent direction

  if (meritDerivZero > 0) {
   this->com->fprintf(stderr,"Reversing the search direction: original direction NOT a descent direction! \n");
    dQ = -1.0*dQ;
    prom = -1.0*prom;
    meritDerivZero = -meritDerivZero;
  }
  this->com->fprintf(stderr,"meritDerivZero = %e, meritZero = %e \n",meritDerivZero, meritZero);
 
  targetMeritDerivZero = eps*fabs(meritZero);  //Compare targetMeritDeriv to meritZero based on Taylor series explanation
 
  if (fabs(meritDerivZero) <= targetMeritDerivZero){ // convergence criterion based on stationarity of full-order residual norm
     convergeFlag = 1;
     return 0;
  }

  meritDerivOld = meritDerivZero;
  double alphaOld = 0.0;  //Left endpoint is zero
  
  while (count<maxIter){
    merit = meritFunction(it, Q, dQ, R, alpha); // evaluate merit function at current alpha
    if (merit>meritZero+c1*alpha*meritDerivZero || (merit>=meritOld && count>0)){
        // alphaLo=alphaOld, alphaHi=alpha
        this->com->fprintf(stderr,"Entering zoom from location 1, alpha = %e \n",alpha);
        alpha = zoom(alphaOld, alpha, meritOld, merit, meritDerivOld, meritDerivZero, meritZero, c1, c2, Q, dQ, R, it);
        return alpha;
    }
    Qnew = Q+alpha*dQ;  // Q at current alpha
    meritDeriv = meritFunctionDeriv(it, Qnew, dQ, R, merit);

    if (fabs(meritDeriv)<=-c2*meritDerivZero){
       this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1*alpha*meritDerivZero);
       this->com->fprintf(stderr,"Returning alpha without zoom, alpha = %e \n",alpha);
       return alpha;  // a valid solution has been found
    }
    if (meritDeriv>=0){
       // alphaLo=alpha, alphaHi=alphaOld
       this->com->fprintf(stderr,"Entering zoom from location 2.  The alphas are: alpha = %e, alphaOld = %e \n",alpha,alphaOld);
       alpha = zoom(alpha, alphaOld, merit, meritOld, meritDeriv, meritDerivZero, meritZero, c1, c2, Q, dQ, R, it);
       return alpha;
    }

    //----- update values of alpha, merit function, and derivative -----//

    alphaOld = alpha;
    meritOld = merit;
    meritDerivOld = meritDeriv;

    //----- choose alpha between current alpha and alpha_max -----//

    alpha = beta*alphaOld;   //increment alpha by a constant factor

    //----- correct if alpha exceeds alpha_max -----

    while (alpha>alpha_max){
      beta=(1.0+beta)/2.0;    //decrease beta
      alpha=beta*alphaOld; //try a new value of alpha
      this->com->fprintf(stderr,"Decreasing beta because alpha_max violated \n");
    }
    ++count;
    
  }

  this->com->fprintf(stderr,"Leaving linesearch because max iterations was exceeded. Conditions NOT satisfied! alpha = %e\n",alpha);  
  return alpha;

}
//-----------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::zoom(double alphaLo, double alphaHi, double meritLo, double meritHi, double meritDerivLo, double meritDerivZero, double meritZero, double c1, double c2, DistSVec<double, dim> Q, DistSVec<double, dim> dQ, DistSVec<double, dim> R, int it){

    // Given a valid bracket, this function finds an alpha within the bracket satisfying the Strong Wolfe conditions
    // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

    // initialize parameters
    int count = 0;
    int maxIter = 40; // max number of zooming steps
    double merit, meritDeriv, alpha;
    DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo());
    double epsilon = 0.1; // Percent reduction needed for quadratic approx before bisection is done
    double tolInterp = 0.05; // Determines if interpolation will be executed-we have had problems with the matrix being ill-conditioned
 
    while (count<maxIter){

      // Perform a quadratic interpolation using slope of the low point IF matrix will be well-conditioned using a heuristic

      if (false) { //(fabs(alphaHi-alphaLo)>tolInterp){
      
            FullM Vdm(3); // Vandermonde matrix
   
            Vdm[0][0] = alphaLo*alphaLo; Vdm[0][1] = alphaLo; Vdm[0][2] = 1.0;
            Vdm[1][0] = 2.0*alphaLo; Vdm[1][1] = 1.0; Vdm[1][2] = 0.0;
            Vdm[2][0] = alphaHi*alphaHi; Vdm[2][1] = alphaHi; Vdm[2][2] = 1.0;
            double *VdmRHS = new double[3];
            VdmRHS[2] = meritHi; VdmRHS[1] = meritDerivLo; VdmRHS[0] = meritLo;
            Vdm.Factor();
            Vdm.ReSolve(VdmRHS); //obtain interpolant coefficients
            alpha = -VdmRHS[1]/(2.0*VdmRHS[0]);//new value from minimizing quadratic model
   
            delete[] VdmRHS;
      }

      // ----- Use bisection if interpolation was bypassed or if there was an insufficient reduction-----//

      if (fabs(alphaHi-alphaLo)<=tolInterp || (alpha-alphaLo)/(alphaHi-alphaLo)<epsilon || (alphaHi-alpha)/(alphaHi-alphaLo)<epsilon)
        alpha=(alphaLo+alphaHi)/2.0; // bisection

      // Evaluate merit function

      merit = meritFunction(it, Q, dQ, R, alpha);

      if (merit>meritZero+c1*alpha*meritDerivZero || merit>=meritLo){
           alphaHi = alpha;
           meritHi = merit;
      }
      else{
           Qnew = Q+alpha*dQ;
           meritDeriv = meritFunctionDeriv(it, Qnew, dQ, R, merit);
           if (fabs(meritDeriv)<=-c2*meritDerivZero){
                this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1*alpha*meritDerivZero);
                this->com->fprintf(stderr,"alpha = %e, # zoom iter = %d\n",alpha, count);

                return alpha;  // alpha satisfies Strong Wolfe conditions: exit loop
           }
           if (meritDeriv*(alphaHi-alphaLo)>=0){
                alphaHi = alphaLo;
                meritHi = meritLo;
           }
           alphaLo = alpha;
           meritLo = merit;
           meritDerivLo = meritDeriv; 
      }
      ++count;
    }
    this->com->fprintf(stderr,"Leaving zoom because max iterations was exceeded. Conditions NOT satisfied! Count = %d, alpha = %e\n",count, alpha);
    return alpha;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::expandVector(Vec<double> &romV, DistSVec<double, dim> &fullV)  {

  fullV = 0.0;
  for (int iVec = 0; iVec < nPod; iVec++)
    fullV += romV[iVec]*pod[iVec];

}


//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::projectVector(VecSet<DistSVec<double, dim> > &leftProj, DistSVec<double,dim> &fullV, Vec<double> &romV)  {

	transMatVecProd(leftProj,fullV,projVectorTmp);

  for (int iVec = 0; iVec < pod.numVectors(); iVec++)
    romV[iVec] = projVectorTmp[iVec];

}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{

  this->com->fprintf(stderr, "*** Warning: Checkfailsafe Not yet implemented\n");
  int failSafeNewton = 0;
  return failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::resetFixesTag()
{
  this->spaceOp->resetTag();
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::computeRedHessianSums(int it, DistSVec<double, dim> &Q)  {

  //mvpfd->evaluate(it, *this->X, *this->A, Q, F);
  
  //for (int iPod = 0; iPod < nPod; iPod++)
  //  mvpfd->apply(pod[iPod], AJ[iPod]);

}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting, DistSVec<double, dim> *R)  {

//  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);

//  if (_pc) {
      
//      this->com->fprintf(stdout, "attempting to apply preconditioner\n");
//      this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc, this->timeState);
//      this->timeState->addToJacobian(*this->A, *_pc, Q);
//      this->spaceOp->applyBCsToJacobian(Q, *_pc);

//      _pc->setup();

//  } else {
  if (R==NULL) R = &F;
 
  mvpfd->evaluate(it, *this->X, *this->A, Q, *R);
  
  for (int iPod = 0; iPod < nPod; iPod++)
    mvpfd->apply(pod[iPod], AJ[iPod]);

  //saveNewtonSystemVectors(totalTimeSteps);
 
  // weight AJ
  if (applyWeighting && (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE)) {
    double weightExp = this->ioData->romOnline.weightingExponent;
    double weightNorm = weightVec->norm();
    weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;
    for (int iVec=0; iVec<nPod; ++iVec) {
      int numLocSub = ((AJ)[iVec]).numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double (*weight)[dim] = weightVec->subData(iSub);
        double (*aj)[dim] = ((AJ)[iVec]).subData(iSub);
        for (int i=0; i<weightVec->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            aj[i][j] = aj[i][j] * weight[i][j];
        }
      }
    }
  }

//  }


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::saveNewtonSystemVectorsAction(const int totalTimeSteps) {
	// only do for PG and Galerkin

  int freq = ioData->output.rom.resjacfrequency;

  if ((freq == 0) || (totalTimeSteps%freq == 0)) {
    DistSVec<double, dim> AJsol(this->domain->getNodeDistInfo()); //CBM--NEED TO CHANGE NAME OF DISTVECTOR
  	AJsol = 0.0;
  	for (int i=0; i<this->nPod; ++i)
  		 AJsol += this->AJ[i] * this->dUromNewtonIt[i]; 

	  // saving 1) residual and 2) this->AJ * this->dUromNewtonIt (for GappyPOD)
	  rom->writeClusteredBinaryVectors(currentCluster, &(this->F), &AJsol);
	}
}


//------------------------------------------------------------------------------

//template<int dim>
//void ImplicitRomTsDesc<dim>::savedUnormAccum() {

//	for (int iPod = 0; iPod < nPod; ++iPod) {
//		dUnormAccum[0][iPod] += fabs(dUromNewtonIt[iPod]);	// 1 norm
//		dUnormAccum[1][iPod] += dUromNewtonIt[iPod] * dUromNewtonIt[iPod];	// 2 norm
//	}

//}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::rstVarImplicitRomTsDesc(IoData &ioData)
{

  mvpfd->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop1(const bool breakloop) {

	return breakloop;
	
}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop2(const bool breakloop) {

	return false;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::updateLeastSquaresWeightingVector() {

  // form the weighting vector for the least-squares system

  int numLocSub = this->domain->getNumLocSub();

  switch (this->ioData->romOnline.weightedLeastSquares) {
    case (NonlinearRomOnlineData::WEIGHTED_LS_FALSE):
      return;
      break;
/*    case (NonlinearRomOnlineData::WEIGHTED_LS_RESIDUAL):
      *(this->weightVec) = *(this->weightFRef);
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_STATE):
      *(this->weightVec) = *(this->weightURef);  
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double (*weight)[dim] = this->weightVec->subData(iSub);
        for (int i=0; i<this->weightVec->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = weight[i][j] - (this->bcData->getInletConservativeState())[j];
        }
      }
      break;*/
    case (NonlinearRomOnlineData::WEIGHTED_LS_CV):
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*weight)[dim] = weightVec->subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = cv[i];
        }
      }
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_BOCOS):
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *ffMask = farFieldMask->subData(iSub); // vector with nonzero entries at farfield nodes
        double *wMask = wallMask->subData(iSub); // vector with nonzero entries at wall nodes
        double (*weight)[dim] = weightVec->subData(iSub);
        for (int i=0; i<weightVec->subSize(iSub); ++i) {
          if (ffMask[i]>0) {
            for (int j=0; j<dim; ++j)
              weight[i][j] = ffWeight;
          } else if (wMask[i]>0) {
            for (int j=0; j<dim; ++j)
              weight[i][j] = wallWeight;
          } else {
            for (int j=0; j<dim; ++j)
               weight[i][j] = interiorWeight;
          }
        }
      }
      break;
    default:
        this->com->fprintf(stderr, "*** Error: Unexpected least-squares weighting method\n");
        exit(-1);
      break;
  }

}
