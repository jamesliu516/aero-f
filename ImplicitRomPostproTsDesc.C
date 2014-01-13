#include <Communicator.h>

#include <cmath>
#include <iostream>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomPostproTsDesc<dim>::ImplicitRomPostproTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), Uinitial(dom->getNodeDistInfo())
{
	this->maxItsNewton = 0;	// never do iterations

  reducedCoordsFile = NULL;

  if (!(strcmp(this->ioData->input.reducedCoords,"")==0)) {
    char *fullReducedCoordsName = new char[strlen(this->ioData->input.prefix) + 1 + strlen(this->ioData->input.reducedCoords) + 1];
    sprintf(fullReducedCoordsName, "%s%s", this->ioData->input.prefix, this->ioData->input.reducedCoords);
    reducedCoordsFile = fopen(fullReducedCoordsName, "r");
    delete [] fullReducedCoordsName;
  }

 if (!reducedCoordsFile)  {
   this->com->fprintf(stderr, "*** Error opening reduced coordinates file: %s\n", reducedCoordsFile);
   exit (-1);
 }

  // avoid calling computeTimeStep for steady
  dt = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeFullResidual(int , DistSVec<double, dim> &) {


	// do nothing
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeAJ(int, DistSVec<double, dim> &)  {

	// do nothing
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  breakloop = true;	// after loop, exit will occur because maxItsNewton = 1
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::checkLocalRomStatus(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  // checks whether the local ROM needs to be modified

  int tmp, _n, closestCluster, nCoords; 
  char switchStr[50], updateStr[50];


  _n = fscanf(reducedCoordsFile, "%d %s %s %d %d %d %d %d", &tmp, switchStr, updateStr, &closestCluster, &nCoords, &tmp, &tmp, &tmp);

  this->com->fprintf(stdout, "%s %s %d %d\n", switchStr, updateStr, closestCluster, nCoords);

  this->clusterSwitch = (strcmp(switchStr,"switch")==0) ? true : false;
  this->updateFreq = (strcmp(updateStr,"update")==0) ? true : false;

  if ((this->rom->nClusters > 1) || (this->basisUpdateFreq>0) || (this->currentCluster == -1))  {
    //rom->closestCenter(U, &closestCluster);

    if (this->rom->nClusters > 1) this->com->fprintf(stdout, " ... using basis number %d\n", closestCluster);

    //updateFreq = ((basisUpdateFreq > 0) && (totalTimeSteps%basisUpdateFreq == 0)) ? true : false;
    //clusterSwitch = (currentCluster != closestCluster) ? true : false;

    if (this->updateFreq || this->clusterSwitch) {
      if (this->clusterSwitch) {
        this->currentCluster = closestCluster;
        this->rom->readClusteredOnlineQuantities(this->currentCluster);  // read state basis, update info, and (if applicable) gnat online matrices
      }

      if (this->ioData->romOnline.basisUpdates) this->rom->updateBasis(this->currentCluster, U);
      //if (this->ioData->romOnline.krylov.include) rom->appendNonStateDataToBasis(currentCluster,"krylov");
      //if (this->ioData->romOnline.sensitivity.include) rom->appendNonStateDataToBasis(currentCluster,"sensitivity");

      this->nPod = this->rom->basis->numVectors();
      this->pod.resize(this->nPod);
      for (int iVec=0; iVec<this->nPod; ++iVec) {
        this->pod[iVec] = (*(this->rom->basis))[iVec];
      }

      //AJ.resize(nPod);
      this->dUromTimeIt.resize(this->nPod);
      //setProblemSize(U);  // defined in derived classes
      //if (clusterSwitch) setReferenceResidual(); // for steady gnat (reference residual is restricted to currently active nodes)

    }
  }

  if (this->nPod != nCoords) {
    this->com->fprintf(stderr, "*** Error: dimension of reduced coordinates (%d) does not match dimension of state ROB (%d)\n"
                       ,nCoords, this->nPod);
    exit(-1);
  }

  double tmp2;
  for (int iPod = 0; iPod < this->nPod; ++iPod) {
    _n = fscanf(reducedCoordsFile, "%le", &tmp2);
    this->dUromTimeIt[iPod] = tmp2;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::postProStep(DistSVec<double, dim> &U, int totalTimeSteps)  {

  DistSVec<double, dim> dU(this->domain->getNodeDistInfo());
	this->expandVector(this->dUromTimeIt, dU); // solution increment in full coordinates
	U += dU;

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomPostproTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// avoid residual calculation

    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitRomPostproTsDesc<dim>::computeTimeStep(int it, double *dtLeft, DistSVec<double,dim> &U, double angle)
{
  double t0 = this->timer->getTime();
  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual, angle);
  int numSubCycles = 1;

  if(it==1) {
    if(this->failSafeFlag == false){
      dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft, &numSubCycles, *this->geoState, *this->X, *this->A, U);
    } else {
      dt = this->timeState->computeTimeStepFailSafe(dtLeft, &numSubCycles);
    }
  }

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*this->refVal->time, numSubCycles);
  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;
}


