#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::ImplicitGnatTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{

  nPodJac = 0;  // set when reading the online matrices
  numResJacMat = this->rom->getNumResJacMat();

	leastSquaresSolver.blockSizeIs(32);

  jactmp = NULL;
  column = NULL;
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::~ImplicitGnatTsDesc() 
{
	if (jactmp) delete [] jactmp;
	if (column) delete [] column;
  
  ResRestrict.reset();
  AJRestrict.reset();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::deleteRestrictedQuantities() {

  ResRestrict.reset();
  AJRestrict.reset();
  
  this->rom->deleteRestrictedQuantities();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy)
{
	// Evaluate residual on full mesh

  if (R==NULL) R=&(this->F);
 
  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, *R, this->timeState, *(this->rom->restrictMapping()));

  if (includeHomotopy) {
    this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Q, *R, (this->rom->restrictMapping())->getRestrictedToOriginLocNode());
  }

  this->spaceOp->applyBCsToResidual(Q, *R);

  if (applyWeighting && (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE)) {
    // weight residual
    double weightExp = this->ioData->romOnline.weightingExponent;
    double weightNorm = this->weightVec->norm();
    weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;
    int numLocSub = R->numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double (*weight)[dim] = this->weightVec->subData(iSub);
      double (*r)[dim] = R->subData(iSub);
      for (int i=0; i<this->weightVec->subSize(iSub); ++i) {
        for (int j=0; j<dim; ++j) {
          //weight[i][j] = pow(abs(weight[i][j])/weightNorm, weightExp);
          r[i][j] = r[i][j] * weight[i][j];
        }
      }
    }
  }

  double t0 = this->timer->getTime();
  (this->rom->restrictMapping())->restriction(*R, *ResRestrict);
  this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting, DistSVec<double, dim> *R)  {

	// Evaluate action of Jacobian on full mesh
  if (R==NULL) R = &this->F;

	this->mvpfd->evaluateRestrict(it, *this->X, *this->A, Q, *R,
			*(this->rom->restrictMapping()));	// very cheap
  
  for (int iPod = 0; iPod < this->nPod; iPod++) {
		this->mvpfd->applyRestrict(this->pod[iPod], this->AJ[iPod],
				*(this->rom->restrictMapping()));
	}

  if (applyWeighting && (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE)) {
    double weightExp = this->ioData->romOnline.weightingExponent;
    double weightNorm = this->weightVec->norm();
    weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      int numLocSub = (this->AJ[iVec]).numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double (*weight)[dim] = this->weightVec->subData(iSub);
        double (*aj)[dim] = (this->AJ[iVec]).subData(iSub);
        for (int i=0; i<this->weightVec->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            aj[i][j] = aj[i][j] * weight[i][j];
        }
      }
    }
  }

	double t0 = this->timer->getTime();
	for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
		(this->rom->restrictMapping())->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
	}
	this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  // Form A * of and distribute

  double t0 = this->timer->getTime();
  transMatMatProd(*(this->rom->getJacMat()), *AJRestrict,jactmp);

  for (int iCol = 0; iCol < leastSquaresSolver.localCols(); ++iCol) {
    const int globalColIdx = leastSquaresSolver.globalColIdx(iCol);
    const int colOffset = globalColIdx * leastSquaresSolver.equationCount();
    for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
      const int globalRowIdx = leastSquaresSolver.globalRowIdx(iRow);
      leastSquaresSolver.matrixEntry(iRow, iCol) = jactmp[globalRowIdx + colOffset];
    }
  }

  // Form B * ResRestrict and distribute
  // NOTE: do not check if current CPU has rhs

  transMatVecProd(*(this->rom->getResMat()), *ResRestrict, column);

  for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
    const int globalRowIdx = leastSquaresSolver.globalRowIdx(iRow);
    leastSquaresSolver.rhsEntry(iRow) = -column[globalRowIdx];
  }
  this->timer->addLinearSystemFormTime(t0);

  // Solve least squares problem
  t0 = this->timer->getTime();
  leastSquaresSolver.solve();
  this->timer->addLinearSystemSolveTime(t0);

  // Update vector: The first nPod rows give the components in the pod basis
  t0 = this->timer->getTime();
  this->dUromNewtonIt = 0.0;
  for (int localIRow = 0; localIRow < leastSquaresSolver.localSolutionRows(); ++localIRow) {
    const int iRow = leastSquaresSolver.globalRhsRowIdx(localIRow);
    this->dUromNewtonIt[iRow] = leastSquaresSolver.rhsEntry(localIRow);
  }
 
  // Consolidate across the cpus
  this->com->globalSum(this->nPod, this->dUromNewtonIt.data());
  res = this->dUromNewtonIt.norm();

  // Convergence criterion
  if (it == 0) {
    this->res0 = res;
    this->target = this->epsNewton * this->res0;
  }

  breakloop = (res == 0.0) || (res <= this->target);
  this->timer->addCheckConvergenceTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::breakloop1(const bool breakloop) {

	return false;
	
}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::breakloop2(const bool breakloop) {

	return breakloop;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitGnatTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &F, double stepLength)  {
	// merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  computeFullResidual(it,newQ,true,&F);

  double merit = 0.0;
  merit += ResRestrict->norm();	// merit function = 1/2 * (norm of full-order residual)^2
  merit *= merit;
//  merit *= 0.5;

  //if (this->ioData->romOnline.lsSolver == 2) {
  //  DistSVec<double, dim>* A_Uerr = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  //  *A_Uerr = 0.0;
  //  int numLocSub = A_Uerr->numLocSub();
//#pragma omp parallel for
  //  for (int iSub=0; iSub<numLocSub; ++iSub) {
  //    double *cv = this->A->subData(iSub); // vector of control volumes
  //    double (*auerr)[dim] = A_Uerr->subData(iSub);
  //    double (*u)[dim] = newQ.subData(iSub);
  //    for (int i=0; i<this->A->subSize(iSub); ++i) {
  //      if (cv[i]>regThresh) {
  //        for (int j=0; j<dim; ++j)
  //          auerr[i][j] = cv[i] * pow(u[i][j] - (this->bcData->getInletConservativeState())[j],2);
  //      }
  //    }
  //  }

  //  DistSVec<double, dim>* ones = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  //  *ones = 1.0;

  //  double regTerm = (*A_Uerr)*(*ones);
  //  regTerm *= regWeight;
  //  merit += regTerm;

  //  delete ones;
  //  delete A_Uerr;
  //}

  return merit;

}


//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {
 
  nPodJac = this->rom->getJacMat()->numVectors();

	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);

  if (jactmp) delete [] jactmp;
  if (column) delete [] column;
	jactmp = new double [nPodJac * this->nPod];
	column = new double [nPodJac];

	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->rom->getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (this->rom->getRestrictedDistInfo())); 
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));  

  if (!this->unsteady) {
    this->com->printf(2, "\nNOTE: For steady GNAT simulations the reported residual is calculated using only the sample mesh,\n");
    this->com->printf(2, "      and is relative to the residual of the initial condition calculated on the same sample mesh.\n");
    this->com->printf(2, "      (This reference residual is re-restricted after every cluster switch for consistency).\n");
 
    this->Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *(this->Uinit) = U;  // needed for computing the restricted residual after each cluster switch
  }

  this->com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if (!this->unsteady && monitorConvergence(it, U))
    return true;

  if (!this->problemType[ProblemData::AERO] && !this->problemType[ProblemData::THERMO] && it >= this->data->maxIts) return true;

  if (this->unsteady)
    if(t >= this->data->maxTime - 0.01 * dt)
      return true;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// only called for steady simulations

  this->data->residual = computeGnatResidualNorm(U);

  if (this->data->residual == 0.0 || this->data->residual < this->data->eps * this->restart->residual || this->data->residual < this->data->epsabs)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitGnatTsDesc<dim>::computeGnatResidualNorm(DistSVec<double,dim>& Q)
{ // spatial only

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, *(this->rom->restrictMapping()));

  this->spaceOp->applyBCsToResidual(Q, this->F);

  double t0 = this->timer->getTime();

  (this->rom->restrictMapping())->restriction(this->F, *ResRestrict); 

  this->timer->addRestrictionTime(t0);

  double res = 0.0;
  res = (*ResRestrict) * (*ResRestrict);

  return sqrt(res);

}


//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::setReferenceResidual()
{
  if (this->Uinit) this->restart->residual = computeGnatResidualNorm(*(this->Uinit));

  this->com->printf(2, "Norm of restricted reference residual = %.12e\n", this->restart->residual);

}


