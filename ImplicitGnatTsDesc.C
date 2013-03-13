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
  
  Uinit = NULL;
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::~ImplicitGnatTsDesc() 
{
	if (jactmp) delete [] jactmp;
	if (column) delete [] column;
  if (Uinit) delete Uinit;
  
  ResRestrict.reset();
  AJRestrict.reset();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q) {

	// Evaluate residual on full mesh

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, *(this->rom->restrictMapping()));

	this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Q,
			this->F, (this->rom->restrictMapping())->getRestrictedToOriginLocNode());

  this->spaceOp->applyBCsToResidual(Q, this->F);

	double t0 = this->timer->getTime();

	(this->rom->restrictMapping())->restriction(this->F, *ResRestrict);

	this->timer->addRestrictionTime(t0);


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action of Jacobian on full mesh

	this->mvpfd->evaluateRestrict(it, *this->X, *this->A, Q, this->F,
			*(this->rom->restrictMapping()));	// very cheap
  
  for (int iPod = 0; iPod < this->nPod; iPod++) {
		this->mvpfd->applyRestrict(this->pod[iPod], this->AJ[iPod],
				*(this->rom->restrictMapping()));
	}

	double t0 = this->timer->getTime();
	for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
		(this->rom->restrictMapping())->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
	}
	this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U)  {

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
  this->dUrom = 0.0;
  for (int localIRow = 0; localIRow < leastSquaresSolver.localSolutionRows(); ++localIRow) {
    const int iRow = leastSquaresSolver.globalRhsRowIdx(localIRow);
    this->dUrom[iRow] = leastSquaresSolver.rhsEntry(localIRow);
  }
 
  // Consolidate across the cpus
  this->com->globalSum(this->nPod, this->dUrom.data());
  res = this->dUrom.norm();

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

  if (!this->problemType[ProblemData::UNSTEADY]) {
    this->com->printf(2, "\nNOTE: For steady GNAT simulations the reported residual is calculated using only the sample mesh,\n");
    this->com->printf(2, "      and is relative to the residual of the initial condition calculated on the same sample mesh.\n");
    this->com->printf(2, "      (This reference residual is re-restricted after every cluster switch for consistency).\n");
 
    Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *Uinit = U;  // needed for computing the restricted residual after each cluster switch
  }

  this->com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if (!this->problemType[ProblemData::UNSTEADY] && monitorConvergence(it, U))
    return true;

  if (!this->problemType[ProblemData::AERO] && !this->problemType[ProblemData::THERMO] && it >= this->data->maxIts) return true;

  if (this->problemType[ProblemData::UNSTEADY] )
    if(t >= this->data->maxTime - 0.01 * dt)
      return true;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGnatTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// only called for steady simulations

  this->data->residual = computeGnatResidualNorm(U);

  if (this->data->residual == 0.0 || this->data->residual < this->data->eps * this->restart->residual)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitGnatTsDesc<dim>::computeGnatResidualNorm(DistSVec<double,dim>& Q)
{ // spatial only

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, *this->R, this->timeState, *(this->rom->restrictMapping()));

  this->spaceOp->applyBCsToResidual(Q, *this->R);

  double t0 = this->timer->getTime();

  (this->rom->restrictMapping())->restriction(*this->R, *ResRestrict); 

  this->timer->addRestrictionTime(t0);

  double res = 0.0;
  res = (*ResRestrict) * (*ResRestrict);

  return sqrt(res);

}


//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::setReferenceResidual()
{
  if (Uinit) this->restart->residual = computeGnatResidualNorm(*Uinit);

  this->com->printf(2, "Norm of restricted reference residual = %.12e\n", this->restart->residual);

}


