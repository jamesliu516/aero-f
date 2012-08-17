#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::ImplicitGnatTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{
/*
	// read in gappy POD matrix for residual
  VecSet<DistSVec<double, dim> > resMatrixFull(0,dom->getNodeDistInfo());
	gnatPrefix = ioData.input.gnatPrefix;
	string fileNameRes;
	ImplicitRomTsDesc<dim>::determineFileName(this->input->resMatrix, ".gappyRes", gnatPrefix, fileNameRes);
  nPodJac = -1;
	dom->readPodBasis(fileNameRes.c_str(), nPodJac, resMatrixFull);

	// read in gappy POD matrix for Jacobian
  VecSet<DistSVec<double, dim> > jacMatrixFull(0,dom->getNodeDistInfo());
	string fileNameJac;
	ImplicitRomTsDesc<dim>::determineFileName(this->input->jacMatrix, ".gappyJac", gnatPrefix, fileNameJac);
	ifstream jacMatFile(fileNameJac.c_str());
	if (jacMatFile.good()) {	// not specified
		numResJacMat = 2;	// same matrix
	}
	else {
		numResJacMat = 1;	// different matrices
	}

	// read in sample nodes
	nSampleNodes = 0;
	string fileNameSample;
	ImplicitRomTsDesc<dim>::determineFileName(this->input->sampleNodes, ".sampledNodes", gnatPrefix, fileNameSample);
	dom->readSampleNodes(sampleNodes, nSampleNodes, fileNameSample.c_str());

	// assume we have nPodJac
	leastSquaresSolver.blockSizeIs(32);
	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
	
	// read in resMatrixFull, jacMatrixFull (temporary) (binary files because in reduced mesh)

	// determine mapping to restricted nodes
	restrictionMapping.reset(new RestrictionMapping<dim>(dom, sampleNodes.begin(), sampleNodes.end()));

	// allocate memory for resMat, jacMat using restrictedDistInfo
	resMat = new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo());

	if (numResJacMat == 1) {
		jacMat = resMat;
	}
	else {
		dom->readPodBasis(fileNameJac.c_str(), nPodJac, jacMatrixFull);
		jacMat = new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo());
	}

	// restrict resMatrixFull and jacMatrixFull to be resMat, jacMat
	for (int i = 0; i < nPodJac; ++i) {
		restrictionMapping->restriction(resMatrixFull[i],(*resMat)[i]);
		if (numResJacMat == 2)
			restrictionMapping->restriction(jacMatrixFull[i],(*jacMat)[i]);
	}

	jactmp = new double [nPodJac * this->nPod];
	column = new double [nPodJac];

	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));
*/
}

template<int dim>
ImplicitGnatTsDesc<dim>::~ImplicitGnatTsDesc() 
{/*
	delete resMat;
	if (numResJacMat == 2)
		delete jacMat;
	if (jactmp) delete [] jactmp;
	if (column) delete [] column;
*/
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q) {

	// Evaluate residual on full mesh
/*
  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, *restrictionMapping);

	this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Q,
			this->F, restrictionMapping->getRestrictedToOriginLocNode());

  this->spaceOp->applyBCsToResidual(Q, this->F);

	double t0 = this->timer->getTime();

	restrictMapping()->restriction(this->F, *ResRestrict);

	this->timer->addRestrictionTime(t0);
*/
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action of Jacobian on full mesh
/*
	this->mvpfd->evaluateRestrict(it, *this->X, *this->A, Q, this->F,
			*restrictionMapping);	// very cheap
  
  for (int iPod = 0; iPod < this->nPod; iPod++) {
		this->mvpfd->applyRestrict(this->pod[iPod], this->AJ[iPod],
				*restrictionMapping);
	}

	double t0 = this->timer->getTime();
	for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
		restrictMapping()->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
	}
	this->timer->addRestrictionTime(t0);
*/
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {
  // Form A * of and distribute
/*
	double t0 = this->timer->getTime();
	transMatMatProd(*jacMat, *AJRestrict,jactmp);

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

	transMatVecProd(*resMat, *ResRestrict, column);

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
*/
}

template<int dim>
bool ImplicitGnatTsDesc<dim>::breakloop1(const bool breakloop) {

	return false;
	
}

template<int dim>
bool ImplicitGnatTsDesc<dim>::breakloop2(const bool breakloop) {

	return breakloop;

}
