#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{

	// NOTE: Amat corresponds to RESIDUAL, Bmat corresponds to JACOBIAN
  nPodJac = ioData.Rob.numROBJac; 
  performRestriction = ioData.Rob.performRestriction == 1 ? true : false; 

	nSampleNodes = 0;

	if (ioData.Rob.sampleNodeFactor != -1) {
		nSampleNodes = static_cast<int>(ceil(double(nPodJac *
						ioData.Rob.sampleNodeFactor)/double(dim)));
	}

	dom->readSampleNodes(sampleNodes, nSampleNodes, this->input->sampleNodes);

	// assume we have nPodJac

	leastSquaresSolver.blockSizeIs(32);
	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
	
	// read in Afull, Bfull (temporary) (binary files because in reduced mesh)
  VecSet<DistSVec<double, dim> > Afull(0,dom->getNodeDistInfo());
  VecSet<DistSVec<double, dim> > Bfull(0,dom->getNodeDistInfo());
	if (*(this->input->aMatrix) != '\0')
		dom->readPodBasis(this->input->aMatrix, nPodJac,Afull);
	if (*(this->input->bMatrix) == '\0') {	// not specified
		numABmat = 1;	// same matrix
	}
	else {
		numABmat = 2;	// different matrices
	}

	// determine mapping to restricted nodes
	restrictionMapping.reset(new RestrictionMapping<dim>(dom, sampleNodes.begin(), sampleNodes.end()));

	// allocate memory for Amat, Bmat using restrictedDistInfo
	Amat = new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo());

	if (numABmat == 1) {
		Bmat = Amat;
	}
	else {
		dom->readPodBasis(this->input->bMatrix, nPodJac,Bfull);
		Bmat = new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo());
	}

	// restrict Afull and Bfull to be Amat, Bmat
	for (int i = 0; i < nPodJac; ++i) {
		restrictionMapping->restriction(Afull[i],(*Amat)[i]);
		if (numABmat == 2)
			restrictionMapping->restriction(Bfull[i],(*Bmat)[i]);
	}

	jactmp = new double [nPodJac * this->nPod];
	column = new double [nPodJac];

	if (performRestriction) {
		AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
		ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));
	}

}

template<int dim>
ImplicitGappyTsDesc<dim>::~ImplicitGappyTsDesc() 
{
	delete Amat;
	if (numABmat == 2)
		delete Bmat;
	if (jactmp) delete [] jactmp;
	if (column) delete [] column;
}
//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q) {

	// Evaluate residual on full mesh

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, *restrictionMapping);

	this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Q,
			this->F, restrictionMapping->getRestrictedToOriginLocNode());

  this->spaceOp->applyBCsToResidual(Q, this->F);

	double t0 = this->timer->getTime();
	// Restrict down
	if (performRestriction) {
		restrictMapping()->restriction(this->F, *ResRestrict);
	}
	this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action of Jacobian on full mesh

	this->mvpfd->evaluateRestrict(it, *this->X, *this->A, Q, this->F,
			*restrictionMapping);	// very cheap
  
  for (int iPod = 0; iPod < this->nPod; iPod++) {
		this->mvpfd->applyRestrict(this->pod[iPod], this->AJ[iPod],
				*restrictionMapping);
	}

	double t0 = this->timer->getTime();
	if (performRestriction) {
		for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
			restrictMapping()->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
		}
	}
	this->timer->addRestrictionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {
  // Form A * of and distribute

	double t0 = this->timer->getTime();
	if (performRestriction) {
		transMatMatProd(*Bmat, *AJRestrict,jactmp);
	}
	else {
		transMatMatProdRestrict(*Bmat, this->AJ,jactmp,restrictionMapping->getRestrictedToOriginLocNode());
	}

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

	if (performRestriction) {
		transMatVecProd(*Amat, *ResRestrict, column);
	}
	else {
		transMatVecProdRestrict(*Amat, this->F, column,restrictionMapping->getRestrictedToOriginLocNode());
	}

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
