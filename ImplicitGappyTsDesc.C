#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)	// all cpus along rows
{

	// NOTE: Amat corresponds to RESIDUAL, Bmat corresponds to JACOBIAN
	dom->readSampleNodes(sampleNodes, nSampleNodes, this->input->sampleNodes);

	// assume we have nPodJac
  nPodJac = ioData.Rob.numROBJac; 

	leastSquaresSolver.blockSizeIs(32);
	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
	
	// read in Afull, Bfull (temporary) (binary files because in reduced mesh)
  VecSet<DistSVec<double, dim> > Afull(0,dom->getNodeDistInfo());
  VecSet<DistSVec<double, dim> > Bfull(0,dom->getNodeDistInfo());
	dom->readPodBasis(this->input->aMatrix, nPodJac,Afull);
	if (*(this->input->bMatrix) == '\0') {	// not specified
		numABmat = 1;	// same matrix
	}
	else {
		numABmat = 2;	// different matrices
	}

	//Amat.reset(new VecSet<DistSVec<double, dim> >(0, dom->getNodeDistInfo()));
	//Bmat.reset(new VecSet<DistSVec<double, dim> >(0, dom->getNodeDistInfo()));
	//dom->readPodBasis(ioData.input.aMatrix, nPodJac,*Amat);
	//dom->readPodBasis(ioData.input.bMatrix, nPodJac,*Bmat);

	// determine mapping to restricted nodes
	restrictionMapping.reset(new RestrictionMapping<dim>(dom, sampleNodes.begin(), sampleNodes.end()));

	// allocate memory for Amat, Bmat using restrictedDistInfo
	//Amat.reset(new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo()));
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

	//AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, dom->getNodeDistInfo()));
	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));

	jactmp = new double [Amat->numVectors() * AJRestrict->numVectors()];
	column = new double [Amat->numVectors()];
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

 ImplicitRomTsDesc<dim>::computeFullResidual(it, Q);

 // Restrict down
 restrictMapping()->restriction(this->F, *ResRestrict);

 //*ResRestrict=this->F;
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action of Jacobian on full mesh

 ImplicitRomTsDesc<dim>::computeAJ(it, Q);

 for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
   restrictMapping()->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
	 //(*AJRestrict)[iPod] = this->AJ[iPod];
 }

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {
  // Form A * of and distribute
	// TODO: don't recreate logic
	transMatMatProd(*Bmat, *AJRestrict,jactmp);
  for (int iCol = 0; iCol < leastSquaresSolver.unknownCount(); ++iCol) {
    const bool hasLocalCol = (leastSquaresSolver.localCpuCol() == leastSquaresSolver.colHostCpu(iCol));
    const int localICol = leastSquaresSolver.localColIdx(iCol);
    for (int iRow = 0; iRow < leastSquaresSolver.equationCount(); ++iRow) {
      const bool hasLocalEntry = hasLocalCol && (leastSquaresSolver.localCpuRow() == leastSquaresSolver.rowHostCpu(iRow));
      if (hasLocalEntry) {
        const int localIRow = leastSquaresSolver.localRowIdx(iRow);
        leastSquaresSolver.matrixEntry(localIRow, localICol) = jactmp[iRow + iCol * leastSquaresSolver.equationCount()];
      }
    }
  }
 
  // Form B * ResRestrict and distribute
  {
    const bool hasLocalRhs = (leastSquaresSolver.localCpuCol() == leastSquaresSolver.rhsRankHostCpu(0));
		transMatVecProd(*Amat, *ResRestrict, column);
    for (int iRow = 0; iRow < leastSquaresSolver.equationCount(); ++iRow) {
      const bool hasLocalEntry = hasLocalRhs && (leastSquaresSolver.localCpuRow() == leastSquaresSolver.rhsRowHostCpu(iRow));
      if (hasLocalEntry) {
        const int localIRow = leastSquaresSolver.localRhsRowIdx(iRow);
        leastSquaresSolver.rhsEntry(localIRow) = -1.0 * column[iRow];
      }
    }
  }

  // Solve least squares problem
  leastSquaresSolver.solve();

  // Update vector: The first nPod rows give the components in the pod basis
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

}
