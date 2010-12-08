//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)	// all cpus along rows
{

	readSampleNodes(ioData.input.sampleNodes);

	// assume we have nPodJac
  nPodJac = ioData.Rob.numROBJac; 
  dom->readPodBasis(ioData.input.aMatrix, nPodJac, *Amat);
  dom->readPodBasis(ioData.input.bMatrix, nPodJac, *Bmat);

	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
	
	// read in Afull, Bfull (temporary) (binary files because in reduced mesh)
  VecSet<DistSVec<double, dim> > Afull(nPodJac,dom->getNodeDistInfo());
  VecSet<DistSVec<double, dim> > Bfull(nPodJac,dom->getNodeDistInfo());
	//dom->readPodBasis(ioData.input.AFile, nPodJac,Afull);
	//dom->readPodBasis(ioData.input.BFile, nPodJac,Bfull);

	// determine mapping to restricted nodes
	restrictionMapping.reset(new RestrictionMapping<dim>(dom, sampleNodes.begin(), sampleNodes.end()));

	// allocate memory for Amat, Bmat using restrictedDistInfo
	Amat.reset(new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo()));
	Bmat.reset(new VecSet<DistSVec<double, dim> >(nPodJac, getRestrictedDistInfo()));

	// restrict Afull and Bfull to be Amat, Bmat
	for (int i = 0; i < nPodJac; ++i) {
		restrictionMapping->restriction(Afull[i],(*Amat)[i]);
		restrictionMapping->restriction(Bfull[i],(*Bmat)[i]);
	}

	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeFullResidual(int it, DistSVec<double,
dim> &Q) {
 // Evaluate residual on full mesh
 ImplicitRomTsDesc<dim>::computeFullResidual(it, Q);

 // Restrict down
 restrictMapping()->restriction(this->F, *ResRestrict);
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action fo Jacobian on full mesh

 this->mvpfd->evaluate(it, *this->X, *this->A, Q, this->F);

 DistSVec<double, dim> AJfull(restrictMapping()->originDistInfo());

 for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
   this->mvpfd->apply(this->pod[iPod], AJfull);
   restrictMapping()->restriction(AJfull, (*AJRestrict)[iPod]);
 }

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::solveNewtonSystem(const int &it, double &res,
bool &breakloop)  {

	// Form A * AJRestrict and distribute

	//TODO FIX THIS
 for (int iCol = 0; iCol < leastSquaresSolver.localCols(); ++ iCol) {
   const int i_y = leastSquaresSolver.globalColIdx(iCol); 
   for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
     const int i_J = leastSquaresSolver.globalRowIdx(iRow);
     leastSquaresSolver.matrixEntry(iRow, iCol) = (*Amat)[i_J] * (*AJRestrict)[i_y];
   }
 }

	// Form B * ResRestrict and distribute
 for (int iRow = 0; iRow < leastSquaresSolver.localRhsRows(); ++iRow) {
   const int i_J = leastSquaresSolver.globalRhsRowIdx(iRow);
   leastSquaresSolver.rhsEntry(iRow) = (*Bmat)[i_J] * (*ResRestrict);
 }

	// Solve least squares problem
 leastSquaresSolver.solve();

 // Compute residual error and update vector
 this->dUrom = 0.0;
 double residualError = 0.0;

 // The first nPod rows give the components in the pod basis
 for (int iRow = 0; iRow < leastSquaresSolver.localSolutionRows(); ++iRow) {
   const int i_y = leastSquaresSolver.globalRhsRowIdx(iRow);
   this->dUrom[i_y] = leastSquaresSolver.rhsEntry(iRow);
 }

 // The sum of squares of the remaining rows give the square of the residual
 // norm
 for (int iRow = leastSquaresSolver.localSolutionRows(); iRow <
leastSquaresSolver.localRhsRows(); ++iRow) {
   const double entry = leastSquaresSolver.rhsEntry(iRow);
   residualError += entry * entry;
 }

 // Consolidate across the cpus
 this->com->globalSum(this->nPod, this->dUrom.data());
 this->com->globalSum(1, &residualError);
 res = sqrt(residualError);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::readSampleNodes(const char *sampleNodeFileName)  {

	// INPUT: sample node file name
	// OUTPUT: nSampleNodes, sampleNodes

	FILE *sampleNodeFile = fopen(sampleNodeFileName, "r");
	fscanf(sampleNodeFile, "%d",&nSampleNodes);	// first entry is the number of sample nodes
	sampleNodes.reserve(nSampleNodes);	// know it will be nSampleNodes long (efficiency)

	int index, currentSampleNode;
	for (int i = 0; i < nSampleNodes; ++i){
		fscanf(sampleNodeFile, "%d",&index);
		fscanf(sampleNodeFile, "%d",&currentSampleNode);
		sampleNodes.push_back(currentSampleNode);
	}
}
