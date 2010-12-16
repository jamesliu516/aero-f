#include <Communicator.h>

#include <cmath>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)	// all cpus along rows
{

	readSampleNodes(ioData.input.sampleNodes);

	// assume we have nPodJac
  nPodJac = ioData.Rob.numROBJac; 

	leastSquaresSolver.blockSizeIs(32);
	leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
	
	// read in Afull, Bfull (temporary) (binary files because in reduced mesh)
  VecSet<DistSVec<double, dim> > Afull(0,dom->getNodeDistInfo());
  VecSet<DistSVec<double, dim> > Bfull(0,dom->getNodeDistInfo());
	dom->readPodBasis(ioData.input.aMatrix, nPodJac,Afull);
	dom->readPodBasis(ioData.input.bMatrix, nPodJac,Bfull);

	//Amat.reset(new VecSet<DistSVec<double, dim> >(0, dom->getNodeDistInfo()));
	//Bmat.reset(new VecSet<DistSVec<double, dim> >(0, dom->getNodeDistInfo()));
	//dom->readPodBasis(ioData.input.aMatrix, nPodJac,*Amat);
	//dom->readPodBasis(ioData.input.bMatrix, nPodJac,*Bmat);

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

	//AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, dom->getNodeDistInfo()));
	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));
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
  for (int iCol = 0; iCol < leastSquaresSolver.unknownCount(); ++iCol) {
    const bool hasLocalCol = (leastSquaresSolver.localCpuCol() == leastSquaresSolver.colHostCpu(iCol));
    const int localICol = leastSquaresSolver.localColIdx(iCol);
    for (int iRow = 0; iRow < leastSquaresSolver.equationCount(); ++iRow) {
      const double entryValue = (*Amat)[iRow] * (*AJRestrict)[iCol];
      const bool hasLocalEntry = hasLocalCol && (leastSquaresSolver.localCpuRow() == leastSquaresSolver.rowHostCpu(iRow));
      if (hasLocalEntry) {
        const int localIRow = leastSquaresSolver.localRowIdx(iRow);
        leastSquaresSolver.matrixEntry(localIRow, localICol) = entryValue;
      }
    }
  }
 
  // Form B * ResRestrict and distribute
  {
    const bool hasLocalRhs = (leastSquaresSolver.localCpuCol() == leastSquaresSolver.rhsRankHostCpu(0));
    for (int iRow = 0; iRow < leastSquaresSolver.equationCount(); ++iRow) {
      const double entryValue = (*Bmat)[iRow] * (*ResRestrict);
      const bool hasLocalEntry = hasLocalRhs && (leastSquaresSolver.localCpuRow() == leastSquaresSolver.rhsRowHostCpu(iRow));
      if (hasLocalEntry) {
        const int localIRow = leastSquaresSolver.localRhsRowIdx(iRow);
        leastSquaresSolver.rhsEntry(localIRow) = -1.0 * entryValue;	// need to make negative
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
		sampleNodes.push_back(currentSampleNode-1);	// reads in the sample node plus one
	}
}
