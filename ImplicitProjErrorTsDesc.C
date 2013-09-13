//------------------------------------------------------------------------------

template<int dim>
ImplicitProjErrorTsDesc<dim>::ImplicitProjErrorTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
	ImplicitRomTsDesc<dim>(ioData, geoSource, dom), Uold(dom->getNodeDistInfo()),
	Unew(dom->getNodeDistInfo()), dUtrue(dom->getNodeDistInfo()),
	dUproj(dom->getNodeDistInfo()), dUerr(dom->getNodeDistInfo()) {

		// assume vectors are orthonormal if not a snapshot basis
		snapshotBasis = (this->ioData->rom.basisType == 0);
		noBasis = (this->ioData->rom.basisType == 2);	// no basis used: just for generating output files
		this->maxItsNewton = 0;	// do 1 iteration

		// jac = phi^T phi

		if (snapshotBasis) {
			this->jac.setNewSize(this->nPod,this->nPod);
			double *result = new double [this->pod.numVectors() * this->pod.numVectors()];
			transMatMatProd(this->pod,this->pod,result);
			for (int iRow = 0; iRow < this->nPod; ++iRow) {
				for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
					this->jac[iRow][iCol] = result[iRow + iCol * this->pod.numVectors()];
				}
			} 
			this->jac.factor();
		}

		this->projVectorTmp = new double [this->nPod];

}

template<int dim>
ImplicitProjErrorTsDesc<dim>::~ImplicitProjErrorTsDesc(){

	if (this->projVectorTmp) delete [] this->projVectorTmp;

}
//------------------------------------------------------------------------------

template<int dim>
void ImplicitProjErrorTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

}

//------------------------------------------------------------------------------ 


template<int dim>
void ImplicitProjErrorTsDesc<dim>::computeFullResidual(int , DistSVec<double, dim> &) {

	// do nothing
}

template<int dim>
void ImplicitProjErrorTsDesc<dim>::computeAJ(int, DistSVec<double, dim> &)  {

	// do nothing
}

template<int dim>
void ImplicitProjErrorTsDesc<dim>::postProStep(DistSVec<double, dim> &U, int totalTimeSteps)  {

	if (totalTimeSteps == 1) {	// first time step
		Uold = U;
	}

	U = Uold;	// take snapshot value as the state

	// compute new state
	bool endOfFile = this->domain->readVectorFromFile(this->input->snapFile, totalTimeSteps, 0, Unew);

	if (noBasis) {
		U = Unew;
		return;
	}
	else {
		dUtrue = Unew - Uold;

		// compute dUproj = phi(phi^Tphi)^{-1}phi^TdUtrue

		this->dUrom = 0.0;
		this->projectVector(this->pod,dUtrue,this->dUrom);	// Phi^T dUtrue
		// output redcoords

		if (snapshotBasis){
			double *x = this->dUrom.data();
			this->jac.reSolve(x);
		}

		// dUtrue
		this->expandVector(this->dUrom, dUproj);
		dUerr = dUtrue - dUproj;

		// error measures
		dUerrnorm = dUerr.norm();
		double dUtruenorm = dUtrue.norm();
		dUrelerr = dUerrnorm/dUtruenorm;

		U += dUproj;
		Uold = Unew;
	}
}

template<int dim>
void ImplicitProjErrorTsDesc<dim>::writeStateRomToDisk(int it, double cpu)  {

	if (!noBasis)
		this->output->writeStateRomToDisk(it, cpu, this->nPod, this->dUrom);

}

template<int dim>
void ImplicitProjErrorTsDesc<dim>::writeErrorToDisk(int it, double cpu)  {
	// write out change in state coordinate

	//RelativeError AbsoluteError
	if (!noBasis) {
		double errorMeasure [2] = {dUrelerr, dUerrnorm};
		this->output->writeErrorToDisk(it, cpu, 2, errorMeasure);
	}

}
