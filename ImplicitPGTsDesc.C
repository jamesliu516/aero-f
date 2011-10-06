//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  residualRef(this->F),
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {

  parallelRom = new ParallelRom<dim>(*dom,this->com);
  parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);  
  lsCoeff = new double*[1];
  lsCoeff[0] = new double[this->nPod];
	this->projVectorTmp = new double [this->nPod];

	lsSolver = this->ioData->rom.lsSolver;
	if (lsSolver == 1){  // normal equations
		jactmp = new double [this->nPod * this->nPod];
		this->jac.setNewSize(this->nPod,this->nPod);
	}
}

//------------------------------------------------------------------------------
template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc(){

    delete [] lsCoeff[0];
    delete [] lsCoeff;  
		if (this->projVectorTmp) delete [] this->projVectorTmp;
		if (jactmp) delete [] jactmp;

}
//-----------------------------------------------------------------------------
template<int dim>
void ImplicitPGTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

	projectVector(this->AJ, this->F, From);
	Vec<double> rhs(this->nPod);
	rhs = -1.0 * From;

	// KTC FIX!
	// saving residual vectors (for GappyPOD)
	//writeBinaryVectorsToDisk1(false, it, 0.0, this->F, Dummy);

	res = rhs*rhs;

	if (res < 0.0){
		fprintf(stderr, "*** negative residual: %e\n", res);
		exit(1);
	}
	res = sqrt(res);

	if (it == 0) {
		this->target = this->epsNewton*res;
		this->res0 = res;
	}

	if (res == 0.0 || res <= this->target) {
		breakloop = true;
		return;	// do not solve the system
	}

	if (lsSolver == 0){	// normal equations

		transMatMatProd(this->AJ,this->AJ,jactmp);	// TODO: make symmetric product
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
			}
		} 

		solveLinearSystem(it, rhs, this->dUrom);
	}

	else if (lsSolver == 1)	{// ScaLAPACK least-squares

		 RefVec<DistSVec<double, dim> > residualRef2(this->F);
		 parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoeff);
		 for (int iPod=0; iPod<this->nPod; ++iPod)
		 this->dUrom[iPod] = -lsCoeff[0][iPod];
	}
}
