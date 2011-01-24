//------------------------------------------------------------------------------

template<int dim>
ImplicitGalerkinTsDesc<dim>::ImplicitGalerkinTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {
  
}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitGalerkinTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		projectVector(this->pod, this->F, From);	// different from PG
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

		// form reduced Jacobian

		this->jac.setNewSize(this->nPod,this->nPod);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol <= iRow; ++iCol) {
				this->jac[iRow][iCol] = this->pod[iRow]*this->AJ[iCol];	// different from PG
				if (iRow > iCol)
					this->jac[iCol][iRow] = this->jac[iRow][iCol];
			}
		} 

    solveLinearSystem(it, rhs, this->dUrom);
}
