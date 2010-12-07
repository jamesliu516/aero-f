//------------------------------------------------------------------------------

template<int dim>
ImplicitBroydenTsDesc<dim>::ImplicitBroydenTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), dFrom(this->nPod) , rhs(this->nPod){
  
  JacSkipNewton = ioData.ts.implicit.newton.JacSkip;
  Fromold.resize(this->nPod);
}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitBroydenTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		if (it==0 && (Git % JacSkipNewton)==0) {
			this->com->fprintf(stderr," ... Computing exact reduced Jacobian \n");
			projectVector(this->AJ, this->F, From);
		}
		else {
			projectVector(this->AJ, this->F, From);
			if (it > 0) {
				dFrom = From-Fromold;
				broydenUpdate(dFrom);
			}
		}
		Fromold = From;
		rhs = -1.0 * From;

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
				this->jac[iRow][iCol] = this->AJ[iRow]*this->AJ[iCol];
				if (iRow > iCol)
					this->jac[iCol][iRow] = this->jac[iRow][iCol];
			}
		} 

    solveLinearSystem(it, rhs, this->dUrom);
}

//------------------------------------------------------------------------------ 

template<int dim>
void ImplicitBroydenTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// do not compute each time

	if (it==0 && (Git % JacSkipNewton)==0) {
		ImplicitRomTsDesc<dim>::computeAJ(it, Q);
	}

}

//------------------------------------------------------------------------------ 

template<int dim>
void ImplicitBroydenTsDesc<dim>::broydenUpdate(Vec<double> &dFrom) {

  // Broyden update of the form B{k+1}=Bk+(yk-Bk*sk)*sk^T/(sk^T*sk) where yk is change in function and sk is step size

  Vec<double> zrom(this->nPod);
  zrom = dFrom;	// zrom = (yk-Bk*sk)
  for (int iPod = 0; iPod < this->nPod; ++iPod) { // KTC: parallelize
    for (int jPod = 0; jPod < this->nPod; ++jPod)
     zrom[iPod] -= this->jac[iPod][jPod]*this->dUrom[jPod];
  }

  double invNormSq = 1.0/ (this->dUrom*this->dUrom);
  zrom *= invNormSq;
  
  for (int iPod = 0; iPod < this->nPod; ++iPod) { // KTC: parallelize
    for (int jPod = 0; jPod < this->nPod; ++jPod)
      this->jac[iPod][jPod] += zrom[iPod]*this->dUrom[jPod];
  }


}
