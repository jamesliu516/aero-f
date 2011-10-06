//------------------------------------------------------------------------------

template<int dim>
ImplicitBroydenTsDesc<dim>::ImplicitBroydenTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), dFrom(this->nPod) , rhs(this->nPod){
  
  JacSkipNewton = ioData.ts.implicit.newton.JacSkip;
  Fromold.resize(this->nPod);

	this->projVectorTmp = new double [this->nPod];
	jactmp = new double [this->nPod * this->nPod];
	this->jac.setNewSize(this->nPod,this->nPod);
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitBroydenTsDesc<dim>::~ImplicitBroydenTsDesc() 
{
	if (this->projVectorTmp) delete [] this->projVectorTmp;
	if (jactmp) delete [] jactmp;
}
//------------------------------------------------------------------------------

template<int dim>
void ImplicitBroydenTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		if (it==0 && (totTimeSteps % JacSkipNewton)==0) {
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
			breakloop = true; return;	// do not solve the system
		}

		// form reduced Jacobian

		// TODO KTC: computing AJ[iRom]AJ[iCol] requires communication after every
		// vec-vec product. We could express 
		transMatMatProd(this->AJ,this->AJ,jactmp);	// TODO: make symmetric product
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
			}
		} 

    solveLinearSystem(it, rhs, this->dUrom);
}

//------------------------------------------------------------------------------ 

template<int dim>
void ImplicitBroydenTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// do not compute each time

	if (it==0 && (totTimeSteps % JacSkipNewton)==0) {
		this->com->fprintf(stderr," ... Computing exact reduced Jacobian \n");
		ImplicitRomTsDesc<dim>::computeAJ(it, Q);
	}

}

//------------------------------------------------------------------------------ 

template<int dim>
void ImplicitBroydenTsDesc<dim>::broydenUpdate(Vec<double> &dFrom) {

  // Broyden update of the form B{k+1}=Bk+(yk-Bk*sk)*sk^T/(sk^T*sk) where yk is change in function and sk is step size

  Vec<double> zrom(this->nPod);
  zrom = dFrom;	// zrom = (yk-Bk*sk)
  for (int iPod = 0; iPod < this->nPod; ++iPod) { // KTC: parallelize rows of Bk
    for (int jPod = 0; jPod < this->nPod; ++jPod)
     zrom[iPod] -= this->jac[iPod][jPod]*this->dUrom[jPod];
  }

  double invNormSq = 1.0/ (this->dUrom*this->dUrom);
  zrom *= invNormSq;
  
  for (int iPod = 0; iPod < this->nPod; ++iPod) { // KTC: parallelize rows of zrom
    for (int jPod = 0; jPod < this->nPod; ++jPod)
      this->jac[iPod][jPod] += zrom[iPod]*this->dUrom[jPod];
  }


}
