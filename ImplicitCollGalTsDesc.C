template<int dim>
ImplicitCollGalTsDesc<dim>::ImplicitCollGalTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitCollLSTsDesc<dim>(ioData, geoSource, dom)
{
	podRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->getRestrictedDistInfo()));
	for (int iPod = 0; iPod < this->nPod; iPod++) { 
		this->restrictMapping()->restriction(this->pod[iPod], (*podRestrict)[iPod]);
		//(*AJRestrict)[iPod] = this->AJ[iPod];
	}
}

template<int dim>
void ImplicitCollGalTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		this->projectVector(*podRestrict, *(this->ResRestrict), this->From);	// different from PG
		this->rhs = -1.0 * this->From;

    res = this->rhs*this->rhs;

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

		transMatMatProd(*podRestrict,*(this->AJRestrict),this->jactmp);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
				this->jac[iRow][iCol] = this->jactmp[iRow + iCol * this->pod.numVectors()];
			}
		} 

    this->solveLinearSystem(it, this->rhs, this->dUrom);
}
