template<int dim>
ImplicitCollGalTsDesc<dim>::ImplicitCollGalTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitGappyTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod)
{
	podRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->getRestrictedDistInfo()));
	for (int iPod = 0; iPod < this->nPod; iPod++) { 
		this->restrictMapping()->restriction(this->pod[iPod], (*podRestrict)[iPod]);
		//(*AJRestrict)[iPod] = this->AJ[iPod];
	}
	this->jac.setNewSize(this->nPod,this->nPod);
	jactmp = new double [this->nPod * this->nPod];
}

template<int dim>
ImplicitCollGalTsDesc<dim>::~ImplicitCollGalTsDesc() 
{

}
template<int dim>
void ImplicitCollGalTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		projectVector(*(this->podRestrict), *(this->ResRestrict), From);	// different from PG
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

		transMatMatProd(*(this->podRestrict),*(this->AJRestrict),this->jactmp);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
				this->jac[iRow][iCol] = this->jactmp[iRow + iCol * this->pod.numVectors()];
			}
		} 

    solveLinearSystem(it, rhs, this->dUrom);
}
