template<int dim>
ImplicitCollLSTsDesc<dim>::ImplicitCollLSTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod)
{
	this->jac.setNewSize(this->nPod,this->nPod);
	jactmp = new double [this->nPod * this->nPod];
	this->projVectorTmp = new double [this->nPod];

	// read in sample nodes
	nSampleNodes = 0;
	dom->readSampleNodes(sampleNodes, nSampleNodes, this->input->sampleNodes);

	// determine mapping to restricted nodes
	restrictionMapping.reset(new RestrictionMapping<dim>(dom, sampleNodes.begin(), sampleNodes.end()));

	AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, getRestrictedDistInfo()));
	ResRestrict.reset(new DistSVec<double, dim> (getRestrictedDistInfo()));
}

template<int dim>
ImplicitCollLSTsDesc<dim>::~ImplicitCollLSTsDesc() 
{
	if (this->projVectorTmp) delete [] this->projVectorTmp;
}


//------------------------------------------------------------------------------

template<int dim>
void ImplicitCollLSTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q) {

	// Evaluate residual on full mesh

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, *restrictionMapping);

	this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Q,
			this->F, restrictionMapping->getRestrictedToOriginLocNode());

  this->spaceOp->applyBCsToResidual(Q, this->F);

	double t0 = this->timer->getTime();

	restrictMapping()->restriction(this->F, *ResRestrict);

	this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCollLSTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

	// Evaluate action of Jacobian on full mesh

	this->mvpfd->evaluateRestrict(it, *this->X, *this->A, Q, this->F,
			*restrictionMapping);	// very cheap
  
  for (int iPod = 0; iPod < this->nPod; iPod++) {
		this->mvpfd->applyRestrict(this->pod[iPod], this->AJ[iPod],
				*restrictionMapping);
	}

	double t0 = this->timer->getTime();
	for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
		restrictMapping()->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
	}
	this->timer->addRestrictionTime(t0);
}

template<int dim>
void ImplicitCollLSTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		this->projectVector(*AJRestrict, *this->ResRestrict, From);	// different from PG
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

		transMatMatProd(*AJRestrict,*AJRestrict,jactmp);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->pod.numVectors()];
			}
		} 

    this->solveLinearSystem(it, rhs, this->dUrom);
}

template<int dim>
bool ImplicitCollLSTsDesc<dim>::breakloop1(const bool breakloop) {

	return false;
	
}

template<int dim>
bool ImplicitCollLSTsDesc<dim>::breakloop2(const bool breakloop) {

	return breakloop;

}
