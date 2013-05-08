//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  residualRef(this->F),
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {

  pc = ImplicitRomTsDesc<dim>::template 
  createPreconditioner<PrecScalar,dim>(this->ioData->ts.implicit.newton.ksp.ns.pc, this->domain);

  parallelRom = new ParallelRom<dim>(*dom,this->com);

  // TODO necessary?
  currentProblemSize = this->nPod;
  parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);  
  lsCoeff = new double*[1];
  lsCoeff[0] = new double[this->nPod];
	this->projVectorTmp = new double [this->nPod];

	lsSolver = this->ioData->romOnline.lsSolver;
	if ((lsSolver == 1) || (lsSolver == 2)){  // normal equations
		jactmp = new double [this->nPod * this->nPod];
		this->jac.setNewSize(this->nPod,this->nPod);
	}

  minRes = -1;
  rhsNormInit = -1;

  this->rom->initializeClusteredOutputs();

  dtInit = this->ioData->romOnline.reducedTimeStep;
  dt = dtInit;

  A_Uinit = NULL;    
  PhiT_A_Uinit = NULL;
  PhiT_A_U = NULL;
  PhiT_A_Phi = NULL;
  A_Phi = NULL;
  regCoeff = this->ioData->romOnline.regCoeff;
  regThresh = this->ioData->romOnline.regThresh;
}

//------------------------------------------------------------------------------
template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc(){

    delete [] lsCoeff[0];
    delete [] lsCoeff;  
		if (this->projVectorTmp) delete [] this->projVectorTmp;
		if (jactmp) delete [] jactmp;
    if (pc) delete pc;
    if (A_Uinit) delete A_Uinit;
    if (PhiT_A_Uinit) delete PhiT_A_Uinit;
    if (PhiT_A_U) delete PhiT_A_U;
    if (PhiT_A_Phi) delete PhiT_A_Phi;
    if (A_Phi) delete A_Phi;
}

//-----------------------------------------------------------------------------
template<int dim>
void ImplicitPGTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

	this->projectVector(this->AJ, this->F, From);
	Vec<double> rhs(this->nPod);
	rhs = -1.0 * From;

	// KTC FIX!
	// saving residual vectors (for GappyPOD)
	//writeBinaryVectorsToDisk1(false, it, 0.0, this->F, Dummy);

	res = rhs*rhs;
  double resReg = 0;

	if (res < 0.0){
		fprintf(stderr, "*** negative residual: %e\n", res);
		exit(1);
	}
	res = sqrt(res);

  if (minRes<=0) minRes = res / this->restart->residual;
  

	if (it == 0) {
		this->target = this->epsNewton*res;
		this->res0 = res;
	}

	if (res == 0.0 || res <= this->target) {
		breakloop = true;
		return;	// do not solve the system
	}

	if (lsSolver == 0){	// ScaLAPACK least-squares

		RefVec<DistSVec<double, dim> > residualRef2(this->F);
		parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoeff);
		for (int iPod=0; iPod<this->nPod; ++iPod)
			this->dUromNewtonIt[iPod] = -lsCoeff[0][iPod];
	}
	else if (lsSolver == 1)	{		// normal equations
		transMatMatProd(this->AJ,this->AJ,jactmp);	// TODO: make symmetric product
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
			}
		} 
		this->solveLinearSystem(it, rhs, this->dUromNewtonIt);
	} 
  else if (lsSolver == 2) {  // regularized normal equations

    if (it==0) minRes = ((res/this->restart->residual)<minRes) ? res/this->restart->residual : minRes;
    double scaling =  minRes * minRes / ( 1.0 + double(totalTimeSteps) );

    this->com->fprintf(stdout, "scaling = %e / %d = %e\n", minRes, totalTimeSteps, scaling); 

    // forming (PhiT * A) * U
    if (PhiT_A_U) delete PhiT_A_U;
    PhiT_A_U = new Vec<double>(this->nPod);
    for (int iVec = 0; iVec < this->nPod; iVec++){
      (*PhiT_A_U)[iVec] = (*A_Phi)[iVec] * U;
    }

    rhs = rhs + (scaling * (*PhiT_A_Uinit - *PhiT_A_U));
    if (rhsNormInit<0) rhsNormInit = rhs.norm();
    if (it==0) {
      double rhsNormPrev = rhs.norm();
      dt = dtInit * (rhsNormInit/rhsNormPrev) ;  
    }
    resReg = sqrt(rhs*rhs);
    Vec<double> dUrom(this->dUromTimeIt);
    dUrom /= dt;
    if (false)
      rhs -= dUrom; 

    transMatMatProd(this->AJ,this->AJ,jactmp);  // TODO: make symmetric product
    for (int iRow = 0; iRow < this->nPod; ++iRow) {
      for (int iCol = 0; iCol < this->nPod; ++iCol) {
        this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod] + (scaling*(*PhiT_A_Phi)[iRow][iCol]);
        if (iRow==iCol && false) this->jac[iRow][iCol] = this->jac[iRow][iCol] + 1.0/dt;
      }
    }
 
    this->solveLinearSystem(it, rhs, this->dUromNewtonIt);    

    this->com->fprintf(stdout, "||(J*Phi)' * R || = %e\n", res);
    this->com->fprintf(stdout, "||(J*Phi)' * R + reg|| = %e\n\n", resReg);
 
    res = resReg;

  }


 


}

//-----------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {

  if (currentProblemSize != this->nPod){
    currentProblemSize = this->nPod;

    parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);

    if (lsCoeff) delete [] lsCoeff[0];
    lsCoeff[0] = new double[this->nPod];

    if (this->projVectorTmp) delete (this->projVectorTmp);
    this->projVectorTmp = new double [this->nPod];

    if ((lsSolver == 1) || (lsSolver == 2)) {
      if (jactmp) delete [] jactmp;
      jactmp = new double [this->nPod * this->nPod];
      this->jac.setNewSize(this->nPod,this->nPod);
    }

    From.resize(this->nPod);
    rhs.resize(this->nPod);
  }

  if (lsSolver == 2) {
    if (A_Uinit==NULL) {  // only needs to be done once
      this->com->fprintf(stdout, " ... forming A * Uinit\n");
      A_Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
      *A_Uinit = U;
      int numLocSub = A_Uinit->numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*au)[dim] = A_Uinit->subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          if (cv[i]>regThresh) {
            for (int j=0; j<dim; ++j)
              au[i][j] *= (regCoeff * cv[i]);
          } else {
            for (int j=0; j<dim; ++j)
              au[i][j] *= 0.0;
          }
        }
      }
    }

    this->com->fprintf(stdout, " ... forming PhiT * (A * Uinit)\n");
    if (PhiT_A_Uinit) delete PhiT_A_Uinit;
    PhiT_A_Uinit = new Vec<double>(this->nPod);
    for (int iVec = 0; iVec < this->nPod; iVec++){
      (*PhiT_A_Uinit)[iVec] = this->pod[iVec] * (*A_Uinit);
    }

    // form PhiT*A (because we need to compute PhiT*A*U at every timestep)

    this->com->fprintf(stdout, " ... forming PhiT * A * Phi\n");

    if (A_Phi) delete A_Phi;
    A_Phi = new VecSet< DistSVec<double, dim> >(this->nPod, this->domain->getNodeDistInfo());
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      (*A_Phi)[iVec] = this->pod[iVec];
      int numLocSub = ((*A_Phi)[iVec]).numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) { 
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*au)[dim] = ((*A_Phi)[iVec]).subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          if (cv[i]>regThresh) {
            for (int j=0; j<dim; ++j)
              au[i][j] *= (regCoeff * cv[i]);
          } else {
            for (int j=0; j<dim; ++j)
              au[i][j] *= 0.0;
          }
        }
      }
    }
    if (PhiT_A_Phi) delete PhiT_A_Phi;
    PhiT_A_Phi = new VecSet< Vec<double> >(this->nPod, this->nPod);

    Vec<double> tmpVec(this->nPod);
    for (int iVec = 0; iVec < this->nPod; ++iVec) {
      for (int jVec = 0; jVec < this->nPod; ++jVec){
        tmpVec[jVec] = this->pod[jVec] * (*A_Phi)[iVec];
      }
      (*PhiT_A_Phi)[iVec] = tmpVec;
    }

  }

}


