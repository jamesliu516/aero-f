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

  A_Uinlet = NULL;    
  PhiT_A_Uinlet = NULL;
  PhiT_A_U = NULL;
  PhiT_A_Phi = NULL;
  A_Phi = NULL;
  regWeightConstant = 0.0;
  Kc = this->ioData->romOnline.constantGain;
  regWeightProportional = 0.0;
  Kp = this->ioData->romOnline.proportionalGain; 
  regWeightIntegral = 0.0;
  dRegWeightIntegral = 0.0;
  Ki = this->ioData->romOnline.integralGain;
  Ki_leak = -1.0 * this->ioData->romOnline.integralLeakGain;
  dKi = 0.0;
  ffError = 0.0;
  ffErrorPrev = 0.0;
  ffErrorTol = this->ioData->romOnline.ffErrorTol;
  
  controlNodeGlobalID = this->ioData->romOnline.controlNodeID - 1;
  controlNodeLocalID = -1;
  controlNodeSubDomain = -1;
  controlNodeCpuNum = -1;

  if (lsSolver == 2) {
    if (controlNodeGlobalID < 0) {
      fprintf(stderr, "*** Error: ControlNodeID must be specified\n");
      exit(-1);
    }    
    int numLocSub = this->domain->getNumLocSub();
    //bool abortOmp = false;
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      #pragma omp flush (abortOmp)
      //if (!abortOmp) {
        controlNodeLocalID = this->domain->getSubDomain()[iSub]->getLocalNodeNum(controlNodeGlobalID);
        if (controlNodeLocalID >= 0) {
          controlNodeSubDomain = iSub;
          controlNodeCpuNum = this->com->cpuNum();
          //abortOmp = true;
          //#pragma omp flush (abortOmp)
          break;
        }
      //}
    }
    this->com->globalMax(1,&controlNodeCpuNum);
  }

}

//------------------------------------------------------------------------------
template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc(){

    delete [] lsCoeff[0];
    delete [] lsCoeff;  
		if (this->projVectorTmp) delete [] this->projVectorTmp;
		if (jactmp) delete [] jactmp;
    if (pc) delete pc;
    if (A_Uinlet) delete A_Uinlet;
    if (PhiT_A_Uinlet) delete PhiT_A_Uinlet;
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

	res = rhs*rhs;
  double resReg = 0;

	if (res < 0.0){
		this->com->fprintf(stderr, "*** negative residual: %e\n", res);
		exit(1);
	}
	res = sqrt(res);

  if (minRes<=0) minRes = res / this->restart->residual; 

	if (it == 0) {
		this->target = this->epsNewton*res;
		this->res0 = res;
	}

	if (res == 0.0 || res <= this->target) {
		this->com->fprintf(stdout, "breakloop=true: res=%e\n", res);
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

    //if (it==0) minRes = ((res/this->restart->residual)<minRes) ? res/this->restart->residual : minRes;
 
    ffError = 0.0;

    if (this->com->cpuNum() == controlNodeCpuNum) {
      double (*u)[dim] = U.subData(controlNodeSubDomain); // vector of control volumes
      for (int j=0; j<dim; ++j) {
         ffError += pow(u[controlNodeLocalID][j] - (this->bcData->getInletConservativeState())[j], 2);
      }
    }

    this->com->globalMax(1,&ffError);

    this->com->fprintf(stdout, " ... Far-field error at control node %d = %e \n", controlNodeGlobalID+1, ffError);

    regWeightConstant = Kc;
    regWeightProportional = (ffError>ffErrorTol) ? Kp*ffError : 0.0 ;
    Ki_leak = (dRegWeightIntegral<0) ? Ki_leak*2.0 : -1.0*(this->ioData->romOnline.integralLeakGain);
    dRegWeightIntegral = (ffError>ffErrorTol) ? Ki*ffError : Ki_leak*ffErrorTol;
    regWeightIntegral += dRegWeightIntegral;
    regWeightIntegral = (regWeightIntegral>0) ? regWeightIntegral : 0;

    this->regWeight = regWeightConstant + regWeightProportional + regWeightIntegral;
    this->regWeight = (this->regWeight>0) ? this->regWeight : 0;

    ffErrorPrev = ffError;

    this->com->fprintf(stdout, " ... Regularization Weighting = %e (P = %e, I = %e)\n",
                       this->regWeight, regWeightProportional, regWeightIntegral); 

    // forming (PhiT * A) * U
    if (PhiT_A_U) delete PhiT_A_U;
    PhiT_A_U = new Vec<double>(this->nPod);
    for (int iVec = 0; iVec < this->nPod; iVec++){
      (*PhiT_A_U)[iVec] = (*A_Phi)[iVec] * U;
    }

    rhs = rhs + (this->regWeight * (*PhiT_A_Uinlet - *PhiT_A_U));
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


   // begin: output convergence information
    this->com->fprintf(stdout, " ... forming A * (U-Uinlet).^2\n");
    DistSVec<double, dim>* A_Uerr = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *A_Uerr = 0.0;
    int numLocSub = A_Uerr->numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = this->A->subData(iSub); // vector of control volumes
      double (*auerr)[dim] = A_Uerr->subData(iSub);
      double (*u)[dim] = U.subData(iSub);
      for (int i=0; i<this->A->subSize(iSub); ++i) {
        if (cv[i]>this->regThresh) {
          for (int j=0; j<dim; ++j)
            auerr[i][j] = cv[i] * pow(u[i][j] - (this->bcData->getInletConservativeState())[j],2);
        }
      }
    }

    DistSVec<double, dim>* ones = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *ones = 1.0;

    this->com->fprintf(stdout, "(||R||)^2 %e\n", (this->F)*(this->F));
    this->com->fprintf(stdout, "(||U - Uinlet||_A)^2 = %e\n", (*A_Uerr) * (*ones));
    this->com->fprintf(stdout, "||(J*Phi)' * R || = %e\n", res);
    this->com->fprintf(stdout, "||(J*Phi)' * R + reg|| = %e\n\n", resReg);

    if (A_Uerr) delete A_Uerr; 
    if (ones) delete ones;
   // end: output convergence information

    res = resReg;

    transMatMatProd(this->AJ,this->AJ,jactmp);  // TODO: make symmetric product
    for (int iRow = 0; iRow < this->nPod; ++iRow) {
      for (int iCol = 0; iCol < this->nPod; ++iCol) {
        this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod] + (this->regWeight*(*PhiT_A_Phi)[iRow][iCol]);
        if (iRow==iCol && false) this->jac[iRow][iCol] = this->jac[iRow][iCol] + 1.0/dt;
      }
    }
 
    this->solveLinearSystem(it, rhs, this->dUromNewtonIt);    

  } 
  
}


//-----------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::applyWeightingToLeastSquaresSystem() {

  // weight the least-squares system

  DistSVec<double, dim> weightVec(this->domain->getNodeDistInfo());
  int numLocSub = this->domain->getNumLocSub();
 
  switch (this->ioData->romOnline.weightedLeastSquares) {
    case (NonlinearRomOnlineData::WEIGHTED_LS_FALSE):
      return;
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_RESIDUAL):
      weightVec = *(this->weightFRef);
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_STATE):
      weightVec = *(this->weightURef);  
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double (*weight)[dim] = weightVec.subData(iSub);
        for (int i=0; i<weightVec.subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = weight[i][j] - (this->bcData->getInletConservativeState())[j];
        }
      }
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_CV):
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*weight)[dim] = weightVec.subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = cv[i];
        }
      }
      break;
    default:
        this->com->fprintf(stderr, "*** Error: Unexpected least-squares weighting method\n");
        exit(-1);
      break;
  } 

  double weightExp = this->ioData->romOnline.weightingExponent;
  double weightNorm = weightVec.norm();
  weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;

  // weight residual
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*weight)[dim] = weightVec.subData(iSub);
    double (*f)[dim] = this->F.subData(iSub);
    for (int i=0; i<weightVec.subSize(iSub); ++i) {
      for (int j=0; j<dim; ++j) {
        weight[i][j] = pow(abs(weight[i][j])/weightNorm, weightExp);
        f[i][j] = f[i][j] * weight[i][j];
      }
    }
  }
    
  // weight AJ
  for (int iVec=0; iVec<this->nPod; ++iVec) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double (*weight)[dim] = weightVec.subData(iSub);
      double (*aj)[dim] = (this->AJ[iVec]).subData(iSub);
      for (int i=0; i<weightVec.subSize(iSub); ++i) {
        for (int j=0; j<dim; ++j)
          aj[i][j] = aj[i][j] * weight[i][j];
      }
    }
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
    if (A_Uinlet==NULL) {  // only needs to be done once
      this->com->fprintf(stdout, " ... forming A * Uinlet\n");
      A_Uinlet = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
      *A_Uinlet = 0.0;
      int numLocSub = A_Uinlet->numLocSub();
      double maxCV = 0;
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*au)[dim] = A_Uinlet->subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          if (cv[i]>this->regThresh) {
            for (int j=0; j<dim; ++j)
              au[i][j] = cv[i] * (this->bcData->getInletConservativeState())[j];
          }
          if (cv[i] > maxCV) maxCV = cv[i];
        }
      }
      this->com->globalMax(1,&maxCV); 
      this->com->fprintf(stdout, " ... regularization CV threshold = %e (largest CV = %e) \n", this->regThresh, maxCV);
    }

   // for (int j=1; j<dim; ++j)
   //    this->com->fprintf(stdout, " ... InletConservativeState = %e \n", (this->bcData->getInletConservativeState())[j]);

    this->com->fprintf(stdout, " ... forming PhiT * (A * Uinlet)\n");
    if (PhiT_A_Uinlet) delete PhiT_A_Uinlet;
    PhiT_A_Uinlet = new Vec<double>(this->nPod);
    for (int iVec = 0; iVec < this->nPod; iVec++){
      (*PhiT_A_Uinlet)[iVec] = this->pod[iVec] * (*A_Uinlet);
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
          if (cv[i]>this->regThresh) {
            for (int j=0; j<dim; ++j)
              au[i][j] *= cv[i];
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


