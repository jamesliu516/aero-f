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
	} else {
    jactmp = NULL;
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

  this->com->fprintf(stdout, "||W*R(U)|| = %1.12e\n", this->F.norm());
  this->com->fprintf(stdout, "|| U || = %1.12e\n", U.norm());
  DistSVec<double, dim> romContribution(this->domain->getNodeDistInfo());
  if (this->Uinit) {
    romContribution = U - *(this->Uinit); 
    this->com->fprintf(stdout, "|| U - Uic || = %1.12e\n", romContribution.norm());
  }

  this->com->fprintf(stdout, "||(W*J*Phi)' * W*R|| = %1.12e\n", res);

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

	if (res == 0.0 || res <= this->target || res <= this->epsAbsResNewton) {
		//this->com->fprintf(stdout, "breakloop=true: res=%e\n", res);
		breakloop = true;
		return;	// do not solve the system
	}

	if (lsSolver == 0){	// ScaLAPACK least-squares

		RefVec<DistSVec<double, dim> > residualRef2(this->F);
		parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoeff);
    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
			this->dUromNewtonIt[iPod] = -lsCoeff[0][iPod];
      //this->com->fprintf(stdout, " ... dUromNewtonIt[%d] = %e \n", iPod, this->dUromNewtonIt[iPod]);
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
    }
    this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %1.12e \n", dUromNewtonItNormSquared);

	} else if (lsSolver == 3){	// Solve Levenberg-Marquardt regularized LS via ScaLAPACK SVD

    // SVD quantities
    VecSet< DistSVec<double, dim> >* U_AJ = new VecSet< DistSVec<double, dim> >(this->nPod, this->domain->getNodeDistInfo());
    Vec<double>* sVals_AJ = new Vec<double>(this->nPod);
    FullM* V_AJ = new FullM(this->nPod);

    parallelRom->parallelSVD(this->AJ, *U_AJ, sVals_AJ->data(), *V_AJ, this->nPod, true);
    // Note: V, not V_transpose

   /* double maxErr = 0.0;
    double avgErr = 0.0;
    for (int iPod = 0; iPod < this->nPod; ++iPod) {
      DistSVec<double, dim> error(this->domain->getNodeDistInfo());
      error  = this->AJ[iPod];
      for (int jPod = 0; jPod < this->nPod; ++jPod)
        error = error - (((*singVals)[jPod]*(*Vtrue)[iPod][jPod])*(*Utrue)[jPod]);
      double errorNorm = error.norm()/((this->AJ[iPod]).norm());
      avgErr += errorNorm;
      if (errorNorm > maxErr)
        maxErr = errorNorm;   
    }
    avgErr /= this->nPod;
  
    this->com->fprintf(stderr, " ... Average error on AJ after SVD = %e\n", avgErr);  
    this->com->fprintf(stderr, " ... Maximum error on AJ after SVD = %e\n", maxErr); */
    
    Vec<double> tmpVec(this->nPod);
    for (int iVec=0; iVec<U_AJ->numVectors(); ++iVec)
      tmpVec[iVec] = (*U_AJ)[iVec] * (-1.0*this->F);

    delete U_AJ;

    if ((*sVals_AJ)[this->nPod-1]>0) {
      this->com->fprintf(stdout, " ... Singular value ratio for AJ = %e \n", (*sVals_AJ)[0]/(*sVals_AJ)[this->nPod-1]);
    } else {
      this->com->fprintf(stdout, " ... AJ is rank deficient! \n");
    }

    double lambdaSquared = pow(this->levenbergMarquardtWeight,2);

    double firstTermNormSquared = 0.0;
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      if ((*sVals_AJ)[iVec]>0) {
        firstTermNormSquared += pow(tmpVec[iVec],2)*pow(lambdaSquared/(pow((*sVals_AJ)[iVec],2) + lambdaSquared),2);
      } else {
        firstTermNormSquared += pow(tmpVec[iVec],2);
      }
    }

    double secondTermNormSquared = 0.0;
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      tmpVec[iVec] = ((*sVals_AJ)[iVec]>0) ? tmpVec[iVec]*(*sVals_AJ)[iVec]/(pow((*sVals_AJ)[iVec],2) + lambdaSquared) : 0.0;
      secondTermNormSquared += pow(tmpVec[iVec],2);
    }

    delete sVals_AJ;

    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      this->dUromNewtonIt[iPod] = 0.0;
      for (int jPod=0; jPod<this->nPod; ++jPod) {
        this->dUromNewtonIt[iPod] += (*V_AJ)[iPod][jPod] * tmpVec[jPod];
      }
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
      //this->com->fprintf(stdout, " ... dUromNewtonIt[%d] = %e \n", iPod, this->dUromNewtonIt[iPod]);
    }
    this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %e \n", dUromNewtonItNormSquared);
    this->com->fprintf(stdout, " ... || Ax - b ||^2 = %1.12e \n", firstTermNormSquared);
    this->com->fprintf(stdout, " ... || x ||^2 = %1.12e \n", secondTermNormSquared);
 
    delete V_AJ;

  } else if (lsSolver == 1)	{		// normal equations
		transMatMatProd(this->AJ,this->AJ,jactmp);	// TODO: make symmetric product
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
			}
		} 
		this->solveLinearSystem(it, rhs, this->dUromNewtonIt);

    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      //this->com->fprintf(stdout, " ... dUromNewtonIt[%d] = %e \n", iPod, this->dUromNewtonIt[iPod]);
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
    }
    this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %e \n", dUromNewtonItNormSquared);

	} else if (lsSolver == 2) {  // regularized normal equations

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
void ImplicitPGTsDesc<dim>::updateLeastSquaresWeightingVector() {

  // form the weighting vector for the least-squares system

  int numLocSub = this->domain->getNumLocSub();
 
  switch (this->ioData->romOnline.weightedLeastSquares) {
    case (NonlinearRomOnlineData::WEIGHTED_LS_FALSE):
      return;
      break;
/*    case (NonlinearRomOnlineData::WEIGHTED_LS_RESIDUAL):
      *(this->weightVec) = *(this->weightFRef);
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_STATE):
      *(this->weightVec) = *(this->weightURef);  
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double (*weight)[dim] = this->weightVec->subData(iSub);
        for (int i=0; i<this->weightVec->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = weight[i][j] - (this->bcData->getInletConservativeState())[j];
        }
      }
      break;*/
    case (NonlinearRomOnlineData::WEIGHTED_LS_CV):
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = this->A->subData(iSub); // vector of control volumes
        double (*weight)[dim] = this->weightVec->subData(iSub);
        for (int i=0; i<this->A->subSize(iSub); ++i) {
          for (int j=0; j<dim; ++j)
            weight[i][j] = cv[i];
        }
      }
      break;
    case (NonlinearRomOnlineData::WEIGHTED_LS_BOCOS):
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *ffMask = this->farFieldMask->subData(iSub); // vector with nonzero entries at farfield nodes
        double (*weight)[dim] = this->weightVec->subData(iSub);
        for (int i=0; i<this->farFieldMask->subSize(iSub); ++i) {
          if (ffMask[i]>0) {
            for (int j=0; j<dim; ++j)
              weight[i][j] = this->ffWeight;
          } else {
            for (int j=0; j<dim; ++j)
               weight[i][j] = 1.0;
          }
        }
      }
      break;
    default:
        this->com->fprintf(stderr, "*** Error: Unexpected least-squares weighting method\n");
        exit(-1);
      break;
  } 
/*
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
*/


}

//-----------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {

  if (currentProblemSize != this->nPod){
    currentProblemSize = this->nPod;

    parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);

    if (lsCoeff) delete [] lsCoeff[0];
    lsCoeff[0] = new double[this->nPod];

    if (this->projVectorTmp) delete [] (this->projVectorTmp);
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

//------------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));

  if (!this->problemType[ProblemData::UNSTEADY]) {
    this->com->printf(2, "\nNOTE: For weighted ROM simulations the reported residual is calculated using a weighted norm,\n");
    this->com->printf(2, "      and is relative to the residual of the initial condition calculated using the same norm.\n");

    this->Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *(this->Uinit) = U;  // needed for computing the restricted residual after each cluster switch
  }

  this->com->printf(2, "\n");
  
}   

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitPGTsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if (!this->problemType[ProblemData::UNSTEADY] && monitorConvergence(it, U))
    return true;

  if (!this->problemType[ProblemData::AERO] && !this->problemType[ProblemData::THERMO] && it >= this->data->maxIts) return true;

  if (this->problemType[ProblemData::UNSTEADY] )
    if(t >= this->data->maxTime - 0.01 * dt)
      return true;

  return false;

}


//------------------------------------------------------------------------------
 
template<int dim>
bool ImplicitPGTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// only called for steady simulations

  this->data->residual = computePGResidualNorm(U);

  if (this->data->residual == 0.0 || this->data->residual < this->data->eps * this->restart->residual || this->data->residual < this->data->epsabs)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitPGTsDesc<dim>::computePGResidualNorm(DistSVec<double,dim>& Q)
{ // spatial only

  this->spaceOp->computeResidual(*this->X, *this->A, Q, this->F, this->timeState);

  this->spaceOp->applyBCsToResidual(Q, this->F);

  // weight residual
  if (this->ioData->romOnline.weightedLeastSquares != NonlinearRomOnlineData::WEIGHTED_LS_FALSE) {
    this->updateLeastSquaresWeightingVector();
    double weightExp = this->ioData->romOnline.weightingExponent;
    double weightNorm = this->weightVec->norm();
    weightNorm = (weightNorm<=0.0) ? 1.0 : weightNorm;
    int numLocSub = Q.numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double (*weight)[dim] = this->weightVec->subData(iSub);
      double (*f)[dim] = this->F.subData(iSub);
      for (int i=0; i<this->weightVec->subSize(iSub); ++i) {
        for (int j=0; j<dim; ++j) {
          //weight[i][j] = pow(abs(weight[i][j])/weightNorm, weightExp);
          f[i][j] = f[i][j] * weight[i][j];
        }
      }
    }
  }

  double res = 0.0;
  res = (this->F) * (this->F);
  return sqrt(res);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::setReferenceResidual()
{
  if (this->Uinit) this->restart->residual = computePGResidualNorm(*(this->Uinit));

  this->com->printf(2, "Norm of reference residual = %.12e\n", this->restart->residual);

}


