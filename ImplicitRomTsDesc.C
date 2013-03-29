#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <NewtonSolver.h>

#include <TsOutput.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::ImplicitRomTsDesc(IoData &_ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(_ioData, geoSource, dom), pod(0, dom->getNodeDistInfo()), F(dom->getNodeDistInfo()), AJ(0, dom->getNodeDistInfo()) {

	ioData = &_ioData;
  tag = 0;

  maxItsNewton = ioData->ts.implicit.newton.maxIts;
  epsNewton = ioData->ts.implicit.newton.eps;  

  //ns = new NewtonSolver<ImplicitRomTsDesc<dim> >(this);
  this->timeState = new DistTimeState<dim>(*ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;
  mvpfd = new MatVecProdFD<dim,dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, *ioData);

  // read Pod Basis
  nPod = ioData->rom.dimension;
	readPodBasis(this->input->podFileState);

	char *snapRefSolFile = this->input->snapRefSolutionFile;
	FILE *inRSFP = fopen(snapRefSolFile, "r");
	if (snapRefSolFile[0] != '\0'){  
		DistSVec<double,dim> referenceSolution(dom->getNodeDistInfo());
		this->com->fprintf(stderr, "Reading reference solution for snapshots in %s\n", snapRefSolFile);
		dom->readVectorFromFile(snapRefSolFile, 0, 0, referenceSolution);
		for (int iPod = 0; iPod < nPod; ++iPod) {
			pod[iPod] -= referenceSolution;
		}
	}

  MemoryPool mp;
  this->mmh = this->createMeshMotionHandler(*ioData, geoSource, &mp);

  AJ.resize(nPod);
  dUrom.resize(nPod);
  UromTotal.resize(nPod);
	UromTotal = 0.0;	// before any time step, it is zero
	dUnormAccum = new Vec<double> [2];
	for (int i = 0 ; i < 2; ++i) {
		dUnormAccum[i].resize(nPod);	// before any time step, it is zero
		dUnormAccum[i] = 0.0;	// before any time step, it is zero
	}

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::~ImplicitRomTsDesc()
{

	// output net dUnormAccum
	for (int iPod = 0; iPod < nPod; ++iPod)
		dUnormAccum[1][iPod] = sqrt(dUnormAccum[1][iPod]);	// to complete 2-norm definition

	const char *dUnormAccumFile = ioData->output.rom.dUnormAccum;
	if (strcmp(dUnormAccumFile, "") != 0)  {
		char *outdUnormAccumFile = new char[strlen(ioData->output.rom.prefix) +
			strlen(dUnormAccumFile)+1];
		if (this->com->cpuNum() ==0) sprintf(outdUnormAccumFile, "%s%s",
				ioData->output.rom.prefix, dUnormAccumFile);
		FILE *writingFile;
		if (this->com->cpuNum() ==0) writingFile = fopen(outdUnormAccumFile, "wt");
		this->com->fprintf(writingFile, "PodCoefficient NetContribution(1norm) NetContribution(2norm)\n");
		for (int i = 0; i < nPod; ++i) { 
			this->com->fprintf(writingFile, "%d %e %e \n", i + 1, dUnormAccum[0][i],dUnormAccum[1][i]);
		}
		delete [] outdUnormAccumFile;
	}

	delete [] dUnormAccum;
  if (tag) delete tag;
}


//------------------------------------------------------------------------------
template<int dim>
int ImplicitRomTsDesc<dim>::solveNonLinearSystem(DistSVec<double, dim> &U, const int totalTimeSteps)  {

	// initializations 

  double t0 = this->timer->getTime();

  int it = 0;
  int fsIt = 0;
	updateGlobalTimeSteps(totalTimeSteps);

  Vec<double> Urom(nPod); // reduced coordinates
  Urom = 0.0; // total solution increment in ROM coordinates (the unknowns for reduced problem)
  DistSVec<double, dim> dUfull(this->domain->getNodeDistInfo());	// solution increment at EACH NEWTON ITERATION in full coordinates
  dUfull = 0.0;	// initial zero increment

  double res;
  bool breakloop = false;
  bool breakloopNow = false;

	// line search variables
  double alpha;
  bool convergeFlag=0;

	postProStep(U,totalTimeSteps);

  for (it = 0; it < maxItsNewton; it++)  {

		double tRes = this->timer->getTime();
    computeFullResidual(it, U);
		computeAJ(it, U);	// skipped some times for Broyden
		this->timer->addResidualTime(tRes);

		solveNewtonSystem(it, res, breakloop);	// 1) check if residual small enough, 2) solve 
			// INPUTS: AJ, F
			// OUTPUTS: dUrom, res, breakloop
		breakloopNow = breakloop1(breakloop);
		if (breakloopNow) break;

// LINE SEARCH
//    // do line search (linesearch exits with alpha=0 and convergenceFlag if convergence criteria is satisfied)
//    alpha = lineSearch(U,dUrom,it,AJ,epsNewton, convergeFlag);
//    if (it > 0 && convergeFlag == 1) break;
//    dUrom *= alpha;
// END LINE SEARCH

		double tSol = this->timer->getTime();
    expandVector(dUrom, dUfull); // solution increment in full coordinates
    Urom += dUrom; // solution increment in reduced coordinates
    UromTotal += dUrom; // solution increment in reduced coordinates
    U += dUfull;
		this->timer->addSolutionIncrementTime(tSol);

		saveNewtonSystemVectors(totalTimeSteps);	// only implemeted for PG rom

    // verify that the solution is physical
    if (this->checkSolution(U)) {
      if (checkFailSafe(U) && fsIt < 5) {
        this->com->fprintf(stderr, "*** Warning: Not yet implemented\n");
        //fprintf(stderr,"*** Warning: Newton solver redoing iteration %d\n", it+1);
        //Q = rhs;
        //--it;
        //++fsIt;
      }
      else{
        this->com->fprintf(stderr, "*** Exiting\n");
        exit(1);
      }
    }
		breakloopNow = breakloop2(breakloop);
		if (breakloopNow) break;
   }	// end Newton loop

	savedUnormAccum();
	if (fsIt > 0 && checkFailSafe(U) == 1)
		resetFixesTag();

  if (it == maxItsNewton && maxItsNewton != 1 && maxItsNewton !=0) {
    this->com->fprintf(stderr, "*** Warning: ROM Newton solver reached %d its", maxItsNewton);
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  this->timer->addFluidSolutionTime(t0);

	// output POD coordinates

  return (maxItsNewton == 0) ? 1 : it;

}
                                                                                                           
//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::solveLinearSystem(int it , Vec<double> &rhs, Vec<double> &dUrom)
{

  double *x = rhs.data();
  FullM myjac(jac);
  myjac.factor();
  myjac.reSolve(x);
  dUrom = rhs;  

  return 0;
}

//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitRomTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q)
{

  //DistSVec<double, dim> F(this->domain->getNodeDistInfo());

  this->spaceOp->computeResidual(*this->X, *this->A, Q, F, this->timeState);

  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);

  this->spaceOp->applyBCsToResidual(Q, F);

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &F, double stepLength)  {
	// merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  computeFullResidual(it,newQ,F);

  double merit = 0.0;
  merit += F.norm();	// merit function = 1/2 * (norm of full-order residual)^2
  merit *= merit;
  merit *= 0.5;

  return merit;

}
//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunctionDeriv(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &p, DistSVec<double, dim> &F)  { 

  double eps = mvpfd->computeEpsilon(Q,p);
  DistSVec<double, dim> newF(this->domain->getNodeDistInfo());

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + eps*p;
  computeFullResidual(it,newQ,newF);

  newF -= F;  // overwrite new Flux with finite difference
  newF *= (1.0/eps);

  double meritDeriv = 0.0;
  meritDeriv = newF*F; // Take inner product
    // merit function = norm of full-order residual
  return meritDeriv;

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::lineSearch(DistSVec<double, dim> &Q, Vec<double> &prom, int it, VecSet<DistSVec<double, dim> > &leftProj, double eps, bool &convergeFlag){
  // This function (along with "zoom"), finds a steplength satisfying the Strong Wolfe conditions in two steps:
  // 1) Find a valid bracket containing a point that satisfies these conditions
  // 2) Perform interpolation/bisection within the bracket to locate this valid point
  // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

  // Parameters for Wolfe conditions

  double c1 = 0.00001;  // for Armijo (sufficient decrease) condition
  double c2 = 0.0001; // relative gradient condition.
  //NOTE: We require 0<c1<c2<1 
  //NOTE: the smaller the value of c2, the closer the steplength will be to a local minimizer

  // Parameters for bracketing step (step 1)

  int maxIter = 100;  // max iterations for finding a bracket
  DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo()); // new state vector initialization
  double beta = 1.2;  // amount to increment alpha each time
  double alpha_max = 50.0;    // bound on step size alpha
  int maxedFlag = 0;  // flag to determine whether or not alpha_max has been reached
  int count = 0;  // counter for number of iterations
  double alpha = 1.0; // initial alpha is just the Newton step
  double merit, meritZero, meritOld;
  double meritDeriv, meritDerivZero,targetMeritDerivZero, meritDerivOld;

  DistSVec<double,dim> dQ(this->domain->getNodeDistInfo());
  DistSVec<double,dim> F(this->domain->getNodeDistInfo());
  Vec<double> From(nPod);
  expandVector(prom, dQ);

  //----- Obtain an initial bracket -----//

  meritZero = meritFunction(it, Q, dQ, F, 0.0);
  meritOld = meritZero;
  meritDerivZero = meritFunctionDeriv(it, Q, dQ, F);

  // Reverse the search direction if it is not a descent direction

  if (meritDerivZero > 0) {
   this->com->fprintf(stderr,"Reversing the search direction: original direction NOT a descent direction! \n");
    dQ = -1.0*dQ;
    prom = -1.0*prom;
    meritDerivZero = -meritDerivZero;
  }
  this->com->fprintf(stderr,"meritDerivZero = %e, meritZero = %e \n",meritDerivZero, meritZero);
 
  targetMeritDerivZero = eps*fabs(meritZero);  //Compare targetMeritDeriv to meritZero based on Taylor series explanation
 
  if (fabs(meritDerivZero) <= targetMeritDerivZero){ // convergence criterion based on stationarity of full-order residual norm
     convergeFlag = 1;
     return 0;
  }

  meritDerivOld = meritDerivZero;
  double alphaOld = 0.0;  //Left endpoint is zero
  
  while (count<maxIter){
    merit = meritFunction(it, Q, dQ, F, alpha); // evaluate merit function at current alpha
    if (merit>meritZero+c1*alpha*meritDerivZero || (merit>=meritOld && count>0)){
        // alphaLo=alphaOld, alphaHi=alpha
        this->com->fprintf(stderr,"Entering zoom from location 1, alpha = %e \n",alpha);
        alpha = zoom(alphaOld, alpha, meritOld, merit, meritDerivOld, meritDerivZero, meritZero, c1, c2, Q, dQ, F, it);
        return alpha;
    }
    Qnew = Q+alpha*dQ;  // Q at current alpha
    meritDeriv = meritFunctionDeriv(it, Qnew, dQ, F);

    if (fabs(meritDeriv)<=-c2*meritDerivZero){
       this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1*alpha*meritDerivZero);
       this->com->fprintf(stderr,"Returning alpha without zoom, alpha = %e \n",alpha);
       return alpha;  // a valid solution has been found
    }
    if (meritDeriv>=0){
       // alphaLo=alpha, alphaHi=alphaOld
       this->com->fprintf(stderr,"Entering zoom from location 2.  The alphas are: alpha = %e, alphaOld = %e \n",alpha,alphaOld);
       alpha = zoom(alpha, alphaOld, merit, meritOld, meritDeriv, meritDerivZero, meritZero, c1, c2, Q, dQ, F, it);
       return alpha;
    }

    //----- update values of alpha, merit function, and derivative -----//

    alphaOld = alpha;
    meritOld = merit;
    meritDerivOld = meritDeriv;

    //----- choose alpha between current alpha and alpha_max -----//

    alpha = beta*alphaOld;   //increment alpha by a constant factor

    //----- correct if alpha exceeds alpha_max -----

    while (alpha>alpha_max){
      beta=(1.0+beta)/2.0;    //decrease beta
      alpha=beta*alphaOld; //try a new value of alpha
      this->com->fprintf(stderr,"Decreasing beta because alpha_max violated \n");
    }
    ++count;
    
  }

  this->com->fprintf(stderr,"Leaving linesearch because max iterations was exceeded. Conditions NOT satisfied! alpha = %e\n",alpha);  
  return alpha;



/*  while (merit <= meritZero+c1*alpha*meritDerivZero && !(merit>=meritOld && count>0)){  // exit when Armijo condition fails
    Qnew = Q+alpha*dQ;
    meritDeriv = meritFunctionDeriv(it, Qnew, dQ, F);    //value of directional derivative at current alpha
    if (meritDeriv>=c2*meritDerivZero)
        break;    //Either: 1) solution obtained, or 2) slope is positive (have a bracket)

    if (merit>meritOld-(alpha-alpha_L)*c2*fabs(meritDerivZero))
 	break;	// a valid bracket has been obtained    

    alpha_L = alpha;  //conditions satisfy left bracket criteria - can update!
    meritOld = merit;   //Save old function
    meritDerivOld = meritDeriv;   //Save old derivative
    alpha = beta*alpha;   //increment alpha by a constant factor

    
    //----- correct if alpha exceeds alpha_max -----
    
    while (alpha>alpha_max){
      beta=(1.0+beta)/2.0;    //decrease beta
      alpha=beta*alpha_L; //try a new value of alpha
      this->com->fprintf(stderr,"Decreasing beta because alpha_max violated \n");
    }
    merit = meritFunction(it, Q, dQ, F, alpha);
    ++count;
    if (count>maxIter){
        this->com->fprintf(stderr,"Maximum # iterations exceeded.  Exiting...\n");
        maxedFlag = 1;  // the constraints are active, so a bracket hasn't been found
        break;
    }
  }
 


  // ----- A valid initial bracket has been obtained----- 
 
  // ----- ZOOM - Perform quadratic approximations to reduce interval ----- 

 
  // using left and right bracket, initialize zoom phase

 
  double meritHi=merit;
  double meritLo=meritOld;
  double meritDerivLo=meritDerivOld; 
  double alpha_Hi=alpha;
  double alpha_Lo=alpha_L;

 
  if (maxedFlag==0) {// only zoom in if the original bracket was valid (not true if alpha_R=alpha_max!
    while ((meritLo>meritZero+c1*alpha_Lo*meritDerivZero || fabs(meritDerivLo)>-c2*meritDerivZero) && fabs(alpha_Hi-alpha_Lo)>tol) {  // continue while strong Wolfe conditions not satisfied or bracket small enough

 
      // Perform a quadratic interpolation using slope of the low point

      FullM Vdm(3); //Vandermonde matrix
      Vdm[0][0] = alpha_Hi*alpha_Hi; Vdm[0][1] = alpha_Hi; Vdm[0][2] = 1.0;
      Vdm[1][0] = 2.0*alpha_Lo; Vdm[1][1] = 1.0; Vdm[1][2] = 0.0;
      Vdm[2][0] = alpha_Lo*alpha_Lo; Vdm[2][1] = alpha_Lo; Vdm[2][2] = 1.0; 
      double *VdmRHS = new double[3];
      VdmRHS[0] = meritHi; VdmRHS[1] = meritDerivLo; VdmRHS[2] = meritLo;
      Vdm.factor();
      Vdm.reSolve(VdmRHS); //obtain interpolant coefficients
      alpha = -VdmRHS[1]/(2.0*VdmRHS[0]);//new value from minimizing quadratic model
      delete[] VdmRHS;
 
      // ----- Test whether or not the interval has been sufficiently reduced -----

 
      this->com->fprintf(stderr,"alpha = %f, alpha_Lo = %f, alpha_Hi = %f\n",alpha,alpha_Lo,alpha_Hi);
      if ((alpha-alpha_Lo)/(alpha_Hi-alpha_Lo)<epsilon || (alpha_Hi-alpha)/(alpha_Hi-alpha_Lo)<epsilon) 
        alpha=(alpha_Lo+alpha_Hi)/2.0; //use bisection if insufficient reduction
         
 
      // ----- Evaluate next high, low points -----

      merit = meritFunction(it, Q, dQ, F, alpha);
 
      // New point is the high point

 
      if (merit>meritZero+c1*alpha*meritDerivZero || merit >=meritLo){
        alpha_Hi=alpha;    //new point is the high point
        meritHi=merit;
      }
      // New point is the low point
      else {  //meritMid<=meritLo and meritMid satisfies Armijo
        Qnew = Q+alpha*dQ;
        meritDerivLo = meritFunctionDeriv(it, Qnew, dQ, F);
        if (meritDerivLo*(alpha_Hi-alpha_Lo)>=0.0){ // depending on slope meritDeriv, change alpha_Hi
          alpha_Hi = alpha_Lo;
          meritHi = meritLo;
        }
        alpha_Lo = alpha;    //update the new low point
        meritLo = merit;
      }
    }
  }
  //------ Loop exits when alpha satisfies both gamma conditions -----
*/
//  return alpha;

}
//-----------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::zoom(double alphaLo, double alphaHi, double meritLo, double meritHi, double meritDerivLo, double meritDerivZero, double meritZero, double c1, double c2, DistSVec<double, dim> Q, DistSVec<double, dim> dQ, DistSVec<double, dim> F, int it){

    // Given a valid bracket, this function finds an alpha within the bracket satisfying the Strong Wolfe conditions
    // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

    // initialize parameters

    int count = 0;
    int maxIter = 40; // max number of zooming steps
    double merit, meritDeriv, alpha;
    DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo());
    double epsilon = 0.1; // Percent reduction needed for quadratic approx before bisection is done
    double tolInterp = 0.05; // Determines if interpolation will be executed-we have had problems with the matrix being ill-conditioned
 
    while (count<maxIter){

      // Perform a quadratic interpolation using slope of the low point IF matrix will be well-conditioned using a heuristic

      if (fabs(alphaHi-alphaLo)>tolInterp){
      
            FullM Vdm(3); // Vandermonde matrix
   
            Vdm[0][0] = alphaLo*alphaLo; Vdm[0][1] = alphaLo; Vdm[0][2] = 1.0;
            Vdm[1][0] = 2.0*alphaLo; Vdm[1][1] = 1.0; Vdm[1][2] = 0.0;
            Vdm[2][0] = alphaHi*alphaHi; Vdm[2][1] = alphaHi; Vdm[2][2] = 1.0;
            double *VdmRHS = new double[3];
            VdmRHS[2] = meritHi; VdmRHS[1] = meritDerivLo; VdmRHS[0] = meritLo;
            Vdm.Factor();
            Vdm.ReSolve(VdmRHS); //obtain interpolant coefficients
            alpha = -VdmRHS[1]/(2.0*VdmRHS[0]);//new value from minimizing quadratic model
   
            delete[] VdmRHS;
      }

      // ----- Use bisection if interpolation was bypassed or if there was an insufficient reduction-----//

      if (fabs(alphaHi-alphaLo)<=tolInterp || (alpha-alphaLo)/(alphaHi-alphaLo)<epsilon || (alphaHi-alpha)/(alphaHi-alphaLo)<epsilon)
        alpha=(alphaLo+alphaHi)/2.0; // bisection

      // Evaluate merit function

      merit = meritFunction(it, Q, dQ, F, alpha);

      if (merit>meritZero+c1*alpha*meritDerivZero || merit>=meritLo){
           alphaHi = alpha;
           meritHi = merit;
      }
      else{
           Qnew = Q+alpha*dQ;
           meritDeriv = meritFunctionDeriv(it, Qnew, dQ, F);
           if (fabs(meritDeriv)<=-c2*meritDerivZero){
                this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1*alpha*meritDerivZero);
                this->com->fprintf(stderr,"alpha = %e, # zoom iter = %d\n",alpha, count);

                return alpha;  // alpha satisfies Strong Wolfe conditions: exit loop
           }
           if (meritDeriv*(alphaHi-alphaLo)>=0){
                alphaHi = alphaLo;
                meritHi = meritLo;
           }
           alphaLo = alpha;
           meritLo = merit;
           meritDerivLo = meritDeriv; 
      }
      ++count;
    }
    this->com->fprintf(stderr,"Leaving zoom because max iterations was exceeded. Conditions NOT satisfied! Count = %d, alpha = %e\n",count, alpha);
    return alpha;
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::expandVector(Vec<double> &romV, DistSVec<double, dim> &fullV)  {

  fullV = 0.0;
  for (int iVec = 0; iVec < nPod; iVec++)
    fullV += romV[iVec]*pod[iVec];
}


//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::projectVector(VecSet<DistSVec<double, dim> > &leftProj, DistSVec<double,dim> &fullV, Vec<double> &romV)  {

	transMatVecProd(leftProj,fullV,projVectorTmp);

  for (int iVec = 0; iVec < pod.numVectors(); iVec++)
    romV[iVec] = projVectorTmp[iVec];
}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{

  this->com->fprintf(stderr, "*** Warning: Checkfailsafe Not yet implemented\n");
  int failSafeNewton = 0;
/*
  if (!failSafeNewton) return 0;

  if (!tag)
    tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *tag);
  this->spaceOp->fix(*tag);
*/
  return failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

  mvpfd->evaluate(it, *this->X, *this->A, Q, F);
  
  for (int iPod = 0; iPod < nPod; iPod++)
    mvpfd->apply(pod[iPod], AJ[iPod]);

}
//------------------------------------------------------------------------------

/*
template<int dim>
int ImplicitRomTsDesc<dim>::solveNonLinearSystem(DistSVec<double, dim> &U, int _it)  {

  double t0 = this->timer->getTime();

  //  write our own newton solver
  int it = 0;
  int fsIt = 0;
  int Git = _it-1;

  DistSVec<double, dim> F(this->domain->getNodeDistInfo());

  Vec<double> From(nPod);
  Vec<double> dFrom(nPod);
  Vec<double> Urom(nPod);
  Vec<double> rhs(nPod);
  projectVector(pod, U, Urom);
  double res0;
  double res;
  double alpha;
  bool convergeFlag=0;

  // AJ is J0*Phi
  // VecSet<DistSVec<double, dim> > AJ(nPod,this->domain->getNodeDistInfo());


  for (it = 0; it < maxItsNewton; it++)  {

    computeFullResidual(it, U, F);

    // Compute reduced Jacobian or do Broyden update for the reduced Jacobian
    
    if (it==0 && (Git % JacSkipNewton)==0) {
      this->com->fprintf(stderr," ... Computing exact reduced Jacobian \n");
      computeAJ(it, U, F, AJ);  // Computes reduced Jacobian J0 and AJ=J0*Phi
      projectVector(AJ, F, From);
    }
    else {
      projectVector(AJ, F, From);
      if (it > 0) {
        dFrom = From-Fromold;
        broydenUpdate(dFrom, dUrom);
      }
    }

    // Using stationary criterion for the full-order residual (meritDerivZero is low) as a convergence criterion.
    // This is consistent with the LINE SEARCH objective function.

    Fromold = From;

    rhs = -1.0 * From;

    // compute search direction

    solveLinearSystem(it, rhs, dUrom);

    // do line search (linesearch exits with alpha=0 and convergenceFlag if convergence criteria is satisfied)

    alpha = lineSearch(U,dUrom,it,AJ,epsNewton, convergeFlag);
    if (it > 0 && convergeFlag == 1) break;
    this->com->fprintf(stderr," \tit = %d, alpha = %e\n",it, alpha);

    dUrom *= alpha;
    Urom += dUrom;
    expandVector(Urom, U);
    // verify that the solution is physical
    if (checkSolution(U)) {
      if (checkFailSafe(U) && fsIt < 5) {
        this->com->fprintf(stderr, "*** Warning: Not yet implemented\n");
        //fprintf(stderr,"*** Warning: Newton solver redoing iteration %d\n", it+1);
        //Q = rhs;
        //--it;
        //++fsIt;
      }
      else{
        this->com->fprintf(stderr, "***Exiting\n");
        exit(1);
      }
    }
   }
    if (fsIt > 0 && checkFailSafe(U) == 1) 
      resetFixesTag();

  if (it == maxItsNewton && maxItsNewton != 1) {
    this->com->fprintf(stderr, "*** Warning: ROM Newton solver reached %d its", maxItsNewton);
  }

  this->timer->addFluidSolutionTime(t0);

  return it;

}
*/

template<int dim>
void ImplicitRomTsDesc<dim>::writeStateRomToDisk(int it, double cpu)  {

	this->output->writeStateRomToDisk(it, cpu, nPod, UromTotal);

}

template<int dim>
void ImplicitRomTsDesc<dim>::saveNewtonSystemVectorsAction(const int totalTimeSteps) {
	// only do for PG and Galerkin

  DistSVec<double, dim> AJsol(this->domain->getNodeDistInfo()); //CBM--NEED TO CHANGE NAME OF DISTVECTOR
	AJsol = 0.0;
	for (int i=0; i<this->nPod; ++i)
		 AJsol += this->AJ[i] * this->dUrom[i]; 

	// saving this->AJ * this->dUrom (for GappyPOD)
	// for now, do not output on last iteration (first argument = false)
	this->writeBinaryVectorsToDiskRom(false, totalTimeSteps, 0.0, &(this->F), &AJsol, &(this->AJ));
	
}
template<int dim>
void ImplicitRomTsDesc<dim>::savedUnormAccum() {
	
	for (int iPod = 0; iPod < nPod; ++iPod) {
		dUnormAccum[0][iPod] += fabs(dUrom[iPod]);	// 1 norm
		dUnormAccum[1][iPod] += dUrom[iPod] * dUrom[iPod];	// 2 norm
	}
}

template<int dim>
void ImplicitRomTsDesc<dim>::readPodBasis(const char *fileName) {

	string fileNameState;
	determineFileName(fileName, ".sampledROBState", ioData->input.gnatPrefix, fileNameState);

	this->domain->readPodBasis(fileNameState.c_str(), nPod,
		pod,this->ioData->rom.basisType == ModelReductionData::SNAPSHOTS);
}

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop1(const bool breakloop) {

	return breakloop;
	
}

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop2(const bool breakloop) {

	return false;

}

template<int dim>
void ImplicitRomTsDesc<dim>::determineFileName(const char *fileNameInput, const char
		*fileNameExtension, const char *prefix, string &fileName) {

	if (strcmp(fileNameInput,"") == 0) {
		fileName = prefix;
	 	fileName +=	fileNameExtension;
	}
	else {
		fileName = fileNameInput;
	}
}

