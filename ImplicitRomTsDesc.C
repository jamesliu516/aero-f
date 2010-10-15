#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <NewtonSolver.h>

#include <TsOutput.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::ImplicitRomTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(ioData, geoSource, dom), pod(0, dom->getNodeDistInfo()), AJ(0, dom->getNodeDistInfo()) {

  tag = 0;

  maxItsNewton = ioData.ts.implicit.newton.maxIts;
  epsNewton = ioData.ts.implicit.newton.eps;  
  JacSkipNewton = ioData.ts.implicit.newton.JacSkip;

  //ns = new NewtonSolver<ImplicitRomTsDesc<dim> >(this);
  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;
  mvpfd = new MatVecProdFD<dim,dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);

  // read Pod Basis
  nPod = ioData.Rob.numROB;
  dom->readPodBasis(ioData.input.podFile, nPod, pod);

  globalSubSet = 0;
  locNodeSet = 0;

  RomSolver = ioData.Rob.romsolver;

  if (RomSolver == 2) {
	// KEVIN FIX!
		dom->readInterpNode(ioData.input.sampleNodes, nIntNodes, globalSubSet, locNodeSet);
		if (ioData.input.aMatrix) dom->readInterpMatrix(ioData.input.aMatrix, dimInterpMat, interpMat1);
		if (ioData.input.bMatrix) dom->readInterpMatrix(ioData.input.bMatrix, dimInterpMat, interpMat2);

		computeRestrictInfo();
  }

  MemoryPool mp;
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

  AJ.resize(nPod);
  dUrom.resize(nPod);
  Fromold.resize(nPod);
  
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::~ImplicitRomTsDesc()
{
  if (tag) delete tag;
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

    computeFunction(it, U, F);

    // Compute reduced Jacobian or do Broyden update for the reduced Jacobian
    
    if (it==0 && (Git % JacSkipNewton)==0) {
      this->com->fprintf(stderr," ... Computing exact reduced Jacobian \n");
      computeJacobian(it, U, F, AJ);  // Computes reduced Jacobian J0 and AJ=J0*Phi
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

//------------------------------------------------------------------------------

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
  Vec<double> Urom(nPod); // reduced coordinates
  Vec<double> rhs(nPod);
  //projectVector(pod, U, Urom);	// KTC change
	Urom = 0.0;	// total solution increment in ROM coordinates (the unknowns for reduced problem)
  DistSVec<double, dim> dUfull(this->domain->getNodeDistInfo());	// solution increment at EACH NEWTON ITERATION in full coordinates
  dUfull = 0.0;	// initial zero increment

  double res0;
  double res;
  double target;
  double alpha;
  bool convergeFlag=0;

  DistSVec<double, dim> Test(this->domain->getNodeDistInfo()); //CBM--NEED TO CHANGE NAME OF DISTVECTOR

  for (it = 0; it < maxItsNewton; it++)  {

    //this->com->fprintf(stderr," --- Newton It # %d\n",it);//CBM

    computeFunction(it, U, F);

    switch (RomSolver) {
      case 0: // Petrov-Galerkin
        computeJacobian(it, U, F, AJ); // KTC: instead, compute QR factorization?
        projectVector(AJ, F, From);
        rhs = -1.0 * From;

        // saving residual vectors (for GappyPOD)
        //writeBinaryVectorsToDisk1(false, _it, 0.0, F, Dummy);

        break;

      case 1: // Broyden
        if (it==0 && (Git % JacSkipNewton)==0) {
          this->com->fprintf(stderr," ... Computing exact reduced Jacobian \n");
          computeJacobian(it, U, F, AJ);
          projectVector(AJ, F, From);
        }
        else {
          projectVector(AJ, F, From);
          if (it > 0) {
            dFrom = From-Fromold;
            broydenUpdate(dFrom, dUrom);
          }
        }
        Fromold = From;
        rhs = -1.0 * From;
        break;

      case 2:  // Gappy POD
        computeJacobianGappy(it, U, F, rhs, AJ);
        break;
    } 
    res = rhs*rhs;

    if (res < 0.0){
      fprintf(stderr, "*** negative residual: %e\n", res);
      exit(1);
    }
    res = sqrt(res);

    if (it == 0) {
      target = epsNewton*res;
      res0 = res;
    }

    if (res == 0.0 || res <= target) break;

    solveLinearSystem(it, rhs, dUrom);
   

// LINE SEARCH
/*
    // do line search (linesearch exits with alpha=0 and convergenceFlag if convergence criteria is satisfied)

    alpha = lineSearch(U,dUrom,it,AJ,epsNewton, convergeFlag);
    if (it > 0 && convergeFlag == 1) break;

    dUrom *= alpha;
*/
// END LINE SEARCH

    expandVector(dUrom, dUfull);	// solution increment in full coordinates
    Urom += dUrom;// solution increment in reduced coordinates
		U += dUfull;

/*
//    if (RomSolver == 0) {
      Test = 0.0;
      for (int i=0; i<nPod; ++i)
         Test += AJ[i] * dUrom[i]; 

      // saving AJ * dUrom (for GappyPOD)
      writeBinaryVectorsToDisk1(false, _it, 0.0, F, Test);
//    }
*/

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
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  this->timer->addFluidSolutionTime(t0);

  return it;

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
void ImplicitRomTsDesc<dim>::computeFunction(int it, DistSVec<double, dim> &Q, 
					  DistSVec<double, dim> &F)
{

  //DistSVec<double, dim> F(this->domain->getNodeDistInfo());

  this->spaceOp->computeResidual(*this->X, *this->A, Q, F, this->timeState);

  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);

  this->spaceOp->applyBCsToResidual(Q, F);


}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &F, double stepLength)  {

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  computeFunction(it,newQ,F);

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
  computeFunction(it,newQ,newF);

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
void ImplicitRomTsDesc<dim>::broydenUpdate(Vec<double> &dFrom, Vec<double> &dUrom) {

  // Broyden update of the form B{k+1}=Bk+(yk-Bk*sk)*sk^T/(sk^T*sk) where yk is change in function and sk is step size

  Vec<double> zrom(nPod);
  zrom = dFrom;	// zrom = (yk-Bk*sk)
  for (int iPod = 0; iPod < nPod; ++iPod) { // KTC: parallelize
    for (int jPod = 0; jPod < nPod; ++jPod)
     zrom[iPod] -= jac[iPod][jPod]*dUrom[jPod];
  }

  double invNormSq = 1.0/ (dUrom*dUrom);
  zrom *= invNormSq;
  
  for (int iPod = 0; iPod < nPod; ++iPod) { // KTC: parallelize
    for (int jPod = 0; jPod < nPod; ++jPod)
      jac[iPod][jPod] += zrom[iPod]*dUrom[jPod];
  }


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

  for (int iVec = 0; iVec < pod.numVectors(); iVec++)
    romV[iVec] = leftProj[iVec]*fullV;
}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::recomputeFunction(DistSVec<double, dim> &Q, Vec<double> &rhsRom)  {

  //DistSVec<double, dim> rhs(this->domain->getNodeDistInfo());

  //this->spaceOp->recomputeRHS(*this->X, Q, rhs);

  //projectVector(pod, rhs, rhsRom);
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
void ImplicitRomTsDesc<dim>::computeJacobian(int it, DistSVec<double, dim> &Q,
                                                 DistSVec<double, dim> &F, VecSet<DistSVec<double, dim> > &AJ)  {

  //DistSVec<double, dim> F(this->domain->getNodeDistInfo());
  //DistSVec<double, dim> jacPart1(this->domain->getNodeDistInfo());

  //Vec<double> jacCol(nPod);

  //expandVector(From, F);

  mvpfd->evaluate(it, *this->X, *this->A, Q, F);
  
  jac.setNewSize(nPod,nPod);

/*  if (galerkin) {

    // populate jac
    for (int iPod = 0; iPod < nPod; iPod++)  {
  
      mvpfd->apply(pod[iPod], jacPart1);
    
      projectVector(pod, jacPart1, jacCol);
    
      for (int iRow = 0; iRow < nPod; iRow++) 
        jac[iRow][iPod] = jacCol[iRow];
    }
  }
  else {
*/


  //VecSet<DistSVec<double, dim> > AJ(nPod,this->domain->getNodeDistInfo());

  // populate jac
  for (int iPod = 0; iPod < nPod; iPod++)  
    mvpfd->apply(pod[iPod], AJ[iPod]);

	// KTC: LS alert!

  for (int iRow = 0; iRow < nPod; ++iRow) {
    for (int iCol = 0; iCol <= iRow; ++iCol) {
      jac[iRow][iCol] = AJ[iRow]*AJ[iCol];
      if (iRow > iCol)
        jac[iCol][iRow] = jac[iRow][iCol];
    }
  } 
  
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::computeJacobianGappy(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F, 
                                                  Vec<double> &rhs, VecSet<DistSVec<double, dim> > &AJ)  {
  
  // The goal of this function is to compute rhs and jac for the newton step
  // nIntNodes = # nodes at which interpolation occurs
  // dimension of interpMat1 = nIntNodes * dim, interpMat2 = nIntNodes * dim

  int debugging = 0;

  mvpfd->evaluate(it, *this->X, *this->A, Q, F);  // prepares mvpfd->apply

  jac.setNewSize(nPod,nPod); // actual Jacobian to be used in Newton iterations

  for (int iPod = 0; iPod < nPod; ++iPod)
    mvpfd->apply(pod[iPod], AJ[iPod]);  // AJ is AJ: the VecSet of DistSVec with AJ[i] = dr/dw*Phi[i]

  FullM AJRestrict(nIntNodes*dim,nPod); // the matrix of dr/dw*Phi[i] restricted to interp nodes
  AJRestrict = 0.0;
  FullM resRestrict(nIntNodes*dim,1); // residual restricted to interp nodes
  resRestrict = 0.0;

  // info for parallel operations
  int numLocSub = this->domain->getNumLocSub();
  int nTotCpus = this->com->size();
  int thisCPU = this->com->cpuNum();
  DistInfo &nodeDistInfo = this->domain->getNodeDistInfo();
  SubDomain** subD = this->domain->getSubDomain();

  // compute AJ and F restricted to the specified nodes
  int oldLocalSub = -1; 
  int localNode,localSub,currentNodeIndex;
  for (int iPod = 0; iPod < nPod; ++iPod) { // loop over all POD vectors
    for (int iMyInterpNode = 0; iMyInterpNode < myNNodeInt; ++iMyInterpNode) { // loop over local interpolation nodes
      currentNodeIndex = myInterpNodes[iMyInterpNode];  // global index of current interpolation node
      localNode = locNodeSet[currentNodeIndex];  // the local node number
      localSub = myLocalSubSet[iMyInterpNode];  // the local subdomain number
      //if (iMyInterpNode == 0 || localSub != oldLocalSub) { // only reload subdomain part of locF, locAJ if needed
      //if (iPod == 0){
      double (*locF)[dim] = F.subData(localSub); // compute local F (should only be loaded for iPod == 0)
      //}
      double (*locAJ)[dim] = AJ[iPod].subData(localSub); // compute local AJ
      //}

      // KTC test: are the CPU #, global Sub #, local Node #, and myNNodeInt correct for a certain global currentNodeIndex?
      if (debugging) {
        fprintf(stderr,"currentNodeIndex %d, CPU %d, globalSub %d, localNode %d, localSub %d, myNNodeInt %d\n",currentNodeIndex,thisCPU,globalSubSet[currentNodeIndex],localNode,localSub, myNNodeInt);
      }
      for (int iDim = 0 ; iDim < dim; ++iDim) {
        // compute restricted residual and jacobian
        if (iPod == 0) resRestrict[currentNodeIndex*dim+iDim][0] = locF[localNode][iDim];  // compute for RHS (only one)
        AJRestrict[currentNodeIndex*dim+iDim][iPod] = locAJ[localNode][iDim]; // fill in the matrix of dr/dw*Phi[i] restricted to interp
      }
      oldLocalSub = localSub;
    }
  }
/*
  delete [] locAJ;
  delete [] locF;
  delete [] subD;
*/

  this->com->globalSum(nIntNodes*dim*nPod,AJRestrict.data()); // ensure all CPUs have the same copy
  this->com->globalSum(nIntNodes*dim,resRestrict.data()); // ensure all CPUs have the same copy

/*
  FullM Arhat(nIntNodes*dim,1);
  Arhat = interpMat*resRestrict;
  FullM rhatArhatmat = resRestrict^Arhat;
  double rhatArhat = rhatArhatmat[0][0];
  this->com->fprintf(stderr,"R^T A R = %e\n",rhatArhat);
*/
  // parallel implementation of AJRestrictTMat1 = AJRestrict^T * interpMat1
  // parallel implementation of AJRestrictTMat2 = AJRestrict^T * interpMat2

  int maxIndex, loadBal, loadBalMod, myMinIndex, myMaxIndex;  // indices used for parallel operations

  FullM AJRestrictTMat1(AJRestrict.numCol(),interpMat1.numCol());
  AJRestrictTMat1 = 0.0;  // initialize to zero

  rowPartition(myMinIndex,myMaxIndex,AJRestrictTMat1.numRow());
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < AJRestrictTMat1.numCol(); ++j){
      for (int k = 0; k < interpMat1.numRow(); ++k) {
        AJRestrictTMat1[i][j] += AJRestrict[k][i]*interpMat1[k][j];
      }
    }
  }
  this->com->globalSum(AJRestrictTMat1.numRow()*AJRestrictTMat1.numCol(),AJRestrictTMat1.data()); // ensure all CPUs have the same copy

  FullM AJRestrictTMat2(AJRestrict.numCol(),interpMat2.numCol());
  AJRestrictTMat2 = 0.0;  // initialize to zero

  rowPartition(myMinIndex,myMaxIndex,AJRestrictTMat2.numRow());
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < AJRestrictTMat2.numCol(); ++j){
      for (int k = 0; k < interpMat2.numRow(); ++k) {
        AJRestrictTMat2[i][j] += AJRestrict[k][i]*interpMat2[k][j];
      }
    }
  }
  this->com->globalSum(AJRestrictTMat2.numRow()*AJRestrictTMat2.numCol(),AJRestrictTMat2.data()); // ensure all CPUs have the same copy

  // compute jac = AJRestrictTMat1*AJRestrict

  jac = 0.0;  // initialize to zero
  rowPartition(myMinIndex,myMaxIndex,jac.numRow(),1);  // symmetric
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j <= i; ++j){
      for (int k = 0; k < AJRestrict.numRow(); ++k) {
        jac[i][j] += AJRestrictTMat1[i][k]*AJRestrict[k][j];
      }
      if (i > j)
        jac[j][i] = jac[i][j];
    }
  }

  this->com->globalSum(jac.numRow()*jac.numCol(),jac.data()); // ensure all CPUs have the same copy
  // compute rhs = -1.0*AJRestrictTMat2*resRestrict

  rowPartition(myMinIndex,myMaxIndex,rhs.size());
  rhs = 0.0; 
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < resRestrict.numRow(); ++j)
      rhs[i] -= AJRestrictTMat2[i][j]*resRestrict[j][0]; // NOTE: the RHS is NEGATIVE!
  }
  this->com->globalSum(rhs.size(),rhs.data()); // ensure all CPUs have the same copy
 
  // KTC test: given AJRestrict, resRestrict, interpMat1, and interpMat2, are the computation of jacobian and rhs correct?
  FullM AJRestrictTMat1test, AJRestrictTMat2test, jacTest,rhsTest;
  if (debugging){
    AJRestrictTMat1test= AJRestrict^interpMat1;
    AJRestrictTMat2test= AJRestrict^interpMat2;
    jacTest = AJRestrictTMat1test*AJRestrict;
    rhsTest = AJRestrictTMat2test*resRestrict;
    rhsTest *= -1.0; 
    this->com->fprintf(stderr,"jac is \n");
    if (thisCPU == 0 ) {jac.print();}
    this->com->fprintf(stderr,"jacTest is \n");
    if (thisCPU == 0 ) {jacTest.print();}
    if (thisCPU == 0 ) {fprintf(stderr, "rhs is \n");}
    for (int j = 0; j < rhs.size(); ++j) { 
      if (thisCPU == 0 ) {fprintf(stderr, "%e \n", rhs[j]);}
    }
    this->com->fprintf(stderr,"rhsTest is \n");
    if (thisCPU == 0 ) {rhsTest.print();}
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::computeRestrictInfo() {

  // computes 1) myNNodeInt: the number of interpolation nodes contained on this CPU
  //          2) myLocalSubSet: the local subdomain numbers
  //          3) myInterpNodes: the global interpolation nodes it contains
  // NOTE: nodes are lumped according to subdomain for efficiency

  // local CPU information

  int numLocSub = this->domain->getNumLocSub();
  DistInfo &nodeDistInfo = this->domain->getNodeDistInfo();
  SubDomain** subD = this->domain->getSubDomain();
  myNNodeInt = 0;
  int globalSubNum;
  int localNode;

  int *myInterpNodesTemp = new int[nIntNodes];
  int *myLocalSubSetTemp = new int[nIntNodes];

  for (int iSub = 0; iSub < numLocSub; ++iSub) { // loop on subdomains on this CPU
    globalSubNum = subD[iSub]->getGlobSubNum(); 
    bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);  // array of locMasterFlag
    for (int iNodeInt = 0; iNodeInt < nIntNodes; ++iNodeInt) {  // loop over global interpolation nodes (zero to nIntNodes-1)
      if (globalSubNum == globalSubSet[iNodeInt]){ // if the current interpolation node is on this subdomain
        localNode = locNodeSet[iNodeInt]; // the local node of interest
        if (locMasterFlag[localNode]){ // if this node is the master on this domain
          myInterpNodesTemp[myNNodeInt] = iNodeInt;
          myLocalSubSetTemp[myNNodeInt]= iSub;
          ++myNNodeInt;
        }
      }
    }
    //delete [] locMasterFlag;
  }

  myInterpNodes = new int[myNNodeInt];
  myLocalSubSet = new int[myNNodeInt];
  for (int i = 0; i < myNNodeInt; ++i) {
    myInterpNodes[i] = myInterpNodesTemp[i];
    myLocalSubSet[i] = myLocalSubSetTemp[i];
  }
/*
  delete [] myInterpNodesTemp;
  delete [] myLocalSubSetTemp;
  delete [] subD;
*/
}
//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::rowPartition(int &myMinIndex, int &myMaxIndex, int nRow, int sym ) {
 
 // this function partitions the rows of a matrix across processors (load balancing)

 long int loadBal;
 int loadBalMod;
 int nTotCpus = this->com->size();
 int thisCPU = this->com->cpuNum();
 switch (sym) {  // is the operation symmetric?
   case 0:  // evenly split rows across processors
     loadBal = nRow/nTotCpus;
     loadBalMod = nRow % nTotCpus;
     myMinIndex = loadBal*thisCPU + (thisCPU < loadBalMod)*thisCPU + (thisCPU >= loadBalMod)*loadBalMod;
     myMaxIndex = loadBal*(thisCPU+1) + ((thisCPU+1) < loadBalMod)*(thisCPU+1) + ((thisCPU +1) >= loadBalMod)*loadBalMod;
   break;
   case 1:  // split rows in a staggered manner
     loadBal = (nRow*nRow)/nTotCpus;
     int *minVal = new int[nTotCpus];
     int *maxVal = new int[nTotCpus];
     minVal[0] = 0;
     maxVal[0] = (int) sqrt(double(loadBal));
     for (int i = 1; i < nTotCpus; ++i){
       minVal[i] = maxVal[i-1];
       maxVal[i] = (int) sqrt(double(loadBal+ minVal[i]*minVal[i]));
     }    
     maxVal[nTotCpus - 1] = nRow;
     myMinIndex = minVal[thisCPU];
     myMaxIndex = maxVal[thisCPU];
     delete [] minVal;
     delete [] maxVal;
   break;
 }
 if (myMaxIndex > nRow || myMinIndex > nRow || myMinIndex > myMaxIndex) {
   fprintf(stderr, "*** Problem with rowPartition!!!");
 }
 // KTC: output minVal, maxVal?
}

