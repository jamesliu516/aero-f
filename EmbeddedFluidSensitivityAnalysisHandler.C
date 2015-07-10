#include <EmbeddedFluidSensitivityAnalysisHandler.h>

#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistVector.h>
#include <MatVecProd.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MemoryPool.h>
#include <cmath>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------

template<int dim>
EmbeddedFluidSensitivityAnalysisHandler<dim>::EmbeddedFluidSensitivityAnalysisHandler
(
  IoData &ioData,
  GeoSource &geoSource,
  Domain *dom
) :
ImplicitEmbeddedCoupledTsDesc<dim>(ioData, geoSource, dom),
domain(dom),
dXdS(dom->getNodeDistInfo()),
dXdSb(dom->getNodeDistInfo()),
Xc(dom->getNodeDistInfo()),
dFdS(dom->getNodeDistInfo()),
dUdS(dom->getNodeDistInfo()),
p(dom->getNodeDistInfo()),
dPdS(dom->getNodeDistInfo()),
Flux(dom->getNodeDistInfo()),
FluxFD(dom->getNodeDistInfo()),
Pin(dom->getFaceDistInfo()),
Uc(dom->getNodeDistInfo())
// Tests
, Xplus(dom->getNodeDistInfo())
, Xminus(dom->getNodeDistInfo())
, dX(dom->getNodeDistInfo())
{

  std::cout<< "EmbeddedFluidSensitivityAnalysisHandler :: \n";

  // Initialize
  step = 0;

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
    Xp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Xm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Ap = new DistVec<double>(domain->getNodeDistInfo());
    Am = new DistVec<double>(domain->getNodeDistInfo());
    Lp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Um = new DistSVec<double,dim>(domain->getNodeDistInfo());
  }  
  else if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    Lp = 0;
    Lm = 0;
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = 0;
    Um = 0;
   }
  else {
    Lp = 0;
    Lm = 0;
    Fp = 0;
    Fm = 0;
  }

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)  
  {
    if (ioData.ts.implicit.mvp == ImplicitData::H2)
    {
      mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, domain, this->geoState);
    }
    else
    {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp, domain, ioData);     
    }
  }
  else 
  {
    if (ioData.ts.implicit.mvp == ImplicitData::H2)
    {
      mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, 0, this->spaceOp, domain, this->geoState);
    }
    else
    {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit,  0,   this->geoState, this->spaceOp, domain, ioData);     
    }
  }
  
  pc = ImplicitEmbeddedTsDesc<dim>::template
    createPreconditioner<dim>(ioData.sa.ksp.pc, domain);
  
  ksp = this->createKrylovSolver(this->getVecInfo(), ioData.sa.ksp, mvp, pc, this->com);

  /*
  MemoryPool mp;
  mvp->exportMemory(&mp);
   pc->exportMemory(&mp);
  */

  std::cout << "IN FSA _fsi "<< this->linRecAtInterface << " " << this->viscSecOrder << " " << this->riemannNormal <<" \n";

  typename MatVecProd<dim,dim>::_fsi fsi = {
    this->distLSS,
   &this->nodeTag,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
    this->Nsbar,
   &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };
  mvp->AttachStructure(fsi);

  length  = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = domain->getNumLocSub();

  dFdS = 0.0;
  dUdS = 0.0;
  p    = 0.0;
  dPdS = 0.0;

  FluxFD = 0.0;
  Flux   = 0.0;

  reynolds0 = ioData.ref.reynolds_mu;
  kenergy0  = ioData.bc.inlet.kenergy;

}

//------------------------------------------------------------------------------

template<int dim>
EmbeddedFluidSensitivityAnalysisHandler<dim>::~EmbeddedFluidSensitivityAnalysisHandler()
{

  if (mvp) delete mvp;

  if (pc) delete pc;

  if (ksp) delete ksp;

}

//------------------------------------------------------------------------------


template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaRestartBcFluxs(IoData &ioData)
{

  double gamma  = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double      R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

// Remark: For internal flows the SA using inlet or outlet 
// Mach, Alpha and Beta should be specified in the input file

  ioData.bc.inlet.mach  = xmach;
  ioData.bc.outlet.mach = xmach;

  ioData.bc.inlet.alpha  = alprad;
  ioData.bc.outlet.alpha = alprad;

  ioData.bc.inlet.beta  = teta;
  ioData.bc.outlet.beta = teta;

  this->com->fprintf(stderr, " FluidSensitivityAnalysis values: \n");
  this->com->fprintf(stderr, " x Mach, alpha, beta = %20.17e  %20.17e %20.17e \n", xmach, alprad, teta);

  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {

     this->com->fprintf(stderr, "NON_DIMENSIONAL type non supported\n");
     exit(1);

  } else if (ioData.problem.mode == ProblemData::DIMENSIONAL) {
    
    //
    // Step 1: Re-scale all the parameters
    //

    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;

    this->com->fprintf(stderr, "-Ref   rho, P = %20.17e %20.17e \n", ioData.ref.rv.density,   ioData.ref.rv.pressure);
    this->com->fprintf(stderr, "-INLET rho, P M = %20.17e %20.17e \n", ioData.bc.inlet.density, ioData.bc.inlet.pressure, ioData.bc.inlet.mach);

    ioData.bc.inlet.density *= ioData.ref.rv.density;
    ioData.bc.inlet.pressure *= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde *= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy *= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps *= ioData.ref.rv.epsilon;
    ioData.bc.outlet.density *= ioData.ref.rv.density;
    ioData.bc.outlet.pressure *= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde *= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy *= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps *= ioData.ref.rv.epsilon;

    ioData.restart.etime *= ioData.ref.rv.time;
    ioData.restart.dt_nm1 *= ioData.ref.rv.time;
    ioData.restart.dt_nm2 *= ioData.ref.rv.time;
    ioData.restart.energy *= ioData.ref.rv.energy;
    ioData.bc.wall.temperature *= ioData.ref.rv.temperature;
    ioData.ts.timestep *= ioData.ref.rv.time;
    ioData.ts.maxTime *= ioData.ref.rv.time;
    ioData.rmesh.vx *= ioData.ref.rv.velocity;
    ioData.rmesh.vy *= ioData.ref.rv.velocity;
    ioData.rmesh.vz *= ioData.ref.rv.velocity;
    ioData.rmesh.ax *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep *= ioData.ref.rv.time;

    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time      *= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ *= ioData.ref.rv.velocity;
    }
    ioData.aero.pressure *= ioData.ref.rv.pressure;
    ioData.forced.timestep *= ioData.ref.rv.time;
    ioData.forced.frequency /= ioData.ref.rv.time;

    ioData.eqs.gravity_x *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth *= ioData.ref.length;

    //
    // Step 2: Restart all values based on new inlet values
    // Step 2.1: Reset values in ioData.ref
    //

    ioData.ref.mach     = ioData.bc.inlet.mach;
    ioData.ref.density  = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;

    double velocity = ioData.ref.mach * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);

    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);

    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
      (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);
    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;
    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->com->fprintf(stderr, "\n\n Reynolds = %e \n\n",ioData.ref.reynolds_mu);

    double dvelocitydMach  = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;

    //
    // Step 2.2: Reset values in ioData.ref.rv
    //

    ioData.ref.rv.mode = RefVal::DIMENSIONAL;
    ioData.ref.rv.density = ioData.ref.density;
    ioData.ref.rv.velocity = velocity;
    ioData.ref.rv.pressure = ioData.ref.density * velocity*velocity;
    ioData.ref.rv.temperature = gamma*(gamma - 1.0) * ioData.ref.mach*ioData.ref.mach * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density);
    ioData.ref.rv.viscosity_mu = viscosity;
    ioData.ref.rv.nutilde = viscosity / ioData.ref.density;
    ioData.ref.rv.kenergy = velocity*velocity;
    ioData.ref.rv.epsilon = velocity*velocity*velocity / ioData.ref.length;
    ioData.ref.rv.time = ioData.ref.length / velocity;    // Problem in RigidMeshMotionHandler
    ioData.ref.rv.force = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.energy = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.power = ioData.ref.density * velocity*velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.tvelocity = velocity / ioData.aero.displacementScaling;
    ioData.ref.rv.tforce = ioData.ref.rv.force / ioData.aero.forceScaling;
    ioData.ref.rv.tpower = ioData.ref.rv.power / ioData.aero.powerScaling;

    ioData.ref.rv.dvelocitydMach = dvelocitydMach;
    ioData.ref.rv.dtimedMach = - ioData.ref.length / (velocity * velocity) * dvelocitydMach;

    ioData.eqs.fluidModel.pmin /= ioData.ref.rv.pressure;

    this->com->fprintf(stderr, "REF+rv rho, P Vel = %20.17e %20.17e\n", ioData.ref.rv.density, ioData.ref.rv.pressure);

    //
    // Step 2.3: Reset values of bc.inlet
    //

    ioData.bc.inlet.density /= ioData.ref.rv.density;
    ioData.bc.inlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps /= ioData.ref.rv.epsilon;

this->com->fprintf(stderr, "NEW INLET rho, P, M= %20.17e %20.17e %20.17e\n", ioData.bc.inlet.density, ioData.bc.inlet.pressure, ioData.bc.inlet.mach);

    //
    // Step 2.4: Reset values of bc.outlet
    //

    ioData.bc.outlet.density /= ioData.ref.rv.density;
    ioData.bc.outlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps /= ioData.ref.rv.epsilon;

    ioData.restart.etime /= ioData.ref.rv.time;
    ioData.restart.dt_nm1 /= ioData.ref.rv.time;
    ioData.restart.dt_nm2 /= ioData.ref.rv.time;
    ioData.restart.energy /= ioData.ref.rv.energy;

    ioData.bc.wall.temperature /= ioData.ref.rv.temperature;
    ioData.linearizedData.stepsize = ioData.ts.timestep;
    ioData.ts.timestep /= ioData.ref.rv.time;             // Problem in RigidRollMeshMotionHandler
    ioData.ts.maxTime /= ioData.ref.rv.time;              // Problem in RigidRollMeshMotionHandler
    ioData.rmesh.vx /= ioData.ref.rv.velocity;            // Problem in RigidMeshMotionHandler
    ioData.rmesh.vy /= ioData.ref.rv.velocity;
    ioData.rmesh.vz /= ioData.ref.rv.velocity;
    ioData.rmesh.ax /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep /= ioData.ref.rv.time;          // Problem in AccMeshMotionHandler
    
    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time     /= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ /= ioData.ref.rv.velocity;
    }
    ioData.aero.pressure /= ioData.ref.rv.pressure;
    ioData.forced.timestep /= ioData.ref.rv.time;         // Problem in ForcedMeshMotionHandler
    ioData.forced.frequency *= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler

    ioData.eqs.gravity_x /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth /= ioData.ref.length;

    double theta_k = 1.0;
    double theta_w = 10.0;
    if (kenergy0 == pow(10.0, -theta_k) * theta_w / reynolds0) {
      ioData.bc.inlet.kenergy = pow(10.0, -theta_k) * theta_w / ioData.ref.reynolds_mu;
      ioData.bc.inlet.eps = ioData.eqs.tc.tm.ke.c_mu * ioData.bc.inlet.kenergy * theta_w;
    }

    Pin = ioData.aero.pressure;

    this->postOp->rstVarPostFcn(ioData);

    this->postOp->rstVar(ioData);

    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->spaceOp->rstVarFet(ioData);

    this->spaceOp->rstFluxFcn(ioData);

    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

    this->rstVarImplicitEmbeddedCoupledTsDesc(ioData);
    
    this->bcData->rstVar(ioData);

    (this->timeState->getData()).rstVar(ioData);

    this->timeState->rstVar(ioData);

    this->output->rstVar(ioData);

    this->restart->rstVar(ioData);

    this->data->rstVar(ioData);

    this->refVal->rstVar(ioData);
    
    this->varFcn->rstVar(ioData);

  }

  // initialize boundary condition fluxes
  this->bcData->initialize(ioData, *this->X);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaSemiAnalytical
(
  IoData &ioData, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistSVec<double,dim> &dF
)
{

  //
  // Error mesage for pointers
  //

  if (Fp == 0) {
    fprintf(stderr, "*** Error: Variable Fp does not exist!\n");
    exit(1);
  }
  if (Fm == 0) {
    fprintf(stderr, "*** Error: Variable Fm does not exist!\n");
    exit(1);
  }


  double dtLeft;
  double eps=ioData.sa.eps;
  double xmachc;
  double alpradc;
  double tetac;

  xmachc  = xmach;
  alpradc = alprad;
  tetac   = teta;

  *Fp = 0.0;
  *Fm = 0.0;
   dF = 0.0;

   *Xp = 0.0;
   *Xm = 0.0;
   *Ap = 0.0;
   *Am = 0.0;

  //
  // Compute the first flux for the FD approach
  //
  Xc  = X;
  *Xp = X + eps*dXdS;
  X   = *Xp;

  xmach  = xmachc  + eps*DFSPAR[0];
  alprad = alpradc + eps*DFSPAR[1];
  teta   = tetac   + eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;
  

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(X, A, U, 
                                 this->Wtemp, this->Wtemp,
				 this->distLSS, this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag,
				 *Fp, 
				 this->riemann, this->riemannNormal,
				 this->Nsbar, 0, this->ghostPoints);

  this->spaceOp->applyBCsToResidual(U, *Fp, this->distLSS);


  //
  // Compute the second flux for the FD approach
  //

  X   = Xc;
  *Xm = X - eps*dXdS;
  X   = *Xm;

  xmach  = xmachc  - eps*DFSPAR[0];
  alprad = alpradc - eps*DFSPAR[1];
  teta   = tetac   - eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(X, A, U, 
                                 this->Wtemp, this->Wtemp,
				 this->distLSS, this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag,
				 *Fm, 
				 this->riemann, this->riemannNormal,
				 this->Nsbar, 0, this->ghostPoints);

  this->spaceOp->applyBCsToResidual(U, *Fm, this->distLSS);

  dF=1.0/(2.0*eps)*((*Fp)-(*Fm));

  //
  // Reset the steady state
  //

  X=Xc;

  xmach  = xmachc;
  alprad = alpradc;
  teta   = tetac;

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  //this->spaceOp->computeGradP(X, A, U);******

  if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    std::cout << "NOT NOW\n";
    this->geoState->updateConfigSA();
    //this->bcData->initializeSA(ioData, X, 0, DFSPAR[0], DFSPAR[1], DFSPAR[2]);
  }

}


//------------------------------------------------------------------------------
/*
template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaAnalytical
(
  IoData &ioData, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistSVec<double,dim> &dFdS
)
{
 
  //
  // Computing the normal, derivative of the normal and of the control volume
  //
  this->geoState->computeDerivatives
  (
    X, dXdS, this->bcData->getVelocityVector(),
    this->bcData->getDerivativeOfVelocityVector(), dAdS
  );

  //
  // Computing the derivatives of the boundary fluxes
  //
  this->bcData->initializeSA
  (
    ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]
  );

  //
  // Computing the partial derivative of the flux with respect to the variables
  //
  this->spaceOp->computeDerivativeOfResidual
  (
    X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS, this->timeState
  );

  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);

}
*/
//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaSetUpLinearSolver
(
  IoData &ioData, 
  DistSVec<double,3> &X, 
  DistVec<double> &A, 
  DistSVec<double,dim> &U, 
  DistSVec<double,dim> &dFdS
)
{

  this->com->fprintf(stderr, "IN setup lin solver \n");

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X); //?

  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A); //?

  this->bcData->update(X); //?

  this->com->fprintf(stderr, "IN setup lin solver: pre spaceOP\n");
  this->spaceOp->computeResidual(X, A, U, 
                                 this->Wtemp, this->Wtemp,
				 this->distLSS, this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag,
				 FluxFD, 
				 this->riemann, this->riemannNormal,
				 this->Nsbar, 0, this->ghostPoints);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY){
    this->timeState->add_dAW_dt(1, *this->geoState, A, U, FluxFD, this->distLSS);
  }

  this->spaceOp->applyBCsToResidual(U, FluxFD, this->distLSS);

  this->com->fprintf(stderr, "IN setup lin solver: pre mvp evaluate\n");
  mvp->evaluate(0, X, A, U, FluxFD);
  this->com->fprintf(stderr, "IN setup lin solver: post mvp evaluate\n");

  DistMat<double,dim> *_pc = dynamic_cast<DistMat<double,dim> *>(pc);

  if (_pc) {

    MatVecProdFD<dim,dim>        *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
    MatVecProdH2<dim,double,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,double,dim> *>(mvp);

    if (mvpfd || mvph2) 
    {
      this->spaceOp->computeJacobian(X, A, U, 
				     this->distLSS, this->nodeTag, 
				     this->riemann, this->riemannNormal, 
				     this->Nsbar,   this->ghostPoints, *_pc, this->timeState);

      this->com->fprintf(stderr, "IN setup lin solver: pc\n");
      if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY){
        this->timeState->addToJacobian(A, *_pc, U);
      }

      this->spaceOp->applyBCsToJacobian(U, *_pc, this->distLSS);

    }

  }

  pc->setup();

  // Computing flux for compatibility correction of the derivative of the flux   
  this->spaceOp->computeResidual(X, A, U, 
                                 this->Wtemp, this->Wtemp,
				 this->distLSS, this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag,
				 Flux, 
				 this->riemann, this->riemannNormal,
				 this->Nsbar, 0, this->ghostPoints);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaLinearSolver
(
  IoData &ioData, 
  DistSVec<double,dim> &dFdS, DistSVec<double,dim> &dUdS
)
{

  dFdS *= (-1.0);
  dUdS = 0.0;

  this->embeddeddQ = 0.0;
  this->embeddedB.ghost() = 0.0;
  this->embeddedB.real()  = dFdS;
 
  ksp->setup(0, 0, this->embeddedB);

  int numberIteration;
  bool istop = false;
  int iter = 0;

  while ((istop == false) && (iter < 100))
  {
    this->com->fprintf(stderr, "iter %d \n", iter);

    numberIteration = ksp->solve(this->embeddedB, this->embeddeddQ);

    dUdS = this->embeddeddQ.real();

    if ((!ioData.sa.excsol) || (numberIteration < ioData.sa.ksp.maxIts))
      istop = true; 
    iter += 1;
  }

  dFdS *= (-1.0);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
int EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR(1)  -  Mach number differential
  // DFSPAR(2)  -  angle of attack differential
  // DFSPAR(3)  -  yaw angle differential

  // Start basic timer
  double MyLocalTimer = -this->timer->getTime();

  this->output->openAsciiFiles();

  // Reseting the configuration control of the geometry datas
  this->geoState->resetConfigSA();

  if (this->com->cpuNum() == 0) {
    outFile = fopen(ioData.sa.sensoutput,"w");
    if (outFile) 
      fclose(outFile);
  }

  xmach  = ioData.sa.machref;
  alprad = ioData.sa.alpharef;
  teta   = ioData.sa.betaref;

  this->com->fprintf(stderr, "xxx Mach, alpha, beta = %20.17e  %20.17e %20.17e \n", xmach, alprad, teta);

  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->computeMeshMetrics();

  this->updateStateVectors(U);  

  fsaSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);

  /*
  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) {

    double tag = 0.0;

    step = 0;
    dXdS = 0.0;
    dXdSb = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 1;

    while (true) {

      // Reading derivative of the overall deformation
      bool readOK = domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdSb);
      if(!readOK) break;

      // Checking if dXdSb has entries different from zero at the interior of the mesh
      this->postOp->checkVec(dXdSb);

      if (dXdSb.norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n *** ERROR *** No Mesh Perturbation \n\n");
        exit(1);
      }

      this->com->fprintf(stderr, "\n ***** Shape variable %d\n", step);

      // Updating the mesh
      dXdS = *this->X;
      mms->solve(dXdSb, dXdS);
      dXdS -= *this->X;
      

      // Check that the mesh perturbation is propagated
      if (dXdS.norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n !!! WARNING !!! No Mesh Perturbation !!!\n\n");
      }

      fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);
  
      fsaComputeSensitivities(ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

      dXdSb = 0.0;

      step = step + 1;

    }

    fsaPrintTextOnScreen("\n ***** Derivatives with respect to the mesh position were computed! \n");

  }
  */

  if (ioData.sa.sensMach == SensitivityAnalysis::ON_SENSITIVITYMACH) {

    //dXdS = 0.0;
    //dAdS = 0.0;
    DFSPAR[0] = 1.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 2;

    fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsaComputeSensitivities(ioData, "Derivatives with respect to the Mach number:", ioData.sa.sensoutput, *this->X, U);

    fsaPrintTextOnScreen("\n ***** Derivatives with respect to the Mach number were computed! \n");

    step = step + 1;
  }

  if (ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) {

    //dXdS = 0.0;
    //dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 1.0;
    DFSPAR[2] = 0.0;
    actvar = 3;

    if (!ioData.sa.angleRad) 
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsaComputeSensitivities(ioData, "Derivatives with respect to the angle of attack:", ioData.sa.sensoutput, *this->X, U);

    fsaPrintTextOnScreen("\n ***** Derivatives with respect to the angle of attack were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
  }

  if (ioData.sa.sensBeta == SensitivityAnalysis::ON_SENSITIVITYBETA) {

    //dXdS = 0.0;
    //dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 1.0;
    actvar = 4;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsaComputeSensitivities(ioData, "Derivatives with respect to the yaw angle:", ioData.sa.sensoutput, *this->X, U);

    fsaPrintTextOnScreen("\n ***** Derivatives with respect to the yaw angle were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
  }
  

  this->output->closeAsciiFiles();

  this->com->barrier();
  MyLocalTimer += this->timer->getTime();

  if (this->com->cpuNum() == 0)
  {
    std::cout << "\n *** FluidSensityAnalysisHandler::fsaHandler >> Exit";
    std::cout << " (" << MyLocalTimer << " s)";
    std::cout << "\n\n";
  }

  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaComputeDerivativesOfFluxAndSolution
(
  IoData &ioData, 
  DistSVec<double,3> &X, 
  DistVec<double> &A, 
  DistSVec<double,dim> &U
)
{

  dFdS = 0.0;

  // Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL ){
    //fsaAnalytical(ioData, X, A, U, dFdS);
  }else{
    fsaSemiAnalytical(ioData, X, A, U, dFdS);
  }

  // Computing the derivative of the fluid variables 
  // with respect to the fsaimization variables
  fsaLinearSolver(ioData, dFdS, dUdS);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidSensitivityAnalysisHandler<dim>::fsaComputeSensitivities
(
 IoData &ioData, 
 const char *mesage, 
 const char *fileName, 
 DistSVec<double,3> &X, 
 DistSVec<double,dim> &U
)
{


  // Computing efforts
  Vec3D F, M;
  F = 0.0;
  M = 0.0;
  //fsaGetEfforts(ioData, X, U, F, M);

  // Computing derivative of the efforts
  Vec3D dFds, dMds;
  
  dFds = 0.0;
  dMds = 0.0;

  /*
  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE )
    fsaGetDerivativeOfEffortsFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, dFds, dMds);
  else
    //fsaGetDerivativeOfEffortsAnalytical(ioData, X, dXdS, U, dUdS, dFds, dMds);
  */

  if ((!ioData.sa.angleRad) && (DFSPAR[1] || DFSPAR[2])) {
    dFds *= acos(-1.0) / 180.0;
    dMds *= acos(-1.0) / 180.0;
  }  

  if (this->com->cpuNum() == 0) {
    outFile = fopen(fileName,"a+");
    if (outFile) {
      this->com->fprintf(outFile,mesage);
      this->com->fprintf(outFile,"\n");
      this->com->fprintf(outFile,"Fx= %16.13e \n",F[0]);
      this->com->fprintf(outFile,"Fy= %16.13e \n",F[1]);
      this->com->fprintf(outFile,"Fz= %16.13e \n",F[2]);
      this->com->fprintf(outFile,"dFx/ds= %16.13e \n",dFds[0]);
      this->com->fprintf(outFile,"dFy/ds= %16.13e \n",dFds[1]);
      this->com->fprintf(outFile,"dFz/ds= %16.13e \n",dFds[2]);
      this->com->fprintf(outFile,"Mx= %16.13e \n",M[0]);
      this->com->fprintf(outFile,"My= %16.13e \n",M[1]);
      this->com->fprintf(outFile,"Mz= %16.13e \n",M[2]);
      this->com->fprintf(outFile,"dMx/ds= %16.13e \n",dMds[0]);
      this->com->fprintf(outFile,"dMy/ds= %16.13e \n",dMds[1]);
      this->com->fprintf(outFile,"dMz/ds= %16.13e \n",dMds[2]);
      this->com->fprintf(outFile,"\n");
      fclose(outFile);
    }
  }
    
  double sboom = 0.0;
  double dSboom = 0.0;

  // This function is simply writing to the disk.
  this->output->writeDerivativeOfForcesToDisk(step, actvar, F, dFds, M, dMds, sboom, dSboom);
  
  //
  // This function is writing to the disk quantities of interest in binary files.
  // The possible quantities of interest include
  // - dUdS
  // - Derivative of Scalar Quantities: Density, Mach, Pressure, Temperature, TotPressure,
  //   NutTurb, EddyViscosity, VelocityScalar
  // - Derivative of Vector Quantities: VelocityVector, Displacement
  dXdS = 0.0;
  this->output->writeBinaryDerivativeOfVectorsToDisk(step+1, actvar, DFSPAR, *this->X, dXdS, U, dUdS, this->timeState);

}

//------------------------------------------------------------------------------
