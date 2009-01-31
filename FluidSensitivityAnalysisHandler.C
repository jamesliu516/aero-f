#include <FluidSensitivityAnalysisHandler.h>

#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistVector.h>
#include <MeshMotionSolver.h>
#include <StructExc.h>
#include <TimeData.h>
#include <MatVecProd.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MemoryPool.h>

#include <math.h>
#include <iomanip.h>

#include <string>

//------------------------------------------------------------------------------

template<int dim>
FluidSensitivityAnalysisHandler<dim>::FluidSensitivityAnalysisHandler(IoData &ioData, GeoSource &geoSource, Domain *dom) :
ImplicitCoupledTsDesc<dim>(ioData, geoSource, dom), domain(dom),
dXdS(dom->getNodeDistInfo()), dXdSb(dom->getNodeDistInfo()), Xc(dom->getNodeDistInfo()), 
dAdS(dom->getNodeDistInfo()), dFdS(dom->getNodeDistInfo()), dUdS(dom->getNodeDistInfo()), 
p(dom->getNodeDistInfo()), dPdS(dom->getNodeDistInfo()), Flux(dom->getNodeDistInfo()), 
FluxFD(dom->getNodeDistInfo()), Pin(dom->getFaceDistInfo()), Uc(dom->getNodeDistInfo())
// Tests
,Xplus(dom->getNodeDistInfo()), Xminus(dom->getNodeDistInfo()), dX(dom->getNodeDistInfo())
{

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
    Xp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Xm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Ap = new DistVec<double>(domain->getNodeDistInfo());
    Am = new DistVec<double>(domain->getNodeDistInfo());
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Um = new DistSVec<double,dim>(domain->getNodeDistInfo());
  }  
  else if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    Xp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Xm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lp = 0;
    Lm = 0;
    Ap = new DistVec<double>(domain->getNodeDistInfo());
    Am = new DistVec<double>(domain->getNodeDistInfo());
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = 0;
    Um = 0;
   }
  else {
    Xp = 0; 
    Xm = 0;
    Lp = 0;
    Lm = 0;
    Ap = 0;
    Am = 0;
    Fp = 0;
    Fm = 0;
  }

  Z = new DistSVec<double,3>(domain->getNodeDistInfo());

  timeDataOpt = this->timeState->getDataOpt();

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)  {
    if (ioData.sa.mvp == SensitivityAnalysis::FD)
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp, domain, ioData, true);
    else if (ioData.sa.mvp == SensitivityAnalysis::H2)
      mvp = new MatVecProdH2<MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, domain, this->geoState);
  }
  else {
    if (ioData.sa.mvp == SensitivityAnalysis::FD)
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit,  0, this->geoState, this->spaceOp, domain, ioData, true);
    else if (ioData.sa.mvp == SensitivityAnalysis::H2)
      mvp = new MatVecProdH2<MatScalar,dim>(ioData, this->varFcn, 0, this->spaceOp, domain, this->geoState);
  }

  pc = ImplicitTsDesc<dim>::template
    createPreconditioner<PrecScalar,dim>(ioData.sa.ksp.ns.pc, domain);

  ksp = createKrylovSolver(this->getVecInfo(), ioData.sa.ksp.ns, mvp, pc, this->com);

  MemoryPool mp;

  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);

  mms = new TetMeshMotionSolver(ioData.dmesh, geoSource.getMatchNodes(), domain, 0);

  length = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = domain->getNumLocSub();

  dXdS=0.0;
  dFdS=0.0;
  dUdS=0.0;
  p=0.0;
  dPdS=0.0;

  FluxFD = 0.0;
  Flux = 0.0;

  reynolds0 = ioData.ref.reynolds_mu;
  kenergy0 = ioData.bc.inlet.kenergy;

}

//------------------------------------------------------------------------------

template<int dim>
FluidSensitivityAnalysisHandler<dim>::~FluidSensitivityAnalysisHandler()
{

  if (mms) delete mms;

  if (timeDataOpt) delete timeDataOpt;

  if (mvp) delete mvp;

  if (pc) delete pc;

  if (ksp) delete ksp;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaRestartBcFluxs(IoData &ioData)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

// Remark: For internal flows the SA using inlet or outlet Mach, Alpha and Beta should be specified in the input file

  ioData.bc.inlet.mach = xmach;
  ioData.bc.outlet.mach = xmach;

  ioData.bc.inlet.alpha = alprad;
  ioData.bc.outlet.alpha = alprad;

  ioData.bc.inlet.beta = teta;
  ioData.bc.outlet.beta = teta;

//  this->com->fprintf(stderr, "\n\n FluidSensitivityAnalysis values: \n\n");
//  this->com->fprintf(stderr, "\n\n Inlet Alpha = %20.17e \n\n",ioData.bc.inlet.alpha);
//  this->com->fprintf(stderr, "\n\n Inlet Beta = %20.17e \n\n",ioData.bc.inlet.beta);
//  this->com->fprintf(stderr, "\n\n Outlet Alpha = %20.17e \n\n",ioData.bc.outlet.alpha);
//  this->com->fprintf(stderr, "\n\n Outlet Beta = %20.17e \n\n",ioData.bc.outlet.beta);
      
  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {

    ioData.ref.mach = xmach;
    ioData.ref.dRe_mudMach = 0.0;
    ioData.ref.dRe_lambdadMach = 0.0;

    if (ioData.sa.densFlag == false) {
      ioData.bc.inlet.density = 1.0;
      ioData.bc.outlet.density = 1.0;
      this->com->fprintf(stderr, "\n\n NonDim Density = %e \n\n",ioData.bc.inlet.density);
    }

    if (ioData.sa.pressFlag == false) {
      ioData.bc.inlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);
      ioData.bc.outlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);

      this->com->fprintf(stderr, "\n\n NonDim Pressure = %e \n\n",ioData.bc.inlet.pressure);
    }

    if (ioData.sa.apressFlag == false) {
      if (ioData.sa.pressFlag == false) {
        ioData.aero.pressure = ioData.bc.inlet.pressure;
        Pin = ioData.aero.pressure;
//        this->com->fprintf(stderr, "\n\n NonDim Internal Pressure = %e \n\n",ioData.aero.pressure);

        this->postOp->rstVarPostFcn(ioData);
        this->postOp->rstVar(ioData);

        this->spaceOp->rstFluxFcn(ioData);
 
        mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

        this->rstVarImplicitCoupledTsDesc(ioData);
      }
    }
  }
  else if (ioData.problem.mode == ProblemData::DIMENSIONAL) {

    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;
    ioData.eqs.fluidModel2.pmin *= ioData.ref.rv.pressure;

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
      ioData.rmesh.vpts[j]->time     *= ioData.ref.rv.time;
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

    ioData.ref.mach = ioData.bc.inlet.mach;
    ioData.ref.density = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;
    double velocity = ioData.ref.mach * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);
    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);
    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
      (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);
    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;
    ioData.ref.reynolds_lambda = -3.0 * ioData.ref.reynolds_mu/2.0;

    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->com->fprintf(stderr, "\n\n Reynolds = %e \n\n",ioData.ref.reynolds_mu);

    double dvelocitydMach = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;
    ioData.ref.dRe_lambdadMach = -3.0 * ioData.ref.dRe_mudMach/2.0;

    ioData.ref.rv.mode = RefVal::DIMENSIONAL;
    ioData.ref.rv.density = ioData.ref.density;
    ioData.ref.rv.velocity = velocity;
    ioData.ref.rv.pressure = ioData.ref.density * velocity*velocity;
    ioData.ref.rv.temperature = gamma*(gamma - 1.0) * ioData.ref.mach*ioData.ref.mach * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density);
    ioData.ref.rv.viscosity_mu = viscosity;
    ioData.ref.rv.viscosity_lambda = -2.0/3.0 * viscosity;
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
    ioData.eqs.fluidModel2.pmin /= ioData.ref.rv.pressure;

    ioData.bc.inlet.density /= ioData.ref.rv.density;
    ioData.bc.inlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps /= ioData.ref.rv.epsilon;
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

    this->rstVarImplicitCoupledTsDesc(ioData);
    
    this->bcData->rstVar(ioData);

    timeDataOpt->rstVar(ioData);

    this->timeState->rstVar(ioData);

    this->output->rstVar(ioData);

    this->restart->rstVar(ioData);

    this->data->rstVar(ioData);

    this->refVal->rstVar(ioData);
    
    this->varFcn->rstVar(ioData);

  }

// initialize boundary condition fluxes
  this->bcData->initialize(ioData, this->varFcn, *this->X);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaGetEfforts(IoData &ioData, 
                                  DistSVec<double,3> &X, DistSVec<double,dim> &U, Vec3D &F, Vec3D &M)
{

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, force, moment;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    F *= 2.0 * this->refVal->length*this->refVal->length / surface;
    M *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
  }
  else {
    F *= this->refVal->force;
    M *= this->refVal->energy;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaGetDerivativeOfEffortsFiniteDifference(IoData &ioData,
                                                          DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &A,
                                                          DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
                                                          Vec3D &dForces, Vec3D &dMoments)
{

//Remark: Error mesage for pointers
  if (Up == 0) {
    fprintf(stderr, "*** Error: Varible Fp does not exist!\n");
    exit(1);
  }
  if (Um == 0) {
    fprintf(stderr, "*** Error: Varible Fm does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Varible Xp does not exist!\n");
    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Varible Xm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Varible Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {
    fprintf(stderr, "*** Error: Varible Am does not exist!\n");
    exit(1);
  }

  double eps=ioData.sa.eps;
  double xmachc;
  double alpradc;
  double tetac;
  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, F, Fplus, Fminus, dF, M, Mp, Mm, dM;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec3D *Fip = new Vec3D[nSurfs];
  Vec3D *Mip = new Vec3D[nSurfs];
  Vec3D *Fvp = new Vec3D[nSurfs];
  Vec3D *Mvp = new Vec3D[nSurfs];

  Vec3D *Fim = new Vec3D[nSurfs];
  Vec3D *Mim = new Vec3D[nSurfs];
  Vec3D *Fvm = new Vec3D[nSurfs];
  Vec3D *Mvm = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  *Up = 0.0;
  *Um = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  this->spaceOp->computeGradP(X, A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  Xc=X;

  *Xp=X+eps*dX;

  X=*Xp;

  *Up=U+eps*dU;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  this->spaceOp->computeGradP(*Xp, *Ap, *Up);

  this->postOp->computeForceAndMoment(x0, *Xp, *Up, 0, Fip, Mip, Fvp, Mvp);

  Fplus = 0.0;
  Mp = 0.0;

  Fplus = Fip[0] + Fvp[0];
  Mp = Mip[0] + Mvp[0];

  X=Xc;

  *Xm=X-eps*dX;

  X=*Xm;

  *Um=U-eps*dU;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  this->spaceOp->computeGradP(*Xm, *Am, *Um);

  this->postOp->computeForceAndMoment(x0, *Xm, *Um, 0, Fim, Mim, Fvm, Mvm);

  Fminus = 0.0;
  Mm = 0.0;

  Fminus = Fim[0] + Fvm[0];
  Mm = Mim[0] + Mvm[0];

  dF=1.0/(2.0*eps)*(Fplus-Fminus);
  dM=1.0/(2.0*eps)*(Mp-Mm);

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  this->spaceOp->computeGradP(X, A, U);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
    dForces=dF;
    dMoments=dM;
  }
  else {
    dF *= this->refVal->force;
    dM *= this->refVal->energy;
    F*=dForce;
    M*=dEnergy;
    dForces=dF+F;
    dMoments=dM+M;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaGetDerivativeOfEffortsAnalytical(IoData &ioData,
                                                          DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                                          DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
                                                          Vec3D &dForces, Vec3D &dMoments)
{


  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, F, dF, M, dM;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec3D *dFi = new Vec3D[nSurfs];
  Vec3D *dMi = new Vec3D[nSurfs];
  Vec3D *dFv = new Vec3D[nSurfs];
  Vec3D *dMv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);

  dF = 0.0;
  dM = 0.0;

  dF = dFi[0] + dFv[0];
  dM = dMi[0] + dMv[0];

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
    dForces=dF;
    dMoments=dM;
  }
  else {
    dF *= this->refVal->force;  
    dM *= this->refVal->energy;
    F *= dForce;
    M *= dEnergy;
    dForces = dF+F;
    dMoments = dM+M;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaGetDerivativeOfLoadFiniteDifference(IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &A,
                                                                                                                            DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
{

//Remark: Error mesage for pointers
  if (Up == 0) {
    fprintf(stderr, "*** Error: Varible Up does not exist!\n");
    exit(1);
  }
  if (Um == 0) {
    fprintf(stderr, "*** Error: Varible Um does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Varible Xp does not exist!\n");
    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Varible Xm does not exist!\n");
    exit(1);
  }
  if (Lp == 0) {
    fprintf(stderr, "*** Error: Varible Lp does not exist!\n");
    exit(1);
  }
  if (Lm == 0) {
    fprintf(stderr, "*** Error: Varible Lm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Varible Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {

    fprintf(stderr, "*** Error: Varible Am does not exist!\n");
    exit(1);
  }

  double eps=ioData.sa.eps;

  double xmachc;
  double alpradc;
  double tetac;
  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  *Up = 0.0;
  *Um = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Lp = 0.0;
  *Lm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  load=0.0;
  dLoad=0.0;

  this->spaceOp->computeGradP(X, A, U);

  this->postOp->computeNodalForce(X, U, Pin, load);

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  Xc=X;

  *Xp=X+eps*dX;

  X=*Xp;

  *Up=U+eps*dU;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  this->spaceOp->computeGradP(*Xp, *Ap, *Up);

  this->postOp->computeNodalForce(*Xp, *Up, Pin, *Lp);

  X=Xc;

  *Xm=X-eps*dX;

  X=*Xm;

  *Um=U-eps*dU;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  this->spaceOp->computeGradP(*Xm, *Am, *Um);

  this->postOp->computeNodalForce(*Xm, *Um, Pin, *Lm);

  dLoad=1.0/(2.0*eps)*((*Lp)-(*Lm));

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  this->spaceOp->computeGradP(X, A, U);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaGetDerivativeOfLoadAnalytical(IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX, 
                                                                                                                  DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  load=0.0;
  dLoad=0.0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeNodalForce(X, U, Pin, load);

  this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  this->postOp->computeDerivativeOfNodalForce(X, dX, U, dU, Pin, DFSPAR, dLoad);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaSemiAnalytical(IoData &ioData, DistSVec<double,3> &X,
                                                                                         DistVec<double> &A,


                                                                                         DistSVec<double,dim> &U,
                                                                                         DistSVec<double,dim> &dF)
{

//Remark: Error mesage for pointers
  if (Fp == 0) {
    fprintf(stderr, "*** Error: Varible Fp does not exist!\n");
    exit(1);
  }
  if (Fm == 0) {
    fprintf(stderr, "*** Error: Varible Fm does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Varible Xp does not exist!\n");

    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Varible Xm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Varible Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {
    fprintf(stderr, "*** Error: Varible Am does not exist!\n");
    exit(1);
  }

  double dtLeft;
  double eps=ioData.sa.eps;
  double xmachc;
  double alpradc;
  double tetac;

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  *Fp = 0.0;
  *Fm = 0.0;
   dF = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  Xc=X;

  *Xp=X+eps*dXdS;

  X=*Xp;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  dtLeft = 0.0;
  computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(*Xp, *Ap, U, *Fp, this->timeState);

  this->spaceOp->applyBCsToResidual(U, *Fp);

  X=Xc;

  *Xm=X-eps*dXdS;

  X=*Xm;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  dtLeft = 0.0;
  computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(*Xm, *Am, U, *Fm, this->timeState);

  this->spaceOp->applyBCsToResidual(U, *Fm);

  dAdS=1.0/(2.0*eps)*((*Ap)-(*Am));

  dF=1.0/(2.0*eps)*((*Fp)-(*Fm));

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  dtLeft = 0.0;
  computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeGradP(X, A, U);

  if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    this->geoState->updateConfigSA();
    this->bcData->initializeSA(ioData, this->varFcn, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);
  }

}


//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaAnalytical(IoData &ioData, DistSVec<double,3> &X,
                                             DistVec<double> &A,
                                             DistSVec<double,dim> &U,
                                             DistSVec<double,dim> &dFdS)
{
 
// Computing the normal, derivative of the normal and of the control volume
    this->geoState->computeDerivatives(X, dXdS, this->bcData->getVelocityVector(), this->bcData->getDerivativeOfVelocityVector(), dAdS);

// Computing the derivatives of the boundary fluxs
    this->bcData->initializeSA(ioData, this->varFcn, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);

// Computing the partial derivative of the flux with respect to the fsaimization variables
    this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS, this->timeState);

    this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);
    
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaSetUpLinearSolver(IoData &ioData, DistSVec<double,3> &X,
                                                                                     DistVec<double> &A,
                                                                                     DistSVec<double,dim> &U,
                                                                                     DistSVec<double,dim> &dFdS)
{

// Prepraring the linear solver

  fsaRestartBcFluxs(ioData);

  this->geoState->reset(X);

  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);

  this->bcData->update(X);

  this->spaceOp->computeResidual(X, A, U, FluxFD, this->timeState);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
    this->timeState->add_dAW_dt(1, *this->geoState, A, U, FluxFD);

  this->spaceOp->applyBCsToResidual(U, FluxFD);

  mvp->evaluate(0, X, A, U, FluxFD);

  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);

  MatVecProdFD<dim,dim> *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
  MatVecProdH2<MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<MatScalar,dim> *>(mvp);

  if (mvpfd || mvph2) {
    this->spaceOp->computeJacobian(X, A, U, *_pc, this->timeState);

    if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
      this->timeState->addToJacobian(A, *_pc, U);

    this->spaceOp->applyBCsToJacobian(U, *_pc);
  }

  pc->setup();

// Computing flux for compatibility correction of the derivative of the flux   
  this->spaceOp->computeResidual(X, A, U, Flux, this->timeState, false);

}

//------------------------------------------------------------------------------


template<int dim>
int FluidSensitivityAnalysisHandler<dim>::fsaLinearSolver(IoData &ioData, DistVec<double> &A, DistSVec<double,3> &X, DistSVec<double,dim> &U, DistSVec<double,dim> &dFdS, DistSVec<double,dim> &dUdS)
{

  int numberIteration;

  dFdS *= (-1.0);

  ksp->setup(0, 0, dFdS);

  numberIteration = ksp->solve(dFdS, dUdS);

  if (ioData.sa.excsol) {
    if (numberIteration < ioData.sa.ksp.ns.maxIts)
      return 1;
    else
      return 0;
  }
  else
    return 1;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaOutput1D(const char *fileName, DistVec<double> &V)
{

outFile = fopen(fileName,"w");


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    for (int j=0; j<V(iSub).size(); ++j)
        fprintf(outFile," %20.17e \n",V(iSub)[j]);
//        fprintf(outFile," %9.6e \n",V(iSub)[j]);


fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaOutput3D(const char *fileName, DistSVec<double,3> &V)
{

outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<3; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
//        fprintf(outFile," %9.6e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");

    }
  }

fclose(outFile);

}

//------------------------------------------------------------------------------


template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaOutputDimD(const char *fileName, DistSVec<double,dim> &V)
{


outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub ; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<dim; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
//        fprintf(outFile," %9.6e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");
    }

  }

fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
int FluidSensitivityAnalysisHandler<dim>::fsaHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR(1)  -  Mach number differential
  // DFSPAR(2)  -  angle of attack differential
  // DFSPAR(3)  -  yaw angle differential

  this->output->openAsciiFiles();

// Reseting the configuration control of the geometry datas
  this->geoState->resetConfigSA();

  if (this->com->cpuNum() == 0) {
    outFile = fopen(ioData.sa.sensoutput,"w");
    if (outFile) 
      fclose(outFile);
  }

  xmach = ioData.sa.machref;
  alprad = ioData.sa.alpharef;
  teta = ioData.sa.betaref;

  double dtLeft = 0.0;
  computeTimeStep(1, &dtLeft, U);

  this->computeMeshMetrics();

  updateStateVectors(U);  

// Setting up the linear solver
  fsaSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) {

    bool stepStop = false;
    double tag;

    step = ioData.sa.si;
    dXdS = 0.0;
    dXdSb = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 1;

    while (!stepStop) {

// Reading derivative of the overall deformation
      domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdSb);

// Checking if dXdSb has entries different from zero at the interior of the mesh
      this->postOp->checkVec(dXdSb);

      if (dXdSb.norm() != 0.0) {

        this->com->fprintf(stderr, "\n ***** Shape variable %d\n", step);

// Updating the mesh
        mms->optSolve(ioData, 0, dXdSb, dXdS, *this->X);

        fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);
  
        fsaComputeSensitivities(ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

        dXdS = 0.0;
        dXdSb = 0.0;

        if ((ioData.sa.sf >= 0) && (ioData.sa.sf == step)) {
          stepStop = true;
          fsaPrintTextOnScreen("\n ***** Derivatives with respect to the mesh position were computed! \n");
        }
        else {
          step = step + 1;
        }

      }
      else {

        stepStop = true;

        fsaPrintTextOnScreen("\n ***** Derivatives with respect to the mesh position were computed! \n");

      }

    }

  }

  if (ioData.sa.sensMach == SensitivityAnalysis::ON_SENSITIVITYMACH) {

    dXdS = 0.0;
    dAdS = 0.0;
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

    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 1.0;
    DFSPAR[2] = 0.0;
    actvar = 3;

    if (!ioData.sa.angleRad) 
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsaComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsaComputeSensitivities(ioData, "Derivatives with respect to the algle of attck:", ioData.sa.sensoutput, *this->X, U);

    fsaPrintTextOnScreen("\n ***** Derivatives with respect to the algle of attck were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
  }

  if (ioData.sa.sensBeta == SensitivityAnalysis::ON_SENSITIVITYBETA) {

    dXdS = 0.0;
    dAdS = 0.0;
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

  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaComputeDerivativesOfFluxAndSolution(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U)
{

  dFdS = 0.0;

  dUdS = 0.0;

// Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL )
    fsaAnalytical(ioData, X, A, U, dFdS);
  else
    fsaSemiAnalytical(ioData, X, A, U, dFdS);

// Computing the derivative of the fluid variables with respect to the fsaimization variables
  int istop = 0;
  while (!istop) {
    istop = fsaLinearSolver(ioData, A, X, U, dFdS, dUdS);
    dFdS *= (-1.0);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidSensitivityAnalysisHandler<dim>::fsaComputeSensitivities(IoData &ioData, const char *mesage, const char *fileName, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

// Computing efforts
  Vec3D F, M;
  fsaGetEfforts(ioData, X, U, F, M);

// Computing derivative of the efforts
  Vec3D dFds, dMds;

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE )
    fsaGetDerivativeOfEffortsFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, dFds, dMds);
  else
    fsaGetDerivativeOfEffortsAnalytical(ioData, X, dXdS, U, dUdS, dFds, dMds);

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

  this->output->writeDerivativeOfForcesToDisk(step, actvar, F, dFds, M, dMds, sboom, dSboom);
  
  this->output->writeBinaryDerivativeOfVectorsToDisk(step+1, actvar, DFSPAR, *this->X, dXdS, U, dUdS, this->timeState);

}

//------------------------------------------------------------------------------
