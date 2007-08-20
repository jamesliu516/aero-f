#include <PostFcn.h>

#include <WallFcn.h>
#include <Vector3D.h>

#include <stdlib.h>
#include <stdio.h>

const double PostFcnEuler::third = 1.0/3.0;
//-----------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	only in the way the constructor creates the PostFcnEuler class
//	and to compute certain values for NS using lame coefficients (lambda and mu)
//------------------------------------------------------------------------------

PostFcn::PostFcn(VarFcn *vf)
{

  varFcn = vf;

}

//------------------------------------------------------------------------------

double PostFcn::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  fprintf(stderr, "*** Warning: computeNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcn::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *X, double *dX, double *V, double *dV, double phi)
{

  fprintf(stderr, "*** Warning: computeDerivativeOfNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcn::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					  Vec3D& n, double d2w[3], double* Vwall, 
					  double* Vface[3], double* Vtet[4])
{

  fprintf(stderr, "*** Warning: computeFaceScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

PostFcnEuler::PostFcnEuler(IoData &iod, VarFcn *vf) : PostFcn(vf)
{

  if (iod.eqs.fluidModel.fluid == FluidModelData::GAS){
    mach = iod.ref.mach;
    pinfty = iod.bc.inlet.pressure;

// Included (MB)
    dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

    if (iod.problem.mode == ProblemData::DIMENSIONAL) {
      dimFlag = true;
    }
    else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
      dimFlag = false;
      if (iod.sa.apressFlag == false) {
        if (iod.sa.pressFlag == false) {
          dPin = dpinfty;
        }
        else {
          dPin = 0.0;
        }      
      }
      else {
        dPin = 0.0;
      }
    }

  } else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID){
    mach = iod.ref.mach;
    double P = vf->getPrefWater();
    double a = vf->getAlphaWater();
    double b = vf->getBetaWater();
    pinfty = (P+a*pow(iod.bc.inlet.density, b));
  }
  depth = iod.bc.hydro.depth;
  gravity = iod.bc.hydro.gravity;
  alpha = iod.bc.hydro.alpha;
  beta  = iod.bc.hydro.beta;
  nGravity[0] = cos(alpha)*cos(beta);
  nGravity[1] = cos(alpha)*sin(beta);
  nGravity[2] = sin(alpha);

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::rstVar(IoData &iod, Communicator *com)
{

  mach = iod.ref.mach;
  pinfty = iod.bc.inlet.pressure;
  dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

  if (iod.problem.mode == ProblemData::DIMENSIONAL) {
    dimFlag = true;
  }
  else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
    dimFlag = false;
    if (iod.sa.apressFlag == false) {
      if (iod.sa.pressFlag == false) {
        dPin = dpinfty;
      }
      else {
        dPin = 0.0;
      }      
    }
    else {
      dPin = 0.0;
    }
  }
  
}

//------------------------------------------------------------------------------

double PostFcnEuler::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  double q = 0.0;
  double n[3];

  if (type == DENSITY)
    q = varFcn->getDensity(V);
  else if (type == MACH)
    q = varFcn->computeMachNumber(V, phi);
  else if (type == WTMACH)
    q = varFcn->computeWtMachNumber(V, phi);
  else if (type == SPEED)
    q = sqrt(varFcn->computeU2(V));
  else if (type == WTSPEED)
    q = sqrt(varFcn->computeWtU2(V));
  else if (type == PRESSURE)
    q = varFcn->getPressure(V, phi);
  else if (type == DIFFPRESSURE)
    q = varFcn->getPressure(V, phi)-pinfty;
  else if (type == TEMPERATURE)
    q = varFcn->computeTemperature(V, phi);
  else if (type == TOTPRESSURE)
    q = varFcn->computeTotalPressure(mach, V, phi);
  else if (type == NUT_TURB)
    q = varFcn->getTurbulentNuTilde(V);
  else if (type == K_TURB)
    q = varFcn->getTurbulentKineticEnergy(V);
  else if (type == EPS_TURB)
    q = varFcn->getTurbulentDissipationRate(V);
  else if(type == HYDROSTATICPRESSURE)
    q = V[0]*gravity*(depth + nGravity[0]*X[0] + nGravity[1]*X[1] + nGravity[2]*X[2]);
  else if(type == HYDRODYNAMICPRESSURE)
    q = varFcn->getPressure(V, phi) - V[0]*gravity*(depth + nGravity[0]*X[0] + nGravity[1]*X[1] + nGravity[2]*X[2]);
  else if(type == PHILEVEL)
    q = phi;

// Included (MB)
  else if (type == VELOCITY_NORM)
    q = varFcn->getVelocityNorm(V);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double q = 0.0;

  if (type == DERIVATIVE_DENSITY)
    q = varFcn->getDensity(dV);
  else if (type == DERIVATIVE_MACH)
    q = varFcn->computeDerivativeOfMachNumber(V, dV, dS[0]);
  else if (type == DERIVATIVE_PRESSURE)
    q = varFcn->getPressure(dV);
  else if (type == DERIVATIVE_TEMPERATURE)
    q = varFcn->computeDerivativeOfTemperature(V, dV);
  else if (type == DERIVATIVE_TOTPRESSURE)
    q = varFcn->computeDerivativeOfTotalPressure(mach, dS[0], V, dV, dS[0]);
  else if (type == DERIVATIVE_NUT_TURB)
    q = varFcn->getTurbulentNuTilde(dV);
  else if (type == DERIVATIVE_VELOCITY_SCALAR)
    q = varFcn->getDerivativeOfVelocityNorm(V, dV);

  return q;

}

//------------------------------------------------------------------------------

void PostFcnEuler::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
				double *Vwall, double *Vface[3], double *Vtet[4], 
		double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv,double *nodalForceWeight, int hydro)
{

  double pcg[3], p[3];
  double pcgin;
  int i;

  if (hydro == 0) {
    for(i=0;i<3;i++)
    pcg[i] = varFcn->getPressure(Vface[i]);
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++)
	pcg[i] = Vface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth);
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i] = varFcn->getPressure(Vface[i]);
      pcg[i] = pcg[i] - Vface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  double w1,w2; // w1 is the main weight
  w1 = nodalForceWeight[0]; 
  w2 = nodalForceWeight[1];
  double totWeight = w1+w2+w2;

  Fi0 = third*( (w1*p[0]+w2*p[1]+w2*p[2])/totWeight )*n;
  Fi1 = third*( (w1*p[1]+w2*p[2]+w2*p[0])/totWeight )*n;
  Fi2 = third*( (w1*p[2]+w2*p[0]+w2*p[1])/totWeight )*n;
  Fv = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
					    double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],  
					    double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double *nodalForceWeight, int hydro)
{

  double pcg[3], p[3];
  double dPcg[3], dP[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]);
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++) {
	pcg[i] = Vface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth);
	dPcg[i] = dVface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth) + Vface[i][0]*gravity * (nGravity[0]*dXface[i][0]+nGravity[1]*dXface[i][1]+nGravity[2]*dXface[i][2]);
     }
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i] = varFcn->getPressure(Vface[i]);
      pcg[i] = pcg[i] - Vface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth);
      dPcg[i] = varFcn->getPressure(dVface[i]);
      dPcg[i] = dPcg[i] - dVface[i][0]*gravity * (nGravity[0]*Xface[i][0]+nGravity[1]*Xface[i][1]+nGravity[2]*Xface[i][2] + depth) - Vface[i][0]*gravity * (nGravity[0]*dXface[i][0]+nGravity[1]*dXface[i][1]+nGravity[2]*dXface[i][2]);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  if (pin) {
    if (dimFlag)
      dPcgin = (-2.0*(*pin)/mach) * dS[0];
    else
      dPcgin = dPin * dS[0];
  }
  else {
    if (hydro == 0)
      dPcgin = dPinfty;
    else
      dPcgin = 0.0;
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  dP[0] = (dPcg[0] - dPcgin) ;
  dP[1] = (dPcg[1] - dPcgin) ;
  dP[2] = (dPcg[2] - dPcgin) ;

  double w1,w2; // w1 is the main weight
  w1 = nodalForceWeight[0]; 
  w2 = nodalForceWeight[1];
  double totWeight = w1+w2+w2;

  dFi0 = third*( (w1*dP[0]+w2*dP[1]+w2*dP[2])/totWeight )*n + third*( (w1*p[0]+w2*p[1]+w2*p[2])/totWeight )*dn;
  dFi1 = third*( (w1*dP[1]+w2*dP[2]+w2*dP[0])/totWeight )*n + third*( (w1*p[1]+w2*p[2]+w2*p[0])/totWeight )*dn;
  dFi2 = third*( (w1*dP[2]+w2*dP[0]+w2*dP[1])/totWeight )*n + third*( (w1*p[2]+w2*p[0]+w2*p[1])/totWeight )*dn;
  dFv = 0.0;

}

//------------------------------------------------------------------------------

double PostFcnEuler::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				      double* Vwall, double* Vface[3], double* Vtet[4])
{

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				      double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcnEuler::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
					  double d2w[3], double* Vwall, double* Vface[3], 
					  double* Vtet[4], double pin)
{

  double p = third * ( varFcn->getPressure(Vface[0]) + varFcn->getPressure(Vface[1]) +
		       varFcn->getPressure(Vface[2]) ) - pin;
  double W = - ndot * p;

  return W;

}

//------------------------------------------------------------------------------

PostFcnNS::PostFcnNS(IoData &iod, VarFcn *vf) 
  : PostFcnEuler(iod, vf), NavierStokesTerm(iod, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, PostFcn::varFcn, viscoFcn);
  else
    wallFcn = 0;

}

//------------------------------------------------------------------------------

PostFcnNS::~PostFcnNS()
{

  if (wallFcn) delete wallFcn;

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnNS::rstVar(IoData &iod, Communicator *com)
{

  PostFcnEuler::rstVar(iod, com);

  NavierStokesTerm::rstVarNS(iod, com);
  
  if (wallFcn)
    wallFcn->rstVar(iod, com);

}

//------------------------------------------------------------------------------

double PostFcnNS::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					    Vec3D& n, double d2w[3], double* Vwall, 
					    double* Vface[3], double* Vtet[4])
{

  double q = 0.0;

  if (type == DELTA_PLUS) {
#if defined(HEAT_FLUX)
    q = computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet) / sqrt(n*n);
#elif defined(SKIN_FRICTION)
    Vec3D t(1.0, 0.0, 0.0);
    Vec3D F = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
    q = 2.0 * t * F / sqrt(n*n);
#else
    if (wallFcn)
      q = wallFcn->computeDeltaPlus(n, d2w, Vwall, Vface);
    else
      fprintf(stderr, "*** Warning: yplus computation not implemented\n");
#endif
  }

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double q = 0.0;

  q = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return q;

}

//------------------------------------------------------------------------------

Vec3D PostFcnNS::computeViscousForce(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				     double* Vwall, double* Vface[3], double* Vtet[4])
{

  Vec3D Fv;

  if (wallFcn)
    Fv = wallFcn->computeForce(n, d2w, Vwall, Vface);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    Fv[0] = tij[0][0] * n[0] + tij[0][1] * n[1] + tij[0][2] * n[2];
    Fv[1] = tij[1][0] * n[0] + tij[1][1] * n[1] + tij[1][2] * n[2];
    Fv[2] = tij[2][0] * n[0] + tij[2][1] * n[1] + tij[2][2] * n[2];
  }

  return -1.0 * Fv;

}

//------------------------------------------------------------------------------

// Included (MB)
Vec3D PostFcnNS::computeDerivativeOfViscousForce(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3],
				     double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  Vec3D dFv;

  if (wallFcn)
    dFv = wallFcn->computeDerivativeOfForce(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double du[4][3], ducg[3];
    computeDerivativeOfVelocity(dVtet, du, ducg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double ddudxj[3][3];
    computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];

    double dooreynolds_lambda = -1.0 / ( reynolds_lambdaNS * reynolds_lambdaNS ) * dRe_lambdadMachNS * dS[0];

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);

    double dmu = dooreynolds_mu * viscoFcn->compute_mu(Tcg) + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dS[0]);

    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dS[0]);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    double dtij[3][3];
    computeDerivativeOfStressTensor(mu, dmu, lambda, dlambda, dudxj, ddudxj, dtij);

    dFv[0] = dtij[0][0] * n[0] + dtij[0][1] * n[1] + dtij[0][2] * n[2] + tij[0][0] * dn[0] + tij[0][1] * dn[1] + tij[0][2] * dn[2];
    dFv[1] = dtij[1][0] * n[0] + dtij[1][1] * n[1] + dtij[1][2] * n[2] + tij[1][0] * dn[0] + tij[1][1] * dn[1] + tij[1][2] * dn[2];
    dFv[2] = dtij[2][0] * n[0] + dtij[2][1] * n[1] + dtij[2][2] * n[2] + tij[2][0] * dn[0] + tij[2][1] * dn[1] + tij[2][2] * dn[2];

  }

  return -1.0 * dFv;

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
			     double *Vwall, double *Vface[3], double *Vtet[4], 
		    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double *nodalForceWeight, int hydro)
{

  PostFcnEuler::computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, nodalForceWeight, hydro);
  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void PostFcnNS::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                            double *Vface[3], double *dVface[3], double *Vtet[4],
                                            double *dVtet[4], double dS[3], double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv,double *nodalForceWeight, int hydro)
{

  PostFcnEuler::computeDerivativeOfForce(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, nodalForceWeight, hydro);

  dFv = computeDerivativeOfViscousForce(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

}

//------------------------------------------------------------------------------

double PostFcnNS::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				   double* Vwall, double* Vface[3], double* Vtet[4])
{

  double hp = 0.0;

  if (wallFcn)
    hp = wallFcn->computeHeatPower(n, d2w, Vwall, Vface);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);
    double qj[3];
    computeHeatFluxVector(kappa, dTdxj, qj);
    hp = qj[0]*n[0] + qj[1]*n[1] + qj[2]*n[2];
  }

  return hp;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				   double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  double dhp = 0.0;

  if (wallFcn)
    dhp = wallFcn->computeDerivativeOfHeatPower(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double ddTdxj[3];
    computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);
    double kappa = ooreynolds * thermalCondFcn->compute(Tcg);
    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];
    double dkappa = dooreynolds_mu * thermalCondFcn->compute(Tcg) + ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcg, dS[0]);
    double qj[3];
    computeHeatFluxVector(kappa, dTdxj, qj);
    double dqj[3];
    computeDerivativeOfHeatFluxVector(kappa, dkappa, dTdxj, ddTdxj, dqj);
    dhp = dqj[0]*n[0] + qj[0]*dn[0] + dqj[1]*n[1] + qj[1]*dn[1] + dqj[2]*n[2] + qj[2]*dn[2];
  }

  return dhp;

}

//------------------------------------------------------------------------------

double PostFcnNS::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
				       double d2w[3], double* Vwall, double* Vface[3], 
				       double* Vtet[4], double pin)
{

  double W = PostFcnEuler::computeInterfaceWork(dp1dxj, n, ndot, d2w, Vwall, Vface, Vtet, pin);

  if (wallFcn)
    W += wallFcn->computeInterfaceWork(n, d2w, Vwall, Vface);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    W += (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * n[0] +
      (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * n[1] +
      (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * n[2];
  }

  return W;

}

//------------------------------------------------------------------------------

PostFcnSA::PostFcnSA(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), SATerm(iod)
{

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnSA::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarSA(iod);

}

//------------------------------------------------------------------------------

double PostFcnSA::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnSA::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------

PostFcnDES::PostFcnDES(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), DESTerm(iod)
{

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnDES::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarDES(iod);

}

//------------------------------------------------------------------------------
                                                                                                
double PostFcnDES::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnDES::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------

PostFcnKE::PostFcnKE(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), KEpsilonTerm(iod)
{

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnKE::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarKE(iod);

}

//------------------------------------------------------------------------------

double PostFcnKE::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY)
    q = computeTurbulentViscosity(V);
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnKE::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY)
    dq = computeDerivativeOfTurbulentViscosity(V, dV, dS[0]);
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------
