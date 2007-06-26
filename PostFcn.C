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

  refRho = iod.ref.rv.density;
  refPressure = iod.ref.rv.pressure;
  refVel = iod.ref.rv.velocity;
}

//------------------------------------------------------------------------------

double PostFcnEuler::computeNodeScalarQuantity(ScalarType type, double *V, double *X, double phi)
{

  double q = 0.0;
  double n[3];

  if (type == DENSITY)  {
    q = varFcn->getDensity(V);
    q *= refRho;
  }
  else if (type == MACH)
    q = varFcn->computeMachNumber(V, phi);
  else if (type == WTMACH)
    q = varFcn->computeWtMachNumber(V, phi);
  else if (type == SPEED)  {
    q = sqrt(varFcn->computeU2(V));
    q *= refVel;
  }
  else if (type == WTSPEED)  {
    q = sqrt(varFcn->computeWtU2(V));
    q *= refVel;
  }
  else if (type == PRESSURE)  {
    q = varFcn->getPressure(V, phi);
    q *= refPressure;
  }
  else if (type == DIFFPRESSURE)  {
    q = varFcn->getPressure(V, phi)-pinfty;
    q *= refPressure;
  }
  else if (type == TEMPERATURE)
    q = varFcn->computeTemperature(V, phi);
  else if (type == TOTPRESSURE)  {
    q = varFcn->computeTotalPressure(mach, V, phi);
    q *= refPressure;
  }
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

double PostFcnEuler::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				      double* Vwall, double* Vface[3], double* Vtet[4])
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

void PostFcnNS::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
			     double *Vwall, double *Vface[3], double *Vtet[4], 
		    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double *nodalForceWeight, int hydro)
{

  PostFcnEuler::computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, nodalForceWeight, hydro);
  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

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

PostFcnDES::PostFcnDES(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), DESTerm(iod)
{

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

PostFcnKE::PostFcnKE(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), KEpsilonTerm(iod)
{

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
