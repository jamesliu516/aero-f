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

  return q;

}

//------------------------------------------------------------------------------
void PostFcnEuler::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                                double *Vwall, double *Vface[3], double *Vtet[4],
                double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro)
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

// ##################
 Vec3D x0, x1, x2, x;
 double temp;
 for(int i = 0; i<3; i++)
 {
  x0[i] = Xface[0][i];
  x1[i] = Xface[1][i];
  x2[i] = Xface[2][i];
 }

// for node 0
  i=0; // i represents node
 x = (7.0/18.0)*(x1 + x2 - 2.0*x0);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi0 = (1.0/6.0*temp)*n;

// for node 1
 i=1;
 x = (7.0/18.0)*(x2 + x0 - 2.0*x1);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi1 = (1.0/6.0*temp)*n;

// for node 2
 i=2;
 x = (7.0/18.0)*(x0 + x1 - 2.0*x2);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi2 = (1.0/6.0*temp)*n;

 Fv = 0.0;

}
//-----------------------------------------------------------------------------------

void PostFcnEuler::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                double *Vwall, double *Vface[3], double *Vtet[4],
                double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro)
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

// ##########################

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 // Computing p_1C, p_0C and p_3C
 i = 0;
 x = x1-xC;
 p_1C = p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xC;
 p_0C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xC;
 p_3C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_1S, p_0S and p_2S
 i = 1;
 x = x1-xS;
 p_1S = p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xS;
 p_0S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x2-xS;
 p_2S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_2P, p_0P and p_3P
 i = 2;
 x = x2-xP;
 p_2P = p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xP;
 p_0P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xP;
 p_3P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0 ;
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0 ;
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0 ;

 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0 ;
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0 ;
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0 ;

 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0 ;

 Fi0 = (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi1 = (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi2 = (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);

 Fv = 0.0;
}

//--------------------------------------------------------------------------------------

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
                    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro)
{

  PostFcnEuler::computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                             double *Vwall, double *Vface[3], double *Vtet[4],
                    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro)
{

  PostFcnEuler::computeForceTransmitted(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro);

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
