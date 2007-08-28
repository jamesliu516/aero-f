#include <WallFcn.h>

#include <IoData.h>
#include <BcDef.h>
#include <VarFcn.h>
#include <ViscoFcn.h>
#include <Vector3D.h>

#include <math.h>

//------------------------------------------------------------------------------

const double WallFcn::third = 1.0 / 3.0;
const double WallFcn::eleventh = 1.0 / 11.0;

//------------------------------------------------------------------------------

WallFcn::WallFcn(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  varFcn(varf), viscoFcn(visf)
{

  gam = varFcn->getGamma();
  vkcst = 0.41;
  reynolds = iod.ref.reynolds_mu;

}

//------------------------------------------------------------------------------

Vec3D WallFcn::computeTangentVector(Vec3D &n, Vec3D &u)
{

  Vec3D un = (u * n) * n;

  Vec3D t = u - un;

  double norm = sqrt(t*t);

  if (norm != 0.0)
    t *= 1.0 / norm;

  return t;

}

//------------------------------------------------------------------------------

void WallFcn::computeFaceValues(double d2wall[3], double *Vwall, double *V[3],
				double &delta, Vec3D &du, double &dT, Vec3D &uw,
				double &rhow, double &muw)
{

  delta = third * (d2wall[0] + d2wall[1] + d2wall[2]);

  rhow = third * ( varFcn->getDensity(V[0]) + 
		   varFcn->getDensity(V[1]) +
		   varFcn->getDensity(V[2]) );

  uw[0] = Vwall[1];
  uw[1] = Vwall[2];
  uw[2] = Vwall[3];

  double Tw = Vwall[4];

  du = third * ( varFcn->getVelocity(V[0]) +
		 varFcn->getVelocity(V[1]) +
		 varFcn->getVelocity(V[2]) ) - uw;

  double T = third * ( varFcn->computeTemperature(V[0]) + 
		       varFcn->computeTemperature(V[1]) +
		       varFcn->computeTemperature(V[2]) );

  dT = T - Tw;

  muw = viscoFcn->compute_mu(T);

}

//------------------------------------------------------------------------------

double WallFcn::computeFrictionVelocity(Vec3D &t, double delta, double rho, Vec3D &u, double mu)
{

  int maxits = 20;
  double eps = 1.e-6;

  double ut = u * t;

  double utau = sqrt(mu * ut / (reynolds * rho * delta));

  int it;
  for (it=0; it<maxits; ++it) {

    double dplus = reynolds * utau * delta * rho / mu;

    if (dplus < 1.0) {
      dplus = 1.0;
      utau = mu / (reynolds * delta * rho);
      break;
    }

    double f = 2.5 * log(1.0 + vkcst*dplus) 
      + 7.8 * (1.0 - exp(-eleventh*dplus) - eleventh*dplus*exp(-0.33*dplus));

    double F = utau * f - ut;

    double res = F * F;

    double target;
    if (it == 0) 
      target = eps * res;

    if (res == 0.0 || res <= target) 
      break;

    double dfdutau = 2.5*vkcst/(1.0 + vkcst*dplus) + 7.8*eleventh* 
      (exp(-eleventh*dplus) - exp(-0.33*dplus) + 0.33*dplus*exp(-0.33*dplus));

    dfdutau *= reynolds * delta * rho / mu;

    utau -= F / (f + utau * dfdutau);

  }

  if (it == maxits) 
    fprintf(stderr, "*** Warning: Newton did not converge on utau\n");

  if (utau <= 0.0)
    fprintf(stderr, "*** Warning: utau=%e\n", utau);

  return utau;

}

//------------------------------------------------------------------------------

double WallFcn::computeFrictionTemperature(double utau, double delta, double rho, 
					   double dT, double mu)
{

  double dplus = reynolds * utau * delta * rho / mu;

  double Tplus;

  if (dplus < 13.2)
    Tplus = prandtl * dplus;
  else
    Tplus = 2.195*log(dplus) + 13.2*prandtl - 5.66;

  return -dT / Tplus;

}

//------------------------------------------------------------------------------

void WallFcn::computeSurfaceTerm(int code, Vec3D &normal, double d2wall[3],
				 double *Vwall, double *V[3], double *term)
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);

  double a = - rhow * utau*utau * norm;

  term[0] = 0.0;
  term[1] = a * t[0];
  term[2] = a * t[1];
  term[3] = a * t[2];
  term[4] = a * (uw * t);

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
    double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);
    term[4] += gam * rhow * utau * Ttau * norm;
  }

  computeWallValues(utau, delta, rhow, du*t, muw, Vwall);

}

//------------------------------------------------------------------------------

Vec3D WallFcn::computeForce(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);
  //double utau = 1.0;

  Vec3D force = - rhow * utau*utau * norm * t;

  return force;

}

//------------------------------------------------------------------------------

double WallFcn::computeHeatPower(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);

  double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);

  double hp = - gam * rhow * utau * Ttau * norm;

  return hp;

}

//------------------------------------------------------------------------------

double WallFcn::computeInterfaceWork(Vec3D& normal, double d2wall[3], double* Vwall, double* V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);

  double W = - rhow * utau*utau * norm * (uw * t);

  return W;

}

//------------------------------------------------------------------------------

double WallFcn::computeDeltaPlus(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);

  return reynolds * utau * delta * rhow / muw;

}

//------------------------------------------------------------------------------

WallFcnSA::WallFcnSA(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  WallFcn(iod, varf, visf)
{

  double cv1 = iod.eqs.tc.tm.sa.cv1;
  cv1_pow3 = cv1*cv1*cv1;

}

//------------------------------------------------------------------------------

void WallFcnSA::computeWallValues(double utau, double delta, double rho, 
				  double ut, double mu, double *V)
{

  const int maxits = 20;
  const double eps = 1.e-6;
  const double coef0 = 0.4 * exp(-0.4*5.5);

  double kuplus = 0.4 * ut / utau;

  /*
  double dplus = reynolds * utau * delta * rho / mu;
  if (dplus > 200.0)
    kuplus = 0.4 * 75.0;
  */

  double mutomu = coef0 * (exp(kuplus) - 1.0 - kuplus - 0.5*kuplus*kuplus);

  if(mutomu<0.0) {
//     fprintf(stderr,"*** Warning: mutomu clipped at wall boundary (mutomu = %e)\n",mutomu);
     mutomu = 0.0;
  }

  double mutildeomu;
  if (mutomu < 3.5)
    mutildeomu = pow(mutomu*cv1_pow3, 0.25);
  else
    mutildeomu = mutomu;

  int it;
  double res, target;
  for (it=0; it<maxits; ++it) {
    double mutildeomu2 = mutildeomu*mutildeomu;
    double mutildeomu3 = mutildeomu2*mutildeomu;
    double f = mutildeomu2*mutildeomu2 - mutomu*(mutildeomu3 + cv1_pow3);
    res = f * f;
    if (it == 0)
      target = eps * res;
    if (res == 0.0 || res <= target)
      break;
    double df = 4.0*mutildeomu3 - 3.0*mutomu*mutildeomu2;
    mutildeomu -= f / df;
  }
  if (it == maxits && target > 1.e-16) {
    fprintf(stderr, "*** Warning: Newton reached %d its on mutilde", maxits);
    fprintf(stderr, " (res=%.2e, target=%.2e)\n", res, target);
  }

  V[5] = mu * mutildeomu;
 
}

//------------------------------------------------------------------------------

WallFcnKE::WallFcnKE(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  WallFcn(iod, varf, visf)
{

  orcmu = 1.0 / sqrt(iod.eqs.tc.tm.ke.c_mu);

}

//------------------------------------------------------------------------------

void WallFcnKE::computeWallValues(double utau, double delta, double rho, 
				  double ut, double mu, double *V)
{

  double nu = mu / rho;
  double dplus = reynolds * utau * delta / nu;
  double dplus2 = dplus * dplus;
  double utau2 = utau * utau;

  double k, eps;

  if (dplus < 10.0) {
    k = orcmu * utau2 * 0.01 * dplus2;
    eps = reynolds * utau2 * utau2 * 0.1 / (vkcst * nu);
    eps *= 0.01 * dplus2 + 0.2*vkcst*orcmu * (1.0 - 0.01 * dplus2);
  }
  else {
    k = orcmu * utau2;
    eps = utau2 * utau / (vkcst * delta);
  }

  V[5] = k;
  V[6] = eps;
 
}

//------------------------------------------------------------------------------
