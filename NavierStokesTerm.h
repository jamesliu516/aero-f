#ifndef _NAVIER_STOKES_TERM_H_
#define _NAVIER_STOKES_TERM_H_

#include <IoData.h>
#include <VarFcn.h>
#include <ViscoFcn.h>
#include <ThermalCondFcn.h>

struct Vec3D;
//-----------------------------------------------------------------------------
//CHANGES_FOR_WATER
// same as in FemEquationTermDesc.C : Lame coefficients behave differently in 
//  gas and water, thus we had to "split" them, and use the Stokes relation 
//   only in the case of gas.
//------------------------------------------------------------------------------

class NavierStokesTerm {

protected:

  static const double third;
  static const double twothird;
  static const double fourth;

  double ooreynolds;
  double ooreynolds_mu;
  double ooreynolds_lambda;

  VarFcn *varFcn;
  ViscoFcn *viscoFcn;
  ThermalCondFcn *thermalCondFcn;

  void computeTemperature(double *, double &);
  void computeVelocity(double *[4], double [4][3], double [3]);
  void computeTemperature(double *[4], double [4], double &);
  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);
  void computeTemperatureGradient(double [4][3], double [4], double [3]);
  void computeStressTensor(double, double, double [3][3], double [3][3]);
  void computeHeatFluxVector(double, double [3], double [3]);

  template<int dim>
  void computeVolumeTermNS(double, double, double, double [3], double [3][3], 
			   double [3], double (*)[dim]);

  template<int neq>
  void computeJacobianVolumeTermNS(double [4][3], double, double,  double, double *[4], 
				   double [4], double (*)[3][neq][neq]);

  void computeSurfaceTermNS(double [4][3], Vec3D &, double *, double *[4], double *);

  template<int neq>
  void computeJacobianSurfaceTermNS(double [4][3], Vec3D &, double *, 
				    double *[4], double (*)[neq][neq]);
  
public:

  NavierStokesTerm(IoData &, VarFcn *);
  ~NavierStokesTerm();

};

//------------------------------------------------------------------------------

inline
NavierStokesTerm::NavierStokesTerm(IoData &iod, VarFcn *vf) : varFcn(vf)
{

  ooreynolds = 1.0/iod.ref.reynolds_mu;
  ooreynolds_mu = 1.0/iod.ref.reynolds_mu;
  ooreynolds_lambda = 1.0/iod.ref.reynolds_lambda;

  viscoFcn = 0;
  thermalCondFcn = 0;

  if (iod.eqs.viscosityModel.type == ViscosityModelData::CONSTANT)
    viscoFcn = new ConstantViscoFcn(iod);
  else if (iod.eqs.viscosityModel.type == ViscosityModelData::SUTHERLAND)
    viscoFcn = new SutherlandViscoFcn(iod);
  else if (iod.eqs.viscosityModel.type == ViscosityModelData::PRANDTL)
    viscoFcn = new PrandtlViscoFcn(iod);
//  else if (iod.eqs.viscosityModel.type == ViscosityModelData::WATER)
//    viscoFcn = new WaterViscoFcn(iod);

  if (iod.eqs.thermalCondModel.type == ThermalCondModelData::CONSTANT_PRANDTL)
    thermalCondFcn = new ConstantPrandtlThermalCondFcn(iod, viscoFcn);
//  else if (iod.eqs.thermalCondModel.type == ThermalCondModelData::WATER)
//    thermalCondFcn = new WaterThermalCondFcn(iod, viscoFcn);

}

//------------------------------------------------------------------------------

inline
NavierStokesTerm::~NavierStokesTerm()
{

  if (viscoFcn) delete viscoFcn;
  if (thermalCondFcn) delete thermalCondFcn;

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeTemperature(double *V, double &T)
{
  T = varFcn->computeTemperature(V);
}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeVelocity(double *V[4], double u[4][3], double ucg[3])
{

  u[0][0] = V[0][1];
  u[0][1] = V[0][2];
  u[0][2] = V[0][3];

  u[1][0] = V[1][1];
  u[1][1] = V[1][2];
  u[1][2] = V[1][3];

  u[2][0] = V[2][1];
  u[2][1] = V[2][2];
  u[2][2] = V[2][3];

  u[3][0] = V[3][1];
  u[3][1] = V[3][2];
  u[3][2] = V[3][3];

  ucg[0] = fourth * (u[0][0] + u[1][0] + u[2][0] + u[3][0]);
  ucg[1] = fourth * (u[0][1] + u[1][1] + u[2][1] + u[3][1]);
  ucg[2] = fourth * (u[0][2] + u[1][2] + u[2][2] + u[3][2]);

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeTemperature(double *V[4], double T[4], double &Tcg)
{

  T[0] = varFcn->computeTemperature(V[0]);
  T[1] = varFcn->computeTemperature(V[1]);
  T[2] = varFcn->computeTemperature(V[2]);
  T[3] = varFcn->computeTemperature(V[3]);

  Tcg = fourth * (T[0] + T[1] + T[2] + T[3]);

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeVelocityGradient(double dp1dxj[4][3], double u[4][3],
					       double dudxj[3][3])
{

  dudxj[0][0] = dp1dxj[0][0]*u[0][0] + dp1dxj[1][0]*u[1][0] + 
    dp1dxj[2][0]*u[2][0] + dp1dxj[3][0]*u[3][0];

  dudxj[0][1] = dp1dxj[0][1]*u[0][0] + dp1dxj[1][1]*u[1][0] + 
    dp1dxj[2][1]*u[2][0] + dp1dxj[3][1]*u[3][0];

  dudxj[0][2] = dp1dxj[0][2]*u[0][0] + dp1dxj[1][2]*u[1][0] + 
    dp1dxj[2][2]*u[2][0] + dp1dxj[3][2]*u[3][0];

  dudxj[1][0] = dp1dxj[0][0]*u[0][1] + dp1dxj[1][0]*u[1][1] + 
    dp1dxj[2][0]*u[2][1] + dp1dxj[3][0]*u[3][1];

  dudxj[1][1] = dp1dxj[0][1]*u[0][1] + dp1dxj[1][1]*u[1][1] + 
    dp1dxj[2][1]*u[2][1] + dp1dxj[3][1]*u[3][1];

  dudxj[1][2] = dp1dxj[0][2]*u[0][1] + dp1dxj[1][2]*u[1][1] + 
    dp1dxj[2][2]*u[2][1] + dp1dxj[3][2]*u[3][1];

  dudxj[2][0] = dp1dxj[0][0]*u[0][2] + dp1dxj[1][0]*u[1][2] + 
    dp1dxj[2][0]*u[2][2] + dp1dxj[3][0]*u[3][2];

  dudxj[2][1] = dp1dxj[0][1]*u[0][2] + dp1dxj[1][1]*u[1][2] + 
    dp1dxj[2][1]*u[2][2] + dp1dxj[3][1]*u[3][2];

  dudxj[2][2] = dp1dxj[0][2]*u[0][2] + dp1dxj[1][2]*u[1][2] + 
    dp1dxj[2][2]*u[2][2] + dp1dxj[3][2]*u[3][2];

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeTemperatureGradient(double dp1dxj[4][3], double T[4],
						  double dTdxj[3])
{

  dTdxj[0] = dp1dxj[0][0]*T[0] + dp1dxj[1][0]*T[1] + 
    dp1dxj[2][0]*T[2] + dp1dxj[3][0]*T[3];

  dTdxj[1] = dp1dxj[0][1]*T[0] + dp1dxj[1][1]*T[1] + 
    dp1dxj[2][1]*T[2] + dp1dxj[3][1]*T[3];

  dTdxj[2] = dp1dxj[0][2]*T[0] + dp1dxj[1][2]*T[1] + 
    dp1dxj[2][2]*T[2] + dp1dxj[3][2]*T[3];

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeStressTensor(double mu, double lambda, double dudxj[3][3], double tij[3][3])
{

  double dij[3][3];

  dij[0][0] = mu * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
  dij[1][1] = mu * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
  dij[2][2] = mu * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
  dij[0][1] = mu * (dudxj[1][0] + dudxj[0][1]);
  dij[0][2] = mu * (dudxj[2][0] + dudxj[0][2]);
  dij[1][2] = mu * (dudxj[2][1] + dudxj[1][2]);
  dij[1][0] = dij[0][1];
  dij[2][0] = dij[0][2];
  dij[2][1] = dij[1][2];


  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2]; 

  tij[0][0] = lambda * div + 2.0 * mu *dudxj[0][0];
  tij[1][1] = lambda * div + 2.0 * mu *dudxj[1][1];
  tij[2][2] = lambda * div + 2.0 * mu *dudxj[2][2];
  tij[0][1] = mu * (dudxj[1][0] + dudxj[0][1]);
  tij[0][2] = mu * (dudxj[2][0] + dudxj[0][2]);
  tij[1][2] = mu * (dudxj[2][1] + dudxj[1][2]);
  tij[1][0] = tij[0][1];
  tij[2][0] = tij[0][2];
  tij[2][1] = tij[1][2];

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeHeatFluxVector(double kappa, double dTdxj[3], double qj[3])
{

  qj[0] = - kappa * dTdxj[0];
  qj[1] = - kappa * dTdxj[1];
  qj[2] = - kappa * dTdxj[2];

}

//------------------------------------------------------------------------------

template<int dim>
inline
void NavierStokesTerm::computeVolumeTermNS(double mu, double lambda, double kappa, double u[3], 
					   double dudxj[3][3], double dTdxj[3], 
					   double (*r)[dim])
{

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double qj[3];
  computeHeatFluxVector(kappa, dTdxj, qj);

  r[0][0] = 0.0;
  r[0][1] = tij[0][0];
  r[0][2] = tij[1][0];
  r[0][3] = tij[2][0];
  r[0][4] = u[0] * tij[0][0] + u[1] * tij[1][0] + u[2] * tij[2][0] - qj[0];

  r[1][0] = 0.0;
  r[1][1] = tij[0][1];
  r[1][2] = tij[1][1];
  r[1][3] = tij[2][1];
  r[1][4] = u[0] * tij[0][1] + u[1] * tij[1][1] + u[2] * tij[2][1] - qj[1]; 

  r[2][0] = 0.0;
  r[2][1] = tij[0][2];
  r[2][2] = tij[1][2];
  r[2][3] = tij[2][2];
  r[2][4] = u[0] * tij[0][2] + u[1] * tij[1][2] + u[2] * tij[2][2] - qj[2];

}

//------------------------------------------------------------------------------

template<int neq>
inline
void NavierStokesTerm::
computeJacobianVolumeTermNS(double dp1dxj[4][3], double mu, double lambda, double kappa, 
			    double *V[4], double T[4], double (*dRdU)[3][neq][neq])
{

  double uu[4][3], ucg[3];
  computeVelocity(V, uu, ucg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, uu, dudxj);

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double txx = tij[0][0];
  double tyy = tij[1][1];
  double tzz = tij[2][2];
  double txy = tij[0][1];
  double txz = tij[0][2];
  double tyz = tij[1][2];

  for (int k=0; k<4; ++k) {
    
    double rho = V[k][0];
    double u = V[k][1];
    double v = V[k][2];
    double w = V[k][3];
    double Temp = T[k];

    double dp1dx = dp1dxj[k][0];
    double dp1dy = dp1dxj[k][1];
    double dp1dz = dp1dxj[k][2];

    double nu = mu / rho;
    double lu = lambda / rho;

//---------------------
    double adtxxdu0 = twothird * nu * (-2.0 * dp1dx * u + dp1dy * v + dp1dz * w);
    double adtxxdu1 = 2.0 * twothird * nu * dp1dx;
    double adtxxdu2 = - twothird * nu * dp1dy;
    double adtxxdu3 = - twothird * nu * dp1dz;

    double adtyydu0 = twothird * nu * (-2.0 * dp1dy * v + dp1dx * u + dp1dz * w);
    double adtyydu1 = - twothird * nu * dp1dx;
    double adtyydu2 = 2.0 * twothird * nu * dp1dy;
    double adtyydu3 = - twothird * nu * dp1dz;

    double adtzzdu0 = twothird * nu * (-2.0 * dp1dz * w + dp1dx * u + dp1dy * v);
    double adtzzdu1 = - twothird * nu * dp1dx;
    double adtzzdu2 = - twothird * nu * dp1dy;
    double adtzzdu3 = 2.0 * twothird * nu * dp1dz;

    double adtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double adtxydu1 = nu * dp1dy;
    double adtxydu2 = nu * dp1dx;

    double adtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double adtxzdu1 = nu * dp1dz;
    double adtxzdu3 = nu * dp1dx;

    double adtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double adtyzdu2 = nu * dp1dz;
    double adtyzdu3 = nu * dp1dy;

//--------------------

    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;

    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;

    double alpha = kappa / rho;
    double c0 = alpha * (Temp - 0.5 * (u*u + v*v + w*w));
    double c1 = alpha * u;
    double c2 = alpha * v;
    double c3 = alpha * w;
    double c4 = - alpha;

    double dqxdu0 = c0 * dp1dx;
    double dqxdu1 = c1 * dp1dx;
    double dqxdu2 = c2 * dp1dx;
    double dqxdu3 = c3 * dp1dx;
    double dqxdu4 = c4 * dp1dx;

    double dqydu0 = c0 * dp1dy;
    double dqydu1 = c1 * dp1dy;
    double dqydu2 = c2 * dp1dy;
    double dqydu3 = c3 * dp1dy;
    double dqydu4 = c4 * dp1dy;

    double dqzdu0 = c0 * dp1dz;
    double dqzdu1 = c1 * dp1dz;
    double dqzdu2 = c2 * dp1dz;
    double dqzdu3 = c3 * dp1dz;
    double dqzdu4 = c4 * dp1dz;

    double beta = 0.25 / rho;

    double dudu0 = - beta * u;
    double dudu1 = beta;

    double dvdu0 = - beta * v;
    double dvdu2 = beta;

    double dwdu0 = - beta * w;
    double dwdu3 = beta;


    //computation of dRdU -> dRdU[k][a][b][c]
// k is the corner of the tetrahedra we are considering
// a is the coordinate direction: x y z
// b is the equation we are considering: 0 for mass conservation
//					 1 for x-momentum
//					 2 for y-momentum
//					 3 for z-momentum
//					 4 for energy
// c is the variable wrt which we derive

    // dRxdU

    dRdU[k][0][0][0] = 0.0;
    dRdU[k][0][0][1] = 0.0;
    dRdU[k][0][0][2] = 0.0;
    dRdU[k][0][0][3] = 0.0;
    dRdU[k][0][0][4] = 0.0;

    dRdU[k][0][1][0] = dtxxdu0;
    dRdU[k][0][1][1] = dtxxdu1;
    dRdU[k][0][1][2] = dtxxdu2;
    dRdU[k][0][1][3] = dtxxdu3;
    dRdU[k][0][1][4] = 0.0;

    dRdU[k][0][2][0] = dtxydu0;
    dRdU[k][0][2][1] = dtxydu1;
    dRdU[k][0][2][2] = dtxydu2;
    dRdU[k][0][2][3] = 0.0;
    dRdU[k][0][2][4] = 0.0;

    dRdU[k][0][3][0] = dtxzdu0;
    dRdU[k][0][3][1] = dtxzdu1;
    dRdU[k][0][3][2] = 0.0;
    dRdU[k][0][3][3] = dtxzdu3;
    dRdU[k][0][3][4] = 0.0;

    dRdU[k][0][4][0] = dudu0 * txx + ucg[0] * dtxxdu0 + dvdu0 * txy + ucg[1] * dtxydu0 + 
                       dwdu0 * txz + ucg[2] * dtxzdu0 - dqxdu0;
    dRdU[k][0][4][1] = dudu1 * txx + ucg[0] * dtxxdu1 +
                       ucg[1] * dtxydu1 + ucg[2] * dtxzdu1 - dqxdu1;
    dRdU[k][0][4][2] = ucg[0] * dtxxdu2 + dvdu2 * txy + ucg[1] * dtxydu2 - dqxdu2;
    dRdU[k][0][4][3] = ucg[0] * dtxxdu3 + dwdu3 * txz + ucg[2] * dtxzdu3 - dqxdu3;
    dRdU[k][0][4][4] = - dqxdu4;

    // dRydU

    dRdU[k][1][0][0] = 0.0;
    dRdU[k][1][0][1] = 0.0;
    dRdU[k][1][0][2] = 0.0;
    dRdU[k][1][0][3] = 0.0;
    dRdU[k][1][0][4] = 0.0;
  
    dRdU[k][1][1][0] = dtxydu0;
    dRdU[k][1][1][1] = dtxydu1;
    dRdU[k][1][1][2] = dtxydu2;
    dRdU[k][1][1][3] = 0.0;
    dRdU[k][1][1][4] = 0.0;

    dRdU[k][1][2][0] = dtyydu0;
    dRdU[k][1][2][1] = dtyydu1;
    dRdU[k][1][2][2] = dtyydu2;
    dRdU[k][1][2][3] = dtyydu3;
    dRdU[k][1][2][4] = 0.0;

    dRdU[k][1][3][0] = dtyzdu0;
    dRdU[k][1][3][1] = 0.0;
    dRdU[k][1][3][2] = dtyzdu2;
    dRdU[k][1][3][3] = dtyzdu3;
    dRdU[k][1][3][4] = 0.0;

    dRdU[k][1][4][0] = dudu0 * txy + ucg[0] * dtxydu0 + dvdu0 * tyy + ucg[1] * dtyydu0 + 
                       dwdu0 * tyz + ucg[2] * dtyzdu0 - dqydu0;
    dRdU[k][1][4][1] = dudu1 * txy + ucg[0] * dtxydu1 + ucg[1] * dtyydu1 - dqydu1;
    dRdU[k][1][4][2] = ucg[0] * dtxydu2 + dvdu2 * tyy + 
                       ucg[1] * dtyydu2 + ucg[2] * dtyzdu2 - dqydu2;
    dRdU[k][1][4][3] = ucg[1] * dtyydu3 + dwdu3 * tyz + ucg[2] * dtyzdu3 - dqydu3;
    dRdU[k][1][4][4] = - dqydu4;

    // dRzdU

    dRdU[k][2][0][0] = 0.0;
    dRdU[k][2][0][1] = 0.0;
    dRdU[k][2][0][2] = 0.0;
    dRdU[k][2][0][3] = 0.0;
    dRdU[k][2][0][4] = 0.0;

    dRdU[k][2][1][0] = dtxzdu0;
    dRdU[k][2][1][1] = dtxzdu1;
    dRdU[k][2][1][2] = 0.0;
    dRdU[k][2][1][3] = dtxzdu3;
    dRdU[k][2][1][4] = 0.0;

    dRdU[k][2][2][0] = dtyzdu0;
    dRdU[k][2][2][1] = 0.0;
    dRdU[k][2][2][2] = dtyzdu2;
    dRdU[k][2][2][3] = dtyzdu3;
    dRdU[k][2][2][4] = 0.0;

    dRdU[k][2][3][0] = dtzzdu0;
    dRdU[k][2][3][1] = dtzzdu1;
    dRdU[k][2][3][2] = dtzzdu2;
    dRdU[k][2][3][3] = dtzzdu3;
    dRdU[k][2][3][4] = 0.0;

    dRdU[k][2][4][0] = dudu0 * txz + ucg[0] * dtxzdu0 + dvdu0 * tyz + ucg[1] * dtyzdu0 + 
                       dwdu0 * tzz + ucg[2] * dtzzdu0 - dqzdu0;
    dRdU[k][2][4][1] = dudu1 * txz + ucg[0] * dtxzdu1 + ucg[2] * dtzzdu1 - dqzdu1;
    dRdU[k][2][4][2] = dvdu2 * tyz + ucg[1] * dtyzdu2 + ucg[2] * dtzzdu2 - dqzdu2;
    dRdU[k][2][4][3] = ucg[0] * dtxzdu3 + ucg[1] * dtyzdu3 + 
                       dwdu3 * tzz + ucg[2] * dtzzdu3 - dqzdu3;
    dRdU[k][2][4][4] = - dqzdu4;

  }

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeSurfaceTermNS(double dp1dxj[4][3], Vec3D &n, 
					    double *Vwall, double *Vtet[4], double *R)
{

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

  R[0] = 0.0;
  R[1] = 0.0;
  R[2] = 0.0;
  R[3] = 0.0;
  R[4] = (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * n[0] +
    (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * n[1] + 
    (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * n[2];

}

//------------------------------------------------------------------------------

template<int neq>
inline
void NavierStokesTerm::computeJacobianSurfaceTermNS(double dp1dxj[4][3], Vec3D &n, 
						    double *Vwall, double *V[4], 
						    double (*dRdU)[neq][neq])
{

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

  for (int k=0; k<4; ++k) {
    
    double rho = V[k][0];
    double u = V[k][1];
    double v = V[k][2];
    double w = V[k][3];

    double dp1dx = dp1dxj[k][0];
    double dp1dy = dp1dxj[k][1];
    double dp1dz = dp1dxj[k][2];

    double nu = mu / rho;
    double lu = lambda / rho;

    double adtxxdu0 = twothird * nu * (-2.0 * dp1dx * u + dp1dy * v + dp1dz * w);
    double adtxxdu1 = 2.0 * twothird * nu * dp1dx;
    double adtxxdu2 = - twothird * nu * dp1dy;
    double adtxxdu3 = - twothird * nu * dp1dz;

    double adtyydu0 = twothird * nu * (-2.0 * dp1dy * v + dp1dx * u + dp1dz * w);
    double adtyydu1 = - twothird * nu * dp1dx;
    double adtyydu2 = 2.0 * twothird * nu * dp1dy;
    double adtyydu3 = - twothird * nu * dp1dz;

    double adtzzdu0 = twothird * nu * (-2.0 * dp1dz * w + dp1dx * u + dp1dy * v);
    double adtzzdu1 = - twothird * nu * dp1dx;
    double adtzzdu2 = - twothird * nu * dp1dy;
    double adtzzdu3 = 2.0 * twothird * nu * dp1dz;



    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;
    
//---check
/*
    if (adtxxdu0!=dtxxdu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu0 = %e and adtxxdu0 = %f differ\n", dtxxdu0, adtxxdu0);
    if (adtxxdu1!=dtxxdu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu1 = %e and adtxxdu1 = %f differ\n", dtxxdu1, adtxxdu1);
    if (adtxxdu2!=dtxxdu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu2 = %e and adtxxdu2 = %f differ\n", dtxxdu2, adtxxdu2);
    if (adtxxdu3!=dtxxdu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu3 = %e and adtxxdu3 = %f differ\n", dtxxdu3, adtxxdu3);
                                                                                                                                                                                                     
                                                                                                                                                                                                     
                                                                                                                                                                                                     
                                                                                                                                                                                                     
    if (adtyydu0!=dtyydu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu0 = %e and adtyydu0 = %f differ\n", dtyydu0, adtyydu0);
    if (adtyydu1!=dtyydu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu1 = %e and adtyydu1 = %f differ\n", dtyydu1, adtyydu1);
    if (adtyydu2!=dtyydu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu2 = %e and adtyydu2 = %f differ\n", dtyydu2, adtyydu2);
    if (adtyydu3!=dtyydu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu3 = %e and adtyydu3 = %f differ\n", dtyydu3, adtyydu3);
                                                                                                                                                                                                     
    if (adtzzdu0!=dtzzdu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu0 = %e and adtzzdu0 = %f differ\n", dtzzdu0, adtzzdu0);
    if (adtzzdu1!=dtzzdu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu1 = %e and adtzzdu1 = %f differ\n", dtzzdu1, adtzzdu1);
    if (adtzzdu2!=dtzzdu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu2 = %e and adtzzdu2 = %f differ\n", dtzzdu2, adtzzdu2);
    if (adtzzdu3!=dtzzdu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu3 = %e and adtzzdu3 = %f differ\n", dtzzdu3, adtzzdu3);
  
*/
    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = 0.0;
    dRdU[k][1][1] = 0.0;
    dRdU[k][1][2] = 0.0;
    dRdU[k][1][3] = 0.0;
    dRdU[k][1][4] = 0.0;

    dRdU[k][2][0] = 0.0;
    dRdU[k][2][1] = 0.0;
    dRdU[k][2][2] = 0.0;
    dRdU[k][2][3] = 0.0;
    dRdU[k][2][4] = 0.0;

    dRdU[k][3][0] = 0.0;
    dRdU[k][3][1] = 0.0;
    dRdU[k][3][2] = 0.0;
    dRdU[k][3][3] = 0.0;
    dRdU[k][3][4] = 0.0;

    dRdU[k][4][0] = (Vwall[1]*dtxxdu0 + Vwall[2]*dtxydu0 + Vwall[3]*dtxzdu0) * n[0] +
      (Vwall[1]*dtxydu0 + Vwall[2]*dtyydu0 + Vwall[3]*dtyzdu0) * n[1] + 
      (Vwall[1]*dtxzdu0 + Vwall[2]*dtyzdu0 + Vwall[3]*dtzzdu0) * n[2];

    dRdU[k][4][1] = (Vwall[1]*dtxxdu1 + Vwall[2]*dtxydu1 + Vwall[3]*dtxzdu1) * n[0] +
      (Vwall[1]*dtxydu1 + Vwall[2]*dtyydu1) * n[1] + 
      (Vwall[1]*dtxzdu1 + Vwall[3]*dtzzdu1) * n[2];

    dRdU[k][4][2] = (Vwall[1]*dtxxdu2 + Vwall[2]*dtxydu2) * n[0] +
      (Vwall[1]*dtxydu2 + Vwall[2]*dtyydu2 + Vwall[3]*dtyzdu2) * n[1] + 
      (Vwall[2]*dtyzdu2 + Vwall[3]*dtzzdu2) * n[2];

    dRdU[k][4][3] = (Vwall[1]*dtxxdu3 + Vwall[3]*dtxzdu3) * n[0] +
      (Vwall[2]*dtyydu3 + Vwall[3]*dtyzdu3) * n[1] + 
      (Vwall[1]*dtxzdu3 + Vwall[2]*dtyzdu3 + Vwall[3]*dtzzdu3) * n[2];

    dRdU[k][4][4] = 0.0;

  }

}

//------------------------------------------------------------------------------

#endif
