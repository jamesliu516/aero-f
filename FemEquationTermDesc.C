#include <FemEquationTermDesc.h>

#include <stdlib.h>
#include <stdio.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
using std::min;
#endif
//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
// as the Stokes' law for gas does not apply anymore, we have to distinguish
// between lambda and mu (the Lame coefficients)
//------------------------------------------------------------------------------

const double NavierStokesTerm::third = 1.0/3.0;
const double NavierStokesTerm::twothird = 2.0/3.0;
const double NavierStokesTerm::fourth = 1.0/4.0;

//------------------------------------------------------------------------------

FemEquationTermNS::FemEquationTermNS(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length; 

}

//------------------------------------------------------------------------------

double FemEquationTermNS::computeViscousTimeStep(double X[3], double *V)
{
  double T;
  computeTemperature(V,T);
  double mul = ooreynolds_mu * viscoFcn->compute_mu(T);
  return mul/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermNS::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol, 
                                          SVec<double,3> &X, int nodeNum[4],  
                                          int volume_id)
{

  bool porousmedia = false; 

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);


  double (*R)[5] = reinterpret_cast<double (*)[5]>(r);
  computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
       porousmedia = true;

       double RR[9], K[9];
       double alpha[3], beta[3];
       double volten = tetVol *(1.0/10.0);

       // non-dimensionalization //

       alpha[0] = it->second->alphax*length/density; 
       alpha[1] = it->second->alphay*length/density;  
       alpha[2] = it->second->alphaz*length/density;

       beta[0]  = it->second->betax*length/(density*velocity);  
       beta[1]  = it->second->betay*length/(density*velocity);   
       beta[2]  = it->second->betaz*length/(density*velocity);

       // transformation matrix // 

       RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
       RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
       RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

       // permittivity matrix  //

       computePermittivityTensor(alpha, beta, ucg, RR, K);

       double SS[4][3];

       SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
       SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
       SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                         
       SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
       SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
       SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                         
       SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
       SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
       SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                         
       SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
       SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
       SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                         
       // FE flux for the porous sink term //
       
       for (int j=0; j<12; ++j) PR[j] = 0.0; // initialize PR
                                                                                                                         
       for (int j=0; j<4; ++j) {
         for (int k=0; k<3; ++k)
           PR[3*j+k] += (K[3*k+0] * SS[j][0] + K[3*k+1] * SS[j][1] + K[3*k+2] * SS[j][2]);
       }
    } 
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      for (int j=0; j<12; ++j) PR[j] = 0.0; 
    }
  }
 
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

  return (porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermNS::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
       porousmedia = true;
       double RR[9], K[9], B[9];
       double alpha[3], beta[3];

       double volten = tetVol *(1.0/10.0);
       double onefourth = 1.0/4.0;
       double onehalf = 1.0/2.0;
  
       double u[4][3], ucg[3];
       computeVelocity(V, u, ucg);

       // non-dimensionalization correction

       alpha[0] = it->second->alphax*length/density;
       alpha[1] = it->second->alphay*length/density;
       alpha[2] = it->second->alphaz*length/density;

       beta[0]  = it->second->betax*length/(density*velocity);
       beta[1]  = it->second->betay*length/(density*velocity);
       beta[2]  = it->second->betaz*length/(density*velocity);

       // transformation matrix
       RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
       RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
       RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

       // permittivity matrix
       computePermittivityTensor(alpha, beta, ucg, RR, K);
      
       // gradient of permittivity matrix
       computeGradPermittivityTensor(alpha, ucg, RR, B);

       double SS[4][3];

       SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
       SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
       SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
       SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
       SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
       SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
       SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
       SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
       SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
       SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
       SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
       SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
       for (int k=0; k<4; ++k) {
         double BB[3];
         BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
         BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
         BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

         double BV[9];
         BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
         BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
         BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
         for (int j=0; j<4; ++j) {
           double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
           double KU[25];
           multiplyBydVdU(v, K, KU, volten);
           double BU[25];
           multiplyBydVdU(v, BV, BU, 1.0);

           if (k == j) {
              dPdU[k][j][0][0] = 0.0; 
              dPdU[k][j][0][1] = 0.0;
              dPdU[k][j][0][2] = 0.0;
              dPdU[k][j][0][3] = 0.0;
              dPdU[k][j][0][4] = 0.0;

              dPdU[k][j][1][0] = KU[5] + BU[5];
              dPdU[k][j][1][1] = KU[6] + BU[6];
              dPdU[k][j][1][2] = KU[7] + BU[7];
              dPdU[k][j][1][3] = KU[8] + BU[8];
              dPdU[k][j][1][4] = KU[9] + BU[9];

              dPdU[k][j][2][0] = KU[10] + BU[10];
              dPdU[k][j][2][1] = KU[11] + BU[11];
              dPdU[k][j][2][2] = KU[12] + BU[12];
              dPdU[k][j][2][3] = KU[13] + BU[13];
              dPdU[k][j][2][4] = KU[14] + BU[14];
          
              dPdU[k][j][3][0] = KU[15] + BU[15];
              dPdU[k][j][3][1] = KU[16] + BU[16];
              dPdU[k][j][3][2] = KU[17] + BU[17];
              dPdU[k][j][3][3] = KU[18] + BU[18];
              dPdU[k][j][3][4] = KU[19] + BU[19];

              dPdU[k][j][4][0] = 0.0; 
              dPdU[k][j][4][1] = 0.0;
              dPdU[k][j][4][2] = 0.0;
              dPdU[k][j][4][3] = 0.0;
              dPdU[k][j][4][4] = 0.0;
          }
          else {
              dPdU[k][j][0][0] = 0.0; 
              dPdU[k][j][0][1] = 0.0;
              dPdU[k][j][0][2] = 0.0;
              dPdU[k][j][0][3] = 0.0;
              dPdU[k][j][0][4] = 0.0;

              dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
              dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
              dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
              dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
              dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];

              dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
              dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
              dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
              dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
              dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
          
              dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
              dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
              dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
              dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
              dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];

              dPdU[k][j][4][0] = 0.0; 
              dPdU[k][j][4][1] = 0.0;
              dPdU[k][j][4][2] = 0.0;
              dPdU[k][j][4][3] = 0.0;
              dPdU[k][j][4][4] = 0.0;
          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      for(int i=0; i<5; ++i)
         for(int j=0; j<5; ++j)
           for(int k=0; k<5; ++k)
             for(int l=0; l<5; ++l)
                dPdU[i][j][k][l] = 0.0;
    }
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{
  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*5*5; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSA::FemEquationTermSA(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcnSA(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;
  
  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.coupling == ImplicitData::STRONG && trip==1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-SA simulation \n");
    trip = 0;
  }


  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermSA::computeViscousTimeStep(double X[3], double *V)
{
  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermSA::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int volume_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip
  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, mul, mutilde);}
    else{
       computeTurbulentViscosity(V, mul, mutilde);
       mut = 0.0; 
    }
  }
  else{
    mut = computeTurbulentViscosity(V, mul, mutilde);
  }


  double mu;
  double kappa;
  double (*R)[6] = reinterpret_cast<double (*)[6]>(r);

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double mu5 = oosigma * (mul + absmutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];

  R[0][5] = mu5 * dnutildedx;
  R[1][5] = mu5 * dnutildedy;
  R[2][5] = mu5 * dnutildedz;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);

  if  (d2wall >= 1.e-15) {
    double chi = max(mutilde/mul, 0.001);
    double chi3 = chi*chi*chi;
    double fv1 = chi3 / (chi3 + cv1_pow3);
    double fv2 = 1.0 + oocv2*chi;
    fv2 = 1.0 / (fv2*fv2*fv2);
    double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
    double oorho = 1.0 / rho;
    double zz = ooreynolds_mu * oovkcst2 * mutilde * oorho * ood2wall2;
    double s12 = dudxj[0][1] - dudxj[1][0];
    double s23 = dudxj[1][2] - dudxj[2][1];
    double s31 = dudxj[2][0] - dudxj[0][2];
    double s = sqrt(s12*s12 + s23*s23 + s31*s31);
    double Stilde = s*fv3 + zz*fv2;
    double rr = min(zz/Stilde, 2.0);
    double rr2 = rr*rr;
    double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
    double gg2 = gg*gg;
    double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);

    double AA = oosigma * cb2 * rho * 
      (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double BB = cb1 * Stilde * absmutilde;
    double CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    S[5] = AA + BB + CC;
  }
  else {
    S[5] = 0.0;
  }

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

      double RR[9], K[9];
      double alpha[3], beta[3];
      double volten = tetVol *(1.0/10.0);
                                                                                                                                       
      // non-dimensionalization correction //
                                                                                                                                       
      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;
                                                                                                                                      
      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);
                                                                                                                                       
      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;
                                                                                                                                       
      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
                                                                                                                                       
      double SS[4][3];
                                                                                                                                       
      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                      
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                       
      // FE flux for the porous sink term //
                                                                                                                                       
      for (int j=0; j<12; ++j) PR[j] = 0.0; // initialize PR
                                                                                                                                       
      for (int j=0; j<4; ++j) {
        for (int k=0; k<3; ++k)
          PR[3*j+k] += (K[3*k+0] * SS[j][0] + K[3*k+1] * SS[j][1] + K[3*k+2] * SS[j][2]);
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid 
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut); 
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R); 
     for (int j=0; j<12; ++j) PR[j] = 0.0;                
   }
  }
  else {
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);
  }
                                                                                                                                       
  return (porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermSA::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  double lambda;
  
  mut = computeTurbulentViscosity(V, mul, mutilde);

  double mu;
  double kappa;

  int k;
  for (k=0; k<4*3*6*6; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*6*6; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][6][6] = reinterpret_cast<double (*)[3][6][6]>(drdu);
  double (*dSdU)[6][6] = reinterpret_cast<double (*)[6][6]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  computeJacobianVolumeTermSA<6,5>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];
             dPdU[k][j][3][5] =  0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];
             dPdU[k][j][3][5] = 0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{

  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  for (int k=0; k<4*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSAmean::FemEquationTermSAmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermSAmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermSAmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);
                                                                                                                                                                   
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
                                                                                                                                                                   
  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                                                                                   
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip
  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, mul, mutilde);}
    else{
       mut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, mul, mutilde);
  }

  double mu, kappa;

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;

          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;

          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						       Vec3D &n, double d2w[4],
						       double *Vwall, double *V[4],
						       double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSAturb::FemEquationTermSAturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermSAturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde; 
  
  computeTurbulentViscosity(V, mul, mutilde);
  
  double (*dRdU)[3][1][1] = reinterpret_cast<double (*)[3][1][1]>(drdu);
  double (*dSdU)[1][1] = reinterpret_cast<double (*)[1][1]>(dsdu);
  computeJacobianVolumeTermSA<1,0>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU);

  return false;

}

//------------------------------------------------------------------------------

FemEquationTermDES::FemEquationTermDES(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcnSA(iod, varFcn, viscoFcn);

  cdes = iod.eqs.tc.tm.des.cdes;

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.coupling == ImplicitData::STRONG && trip == 1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-DES simulation \n");
    trip = 0;
  }

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermDES::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int volume_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);
                                                                                                   
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip
  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, mul, mutilde);}
    else{
       computeTurbulentViscosity(V, mul, mutilde);
       mut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, mul, mutilde);
  }

  double mu;
  double kappa;

  double (*R)[6] = reinterpret_cast<double (*)[6]>(r);

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double mu5 = oosigma * (mul + absmutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];

  R[0][5] = mu5 * dnutildedx;
  R[1][5] = mu5 * dnutildedy;
  R[2][5] = mu5 * dnutildedz;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

  double maxl,sidel;
  maxl=-1.0;

  for (int i=0; i<4; i++) {
    for (int j=i+1; j<4; j++){
      sidel=sqrt((X[nodeNum[i]][0]-X[nodeNum[j]][0])*(X[nodeNum[i]][0]-X[nodeNum[j]][0]) +
	     (X[nodeNum[i]][1]-X[nodeNum[j]][1])*(X[nodeNum[i]][1]-X[nodeNum[j]][1]) +
	      (X[nodeNum[i]][2]-X[nodeNum[j]][2])*(X[nodeNum[i]][2]-X[nodeNum[j]][2]));
      maxl = max(maxl,sidel);
    }
  }

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  d2wall = min(d2wall,cdes*maxl);

  if  (d2wall >= 1.e-15) {
    double chi = max(mutilde/mul, 0.001);
    double chi3 = chi*chi*chi;
    double fv1 = chi3 / (chi3 + cv1_pow3);
    double fv2 = 1.0 + oocv2*chi;
    fv2 = 1.0 / (fv2*fv2*fv2);
    double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
    double oorho = 1.0 / rho;
    double zz = ooreynolds_mu * oovkcst2 * mutilde * oorho * ood2wall2;
    double s12 = dudxj[0][1] - dudxj[1][0];
    double s23 = dudxj[1][2] - dudxj[2][1];
    double s31 = dudxj[2][0] - dudxj[0][2];
    double s = sqrt(s12*s12 + s23*s23 + s31*s31);
    double Stilde = s*fv3 + zz*fv2;
    double rr = min(zz/Stilde, 2.0);
    double rr2 = rr*rr;
    double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
    double gg2 = gg*gg;
    double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);

    double AA = oosigma * cb2 * rho * 
      (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double BB = cb1 * Stilde * absmutilde;
    double CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    S[5] = AA + BB + CC;
  }
  else {
    S[5] = 0.0;
  }

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

      double RR[9], K[9];
      double alpha[3], beta[3];
      double volten = tetVol *(1.0/10.0);
                                                                                                                                       
      // non-dimensionalization correction //
                                                                                                                                       
      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;
                                                                                                                                      
      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);
                                                                                                                                       
      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;
                                                                                                                                       
      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
                                                                                                                                       
      double SS[4][3];
                                                                                                                                       
      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                      
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                       
      // FE flux for the porous sink term //
                                                                                                                                       
      for (int j=0; j<12; ++j) PR[j] = 0.0; // initialize PR
                                                                                                                                       
      for (int j=0; j<4; ++j) {
        for (int k=0; k<3; ++k)
          PR[3*j+k] += (K[3*k+0] * SS[j][0] + K[3*k+1] * SS[j][1] + K[3*k+2] * SS[j][2]);
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid 
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut); 
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R); 
     for (int j=0; j<12; ++j) PR[j] = 0.0;                
   }
  }
  else {
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);
  }
                                                                                                                                       
  return (porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                   
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  
  mut = computeTurbulentViscosity(V, mul, mutilde);
                                                                                                   
  double mu = ooreynolds_mu * (mul + mut);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  double kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
                                                                                                   
  int k;
  for (k=0; k<4*3*6*6; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*6*6; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][6][6] = reinterpret_cast<double (*)[3][6][6]>(drdu);
  double (*dSdU)[6][6] = reinterpret_cast<double (*)[6][6]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  computeJacobianVolumeTermDES<6,5>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU, X, nodeNum);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];
             dPdU[k][j][3][5] =  0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];
             dPdU[k][j][3][5] = 0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{

  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  for (int k=0; k<4*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermDESmean::FemEquationTermDESmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermDESmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermDESmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);
                                                                                                                                                                   
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
                                                                                                                                                                   
  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                                                                                   
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip

  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, mul, mutilde);}
    else{
       mut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, mul, mutilde);
  }

  double mu, kappa;

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;

          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						       Vec3D &n, double d2w[4],
						       double *Vwall, double *V[4],
						       double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermDESturb::FemEquationTermDESturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermDESturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                  
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  
  computeTurbulentViscosity(V, mul, mutilde);
  
  double (*dRdU)[3][1][1] = reinterpret_cast<double (*)[3][1][1]>(drdu);
  double (*dSdU)[1][1] = reinterpret_cast<double (*)[1][1]>(dsdu);
  computeJacobianVolumeTermDES<1,0>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU, X, nodeNum);

  return false;

}

//------------------------------------------------------------------------------

FemEquationTermKE::FemEquationTermKE(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  wallFcn = new WallFcnKE(iod, varFcn, viscoFcn);

  x0 = iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.coupling == ImplicitData::STRONG && trip == 1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-KEpsilon simulation \n");
    trip = 0;
  }

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermKE::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X, int nodeNum[4], int volume_id)
{
  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double rhok, rhoeps;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip

  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, rhok, rhoeps);}
    else{
       computeTurbulentViscosity(V, rhok, rhoeps);
       mut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, rhok, rhoeps);
  }

  double mu;
  double kappa;

  double (*R)[7] = reinterpret_cast<double (*)[7]>(r);

  double muk = ooreynolds_mu * (mul + sigma_k * mut);
  double mueps = ooreynolds_mu * (mul + sigma_eps * mut);

  R[0][5] = muk * (dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
		   dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5]);
  R[0][6] = mueps * (dp1dxj[0][0]*V[0][6] + dp1dxj[1][0]*V[1][6] + 
		     dp1dxj[2][0]*V[2][6] + dp1dxj[3][0]*V[3][6]);

  R[1][5] = muk * (dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
		   dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5]);
  R[1][6] = mueps * (dp1dxj[0][1]*V[0][6] + dp1dxj[1][1]*V[1][6] + 
		     dp1dxj[2][1]*V[2][6] + dp1dxj[3][1]*V[3][6]);

  R[2][5] = muk * (dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
		   dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5]);
  R[2][6] = mueps * (dp1dxj[0][2]*V[0][6] + dp1dxj[1][2]*V[1][6] + 
		     dp1dxj[2][2]*V[2][6] + dp1dxj[3][2]*V[3][6]);

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2];
  double div2 = dudxj[0][0]*dudxj[0][0] + dudxj[1][1]*dudxj[1][1] + dudxj[2][2]*dudxj[2][2];
  double a = dudxj[0][1] + dudxj[1][0];
  double b = dudxj[0][2] + dudxj[2][0];
  double c = dudxj[1][2] + dudxj[2][1];
  double prod = ooreynolds_mu * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
    - twothird * rhok * div;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = - rhoeps + prod;
  S[6] = (sigma_eps1 * rhoeps * prod - sigma_eps2 * rhoeps*rhoeps) / rhok;

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

      double RR[9], K[9];
      double alpha[3], beta[3];
      double volten = tetVol *(1.0/10.0);
                                                                                                                                       
      // non-dimensionalization correction //
                                                                                                                                       
      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;
                                                                                                                                      
      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);
                                                                                                                                       
      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;
                                                                                                                                       
      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
                                                                                                                                       
      double SS[4][3];
                                                                                                                                       
      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                       
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                      
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                       
      // FE flux for the porous sink term //
                                                                                                                                       
      for (int j=0; j<12; ++j) PR[j] = 0.0; // initialize PR
                                                                                                                                       
      for (int j=0; j<4; ++j) {
        for (int k=0; k<3; ++k)
          PR[3*j+k] += (K[3*k+0] * SS[j][0] + K[3*k+1] * SS[j][1] + K[3*k+2] * SS[j][2]);
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid 
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut); 
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R); 
     for (int j=0; j<12; ++j) PR[j] = 0.0;                
   }
  }
  else {
     mu = ooreynolds_mu * (mul + mut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
     computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);
  }

  return(porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mul = viscoFcn->compute_mu(Tcg);
  double rhok, rhoeps;
  double mut;
  double lambda;

  mut = computeTurbulentViscosity(V, rhok, rhoeps);

  double mu;
  double kappa;

  int k;
  for (k=0; k<4*3*7*7; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*7*7; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][7][7] = reinterpret_cast<double (*)[3][7][7]>(drdu);
  double (*dSdU)[7][7] = reinterpret_cast<double (*)[7][7]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  computeJacobianVolumeTermKE<7,5>(dp1dxj, mul, mut, rhok, rhoeps, V, dRdU, dSdU);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];
             dPdU[k][j][3][5] =  0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;
             dPdU[k][j][0][5] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];
             dPdU[k][j][1][5] = 0.0;

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
             dPdU[k][j][2][5] = 0.0;
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];
             dPdU[k][j][3][5] = 0.0;

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;
             dPdU[k][j][4][5] = 0.0;

             dPdU[k][j][5][0] = 0.0; 
             dPdU[k][j][5][1] = 0.0;
             dPdU[k][j][5][2] = 0.0;
             dPdU[k][j][5][3] = 0.0;
             dPdU[k][j][5][4] = 0.0;
             dPdU[k][j][5][5] = 0.0;
          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{
  
  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;
  R[6] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*7*7; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[7][7] = reinterpret_cast<double (*)[7][7]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermKEmean::FemEquationTermKEmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

  wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

double FemEquationTermKEmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

bool FemEquationTermKEmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);
                                                                                                                                                                   
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
                                                                                                                                                                   
  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                                                                                   
  double mul = viscoFcn->compute_mu(Tcg);
  double rhok, rhoeps;
  double mut;
  double lambda;

  // Applying the laminar-turbulent trip


  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
       X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
       X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
       X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
       X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
       mut = computeTurbulentViscosity(V, rhok, rhoeps);}
    else{
       mut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, rhok, rhoeps);
  }

  double mu, kappa;

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(volume_id>0) {
    map<int,VolumeData *>::iterator it = volInfo.find(volume_id);
    if(it!=  volInfo.end()) {     // if porous media with volume_id has been defined in the input file
      porousmedia = true;
      double cmu = 0.09;
      double coeff = 1.2247*pow(cmu,0.25);
      double Idr = it->second->idr;                    // average turbulence intensity
      double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
      double vel = sqrt(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2]);

      mut = coeff*Idr*Ldr*vel;
      mu  = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      double RR[9], K[9], B[9];
      double alpha[3], beta[3];

      double volten = tetVol *(1.0/10.0);
      double onefourth = 1.0/4.0;
      double onehalf = 1.0/2.0;
  
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);

      // non-dimensionalization correction //

      alpha[0] = it->second->alphax*length/density;
      alpha[1] = it->second->alphay*length/density;
      alpha[2] = it->second->alphaz*length/density;

      beta[0]  = it->second->betax*length/(density*velocity);
      beta[1]  = it->second->betay*length/(density*velocity);
      beta[2]  = it->second->betaz*length/(density*velocity);

      // transformation matrix //

      RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
      RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
      RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

      // permittivity matrix //

      computePermittivityTensor(alpha, beta, ucg, RR, K);
      
      // gradient of permittivity matrix //

      computeGradPermittivityTensor(alpha, ucg, RR, B);

      double SS[4][3];

      SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
      SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
      SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));
                                                                                                                                                       
      SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
      SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
      SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));
                                                                                                                                                         
      SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
      SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
      SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));
                                                                                                                                                         
      SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
      SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
      SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));
                                                                                                                                                         
      for (int k=0; k<4; ++k) {
        double BB[3];
        BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
        BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
        BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

        double BV[9];
        BV[0] = ucg[0]*BB[0]; BV[1] = ucg[1]*BB[0]; BV[2] = ucg[2]*BB[0]; 
        BV[3] = ucg[0]*BB[1]; BV[4] = ucg[1]*BB[1]; BV[5] = ucg[2]*BB[1]; 
        BV[6] = ucg[0]*BB[2]; BV[7] = ucg[1]*BB[2]; BV[8] = ucg[2]*BB[2]; 
         
        for (int j=0; j<4; ++j) {
          double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};
          double KU[25];
          multiplyBydVdU(v, K, KU, volten);
          double BU[25];
          multiplyBydVdU(v, BV, BU, 1.0);

          if (k == j) {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = KU[5] + BU[5];
             dPdU[k][j][1][1] = KU[6] + BU[6];
             dPdU[k][j][1][2] = KU[7] + BU[7];
             dPdU[k][j][1][3] = KU[8] + BU[8];
             dPdU[k][j][1][4] = KU[9] + BU[9];

             dPdU[k][j][2][0] = KU[10] + BU[10];
             dPdU[k][j][2][1] = KU[11] + BU[11];
             dPdU[k][j][2][2] = KU[12] + BU[12];
             dPdU[k][j][2][3] = KU[13] + BU[13];
             dPdU[k][j][2][4] = KU[14] + BU[14];
          
             dPdU[k][j][3][0] = KU[15] + BU[15];
             dPdU[k][j][3][1] = KU[16] + BU[16];
             dPdU[k][j][3][2] = KU[17] + BU[17];
             dPdU[k][j][3][3] = KU[18] + BU[18];
             dPdU[k][j][3][4] = KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;

          }
          else {
             dPdU[k][j][0][0] = 0.0; 
             dPdU[k][j][0][1] = 0.0;
             dPdU[k][j][0][2] = 0.0;
             dPdU[k][j][0][3] = 0.0;
             dPdU[k][j][0][4] = 0.0;

             dPdU[k][j][1][0] = 0.5*KU[5] + BU[5];
             dPdU[k][j][1][1] = 0.5*KU[6] + BU[6];
             dPdU[k][j][1][2] = 0.5*KU[7] + BU[7];
             dPdU[k][j][1][3] = 0.5*KU[8] + BU[8];
             dPdU[k][j][1][4] = 0.5*KU[9] + BU[9];

             dPdU[k][j][2][0] = 0.5*KU[10] + BU[10];
             dPdU[k][j][2][1] = 0.5*KU[11] + BU[11];
             dPdU[k][j][2][2] = 0.5*KU[12] + BU[12];
             dPdU[k][j][2][3] = 0.5*KU[13] + BU[13];
             dPdU[k][j][2][4] = 0.5*KU[14] + BU[14];
          
             dPdU[k][j][3][0] = 0.5*KU[15] + BU[15];
             dPdU[k][j][3][1] = 0.5*KU[16] + BU[16];
             dPdU[k][j][3][2] = 0.5*KU[17] + BU[17];
             dPdU[k][j][3][3] = 0.5*KU[18] + BU[18];
             dPdU[k][j][3][4] = 0.5*KU[19] + BU[19];

             dPdU[k][j][4][0] = 0.0; 
             dPdU[k][j][4][1] = 0.0;
             dPdU[k][j][4][2] = 0.0;
             dPdU[k][j][4][3] = 0.0;
             dPdU[k][j][4][4] = 0.0;

          }
        }
      }
    }
    else { // if porous media with volume_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermKEmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermKEturb::FemEquationTermKEturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), volInfo(iod.porousmedia.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermKEturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int volume_id)
{

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mul = viscoFcn->compute_mu(Tcg);

  double rhok, rhoeps;
  double mut;

  mut = computeTurbulentViscosity(V, rhok, rhoeps);
 
  double (*dRdU)[3][2][2] = reinterpret_cast<double (*)[3][2][2]>(drdu);
  double (*dSdU)[2][2] = reinterpret_cast<double (*)[2][2]>(dsdu);
  computeJacobianVolumeTermKE<2,0>(dp1dxj, mul, mut, rhok, rhoeps, V, dRdU, dSdU);

  return false;

}

//------------------------------------------------------------------------------
