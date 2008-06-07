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
  NavierStokesTerm(iod, vf), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length; 

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::GAS)
     completeJac = true;
   else
     completeJac = false;
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

// Included (MB)
double FemEquationTermNS::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;
  double mul = ooreynolds_mu * viscoFcn->compute_mu(T);
  double dmul = dooreynolds_mu * viscoFcn->compute_mu(T) + ooreynolds_mu * viscoFcn->compute_muDerivative(T,dT,dMach);
  return (dmul*V[0]-mul*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermNS::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol, 
                                          SVec<double,3> &X, int nodeNum[4],  
                                          int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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

// Included (MB)
bool FemEquationTermNS::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false; 

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  double dooreynolds_lambda = -1.0 / ( reynolds_lambdaNS * reynolds_lambdaNS ) * dRe_lambdadMachNS * dMach;

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);

  double dmu = dooreynolds_mu * viscoFcn->compute_mu(Tcg) + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
  double dkappa = dooreynolds_mu * thermalCondFcn->compute(Tcg) + ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcg, dMach);

  double (*dR)[5] = reinterpret_cast<double (*)[5]>(dr);
  computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
       porousmedia = true;
       fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
       exit(1);
    } 
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
      for (int j=0; j<12; ++j) dPR[j] = 0.0; 
    }
  }
 
  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;

  return (porousmedia);

}

//------------------------------------------------------------------------------

// This function was modified to account for the derivative of the mu with respect to the conservative variables
// Included (MB*)
bool FemEquationTermNS::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  double dmu[4][5];

  double dlambda[4][5];
  
  double dkappa[4][5];

  double dTcgdu0 = 0.0;
  double dTcgdu1 = 0.0;
  double dTcgdu2 = 0.0;
  double dTcgdu3 = 0.0;
  double dTcgdu4 = 0.0;
  
  double dMach = 0.0;

  for (int k=0; k<4; ++k) {

    if (completeJac) { 
      dTcgdu0 = - 0.25 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
      dTcgdu1 = - 0.25 / V[k][0] * V[k][1];
      dTcgdu2 = - 0.25 / V[k][0] * V[k][2];
      dTcgdu3 = - 0.25 / V[k][0] * V[k][3];
      dTcgdu4 = 0.25 / V[k][0];
    }

    dkappa[k][0] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu0, dMach);
    dkappa[k][1] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu1, dMach);
    dkappa[k][2] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu2, dMach);
    dkappa[k][3] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu3, dMach); 
    dkappa[k][4] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu4, dMach); 

    dmu[k][0] = ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcgdu0, dMach);
    dmu[k][1] = ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcgdu1, dMach);
    dmu[k][2] = ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcgdu2, dMach);
    dmu[k][3] = ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcgdu3, dMach); 
    dmu[k][4] = ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcgdu4, dMach); 

    dlambda[k][0] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][0], dMach);
    dlambda[k][1] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][1], dMach);
    dlambda[k][2] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][2], dMach);
    dlambda[k][3] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][3], dMach); 
    dlambda[k][4] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][4], dMach); 

  }
  
  computeJacobianVolumeTermNS(dp1dxj, mu, dmu, lambda, dlambda, kappa, dkappa, V, T, dRdU);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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

// Original
/*
bool FemEquationTermNS::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
      for(int i=0; i<5; ++i)
         for(int j=0; j<5; ++j)
           for(int k=0; k<5; ++k)
             for(int l=0; l<5; ++l)
                dPdU[i][j][k][l] = 0.0;
    }
  }

  return (porousmedia);

}
*/

//------------------------------------------------------------------------------

void FemEquationTermNS::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{
  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermNS::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

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

// Included (MB)
void FemEquationTermNS::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code,
					   Vec3D &n, Vec3D &dn, double d2w[4],
					   double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{

  computeDerivativeOfSurfaceTermNS(dp1dxj, ddp1dxj, n, dn, Vwall, dVwall, V, dV, dMach, dR);

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
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::GAS)
     completeJac = true;
   else
     completeJac = false;

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

// Included (MB)
double FemEquationTermSA::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{
  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermSA::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid 
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

// Included (MB)
bool FemEquationTermSA::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double dmul = viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);

// Test
  double mutilde = 0.0;
  double dmutilde = 0.0;

  double mut;
  double dmut;
  double lambda;
  double dlambda;
  
  // Applying the laminar-turbulent trip
  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
    	X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
    	X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
    	(X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
    	X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
    	X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)) {
    	mut = computeTurbulentViscosity(V, mul, mutilde);
    	dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
    }
    else {
        computeTurbulentViscosity(V, mul, mutilde);
        computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
    	mut = 0.0; 
    	dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul, mutilde);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;
  double dooreynolds_lambda = -1.0 / ( reynolds_lambdaNS * reynolds_lambdaNS ) * dRe_lambdadMachNS * dMach;
  
  double mu;
  double dmu;
  double kappa;
  double dkappa;

  double (*dR)[6] = reinterpret_cast<double (*)[6]>(dr);

  double absmutilde = fabs(mutilde);
  double dabsmutilde;
  if  (mutilde != 0.0)
    dabsmutilde = ( fabs(mutilde) / mutilde ) * dmutilde;
  else {
    fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varible mutilde is zero *****\n");
    //exit(1);
    dabsmutilde = 0.0;
  }

  if (mutilde == 0.0) {
    fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varibles in the function max are equal *****\n");
    //exit(1);
  }
  double maxmutilde = max(mutilde, 0.0);
  double dmaxmutilde;
  if ( maxmutilde == 0.0 )  {
    dmaxmutilde = 0.0;
  }
  else {
    dmaxmutilde = dmutilde;
  }
  double mu5 = oosigma * (mul + absmutilde);
  double dmu5 = oosigma * (dmul + dabsmutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] +
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double ddnutildedx = ddp1dxj[0][0]*V[0][5] + dp1dxj[0][0]*dV[0][5] + ddp1dxj[1][0]*V[1][5] + dp1dxj[1][0]*dV[1][5] +
    ddp1dxj[2][0]*V[2][5] + dp1dxj[2][0]*dV[2][5] + ddp1dxj[3][0]*V[3][5] + dp1dxj[3][0]*dV[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] +
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double ddnutildedy = ddp1dxj[0][1]*V[0][5] + dp1dxj[0][1]*dV[0][5] + ddp1dxj[1][1]*V[1][5] + dp1dxj[1][1]*dV[1][5] +
    ddp1dxj[2][1]*V[2][5] + dp1dxj[2][1]*dV[2][5] + ddp1dxj[3][1]*V[3][5] + dp1dxj[3][1]*dV[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] +
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];
  double ddnutildedz = ddp1dxj[0][2]*V[0][5] + dp1dxj[0][2]*dV[0][5] + ddp1dxj[1][2]*V[1][5] + dp1dxj[1][2]*dV[1][5] +
    ddp1dxj[2][2]*V[2][5] + dp1dxj[2][2]*dV[2][5] + ddp1dxj[3][2]*V[3][5] + dp1dxj[3][2]*dV[3][5];

  dR[0][5] = dmu5 * dnutildedx + mu5 * ddnutildedx;
  dR[1][5] = dmu5 * dnutildedy + mu5 * ddnutildedy;
  dR[2][5] = dmu5 * dnutildedz + mu5 * ddnutildedz;

  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  if (d2wall < 1.e-15) {
    if (mutilde/mul == 0.001) {
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varibles in the function max are equal *****\n");
      //exit(1);
    }
    double chi = max(mutilde/mul, 0.001);
    double dchi;
    if (chi == 0.001)
      dchi = 0.0;
    else
      dchi = ( dmutilde * mul - mutilde * dmul ) / ( mul * mul );
    double chi3 = chi*chi*chi;
    double fv1 = chi3 / (chi3 + cv1_pow3);
    double dfv1 = ( 3.0*chi*chi*dchi*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    double fv2 = 1.0 + oocv2*chi;
    double dfv2 = oocv2*dchi;
    dfv2 = -3.0 / (fv2*fv2*fv2*fv2)*dfv2;
    fv2 = 1.0 / (fv2*fv2*fv2);
    double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
    double dfv3 = ( ( dchi*fv1 + chi*dfv1 ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi ) / ( chi * chi );
    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
    double drho = 0.25 * (dV[0][0] + dV[1][0] + dV[2][0] + dV[3][0]);
    double oorho = 1.0 / rho;
    double doorho = -1.0 / ( rho * rho ) * drho;
    double zz = ooreynolds * oovkcst2 * mutilde * oorho * ood2wall2;
    double dzz = dooreynolds_mu * oovkcst2 * mutilde * oorho * ood2wall2 + ooreynolds_mu * oovkcst2 * dmutilde * oorho * ood2wall2 + ooreynolds * oovkcst2 * mutilde * doorho * ood2wall2;
    double s12 = dudxj[0][1] - dudxj[1][0];
    double ds12 = ddudxj[0][1] - ddudxj[1][0];
    double s23 = dudxj[1][2] - dudxj[2][1];
    double ds23 = ddudxj[1][2] - ddudxj[2][1];
    double s31 = dudxj[2][0] - dudxj[0][2];
    double ds31 = ddudxj[2][0] - ddudxj[0][2];
    double s = sqrt(s12*s12 + s23*s23 + s31*s31);
    double ds = 1.0 / ( 2.0*s ) * (2.0*s12*ds12 + 2.0*s23*ds23 + 2.0*s31*ds31);
    double Stilde = s*fv3 + zz*fv2;
    double dStilde = ds*fv3 + s*dfv3 + dzz*fv2 + zz*dfv2;
    double rr = min(zz/Stilde, 2.0);
    double drr;
    if (rr==2.0)
      drr = 0.0;
    else
      drr = ( dzz * Stilde - zz*dStilde ) / ( Stilde * Stilde );
    double rr2 = rr*rr;
    double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
    double dgg = drr + cw2 * (6.0*rr*rr2*rr2*drr - drr);
    double gg2 = gg*gg;
    double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);
    double dfw = opcw3_pow * dgg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg;

//  double AA = oosigma * cb2 * rho *
//    (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double dAA = oosigma * cb2 * drho * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz) +
                            oosigma * cb2 * rho * (2.0*dnutildedx*ddnutildedx + 2.0*dnutildedy*ddnutildedy + 2.0*dnutildedz*ddnutildedz);
//  double BB = cb1 * Stilde * absmutilde;
    double dBB = cb1 * dStilde * absmutilde + cb1 * Stilde * dabsmutilde;
//  double CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    double dCC = - cw1 * dfw * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde * ood2wall2;
    dS[5] = dAA + dBB + dCC;
  }
  else {
    dS[5] = 0.0;
  }


  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      porousmedia = true;
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
      exit(1);
    }
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid 
      mu = ooreynolds_mu * (mul + mut);
      dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      dkappa = dooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut) + ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcg, dMach) + alpha * dmut);
      computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
      for (int j=0; j<12; ++j) dPR[j] = 0.0;                
   }
  }
  else {
    mu = ooreynolds_mu * (mul + mut);
    dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);
    lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
    dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
    kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
    dkappa = dooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut) + ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcg, dMach) + alpha * dmut);
    computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
  }

}

//------------------------------------------------------------------------------

// This function was modified to account for the derivative of the mu with respect to the conservative variables
// Included (MB*)
bool FemEquationTermSA::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde = 0.0;
  double mut;
  double lambda;

  mut = computeTurbulentViscosity(V, mul, mutilde);

  double mu;
  double kappa;

  mu = ooreynolds_mu * (mul + mut);
  lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
  kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);

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

  double dmul[4][6];
  double dmutilde[4][6];
  double dmut[4][6];
  double dmu[4][6];
  double dlambda[4][6];
  double dkappa[4][6];

  double dTcgdu0 = 0.0;
  double dTcgdu1 = 0.0;
  double dTcgdu2 = 0.0;
  double dTcgdu3 = 0.0;
  double dTcgdu4 = 0.0;
  double dTcgdu5 = 0.0;
  
  double dMach = 0.0;

  for (int k=0; k<4; ++k) {

    if (completeJac) {
      dTcgdu0 = - 0.25 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
      dTcgdu1 = - 0.25 / V[k][0] * V[k][1];
      dTcgdu2 = - 0.25 / V[k][0] * V[k][2];
      dTcgdu3 = - 0.25 / V[k][0] * V[k][3];
      dTcgdu4 = 0.25 / V[k][0];
      dTcgdu5 = 0.0;
    }

    dmul[k][0] = viscoFcn->compute_muDerivative(Tcg, dTcgdu0, dMach);
    dmul[k][1] = viscoFcn->compute_muDerivative(Tcg, dTcgdu1, dMach);
    dmul[k][2] = viscoFcn->compute_muDerivative(Tcg, dTcgdu2, dMach);
    dmul[k][3] = viscoFcn->compute_muDerivative(Tcg, dTcgdu3, dMach); 
    dmul[k][4] = viscoFcn->compute_muDerivative(Tcg, dTcgdu4, dMach); 
    dmul[k][5] = viscoFcn->compute_muDerivative(Tcg, dTcgdu5, dMach); 

    dmutilde[k][0] = 0.0;
    dmutilde[k][1] = 0.0;
    dmutilde[k][2] = 0.0;
    dmutilde[k][3] = 0.0;
    dmutilde[k][4] = 0.0;
    if (completeJac) 
      dmutilde[k][5] = 0.25;
    else
      dmutilde[k][5] = 0.0;

    dmut[k][0] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][0], dmutilde[k][0]);
    dmut[k][1] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][1], dmutilde[k][1]);
    dmut[k][2] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][2], dmutilde[k][2]);
    dmut[k][3] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][3], dmutilde[k][3]);
    dmut[k][4] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][4], dmutilde[k][4]);
    dmut[k][5] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][5], dmutilde[k][5]);

    dkappa[k][0] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu0, dMach) + alpha * dmut[k][0]);
    dkappa[k][1] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu1, dMach) + alpha * dmut[k][1]);
    dkappa[k][2] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu2, dMach) + alpha * dmut[k][2]);
    dkappa[k][3] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu3, dMach) + alpha * dmut[k][3]); 
    dkappa[k][4] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu4, dMach) + alpha * dmut[k][4]); 
    dkappa[k][5] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu5, dMach) + alpha * dmut[k][5]); 

    dmu[k][0] = ooreynolds_mu * (dmul[k][0] + dmut[k][0]);
    dmu[k][1] = ooreynolds_mu * (dmul[k][1] + dmut[k][1]);
    dmu[k][2] = ooreynolds_mu * (dmul[k][2] + dmut[k][2]);
    dmu[k][3] = ooreynolds_mu * (dmul[k][3] + dmut[k][3]); 
    dmu[k][4] = ooreynolds_mu * (dmul[k][4] + dmut[k][4]); 
    dmu[k][5] = ooreynolds_mu * (dmul[k][5] + dmut[k][5]); 

    dlambda[k][0] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][0], dMach);
    dlambda[k][1] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][1], dMach);
    dlambda[k][2] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][2], dMach);
    dlambda[k][3] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][3], dMach); 
    dlambda[k][4] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][4], dMach); 
    dlambda[k][5] = ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu[k][5], dMach); 

  }
  
  double s = sqrt((dudxj[0][1] - dudxj[1][0])*(dudxj[0][1] - dudxj[1][0]) + (dudxj[1][2] - dudxj[2][1])*(dudxj[1][2] - dudxj[2][1]) + (dudxj[2][0] - dudxj[0][2])*(dudxj[2][0] - dudxj[0][2]));
  if (s != 0.0)
    computeJacobianVolumeTermSA<6,5>(dp1dxj, d2w, dudxj, mul, dmul, mutilde, dmutilde, V, dRdU, dSdU);
  else
    computeJacobianVolumeTermSA<6,5>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU);
     
  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
      mu = ooreynolds_mu * (mul + mut);
      lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
      kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
      computeJacobianVolumeTermNS(dp1dxj, mu, dmu, lambda, dlambda, kappa, dkappa, V, T, dRdU);

      for(int i=0; i<5; ++i)
        for(int j=0; j<5; ++j)
          for(int k=0; k<5; ++k)
            for(int l=0; l<5; ++l)
               dPdU[i][j][k][l] = 0.0;
    }
  }
  else {
    computeJacobianVolumeTermNS(dp1dxj, mu, dmu, lambda, dlambda, kappa, dkappa, V, T, dRdU);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

// Original
/*
bool FemEquationTermSA::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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
*/

//------------------------------------------------------------------------------

void FemEquationTermSA::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{

  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermSA::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;

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

// Included (MB)
void FemEquationTermSA::computeBCsJacobianWallValues(int code, Vec3D &n, double d2w[3], double *Vwall, double *dVwall, double *V[3])
{

  wallFcn->computeBCsJacobianWallValues<6>(code, n, d2w, Vwall, dVwall, V);

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

// Included (MB)
void FemEquationTermSA::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code,
					   Vec3D &n, Vec3D &dn, double d2w[4],
					   double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{

  computeDerivativeOfSurfaceTermNS(dp1dxj, ddp1dxj, n, dn, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;

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
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
double FemEquationTermSAmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermSAmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermSAturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::GAS)
     completeJac = true;
   else
     completeJac = false;

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

// Included (MB)
double FemEquationTermDES::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid 
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

// Included (MB)
bool FemEquationTermDES::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
  exit(1);

  return (porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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

// Included (MB)
void FemEquationTermDES::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{
  fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to FemEquationTermDES::computeDerivativeOfSurfaceTerm is not implemented *****\n");
  exit(1);
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

// Included (MB)
void FemEquationTermDES::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code, Vec3D &n, Vec3D &dn, double d2w[4], double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{
  fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to FemEquationTermDES::computeDerivativeOfSurfaceTerm is not implemented *****\n");
  exit(1);
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
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
double FemEquationTermDESmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermDESmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermDESturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::GAS)
     completeJac = true;
   else
     completeJac = false;

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

// Included (MB)
double FemEquationTermKE::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V);
       dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V);
    dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid 
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

// Included (MB)
bool FemEquationTermKE::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double dmul = viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double rhok, rhoeps;
  double mut;
  double lambda;

  double drhok, drhoeps;
  double dmut;
  double dlambda;

  // Applying the laminar-turbulent trip
  if(trip){

    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
    	X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
    	X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
    	(X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
    	X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
    	X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)) {
        mut = computeTurbulentViscosity(V, rhok, rhoeps);
        dmut = computeDerivativeOfTurbulentViscosity(V, dV, drhok, drhoeps, dMach);     
    }
    else{
        computeTurbulentViscosity(V, rhok, rhoeps);
        computeDerivativeOfTurbulentViscosity(V, dV, drhok, drhoeps, dMach);     
        mut = 0.0; 
        dmut = 0.0;
    }
  }
  else{
    mut = computeTurbulentViscosity(V, rhok, rhoeps);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, drhok, drhoeps, dMach);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  double dooreynolds_lambda = -1.0 / ( reynolds_lambdaNS * reynolds_lambdaNS ) * dRe_lambdadMachNS * dMach;

  double mu;
  double dmu;
  double kappa;
  double dkappa;

  double (*dR)[7] = reinterpret_cast<double (*)[7]>(dr);

  double muk = ooreynolds_mu * (mul + sigma_k * mut);
  double dmuk = dooreynolds_mu * (mul + sigma_k * mut) + ooreynolds_mu * (dmul + sigma_k * dmut);
  double mueps = ooreynolds_mu * (mul + sigma_eps * mut);
  double dmueps = dooreynolds_mu * (mul + sigma_eps * mut) + ooreynolds_mu * (dmul + sigma_eps * dmut);

  dR[0][5] = dmuk * (dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] +
		   dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5]) +
           muk * (ddp1dxj[0][0]*V[0][5] + dp1dxj[0][0]*dV[0][5] + ddp1dxj[1][0]*V[1][5] + dp1dxj[1][0]*dV[1][5] +
		   ddp1dxj[2][0]*V[2][5] + dp1dxj[2][0]*dV[2][5] + ddp1dxj[3][0]*V[3][5] + dp1dxj[3][0]*dV[3][5]);
  dR[0][6] = dmueps * (dp1dxj[0][0]*V[0][6] + dp1dxj[1][0]*V[1][6] +
		     dp1dxj[2][0]*V[2][6] + dp1dxj[3][0]*V[3][6]) +
             mueps * (ddp1dxj[0][0]*V[0][6] + dp1dxj[0][0]*dV[0][6] + ddp1dxj[1][0]*V[1][6] + dp1dxj[1][0]*dV[1][6] +
		     ddp1dxj[2][0]*V[2][6] + dp1dxj[2][0]*dV[2][6] + ddp1dxj[3][0]*V[3][6] + dp1dxj[3][0]*dV[3][6]);

  dR[1][5] = dmuk * (dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] +
		   dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5]) +
           muk * (ddp1dxj[0][1]*V[0][5] + dp1dxj[0][1]*dV[0][5] + ddp1dxj[1][1]*V[1][5] + dp1dxj[1][1]*dV[1][5] +
		   ddp1dxj[2][1]*V[2][5] + dp1dxj[2][1]*dV[2][5] + ddp1dxj[3][1]*V[3][5] + dp1dxj[3][1]*dV[3][5]);
  dR[1][6] = dmueps * (dp1dxj[0][1]*V[0][6] + dp1dxj[1][1]*V[1][6] +
		     dp1dxj[2][1]*V[2][6] + dp1dxj[3][1]*V[3][6]) +
             mueps * (ddp1dxj[0][1]*V[0][6] + dp1dxj[0][1]*dV[0][6] + ddp1dxj[1][1]*V[1][6] + dp1dxj[1][1]*dV[1][6] +
		     ddp1dxj[2][1]*V[2][6] + dp1dxj[2][1]*dV[2][6] + ddp1dxj[3][1]*V[3][6] + dp1dxj[3][1]*dV[3][6]);

  dR[2][5] = dmuk * (dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] +
		   dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5]) +
           muk * (ddp1dxj[0][2]*V[0][5] + dp1dxj[0][2]*dV[0][5] + ddp1dxj[1][2]*V[1][5] + dp1dxj[1][2]*dV[1][5] +
		   ddp1dxj[2][2]*V[2][5] + dp1dxj[2][2]*dV[2][5] + ddp1dxj[3][2]*V[3][5] + dp1dxj[3][2]*dV[3][5]);
  dR[2][6] = dmueps * (dp1dxj[0][2]*V[0][6] + dp1dxj[1][2]*V[1][6] +
		     dp1dxj[2][2]*V[2][6] + dp1dxj[3][2]*V[3][6]) +
             mueps * (ddp1dxj[0][2]*V[0][6] + dp1dxj[0][2]*dV[0][6] + ddp1dxj[1][2]*V[1][6] + dp1dxj[1][2]*dV[1][6] +
		     ddp1dxj[2][2]*V[2][6] + dp1dxj[2][2]*dV[2][6] + ddp1dxj[3][2]*V[3][6] + dp1dxj[3][2]*dV[3][6]);

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2];
  double ddiv = ddudxj[0][0] + ddudxj[1][1] + ddudxj[2][2];
  double div2 = dudxj[0][0]*dudxj[0][0] + dudxj[1][1]*dudxj[1][1] + dudxj[2][2]*dudxj[2][2];
  double ddiv2 = 2.0*dudxj[0][0]*ddudxj[0][0] + 2.0*dudxj[1][1]*ddudxj[1][1] + 2.0*dudxj[2][2]*ddudxj[2][2];
  double a = dudxj[0][1] + dudxj[1][0];
  double da = ddudxj[0][1] + ddudxj[1][0];
  double b = dudxj[0][2] + dudxj[2][0];
  double db = ddudxj[0][2] + ddudxj[2][0];
  double c = dudxj[1][2] + dudxj[2][1];
  double dc = ddudxj[1][2] + ddudxj[2][1];
  double prod = ooreynolds * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
    - twothird * rhok * div;

  double dprod = dooreynolds_mu * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
                 - twothird * rhok * div + ooreynolds * dmut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c) +
                 ooreynolds_mu * mut * (2.0 * ddiv2 - twothird * 2.0*div*ddiv + 2.0*a*da + 2.0*b*db + 2.0*c*dc)
                 - twothird * drhok * div - twothird * rhok * ddiv;

  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;
  dS[5] = - drhoeps + dprod;
  dS[6] = ( (sigma_eps1 * drhoeps * prod + sigma_eps1 * rhoeps * dprod  - sigma_eps2 * 2.0*rhoeps*drhoeps) * rhok -  (sigma_eps1 * rhoeps * prod - sigma_eps2 * rhoeps*rhoeps) * drhok ) / ( rhok * rhok );

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      porousmedia = true;
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
      exit(1);

    }
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid 
     mu = ooreynolds_mu * (mul + mut);
     dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);												    
     dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
     dkappa = dooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut) + ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcg, dMach) + alpha * dmut);
     computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
     for (int j=0; j<12; ++j) dPR[j] = 0.0;                
   }
  }
  else {
     mu = ooreynolds_mu * (mul + mut);
     dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);
     lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);
     dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
     kappa = ooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut);
     dkappa = dooreynolds_mu * (thermalCondFcn->compute(Tcg) + alpha * mut) + ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcg, dMach) + alpha * dmut);
     computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
  }

  return(porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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

// Included (MB)
void FemEquationTermKE::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;
  dR[6] = 0.0;

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
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
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

// Included (MB)
double FemEquationTermKEmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V);
       dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V);
    dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermKEmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
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
    else { // if porous media with material_id has NOT been defined in the input file -> treat as standard fluid
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
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermKEturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
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
