#include <VarFcn.h>

#include <DistVector.h>
#include <Vector3D.h>


#include<math.h>
#include<complex.h>
typedef complex<double> bcomp;


//------------------------------------------------------------------------------

VarFcn::VarFcn(IoData &iod)
{
  pmin  = iod.eqs.fluidModel.pmin;
  pminp = iod.eqs.fluidModel2.pmin;
  verif_clipping = false;
  node_change = true;

  gravity = 0.0;
  ngravity[0] = cos(iod.bc.hydro.alpha)*cos(iod.bc.hydro.beta);
  ngravity[1] = cos(iod.bc.hydro.alpha)*sin(iod.bc.hydro.beta);
  ngravity[2] = sin(iod.bc.hydro.alpha);
  if (iod.bc.hydro.type == BcsHydroData::GRAVITY)
    gravity = iod.bc.hydro.gravity;

  meshVel = 0.0;
}

//------------------------------------------------------------------------------
void VarFcn::computeNewPrimitiveLiquidEuler(double Pr, double a, double b, double *U, double *V)
{
  V[0] = Pr+a*pow(U[0],b);
  V[1] = U[1];
  V[2] = U[2];
  V[3] = U[3];
  V[4] = U[4];
}
void VarFcn::computeOldPrimitiveLiquidEuler(double Pr, double a, double b, double *U, double *V)
{
  V[0] = pow((U[0]-Pr)/a,b);
  V[1] = U[1];
  V[2] = U[2];
  V[3] = U[3];
  V[4] = U[4];
}

void VarFcn::conservativeToPrimitiveGasEuler(double g, double Ps, double *U, double *V)
{
                         
  V[0] = U[0];
           
  double invRho = 1.0 / U[0];
  if (isnan(invRho)){
        fprintf(stderr, "ERROR*** conservativeToPrimitive invRho = 1 / %f\n", U[0]);
        exit(1);
  }
                                                                     
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
  if (isnan(V[1])){
        fprintf(stderr, "ERROR*** conservativeToPrimitive u = %e / %e\n", U[1], U[0]);
        exit(1);
  }
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                    
  V[4] = (g-1.0) * (U[4] - 0.5 * U[0] * vel2) - g*Ps;
                                                                     
}
                                                                   
//------------------------------------------------------------------------------

void VarFcn::primitiveToConservativeGasEuler(double g, double invg1, double Ps,
                                             double *V, double *U)
{
                                                           
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                                                 
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = (V[4]+g*Ps) * invg1 + 0.5 * V[0] * vel2;
                                            
}

//------------------------------------------------------------------------------
//FLAG_ARL:  so be careful for multiphase flows. - Any problem with the norm of n??
//           or with the stiffened gas formulation???
void VarFcn::primitiveToCharacteristicVariationsGasEuler(double n[3],
                                 double g, double Ps, 
                                 double *V, double *dV, double *dW)
{
  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double ooc2 = V[0]/(g*(V[4]+Ps));
  double oorhoc = sqrt(ooc2)/V[0];
  double coeff1 = dV[0] - ooc2*dV[4];

  dW[0] = n[0]*coeff1 + n[2]*dV[2] - n[1]*dV[3];
  dW[1] = n[1]*coeff1 + n[0]*dV[3] - n[2]*dV[1];
  dW[2] = n[2]*coeff1 + n[1]*dV[1] - n[0]*dV[2];
  dW[3] = dVn + oorhoc*dV[4];
  dW[4] =-dVn + oorhoc*dV[4];
}
                                                                                                    
//-----------------------------------------------------------------------------
//FLAG_ARL:  so be careful for multiphase flows.
void VarFcn::characteristicToPrimitiveVariationsGasEuler(double n[3],
                                 double g, double Ps, 
                                 double *V, double *dW, double *dV)
{
  double sum = dW[3]+dW[4]; 
  double diff = dW[3]-dW[4];
  double c = sqrt(g*(V[4]+Ps)/V[0]);

  dV[0] = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2] + 0.5*V[0]*sum/c;
  dV[1] = n[1]*dV[2]-n[2]*dV[1] + 0.5*n[0]*diff;
  dV[1] = n[2]*dV[0]-n[0]*dV[2] + 0.5*n[1]*diff;
  dV[1] = n[0]*dV[1]-n[1]*dV[0] + 0.5*n[2]*diff;
  dV[4] = 0.5*V[0]*c*sum;
                                                                                                    
}
                                                                                                    
//-----------------------------------------------------------------------------

void VarFcn::computedVdUGasEuler(double g1, double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  dVdU[20] = g1 * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dVdU[21] = -g1 * V[1];
  dVdU[22] = -g1 * V[2];
  dVdU[23] = -g1 * V[3];
  dVdU[24] = g1;
}

//------------------------------------------------------------------------------

void VarFcn::computedUdVGasEuler(double invg1, double *V, double *dUdV)
{

  dUdV[0]  = 1.0;
  dUdV[5]  = V[1];
  dUdV[6]  = V[0];
  dUdV[10] = V[2];
  dUdV[12] = V[0];
  dUdV[15] = V[3];
  dUdV[18] = V[0];
  dUdV[20] = 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[21] = V[0] * V[1];
  dUdV[22] = V[0] * V[2];
  dUdV[23] = V[0] * V[3];
  dUdV[24] = invg1;

}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUGasEuler(double g1, double *V, double *vec, double *res)
{
  double dVdU[25];
  computedVdUGasEuler(g1, V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*vec[4];

}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUGasEuler(double g1, double *V, bcomp *vec, bcomp *res)
{
  double dVdU[25];
  computedVdUGasEuler(g1, V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*
vec[4];

}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUTGasEuler(double g1,  double *V, double *vec, double *res)
{
  double dVdU[25];
  computedVdUGasEuler(g1, V, dVdU);

  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];

}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUTGasEuler(double g1, double *V, bcomp *vec, bcomp *res)
{
  double dVdU[25];
  computedVdUGasEuler(g1, V, dVdU);
                                                                                                                                                                                                     
  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];
}

//------------------------------------------------------------------------------

void VarFcn::preMultiplyBydUdVGasEuler(double invg1, double *V, double *mat, double *res)
{
  double dUdV[25];
  computedUdVGasEuler(invg1, V, dUdV);

  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[5]*mat[0]+dUdV[6]*mat[5];
  res[6] = dUdV[5]*mat[1]+dUdV[6]*mat[6];
  res[7] = dUdV[5]*mat[2]+dUdV[6]*mat[7];
  res[8] = dUdV[5]*mat[3]+dUdV[6]*mat[8];
  res[9] = dUdV[5]*mat[4]+dUdV[6]*mat[9];
  res[10] = dUdV[10]*mat[0]+dUdV[12]*mat[10];
  res[11] = dUdV[10]*mat[1]+dUdV[12]*mat[11];
  res[12] = dUdV[10]*mat[2]+dUdV[12]*mat[12];
  res[13] = dUdV[10]*mat[3]+dUdV[12]*mat[13];
  res[14] = dUdV[10]*mat[4]+dUdV[12]*mat[14];
  res[15] = dUdV[15]*mat[0]+dUdV[18]*mat[15];
  res[16] = dUdV[15]*mat[1]+dUdV[18]*mat[16];
  res[17] = dUdV[15]*mat[2]+dUdV[18]*mat[17];
  res[18] = dUdV[15]*mat[3]+dUdV[18]*mat[18];
  res[19] = dUdV[15]*mat[4]+dUdV[18]*mat[19];
  res[20] = dUdV[20]*mat[0]+dUdV[21]*mat[5]+dUdV[22]*mat[10]+dUdV[23]*mat[15]+dUdV[24]*mat[20];
  res[21] = dUdV[20]*mat[1]+dUdV[21]*mat[6]+dUdV[22]*mat[11]+dUdV[23]*mat[16]+dUdV[24]*mat[21];
  res[22] = dUdV[20]*mat[2]+dUdV[21]*mat[7]+dUdV[22]*mat[12]+dUdV[23]*mat[17]+dUdV[24]*mat[22];
  res[23] = dUdV[20]*mat[3]+dUdV[21]*mat[8]+dUdV[22]*mat[13]+dUdV[23]*mat[18]+dUdV[24]*mat[23];
  res[24] = dUdV[20]*mat[4]+dUdV[21]*mat[9]+dUdV[22]*mat[14]+dUdV[23]*mat[19]+dUdV[24]*mat[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydVdUGasEuler(double g1, double *V, double *mat, double *res)
{
  double dVdU[25];
  computedVdUGasEuler(g1, V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[5]+mat[2]*dVdU[10]+mat[3]*dVdU[15]+mat[4]*dVdU[20];
  res[1] = mat[1]*dVdU[6]+mat[4]*dVdU[21];
  res[2] = mat[2]*dVdU[12]+mat[4]*dVdU[22];
  res[3] = mat[3]*dVdU[18]+mat[4]*dVdU[23];
  res[4] = mat[4]*dVdU[24];
  res[5] = mat[5]*dVdU[0]+mat[6]*dVdU[5]+mat[7]*dVdU[10]+mat[8]*dVdU[15]+mat[9]*dVdU[20];
  res[6] = mat[6]*dVdU[6]+mat[9]*dVdU[21];
  res[7] = mat[7]*dVdU[12]+mat[9]*dVdU[22];
  res[8] = mat[8]*dVdU[18]+mat[9]*dVdU[23];
  res[9] = mat[9]*dVdU[24];
  res[10] = mat[10]*dVdU[0]+mat[11]*dVdU[5]+mat[12]*dVdU[10]+mat[13]*dVdU[15]+mat[14]*dVdU[20];  res[11] = mat[11]*dVdU[6]+mat[14]*dVdU[21];
  res[12] = mat[12]*dVdU[12]+mat[14]*dVdU[22];
  res[13] = mat[13]*dVdU[18]+mat[14]*dVdU[23];
  res[14] = mat[14]*dVdU[24];
  res[15] = mat[15]*dVdU[0]+mat[16]*dVdU[5]+mat[17]*dVdU[10]+mat[18]*dVdU[15]+mat[19]*dVdU[20];  res[16] = mat[16]*dVdU[6]+mat[19]*dVdU[21];
  res[17] = mat[17]*dVdU[12]+mat[19]*dVdU[22];
  res[18] = mat[18]*dVdU[18]+mat[19]*dVdU[23];
  res[19] = mat[19]*dVdU[24];
  res[20] = mat[20]*dVdU[0]+mat[21]*dVdU[5]+mat[22]*dVdU[10]+mat[23]*dVdU[15]+mat[24]*dVdU[20];  res[21] = mat[21]*dVdU[6]+mat[24]*dVdU[21];
  res[22] = mat[22]*dVdU[12]+mat[24]*dVdU[22];
  res[23] = mat[23]*dVdU[18]+mat[24]*dVdU[23];
  res[24] = mat[24]*dVdU[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVGasEuler(double invg1, double *V, double *mat, double *res)
{
  double dUdV[25];
  computedUdVGasEuler(invg1, V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVGasEuler(double invg1, double *V, bcomp *mat, bcomp *res)
{
  double dUdV[25];
  computedUdVGasEuler(invg1, V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];
}

//------------------------------------------------------------------------------

void VarFcn::extrapolatePrimitiveGasEuler(double un, double c, double *Vb,
                                         double *Vinter, double *V)
{

  if (un == 0.0){
    V[0] = Vb[0];
    V[1] = Vb[1];
    V[2] = Vb[2];
    V[3] = Vb[3];
    V[4] = Vinter[4];
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        V[0] = Vb[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }else{                 //    SUBSONIC
        V[0] = Vb[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vinter[4];
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
        V[0] = Vinter[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }else{                //     SUBSONIC
        V[0] = Vinter[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vb[4];
      }
    }
  }
}

//------------------------------------------------------------------------------

void VarFcn::extrapolateCharacteristicGasEuler(double gam, double Ps, double n[3],
                                        double un, double c, double *Vb,
                                        double *dV)
{
/* routine computes boundary conditions using characteristic methods
 * and assuming that values are small perturbations of values at infinity
 * initially dV contains perturbations of primitive variables to be extrapolated
 * at return, dV contains perturbations of primitive variables at the boundary
 */

  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double ooc2 = 1.0/(c*c);
  double oorhoc = sqrt(ooc2)/Vb[0];
  double coeff1 = dV[0] - ooc2*dV[4];

// step 1: primitive to characteristic variations
  double dW[5];
  dW[0] = n[0]*coeff1 + n[2]*dV[2] - n[1]*dV[3];
  dW[1] = n[1]*coeff1 + n[0]*dV[3] - n[2]*dV[1];
  dW[2] = n[2]*coeff1 + n[1]*dV[1] - n[0]*dV[2];
  dW[3] = dVn + oorhoc*dV[4];
  dW[4] =-dVn + oorhoc*dV[4];

// step 2: choose variations to be extrapolated
//         if incoming characteristic, then the
//         characteristic is set to 0.0
//         else, there is nothing to do
  if (un == 0.0){
    dW[0] = 0.0;
    dW[1] = 0.0;
    dW[2] = 0.0;
    dW[3] = 0.0;
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
        dW[4] = 0.0;
      }else{                 //    SUBSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
      }else{                //     SUBSONIC
        dW[4] = 0.0;
      }
    }
  }

// step 3: characteristic to primitive variations
  double sum = dW[3]+dW[4];
  double diff = dW[3]-dW[4];

  dV[0] = dW[0]*n[0]+dW[1]*n[1]+dW[2]*n[2] + 0.5*Vb[0]*sum/c;
  dV[1] = n[1]*dW[2]-n[2]*dW[1] + 0.5*n[0]*diff;
  dV[2] = n[2]*dW[0]-n[0]*dW[2] + 0.5*n[1]*diff;
  dV[3] = n[0]*dW[1]-n[1]*dW[0] + 0.5*n[2]*diff;
  dV[4] = 0.5*Vb[0]*c*sum;

}

//------------------------------------------------------------------------------

void VarFcn::extrapolatePrimitiveGasEulerPUT(double un, double c, double *Vb,
                                         double *Vinter, double *V)
{
  exit(1);
  double gam1 = getGamma1();
  if (getGamma1() == 0.0){
    fprintf(stderr, "extrapolatePrimitiveGasEulerPUT cannot be used. Exiting.\n");
    exit(1);
  }

  double Pb = Vb[4];
  double Rb = Vb[0];
  double Tb = Pb/(Rb*gam1);
  
  double P = Vinter[4];
  double R = Vinter[0]; 
  double T = P/(gam1*R);

  if (un == 0.0){
    V[0] = P/(gam1*Tb);
    V[1] = Vb[1];
    V[2] = Vb[2];
    V[3] = Vb[3];
    V[4] = P;
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        V[0] = Vb[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }else{                 //    SUBSONIC
        V[0] = P/(gam1*Tb);
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = P;
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
        V[0] = Vinter[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }else{                //     SUBSONIC
        V[0] = Pb/(gam1*T);;
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Pb;
      }
    }
  }
}
//---------------------------------------------------------------------------

int VarFcn::VerificationGasEuler(int glob, double pmin, double gam, double Pstiff,
	     double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  if(V[4]<pmin){
    if (verif_clipping)
      fprintf(stdout, "clip pressure[%d] in gas from %e to %e\n", glob, V[4], pmin);
    V[4] = pmin;
    U[4] = 0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
          +(pmin+gam*Pstiff)/(gam-1.0);
    return 1;
  }
  return 0;

}

//---------------------------------------------------------------------------

void VarFcn::conservativeToPrimitiveGasSA(double g1, double *U, double *V)
{
  V[0] = U[0];

  double invRho = 1.0 / U[0];

  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
  V[4] = g1 * (U[4] - 0.5 * U[0] * vel2);
  V[5] = U[5] * invRho;
}
//------------------------------------------------------------------------------

void VarFcn::primitiveToConservativeGasSA(double *V, double *U)
{
  fprintf(stderr, "*** Error: primitiveToConservative not implemented\n");
}
//------------------------------------------------------------------------------

void VarFcn::computedUdVGasSA(double invg1, double *V, double *dUdV)
{
  dUdV[0]  = 1.0;
  dUdV[6]  = V[1];
  dUdV[7]  = V[0];
  dUdV[12] = V[2];
  dUdV[14] = V[0];
  dUdV[18] = V[3];
  dUdV[21] = V[0];
  dUdV[24] = 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[25] = V[0] * V[1];
  dUdV[26] = V[0] * V[2];
  dUdV[27] = V[0] * V[3];
  dUdV[28] = invg1;
  dUdV[30] = V[5];
  dUdV[35] = V[0];
}
//------------------------------------------------------------------------------

void VarFcn::computedVdUGasSA(double g1, double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[6]  = -invrho * V[1];
  dVdU[7]  = invrho;
  dVdU[12] = -invrho * V[2];
  dVdU[14] = invrho;
  dVdU[18] = -invrho * V[3];
  dVdU[21] = invrho;
  dVdU[24] = g1 * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dVdU[25] = -g1 * V[1];
  dVdU[26] = -g1 * V[2];
  dVdU[27] = -g1 * V[3];
  dVdU[28] = g1;
  dVdU[30] = -invrho * V[5];
  dVdU[35] = invrho;
}
//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUGasSA(double *V, double *vec, double *res)
{
  fprintf(stderr, "*** Error: multiplyBydVdU not implemented\n");
}
//------------------------------------------------------------------------------

void VarFcn::preMultiplyBydUdVGasSA(double invg1, double *V, double *mat, double *res)
{
  double dUdV[36];
  computedUdVGasSA(invg1, V, dUdV);
  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[0]*mat[5];
  res[6] = dUdV[6]*mat[0]+dUdV[7]*mat[6];
  res[7] = dUdV[6]*mat[1]+dUdV[7]*mat[7];
  res[8] = dUdV[6]*mat[2]+dUdV[7]*mat[8];
  res[9] = dUdV[6]*mat[3]+dUdV[7]*mat[9];
  res[10] = dUdV[6]*mat[4]+dUdV[7]*mat[10];
  res[11] = dUdV[6]*mat[5]+dUdV[7]*mat[11];
  res[12] = dUdV[12]*mat[0]+dUdV[14]*mat[12];
  res[13] = dUdV[12]*mat[1]+dUdV[14]*mat[13];
  res[14] = dUdV[12]*mat[2]+dUdV[14]*mat[14];
  res[15] = dUdV[12]*mat[3]+dUdV[14]*mat[15];
  res[16] = dUdV[12]*mat[4]+dUdV[14]*mat[16];
  res[17] = dUdV[12]*mat[5]+dUdV[14]*mat[17];
  res[18] = dUdV[18]*mat[0]+dUdV[21]*mat[18];
  res[19] = dUdV[18]*mat[1]+dUdV[21]*mat[19];
  res[20] = dUdV[18]*mat[2]+dUdV[21]*mat[20];
  res[21] = dUdV[18]*mat[3]+dUdV[21]*mat[21];
  res[22] = dUdV[18]*mat[4]+dUdV[21]*mat[22];
  res[23] = dUdV[18]*mat[5]+dUdV[21]*mat[23];
  res[24] = dUdV[24]*mat[0]+dUdV[25]*mat[6]+dUdV[26]*mat[12]+dUdV[27]*mat[18]+dUdV[28]*mat[24];
  res[25] = dUdV[24]*mat[1]+dUdV[25]*mat[7]+dUdV[26]*mat[13]+dUdV[27]*mat[19]+dUdV[28]*mat[25];
  res[26] = dUdV[24]*mat[2]+dUdV[25]*mat[8]+dUdV[26]*mat[14]+dUdV[27]*mat[20]+dUdV[28]*mat[26];
  res[27] = dUdV[24]*mat[3]+dUdV[25]*mat[9]+dUdV[26]*mat[15]+dUdV[27]*mat[21]+dUdV[28]*mat[27];
  res[28] = dUdV[24]*mat[4]+dUdV[25]*mat[10]+dUdV[26]*mat[16]+dUdV[27]*mat[22]+dUdV[28]*mat[28];
  res[29] = dUdV[24]*mat[5]+dUdV[25]*mat[11]+dUdV[26]*mat[17]+dUdV[27]*mat[23]+dUdV[28]*mat[29];
  res[30] = dUdV[30]*mat[0]+dUdV[35]*mat[30];
  res[31] = dUdV[30]*mat[1]+dUdV[35]*mat[31];
  res[32] = dUdV[30]*mat[2]+dUdV[35]*mat[32];
  res[33] = dUdV[30]*mat[3]+dUdV[35]*mat[33];
  res[34] = dUdV[30]*mat[4]+dUdV[35]*mat[34];
  res[35] = dUdV[30]*mat[5]+dUdV[35]*mat[35];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydVdUGasSA(double g1, double *V, double *mat, double *res)
{
  double dVdU[36];
  computedVdUGasSA(g1, V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[6]+mat[2]*dVdU[12]+
    mat[3]*dVdU[18]+mat[4]*dVdU[24]+mat[5]*dVdU[30];
  res[1] = mat[1]*dVdU[7]+mat[4]*dVdU[25];
  res[2] = mat[2]*dVdU[14]+mat[4]*dVdU[26];
  res[3] = mat[3]*dVdU[21]+mat[4]*dVdU[27];
  res[4] = mat[4]*dVdU[28];
  res[5] = mat[5]*dVdU[35];
  res[6] = mat[6]*dVdU[0]+mat[7]*dVdU[6]+mat[8]*dVdU[12]+
    mat[9]*dVdU[18]+mat[10]*dVdU[24]+mat[11]*dVdU[30];
  res[7] = mat[7]*dVdU[7]+mat[10]*dVdU[25];
  res[8] = mat[8]*dVdU[14]+mat[10]*dVdU[26];
  res[9] = mat[9]*dVdU[21]+mat[10]*dVdU[27];
  res[10] = mat[10]*dVdU[28];
  res[11] = mat[11]*dVdU[35];
  res[12] = mat[12]*dVdU[0]+mat[13]*dVdU[6]+mat[14]*dVdU[12]+
    mat[15]*dVdU[18]+mat[16]*dVdU[24]+mat[17]*dVdU[30];
  res[13] = mat[13]*dVdU[7]+mat[16]*dVdU[25];
  res[14] = mat[14]*dVdU[14]+mat[16]*dVdU[26];
  res[15] = mat[15]*dVdU[21]+mat[16]*dVdU[27];
  res[16] = mat[16]*dVdU[28];
  res[17] = mat[17]*dVdU[35];
  res[18] = mat[18]*dVdU[0]+mat[19]*dVdU[6]+mat[20]*dVdU[12]+
    mat[21]*dVdU[18]+mat[22]*dVdU[24]+mat[23]*dVdU[30];
  res[19] = mat[19]*dVdU[7]+mat[22]*dVdU[25];
  res[20] = mat[20]*dVdU[14]+mat[22]*dVdU[26];
  res[21] = mat[21]*dVdU[21]+mat[22]*dVdU[27];
  res[22] = mat[22]*dVdU[28];
  res[23] = mat[23]*dVdU[35];
  res[24] = mat[24]*dVdU[0]+mat[25]*dVdU[6]+mat[26]*dVdU[12]+
    mat[27]*dVdU[18]+mat[28]*dVdU[24]+mat[29]*dVdU[30];
  res[25] = mat[25]*dVdU[7]+mat[28]*dVdU[25];
  res[26] = mat[26]*dVdU[14]+mat[28]*dVdU[26];
  res[27] = mat[27]*dVdU[21]+mat[28]*dVdU[27];
  res[28] = mat[28]*dVdU[28];
  res[29] = mat[29]*dVdU[35];
  res[30] = mat[30]*dVdU[0]+mat[31]*dVdU[6]+mat[32]*dVdU[12]+
    mat[33]*dVdU[18]+mat[34]*dVdU[24]+mat[35]*dVdU[30];
  res[31] = mat[31]*dVdU[7]+mat[34]*dVdU[25];
  res[32] = mat[32]*dVdU[14]+mat[34]*dVdU[26];
  res[33] = mat[33]*dVdU[21]+mat[34]*dVdU[27];
  res[34] = mat[34]*dVdU[28];
  res[35] = mat[35]*dVdU[35];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVGasSA(double invg1, double *V, double *mat, double *res)
{
  double dUdV[36];
  computedUdVGasSA(invg1, V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[6]+mat[2]*dUdV[12]+
    mat[3]*dUdV[18]+mat[4]*dUdV[24]+mat[5]*dUdV[30];
  res[1] = mat[1]*dUdV[7]+mat[4]*dUdV[25];
  res[2] = mat[2]*dUdV[14]+mat[4]*dUdV[26];
  res[3] = mat[3]*dUdV[21]+mat[4]*dUdV[27];
  res[4] = mat[4]*dUdV[28];
  res[5] = mat[5]*dUdV[35];
  res[6] = mat[6]*dUdV[0]+mat[7]*dUdV[6]+mat[8]*dUdV[12]+
    mat[9]*dUdV[18]+mat[10]*dUdV[24]+mat[11]*dUdV[30];
  res[7] = mat[7]*dUdV[7]+mat[10]*dUdV[25];
  res[8] = mat[8]*dUdV[14]+mat[10]*dUdV[26];
  res[9] = mat[9]*dUdV[21]+mat[10]*dUdV[27];
  res[10] = mat[10]*dUdV[28];
  res[11] = mat[11]*dUdV[35];
  res[12] = mat[12]*dUdV[0]+mat[13]*dUdV[6]+mat[14]*dUdV[12]+
    mat[15]*dUdV[18]+mat[16]*dUdV[24]+mat[17]*dUdV[30];
  res[13] = mat[13]*dUdV[7]+mat[16]*dUdV[25];
  res[14] = mat[14]*dUdV[14]+mat[16]*dUdV[26];
  res[15] = mat[15]*dUdV[21]+mat[16]*dUdV[27];
  res[16] = mat[16]*dUdV[28];
  res[17] = mat[17]*dUdV[35];
  res[18] = mat[18]*dUdV[0]+mat[19]*dUdV[6]+mat[20]*dUdV[12]+
    mat[21]*dUdV[18]+mat[22]*dUdV[24]+mat[23]*dUdV[30];
  res[19] = mat[19]*dUdV[7]+mat[22]*dUdV[25];
  res[20] = mat[20]*dUdV[14]+mat[22]*dUdV[26];
  res[21] = mat[21]*dUdV[21]+mat[22]*dUdV[27];
  res[22] = mat[22]*dUdV[28];
  res[23] = mat[23]*dUdV[35];
  res[24] = mat[24]*dUdV[0]+mat[25]*dUdV[6]+mat[26]*dUdV[12]+
    mat[27]*dUdV[18]+mat[28]*dUdV[24]+mat[29]*dUdV[30];
  res[25] = mat[25]*dUdV[7]+mat[28]*dUdV[25];
  res[26] = mat[26]*dUdV[14]+mat[28]*dUdV[26];
  res[27] = mat[27]*dUdV[21]+mat[28]*dUdV[27];
  res[28] = mat[28]*dUdV[28];
  res[29] = mat[29]*dUdV[35];
  res[30] = mat[30]*dUdV[0]+mat[31]*dUdV[6]+mat[32]*dUdV[12]+
    mat[33]*dUdV[18]+mat[34]*dUdV[24]+mat[35]*dUdV[30];
  res[31] = mat[31]*dUdV[7]+mat[34]*dUdV[25];
  res[32] = mat[32]*dUdV[14]+mat[34]*dUdV[26];
  res[33] = mat[33]*dUdV[21]+mat[34]*dUdV[27];
  res[34] = mat[34]*dUdV[28];
  res[35] = mat[35]*dUdV[35];
  //fprintf(stderr, "*** Error: postMultiplyBydUdV not implemented\n"); //ARL
}
//---------------------------------------------------------------------------

int VarFcn::VerificationGasSA(int glob, double pmin, double gam1, double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  if(V[4]<pmin){
    V[4] = pmin;
    U[4] = 0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
          +pmin/gam1;
    return 1;
  }
  return 0;

}

//------------------------------------------------------------------------------

void VarFcn::conservativeToPrimitiveGasKE(double g1, double *U, double *V)
{
  V[0] = U[0];
  double invRho = 1.0 / U[0];
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
  V[4] = g1 * (U[4] - 0.5 * U[0] * vel2);
  V[5] = U[5] * invRho;
  V[6] = U[6] * invRho;
}

//------------------------------------------------------------------------------

void VarFcn::primitiveToConservativeGasKE(double *V, double *U)
{
  fprintf(stderr, "*** Error: primitiveToConservative not implemented\n");
}
//------------------------------------------------------------------------------

void VarFcn::computedUdVGasKE(double invg1, double *V, double *dUdV)
{
  dUdV[0]  = 1.0;
  dUdV[7]  = V[1];
  dUdV[8]  = V[0];
  dUdV[14] = V[2];
  dUdV[16] = V[0];
  dUdV[21] = V[3];
  dUdV[24] = V[0];
  dUdV[28] = 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[29] = V[0] * V[1];
  dUdV[30] = V[0] * V[2];
  dUdV[31] = V[0] * V[3];
  dUdV[32] = invg1;
  dUdV[35] = V[5];
  dUdV[40] = V[0];
  dUdV[42] = V[6];
  dUdV[48] = V[0];
}

//------------------------------------------------------------------------------

void VarFcn::computedVdUGasKE(double g1, double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[7]  = -invrho * V[1];
  dVdU[8]  = invrho;
  dVdU[14] = -invrho * V[2];
  dVdU[16] = invrho;
  dVdU[21] = -invrho * V[3];
  dVdU[24] = invrho;
  dVdU[28] = g1 * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dVdU[29] = -g1 * V[1];
  dVdU[30] = -g1 * V[2];
  dVdU[31] = -g1 * V[3];
  dVdU[32] = g1;
  dVdU[35] = -invrho * V[5];
  dVdU[40] = invrho;
  dVdU[42] = -invrho * V[6];
  dVdU[48] = invrho;
}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUGasKE(double *V, double *vec, double *res)
{
  fprintf(stderr, "*** Error: multiplyBydVdU not implemented\n");
}

//------------------------------------------------------------------------------

void VarFcn::preMultiplyBydUdVGasKE(double invg1, double *V, double *mat, double *res)
{
  double dUdV[49];
  computedUdVGasKE(invg1, V, dUdV);

  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[0]*mat[5];
  res[6] = dUdV[0]*mat[6];
  res[7] = dUdV[7]*mat[0]+dUdV[8]*mat[7];
  res[8] = dUdV[7]*mat[1]+dUdV[8]*mat[8];
  res[9] = dUdV[7]*mat[2]+dUdV[8]*mat[9];
  res[10] = dUdV[7]*mat[3]+dUdV[8]*mat[10];
  res[11] = dUdV[7]*mat[4]+dUdV[8]*mat[11];
  res[12] = dUdV[7]*mat[5]+dUdV[8]*mat[12];
  res[13] = dUdV[7]*mat[6]+dUdV[8]*mat[13];
  res[14] = dUdV[14]*mat[0]+dUdV[16]*mat[14];
  res[15] = dUdV[14]*mat[1]+dUdV[16]*mat[15];
  res[16] = dUdV[14]*mat[2]+dUdV[16]*mat[16];
  res[17] = dUdV[14]*mat[3]+dUdV[16]*mat[17];
  res[18] = dUdV[14]*mat[4]+dUdV[16]*mat[18];
  res[19] = dUdV[14]*mat[5]+dUdV[16]*mat[19];
  res[20] = dUdV[14]*mat[6]+dUdV[16]*mat[20];
  res[21] = dUdV[21]*mat[0]+dUdV[24]*mat[21];
  res[22] = dUdV[21]*mat[1]+dUdV[24]*mat[22];
  res[23] = dUdV[21]*mat[2]+dUdV[24]*mat[23];
  res[24] = dUdV[21]*mat[3]+dUdV[24]*mat[24];
  res[25] = dUdV[21]*mat[4]+dUdV[24]*mat[25];
  res[26] = dUdV[21]*mat[5]+dUdV[24]*mat[26];
  res[27] = dUdV[21]*mat[6]+dUdV[24]*mat[27];
  res[28] = dUdV[28]*mat[0]+dUdV[29]*mat[7]+dUdV[30]*mat[14]+dUdV[31]*mat[21]+dUdV[32]*mat[28];
  res[29] = dUdV[28]*mat[1]+dUdV[29]*mat[8]+dUdV[30]*mat[15]+dUdV[31]*mat[22]+dUdV[32]*mat[29];
  res[30] = dUdV[28]*mat[2]+dUdV[29]*mat[9]+dUdV[30]*mat[16]+dUdV[31]*mat[23]+dUdV[32]*mat[30];
  res[31] = dUdV[28]*mat[3]+dUdV[29]*mat[10]+dUdV[30]*mat[17]+dUdV[31]*mat[24]+dUdV[32]*mat[31];
  res[32] = dUdV[28]*mat[4]+dUdV[29]*mat[11]+dUdV[30]*mat[18]+dUdV[31]*mat[25]+dUdV[32]*mat[32];
  res[33] = dUdV[28]*mat[5]+dUdV[29]*mat[12]+dUdV[30]*mat[19]+dUdV[31]*mat[26]+dUdV[32]*mat[33];
  res[34] = dUdV[28]*mat[6]+dUdV[29]*mat[13]+dUdV[30]*mat[20]+dUdV[31]*mat[27]+dUdV[32]*mat[34];
  res[35] = dUdV[35]*mat[0]+dUdV[40]*mat[35];
  res[36] = dUdV[35]*mat[1]+dUdV[40]*mat[36];
  res[37] = dUdV[35]*mat[2]+dUdV[40]*mat[37];
  res[38] = dUdV[35]*mat[3]+dUdV[40]*mat[38];
  res[39] = dUdV[35]*mat[4]+dUdV[40]*mat[39];
  res[40] = dUdV[35]*mat[5]+dUdV[40]*mat[40];
  res[41] = dUdV[35]*mat[6]+dUdV[40]*mat[41];
  res[42] = dUdV[42]*mat[0]+dUdV[48]*mat[42];
  res[43] = dUdV[42]*mat[1]+dUdV[48]*mat[43];
  res[44] = dUdV[42]*mat[2]+dUdV[48]*mat[44];
  res[45] = dUdV[42]*mat[3]+dUdV[48]*mat[45];
  res[46] = dUdV[42]*mat[4]+dUdV[48]*mat[46];
  res[47] = dUdV[42]*mat[5]+dUdV[48]*mat[47];
  res[48] = dUdV[42]*mat[6]+dUdV[48]*mat[48];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydVdUGasKE(double g1, double *V, double *mat, double *res)
{
  double dVdU[49];
  computedVdUGasKE(g1, V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[7]+mat[2]*dVdU[14]+mat[3]*dVdU[21]+
    mat[4]*dVdU[28]+mat[5]*dVdU[35]+mat[6]*dVdU[42];
  res[1] = mat[1]*dVdU[8]+mat[4]*dVdU[29];
  res[2] = mat[2]*dVdU[16]+mat[4]*dVdU[30];
  res[3] = mat[3]*dVdU[24]+mat[4]*dVdU[31];
  res[4] = mat[4]*dVdU[32];
  res[5] = mat[5]*dVdU[40];
  res[6] = mat[6]*dVdU[48];
  res[7] = mat[7]*dVdU[0]+mat[8]*dVdU[7]+mat[9]*dVdU[14]+mat[10]*dVdU[21]+
    mat[11]*dVdU[28]+mat[12]*dVdU[35]+mat[13]*dVdU[42];
  res[8] = mat[8]*dVdU[8]+mat[11]*dVdU[29];
  res[9] = mat[9]*dVdU[16]+mat[11]*dVdU[30];
  res[10] = mat[10]*dVdU[24]+mat[11]*dVdU[31];
  res[11] = mat[11]*dVdU[32];
  res[12] = mat[12]*dVdU[40];
  res[13] = mat[13]*dVdU[48];
  res[14] = mat[14]*dVdU[0]+mat[15]*dVdU[7]+mat[16]*dVdU[14]+mat[17]*dVdU[21]+
    mat[18]*dVdU[28]+mat[19]*dVdU[35]+mat[20]*dVdU[42];
  res[15] = mat[15]*dVdU[8]+mat[18]*dVdU[29];
  res[16] = mat[16]*dVdU[16]+mat[18]*dVdU[30];
  res[17] = mat[17]*dVdU[24]+mat[18]*dVdU[31];
  res[18] = mat[18]*dVdU[32];
  res[19] = mat[19]*dVdU[40];
  res[20] = mat[20]*dVdU[48];
  res[21] = mat[21]*dVdU[0]+mat[22]*dVdU[7]+mat[23]*dVdU[14]+mat[24]*dVdU[21]+
    mat[25]*dVdU[28]+mat[26]*dVdU[35]+mat[27]*dVdU[42];
  res[22] = mat[22]*dVdU[8]+mat[25]*dVdU[29];
  res[23] = mat[23]*dVdU[16]+mat[25]*dVdU[30];
  res[24] = mat[24]*dVdU[24]+mat[25]*dVdU[31];
  res[25] = mat[25]*dVdU[32];
  res[26] = mat[26]*dVdU[40];
  res[27] = mat[27]*dVdU[48];
  res[28] = mat[28]*dVdU[0]+mat[29]*dVdU[7]+mat[30]*dVdU[14]+mat[31]*dVdU[21]+
    mat[32]*dVdU[28]+mat[33]*dVdU[35]+mat[34]*dVdU[42];
  res[29] = mat[29]*dVdU[8]+mat[32]*dVdU[29];
  res[30] = mat[30]*dVdU[16]+mat[32]*dVdU[30];
  res[31] = mat[31]*dVdU[24]+mat[32]*dVdU[31];
  res[32] = mat[32]*dVdU[32];
  res[33] = mat[33]*dVdU[40];
  res[34] = mat[34]*dVdU[48];
  res[35] = mat[35]*dVdU[0]+mat[36]*dVdU[7]+mat[37]*dVdU[14]+mat[38]*dVdU[21]+
    mat[39]*dVdU[28]+mat[40]*dVdU[35]+mat[41]*dVdU[42];
  res[36] = mat[36]*dVdU[8]+mat[39]*dVdU[29];
  res[37] = mat[37]*dVdU[16]+mat[39]*dVdU[30];
  res[38] = mat[38]*dVdU[24]+mat[39]*dVdU[31];
  res[39] = mat[39]*dVdU[32];
  res[40] = mat[40]*dVdU[40];
  res[41] = mat[41]*dVdU[48];
  res[42] = mat[42]*dVdU[0]+mat[43]*dVdU[7]+mat[44]*dVdU[14]+mat[45]*dVdU[21]+
    mat[46]*dVdU[28]+mat[47]*dVdU[35]+mat[48]*dVdU[42];
  res[43] = mat[43]*dVdU[8]+mat[46]*dVdU[29];
  res[44] = mat[44]*dVdU[16]+mat[46]*dVdU[30];
  res[45] = mat[45]*dVdU[24]+mat[46]*dVdU[31];
  res[46] = mat[46]*dVdU[32];
  res[47] = mat[47]*dVdU[40];
  res[48] = mat[48]*dVdU[48];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVGasKE(double invg1, double *V, double *mat, double *res)
{
  double dUdV[49];
  computedUdVGasKE(invg1, V, dUdV);

  for (int i=0; i<7; i++){
    res[0+7*i] = mat[0+7*i]*dUdV[0]+mat[1+7*i]*dUdV[7]+mat[2+7*i]*dUdV[14]+
      mat[3+7*i]*dUdV[21]+mat[4+7*i]*dUdV[28]+mat[5+7*i]*dUdV[35]+mat[6+7*i]*dUdV[42];
    res[1+7*i] = mat[1+7*i]*dUdV[8]+mat[4+7*i]*dUdV[29];
    res[2+7*i] = mat[2+7*i]*dUdV[16]+mat[4+7*i]*dUdV[30];
    res[3+7*i] = mat[3+7*i]*dUdV[24]+mat[4+7*i]*dUdV[31];
    res[4+7*i] = mat[4+7*i]*dUdV[32];
    res[5+7*i] = mat[5+7*i]*dUdV[40];
    res[6+7*i] = mat[6+7*i]*dUdV[48];
  }

  //fprintf(stderr, "*** Error: postMultiplyBydUdV not implemented\n"); //ARL
}

//---------------------------------------------------------------------------

int VarFcn::VerificationGasKE(int glob, double pmin, double gam1, double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  if(V[4]<pmin){
    V[4] = pmin;
    U[4] = 0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
          +pmin/gam1;
    return 1;
  }
  return 0;

}




//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//                                                                            //
//                  BAROTROPIC LIQUID                                         //
//                                                                            //
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void VarFcn::conservativeToPrimitiveLiquidEuler(double invCv, double *U, double *V)
{
  V[0] = U[0];
                                       
  double invRho = 1.0 / U[0];
  if (isnan(invRho)){
        fprintf(stderr, "ERROR*** conservativeToPrimitive\n");
        exit(1);
  }
   
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
      
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
    
  V[4] = invRho * invCv * (U[4] - 0.5 * U[0] * vel2);
}
      
//------------------------------------------------------------------------------

void VarFcn::primitiveToConservativeLiquidEuler(double Cv, double *V, double *U)
{
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                                            
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = V[0] * Cv * V[4] + 0.5 * V[0] * vel2;
}

//-----------------------------------------------------------------------------

void VarFcn::primitiveToCharacteristicVariationsLiquidEuler(double n[3],
                                 double Cv, double a, double b, double Pr,
                                 double *V, double *dV, double *dW)
{

  double Ptemp = -(Pr+a*pow(V[0],b))/pow(V[0],2);
  double c = computeSoundSpeed(V)/V[0];
  double coeff1 = Ptemp*dV[0]+ Cv*dV[4];
  double dVn = n[0]*dV[1]+n[1]*dV[2]+n[2]*dV[3];

  dW[0] = coeff1*n[0] - n[2]*dV[2] + n[1]*dV[3];
  dW[1] = coeff1*n[1] - n[0]*dV[3] + n[2]*dV[1];
  dW[2] = coeff1*n[2] - n[1]*dV[1] + n[0]*dV[2];
  dW[3] = 0.5*(c*dV[0]+dVn);
  dW[4] = 0.5*(c*dV[0]-dVn);
}

//-----------------------------------------------------------------------------

void VarFcn::characteristicToPrimitiveVariationsLiquidEuler(double n[3],
                                 double Cv, double a, double b, double Pr,
                                 double *V, double *dW, double *dV)
{
  double ooCv = 1.0/Cv;
  double c = computeSoundSpeed(V);
  double Ptemp = (Pr+a*pow(V[0],b))/(V[0]*Cv*c);
  double rho = V[0]/c;

  dV[0] = rho*(dW[3]+dW[4]);
  dV[1] = n[2]*dV[1]-n[1]*dV[2] + n[0]*(dW[3]-dW[4]);
  dV[1] = n[0]*dV[2]-n[2]*dV[0] + n[1]*(dW[3]-dW[4]);
  dV[1] = n[1]*dV[0]-n[0]*dV[1] + n[2]*(dW[3]-dW[4]);
  dV[4] = ooCv*(n[0]*dW[0]+n[1]*dW[1]+n[2]*dW[2]) + Ptemp*(dW[3]+dW[4]);

}

//-----------------------------------------------------------------------------

void VarFcn::computedVdULiquidEuler(double invCv, double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  double invrhoCv = invrho*invCv;
  double Cv = 1.0/invCv;
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  dVdU[20] = invrhoCv * ( -Cv*V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) );
  dVdU[21] = -invrhoCv * V[1];
  dVdU[22] = -invrhoCv * V[2];
  dVdU[23] = -invrhoCv * V[3];
  dVdU[24] = invrhoCv;
}

//------------------------------------------------------------------------------

void VarFcn::computedUdVLiquidEuler(double Cv, double *V, double *dUdV)
{
  dUdV[0]  = 1.0;
  dUdV[5]  = V[1];
  dUdV[6]  = V[0];
  dUdV[10] = V[2];
  dUdV[12] = V[0];
  dUdV[15] = V[3];
  dUdV[18] = V[0];
  dUdV[20] = Cv * V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[21] = V[0] * V[1];
  dUdV[22] = V[0] * V[2];
  dUdV[23] = V[0] * V[3];
  dUdV[24] = V[0] * Cv;
}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdULiquidEuler(double invCv, double *V, double *vec, double *res)
{
  double dVdU[25];
  computedVdULiquidEuler(invCv, V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*vec[4];
}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdULiquidEuler(double invCv, double *V, bcomp *vec, bcomp *res)
{
  double dVdU[25];
  computedVdULiquidEuler(invCv, V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*
vec[4];
}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUTLiquidEuler(double invCv, double *V, double *vec, double *res)
{
  double dVdU[25];
  computedVdULiquidEuler(invCv, V, dVdU);

  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];
}

//------------------------------------------------------------------------------

void VarFcn::multiplyBydVdUTLiquidEuler(double invCv, double *V, bcomp *vec, bcomp *res)
{
  double dVdU[25];
  computedVdULiquidEuler(invCv, V, dVdU);

  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];
}

//------------------------------------------------------------------------------

void VarFcn::preMultiplyBydUdVLiquidEuler(double Cv, double *V, double *mat, double *res)
{
  double dUdV[25];
  computedUdVLiquidEuler(Cv, V, dUdV);

  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[5]*mat[0]+dUdV[6]*mat[5];
  res[6] = dUdV[5]*mat[1]+dUdV[6]*mat[6];
  res[7] = dUdV[5]*mat[2]+dUdV[6]*mat[7];
  res[8] = dUdV[5]*mat[3]+dUdV[6]*mat[8];
  res[9] = dUdV[5]*mat[4]+dUdV[6]*mat[9];
  res[10] = dUdV[10]*mat[0]+dUdV[12]*mat[10];
  res[11] = dUdV[10]*mat[1]+dUdV[12]*mat[11];
  res[12] = dUdV[10]*mat[2]+dUdV[12]*mat[12];
  res[13] = dUdV[10]*mat[3]+dUdV[12]*mat[13];
  res[14] = dUdV[10]*mat[4]+dUdV[12]*mat[14];
  res[15] = dUdV[15]*mat[0]+dUdV[18]*mat[15];
  res[16] = dUdV[15]*mat[1]+dUdV[18]*mat[16];
  res[17] = dUdV[15]*mat[2]+dUdV[18]*mat[17];
  res[18] = dUdV[15]*mat[3]+dUdV[18]*mat[18];
  res[19] = dUdV[15]*mat[4]+dUdV[18]*mat[19];
  res[20] = dUdV[20]*mat[0]+dUdV[21]*mat[5]+dUdV[22]*mat[10]+dUdV[23]*mat[15]+dUdV[24]*mat[20];
  res[21] = dUdV[20]*mat[1]+dUdV[21]*mat[6]+dUdV[22]*mat[11]+dUdV[23]*mat[16]+dUdV[24]*mat[21];
  res[22] = dUdV[20]*mat[2]+dUdV[21]*mat[7]+dUdV[22]*mat[12]+dUdV[23]*mat[17]+dUdV[24]*mat[22];
  res[23] = dUdV[20]*mat[3]+dUdV[21]*mat[8]+dUdV[22]*mat[13]+dUdV[23]*mat[18]+dUdV[24]*mat[23];
  res[24] = dUdV[20]*mat[4]+dUdV[21]*mat[9]+dUdV[22]*mat[14]+dUdV[23]*mat[19]+dUdV[24]*mat[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydVdULiquidEuler(double invCv, double *V, double *mat, double *res)
{
  double dVdU[25];
  computedVdULiquidEuler(invCv, V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[5]+mat[2]*dVdU[10]+mat[3]*dVdU[15]+mat[4]*dVdU[20];
  res[1] = mat[1]*dVdU[6]+mat[4]*dVdU[21];
  res[2] = mat[2]*dVdU[12]+mat[4]*dVdU[22];
  res[3] = mat[3]*dVdU[18]+mat[4]*dVdU[23];
  res[4] = mat[4]*dVdU[24];
  res[5] = mat[5]*dVdU[0]+mat[6]*dVdU[5]+mat[7]*dVdU[10]+mat[8]*dVdU[15]+mat[9]*dVdU[20];
  res[6] = mat[6]*dVdU[6]+mat[9]*dVdU[21];
  res[7] = mat[7]*dVdU[12]+mat[9]*dVdU[22];
  res[8] = mat[8]*dVdU[18]+mat[9]*dVdU[23];
  res[9] = mat[9]*dVdU[24];
  res[10] = mat[10]*dVdU[0]+mat[11]*dVdU[5]+mat[12]*dVdU[10]+mat[13]*dVdU[15]+mat[14]*dVdU[20];  res[11] = mat[11]*dVdU[6]+mat[14]*dVdU[21];
  res[12] = mat[12]*dVdU[12]+mat[14]*dVdU[22];
  res[13] = mat[13]*dVdU[18]+mat[14]*dVdU[23];
  res[14] = mat[14]*dVdU[24];
  res[15] = mat[15]*dVdU[0]+mat[16]*dVdU[5]+mat[17]*dVdU[10]+mat[18]*dVdU[15]+mat[19]*dVdU[20];  res[16] = mat[16]*dVdU[6]+mat[19]*dVdU[21];
  res[17] = mat[17]*dVdU[12]+mat[19]*dVdU[22];
  res[18] = mat[18]*dVdU[18]+mat[19]*dVdU[23];
  res[19] = mat[19]*dVdU[24];
  res[20] = mat[20]*dVdU[0]+mat[21]*dVdU[5]+mat[22]*dVdU[10]+mat[23]*dVdU[15]+mat[24]*dVdU[20];  res[21] = mat[21]*dVdU[6]+mat[24]*dVdU[21];
  res[22] = mat[22]*dVdU[12]+mat[24]*dVdU[22];
  res[23] = mat[23]*dVdU[18]+mat[24]*dVdU[23];
  res[24] = mat[24]*dVdU[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVLiquidEuler(double Cv, double *V, double *mat, double *res)
{
  double dUdV[25];
  computedUdVLiquidEuler(Cv, V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];
}

//------------------------------------------------------------------------------

void VarFcn::postMultiplyBydUdVLiquidEuler(double Cv, double *V, bcomp *mat, bcomp *res)
{
  double dUdV[25];
  computedUdVLiquidEuler(Cv, V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];
}
   
//------------------------------------------------------------------------------
void VarFcn::extrapolatePrimitiveLiquidEuler(double un, double c, double *Vb,
                                            double *Vinter, double *V)
{
  if (un == 0.0){
    V[0] = Vinter[0];
    V[1] = Vb[1];
    V[2] = Vb[2];
    V[3] = Vb[3];
    V[4] = Vb[4];
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        V[0] = Vb[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }else{                 //    SUBSONIC
        V[0] = Vinter[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
        V[0] = Vinter[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }else{                //     SUBSONIC
        V[0] = Vb[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }
    }
  }
}

//---------------------------------------------------------------------------

void VarFcn::extrapolateCharacteristicLiquidEuler(double cv, double Pr,
					double a, double b,
					double n[3], double un, double c,
                                        double *Vb, double *dV)
{
//cf Research/latex/notes/matrices, and look at L and L^{-1}
/* routine computes boundary conditions using characteristic methods
 * and assuming that values are small perturbations of values at infinity
 * initially dV contains perturbations of primitive variables to be extrapolated
 * at return, dV contains perturbations of primitive variables at the boundary
 */

  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double rhooc = Vb[0]/c;
  double P = Pr +a*pow(Vb[0],b);
  double coeff1 = -P*dV[0]/Vb[0] + Vb[0]*cv*dV[4];

// step 1: primitive to characteristic variations
  double dW[5];
  dW[0] = coeff1*n[0] - Vb[0]*(n[2]*dV[2]-n[1]*dV[3]);
  dW[1] = coeff1*n[1] - Vb[0]*(n[0]*dV[3]-n[2]*dV[1]);
  dW[2] = coeff1*n[2] - Vb[0]*(n[1]*dV[0]-n[0]*dV[2]);
  dW[3] = 0.5*(dV[0] + rhooc*dVn);
  dW[4] = 0.5*(dV[0] - rhooc*dVn);

// step 2: choose variations to be extrapolated
//         if incoming characteristic, then the
//         characteristic is set to 0.0
//         else, there is nothing to do
  if (un == 0.0){
    dW[0] = 0.0;
    dW[1] = 0.0;
    dW[2] = 0.0;
    dW[3] = 0.0;
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
        dW[4] = 0.0;
      }else{                 //    SUBSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
      }else{                //     SUBSONIC
        dW[4] = 0.0;
      }
    }
  }

// step 3: characteristic to primitive variations
  double sum  = dW[3]+dW[4];
  double diff = dW[3]-dW[4];
  double corho = 1.0/rhooc;
  double oorho = 1.0/Vb[0];

  dV[0] = sum;
  dV[1] = oorho*(n[2]*dW[1]-n[1]*dW[2]) + corho*n[0]*diff;
  dV[2] = oorho*(n[0]*dW[2]-n[2]*dW[0]) + corho*n[1]*diff;
  dV[3] = oorho*(n[1]*dW[0]-n[0]*dW[1]) + corho*n[2]*diff;
  dV[4] = ( dW[0]*n[0]+dW[1]*n[1]+dW[2]*n[2] + P*sum*oorho)*oorho/cv;

}
//---------------------------------------------------------------------------
int VarFcn::VerificationLiquidEuler(int glob, double pmin, double a, double b, double Pr,
	     double oocv, double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  double locPressure = Pr+a*pow(V[0],b);
  if(locPressure<pmin){
    if(verif_clipping)
      fprintf(stdout, "clip pressure[%d] in tait from %e to %e\n", glob, locPressure, pmin);
    V[0] = pow((pmin-Pr)/a, 1.0/b);
    U[0] = V[0];
    U[1] = V[0]*V[1];
    U[2] = V[0]*V[2];
    U[3] = V[0]*V[3];
    U[4] = V[0]*V[4]/oocv+0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    return 1;
  }
  return 0;

}

//------------------------------------------------------------------------------
