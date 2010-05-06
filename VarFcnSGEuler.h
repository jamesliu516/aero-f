#ifndef _VAR_FCN_SGEULER_H
#define _VAR_FCN_SGEULER_H

#include <VarFcnBase.h>

//--------------------------------------------------------------------------
// This class is the VarFcn class for the Stiffened Gas EOS in Euler 
// Equations. Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state. 
//
// lay-out of the base class is:
//  - 1 -  Transformation Operators
//  - 2 -  General Functions
//  - 3 -  Equations of State Parameters
//  - 4 -  EOS related functions
//
//--------------------------------------------------------------------------
//
// EOS: Pressure = (gam - 1)*Density*e - gam*Pc
// where
//   e  : internal energy per unit mass.
//   Pc : pressure constant.
//
//--------------------------------------------------------------------------
class VarFcnSGEuler : public VarFcnBase {

private:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;
  double dPstiff;

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);
  int verification(int glob, double *U, double *V);

public:
  VarFcnSGEuler(FluidModelData &data);
  ~VarFcnSGEuler() { delete [] pname; }

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *);
  void primitiveToConservativeDerivative(double *, double *, double *, double *);
  
  void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V);
  void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV);
  void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW);
  void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV);

  //----- General Functions -----//
  double checkPressure(double *V) const { 
    return V[4]+Pstiff; 
  }
  double computeTemperature(double *V) const {
    if (isnan(1.0/V[0])) {
      fprintf(stderr, "ERROR*** computeTemp\n");
      exit(1);
    }
    return invgam1 * (V[4]+gam*Pstiff) / V[0];
  }
  double computeRhoEnergy(double *V) const {
    return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V) const { 
    return invgam1 * (V[4]+gam*Pstiff); 
  }
  double computeSoundSpeed(double *V) const { 
    return sqrt(gam * (V[4]+Pstiff) / V[0]); 
  }
  double computeSoundSpeed(double density, double entropy) const {
    double c2 = gam * entropy*pow(density,gam-1.0);
    if(c2>0) return sqrt(c2);
    return 0.0; 
  }
  double computeEntropy(double density, double pressure) const {
    return (pressure+Pstiff)/pow(density,gam); 
  }
  double computeIsentropicPressure(double entropy, double density) const {
    return entropy*pow(density,gam)-Pstiff; 
  }
  double computePressureCoefficient(double *V, double pinfty, double mach, bool dimFlag) const {
    if (dimFlag)
      return 2.0 * (V[4] - pinfty);
    else
      return 2.0 * (V[4] - 1.0/(gam*mach*mach));
  }
  double computeTotalPressure(double machr, double* V) const {
    double mach = computeMachNumber(V);
    double machr2 = machr*machr;
    double popr = V[4]*gam*machr2;
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    return V[4]*pow(opmach, gam*invgam1);
  }
  double computeDerivativeOfTemperature(double *V, double *dV) const {
    return invgam1 * ( dV[4]*V[0] - V[4]*dV[0] ) / ( V[0]*V[0] );
  }
  double computeDerivativeOfMachNumber(double *V, double *dV, double dMach) const {
    return 1/(2.0*sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)))) * ( ( (2.0*(V[1]*dV[1] + V[2]*dV[2] + V[3]*dV[3]) * V[0] + (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * dV[0]) * (V[4]+Pstiff) - (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] * (dV[4] + dPstiff*dMach) ) / ( (V[4]+Pstiff) * (V[4]+Pstiff) ) );
  }
  double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach) const {
    return 1.0/( 2.0*sqrt(gam * (V[4]+Pstiff) / V[0]) ) * gam * ( (dV[4]+dPstiff*dMach) * V[0] - (V[4]+Pstiff) * dV[0] ) / ( V[0] * V[0] );
  }
  double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dMach) const {
    double mach = computeMachNumber(V);
    double dmach = computeDerivativeOfMachNumber(V, dV, dMach);
    double machr2 = machr*machr;
    double dmachr2 = 2.0*machr*dmachr;
    double popr = V[4]*gam*machr2;
    double dpopr = dV[4]*gam*machr2 + V[4]*gam*dmachr2;
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    double dopmach = gam1*mach*dmach;
    return dV[4]*pow(opmach, gam*invgam1) + V[4]*gam*invgam1*pow(opmach, (gam*invgam1-1))*dopmach;
  }
  void rstVar(IoData &iod) {
    dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));
    rV(iod);
  }

  //----- Equation of State Parameters -----//
  double getGamma()                           const {return gam;}
  double getGamma1()                          const {return gam1;}
  double getPressureConstant()                const {return Pstiff;}
  double getPressure(double *V)               const {return V[4];}
  double getDerivativeOfPressureConstant()    const {return dPstiff;}

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnSGEuler::VarFcnSGEuler(FluidModelData &data) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::GAS){
    fprintf(stderr, "*** Error: FluidModelData is not of type GAS\n");
    exit(1);
  }

  if (data.gasModel.type == GasModelData::IDEAL)
    type = PERFECTGAS;
  else if(data.gasModel.type == GasModelData::STIFFENED)
    type = STIFFENEDGAS;
  else
    fprintf(stdout, "*** Error: VarFcnSGEuler::type is undefined since data.gasModel.type = %d\n", data.gasModel.type);

  gam = data.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = data.gasModel.pressureConstant;
  dPstiff = 0.0;//iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));

  pname = new const char*[5];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::conservativeToPrimitive(double *U, double *V)
{

  V[0] = U[0];

  double invRho = 1.0 / U[0];

  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];

  V[4] = (gam-1.0) * (U[4] - 0.5 * U[0] * vel2) - gam*Pstiff;

}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV)
{

  dV[0] = dU[0];

  double invRho = 1.0 / V[0];

  dV[1] = ( dU[1]  - dV[0] * V[1] ) * invRho;
  dV[2] = ( dU[2]  - dV[0] * V[2] ) * invRho;
  dV[3] = ( dU[3]  - dV[0] * V[3] ) * invRho;

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];

  double dvel2 = 2.0 * V[1] * dV[1] + 2.0 * V[2] * dV[2] + 2.0 * V[3] * dV[3];

  dV[4] = (gam-1.0) * (dU[4] - 0.5 * dU[0] * vel2  - 0.5 * U[0] * dvel2) - gam*dPstiff;

}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::primitiveToConservative(double *V, double *U)
{

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];

  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = (V[4]+gam*Pstiff) * invgam1 + 0.5 * V[0] * vel2;

}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU)
{

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
  double dvel2 = 2.0 * V[1] * dV[1] + 2.0 * V[2] * dV[2] + 2.0 * V[3] * dV[3];

  dU[0] = dV[0];
  dU[1] = dV[0] * V[1] + V[0] * dV[1];
  dU[2] = dV[0] * V[2] + V[0] * dV[2];
  dU[3] = dV[0] * V[3] + V[0] * dV[3];
  dU[4] = (dV[4]+gam*dPstiff) * invgam1 + 0.5 * dV[0] * vel2 + 0.5 * V[0] * dvel2;

}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::extrapolatePrimitive(double un, double c, double *Vb,
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
inline
void VarFcnSGEuler::extrapolateCharacteristic(double n[3], double un, double c, 
                                              double *Vb, double *dV)
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
inline
void VarFcnSGEuler::primitiveToCharacteristicVariations(double n[3], double *V, 
                                                        double *dV, double *dW)
{
  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double ooc2 = V[0]/(gam*(V[4]+Pstiff));
  double oorhoc = sqrt(ooc2)/V[0];
  double coeff1 = dV[0] - ooc2*dV[4];

  dW[0] = n[0]*coeff1 + n[2]*dV[2] - n[1]*dV[3];
  dW[1] = n[1]*coeff1 + n[0]*dV[3] - n[2]*dV[1];
  dW[2] = n[2]*coeff1 + n[1]*dV[1] - n[0]*dV[2];
  dW[3] = dVn + oorhoc*dV[4];
  dW[4] =-dVn + oorhoc*dV[4];
}
//-----------------------------------------------------------------------------
inline
void VarFcnSGEuler::characteristicToPrimitiveVariations(double n[3], double *V, 
                                                        double *dW, double *dV)
{
  double sum = dW[3]+dW[4];
  double diff = dW[3]-dW[4];
  double c = sqrt(gam*(V[4]+Pstiff)/V[0]);

  dV[0] = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2] + 0.5*V[0]*sum/c;
  dV[1] = n[1]*dV[2]-n[2]*dV[1] + 0.5*n[0]*diff;
  dV[2] = n[2]*dV[0]-n[0]*dV[2] + 0.5*n[1]*diff;
  dV[3] = n[0]*dV[1]-n[1]*dV[0] + 0.5*n[2]*diff;
  dV[4] = 0.5*V[0]*c*sum;

}
//------------------------------------------------------------------------------
inline
int VarFcnSGEuler::verification(int glob, double *U, double *V)
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
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::computedVdU(double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  dVdU[20] = gam1 * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dVdU[21] = -gam1 * V[1];
  dVdU[22] = -gam1 * V[2];
  dVdU[23] = -gam1 * V[3];
  dVdU[24] = gam1;
}
//------------------------------------------------------------------------------
inline
void VarFcnSGEuler::computedUdV(double *V, double *dUdV)
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
  dUdV[24] = invgam1;
}
//------------------------------------------------------------------------------

#endif