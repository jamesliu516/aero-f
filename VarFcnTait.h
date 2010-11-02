#ifndef _VAR_FCN_TAIT_H
#define _VAR_FCN_TAIT_H

#include <VarFcnBase.h>

//--------------------------------------------------------------------------
// This class is the VarFcn class for the Tait EOS.
// Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state, since it is assumed that the Tait EOS that must be used at this 
// point.
//
// lay-out of the base class is:
//  - 1 -  Transformation Operators
//  - 2 -  General Functions
//  - 3 -  Equations of State Parameters
//  - 4 -  EOS related functions
//
//--------------------------------------------------------------------------
//
// EOS: Pressure = p + a*Density^b
//
//--------------------------------------------------------------------------
class VarFcnTait : public VarFcnBase {

private:
  double Cv_;
  double invCv_;
  double a_;
  double b_;
  double p_;

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);
  int verification(int glob, double *U, double *V);

public:
  VarFcnTait(FluidModelData &data);
  ~VarFcnTait() { delete [] pname; }

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);
  
  void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V);
  void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV);
  void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW);
  void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV);

  //----- General Functions -----//
  double getPressure(double *V) const{ return p_ + a_*pow(V[0],b_); }

  void setPressure(const double p, double *V){V[0] = pow( (p-p_)/a_ , 1.0/b_); }
  void setPressure(double *V, double *Vorig) {V[0] = Vorig[0];}

  //checks that the Euler equations are still hyperbolic
  double checkPressure(double *V) const{ return getPressure(V); }

  double computeTemperature(double *V) const{ return V[4]; }
  double computeRhoEnergy(double *V)   const{
    return computeRhoEpsilon(V) + 0.5 * V[0] * getVelocitySquare(V);
  }
  //computes internal energy (=rho*e-0.5*rho*u^2)
  double computeRhoEpsilon(double *V)  const{ return V[0] * Cv_ * V[4]; }
  double computeSoundSpeed(double *V)  const{ 
    double c2 = a_ * b_ * pow(V[0], b_ - 1.0);
    if (c2>0) return sqrt(c2);
    return 0.0;
  }
  double computeSoundSpeed(const double density, const double entropy) const{
    double c2 = a_ * b_ * pow(density, b_ - 1.0);
    if (c2>0) return sqrt(c2);
    return 0.0;
  }
  //double computeEntropy(const double density, const double pressure) const{ }
  //double computeIsentropicPressure(const double entropy, const double density) const{ }

  //----- Equation of State Parameters -----//
  double getCv()         const{ return Cv_; }
  double getAlphaWater() const{ return a_; }
  double getBetaWater()  const{ return b_; }
  double getPrefWater()  const{ return p_; }

  //----- EOS related functions -----//

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnTait::VarFcnTait(FluidModelData &data) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::LIQUID){
    fprintf(stderr, "*** Error: FluidModelData is not of type Tait\n");
    exit(1);
  }

  type = TAIT;
 
  Cv_    = 1.0;
  invCv_ = 1.0/Cv_;
  a_     = data.liquidModel.alpha;
  b_     = data.liquidModel.beta;
  p_     = data.liquidModel.Pref;

  pname = new const char*[5];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature";
}
//------------------------------------------------------------------------------
inline
void VarFcnTait::conservativeToPrimitive(double *U, double *V){

  V[0] = U[0];

  double invRho = 1.0 / U[0];
   
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
      
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
    
  V[4] = invRho * invCv_ * (U[4] - 0.5 * U[0] * vel2);

}
//------------------------------------------------------------------------------
inline
int VarFcnTait::verification(int glob, double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  double locPressure = p_ +a_*pow(V[0],b_);
  if(locPressure<pmin){
    if(verif_clipping)
      fprintf(stdout, "clip pressure[%d] in tait from %e to %e\n", glob, locPressure, pmin);
    V[0] = pow((pmin-p_)/a_, 1.0/b_);
    U[0] = V[0];
    U[1] = V[0]*V[1];
    U[2] = V[0]*V[2];
    U[3] = V[0]*V[3];
    U[4] = V[0]*V[4]*Cv_+0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    return 1;
  }
  return 0;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::primitiveToConservative(double *V, double *U) {
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                                            
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = V[0] * Cv_ * V[4] + 0.5 * V[0] * vel2;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::extrapolatePrimitive(double un, double c, double *Vb,
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
//------------------------------------------------------------------------------
inline
void VarFcnTait::extrapolateCharacteristic(double n[3], double un, double c,
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
  double P = p_ +a_*pow(Vb[0],b_);
  double coeff1 = -P*dV[0]/Vb[0] + Vb[0]*Cv_*dV[4];

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
  dV[4] = ( dW[0]*n[0]+dW[1]*n[1]+dW[2]*n[2] + P*sum*oorho)*oorho/Cv_;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::primitiveToCharacteristicVariations(double n[3], double *V, 
                                                     double *dV, double *dW)
{

  double Ptemp = -(p_+a_*pow(V[0],b_))/pow(V[0],2);
  double c = computeSoundSpeed(V)/V[0];
  double coeff1 = Ptemp*dV[0]+ Cv_*dV[4];
  double dVn = n[0]*dV[1]+n[1]*dV[2]+n[2]*dV[3];

  dW[0] = coeff1*n[0] - n[2]*dV[2] + n[1]*dV[3];
  dW[1] = coeff1*n[1] - n[0]*dV[3] + n[2]*dV[1];
  dW[2] = coeff1*n[2] - n[1]*dV[1] + n[0]*dV[2];
  dW[3] = 0.5*(c*dV[0]+dVn);
  dW[4] = 0.5*(c*dV[0]-dVn);
}
//------------------------------------------------------------------------------
inline
void VarFcnTait::characteristicToPrimitiveVariations(double n[3], double *V, 
                                                     double *dW, double *dV)
{
  double ooCv = 1.0/Cv_;
  double c = computeSoundSpeed(V);
  double Ptemp = (p_+a_*pow(V[0],b_))/(V[0]*Cv_*c);
  double rho = V[0]/c;

  dV[0] = rho*(dW[3]+dW[4]);
  dV[1] = n[2]*dV[1]-n[1]*dV[2] + n[0]*(dW[3]-dW[4]);
  dV[2] = n[0]*dV[2]-n[2]*dV[0] + n[1]*(dW[3]-dW[4]);
  dV[3] = n[1]*dV[0]-n[0]*dV[1] + n[2]*(dW[3]-dW[4]);
  dV[4] = ooCv*(n[0]*dW[0]+n[1]*dW[1]+n[2]*dW[2]) + Ptemp*(dW[3]+dW[4]);

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::computedVdU(double *V, double *dVdU) {

  double invrho = 1.0 / V[0];
  double invrhoCv = invrho*invCv_;
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  dVdU[20] = invrhoCv * ( -Cv_*V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) );
  dVdU[21] = -invrhoCv * V[1];
  dVdU[22] = -invrhoCv * V[2];
  dVdU[23] = -invrhoCv * V[3];
  dVdU[24] = invrhoCv;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::computedUdV(double *V, double *dUdV) {

  dUdV[0]  = 1.0;
  dUdV[5]  = V[1];
  dUdV[6]  = V[0];
  dUdV[10] = V[2];
  dUdV[12] = V[0];
  dUdV[15] = V[3];
  dUdV[18] = V[0];
  dUdV[20] = Cv_ * V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[21] = V[0] * V[1];
  dUdV[22] = V[0] * V[2];
  dUdV[23] = V[0] * V[3];
  dUdV[24] = V[0] * Cv_;

}
//------------------------------------------------------------------------------

#endif
