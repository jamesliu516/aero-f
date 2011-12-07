#ifndef _ONE_DIMENSIONAL_SOURCE_TERM_H_
#define _ONE_DIMENSIONAL_SOURCE_TERM_H_

#include <cmath>

#include "VarFcn.h"

//------------------------------------------------------------------------------
// This class computes the source terms encountered in the various forms of
// the one-dimensional Euler equations and LevelSet equation.
// See OneDimensionalSolver.h for more details on the various forms of the
// one-dimensional Euler equations.
// In addition, depending whether the cartesian, cylindrical or spherical
// coordinates are considered, the source terms change.
//------------------------------------------------------------------------------

class OneDimensionalSourceTermBase {

protected:
  VarFcn *varFcn;

public:
  OneDimensionalSourceTermBase(VarFcn *vf) : varFcn(vf) {}
  virtual ~OneDimensionalSourceTermBase() { varFcn = 0; }

  virtual void computeSourceTerm(double *Vm,double* Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) = 0;
  virtual void computeLevelSetSourceTerm(double *phim, double *Vm,double *phip, double *Vp,double x,  double *Ym, double *Yp, double *res, int fluidId = 0) { exit(-1); }

};

//------------------------------------------------------------------------------

class CartesianOneDSourceTerm : public OneDimensionalSourceTermBase {

public:
  CartesianOneDSourceTerm(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~CartesianOneDSourceTerm() { varFcn = 0; }

  void computeSourceTerm(double *Vm,double* Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) { 
    for(int i=0; i<5; i++) res[i] = 0.0;
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    *res = 0.0;
  }

};

//------------------------------------------------------------------------------

class CylindricalOneDSourceTerm : public OneDimensionalSourceTermBase {

public:
  CylindricalOneDSourceTerm(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~CylindricalOneDSourceTerm() { varFcn = 0; }

  void computeSourceTerm(double *Vm,double* Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    if(Ym[0]>0.0){
      double logTerm = log(Yp[0]/Ym[0]);
      res[0] = Vm[0]*Vm[1]*logTerm;
      res[1] = Vm[0]*Vm[1]*Vm[1]*logTerm;
      res[4] = (this->varFcn->computeRhoEnergy(Vm,fluidId)+this->varFcn->getPressure(Vm,fluidId))*Vm[1]*logTerm;
    }else{ //at center, assuming a linear variation of velocity from 0 to XXX while maintaining conservation
      double velocityGrad = 2.0*Vm[1]; //actually velocity gradient * Yp[0]
      res[0] = Vm[0]*velocityGrad;
      res[1] = Vm[0]*velocityGrad*velocityGrad/2.0;
      res[4] = (this->varFcn->computeRhoEnergy(Vm,fluidId)+this->varFcn->getPressure(Vm,fluidId))*velocityGrad;
    }
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    if(Ym[0]>0.0) *res = phi[0]*V[1]*log(Yp[0]/Ym[0]);
    else          *res = phi[0]*2.0*V[1];
  }

};

//------------------------------------------------------------------------------
inline double sph_int_oned2(double q, double r0,double r1,double a, double b,double r) {

  return q/(r0-r1)*(log(r)*(a*r0-b*r1)+r*(b-a));
}
 
inline double sph_int_oned(double q, double r0,double r1,double a, double b, double ra) {
  if (ra > 0.0) {
    return sph_int_oned2(q,r0,r1,a,b,r1)-sph_int_oned2(q,r0,r1,a,b,ra);
  } else {
    return sph_int_oned2(q,r0,r1,a,b,-ra)-sph_int_oned2(q,r0,r1,a,b,r0);
  }
}

class SphericalOneDSourceTerm : public OneDimensionalSourceTermBase {

public:
  SphericalOneDSourceTerm(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~SphericalOneDSourceTerm() { varFcn = 0; } 

  void computeSourceTerm(double *Vm,double* Vp,double x,double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    double A[25];
    if(Ym[0]>0.0 || x > 0.0){
      if (fluidId == fluidId2) {
        res[0] = sph_int_oned(2.0,Ym[0],Yp[0],Vm[0]*Vm[1],Vp[0]*Vp[1],x);
        res[1] = sph_int_oned(2.0,Ym[0],Yp[0],Vm[0]*Vm[1]*Vm[1],Vp[0]*Vp[1]*Vp[1],x);
        res[4] = sph_int_oned(2.0,Ym[0],Yp[0],(this->varFcn->computeRhoEnergy(Vm,fluidId)+this->varFcn->getPressure(Vm,fluidId))*Vm[1],(this->varFcn->computeRhoEnergy(Vp,fluidId2)+this->varFcn->getPressure(Vp,fluidId2))*Vp[1],x);
      }
      else if (x > 0.0) {
        res[0] = sph_int_oned(2.0,Ym[0],Yp[0],Vp[0]*Vp[1],Vp[0]*Vp[1],x);
        res[1] = sph_int_oned(2.0,Ym[0],Yp[0],Vp[0]*Vp[1]*Vp[1],Vp[0]*Vp[1]*Vp[1],x);
        res[4] = sph_int_oned(2.0,Ym[0],Yp[0],(this->varFcn->computeRhoEnergy(Vp,fluidId2)+this->varFcn->getPressure(Vp,fluidId2))*Vp[1],(this->varFcn->computeRhoEnergy(Vp,fluidId2)+this->varFcn->getPressure(Vp,fluidId2))*Vp[1],x);
      }
      else {
        res[0] = sph_int_oned(2.0,Ym[0],Yp[0],Vm[0]*Vm[1],Vm[0]*Vm[1],x);
        res[1] = sph_int_oned(2.0,Ym[0],Yp[0],Vm[0]*Vm[1]*Vm[1],Vm[0]*Vm[1]*Vm[1],x);
        res[4] = sph_int_oned(2.0,Ym[0],Yp[0],(this->varFcn->computeRhoEnergy(Vm,fluidId)+this->varFcn->getPressure(Vm,fluidId))*Vm[1],(this->varFcn->computeRhoEnergy(Vm,fluidId)+this->varFcn->getPressure(Vm,fluidId))*Vm[1],x);      
      } 
    }else{ //at center, assuming a linear variation of velocity from 0 to XXX while maintaining conservation
      res[0] = 2.0*Vp[0]*Vp[1]*(-x)/Yp[0];
      res[1] = 2.0*Vp[0]*Vp[1]*Vp[1]*(-x)/Yp[0];//.0*Vm[0]*velocityGrad*velocityGrad/2.0;
      res[4] = 2.0*(this->varFcn->computeRhoEnergy(Vp,fluidId2)+this->varFcn->getPressure(Vp,fluidId2))*Vp[1]*(-x)/Yp[0];
    }
  }
  void computeLevelSetSourceTerm(double *phim, double *Vm,double *phip, double *Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0) {

    if(Ym[0]>0.0 || x > 0.0){
      *res = sph_int_oned(2.0,Ym[0],Yp[0],phim[0]*Vm[1],phip[0]*Vp[1],x);
    } else {
      *res = 2.0*phip[0]*Vp[1]*(-x)/Yp[0];
    }

    //if(Ym[0]>0.0) *res = 2.0*phi[0]*V[1]*log(Yp[0]/Ym[0]);
    //else          *res = 2.0*phi[0]*2.0*V[1];
  }

};

//------------------------------------------------------------------------------

class CylindricalOneDSourceTerm2 : public OneDimensionalSourceTermBase {

public:
  CylindricalOneDSourceTerm2(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~CylindricalOneDSourceTerm2() { varFcn = 0; } 

  void computeSourceTerm(double *Vm,double* Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    //res[1] = -this->varFcn->getPressure(V,fluidId)*(Yp[0]*Yp[0]-Ym[0]*Ym[0]);
    // TODO
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    *res = 0;
  }

};

//------------------------------------------------------------------------------

class SphericalOneDSourceTerm2 : public OneDimensionalSourceTermBase {

public:
  SphericalOneDSourceTerm2(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~SphericalOneDSourceTerm2() { varFcn = 0; } 

  void computeSourceTerm(double *Vm,double* Vp,double x, double *Ym, double *Yp, double *res, int fluidId = 0,int fluidId2 = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    res[1] = (-this->varFcn->getPressure(Vp,fluidId)*Yp[0]*Yp[0]+this->varFcn->getPressure(Vp,fluidId)*Ym[0]*Ym[0] );
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    *res = 0;
  }

};

//------------------------------------------------------------------------------

#endif
