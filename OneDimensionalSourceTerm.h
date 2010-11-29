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

  virtual void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) = 0;
  virtual void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) = 0;

};

//------------------------------------------------------------------------------

class CartesianOneDSourceTerm : public OneDimensionalSourceTermBase {

public:
  CartesianOneDSourceTerm(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~CartesianOneDSourceTerm() { varFcn = 0; }

  void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) { 
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

  void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    if(Ym[0]>0.0){
      double logTerm = log(Yp[0]/Ym[0]);
      res[0] = V[0]*V[1]*logTerm;
      res[1] = V[0]*V[1]*V[1]*logTerm;
      res[4] = (this->varFcn->computeRhoEnergy(V,fluidId)+this->varFcn->getPressure(V,fluidId))*V[1]*logTerm;
    }else{ //at center, assuming a linear variation of velocity from 0 to XXX while maintaining conservation
      double velocityGrad = 2.0*V[1]; //actually velocity gradient * Yp[0]
      res[0] = V[0]*velocityGrad;
      res[1] = V[0]*velocityGrad*velocityGrad/2.0;
      res[4] = (this->varFcn->computeRhoEnergy(V,fluidId)+this->varFcn->getPressure(V,fluidId))*velocityGrad;
    }
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    if(Ym[0]>0.0) *res = phi[0]*V[1]*log(Yp[0]/Ym[0]);
    else          *res = phi[0]*2.0*V[1];
  }

};

//------------------------------------------------------------------------------

class SphericalOneDSourceTerm : public OneDimensionalSourceTermBase {

public:
  SphericalOneDSourceTerm(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~SphericalOneDSourceTerm() { varFcn = 0; } 

  void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    if(Ym[0]>0.0){
      double logTerm = log(Yp[0]/Ym[0]);
      res[0] = 2.0*V[0]*V[1]*logTerm;
      res[1] = 2.0*V[0]*V[1]*V[1]*logTerm;
      res[4] = 2.0*(this->varFcn->computeRhoEnergy(V,fluidId)+this->varFcn->getPressure(V,fluidId))*V[1]*logTerm;
    }else{ //at center, assuming a linear variation of velocity from 0 to XXX while maintaining conservation
      double velocityGrad = 2.0*V[1]; //actually velocity gradient * Yp[0]
      res[0] = 2.0*V[0]*velocityGrad;
      res[1] = 2.0*V[0]*velocityGrad*velocityGrad/2.0;
      res[4] = 2.0*(this->varFcn->computeRhoEnergy(V,fluidId)+this->varFcn->getPressure(V,fluidId))*velocityGrad;
    }
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    if(Ym[0]>0.0) *res = 2.0*phi[0]*V[1]*log(Yp[0]/Ym[0]);
    else          *res = 2.0*phi[0]*2.0*V[1];
  }

};

//------------------------------------------------------------------------------

class CylindricalOneDSourceTerm2 : public OneDimensionalSourceTermBase {

public:
  CylindricalOneDSourceTerm2(VarFcn *vf) : OneDimensionalSourceTermBase(vf) {}
  ~CylindricalOneDSourceTerm2() { varFcn = 0; } 

  void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
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

  void computeSourceTerm(double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    for(int i=0; i<5; i++) res[i] = 0.0;
    res[1] = -this->varFcn->getPressure(V,fluidId)*(Yp[0]*Yp[0]-Ym[0]*Ym[0]);
  }
  void computeLevelSetSourceTerm(double *phi, double *V, double *Ym, double *Yp, double *res, int fluidId = 0) {
    *res = 0;
  }

};

//------------------------------------------------------------------------------

#endif
