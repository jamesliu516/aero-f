#ifndef _FLUX_FCN_BASE_H_
#define _FLUX_FCN_BASE_H_

#include <VarFcnBase.h>

//----------------------------------------------------------------------------------------

class FluxFcnBase {

public:

  enum Type {CONSERVATIVE = 0, PRIMITIVE = 1} typeJac;

protected:
  VarFcnBase *vf;

public:
  FluxFcnBase(VarFcnBase *varFcn, Type tp);
  virtual ~FluxFcnBase() { vf = 0; }

  virtual void compute(double, double, double *, double, double *, double *, double *, bool) {}
  virtual void computeJacobian(double, double, double *, double, double *, double *, double *, bool) {}
  virtual void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool) {}

  VarFcnBase* getVarFcnBase() const { return vf; }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *vl, double *dvl, double *vr, double *dvr, 
    double dmach, double *f, double *df
    , bool useLimiter = false
  ) 
  {
    std::cout << "\n !!! FluxFcnBase::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *v, double *ub, double *dub, double *f, double *df
  ) 
  {
    std::cout << "\n !!! FluxFcnBase::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};


//----------------------------------------------------------------------------------------

inline
FluxFcnBase::FluxFcnBase(VarFcnBase *varFcn,Type tp) : vf(varFcn) {
 
  typeJac = tp;
  
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

template<int dim>
class FluxFcnFD : public FluxFcnBase {

public:
  FluxFcnFD(VarFcnBase *vf,Type tp) : FluxFcnBase(vf,tp) {}
  ~FluxFcnFD() { vf = 0; }

  virtual void compute(double, double, double *, double, double *, double *, double *, bool){} 
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);
  
};

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobian(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *jacL, bool useLimiter)
{

  const double eps0 = 1.e-6;

  double Veps[dim], flux[dim], fluxeps[dim], dfdVL[dim*dim];

  compute(length, irey, normal, normalVel, VL, VR, flux, useLimiter); 

  int k;
  for (k=0; k<dim; ++k)
    Veps[k] = VL[k];

  for (k=0; k<dim; ++k) {

    double eps;
    if (fabs(Veps[k]) < 1.e-10)
      eps = eps0;
    else
      eps = eps0 * Veps[k];

    double inveps = 1.0 / eps;

    Veps[k] += eps;

    if (k != 0)
      Veps[k-1] = VL[k-1];

    compute(length, irey, normal, normalVel, Veps, VR, fluxeps, useLimiter);

    for (int j=0; j<dim; ++j) 
      dfdVL[dim*j + k] = (fluxeps[j]-flux[j]) * inveps;

  }

  if (typeJac == CONSERVATIVE)
    vf->postMultiplyBydVdU(VL, dfdVL, jacL);
  else{
    for (k=0; k<dim*dim; ++k) 
      jacL[k] = dfdVL[k];
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobians(double length, double irey, double *normal, double normalVel, double *VL,
                                      double *VR, double *jacL, double *jacR, bool useLimiter)
{

  double n[3] = {normal[0], normal[1], normal[2]};

  computeJacobian(length, irey, n, normalVel, VL, VR, jacL, useLimiter);

  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  computeJacobian(length, irey, n, -normalVel, VR, VL, jacR, useLimiter);
  for (int k=0; k<dim*dim; ++k) jacR[k] = -jacR[k];

}

//------------------------------------------------------------------------------

#endif
