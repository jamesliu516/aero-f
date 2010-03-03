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

  virtual void compute(double, double, double *, double, double *, double *, double *) {}
  virtual void computeJacobian(double, double, double *, double, double *, double *, double *) {}
  virtual void computeJacobians(double, double, double *, double, double *, double *, double *, double *) {}

// Included (MB)
  virtual void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *vl, double *dvl, double *vr, double *dvr, double dmach, double *f, double *df) {}
  virtual void computeDerivative(double ire, double dIre, double *n, double *dn, double nv, double dnv, double *v, double *ub, double *dub, double *f, double *df) {}

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

  virtual void compute(double, double, double *, double, double *, double *, double *){} 
  void computeJacobian(double, double, double *, double, double *, double *, double *);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *);
  
};

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobian(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *jacL)
{

  const double eps0 = 1.e-6;

  double Veps[dim], flux[dim], fluxeps[dim], dfdVL[dim*dim];

  compute(length, irey, normal, normalVel, VL, VR, flux); 

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

    compute(length, irey, normal, normalVel, Veps, VR, fluxeps);

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
                                      double *VR, double *jacL, double *jacR)
{

  double n[3] = {normal[0], normal[1], normal[2]};

  computeJacobian(length, irey, n, normalVel, VL, VR, jacL);

  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  computeJacobian(length, irey, n, -normalVel, VR, VL, jacR);
  for (int k=0; k<dim*dim; ++k) jacR[k] = -jacR[k];

}

//------------------------------------------------------------------------------

#endif
