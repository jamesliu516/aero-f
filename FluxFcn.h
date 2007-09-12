#ifndef _FLUX_FCN_H_
#define _FLUX_FCN_H_

#include <VarFcn.h>
#include <LinkF77.h>
                                                                                                 
#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif
/*
extern "C" {
                                                                                                 
   void F77NAME(eriemanngw) (const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&);
   void F77NAME(eriemanngg) (const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&);
   void F77NAME(eriemannww) (const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&, const double&,
														 const double&);
};*/
//----------------------------------------------------------------------------------------
//CHANGES_FOR_WATER
// the class FluxFcn and its subclasses have been modified because the computation of
// the fluxes and their jacobians is different.
//----------------------------------------------------------------------------------------

class FluxFcn {

public:

  enum Type {CONSERVATIVE = 0, PRIMITIVE = 1} type;

protected:
  VarFcn *vf;

public:
  FluxFcn(VarFcn *varFcn, Type tp);
  ~FluxFcn() { if (vf) delete vf; }

  virtual void compute(double, double, double *, double, double *, double *, double *, int = 1) = 0;
  virtual void computeJacobian(double, double *, double, double *, double *, double *, int = 1) {}
  virtual void computeJacobians(double, double *, double, double *, double *, double *, double *, int = 1) {}

  VarFcn *getVarFcn() { return vf; }
/*  void compute_riemann(double, double, double, double, double, double, 
                       double &, double &, double &, double &,
											 double, double, double, double);
  void compute_riemann(double, double, double, double, double, double, 
                       double &, double &, double &, double &,
											 double, double, double, double, int);
	void compute_riemann(double, double, double, double, double, double,
											 double &, double &, double &, double &,
											 double, double, double, double, double, double);
	*/										 
};


//----------------------------------------------------------------------------------------

inline
FluxFcn::FluxFcn(VarFcn *varFcn,Type tp) : vf(varFcn) {
 
  type = tp;
  
}
/*
//------------------------------------------------------------------------------
inline 
void FluxFcn::compute_riemann(double DL, double UL, double PL,
                              double DR, double UR, double PR,
                              double &PI, double &UI, double &RIL, double &RIR,
                              double alpha, double beta, double pref, double gam)
{
   F77NAME(eriemanngw) (DL, UL, PL, DR, UR, PR, PI, UI, RIL, RIR,alpha,beta,pref,gam);
}

//------------------------------------------------------------------------------
inline 
void FluxFcn::compute_riemann(double DL, double UL, double PL,
                              double DR, double UR, double PR,
                              double &PI, double &UI, double &RIL, double &RIR,
			      double gamL, double prefL, double gamR, double prefR, int i)
{
   F77NAME(eriemanngg) (DL, UL, PL, DR, UR, PR, PI, UI, RIL, RIR,gamL,prefL,gamR,prefR);
}

//------------------------------------------------------------------------------
inline 
void FluxFcn::compute_riemann(double DL, double UL, double PL,
                              double DR, double UR, double PR,
                              double &PI, double &UI, double &RIL, double &RIR,
														  double alphal, double betal, double prefl,
											  			double alphar, double betar, double prefr)
{
   F77NAME(eriemannww) (DL, UL, PL, DR, UR, PR, PI, UI, RIL, RIR,
                        alphal,betal,prefl,alphar,betar,prefr);
}
*/
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

template<int dim>
class FluxFcnFD : public FluxFcn {

public:
  FluxFcnFD(VarFcn *vf,Type tp) : FluxFcn(vf,tp) {}
  ~FluxFcnFD() {}

  virtual void compute(double, double, double *, double, double *, double *, double *, int) = 0;
  void computeJacobian(double, double *, double, double *, double *, double *, int);
  void computeJacobians(double, double *, double, double *, double *, double *, double *, int);
  
};

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobian(double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *jacL, int flag)
{

	double dflag = double(flag);
  const double eps0 = 1.e-6;

  double Veps[dim], flux[dim], fluxeps[dim], dfdVL[dim*dim];

  compute(1.0, irey, normal, normalVel, VL, VR, flux, flag); //1.0 for length

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

    compute(1.0, irey, normal, normalVel, Veps, VR, fluxeps, flag);

    for (int j=0; j<dim; ++j) 
      dfdVL[dim*j + k] = (fluxeps[j]-flux[j]) * inveps;

  }

  if (type == CONSERVATIVE){
    vf->postMultiplyBydVdU(VL, dfdVL, jacL, dflag);
    //vf->postMultiplyBydVdU(VL, dfdVL, jacL);
  }
  else{
    for (k=0; k<dim*dim; ++k) 
      jacL[k] = dfdVL[k];
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobians(double irey, double *normal, double normalVel, double *VL,
                                      double *VR, double *jacL, double *jacR, int flag)
{

  double n[3] = {normal[0], normal[1], normal[2]};

  computeJacobian(irey, n, normalVel, VL, VR, jacL, flag);

  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  computeJacobian(irey, n, -normalVel, VR, VL, jacR, flag);
  for (int k=0; k<dim*dim; ++k) jacR[k] = -jacR[k];

}

//------------------------------------------------------------------------------

#endif
