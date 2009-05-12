#ifndef _LOCAL_RIEMANN_H
#define _LOCAL_RIEMANN_H

#include <LinkF77.h>

class VarFcn;

extern "C" {
  void F77NAME(eriemanngw) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double &,
                            const double &, const double&, const double&,
                            const double &, const double&);
  void F77NAME(eriemanngg) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double &,
                            const double &, const double&, const double&,
                            const double &, const double&);
  void F77NAME(eriemannww) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&);
};
//------------------------------------------------------------------------------

class LocalRiemann {
//this is a virtual class.
//only subclasses will be created:
//    -LocalRiemannGfmpGasGas
//    -LocalRiemannGfmpTaitTait
//    -LocalRiemannGfmparGasGas
//    -LocalRiemannGfmparGasTait
//    -LocalRiemannGfmparTaitTait

public:
  LocalRiemann() {}
  ~LocalRiemann() {}

  // multiphase Riemann problem
  virtual void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj, int it){}

  virtual void eriemann(double rhol, double ul, double pl, 
                        double rhor, double ur, double pr, 
                        double &pi, double &ui,  
                        double &rhoil, double &rhoir,
                        VarFcn *vf){}
  // FSI-type Riemann problem
  virtual void computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei, 
                            double &weighti, int it){}
  // FSI-type Riemann problem
  virtual void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei, 
                            double &weighti, int it){}

protected:
  virtual void solve2x2System(double *mat, double *rhs, double *res);

  void riemannInvariant(VarFcn *vf, double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp);
  void shockJWL(VarFcn *vf, double phi, double omega,
                double omp1oom, double frho, double frhoi, 
                double frhopi,
                double v, double u, double p, double vi,
                double &ui, double &pi, double &dui, double &dpi);

  void shockGAS(double phi, double gamogam1,
                double pref,
                double v, double u, double p, 
                double vi, double &ui, double &pi,
                double &dui, double &dpi);
  void riemannInvariantGAS(VarFcn *vf, double phi,
                   double gam, double gam1, double pref, double c1, 
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp);
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
inline
void LocalRiemann::solve2x2System(double *mat, double *rhs, double *res)
{
  double determinant = mat[0]*mat[3]-mat[1]*mat[2];
  double eps = 1.0e-15;
  double norm = 1.0;
  if(fabs(determinant)>eps*norm){
    res[0] = ( mat[3]*rhs[0]-mat[1]*rhs[1])/determinant;
    res[1] = (-mat[2]*rhs[0]+mat[0]*rhs[1])/determinant;
  }else{
    fprintf(stdout, "zero-determinant\n");
    res[0] = 0.0;
    res[1] = 0.0;
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
inline
void LocalRiemann::riemannInvariant(VarFcn *vf, double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp){
//compute integrals using ODEs. Integrand are not approximated.
//integrate 
  //fprintf(stdout, "riemannInvariant\n");
  int N = 1000;
  double t  = v1;
  double dt = (v-v1)/N;
  double c;
  u = u1; p = p1;
  double V[5] = { 1.0/v1, u1, 0.0, 0.0, p1 };
  bool continueCondition = true;
  //fprintf(stdout, "begin loop\n");
  while(continueCondition){
    c = vf->computeSoundSpeed(V,phi);
    u -= phi*c/t*dt;
    p -= c*c/(t*t)*dt;
    t  += dt;
    continueCondition = (t<v-dt/2.0);
    V[0] = 1.0/t; V[1] = u; V[4] = p;
  }
  //fprintf(stdout, "end loop\n");
  c = vf->computeSoundSpeed(V,phi);
  du = -phi*c/v;
  dp = -c*c/(v*v);
  //fprintf(stdout, "end riemannInvariant\n");
}
inline
void LocalRiemann::shockJWL(VarFcn *vf, double phi, double omega,
                   double omp1oom, double frho, double frhoi, 
                   double frhopi,
                   double v, double u, double p, double vi,
                   double &ui, double &pi, double &dui, double &dpi){
//phi = -1 => left
//phi = +1 => right
  double den=omp1oom*vi-0.5*(vi+v);

  pi = (omp1oom*v*p - 0.5*(vi+v)*p
     + (frhoi*vi-frho*v)/omega
       )/den;

  ui = u + phi * sqrt(-(pi-p)*(vi - v));

  dpi = ((frhoi - frhopi/vi)/omega
      - 0.5*p - (omp1oom - 0.5)*pi)/den;

  dui = -0.5*((vi-v)*dpi+pi-p)/(ui-u);

}
//---------------------------------------------------------------------------
inline
void LocalRiemann::riemannInvariantGAS(VarFcn *vf, double phi,
                   double gam, double gam1, double pref, double c1,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp){
  p  = (p1+pref)*pow(v1/v,gam)-pref;
  double V[5] = {1.0/v, u, 0.0, 0.0, p};
  double c = vf->computeSoundSpeed(V,phi);
  u  = u1 - phi*2.0/gam1*(c1 - c);

  dp = -c*c/(v*v);
  du = -phi*c/v;
}
inline
void LocalRiemann::shockGAS(double phi, double gamogam1,
                   double pref,
                   double v, double u, double p, 
                   double vi, double &ui, double &pi,
                   double &dui, double &dpi){
  double den=gamogam1*vi-0.5*(vi+v);

  pi = (gamogam1*v - 0.5*(vi+v))*(p+pref)/den - pref;

  ui = u + phi * sqrt(-(pi-p)*(vi - v));

  dpi = ( -0.5*(p+pref) - (gamogam1 - 0.5)*(pi+pref))/den;

  dui = -0.5*((vi-v)*dpi+pi-p)/(ui-u);
}

#endif

