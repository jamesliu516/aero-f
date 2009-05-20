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
//    -LocalRiemannGfmparJWLJWL
//    -LocalRiemannGfmparJWLGas

public:
  LocalRiemann() {}
  ~LocalRiemann() {}

  virtual void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj, 
                            double dx[3], int it){}

  virtual void eriemann(double rhol, double ul, double pl, 
                        double rhor, double ur, double pr, 
                        double &pi, double &ui,  
                        double &rhoil, double &rhoir,
                        VarFcn *vf){}

protected:
  virtual void solve2x2System(double *mat, double *rhs, double *res);

  void riemannInvariant2(VarFcn *vf, double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp);
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

  void updatePhaseChangingNodeValues(double * const dx, 
                                     double * const Wi, double * const Wj,
                                     double &weighti, double *rupdatei, 
                                     double &weightj, double *rupdatej);
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
inline
void LocalRiemann::updatePhaseChangingNodeValues(
                      double * const dx, double * const Wi, double * const Wj,
                      double &weighti, double *rupdatei, 
                      double &weightj, double *rupdatej)
{
// In multiphase flow:
// From one iteration to the next, some nodes change phases.
// The state values at these nodes are not relevant to the thermodynamics
//     of the fluid they belong to at the end of the iteration and thus
//     must somehow be replaced or "updated".
// Appropriate state values are given by the interfacial states of the
//     solution of the two-phase Riemann problem.
// In three dimensions, there are several interfacial states to consider.
// The present routine uses all interfacial states that are upwind of
//    the node that may need to be "updated" and weighs them according
//    to their upwind position.

  int dim = 5;

  double temp = 0.0;
  double normdx2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
  double normWi2 = Wi[1]*Wi[1]+Wi[2]*Wi[2]+Wi[3]*Wi[3];
  double normWj2 = Wj[1]*Wj[1]+Wj[2]*Wj[2]+Wj[3]*Wj[3];

  if(normdx2 > 0.0 && normWj2 > 0.0)
    temp = -(Wj[1]*dx[0]+Wj[2]*dx[1]+Wj[3]*dx[2])/sqrt(normdx2*normWj2);
    if (temp > 0.0){
      weighti += temp;
    for (int k=0; k<dim; k++)
      rupdatei[k] += temp*Wj[k];
  }
  temp = 0.0;
  if(normdx2 > 0.0 && normWi2 > 0.0)
    temp = (Wi[1]*dx[0]+Wi[2]*dx[1]+Wi[3]*dx[2])/sqrt(normdx2*normWi2);
  if(temp > 0.0){ // for update of node j
    weightj += temp;
    for (int k=0; k<dim; k++)
      rupdatej[k] += temp*Wi[k];
  }

}
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
    fprintf(stdout, "zero-determinant (mat = [%e , %e ; %e , %e] = %e)\n", mat[0],mat[3],mat[1],mat[2],determinant);
    res[0] = 0.0;
    res[1] = 0.0;
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline
void LocalRiemann::riemannInvariant2(VarFcn *vf, double phi,
                   double v1, double u1, double p1,
                   double v,  double &u, double &p,
                   double &du, double &dp){
//compute integrals using ODEs. Integrand are not approximated.
//integrate using 2nd order integration
  int N  = 500;
  double dt = (v-v1)/N;
  double t = v1; u = u1; p = p1;
  double V[5] = { 1.0/v1, u1, 0.0, 0.0, p1 };
  double c = vf->computeSoundSpeed(V,phi);
  bool continueCondition = true;

  while (continueCondition) {
    V[1] = u - phi*c/t*dt/2;
    V[4] = p - c*c/(t*t)*dt/2;
    t  += dt/2;
    V[0] = 1.0/t;
    if(vf->checkPressure(V,phi) < 0.0) break;
    c = vf->computeSoundSpeed(V,phi);

    u -= phi*c/t*dt;
    p -= c*c/(t*t)*dt;
    t  += dt/2;

    V[0] = 1.0/t; V[1] = u; V[4] = p;
    if(vf->checkPressure(V,phi) < 0.0) break;
    c = vf->computeSoundSpeed(V,phi);

    continueCondition = (t<v-dt/4.0);

  }
  du = -phi*c/v;
  dp = -c*c/(v*v);
  if(vf->checkPressure(V,phi)<0.0){
    u=0.0; p=0.0; du=0.0; dp=0.0;
  }

}

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

