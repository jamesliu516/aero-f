#ifndef _LOCAL_RIEMANN_H
#define _LOCAL_RIEMANN_H

#include <LinkF77.h>
//#include <SparseGrid.h>
#include <IoData.h>
#include <assert.h>

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
/* Depending on the EOS considered, the resolution of the Riemann problem
 * can be done in different manners (cf Toro as well as Quartapelle).
 */

protected:
  VarFcn *vf_;
  double invRhoRef;
  double densityRef;
// to see if an integral computed on the fly or if a tabulated value
// is used for the computation of quantities related to the Riemann invariants.

public:
  LocalRiemann() {densityRef=1.63;}
  LocalRiemann(VarFcn *vf) {vf_ = vf; densityRef=1.63;}
  ~LocalRiemann() {vf_=0;}

  void setReferenceDensity(const double density){ densityRef=density; }

  virtual void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj, 
                            double dx[3], int it){}

  virtual void eriemann(double rhol, double ul, double pl, 
                        double rhor, double ur, double pr, 
                        double &pi, double &ui,  
                        double &rhoil, double &rhoir){}

  void riemannInvariantGeneral1stOrder(double *in, double *res, double *phi);
  void riemannInvariantGeneral2ndOrder(double *in, double *res, double *phi);
protected:
  virtual bool solve2x2System(double *mat, double *rhs, double *res);

  // functions used to determine the Riemann invariants as needed
  // in Quartapelle's algorithm to compute the solution of the Riemann problem

  // valid for General EOS
  virtual void riemannInvariantGeneralTabulation(double *in, double *res);

  // valid for JWL only
  void rarefactionJWL(double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp,
                   MultiFluidData::RiemannComputation type = MultiFluidData::RK2, int flag = 0);
  void rarefactionJWL2ndOrder(double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp);
  void rarefactionJWL1stOrder(double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp);
  void shockJWL(double phi, double omega,
                double omp1oom, double frho, double frhoi, 
                double frhopi,
                double v, double u, double p, double vi,
                double &ui, double &pi, double &dui, double &dpi);

  // valid for Gas (Perfect and Stiffened)
  void shockGAS(double phi, double gamogam1,
                double pref,
                double v, double u, double p, 
                double vi, double &ui, double &pi,
                double &dui, double &dpi);
  void rarefactionGAS(double phi,
                   double gam, double gam1, double pref, double c1, 
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag = 0);


  // function used for the multiphase flow algorithm to update
  // nodes that change phases
  void updatePhaseChangingNodeValues(double * const dx, 
                                     double * const Wi, double * const Wj,
                                     double &weighti, double *rupdatei, 
                                     double &weightj, double *rupdatej);
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline
void LocalRiemann::rarefactionJWL(double phi,
                   double v1, double u1, double p1,
                   double v,  double &u, double &p,
                   double &du, double &dp, 
                   MultiFluidData::RiemannComputation type, int flag){

  //fprintf(stdout, "*** rarefactionJWL %e %e %e %e\n",1.0/v1,u1,p1,1.0/v);
  double entropy = vf_->computeEntropy(1.0/v1,p1, phi);
  double in[2] = {1.0/v1, entropy};
  double res1[1] = {0.0};
  if(type == MultiFluidData::FE){
    riemannInvariantGeneral1stOrder(in,res1,&phi);
  }else if(type == MultiFluidData::RK2){
    riemannInvariantGeneral2ndOrder(in,res1,&phi);
  }else if(type == MultiFluidData::TABULATION2){
    riemannInvariantGeneralTabulation(in,res1);
  }
  in[0] = 1.0/v;
  double res2[1] = {0.0};
  if(type == MultiFluidData::FE)
    riemannInvariantGeneral1stOrder(in,res2,&phi);
  else if(type == MultiFluidData::RK2)
    riemannInvariantGeneral2ndOrder(in,res2,&phi);
  else if(type == MultiFluidData::TABULATION2)
    riemannInvariantGeneralTabulation(in,res2);

  u = u1 - phi*(res2[0]-res1[0]);
  p = vf_->computeIsentropicPressure(entropy, 1.0/v, phi);
  //fprintf(stdout, "*** rarefactionJWL2 %e %e %e %e\n", u1, res2[0], res1[0], p);
  double c = vf_->computeSoundSpeed(1.0/v, entropy, phi);
  du = -phi*c/v;
  dp = -c*c/(v*v);
  if (flag>0 && c<= 0.0) fprintf(stdout, "*** rarefactionJWL returns c=%e, u=%e, p=%e, du=%e, dp=%e\n", c,u,p,du,dp);

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline
void LocalRiemann::riemannInvariantGeneral1stOrder(double *in, double *res,
                                                   double *phi){
// in contains density and pressure
// res is the output result and contains the variation of velocity
  res[0] = 0.0;
  int N  = 5000;
  double density = densityRef; double entropy = in[1];
  //fprintf(stdout, "densityRef = %e and entropy = %e\n", densityRef, entropy);
  double ddensity = (in[0] - densityRef)/N;
  double c = vf_->computeSoundSpeed(density,entropy,*phi);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    res[0] -= c/density*ddensity;
    density  += ddensity;
    if(ddensity>0.0) density = density>in[0] ? in[0] : density;
    else             density = density<in[0] ? in[0] : density;
    //fprintf(stdout, "densityRef = %e and entropy = %e - density = %e\n", densityRef, entropy, density);
    c = vf_->computeSoundSpeed(density,entropy,*phi);
    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/2.0);
    else continueCondition = (density>in[0]-ddensity/2.0);
    it++;
    if(it==N) break;
  }

}

//----------------------------------------------------------------------------

inline
void LocalRiemann::riemannInvariantGeneral2ndOrder(double *in, double *res,
                                                   double *phi){
// in contains density and entropy
// res is the output result and contains the variation of velocity
  res[0] = 0.0;
  int N  = 20000;
  double density = densityRef; double entropy = in[1];
  double ddensity = (in[0] - densityRef)/N;
  double c = vf_->computeSoundSpeed(density,entropy,*phi);
  //fprintf(stdout, "riemannInvariantGeneral2ndOrder-- %e %e %e %e\n", density, entropy, in[0], ddensity);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    density += 0.5*ddensity;
    c = vf_->computeSoundSpeed(density,entropy,*phi);
    res[0] -= c/density*ddensity;
    density += 0.5*ddensity;

    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/2.0);
    else continueCondition = (density>in[0]-ddensity/2.0);
    it++;
    if(it==N) break;

    /* RK-type
    //advance by first half density-step
    res[0] -= c/density*ddensity/2.0;

    density += ddensity;
    if(ddensity>0.0) density = density>in[0] ? in[0] : density;
    else             density = density<in[0] ? in[0] : density;
    //advance by second half density-step
    c = vf_->computeSoundSpeed(density,entropy,*phi);
    res[0] -= c/density*ddensity/2.0;

    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/4.0);
    else continueCondition = (density>in[0]-ddensity/4.0);
    fprintf(stdout, "adim = %e %e %e %e\n", density, ddensity, c, res[0]);
    it++;
    if(it==N) break;*/
  }

}

//----------------------------------------------------------------------------

inline
void LocalRiemann::riemannInvariantGeneralTabulation(double *in, double *res){
  fprintf(stderr, "*** Error: tabulation of the Riemann invariant is only available for Gas-JWL simulation\n");
  exit(1);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// OLD JWL relation routines (for rarefaction only)
inline
void LocalRiemann::rarefactionJWL2ndOrder(double phi,
                   double v1, double u1, double p1,
                   double v,  double &u, double &p,
                   double &du, double &dp){
//compute integrals using ODEs. Integrand are not approximated.
//integrate using 2nd order integration
  int N  = 500;
  double dt = (v-v1)/N;
  double t = v1; u = u1; p = p1;
  double V[5] = { 1.0/v1, u1, 0.0, 0.0, p1 };
  double c = vf_->computeSoundSpeed(V,phi);
  bool continueCondition = true;

  while (continueCondition) {
    V[1] = u - phi*c/t*dt/2;
    V[4] = p - c*c/(t*t)*dt/2;
    t  += dt/2;
    V[0] = 1.0/t;
    if(vf_->checkPressure(V,phi) < 0.0) break;
    c = vf_->computeSoundSpeed(V,phi);

    u -= phi*c/t*dt;
    p -= c*c/(t*t)*dt;
    t  += dt/2;

    V[0] = 1.0/t; V[1] = u; V[4] = p;
    if(vf_->checkPressure(V,phi) < 0.0) break;
    c = vf_->computeSoundSpeed(V,phi);

    continueCondition = (t<v-dt/4.0);

  }
  du = -phi*c/v;
  dp = -c*c/(v*v);
  if(vf_->checkPressure(V,phi)<0.0){
    u=0.0; p=0.0; du=0.0; dp=0.0;
  }

}

inline
void LocalRiemann::rarefactionJWL1stOrder(double phi,
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
    c = vf_->computeSoundSpeed(V,phi);
    u -= phi*c/t*dt;
    p -= c*c/(t*t)*dt;
    t  += dt;
    continueCondition = (t<v-dt/2.0);
    V[0] = 1.0/t; V[1] = u; V[4] = p;
  }
  //fprintf(stdout, "end loop\n");
  c = vf_->computeSoundSpeed(V,phi);
  du = -phi*c/v;
  dp = -c*c/(v*v);
  //fprintf(stdout, "end riemannInvariant\n");
}

inline
void LocalRiemann::shockJWL(double phi, double omega,
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
// GAS relation routines
inline
void LocalRiemann::rarefactionGAS(double phi,
                   double gam, double gam1, double pref, double c1,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag){

  double ppref  = (p1+pref)*pow(v1/v,gam);
  p = ppref - pref;
  double c = sqrt(gam*ppref*v);

  u  = u1 - phi*2.0/gam1*(c1 - c);

  dp = -c*c/(v*v);
  du = -phi*c/v;

  if(flag>0){
    fprintf(stdout, "%e %e %e %e %e %e %e %e %e\n", phi,gam,gam1,pref,c1,1.0/v1,u1,p1,1.0/v);
    fprintf(stdout, "p = %e\n", p);
    fprintf(stdout, "c = %e\n", c);
    fprintf(stdout, "u = %e\n", u);
    fprintf(stdout, "dp = %e\n", dp);
    fprintf(stdout, "du = %e\n", du);
    exit(1);
  }
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
bool LocalRiemann::solve2x2System(double *mat, double *rhs, double *res)
{
  double determinant = mat[0]*mat[3]-mat[1]*mat[2];
  double eps = 1.0e-15;
  double norm = fmax(fabs(mat[0])+fabs(mat[2]), fabs(mat[1])+fabs(mat[3]));
  if(fabs(determinant)>eps*norm){
    res[0] = ( mat[3]*rhs[0]-mat[1]*rhs[1])/determinant;
    res[1] = (-mat[2]*rhs[0]+mat[0]*rhs[1])/determinant;
    return true;
  }else{
    fprintf(stdout, "zero-determinant (mat = [%e , %e ; %e , %e] = %e)\n", mat[0],mat[3],mat[1],mat[2],determinant);
    res[0] = 0.0;
    res[1] = 0.0;
    return false;
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#endif

