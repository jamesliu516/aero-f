#ifndef _LOCAL_RIEMANN_H
#define _LOCAL_RIEMANN_H

#include <IoData.h>
#include "VarFcn.h"

//----------------------------------------------------------------------------
// Virtual base class to provide nodal values to compute fluxes at the interface
// between two fluids
class LocalRiemann {

protected:
  VarFcn *vf_;
  int fluid1, fluid2;  //             fluid1 ~ phi>0           fluid2 ~ phi<0
                       //          ~ "outside" ~ "right"    ~ "inside" ~ "left"
                       // GasGas            Gas1                    Gas2
                       // GasTait           Tait                    Gas
                       // TaitTait          Tait1                   Tait2	 
                       // JWLJWL            JWL1                    JWL2
                       // GasJWL            Gas                     JWL
                       // FluidStruct       Gas                     N/A
public:
  LocalRiemann()           { vf_ = 0; fluid1 = 0; fluid2 = 1;}
  LocalRiemann(VarFcn *vf) { vf_ = vf; fluid1 = 0; fluid2 = 1;}
  virtual ~LocalRiemann()  { vf_ = 0; }

  // multiphase Riemann problem
  virtual void storePreviousPrimitive(double *V, int ID, double *storeV, double &weight){}
  virtual void updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight){}
  virtual void computeRiemannSolution(double *Vi, double *Vj,
                            int IDi, int IDj, double *nphi,
                            double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj, 
                            double dx[3], int it) {} 

  virtual void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {} //KW: to be removed
  virtual void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {}
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// class used when the original GFMP is used (no need to solve a two-phase
//                                            Riemann problem)
class LocalRiemannGfmp : public LocalRiemann {

public:
  LocalRiemannGfmp() : LocalRiemann() {}
  LocalRiemannGfmp(VarFcn *vf) : LocalRiemann(vf) {}
  virtual ~LocalRiemannGfmp() { vf_ = 0; }

  void storePreviousPrimitive(double *V, int ID, double *storeV, double &weight){/*nothing to do for GFMP*/}
  void updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight){///*nothing to do for GFMP*/}
    if(ID != IDn) fprintf(stdout, "node changes phase!\n");
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// class used when the exact two-phase Riemann problem is solved
class LocalRiemannGfmpar : public LocalRiemann {

  MultiFluidData::TypePhaseChange phaseChangeType_;

public:
  LocalRiemannGfmpar() : LocalRiemann() {}
  LocalRiemannGfmpar(VarFcn *vf, MultiFluidData::TypePhaseChange phaseChangeType) : LocalRiemann(vf), phaseChangeType_(phaseChangeType) {}
  virtual ~LocalRiemannGfmpar() { vf_ = 0; }

  void storePreviousPrimitive(double *V, int ID, double *storeV, double &weight);
  void updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight){
    if(ID == IDn) return; /* node does not change phase: nothing to do*/
    if(weight<=0.0){ fprintf(stdout, "*** Error: negative weight in LocalRiemannGfmpar::updatePhaseChange\n");
                     exit(1); }
    for(int k=0; k<5; k++) V[k] = newV[k]/weight;
  }

protected:
  // following functions used when Riemann problem formulated as a system of nonlinear equations.
  bool solve2x2System(double *mat, double *rhs, double *res);
  // valid for JWL phase
  void rarefactionJWL(double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp,
                   MultiFluidData::RiemannComputation type = MultiFluidData::RK2, int flag = 0);
  virtual void riemannInvariantGeneralTabulation(double *in, double *res);
  void riemannInvariantGeneral1stOrder(double *in, double *res, double *phi);
  void riemannInvariantGeneral2ndOrder(double *in, double *res, double *phi);

  void shockJWL(double phi, double omega,
                double omp1oom, double frho, double frhoi, 
                double frhopi,
                double v, double u, double p, double vi,
                double &ui, double &pi, double &dui, double &dpi);

  // valid for Gas (Perfect and Stiffened) phase
  void rarefactionGAS(double phi,
                   double gam, double gam1, double pref, double c1, 
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag = 0);
  void shockGAS(double phi, double gamogam1,
                double pref,
                double v, double u, double p, 
                double vi, double &ui, double &pi,
                double &dui, double &dpi);

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
void LocalRiemannGfmpar::storePreviousPrimitive(double *V, int ID, double *storeV, double &weight)
{
  if(phaseChangeType_ == MultiFluidData::RIEMANN_SOLUTION) return;
  else if(phaseChangeType_ == MultiFluidData::EXTRAPOLATION){
  }
  else if(phaseChangeType_ == MultiFluidData::ASIS) return;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmpar::solve2x2System(double *mat, double *rhs, double *res)
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

inline
void LocalRiemannGfmpar::rarefactionJWL(double phi,
                   double v1, double u1, double p1,
                   double v,  double &u, double &p,
                   double &du, double &dp, 
                   MultiFluidData::RiemannComputation type, int flag){

  int myFluidId = (phi>=0) ? fluid1 : fluid2;
  double entropy = vf_->computeEntropy(1.0/v1,p1, myFluidId);
  double *in = 0;
  double res1[1] = {0.0};
  if(type == MultiFluidData::FE){
    in = new double[3];
    in[0] = 1.0/v1; in[1] = entropy; in[2] = 1.0/v;
    riemannInvariantGeneral1stOrder(in,res1,&phi);
  }else if(type == MultiFluidData::RK2){
    in = new double[3];
    in[0] = 1.0/v1; in[1] = entropy; in[2] = 1.0/v;
    riemannInvariantGeneral2ndOrder(in,res1,&phi);
  }else if(type == MultiFluidData::TABULATION2){
    in = new double[2];
    in[0] = 1.0/v1; in[1] = entropy;
    riemannInvariantGeneralTabulation(in,res1);
  }

  double res2[1] = {0.0};
  if(type == MultiFluidData::TABULATION2){
    in[0] = 1.0/v;
    riemannInvariantGeneralTabulation(in,res2);
  }

  u = u1 - phi*(res2[0]-res1[0]);
  p = vf_->computeIsentropicPressure(entropy, 1.0/v, myFluidId);
  double c = vf_->computeSoundSpeed(1.0/v, entropy, myFluidId);
  du = -phi*c/v;
  dp = -c*c/(v*v);
  if (flag>0 && c<= 0.0) fprintf(stdout, "*** rarefactionJWL returns c=%e, u=%e, p=%e, du=%e, dp=%e\n", c,u,p,du,dp);

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::riemannInvariantGeneral1stOrder(double *in, double *res,
                                                   double *phi){
// in contains density and pressure and density
// res is the output result and contains the variation of velocity
// integrates an ODE with first order integration

  int myFluidId = (*phi>=0) ? fluid1 : fluid2;
  res[0] = 0.0;
  int N  = 5000;
  double density = in[2]; double entropy = in[1];
  double ddensity = (in[0] - in[2])/N;
  double c = vf_->computeSoundSpeed(density,entropy,myFluidId);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    res[0] -= c/density*ddensity;
    density  += ddensity;
    if(ddensity>0.0) density = density>in[0] ? in[0] : density;
    else             density = density<in[0] ? in[0] : density;
    c = vf_->computeSoundSpeed(density,entropy,myFluidId);
    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/2.0);
    else continueCondition = (density>in[0]-ddensity/2.0);
    it++;
    if(it==N) break;
  }

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::riemannInvariantGeneral2ndOrder(double *in, double *res,
                                                   double *phi){
// in contains density and entropy and density
// res is the output result and contains the variation of velocity
// integrates an ODE with second order integration (Mid-Point Rule)
  int myFluidId = (*phi>=0) ? fluid1 : fluid2;
  res[0] = 0.0;
  int N  = 2000;
  double density = in[2]; double entropy = in[1];
  double ddensity = (in[0] - in[2])/N;
  double c = vf_->computeSoundSpeed(density,entropy,myFluidId);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    density += 0.5*ddensity;
    c = vf_->computeSoundSpeed(density,entropy,myFluidId);
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
    c = vf_->computeSoundSpeed(density,entropy,myFluidId);
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
void LocalRiemannGfmpar::riemannInvariantGeneralTabulation(double *in, double *res){
  fprintf(stderr, "*** Error: tabulation of the Riemann invariant is only available for Gas-JWL simulation\n");
  exit(1);
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::shockJWL(double phi, double omega,
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
void LocalRiemannGfmpar::rarefactionGAS(double phi,
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

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::shockGAS(double phi, double gamogam1,
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
void LocalRiemannGfmpar::updatePhaseChangingNodeValues( double * const dx, 
                      double * const Wi, double * const Wj,
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
//----------------------------------------------------------------------------
#endif

