#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include "LinkF77.h"
#include <LocalRiemann.h>
#include <VarFcn.h>
#include "IoData.h"
#include "SparseGridCluster.h"
#include <cmath>
#include "ImplicitRiemann.h"
#include "DenseMatrixOps.h"
#include "DebugTools.h"
//----------------------------------------------------------------------------
// First the derived classes of LocalRiemannGfmp (no exact Riemann problem)
// Second the derived classes of LocalRiemannGfmpar (with exact Riemann prob)
// Third the derived classes of LocalRiemannFluidStructure
//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

class LocalRiemannGfmpGasGas : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpGasGas(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpGasGas() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpGasGas!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n");} 

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n");}

private:
  LocalRiemannGfmpGasGas();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasGas::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpTaitTait : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpTaitTait(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpTaitTait() { vf_ = 0; }

void computeRiemannSolution(double *Vi, double *Vj,
                            int IDi, int IDj, double *nphi,
                            double *initWi, double *initWj,
                            double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpTaitTait!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n");}

private:
  LocalRiemannGfmpTaitTait();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  double a1 = vf_->getAlphaWater(this->fluid1);
  double b1 = vf_->getBetaWater(this->fluid1);
  double p1 = vf_->getPrefWater(this->fluid1);
  double a2 = vf_->getAlphaWater(this->fluid2);
  double b2 = vf_->getBetaWater(this->fluid2);
  double p2 = vf_->getPrefWater(this->fluid2);

  double temp = 0;
	
  if(IDi==fluid1){
    temp = p2+a2*pow(Vj[0],b2);
    Wi[0] = pow((temp-p1)/a1,1.0/b1);
    temp = p2+a2*pow(Vj[5],b2);
    Wi[5] = pow((temp-p1)/a1,1.0/b1);
    temp = p1+a1*pow(Vi[0],b1);
    Wj[0] = pow((temp-p2)/a2,1.0/b2);
    temp = p1+a1*pow(Vi[5],b1);
    Wj[5] = pow((temp-p2)/a2,1.0/b2);
  }else{
    temp = p2+a2*pow(Vi[0],b2);
    Wj[0] = pow((temp-p1)/a1,1.0/b1);
    temp = p2+a2*pow(Vi[5],b2);
    Wj[5] = pow((temp-p1)/a1,1.0/b1);
    temp = p1+a1*pow(Vj[0],b1);
    Wi[0] = pow((temp-p2)/a2,1.0/b2);
    temp = p1+a1*pow(Vj[5],b1);
    Wi[5] = pow((temp-p2)/a2,1.0/b2);
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpJWLJWL : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpJWLJWL(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpJWLJWL() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpJWLJWL!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n");}

private:
  LocalRiemannGfmpJWLJWL();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpGasJWL : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpGasJWL(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpGasJWL() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpGasJWL!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n");}

private:
  LocalRiemannGfmpGasJWL();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasJWL::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  //Fedkiw's GFM -- extrapolation of density in the ghost cell
  //fprintf(stderr, "using density to extrapolate in ghost cell\n");
  //Wi[0] = Vi[0]; Wi[5] = Vi[5];
  //Wj[0] = Vj[0]; Wj[5] = Vj[5];

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasGas : public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparGasGas(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparGasGas() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){

  F77NAME(eriemanngg)(  rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,
                        vf_->getGamma(fluid2), vf_->getPressureConstant(fluid2), 
                        vf_->getGamma(fluid1), vf_->getPressureConstant(fluid1)); 
  }

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n");}

private:
  LocalRiemannGfmparGasGas();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasGas::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, 
                double &weighti, double &weightj,
                double dx[3], int it)
{
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf_->getGamma(fluid1);
  double pref1 = vf_->getPressureConstant(fluid1);
  double gam2  = vf_->getGamma(fluid2);
  double pref2 = vf_->getPressureConstant(fluid2);

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}

inline 
void LocalRiemannGfmparGasGas::computeRiemannJacobian(double *Vi, double *Vj,
				                      int IDi, int IDj, double *nphi,
				                      double *Wi, double *Wj,
				                      double dx[3],int it,
						      double* dWidWi,double*  dWidWj,
						      double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf_->getGamma(fluid1);
  double pref1 = vf_->getPressureConstant(fluid1);
  double gam2  = vf_->getGamma(fluid2);
  double pref2 = vf_->getPressureConstant(fluid2);

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    ImplicitRiemann::computeGasGasJacobian(Wi[4], gam2,pref2,P_2,R_2, gam1, pref1, P_1,R_1, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3 );

    /*const double eps = 1.0e-5;
    double r2p = R_2+eps;
    double u2p = U_2+eps;
    double p2p = P_2+eps;
    double r1p = R_1+eps;
    double u1p = U_1+eps;
    double p1p = P_1+eps;
    double vn[4];
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    F77NAME(eriemanngg)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<< "dUj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpj/duj=" << dWjdWj3[1] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<<  "dpi/duj=" << dWidWj3[1] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<<"dpj/dPj=" << dWjdWj3[2] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<< "duj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<<"dPj/dPj=" <<  dWjdWj3[8] << " " << dWidWj3[8] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dui/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[2]-R_i2)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpi/dUi=" <<dWidWi3[1] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUi/dUi=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPi/dUi=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpj/dUi=" << dWjdWi3[1] << " " << (vn[2]-R_i2)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpi/dPi=" <<dWidWi3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUi/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPi/dPi=" << dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[2]-R_i2)/eps << std::endl;*/
  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    ImplicitRiemann::computeGasGasJacobian(Wi[4], gam2, pref2, P_2,R_2,gam1,pref1,P_1,R_1, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );
    /*const double eps = 1.0e-5;
    double r2p = R_2+eps;
    double u2p = U_2+eps;
    double p2p = P_2+eps;
    double r1p = R_1+eps;
    double u1p = U_1+eps;
    double p1p = P_1+eps;
    double vn[4];
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    F77NAME(eriemanngg)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/eps << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/eps << std::endl;
    
    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],gam2,pref2,gam1,pref1);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/eps << std::endl;*/
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasTait: public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparGasTait(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparGasTait() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj, 
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n");}

private:
  LocalRiemannGfmparGasTait();
};

//----------------------------------------------------------------------------

inline void solveSGTait(double Rg,double Ug,double Pg, 
			double Rw,double Uw,double Pw,
			double &Pi,double &Ui, 
			double &Rig, double &Riw,
			double alpha,double beta,
			double pref, double gamma,
			double Pinf) {

  double Q,f,m,n,dQ,df,g,dg,db;
  double ag = sqrt(gamma/Rg*(Pg+Pinf));
  double aw = sqrt(alpha*beta*pow(Rw,beta-1.0)),b;
  Pi = sqrt(Pw*Pg);
  Riw = pow((Pi-pref)/alpha,1.0/beta);
  double dpdrho = alpha*beta*pow(Riw,beta-1.0);
  
  int k = 0;
  while (++k < 100) {

    // Gas relations
    if (Pi > Pg) {
      // Shock
      m = 2.0/((gamma+1.0)*Rg);
      n = (Pg+Pinf)*(gamma-1.0)/(gamma+1.0);
      Q = sqrt( (Pi+Pinf+n)/m );
      dQ = 0.5/(Q*m);
      f = (Pi-Pg)/Q;
      df = 1.0/Q - (Pi-Pg)/(Q*Q)*dQ;
    } else {

      // Rarefaction
      f = 2.0*ag/(gamma-1.0)*( pow((Pi+Pinf)/(Pg+Pinf), (gamma-1.0)/(2.0*gamma) )-1.0 );
      df = 2.0*ag/(gamma-1.0)*((gamma-1.0)/(2.0*gamma)/(Pg+Pinf)*pow((Pi+Pinf)/(Pg+Pinf), (gamma-1.0)/(2.0*gamma)-1.0 ));
    }

    // Liquid relations
    if (Pi > Pw) {

      g = max(1.0e-8,sqrt( alpha*(pow(Riw,beta)-pow(Rw,beta))*(Riw-Rw)/(Riw*Rw) ));
      dg = 0.5/g*(alpha*( beta*pow(Riw,beta-1.0))*(Riw-Rw)/(Riw*Rw) + 
		  (pow(Riw,beta)-pow(Rw,beta))/(Riw*Rw) - 
		  (pow(Riw,beta)-pow(Rw,beta))*(Riw-Rw)/(Riw*Riw*Rw));
    } else {

      g = 2.0*aw/(beta-1.0)*( pow(Riw/Rw, (beta-1.0)*0.5) - 1.0);
      dg = 2.0*aw/(beta-1.0)*(beta-1.0)*0.5/Rw*pow(Riw/Rw, (beta-1.0)*0.5-1.0);
    }

    //b = f+g+Uw-Ug;
    b = f+g-Uw+Ug;
    //std::cout << Pi << " " << b << std::endl;
    db = df+dg/dpdrho;   
    if (fabs(b/db) < 1.0e-6*Pi || fabs(b) < 1.0e-8)
      break;
    Pi -= b/db;
    if (Pi < 1.0e-10)
      Pi = 1.0e-10;

    Riw = pow((Pi-pref)/alpha,1.0/beta);
    dpdrho = alpha*beta*pow(Riw,beta-1.0);
    
  }

  if (k == 100) {
    std::cout << "No convergence in Newton iteration for gas-tait ERS" << std::endl;
    std::cout << "Rg = " << Rg << " Ug = " << Ug << " Pg = " << Pg << std::endl;
    std::cout << "Rw = " << Rw << " Uw = " << Uw << " Pw = " << Pw << std::endl;
    std::cout << "Pi = " << Pi << " Riw = " << Riw << " b = " << b << " db = " << db << std::endl;
    
  }

  Ui = 0.5*(Uw+Ug)+0.5*(f-g);
  //Ui = 0.5*(Uw+Ug)+0.5*(g-f);
  if (Pi > Pg) {
    double h = (gamma-1.0)/(gamma+1.0);
    double j = (Pi+Pinf)/(Pg+Pinf);
    Rig = Rg*(j+h)/(j*h+1.0);
  } else {
    Rig = Rg*pow( (Pi+Pinf)/(Pg+Pinf), 1.0/gamma);
  }

  //std::cout << "Rg = " << Rg << " Ug = " << Ug << " Pg = " << Pg << std::endl;
  //std::cout << "Rw = " << Rw << " Uw = " << Uw << " Pw = " << Pw << std::endl;
  //std::cout << "Pi = " << Pi << " Riw = " << Riw << " Rig = " << Rig << std::endl;
}

inline
void LocalRiemannGfmparGasTait::computeRiemannSolution(double *Vi, double *Vj,
	  int IDi, int IDj, double *nphi,
          double *initWi, double *initWj,
	  double *Wi, double *Wj,
          double *rupdatei, double *rupdatej, double &weighti, double &weightj,
          double dx[3], int it)
{
  int dim = 5;

  double alpha   = vf_->getAlphaWater(fluid2);
  double beta    = vf_->getBetaWater(fluid2);
  double pref    = vf_->getPrefWater(fluid2);
  double gam     = vf_->getGamma(fluid1);
  double Pinf    = vf_->getPressureConstant(fluid1);

  double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (IDi==fluid2) {
    // cell j is gas
    // cell i is tait
    R_g  = Vj[0];     R_w  = Vi[0];
    U_g  = vnj;       U_w  = vni;
    P_g  = vf_->getPressure(Vj, IDj);
    P_w  = vf_->getPressure(Vi, IDi);


    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    solveSGTait(R_g,U_g,P_g, 
		R_w,U_w,P_w,
		P_i,U_i, 
		R_il,R_ir,
		alpha,beta,
		pref,gam,
		Pinf);

    Wi[0]  = R_ir;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];
    T_w  = vf_->computeTemperature(Wi, IDj); //KW: need an explanation
    //T_w  = vf_->computeTemperature(Vi, IDi);
    Wi[4]  = T_w;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_il;                      Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];        Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];        Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];        Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                       Wj[dim+4]  = P_i;

  }else{
    // cell j is tait
    // cell i is gas
    R_g  = Vi[0];     R_w  = Vj[0];
    U_g  = vni;       U_w  = vnj;
    P_g  = vf_->getPressure(Vi, IDi);
    P_w  = vf_->getPressure(Vj, IDj);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    solveSGTait(R_g,U_g,P_g, 
		R_w,U_w,P_w,
		P_i,U_i, 
		R_il,R_ir,
		alpha,beta,
		pref,gam,
		Pinf);

    //std::cout << "P_i = " << P_i << " " << R_il << " " << R_ir << std::endl;

    Wi[0]  = R_il;                      Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];        Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];        Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];        Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                       Wi[dim+4]  = P_i;

    Wj[0]  = R_ir;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
    T_w  = vf_->computeTemperature(Wj, IDi);
    //T_w  = vf_->computeTemperature(Vj, IDj);
    Wj[4]  = T_w;                     Wj[dim+4]  = Wj[4];

  }
  
// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}
//----------------------------------------------------------------------------

inline 
void LocalRiemannGfmparGasTait::computeRiemannJacobian(double *Vi, double *Vj,
						       int IDi, int IDj, double *nphi,
						       double *Wi, double *Wj,
						       double dx[3],int it,
						       double* dWidWi,double*  dWidWj,
						       double* dWjdWi,double*  dWjdWj) {
  int dim = 5;

  double alpha   = vf_->getAlphaWater(fluid2);
  double beta    = vf_->getBetaWater(fluid2);
  double pref    = vf_->getPrefWater(fluid2);
  double gam     = vf_->getGamma(fluid1);
  double Pinf    = vf_->getPressureConstant(fluid1);

  double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir,R_i1,R_i2;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];
  double dTdrho,dTdp;

  if (IDi==fluid2) {
    // cell j is gas
    // cell i is tait
    R_g  = Vj[0];     R_w  = Vi[0];
    U_g  = vnj;       U_w  = vni;
    P_g  = vf_->getPressure(Vj, IDj);
    P_w  = vf_->getPressure(Vi, IDi);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);

    ImplicitRiemann::computeGasTaitJacobian(Wj[4], gam, Pinf, P_g, R_g, alpha,
					    beta, pref, P_w, R_w, 
					    dWjdWj3, dWjdWi3, dWidWi3, dWidWj3);

    /*dTdrho = 1.0/(Wj[0]*(vf_->getGamma(IDj)-1.0));
    dTdp = -vf_->computeTemperature(Wj,IDj)/Wj[0];

    dWidWi3[6] = dTdrho*dWjdWi3[0]+dTdp*dWjdWi[2];
    dWidWi3[7] = dTdrho*dWjdWi3[3]+dTdp*dWjdWi[5];
    dWidWi3[8] = dTdrho*dWjdWi3[6]+dTdp*dWjdWi[8];

    dWidWj3[6] = dTdrho*dWjdWj3[0]+dTdp*dWjdWj[2];
    dWidWj3[7] = dTdrho*dWjdWj3[3]+dTdp*dWjdWj[5];
    dWidWj3[8] = dTdrho*dWjdWj3[6]+dTdp*dWjdWj[8];
    */   
    /* alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);  
    solveSGTait(r2p,U_g,P_g, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(r2p,U_g,P_g,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout << "H1" << std::endl;
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_g) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_g) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_g) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_g) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,u2p,P_g, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,u2p,P_g,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_g) << std::endl;
    std::cout<<  "dUi/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_g) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_g) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_g) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,p2p, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,p2p,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_g) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_g) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_g) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_g) << std::endl;
    
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,P_g, 
		r1p,U_w,p1p,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,r1p,U_w,p1p,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_w) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_w) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_w) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_w) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,P_g, 
		R_w,u1p,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,u1p,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_w) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_w) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_w) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_w) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,p1p,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/eps << std::endl;
    */
  }else{
    // cell j is tait
    // cell i is gas
    R_g  = Vi[0];     R_w  = Vj[0];
    U_g  = vni;       U_w  = vnj;
    P_g  = vf_->getPressure(Vi, IDi);
    P_w  = vf_->getPressure(Vj, IDj);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);

    ImplicitRiemann::computeGasTaitJacobian(Wi[4], gam, 0.0, P_g, R_g, alpha,
					    beta, pref, P_w, R_w, 
					    dWidWi3, dWidWj3, dWjdWj3, dWjdWi3);

    /*dTdrho = 1.0/(Wi[0]*(vf_->getGamma(IDi)-1.0));
    dTdp = -vf_->computeTemperature(Wi,IDi)/Wi[0];

    dWjdWj3[6] = dTdrho*dWidWj3[0]+dTdp*dWidWj[2];
    dWjdWj3[7] = dTdrho*dWidWj3[3]+dTdp*dWidWj[5];
    dWjdWj3[8] = dTdrho*dWidWj3[6]+dTdp*dWidWj[8];

    dWjdWi3[6] = dTdrho*dWidWi3[0]+dTdp*dWidWi[2];
    dWjdWi3[7] = dTdrho*dWidWi3[3]+dTdp*dWidWi[5];
    dWjdWi3[8] = dTdrho*dWidWi3[6]+dTdp*dWidWi[8];
    */
    /*const double eps = 1.0+1.0e-7,eps2=1.0e-6;
    double r2p = R_g*eps+eps2;
    double u2p = U_g*eps+eps2;
    double p2p = P_g*eps+eps2;
    double r1p = R_w*eps+eps2;
    double u1p = U_w*eps+eps2;
    double p1p = pref+alpha*pow(r1p,beta);
    double vn[4];   
    solveSGTait(R_g,U_g,P_g, 
		R_w,U_w,P_w,
		P_i,U_i, 
		R_i2,R_i1,
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_i2,R_i1,alpha,beta,pref,gam);
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);  
    solveSGTait(r2p,U_g,P_g, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(r2p,U_g,P_g,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout << "H1" << std::endl;
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_g) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_g) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_g) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_g) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,u2p,P_g, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,u2p,P_g,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<computeTaitJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
					    double* Vi, double* Vj,
					    double* Wi, double* Wj,
					    double* jacii,double* jacij,
					    double* jacjj,double* jacji)< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_g) << std::endl;
    std::cout<<  "dUi/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_g) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_g) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_g) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,p2p, 
		R_w,U_w,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,p2p,R_w,U_w,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_g) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_g) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_g) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_g) << std::endl;
    
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,P_g, 
		r1p,U_w,p1p,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,r1p,U_w,p1p,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_w) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_w) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_w) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_w) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    solveSGTait(R_g,U_g,P_g, 
		R_w,u1p,P_w,
		vn[0],vn[1], 
		vn[2],vn[3],
		alpha,beta,
		pref,gam,
		0.0);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,u1p,P_w,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_w) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_w) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_w) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_w) << std::endl;
    
    alpha   = vf_->getAlphaWater(fluid1);
    beta    = vf_->getBetaWater(fluid1);
    pref    = vf_->getPrefWater(fluid1);
    gam     = vf_->getGamma(fluid2);
    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,p1p,vn[0],vn[1],vn[2],vn[3],alpha,beta,pref,gam);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/eps << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/eps << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/eps << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/eps << std::endl;
    */
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------

class LocalRiemannGfmparTaitTait: public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparTaitTait(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparTaitTait() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n");}

private:
  LocalRiemannGfmparTaitTait();
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  int dim = 5;

  double alpha1   = vf_->getAlphaWater(fluid1);
  double beta1    = vf_->getBetaWater(fluid1);
  double pref1    = vf_->getPrefWater(fluid1);
  double alpha2   = vf_->getAlphaWater(fluid2);
  double beta2    = vf_->getBetaWater(fluid2);
  double pref2    = vf_->getPrefWater(fluid2);

  //double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_1, P_2, R_1, R_2, U_1, U_2, T_1, T_2;
  double P_i, U_i, R_i1, R_i2;
  
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};
  
  if (IDi==fluid1) {
    // cell j is tait2
    // cell i is tait1
    R_1  = Vi[0];     R_2  = Vj[0];
    U_1  = vni;       U_2  = vnj;
    P_1  = vf_->getPressure(Vi, IDi);
    P_2  = vf_->getPressure(Vj, IDj);
    T_1  = vf_->computeTemperature(Vi, IDi);
    T_2  = vf_->computeTemperature(Vj, IDj);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);
    
    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = T_1;                     Wi[dim+4]  = Wi[4];
    //Wi[4]  = T_2;                     Wi[dim+4]  = Wi[4];
    
    Wj[0]  = R_i2;                      Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];        Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];        Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];        Wj[dim+3]  = Wj[3];
    Wj[4]  = T_2;                       Wj[dim+4]  = Wj[4];
    //Wj[4]  = T_1;                       Wj[dim+4]  = Wj[4];
    
  }else{
    // cell j is tait1
    // cell i is tait2
    R_1  = Vj[0];     R_2  = Vi[0];
    U_1  = vnj;       U_2  = vni;
    P_1  = vf_->getPressure(Vj, IDj);
    P_2  = vf_->getPressure(Vi, IDi);
    T_1  = vf_->computeTemperature(Vj, IDj);
    T_2  = vf_->computeTemperature(Vi, IDi);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);
    
    Wi[0]  = R_i2;                      Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];        Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];        Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];        Wi[dim+3]  = Wi[3];
    Wi[4]  = T_2;                       Wi[dim+4]  = Wi[4];
    //Wi[4]  = T_1;                       Wi[dim+4]  = Wi[4];
    
    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = T_1;                     Wj[dim+4]  = Wj[4];
    //Wj[4]  = T_2;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}

inline 
void LocalRiemannGfmparTaitTait::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  
  double alpha1   = vf_->getAlphaWater(fluid1);
  double beta1    = vf_->getBetaWater(fluid1);
  double pref1    = vf_->getPrefWater(fluid1);
  double alpha2   = vf_->getAlphaWater(fluid2);
  double beta2    = vf_->getBetaWater(fluid2);
  double pref2    = vf_->getPrefWater(fluid2);

  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);

    ImplicitRiemann::computeTaitTaitJacobian(P_i, alpha2,beta2,pref2,P_2,R_2, alpha1, beta1,pref1, P_1,R_1, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3 );
    
    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);//R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/duj=" << dWjdWj3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "duj/duj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpi/duj=" << dWidWj3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<<"dpj/dPj=" << dWjdWj3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "duj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPj/dPj=" <<  dWjdWj3[8] << " " << dWidWj3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "dui/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dUi=" <<dWidWi3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUi/dUi=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPi/dUi=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpj/dUi=" << dWjdWi3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dPi=" <<dWidWi3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUi/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPi/dPi=" << dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
  }else{
    // cell i is fluid2
    // cell j is fluid1
    //std::cout << "ji" << std::endl << std::endl ;
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);

    ImplicitRiemann::computeTaitTaitJacobian(P_i, alpha2,beta2,pref2, P_2,R_2,alpha1, beta1,pref1,P_1,R_1, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );

    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,alpha2,beta2,pref2,alpha1,beta1,pref1);
    F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "dui/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
    
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparJWLJWL : public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparJWLJWL(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparJWLJWL() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){ 
    eriemannjj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir); }

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n");}

private:
  LocalRiemannGfmparJWLJWL();
  void eriemannjj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir);
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it)
{

  bool computeRiemannSolutionJWLJWLimplemented = false;
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1);

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1);

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}
//----------------------------------------------------------------------------

inline 
void LocalRiemannGfmparJWLJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    ImplicitRiemann::computeJwlJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );
    
    /*const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);//R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/duj=" << dWjdWj3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "duj/duj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpi/duj=" << dWidWj3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<<"dpj/dPj=" << dWjdWj3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "duj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPj/dPj=" <<  dWjdWj3[8] << " " << dWidWj3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "dui/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dUi=" <<dWidWi3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUi/dUi=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPi/dUi=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpj/dUi=" << dWjdWi3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dPi=" <<dWidWi3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUi/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPi/dPi=" << dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;*/
  }else{
    // cell i is fluid2
    // cell j is fluid1

    ImplicitRiemann::computeJwlJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );

    /*const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,alpha2,beta2,pref2,alpha1,beta1,pref1);
    F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "dui/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
    */
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparJWLJWL::eriemannjj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir){

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-6;
  int MaxIts = 100;
  int it = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  double omegal = vf_->getOmega(fluid2);
  double omegar = vf_->getOmega(fluid1);
  double omp1ooml = (omegal+1.0)/omegal;
  double omp1oomr = (omegar+1.0)/omegar;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhor = vf_->computeFrho(1.0/vr,fluid2);
  double frhoil = frhol;
  double frhoir = frhor;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid1);
  double frhopir = vf_->computeFrhop(1.0/vr,fluid1);
//check vacuum ?

//start newton iteration loop
  while(!convergence){
    //fprintf(stdout, "%e %e %e %e\n", 1.0/vil, 1.0/vir, ui, pi);

  //compute left  term (shock or rarefaction)
    if( vil < vl){
      //fprintf(stdout, "leftshock\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      //fprintf(stdout, "leftraref\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, 
                     MultiFluidData::RK2, 1);
    }
  //compute right term (shock or rarefaction)
    if( vir < vr){
      //fprintf(stdout, "rightshock\n");
      frhoir  = vf_->computeFrho(1.0/vir,fluid1);
      frhopir = vf_->computeFrhop(1.0/vir,fluid1);
      shockJWL(1.0, omegar, omp1oomr, frhor, frhoir, frhopir, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      //fprintf(stdout, "rightraref\n");
      rarefactionJWL(1.0, vr, ur, pr, vir, uir, pir, duir, dpir, 
                     MultiFluidData::RK2, 1);
    }

    //fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    //fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    //fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    //fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);

  //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    solve2x2System(jacobian,function,increment);

    //fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

  //update values and check bounds
    //fprintf(stdout, "1 -- vil = %e and vir = %e\n", vil, vir);
    //fprintf(stdout, "11-- vil = %e and vir = %e\n", vil-increment[0], vir-increment[1]);
    if(vil - increment[0] < 0.0)
      increment[0] = 0.5*vil;
    if(vir - increment[1] < 0.0)
      increment[1] = 0.5*vir;
    
    vil -= increment[0];
    vir -= increment[1];
    //fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp)
        vil = 0.5*(vil+increment[0]+temp);
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      double temp = omegar*vr/(omegar+2.0);
      if(vir<temp)
        vir = 0.5*(vir+increment[1]+temp);
    }
    //fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop
  if(convergence){
    //fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
    rhoil = 1.0/vil;
    rhoir = 1.0/vir;
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);
  }else{
    fprintf(stdout, "riemann solver did not converged\n");
    exit(1);
  }


}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasJWL : public LocalRiemannGfmpar {

private:
  MultiFluidData::RiemannComputation riemannComputationType_;
  SparseGridCluster *sgCluster_;

public:
  LocalRiemannGfmparGasJWL(VarFcn *vf, int tag1, int tag2, SparseGridCluster *sgCluster, 
                           MultiFluidData::RiemannComputation riemannComputation,
                           MultiFluidData::TypePhaseChange typePhaseChange = MultiFluidData::RIEMANN_SOLUTION) : 
  LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {
    riemannComputationType_ = riemannComputation;
    sgCluster_ = sgCluster;
  }
  ~LocalRiemannGfmparGasJWL(){ vf_ = 0; sgCluster_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){
    eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,-1.0,-1.0); }
  void eriemanngj_wrapper(double *in, double *res, double *para);
  void riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para);
  void riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  void computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n");}

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n");}


private:
  LocalRiemannGfmparGasJWL();

protected:
  void eriemanngj_selector(double rhol, double ul, double pl, 
                           double rhor, double ur, double pr, 
                           double &pi, double &ui, double &rhoil, double &rhoir,
                           double initrhol, double initrhor);
  bool eriemanngj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  double initrhol, double initrhor);

  int riemannInvariantGeneralTabulation(double *in, double *res);
  bool vacuum(const double rhol, const double ul, const double pl,
              const double rhor, const double ur, const double pr,
              double vacuumValues[6]);
  double jwlZeroSoundSpeedJwlDensity(const double density, 
                                     const double pressure);
  double sgZeroDensityPJwlDensity(const double density, 
                                  const double pressure,
                                  const double rho_c0);
  double pressureEqGasDensity(const double gasDensity, 
                              const double gasPressure,
                              const double jwlDensity, 
                              const double jwlPressure,
                              const double interfacialJwlDensity);
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it)
{

  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    //eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,initWj[0],initWi[0]); 
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1); 
    initWi[0] = R_i1;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i2;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    //eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,initWi[0],initWj[0]); 
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1); 
    initWi[0] = R_i2;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i1;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}
inline 
void LocalRiemannGfmparGasJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      rarefactionJWLderivs(-1.0, 1.0/Vj[0], vnj, Vj[4], 1.0/Wj[0] , dVdv, sgCluster_ );
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, &dVdv[0]  );
    } else {
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, NULL);
    }
  
    dWidWi3[1] *= -1.0;
    dWidWi3[3] *= -1.0;
    dWidWi3[5] *= -1.0;
    dWidWi3[7] *= -1.0;
    dWidWj3[1] *= -1.0;
    dWidWj3[3] *= -1.0;
    dWidWj3[5] *= -1.0;
    dWidWj3[7] *= -1.0;
    dWjdWi3[1] *= -1.0;
    dWjdWi3[3] *= -1.0;
    dWjdWi3[5] *= -1.0;
    dWjdWi3[7] *= -1.0;
    dWjdWj3[1] *= -1.0;
    dWjdWj3[3] *= -1.0;
    dWjdWj3[5] *= -1.0;
    dWjdWj3[7] *= -1.0;

    
    /*R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
    //                    alpha2,beta2,pref2,alpha1,beta1,pref1);//R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1);
    eriemanngj_selector(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    //F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/duj=" << dWjdWj3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "duj/duj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpi/duj=" << dWidWj3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<<"dpj/dPj=" << dWjdWj3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "duj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPj/dPj=" <<  dWjdWj3[8] << " " << dWidWj3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "dui/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dUi=" <<dWidWi3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUi/dUi=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPi/dUi=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpj/dUi=" << dWjdWi3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dPi=" <<dWidWi3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUi/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPi/dPi=" << dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;*/
  }else{
    // cell i is fluid2
    // cell j is fluid1

    //std::cout << "ji" << std::endl << std::endl ;

    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      rarefactionJWLderivs(-1.0, 1.0/Vi[0], vni, Vi[4], 1.0/Wi[0] , dVdv, sgCluster_ );
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, &dVdv[0]  );
    }
    else {
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, NULL );
    }
    
    dWidWi3[1] *= -1.0;
    dWidWi3[3] *= -1.0;
    dWidWi3[5] *= -1.0;
    dWidWi3[7] *= -1.0;
    dWidWj3[1] *= -1.0;
    dWidWj3[3] *= -1.0;
    dWidWj3[5] *= -1.0;
    dWidWj3[7] *= -1.0;
    dWjdWi3[1] *= -1.0;
    dWjdWi3[3] *= -1.0;
    dWjdWi3[5] *= -1.0;
    dWjdWi3[7] *= -1.0;
    dWjdWj3[1] *= -1.0;
    dWjdWj3[3] *= -1.0;
    dWjdWj3[5] *= -1.0;
    dWjdWj3[7] *= -1.0;

    /*R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);
    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1);
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    //F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "dui/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
    */
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::eriemanngj_wrapper(
                               double *in, double *res, double *para)
{

  double dummy1, dummy2;
  eriemanngj(in[0], 0.0, in[1], in[2], in[4], in[3], dummy1, dummy2, res[0], res[1], -1.0, -1.0);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::eriemanngj_selector(
                               double rhol, double ul, double pl, 
                               double rhor, double ur, double pr, 
                               double &pi, double &ui,  
                               double &rhoil, double &rhoir,
                               double initrhol, double initrhor)
{

  if(riemannComputationType_==MultiFluidData::TABULATION5){
    double *in  = new double[5]; in[0]=rhol;in[1]=pl;in[2]=rhor;in[3]=pr;in[4]=ur-ul;
    double *res = new double[2]; res[0]=0.0;res[1]=0.0;
    sgCluster_->interpolate(1,&in,&res);
    rhoil = fmax(res[0],0.0); rhoir = fmax(res[1],0.0);
    double d[2]; //dummy variable
    double uir, pir, uil, pil;

    double omegal = vf_->getOmega(fluid2);
    double omp1ooml = (omegal+1.0)/omegal;
    double frhol = vf_->computeFrho(rhol,fluid2);
    double gamr = vf_->getGamma(fluid1);
    double prefr = vf_->getPressureConstant(fluid1);
    double gam1r = vf_->getGamma(fluid1)-1.0;
    double gamogam1r = gamr/gam1r;
    double Vr[5] = { rhor, ur, 0.0, 0.0, pr };
    double cr = vf_->computeSoundSpeed(Vr,fluid1);

    if( rhoil > rhol){
      double frhoil  = vf_->computeFrho(rhoil,fluid2);
      double frhopil = vf_->computeFrhop(rhoil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1]);
    }else
      rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1], riemannComputationType_,1);
    if( rhoir > rhor)
      shockGAS(1.0, gamogam1r, prefr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);
    else
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);

    ui = 0.5*(uil+uir);
    pi = 0.5*(pil+pir);
    // not checking for vacuum!
  }
  else
    eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,initrhol,initrhor);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::eriemanngj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          double initrhol, double initrhor){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0
  int verbose = -1;
  if(verbose>0){
    fprintf(stdout, "---- new Riemann ----\n");
    fprintf(stdout, "initial rhoil, rhoir = %e %e\n", rhol, rhor);
    fprintf(stdout, "initial vil,   vir   = %e %e\n", 1.0/rhol, 1.0/rhor);
  }

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-3;
  int MaxIts = 100;
  int it = 0;
  double relaxationFactorJwl = 0.25;//1.0; //0.85; // must be between 0 and 1
  double relaxationFactorGas = 0.25;//1.0; //0.85; // must be between 0 and 1
  int count = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  vil = initrhol>0.0 ? 1.0/initrhol : vl;
  vir = initrhor>0.0 ? 1.0/initrhor : vr;

  double omegal = vf_->getOmega(fluid2);
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid2);


  double gamr = vf_->getGamma(fluid1);
  double prefr = vf_->getPressureConstant(fluid1);
  double gam1r = vf_->getGamma(fluid1)-1.0;
  double gamogam1r = gamr/gam1r;
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,fluid1);
  double pastiterates[100][2];

//check vacuum
  if(verbose>4) fprintf(stdout, "checking vacuum possibilities\n");
  double vacuumValues[6]; /* rhoil, uil, pil, rhoir, uir, pir */
  vacuumValues[0] = -1.0; // positive if proper vacuum values are computed
  bool checkVacuumValues = false;
  if(checkVacuumValues){
    if(vacuum(rhol,ul,pl,rhor,ur,pr,vacuumValues)){
      if(verbose>-1){
        fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", vacuumValues[0], vacuumValues[3]);
        fprintf(stdout, "uil_vac   = %e and uir_vac   = %e\n", vacuumValues[1], vacuumValues[4]);
        fprintf(stdout, "pil_vac   = %e and pir_vac   = %e\n", vacuumValues[2], vacuumValues[5]);
      }
      rhoil = vacuumValues[0];
      rhoir = vacuumValues[3];
      ui    = 0.5*(vacuumValues[1]+vacuumValues[4]);
      pi    = 0.5*(vacuumValues[2]+vacuumValues[5]);
      return true;
    }
    if(verbose>4) fprintf(stdout, "checking vacuum possibilities -- DONE\n");
  }else{
    if(verbose>4) fprintf(stdout, "no checking of vacuum possibilities\n");
  }

//start newton iteration loop
  while(!convergence){
    if(verbose>0) fprintf(stdout, "\n");
    int status = 1;

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      if(verbose>0) fprintf(stdout, "shockJWL\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      if(verbose>0) fprintf(stdout, "rarefactionJWL\n");
      status = rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,1);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      if(verbose>0) fprintf(stdout, "shockGAS\n");
      shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      if(verbose>0) fprintf(stdout, "rarefactionGAS\n");
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    if(verbose>1){
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
    }

    pastiterates[it][0] = vil;
    pastiterates[it][1] = vir;

    if (!status) {
      fprintf(stdout, "$$$$ status = 0\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      for (int kk = 0; kk <= it; ++kk) {
	fprintf(stdout, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
      }
      exit(1);
    }

    //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    bool solved = solve2x2System(jacobian,function,increment);
    if(!solved){
      fprintf(stdout, "$$$$\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      for (int kk = 0; kk <= it; ++kk) {
	fprintf(stdout, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
      }
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir,1);
    }

    if(verbose>2) fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

    //update values and check bounds

    if(verbose>3) fprintf(stdout, "increment/v = %e %e\n", increment[0]/vil, increment[1]/vir);

    // prevent large increases
    if(-increment[0]>2.0*vil) increment[0] = -2.0*vil;
    if(-increment[1]>2.0*vir) increment[1] = -2.0*vir;
    // prevent large decreases
    if(increment[0]>0.5*vil) increment[0] = 0.5*vil;
    if(increment[1]>0.5*vir) increment[1] = 0.5*vir;

    increment[0] *= relaxationFactorJwl;
    increment[1] *= relaxationFactorGas;
    
    vil -= increment[0];
    vir -= increment[1];
    if(verbose>2) fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vil)/increment[0];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      double temp = (gamr-1.0)/(gamr+1.0)*vr;
      if(vir<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vir)/increment[1];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if(verbose>2) fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    //check - in case of rarefaction at next iteration, 1.0/rhoil may not be above a certain value
    if(vacuumValues[0] > 0.0 && vil>1.0/vacuumValues[0]){ // vacuumValues is negative if it does not contain any proper value(see declaration and definition above)
      vil += increment[0];
      increment[0] = -0.5*(1.0/vacuumValues[0] - vil);
      vil -= increment[0];
    }
    if(verbose>2) fprintf(stdout, "4 -- vil = %e and vir = %e\n", vil, vir);
    if(verbose>0) fprintf(stdout, "rhoil = %e and rhoir = %e\n", 1.0/vil, 1.0/vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop

  if( vil < vl){
    frhoil  = vf_->computeFrho(1.0/vil,fluid2);
    frhopil = vf_->computeFrhop(1.0/vil,fluid2);
    shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
  }else rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,0);
  
  if( vir < vr) shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
  else rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);

  rhoil = 1.0/vil;
  rhoir = 1.0/vir;
  ui    = 0.5*(uil+uir);
  pi    = 0.5*(pil+pir);

  //std::cout << "JWL results: rhoil = " << 1.0/vil << " rhoir = "  << 1.0/vir << " ui = " << ui << " pi = " << pi << std::endl;

  //double Vpp[5] = {rhoil,0,0,0,pi};
  //double cll = vf_->computeSoundSpeed(Vpp,fluid2);

  if(convergence){
    if(verbose>-1) fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
  }else{
    fprintf(stderr, "riemann solver did not converged\n");
    fprintf(stderr, "Warning: solution will be state given by vacuum\n");
    fprintf(stderr, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stderr, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stderr, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stderr, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
    fprintf(stderr, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stderr, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stderr, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stderr, "dpil = %e and dpir = %e\n", dpil, dpir);
    for (int kk = 0; kk <= it; ++kk) {
      fprintf(stderr, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
    }

    fflush(stderr);

    rhoil = vacuumValues[0];
    rhoir = vacuumValues[3];
    uil   = vacuumValues[1];
    uir   = vacuumValues[4];
    pil   = vacuumValues[2];
    pir   = vacuumValues[5];
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);
    if(verbose>-1) fprintf(stdout, "Warning: uil = %e and uir = %e\n", uil, uir);
  }

  if(verbose>-1){
    fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
    fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
  }

  if(convergence) return true;
  else            return false;


}
//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::vacuum(const double rhol, const double ul, const double pl,
                                      const double rhor, const double ur, const double pr,
                                      double vacuumValues[6]){
// notation: JWL on the left and SG on the right
// remember vacuum can occur only between two rarefaction waves, thus decrease of densities

// 1st step: find JWL-density for which there is loss of positivity of c^2 in JWL gas
  double min1 = jwlZeroSoundSpeedJwlDensity(rhol,pl); // returns -1 if none found

// 2nd step: find JWL-density for which the SG-density would become zero when
//           expressing equality of pressures on both sides of interface
//           equivalent to JWL-density for JWL-pressure is below the
//           lowest value of the SG-pressure
  double min2 = sgZeroDensityPJwlDensity(rhol,pl,min1);    // returns -1 if none found

// 3rd step: find max3, JWL-density for which the SG-density would become zero when
//           expressing equality of velocities on both sides of interface
// WARNING: not done because it would be too costly. Way things are computed
//          defines a zero sg-density for JWL-density above max3

  //fprintf(stdout, "min1 = %e - min2 = %e\n", min1, min2);

// compute SG-density corresponding to JWL-density-bound = max(0,min1,min2)
  double rhoil_vac, rhoir_vac;
  if(min1 < 0 && min2 < 0){
    rhoil_vac = 1.0e-14; // 0.0
    rhoir_vac = pressureEqGasDensity(rhor,pr,rhol,pl,0.0);
  }else if(min1 > 0 && min2 < 0){
    rhoil_vac = min1;
    rhoir_vac = pressureEqGasDensity(rhor,pr,rhol,pl,min1);
  }else if(min1 < 0 && min2 > 0){
    rhoil_vac = min2;
    rhoir_vac = 1.0e-14; // 0.0
  }else{
    rhoil_vac = min1>min2 ? min1 : min2;
    rhoir_vac = min1>min2 ? pressureEqGasDensity(rhor,pr,rhol,pl,rhoil_vac) : 1.0e-14;
  }
  //fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", rhoil_vac, rhoir_vac);

  double uil,pil,duil,dpil,uir,pir,duir,dpir;
  rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil_vac, uil, pil, duil, dpil);
  double Vr[5] = { rhor, ur, 0.0, 0.0, pr};
  double cr = vf_->computeSoundSpeed(Vr,fluid1);
  rarefactionGAS(1.0, vf_->getGamma(fluid1), vf_->getGamma(fluid1)-1.0, vf_->getPressureConstant(fluid1), cr, 1.0/rhor, ur, pr, 1.0/rhoir_vac, uir, pir, duir, dpir);
  //fprintf(stdout, "uil_vac = %e and uir_vac = %e\n", uil,uir);

  vacuumValues[0] = rhoil_vac;
  vacuumValues[1] = uil;
  vacuumValues[2] = pil;
  vacuumValues[3] = rhoir_vac;
  vacuumValues[4] = uir;
  vacuumValues[5] = pir;

  if(uil<uir) return true;
  return false; //vacuumValues[0] then contains the lower bound for JWL-density
}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::jwlZeroSoundSpeedJwlDensity(const double density, const double pressure)
{

    bool convergence = false;
    double tol = 1.0e-4;
    int it = 0, maxIt = 30;
    double relaxation = 1.0;

    double entropy = vf_->computeEntropy(density,pressure,fluid2);
    if(entropy >= 0.0) return -1.0;

    double xn = density;
    double xnm1 = xn;
    double dx, fn, dfn;
    double lowerBound = 0.0; // this equation may have more than one root (xn = 0 is always root)
                             // but we want the larger positive root below density
                             // so we use this lowerBound to reduce the search domain.

    while(!convergence){
        fn = (vf_->getOmega(fluid2)+1.0)*entropy*pow(xn,vf_->getOmega(fluid2)) + vf_->computeExponentials2(xn,fluid2);
        if(fn<0.0) lowerBound = xn;
        dfn = entropy*(vf_->getOmega(fluid2)+1.0)*vf_->getOmega(fluid2)*pow(xn,vf_->getOmega(fluid2)-1) + vf_->computeDerivativeOfExponentials2(xn,fluid2);
        //fprintf(stdout, "xn = %e - fn = %e - dfn = %e\n", xn, fn, dfn);
        if(dfn!=0) dx = -relaxation*fn/dfn;
        else{ dx = -0.75*dx; continue; }

        //check lower and upper bounds
        if(xn+dx>density)    dx = 0.25*relaxation*(density-xn);
        if(xn+dx<lowerBound) dx = 0.25*relaxation*(lowerBound-xn);

        if(fn<0 && dfn<0){
          dx = 0.25*(xn-xnm1);
          xn = xnm1;
        }

        if(fabs(2.0*dx/(2*xn+dx))<tol) convergence = true;

        xnm1 = xn;
        xn += dx;
        it++;
        if(it>maxIt) break;
    }

    if(!convergence){
        //if(it>maxIt) fprintf(stdout, "vacuumJwlDensity::maximum number of iteration %d reached\n", it);
        //if(dfn==0)   fprintf(stdout, "vacuumJwlDensity::zero derivative after %d iterations\n", it);
        xn = -1.0; // non-convergence value
    }//else fprintf(stdout, "vacuumJwlDensity::convergence after %d iterations\n", it);
    //fprintf(stdout, "vacuumJwlDensity::fn = %e, dfn = %e, dx = %e, xn = %e\n", fn, dfn, dx, xn);
    return xn;
}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::pressureEqGasDensity(const double gasDensity, const double gasPressure,
                                                      const double jwlDensity, const double jwlPressure,
                                                      const double interfacialJwlDensity)
{
  if(interfacialJwlDensity == 0)
    return gasDensity/pow(1+gasPressure/vf_->getPressureConstant(fluid1),1.0/vf_->getGamma(fluid1));

  double jwlEntropy = vf_->computeEntropy(jwlDensity, jwlPressure, fluid2);
  double gasEntropy = vf_->computeEntropy(gasDensity, gasPressure, fluid1);
  return pow((jwlEntropy*pow(interfacialJwlDensity,vf_->getOmega(fluid2)+1.0)+vf_->computeExponentials(interfacialJwlDensity,fluid2)+vf_->getPressureConstant(fluid1))/gasEntropy,1.0/vf_->getGamma(fluid1));

}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::sgZeroDensityPJwlDensity(const double density, const double pressure,
                                                          const double rho_c0)
{

  int verbose=0;
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - density=%e and pressure=%e and rho_c0=%e\n", density, pressure, rho_c0);
  double entropy = vf_->computeEntropy(density,pressure,fluid2);

  double fn = entropy*pow(rho_c0>0.0 ? rho_c0 : 1.e-14,vf_->getOmega(fluid2)+1.0) + vf_->computeExponentials(rho_c0>0.0 ? rho_c0 : 1.e-14,fluid2) + vf_->getPressureConstant(fluid1);
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - fn(max(rho_c0,0)) = %e\n", fn);
  if(fn>0.0) return -1.0;

  bool convergence = false;
  double tol = 1.0e-4;
  int it =0, maxIt = 30;
  double relaxation = 1.0;


  double xn = density;
  double xnm1 = xn;
  double dx=0, dfn;

  while(!convergence){
    fn = entropy*pow(xn,vf_->getOmega(fluid2)+1.0) + vf_->computeExponentials(xn,fluid2) + vf_->getPressureConstant(fluid1);
    dfn = entropy*(vf_->getOmega(fluid2)+1.0)*pow(xn,vf_->getOmega(fluid2)) + vf_->computeDerivativeOfExponentials(xn,fluid2);
    //fprintf(stdout, "it = %d - xn = %e - fn = %e - dfn = %e\n", it, xn, fn, dfn);
    if(dfn>0) dx = -relaxation*fn/dfn;
    
    if(xn+dx<0)       dx = -0.25*relaxation*xn;
    if(xn+dx>density) dx = 0.25*relaxation*(density-xn);

    if(fn<0 && dfn<0){
      dx = 0.5*(xn-xnm1);
      xn = xnm1;
    }

    if(fabs(2.0*dx/(2*xn+dx))<tol) convergence = true;

    xnm1 = xn;
    xn += dx;
    it++;
    if(it>maxIt) break;
  }

  if(!convergence){
    xn = -1.0;
  }
  return xn;

}

//----------------------------------------------------------------------------
inline
int LocalRiemannGfmparGasJWL::riemannInvariantGeneralTabulation(double *in, 
                                                                 double *res){
  //fprintf(stdout, "in-value = %e %e\n", in[0],in[1]);
  return sgCluster_->interpolate(1,&in,&res);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para){

  fprintf(stdout, "in[0] = %e\n", in[0]);
  fprintf(stdout, "in[1] = %e\n", in[1]);
  fprintf(stdout, "para[1] = %e\n", para[1]);
  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral1stOrder(locin, res, &(para[0]));

}
//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para){

  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral2ndOrder(locin, res, &(para[0]));

}

//---------------------------------------------------------------------------
// Tait-JWL
//---------------------------------------------------------------------------
class LocalRiemannGfmparTaitJWL : public LocalRiemannGfmpar {

private:
  MultiFluidData::RiemannComputation riemannComputationType_;
  SparseGridCluster *sgCluster_;

  double sign;

public:
  LocalRiemannGfmparTaitJWL(VarFcn *vf, int tag1, int tag2, SparseGridCluster *sgCluster, 
                           MultiFluidData::RiemannComputation riemannComputation,
			    MultiFluidData::TypePhaseChange typePhaseChange = MultiFluidData::RIEMANN_SOLUTION, 
			    double sgn = 1.0) : 
    LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) , sign(sgn) {
    riemannComputationType_ = riemannComputation;
    sgCluster_ = sgCluster;
  }
  ~LocalRiemannGfmparTaitJWL(){ vf_ = 0; sgCluster_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){
    eriemanntj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,-1.0,-1.0); }
  void eriemanntj_wrapper(double *in, double *res, double *para);
  void riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para);
  void riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para);

private:
  LocalRiemannGfmparTaitJWL();

protected:
  bool eriemanntj_selector(double rhol, double ul, double pl, 
                           double rhor, double ur, double pr, 
                           double &pi, double &ui, double &rhoil, double &rhoir,
                           double initrhol, double initrhor);
  bool eriemanntj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  double initrhol, double initrhor);

  int riemannInvariantGeneralTabulation(double *in, double *res);
};

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparTaitJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it)
{

  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1); 
    initWi[0] = R_i1;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i2;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4] = P_i;
    Wi[4]  = vf_->computeTemperature(Wi, IDj);         Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1); 
    initWi[0] = R_i2;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i1;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;
    Wj[4]  = vf_->computeTemperature(Wj, IDi);                 Wj[dim+4]  = Wj[4];
  }

// METHOD 2 : combine averaging and direction of flow
  if (it == 1)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

}
inline 
void LocalRiemannGfmparTaitJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      rarefactionJWLderivs(-1.0, 1.0/Vj[0], vnj, Vj[4], 1.0/Wj[0] , dVdv, sgCluster_ );
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, &dVdv[0] );
    } else
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3,NULL);

    dWidWi3[1] *= -1.0*sign;
    dWidWi3[3] *= -1.0*sign;
    dWidWi3[5] *= -1.0*sign;
    dWidWi3[7] *= -1.0*sign;
    dWidWj3[1] *= -1.0*sign;
    dWidWj3[3] *= -1.0*sign;
    dWidWj3[5] *= -1.0*sign;
    dWidWj3[7] *= -1.0*sign;
    dWjdWi3[1] *= -1.0*sign;
    dWjdWi3[3] *= -1.0*sign;
    dWjdWi3[5] *= -1.0*sign;
    dWjdWi3[7] *= -1.0*sign;
    dWjdWj3[1] *= -1.0*sign;
    dWjdWj3[3] *= -1.0*sign;
    dWjdWj3[5] *= -1.0*sign;
    dWjdWj3[7] *= -1.0*sign;

    
    /*R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
    //                    alpha2,beta2,pref2,alpha1,beta1,pref1);//R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1);
    eriemanntj_selector(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    //F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/duj=" << dWjdWj3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "duj/duj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPj/duj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpi/duj=" << dWidWj3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<<"dpj/dPj=" << dWjdWj3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "duj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPj/dPj=" <<  dWjdWj3[8] << " " << dWidWj3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "dui/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dUi=" <<dWidWi3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUi/dUi=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPi/dUi=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpj/dUi=" << dWjdWi3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dPi=" <<dWidWi3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUi/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPi/dPi=" << dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
*/
  }else{
    // cell i is fluid2
    // cell j is fluid1

    //std::cout << "ji" << std::endl << std::endl ;
    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      rarefactionJWLderivs(-1.0, 1.0/Vi[0], vni, Vi[4], 1.0/Wi[0] , dVdv, sgCluster_ );
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, &dVdv[0] );
    } else
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3,NULL);

    dWidWi3[1] *= -1.0*sign;
    dWidWi3[3] *= -1.0*sign;
    dWidWi3[5] *= -1.0*sign;
    dWidWi3[7] *= -1.0*sign;
    dWidWj3[1] *= -1.0*sign;
    dWidWj3[3] *= -1.0*sign;
    dWidWj3[5] *= -1.0*sign;
    dWidWj3[7] *= -1.0*sign;
    dWjdWi3[1] *= -1.0*sign;
    dWjdWi3[3] *= -1.0*sign;
    dWjdWi3[5] *= -1.0*sign;
    dWjdWi3[7] *= -1.0*sign;
    dWjdWj3[1] *= -1.0*sign;
    dWjdWj3[3] *= -1.0*sign;
    dWjdWj3[5] *= -1.0*sign;
    dWjdWj3[7] *= -1.0*sign;

    /*R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);
    const double eps = 1.0+1.0e-7;
    double r2p = R_2*eps+1.0e-8;
    double u2p = U_2*eps+1.0e-8;
    double p2p = P_2*eps+1.0e-8;
    double r1p = R_1*eps+1.0e-8;
    double u1p = U_1*eps+1.0e-8;
    double p1p = P_1*eps+1.0e-8;
    double vn[4];
    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1);
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    //F77NAME(eriemannww)(r2p,U_2,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    std::cout<< "dpi/dpi=" << dWidWi3[0] << " " << (vn[2]-R_i2)/(r2p-R_2) << std::endl;
    std::cout<< "dUi/dpi=" << dWidWi3[3] << " " << dWjdWi3[3] << " " << (vn[1]-U_i)/(r2p-R_2) << std::endl;
    std::cout<< "dPi/dpi=" << dWidWi3[6] << " " << dWjdWi3[6] << " " << (vn[0]-P_i)/(r2p-R_2) << std::endl;
    std::cout<< "dpj/dpi=" << dWjdWi3[0] << " " << (vn[3]-R_i1)/(r2p-R_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,u2p,P_2,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpi/dui=" << dWidWi3[1] << " " << (vn[2]-R_i2)/(u2p-U_2) << std::endl;
    std::cout<<  "dui/dui=" << dWidWi3[4] << " " << dWjdWi3[4] << " " <<(vn[1]-U_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dPi/dui=" << dWidWi3[7] << " " << dWjdWi3[7] << " " <<(vn[0]-P_i)/(u2p-U_2) << std::endl;
    std::cout<<  "dpj/dui=" << dWjdWi3[1] << " " << (vn[3]-R_i1)/(u2p-U_2) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,p2p,R_1,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<<"dpi/dPi=" << dWidWi3[2] << " " << (vn[2]-R_i2)/(p2p-P_2) << std::endl;
    std::cout<< "dui/dPi=" << dWidWi3[5] << " " << dWjdWi3[5] << " " <<(vn[1]-U_i)/(p2p-P_2) << std::endl;
    std::cout<<"dPi/dPi=" <<  dWidWi3[8] << " " << dWjdWi3[8] << " " <<(vn[0]-P_i)/(p2p-P_2) << std::endl;
    std::cout<< "dpj/dPi=" << dWjdWi3[2] << " " << (vn[3]-R_i1)/(p2p-P_2) << std::endl;
    
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,r1p,U_1,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dpj=" << dWjdWj3[0] << " " << (vn[3]-R_i1)/(r1p-R_1) << std::endl;
    std::cout<< "duj/dpj=" << dWjdWj3[3] << " " << dWidWj3[3] << " " <<(vn[1]-U_i)/(r1p-R_1) << std::endl;
    std::cout<< "dPj/dpj=" << dWjdWj3[6] << " " << dWidWj3[6] << " " <<(vn[0]-P_i)/(r1p-R_1) << std::endl;
    std::cout<< "dpi/dpj=" << dWidWj3[0] << " " << (vn[2]-R_i2)/(r1p-R_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,R_1,u1p,P_1,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dUj=" <<dWjdWj3[1] << " " << (vn[3]-R_i1)/(u1p-U_1) << std::endl;
    std::cout<< "dUj/dUj=" << dWjdWj3[4] << " " << dWidWj3[4] << " " << (vn[1]-U_i)/(u1p-U_1) << std::endl;
    std::cout<< "dPj/dUj=" << dWjdWj3[7] << " " << dWidWj3[7] << " " << (vn[0]-P_i)/(u1p-U_1) << std::endl;
    std::cout<< "dpi/dUj=" << dWidWj3[1] << " " << (vn[2]-R_i2)/(u1p-U_1) << std::endl;
    
    //F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],alpha2,beta2,pref2,alpha1,beta1,pref1);
    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,p1p,vn[0],vn[1],vn[2],vn[3],R_2,R_1);
    std::cout<< "dpj/dPj=" <<dWjdWj3[2] << " " << (vn[3]-R_i1)/(p1p-P_1) << std::endl;
    std::cout<< "dUj/dPj=" << dWjdWj3[5] << " " << dWidWj3[5] << " " << (vn[1]-U_i)/(p1p-P_1) << std::endl;
    std::cout<< "dPj/dPj=" << dWjdWj3[8] << " " << dWidWj3[8] << " " << (vn[0]-P_i)/(p1p-P_1) << std::endl;
    std::cout<< "dpi/dPj=" << dWidWj3[2] << " " << (vn[2]-R_i2)/(p1p-P_1) << std::endl;
    */
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::eriemanntj_wrapper(
                               double *in, double *res, double *para)
{

  double dummy1, dummy2;
  eriemanntj(in[0], 0.0, in[1], in[2], in[4], in[3], dummy1, dummy2, res[0], res[1], -1.0, -1.0);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparTaitJWL::eriemanntj_selector(
                               double rhol, double ul, double pl, 
                               double rhor, double ur, double pr, 
                               double &pi, double &ui,  
                               double &rhoil, double &rhoir,
                               double initrhol, double initrhor)
{

  if(riemannComputationType_==MultiFluidData::TABULATION5){
    double *in  = new double[5]; in[0]=rhol;in[1]=pl;in[2]=rhor;in[3]=pr;in[4]=ur-ul;
    double *res = new double[2]; res[0]=0.0;res[1]=0.0;
    sgCluster_->interpolate(1,&in,&res);
    rhoil = fmax(res[0],0.0); rhoir = fmax(res[1],0.0);
    double d[2]; //dummy variable
    double uir, pir, uil, pil;

    double omegal = vf_->getOmega(fluid2);
    double omp1ooml = (omegal+1.0)/omegal;
    double frhol = vf_->computeFrho(rhol,fluid2);
    double alphar = vf_->getAlphaWater(fluid1);
    double betar = vf_->getBetaWater(fluid1);
    double pinfr = vf_->getPrefWater(fluid1);
    double Vr[5] = { rhor, ur, 0.0, 0.0, pr };
    double cr = vf_->computeSoundSpeed(Vr,fluid1);

    if( rhoil > rhol){
      double frhoil  = vf_->computeFrho(rhoil,fluid2);
      double frhopil = vf_->computeFrhop(rhoil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1]);
    }else
      rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1], riemannComputationType_,1);
    if( rhoir > rhor)
      shockTAIT(1.0, alphar,betar,pinfr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1],0);
    else
      rarefactionTAIT(1.0, alphar,betar,pinfr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1],0);

    ui = 0.5*(uil+uir);
    pi = 0.5*(pil+pir);
    // not checking for vacuum!
  }
  else
    eriemanntj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,initrhol,initrhor);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparTaitJWL::eriemanntj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          double initrhol, double initrhor){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0
  int verbose = -1;
  if(verbose>0){
    fprintf(stdout, "---- new Riemann ----\n");
    fprintf(stdout, "initial rhoil, rhoir = %e %e\n", rhol, rhor);
    fprintf(stdout, "initial vil,   vir   = %e %e\n", 1.0/rhol, 1.0/rhor);
  }

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-3;
  int MaxIts = 100;
  int it = 0;
  double relaxationFactorJwl = 1.0; //0.85; // must be between 0 and 1
  double relaxationFactorGas = 1.0; //0.85; // must be between 0 and 1
  int count = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  double pmin = max( vf_->getPmin(fluid1), vf_->getPmin(fluid2) );
  vil = initrhol>0.0 ? 1.0/initrhol : vl;
  vir = initrhor>0.0 ? 1.0/initrhor : vr;

  double omegal = vf_->getOmega(fluid2);
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid2);


  double alphar = vf_->getAlphaWater(fluid1);
  double betar = vf_->getBetaWater(fluid1);
  double pinfr = vf_->getPrefWater(fluid1);
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,fluid1);

//check vacuum
/*  if(verbose>4) fprintf(stdout, "checking vacuum possibilities\n");
  double vacuumValues[6]; // rhoil, uil, pil, rhoir, uir, pir 
  vacuumValues[0] = -1.0; // positive if proper vacuum values are computed
  bool checkVacuumValues = false;
  if(checkVacuumValues){
    if(vacuum(rhol,ul,pl,rhor,ur,pr,vacuumValues)){
      if(verbose>-1){
        fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", vacuumValues[0], vacuumValues[3]);
        fprintf(stdout, "uil_vac   = %e and uir_vac   = %e\n", vacuumValues[1], vacuumValues[4]);
        fprintf(stdout, "pil_vac   = %e and pir_vac   = %e\n", vacuumValues[2], vacuumValues[5]);
      }
      rhoil = vacuumValues[0];
      rhoir = vacuumValues[3];
      ui    = 0.5*(vacuumValues[1]+vacuumValues[4]);
      pi    = 0.5*(vacuumValues[2]+vacuumValues[5]);
      return true;
    }
    if(verbose>4) fprintf(stdout, "checking vacuum possibilities -- DONE\n");
  }else{
    if(verbose>4) fprintf(stdout, "no checking of vacuum possibilities\n");
  }
*/
//start newton iteration loop
  while(!convergence){
    if(verbose>0) fprintf(stdout, "\n");

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      if(verbose>0) fprintf(stdout, "shockJWL\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      if(verbose>0) fprintf(stdout, "rarefactionJWL\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,1);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      if(verbose>0) fprintf(stdout, "shockTAIT\n");
      shockTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
    }
    else{
      if(verbose>0) fprintf(stdout, "rarefactionTAIT\n");
      rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
    }

    if(verbose>1){
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
    }

    //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    bool solved = solve2x2System(jacobian,function,increment);
    if(!solved){
      fprintf(stdout, "$$$$\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,1);
    }

    if(verbose>2) fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

    //update values and check bounds

    if(verbose>3) fprintf(stdout, "increment/v = %e %e\n", increment[0]/vil, increment[1]/vir);

    // prevent large increases
    if(-increment[0]>2.0*vil) increment[0] = -2.0*vil;
    if(-increment[1]>2.0*vir) increment[1] = -2.0*vir;
    // prevent large decreases
    if(increment[0]>0.5*vil) increment[0] = 0.5*vil;
    if(increment[1]>0.5*vir) increment[1] = 0.5*vir;

    increment[0] *= relaxationFactorJwl;
    increment[1] *= relaxationFactorGas;
    
    vil -= increment[0];
    vir -= increment[1];
    if(verbose>2) fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vil)/increment[0];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      /*fprintf(stdout,"??\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      double pit = pinfr+alphar*pow(vir, -betar);
      fprintf(stdout,"pit = %e pmin = %e\n",pit,pmin);
      if (pit < pmin) {

	vir = pow( (pmin-pinfr) / alphar, betar);
      }

      double Pt,Ut,Rigt,Riwt;
      solveSGTait(rhol,ul,pl, 
		  rhor,ur,pr,
		  Pt,Ut, 
		  Rigt, Riwt,
		  alphar,betar,
		  pinfr, omegal+1.0,0.0);
      fprintf(stdout, "SGTAIT: rhoil  = %e and rhoir  = %e\n", Rigt, Riwt);
      fprintf(stdout, "SGTAIT: uil  = %e and uir  = %e\n", Ut, Ut);
      fprintf(stdout, "SGTAIT: pil  = %e and pir  = %e\n", Pt,Pt);*/
      /*double temp = (gamr-1.0)/(gamr+1.0)*vr;
      if(vir<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vir)/increment[1];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
	}*/
      }
    if(verbose>2) fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    //check - in case of rarefaction at next iteration, 1.0/rhoil may not be above a certain value
    /*if(vacuumValues[0] > 0.0 && vil>1.0/vacuumValues[0]){ // vacuumValues is negative if it does not contain any proper value(see declaration and definition above)
      vil += increment[0];
      increment[0] = -0.5*(1.0/vacuumValues[0] - vil);
      vil -= increment[0];
      }*/
    if(verbose>2) fprintf(stdout, "4 -- vil = %e and vir = %e\n", vil, vir);
    if(verbose>0) fprintf(stdout, "rhoil = %e and rhoir = %e\n", 1.0/vil, 1.0/vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop

  if( vil < vl){
    frhoil  = vf_->computeFrho(1.0/vil,fluid2);
    frhopil = vf_->computeFrhop(1.0/vil,fluid2);
    shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
  }else rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,0);
  
  if( vir < vr) shockTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
  else rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);

  rhoil = 1.0/vil;
  rhoir = 1.0/vir;
  ui    = 0.5*(uil+uir);
  pi    = 0.5*(pil+pir);

  if(convergence){
    if(verbose>-1) fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
  }else{
    if(verbose>-1) fprintf(stdout, "riemann solver did not converged\n");
    if(verbose>-1) fprintf(stdout, "Warning: solution will be state given by vacuum\n");
    /*    rhoil = vacuumValues[0];
    rhoir = vacuumValues[3];
    uil   = vacuumValues[1];
    uir   = vacuumValues[4];
    pil   = vacuumValues[2];
    pir   = vacuumValues[5];
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);*/
    if(verbose>-1) fprintf(stdout, "Warning: uil = %e and uir = %e\n", uil, uir);
  }

  //if (vir < vr) {
  /*fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
    fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);

    double Pt,Ut,Rigt,Riwt;
  solveSGTait(rhol,ul,pl, 
	      rhor,ur,pr,
	      Pt,Ut, 
	      Rigt, Riwt,
	      alphar,betar,
	      pinfr, omegal+1.0,0.0);
    fprintf(stdout, "SGTAIT: rhoil  = %e and rhoir  = %e\n", Rigt, Riwt);
    fprintf(stdout, "SGTAIT: uil  = %e and uir  = %e\n", Ut, Ut);
    fprintf(stdout, "SGTAIT: pil  = %e and pir  = %e\n", Pt,Pt);*/
    //}
  
  if(verbose>-1){
    fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
    fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
  }

  if(convergence) return true;
  else            return false;


}
//----------------------------------------------------------------------------
inline
int LocalRiemannGfmparTaitJWL::riemannInvariantGeneralTabulation(double *in, 
                                                                 double *res){
  //fprintf(stdout, "in-value = %e %e\n", in[0],in[1]);
  return sgCluster_->interpolate(1,&in,&res);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para){

  fprintf(stdout, "in[0] = %e\n", in[0]);
  fprintf(stdout, "in[1] = %e\n", in[1]);
  fprintf(stdout, "para[1] = %e\n", para[1]);
  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral1stOrder(locin, res, &(para[0]));

}
//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para){

  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral2ndOrder(locin, res, &(para[0]));

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Adam 2010.08.17
// This Local Riemann is now templated by dim. It can now handle Turbulent simulations
template<int dim>
class LocalRiemannFluidStructure : public LocalRiemann {

public:
  LocalRiemannFluidStructure() : LocalRiemann() {fluid1 = fluid2 = 0;}
  LocalRiemannFluidStructure(VarFcn *vf) : LocalRiemann(vf,0,0) {fluid1 = fluid2 = 0;}
  virtual ~LocalRiemannFluidStructure() { vf_ = 0; }

  void computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it, int Id = 0);
  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it);//TODO:not needed!

  void computeRiemannJacobian(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it, double* WstardU,int Id = 0);

  // Multi-Phase Riemann solvers (implemented here just to stop compiler's complaining...)
  void computeRiemannSolution(double *Vi, double *Vj,
                            int IDi, int IDj, double *nphi,
                            double *initWi, double *initWj,
                            double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj,
                            double dx[3], int it) {
    fprintf(stderr,"ERROR: Should not call the two-phase Riemann solver in LocalRiemannFluidStructure!\n");}

  void computeRiemannJacobian(double *Vi, double *Vj,
                                      int IDi, int IDj, double *nphi,
                                      double *Wi, double *Wj,
                                      double dx[3],int it,
                                      double* dWidWi,double*  dWidWj,
                                      double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: Should not call the two-phase Riemann solver in LocalRiemannFluidStructure!\n");}


private:
  void eriemannfs(double rhol, double ul, double pl,
                  double &rhoi, double ui, double &pi,
                  VarFcn *vf, int Id = 0); //note: ui shouldn't be changed. so the value (instead of reference) is used.

  void eriemannfs_grad(double rho, double u, double p,
                       double &rhoi, double ui, double &pi,
                       VarFcn *vf, double* dWdWi,int Id); //Caution: "ui" will not be modified!
  
  template <int d>
  friend class FSJac;
};

//------------------------------------------------------------------------------

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej,
                            double &weightj, int it, int Id)
{
  // Commented by Adam on 2010.08.17 because it has to handle dim > 5
  //int dim = 5;

  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver

  //---------------------------------------------------------------

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};


  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi,Id);
  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id); //caution: U_i will not be modified!

  Wstar[0]  = R_i;
  Wstar[1]  = vti[0]+U_i*nphi[0];
  Wstar[2]  = vti[1]+U_i*nphi[1];
  Wstar[3]  = vti[2]+U_i*nphi[2];
  Wstar[4]  = P_i;
  if(dim == 6)
    {
      Wstar[5]  = 0.0;// Boundary Condition: nuTilde = 0
    }
  else if(dim == 7)  // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[5] = 0.0;
      Wstar[6] = 0.0;
    }


  //---------------------------------------------------------------

  double Vi0[dim];
  for(int i=0; i<dim; i++)
    Vi0[i] = Vi[dim+i];

  vni = Vi0[1]*nphi[0]+Vi0[2]*nphi[1]+Vi0[3]*nphi[2];
  vti[0] = Vi0[1] - vni*nphi[0];
  vti[1] = Vi0[2] - vni*nphi[1];
  vti[2] = Vi0[3] - vni*nphi[2];

  R_1 = Vi0[0];
  U_1 = vni;
  P_1 = vf->getPressure(Vi0,Id);

  // U_i is the same.
  eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id); //caution: U_i will not be modified!

  Wstar[dim]    = R_i;
  Wstar[dim+1]  = vti[0]+U_i*nphi[0];
  Wstar[dim+2]  = vti[1]+U_i*nphi[1];
  Wstar[dim+3]  = vti[2]+U_i*nphi[2];
  Wstar[dim+4]  = P_i;
  if(dim == 6)
    {
      Wstar[dim+5]  = 0.0; // Boundary Condition: nuTilde = 0
    }
  else if(dim == 7) // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[dim+5]  = 0.0; 
      Wstar[dim+6]  = 0.0; 
    }

  //-----------------------------------------------------------------
/*
  if(it==1){
    weightj += 1.0;
    for (int k=0; k<dim; k++)
      rupdatej[k] += Wstar[k];  //TODO: rupdate is never used for FSI. (only used for MPF)
  }
*/
}

template<int dim>
class FSJac {
 
  int Id;
  VarFcn* vf;
  double U_i;
  LocalRiemannFluidStructure<dim>* ls; 
  public:
  FSJac(LocalRiemannFluidStructure<dim>* _ls,VarFcn* _vf, int _Id,double _ui) : vf(_vf), Id(_Id),ls(_ls),U_i(_ui) {}
  void Compute(const double u[3],double f[3]) const {

    ls->eriemannfs(u[0],u[1],u[2],f[0],U_i,f[2],vf,Id); //caution: U_i will not be modified!
  }
};

//------------------------------------------------------------------------------
template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannJacobian(double *Vi, double *Vstar,
                                      double *nphi, VarFcn *vf,
                                      double *Wstar, double *rupdatei,
                                      double &weighti, int it, double* dWstardU,int Id) {

  double P_1, U_1,R_1; // pass to 1D-FSI Riemann solver
  double P_i=Wstar[4], U_i, R_i=Wstar[0]; // solution given by 1D-FSI Riemann solver
  //---------------------------------------------------------------

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  double dWdW[9]={1,0,0,0,1,0,0,0,1};

  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi,Id);
  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  eriemannfs_grad(R_1,U_1,P_1,R_i,U_i,P_i,vf,dWdW,Id); //caution: U_i will not be modified!

  // Checking jacobian
  //double u[3] = {R_1,U_1,P_1};
  //double f[3] = {R_i,U_i,P_i};
  //DebugTools::CheckJacobian<3>(dWdW, u,f,FSJac<dim>(this,vf,Id,U_i), "FSJacobian\n-------------------------\n"); 

  if(dim == 6)
    {
      Wstar[5]  = 0.0;// Boundary Condition: nuTilde = 0
    }
  else if(dim == 7) // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[5] = 0.0;
      Wstar[6] = 0.0;
    }
 
  memset(dWstardU, 0, sizeof(double)*dim*dim);
  dWstardU[0] = dWdW[0];
  dWstardU[dim*4+4] = dWdW[8];

  for (int i = 0; i < 3; ++i) {
    dWstardU[(i+1)] = nphi[i]*dWdW[1];
    dWstardU[(i+1)+4*dim] = nphi[i]*dWdW[7];
  }

  dWstardU[4] = dWdW[2];
  dWstardU[dim*4] = dWdW[6]; 
}

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej,
                            double &weightj, int it)
{
  // Adam 2010.08.18
  // This function doesn't seem to be used anymore.
  // To be removed in a couple of months
  fprintf(stderr,"Oh Sorry ! Please uncomment the function (LocalRiemmannDesc.h:1596). I thought it wasn't needed anymore\n");
  exit(-1);
  /*
  // Commented by Adam on 2010.08.17 because it has to handle dim > 5
  //  int dim = 5;

  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};


  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi);

  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,tag); //caution: U_i will not be modified!

  Wstar[0]  = R_i;                     Wstar[dim]    = Wstar[0];
  Wstar[1]  = vti[0]+U_i*nphi[0];      Wstar[dim+1]  = Wstar[1];
  Wstar[2]  = vti[1]+U_i*nphi[1];      Wstar[dim+2]  = Wstar[2];
  Wstar[3]  = vti[2]+U_i*nphi[2];      Wstar[dim+3]  = Wstar[3];
  Wstar[4]  = P_i;                     Wstar[dim+4]  = Wstar[4];
  if(dim > 5)
    {
      Wstar[5]  = 0.0;                 Wstar[dim+5]  = 0.0; // Boundary Condition: nuTilde = 0
    }

  if(it==1){
    weightj += 1.0;
    for (int k=0; k<dim; k++)
      rupdatej[k] += Wstar[k];
  }
  */
}

//------------------------------------------------------------------------------

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::eriemannfs(double rho, double u, double p,
                                            double &rhoi, double ui, double &pi,
                                            VarFcn *vf, int Id) //Caution: "ui" will not be modified!
{
  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double gamma = vf->getGamma(Id);
  double pref  = vf->getPressureConstant(Id);

  if(u==ui){ // contact
    rhoi = rho;
    pi   = p;
    return;
  }

  if(ui<u){ // rarefaction
 /*   double power = 2*gamma/(gamma-1.0);
    double a = sqrt(gamma*(p+pref)/rho);
    double pbar = p + pref;
    pi = pbar*pow(0.5*(gamma-1.0)*(ui-u)/a + 1.0,power)-pref;
    rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma); */

    double power = gamma/(gamma-1.0);
    double a = sqrt(gamma*(p+pref)/rho);
    double pbar = p + pref;
    double dee = 0.5*(gamma-1.0)*(ui-u)/a + 1.0;
    pi = pbar*pow(dee*dee,power)-pref;
    rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
  }
  else{ // shock
    double temp = ((gamma+1.0)*rho*(ui-u)*(ui-u))/2.0;
    pi = p + 0.5*temp + sqrt(0.25*temp*temp + 2.0*gamma*temp*(p+pref)/(gamma+1.0));
    temp = (gamma-1.0)/(gamma+1.0);
    double pstarbar = pi + pref;
    double pbar = p + pref;
    rhoi = rho*(pstarbar/pbar+temp)/(temp*pstarbar/pbar+1);
  }
}

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::eriemannfs_grad(double rho, double u, double p,
                                                 double &rhoi, double ui, double &pi,
                                                 VarFcn *vf, double* dWidWi,int Id) //Caution: "ui" will not be modified!
{

  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double gamma = vf->getGamma(Id);
  double pref  = vf->getPressureConstant(Id);
  memset(dWidWi, 0,sizeof(double)*9);
  if(u==ui){ // contact
    dWidWi[0] = 1.0; 
    dWidWi[8] = 1.0;
    return;
  }

  double q = (gamma-1.0)/(gamma+1.0);
  if(ui<u){ // rarefaction
    double power = 2*gamma/(gamma-1.0);
    double a = sqrt(gamma*(p+pref)/rho);
    double pbar = p + pref;

    double s = 0.5*(gamma-1.0)*(ui-u)/a + 1.0;
//    if (s < 0.0) {
//      fprintf(stderr,"Warning: s (%lf) in fs_grad is < 0!\n",s);
//    }
    double eta = pbar*power*pow(s*s,q*0.5);
    double xi = eta*(-0.5/(a*a)*(gamma-1.0)*(ui-u));
 
    double dadp = 0.5/a*(gamma/rho), dadrho = -0.5*a/rho;
     
    // dpi/dp
    dWidWi[8] = (pi+pref)/pbar+xi*dadp;
    // dpidrho
    dWidWi[6] = xi*dadrho;
    // dpidu
    dWidWi[7] = eta*(-0.5/a*(gamma-1.0));

    double mu = rho/gamma*pow((pi+pref)/pbar, (1.0-gamma)/gamma); 
    dWidWi[2] = mu*(1.0/pbar*dWidWi[8]-(pi+pref)/(pbar*pbar));
    dWidWi[1] = mu*dWidWi[7]/pbar;
    dWidWi[0] = rhoi/rho+mu/pbar*dWidWi[6];
  }
  else{ // shock
    double power = 2*gamma/(gamma+1.0);
    double t = ((gamma+1)*rho*(ui-u)*(ui-u))/2.0;
    double pstarbar = pi + pref;
    double pbar = p + pref;
 
    double dtdrho = t/rho, dtdu = -(gamma+1.0)*rho*(ui-u);
    double xi = sqrt(0.5*t*t+power*t*pbar);
    double eta = 0.5+0.5/xi*(0.5*t+power*pbar);
    dWidWi[8] = 1.0+0.5/xi*power*t;
    dWidWi[7] = eta*dtdu;
    dWidWi[6] = eta*dtdrho;

    double s = q*pstarbar/pbar+1.0;
    double deriv = 1.0/(pbar*s)-(pstarbar/pbar+q)/(s*s)*(q/pbar);
    double deriv2 = -pstarbar/(pbar*pbar*s)+(pstarbar/pbar+q)*(q*pstarbar/(pbar*pbar))/(s*s);
    dWidWi[0] = rhoi/rho+rho*deriv*dWidWi[6];
    dWidWi[1] = rho*deriv*dWidWi[7];
    dWidWi[2] = rho*(deriv*dWidWi[8]+deriv2);
  }
}

//------------------------------------------------------------------------------
/*
inline
void LocalRiemannFluidStructure::eriemannfs(double rho, double u, double p,
                                            double &rhoi, double ui, double &pi,
                                            VarFcn *vf, int tag) //Caution: "ui" will not be modified!
{

  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double gamma = (tag<=0) ? vf->getGammabis() : vf->getGamma();
  double pref  = (tag<=0) ? vf->getPressureConstantbis() : vf->getPressureConstant();

  if(u==ui){ // contact
    rhoi = rho;
    pi   = p;
    return;
  }

  if(ui<u){ // rarefaction
    double power = 2*gamma/(gamma-1.0);
    double a = sqrt(gamma*(p+pref)/rho);
    double pbar = p + pref;
    pi = pbar*pow(0.5*(gamma-1.0)*(ui-u)/a + 1,power)-pref;
    rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
    if (pi>p) {fprintf(stderr,"ERROR: Wrong solution to FS Riemann problem. Aborting.\n"); exit(-1);}
  }
  else{ // shock
    double temp = 2.0/((gamma+1)*rho*(ui-u)*(ui-u));
    pi = p + 0.5*(1.0+sqrt(8.0*gamma*temp*(p+pref)/(gamma+1.0)+1.0))/temp;
    temp = (gamma-1.0)/(gamma+1.0);
    double pstarbar = pi + pref;
    double pbar = p + pref;
    rhoi = rho*(pstarbar/pbar+temp)/(temp*pstarbar/pbar+1);
    if (pi<p) {fprintf(stderr,"ERROR: Wrong solution to FS Riemann problem. Aborting.\n"); exit(-1);}
  }

}
*/
//----------------------------------------------------------------------------
#endif
