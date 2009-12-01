#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include <LocalRiemann.h>
#include <VarFcn.h>
#include "IoData.h"
#include "SparseGrid.h"
#include <cmath>

//----------------------------------------------------------------------------

class LocalRiemannGfmpGasGas : public LocalRiemann {

public:
  LocalRiemannGfmpGasGas(VarFcn *vf);
  ~LocalRiemannGfmpGasGas();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpGasGas::LocalRiemannGfmpGasGas(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasGas::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  if(Phii>=0.0){
    epsi =  1;
    epsj = -1;
  }else{
    epsi = -1;
    epsj =  1;
  }


  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpTaitTait : public LocalRiemann {

public:
  LocalRiemannGfmpTaitTait(VarFcn *vf);
  ~LocalRiemannGfmpTaitTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpTaitTait::LocalRiemannGfmpTaitTait(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  double a1 = vf_->getAlphaWater();
  double b1 = vf_->getBetaWater();
  double p1 = vf_->getPrefWater();
  double a2 = vf_->getAlphaWaterbis();
  double b2 = vf_->getBetaWaterbis();
  double p2 = vf_->getPrefWaterbis();

  double temp = 0;
	
  if(Phii>=0.0){
    epsi =  1;
    epsj = -1;
    temp = p2+a2*pow(Vj[0],b2);
    Wi[0] = pow((temp-p1)/a1,1.0/b1);
    temp = p2+a2*pow(Vj[5],b2);
    Wi[5] = pow((temp-p1)/a1,1.0/b1);
    temp = p1+a1*pow(Vi[0],b1);
    Wj[0] = pow((temp-p2)/a2,1.0/b2);
    temp = p1+a1*pow(Vi[5],b1);
    Wj[5] = pow((temp-p2)/a2,1.0/b2);
  }else{
    epsi = -1;
    epsj =  1;
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

class LocalRiemannGfmpJWLJWL : public LocalRiemann {

public:
  LocalRiemannGfmpJWLJWL(VarFcn *vf);
  ~LocalRiemannGfmpJWLJWL();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpJWLJWL::LocalRiemannGfmpJWLJWL(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  if(Phii>=0.0){
    epsi =  1;
    epsj = -1;
  }else{
    epsi = -1;
    epsj =  1;
  }


  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpGasJWL : public LocalRiemann {

public:
  LocalRiemannGfmpGasJWL(VarFcn *vf);
  ~LocalRiemannGfmpGasJWL();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpGasJWL::LocalRiemannGfmpGasJWL(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasJWL::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  if(Phii>=0.0){
    epsi =  1;
    epsj = -1;
  }else{
    epsi = -1;
    epsj =  1;
  }


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

class LocalRiemannGfmparGasGas : public LocalRiemann {

public:
  LocalRiemannGfmparGasGas(VarFcn *vf);
  ~LocalRiemannGfmparGasGas();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it);
void eriemann(double rhol, double ul, double pl, 
              double rhor, double ur, double pr, 
              double &pi, double &ui, double &rhoil, double &rhoir){
  F77NAME(eriemanngg)(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,vf_->getGammabis(),
             vf_->getPressureConstantbis(), vf_->getGamma(), vf_->getPressureConstant()); 
}

};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasGas::LocalRiemannGfmparGasGas(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasGas::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi,
		int &epsi, int &epsj,	double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, 
                double &weighti, double &weightj,
                double dx[3], int it)
{
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf_->getGamma();
  double pref1 = vf_->getPressureConstant();
  double gam2  = vf_->getGammabis();
  double pref2 = vf_->getPressureConstantbis();

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (Phii >= 0.0) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, Phij);
    P_1  = vf_->getPressure(Vi, Phii);

    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    epsi = 1;
    epsj = -1;

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
    P_2  = vf_->getPressure(Vi, Phii);
    P_1  = vf_->getPressure(Vj, Phij);

    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1);
    epsi = -1;
    epsj = 1;

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
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasTait: public LocalRiemann {

public:
  LocalRiemannGfmparGasTait(VarFcn *vf);
  ~LocalRiemannGfmparGasTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, 
                            double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasTait::LocalRiemannGfmparGasTait(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasTait::computeRiemannSolution(double *Vi, double *Vj,
	  double Phii, double Phij, double *nphi,
	  int &epsi, int &epsj, double *Wi, double *Wj,
          double *rupdatei, double *rupdatej, double &weighti, double &weightj,
          double dx[3], int it)
{
  int dim = 5;

  double alpha   = vf_->getAlphaWater();
  double beta    = vf_->getBetaWater();
  double pref    = vf_->getPrefWater();
  double gam     = vf_->getGamma();

  double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (Phii >= 0.0) {
    // cell j is gas
    // cell i is tait
    R_g  = Vj[0];     R_w  = Vi[0];
    U_g  = vnj;       U_w  = vni;
    P_g  = vf_->getPressure(Vj, Phij);
    P_w  = vf_->getPressure(Vi, Phii);

    F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    epsi = 1;
    epsj = -1;

    Wi[0]  = R_ir;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];
    T_w  = vf_->computeTemperature(Wi, Phij);
    //T_w  = vf_->computeTemperature(Vi, Phii);
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
    P_g  = vf_->getPressure(Vi, Phii);
    P_w  = vf_->getPressure(Vj, Phij);

    F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    epsi = -1;
    epsj = 1;

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
    T_w  = vf_->computeTemperature(Wj, Phii);
    //T_w  = vf_->computeTemperature(Vj, Phij);
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
//----------------------------------------------------------------------------

class LocalRiemannGfmparTaitTait: public LocalRiemann {

public:
  LocalRiemannGfmparTaitTait(VarFcn *vf);
  ~LocalRiemannGfmparTaitTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparTaitTait::LocalRiemannGfmparTaitTait(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it)
{

  int dim = 5;

  double alpha1   = vf_->getAlphaWater();
  double beta1    = vf_->getBetaWater();
  double pref1    = vf_->getPrefWater();
  double alpha2   = vf_->getAlphaWaterbis();
  double beta2    = vf_->getBetaWaterbis();
  double pref2    = vf_->getPrefWaterbis();

  //double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_1, P_2, R_1, R_2, U_1, U_2, T_1, T_2;
  double P_i, U_i, R_i1, R_i2;
  
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};
  
  if (Phii >= 0.0) {
    // cell j is tait2
    // cell i is tait1
    R_1  = Vi[0];     R_2  = Vj[0];
    U_1  = vni;       U_2  = vnj;
    P_1  = vf_->getPressure(Vi, Phii);
    P_2  = vf_->getPressure(Vj, Phij);
    T_1  = vf_->computeTemperature(Vi, Phii);
    T_2  = vf_->computeTemperature(Vj, Phij);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);
    epsi = 1;
    epsj = -1;
    
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
    P_1  = vf_->getPressure(Vj, Phij);
    P_2  = vf_->getPressure(Vi, Phii);
    T_1  = vf_->computeTemperature(Vj, Phij);
    T_2  = vf_->computeTemperature(Vi, Phii);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1);
    epsi = -1;
    epsj = 1;
    
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
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparJWLJWL : public LocalRiemann {

public:
  LocalRiemannGfmparJWLJWL(VarFcn *vf);
  ~LocalRiemannGfmparJWLJWL();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){ 
    /*eriemannjj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir);*/ }

private:
  void eriemannjj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparJWLJWL::LocalRiemannGfmparJWLJWL(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi,
		int &epsi, int &epsj,	double *Wi, double *Wj,
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

  if (Phii >= 0.0) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, Phij);
    P_1  = vf_->getPressure(Vi, Phii);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1);

    epsi = 1;
    epsj = -1;

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
    P_2  = vf_->getPressure(Vi, Phii);
    P_1  = vf_->getPressure(Vj, Phij);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1);

    epsi = -1;
    epsj = 1;

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
  double omegal = vf_->getOmegabis();
  double omegar = vf_->getOmega();
  double omp1ooml = (omegal+1.0)/omegal;
  double omp1oomr = (omegar+1.0)/omegar;
  double frhol = vf_->computeFrho(1.0/vl,-1.0);
  double frhor = vf_->computeFrho(1.0/vr,-1.0);
  double frhoil = frhol;
  double frhoir = frhor;
  double frhopil = vf_->computeFrhop(1.0/vl,1.0);
  double frhopir = vf_->computeFrhop(1.0/vr,1.0);
//check vacuum ?

//start newton iteration loop
  while(!convergence){
    //fprintf(stdout, "%e %e %e %e\n", 1.0/vil, 1.0/vir, ui, pi);

  //compute left  term (shock or rarefaction)
    if( vil < vl){
      //fprintf(stdout, "leftshock\n");
      frhoil  = vf_->computeFrho(1.0/vil,-1.0);
      frhopil = vf_->computeFrhop(1.0/vil,-1.0);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      //fprintf(stdout, "leftraref\n");
      rarefactionJWL2ndOrder(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil);
    }
  //compute right term (shock or rarefaction)
    if( vir < vr){
      //fprintf(stdout, "rightshock\n");
      frhoir  = vf_->computeFrho(1.0/vir,1.0);
      frhopir = vf_->computeFrhop(1.0/vir,1.0);
      shockJWL(1.0, omegar, omp1oomr, frhor, frhoir, frhopir, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      //fprintf(stdout, "rightraref\n");
      rarefactionJWL2ndOrder(1.0, vr, ur, pr, vir, uir, pir, duir, dpir);
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

class LocalRiemannGfmparGasJWL : public LocalRiemann {

  MultiFluidData::RiemannComputation riemannComputationType;
  SparseGrid *tabulation;

public:
  LocalRiemannGfmparGasJWL(VarFcn *vf, IoData &iod);
  ~LocalRiemannGfmparGasJWL(){ delete tabulation; }

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it);
  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){
    /*eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir);*/ }
  void eriemanngj_wrapper(double *in, double *res, double *phi);

private:
  bool eriemanngj_selector(double rhol, double ul, double pl, 
                           double rhor, double ur, double pr, 
                           double &pi, double &ui, double &rhoil, double &rhoir);
  bool eriemanngj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir);

  void riemannInvariantGeneralTabulation(double *in, double *res);
  bool vacuum(const double rhol, const double ul, const double pl,
              const double rhor, const double ur, const double pr,
              double vacuumValues[6]);
  double jwlZeroSoundSpeedJwlDensity(const double density, const double pressure);
  double sgZeroDensityPJwlDensity(const double density, const double pressure,
                                  const double rho_c0);
  double pressureEqGasDensity(const double gasDensity, const double gasPressure,
                              const double jwlDensity, const double jwlPressure,
                              const double interfacialJwlDensity);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasJWL::LocalRiemannGfmparGasJWL(VarFcn *vf, IoData &iod) : 
LocalRiemann()
{
  vf_ = vf;
  riemannComputationType = iod.mf.riemannComputation;
  if(riemannComputationType == MultiFluidData::TABULATION2){
    double *refIn = new double[2]; double *refOut = new double[1];
    refIn[0] = iod.ref.rv.density;
    //refIn[1] = iod.ref.rv.entropy;
    refIn[1] = pow(iod.ref.rv.density,-iod.eqs.fluidModel2.jwlModel.omega)*iod.ref.rv.velocity*iod.ref.rv.velocity;
    refOut[0] = iod.ref.rv.velocity;
    fprintf(stdout, "adim SparseGrid %e %e %e\n", refIn[0], refIn[1], refOut[0]);

    tabulation = new SparseGrid;
    tabulation->readFromFile(refIn, refOut);
    //tabulation->scaleGrid(refIn, refOut);
  }else if(riemannComputationType == MultiFluidData::TABULATION5){

    double *refIn = new double[5]; double *refOut = new double[2];

    refIn[0] = iod.ref.rv.density;
    refIn[1] = iod.ref.rv.pressure;
    refIn[2] = iod.ref.rv.density;
    refIn[3] = iod.ref.rv.pressure;
    refIn[4] = iod.ref.rv.velocity;

    refOut[0] = iod.ref.rv.density;
    refOut[1] = iod.ref.rv.density;
    //fprintf(stdout, "adim SparseGrid %e %e %e %e %e %e %e\n", refIn[0],refIn[1],refIn[2],refIn[3],refIn[4],refOut[0],refOut[1]);
    tabulation = new SparseGrid;
    tabulation->readFromFile(refIn, refOut);
    //tabulation->scaleGrid(refIn, refOut);
  }
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi,
		int &epsi, int &epsj,	double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it)
{

  bool computeRiemannSolutionGasJWLimplemented = false;
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (Phii >= 0.0) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, Phij);
    P_1  = vf_->getPressure(Vi, Phii);

    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1); 
    epsi = 1;
    epsj = -1;

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
    P_2  = vf_->getPressure(Vi, Phii);
    P_1  = vf_->getPressure(Vj, Phij);

    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1); 
    epsi = -1;
    epsj = 1;

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
void LocalRiemannGfmparGasJWL::eriemanngj_wrapper(double *in, double *res, double *phi)
{

  double dummy1, dummy2, dummy3;
  eriemanngj(in[0], 0.0, in[1], in[2], in[4], in[3], dummy1, dummy2, res[0], res[1]);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::eriemanngj_selector(double rhol, double ul, double pl, 
                                                   double rhor, double ur, double pr, 
                                                   double &pi, double &ui,  
                                                   double &rhoil, double &rhoir)
{

  if(riemannComputationType==MultiFluidData::TABULATION5){
    double *in  = new double[5]; in[0]=rhol;in[1]=pl;in[2]=rhor;in[3]=pr;in[4]=ur-ul;
    double *res = new double[2]; res[0]=0.0;res[1]=0.0;
    tabulation->interpolate(1,&in,&res);
    rhoil = fmax(res[0],0.0); rhoir = fmax(res[1],0.0);
    double d[2]; //dummy variable
    double uir, pir, uil, pil;

    double omegal = vf_->getOmega();
    double omp1ooml = (omegal+1.0)/omegal;
    double frhol = vf_->computeFrho(rhol,-1.0);
    double gamr = vf_->getGamma();
    double prefr = vf_->getPressureConstant();
    double gam1r = vf_->getGamma()-1.0;
    double gamogam1r = gamr/gam1r;
    double Vr[5] = { rhor, ur, 0.0, 0.0, pr };
    double cr = vf_->computeSoundSpeed(Vr,1.0);

    if( rhoil > rhol){
      double frhoil  = vf_->computeFrho(rhoil,-1.0);
      double frhopil = vf_->computeFrhop(rhoil,-1.0);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1]);
    }else
      rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1], riemannComputationType,1);
    if( rhoir > rhor)
      shockGAS(1.0, gamogam1r, prefr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);
    else
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);

    ui = 0.5*(uil+uir);
    pi = 0.5*(pil+pir);
    // not checking for vacuum!
  }
  else
    eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::eriemanngj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir){
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
  double eps = 1.e-6;
  int MaxIts = 400;
  int it = 0;
  double relaxationFactorJwl = 0.85; // must be between 0 and 1
  double relaxationFactorGas = 0.85; // must be between 0 and 1
  int count = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;

  double omegal = vf_->getOmega();
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,-1.0);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,1.0);


  double gamr = vf_->getGamma();
  double prefr = vf_->getPressureConstant();
  double gam1r = vf_->getGamma()-1.0;
  double gamogam1r = gamr/gam1r;
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,1.0);

//check vacuum
  if(verbose>4) fprintf(stdout, "checking vacuum possibilities\n");
  double vacuumValues[6]; /* rhoil, uil, pil, rhoir, uir, pir */
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

//start newton iteration loop
  while(!convergence){
    if(verbose>0) fprintf(stdout, "\n");

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      if(verbose>0) fprintf(stdout, "shockJWL\n");
      frhoil  = vf_->computeFrho(1.0/vil,-1.0);
      frhopil = vf_->computeFrhop(1.0/vil,-1.0);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      if(verbose>0) fprintf(stdout, "rarefactionJWL\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType,1);
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
    if(vil>1.0/vacuumValues[0]){
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
    frhoil  = vf_->computeFrho(1.0/vil,-1.0);
    frhopil = vf_->computeFrhop(1.0/vil,-1.0);
    shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
  }else rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType,0);
  
  if( vir < vr) shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
  else rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);

  rhoil = 1.0/vil;
  rhoir = 1.0/vir;
  ui    = 0.5*(uil+uir);
  pi    = 0.5*(pil+pir);

  if(convergence){
    if(verbose>-1) fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
  }else{
    if(verbose>-1) fprintf(stdout, "riemann solver did not converged\n");
    if(verbose>-1) fprintf(stdout, "Warning: solution will be state given by vacuum\n");
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
  double cr = vf_->computeSoundSpeed(Vr,1.0);
  rarefactionGAS(1.0, vf_->getGamma(), vf_->getGamma()-1.0, vf_->getPressureConstant(), cr, 1.0/rhor, ur, pr, 1.0/rhoir_vac, uir, pir, duir, dpir);
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
    double tol = 1.0e-6;
    int it = 0, maxIt = 100;
    double relaxation = 1.0;

    double entropy = vf_->computeEntropy(density,pressure,-1.0);

    double xn = density;
    double xnm1 = xn;
    double dx, fn, dfn;
    double lowerBound = 0.0; // this equation may have more than one root (xn = 0 is always root)
                             // but we want the larger positive root below density
                             // so we use this lowerBound to reduce the search domain.

    while(!convergence){
        fn = (vf_->getOmega()+1.0)*entropy*pow(xn,vf_->getOmega()) + vf_->computeExponentials2(xn);
        if(fn<0.0) lowerBound = xn;
        dfn = entropy*(vf_->getOmega()+1.0)*vf_->getOmega()*pow(xn,vf_->getOmega()-1) + vf_->computeDerivativeOfExponentials2(xn);
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
    return gasDensity/pow(1+gasPressure/vf_->getPressureConstant(),1.0/vf_->getGamma());

  double jwlEntropy = vf_->computeEntropy(jwlDensity, jwlPressure, -1.0);
  double gasEntropy = vf_->computeEntropy(gasDensity, gasPressure, 1.0);
  return pow((jwlEntropy*pow(interfacialJwlDensity,vf_->getOmega()+1.0)+vf_->computeExponentials(interfacialJwlDensity)+vf_->getPressureConstant())/gasEntropy,1.0/vf_->getGamma());

}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::sgZeroDensityPJwlDensity(const double density, const double pressure,
                                                          const double rho_c0)
{

  int verbose=0;
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - density=%e and pressure=%e and rho_c0=%e\n", density, pressure, rho_c0);
  double entropy = vf_->computeEntropy(density,pressure,-1.0);

  double fn = entropy*pow(rho_c0>0.0 ? rho_c0 : 1.e-14,vf_->getOmega()+1.0) + vf_->computeExponentials(rho_c0>0.0 ? rho_c0 : 1.e-14) + vf_->getPressureConstant();
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - fn(max(rho_c0,0)) = %e\n", fn);
  if(fn>0.0) return -1.0;

  bool convergence = false;
  double tol = 1.0e-6;
  int it =0, maxIt = 100;
  double relaxation = 1.0;


  double xn = density;
  double xnm1 = xn;
  double dx=0, dfn;

  while(!convergence){
    fn = entropy*pow(xn,vf_->getOmega()+1.0) + vf_->computeExponentials(xn) + vf_->getPressureConstant();
    dfn = entropy*(vf_->getOmega()+1.0)*pow(xn,vf_->getOmega()) + vf_->computeDerivativeOfExponentials(xn);
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
void LocalRiemannGfmparGasJWL::riemannInvariantGeneralTabulation(double *in, 
                                                                 double *res){
  //fprintf(stdout, "in-value = %e %e\n", in[0],in[1]);
  tabulation->interpolate(1,&in,&res);

}

//----------------------------------------------------------------------------

#endif
