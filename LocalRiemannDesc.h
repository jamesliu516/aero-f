#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include <LocalRiemann.h>
#include <VarFcn.h>
#include <math.h>

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

public:
  LocalRiemannGfmparGasJWL(VarFcn *vf);
  ~LocalRiemannGfmparGasJWL();

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

private:
  void eriemanngj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasJWL::LocalRiemannGfmparGasJWL(VarFcn *vf) : LocalRiemann()
{
  vf_ = vf;
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

    eriemanngj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1); 
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

    eriemanngj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1); 
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
void LocalRiemannGfmparGasJWL::eriemanngj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-8;
  int MaxIts = 400;
  int it = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;

  double omegal = vf_->getOmega();
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,-1.0);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,1.0);

  fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
  fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);

  double gamr = vf_->getGamma();
  double prefr = vf_->getPressureConstant();
  double gam1r = vf_->getGamma()-1.0;
  double gamogam1r = gamr/gam1r;
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,1.0);
//check vacuum ?

//start newton iteration loop
  while(!convergence){

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      fprintf(stdout, "shock\n");
      frhoil  = vf_->computeFrho(1.0/vil,-1.0);
      frhopil = vf_->computeFrhop(1.0/vil,-1.0);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      fprintf(stdout, "rarefaction\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }

    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);

  //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    solve2x2System(jacobian,function,increment);
    increment[0] /= 5;
    increment[1] /= 5;

    fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

  //update values and check bounds
    fprintf(stdout, "1 -- vil = %e and vir = %e\n", vil, vir);
    fprintf(stdout, "11-- vil = %e and vir = %e\n", vil-increment[0], vir-increment[1]);
    if(vil - increment[0] < 0.0)
      increment[0] = 0.5*vil;
    if(vir - increment[1] < 0.0)
      increment[1] = 0.5*vir;
    
    vil -= increment[0];
    vir -= increment[1];
    fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp)
        vil = 0.5*(vil+increment[0]+temp);
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      double temp = (gamr-1.0)/(gamr+1.0)*vr;
      if(vir<temp)
        vir = 0.5*(vir+increment[1]+temp);
    }
    fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop
  if(convergence){
    fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
    rhoil = 1.0/vil;
    rhoir = 1.0/vir;
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);
  }else{
    fprintf(stdout, "riemann solver did not converged\n");
    exit(1);
  }
  fprintf(stdout, "rhoil, rhoir, ui, pi = %e %e %e %e\n", rhoil,rhoir, ui, pi);
  exit(1);


}
//----------------------------------------------------------------------------


#endif
