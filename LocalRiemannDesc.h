#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include <LocalRiemann.h>
#include <VarFcn.h>
#include <math.h>

//----------------------------------------------------------------------------

class LocalRiemannGfmpGasGas : public LocalRiemann {

public:
  LocalRiemannGfmpGasGas();
  ~LocalRiemannGfmpGasGas();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi, VarFcn *vf,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpGasGas::LocalRiemannGfmpGasGas() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasGas::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi, VarFcn *vf,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
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
  LocalRiemannGfmpTaitTait();
  ~LocalRiemannGfmpTaitTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpTaitTait::LocalRiemannGfmpTaitTait() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi, VarFcn *vf,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  double a1 = vf->getAlphaWater();
  double b1 = vf->getBetaWater();
  double p1 = vf->getPrefWater();
  double a2 = vf->getAlphaWaterbis();
  double b2 = vf->getBetaWaterbis();
  double p2 = vf->getPrefWaterbis();

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
  LocalRiemannGfmpJWLJWL();
  ~LocalRiemannGfmpJWLJWL();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi, VarFcn *vf,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpJWLJWL::LocalRiemannGfmpJWLJWL() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi, VarFcn *vf,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
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
  LocalRiemannGfmpGasJWL();
  ~LocalRiemannGfmpGasJWL();

  void computeRiemannSolution(double *Vi, double *Vj,
                              double Phii, double Phij, double *nphi, VarFcn *vf,
                              int &epsi, int &epsj, double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmpGasJWL::LocalRiemannGfmpGasJWL() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpGasJWL::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi, VarFcn *vf,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
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
  LocalRiemannGfmparGasGas();
  ~LocalRiemannGfmparGasGas();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);
void eriemann(double rhol, double ul, double pl, 
              double rhor, double ur, double pr, 
              double &pi, double &ui, double &rhoil, double &rhoir,
              VarFcn *vf){ 
  F77NAME(eriemanngg)(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,vf->getGammabis(),
                vf->getPressureConstantbis(), vf->getGamma(), vf->getPressureConstant()); 
}

};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasGas::LocalRiemannGfmparGasGas() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasGas::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi, VarFcn *vf, 
		int &epsi, int &epsj,	double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
{
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf->getGamma();
  double pref1 = vf->getPressureConstant();
  double gam2  = vf->getGammabis();
  double pref2 = vf->getPressureConstantbis();

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  if (Phii >= 0.0) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf->getPressure(Vj, Phij);
    P_1  = vf->getPressure(Vi, Phii);

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
    P_2  = vf->getPressure(Vi, Phii);
    P_1  = vf->getPressure(Vj, Phij);

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

  if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasTait: public LocalRiemann {

public:
  LocalRiemannGfmparGasTait();
  ~LocalRiemannGfmparGasTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasTait::LocalRiemannGfmparGasTait() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasTait::computeRiemannSolution(double *Vi, double *Vj,
	  double Phii, double Phij, double *nphi, VarFcn *vf,
	  int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
{
  int dim = 5;

  double alpha   = vf->getAlphaWater();
  double beta    = vf->getBetaWater();
  double pref    = vf->getPrefWater();
  double gam     = vf->getGamma();

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
    P_g  = vf->getPressure(Vj, Phij);
    P_w  = vf->getPressure(Vi, Phii);

    F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    epsi = 1;
    epsj = -1;

    Wi[0]  = R_ir;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];
    T_w  = vf->computeTemperature(Wi, Phij);
    //T_w  = vf->computeTemperature(Vi, Phii);
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
    P_g  = vf->getPressure(Vi, Phii);
    P_w  = vf->getPressure(Vj, Phij);

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
    T_w  = vf->computeTemperature(Wj, Phii);
    //T_w  = vf->computeTemperature(Vj, Phij);
    Wj[4]  = T_w;                     Wj[dim+4]  = Wj[4];

  }
  
  if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparTaitTait: public LocalRiemann {

public:
  LocalRiemannGfmparTaitTait();
  ~LocalRiemannGfmparTaitTait();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparTaitTait::LocalRiemannGfmparTaitTait() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    double Phii, double Phij, double *nphi, VarFcn *vf,
    int &epsi, int &epsj, double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
{

  int dim = 5;

  double alpha1   = vf->getAlphaWater();
  double beta1    = vf->getBetaWater();
  double pref1    = vf->getPrefWater();
  double alpha2   = vf->getAlphaWaterbis();
  double beta2    = vf->getBetaWaterbis();
  double pref2    = vf->getPrefWaterbis();

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
    P_1  = vf->getPressure(Vi, Phii);
    P_2  = vf->getPressure(Vj, Phij);
    T_1  = vf->computeTemperature(Vi, Phii);
    T_2  = vf->computeTemperature(Vj, Phij);
    
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
    P_1  = vf->getPressure(Vj, Phij);
    P_2  = vf->getPressure(Vi, Phii);
    T_1  = vf->computeTemperature(Vj, Phij);
    T_2  = vf->computeTemperature(Vi, Phii);
    
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

if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparJWLJWL : public LocalRiemann {

public:
  LocalRiemannGfmparJWLJWL();
  ~LocalRiemannGfmparJWLJWL();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);

  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir,
                VarFcn *vf){ /*eriemannjj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,vf);*/ }

private:
  void eriemannjj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  VarFcn *vf);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparJWLJWL::LocalRiemannGfmparJWLJWL() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi, VarFcn *vf, 
		int &epsi, int &epsj,	double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
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
    P_2  = vf->getPressure(Vj, Phij);
    P_1  = vf->getPressure(Vi, Phii);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,vf); // vf->omegap is for left, and vf->omega is for right

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
    P_2  = vf->getPressure(Vi, Phii);
    P_1  = vf->getPressure(Vj, Phij);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,vf);

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

  if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparJWLJWL::eriemannjj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          VarFcn *vf){

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
  double omegal = vf->getOmegabis();
  double omegar = vf->getOmega();
  double omp1ooml = (omegal+1.0)/omegal;
  double omp1oomr = (omegar+1.0)/omegar;
  double frhol = vf->computeFrho(1.0/vl,-1.0);
  double frhor = vf->computeFrho(1.0/vr,-1.0);
  double frhoil = frhol;
  double frhoir = frhor;
  double frhopil = vf->computeFrhop(1.0/vl,1.0);
  double frhopir = vf->computeFrhop(1.0/vr,1.0);
//check vacuum ?

//start newton iteration loop
  while(!convergence){
    //fprintf(stdout, "%e %e %e %e\n", 1.0/vil, 1.0/vir, ui, pi);

  //compute left  term (shock or rarefaction)
    if( vil < vl){
      //fprintf(stdout, "leftshock\n");
      frhoil  = vf->computeFrho(1.0/vil,-1.0);
      frhopil = vf->computeFrhop(1.0/vil,-1.0);
      shockJWL(vf, -1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      //fprintf(stdout, "leftraref\n");
      riemannInvariant(vf, -1.0, vl, ul, pl, vil, uil, pil, duil, dpil);
    }
  //compute right term (shock or rarefaction)
    if( vir < vr){
      //fprintf(stdout, "rightshock\n");
      frhoir  = vf->computeFrho(1.0/vir,1.0);
      frhopir = vf->computeFrhop(1.0/vir,1.0);
      shockJWL(vf, 1.0, omegar, omp1oomr, frhor, frhoir, frhopir, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      //fprintf(stdout, "rightraref\n");
      riemannInvariant(vf, 1.0, vr, ur, pr, vir, uir, pir, duir, dpir);
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
  LocalRiemannGfmparGasJWL();
  ~LocalRiemannGfmparGasJWL();

void computeRiemannSolution(double *Vi, double *Vj,
                            double Phii, double Phij, double *nphi, VarFcn *vf,
                            int &epsi, int &epsj, double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj, int it);
  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir,
                VarFcn *vf){ /*eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,vf);*/ }

private:
  void eriemanngj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  VarFcn *vf);
};

//----------------------------------------------------------------------------

inline
LocalRiemannGfmparGasJWL::LocalRiemannGfmparGasJWL() : LocalRiemann()
{

}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmparGasJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	double Phii, double Phij, double *nphi, VarFcn *vf, 
		int &epsi, int &epsj,	double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it)
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
    P_2  = vf->getPressure(Vj, Phij);
    P_1  = vf->getPressure(Vi, Phii);

    eriemanngj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,vf); 
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
    P_2  = vf->getPressure(Vi, Phii);
    P_1  = vf->getPressure(Vj, Phij);

    eriemanngj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,vf); 
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

  if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::eriemanngj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          VarFcn *vf){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-4;
  int MaxIts = 400;
  int it = 0;

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;

  double omegal = vf->getOmega();
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf->computeFrho(1.0/vl,-1.0);
  double frhoil = frhol;
  double frhopil = vf->computeFrhop(1.0/vl,1.0);

  //fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
  //fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);

  double gamr = vf->getGamma();
  double prefr = vf->getPressureConstant();
  double gam1r = vf->getGamma()-1.0;
  double gamogam1r = gamr/gam1r;
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf->computeSoundSpeed(Vr,1.0);
//check vacuum ?

//start newton iteration loop
  while(!convergence){

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      //fprintf(stdout, "shock\n");
      frhoil  = vf->computeFrho(1.0/vil,-1.0);
      frhopil = vf->computeFrhop(1.0/vil,-1.0);
      shockJWL(vf, -1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      //fprintf(stdout, "rarefaction\n");
      riemannInvariant(vf, -1.0, vl, ul, pl, vil, uil, pil, duil, dpil);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      riemannInvariantGAS(vf, 1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);
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
    increment[0] /= 5;
    increment[1] /= 5;

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
      double temp = (gamr-1.0)/(gamr+1.0)*vr;
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
  //fprintf(stdout, "rhoil, rhoir, ui, pi = %e %e %e %e\n", rhoil,rhoir, ui, pi);


}
//----------------------------------------------------------------------------


class LocalRiemannFluidStructure : public LocalRiemann {

public:
  LocalRiemannFluidStructure();
  ~LocalRiemannFluidStructure();

void computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei, 
                            double &weighti, int it);
void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei, 
                            double &weighti, int it);

private:
  void eriemannfs(double rhol, double ul, double pl,
                  double &rhoi, double ui, double &pi,
                  VarFcn *vf); //note: ui shouldn't be changed. so the value (instead of reference) is used.
  void eriemannfs(double rhol, double ul, double pl,
                  double &rhoi, double ui, double &pi,
                  VarFcn *vf, int tag); //note: ui shouldn't be changed. so the value (instead of reference) is used.
};

//------------------------------------------------------------------------------

inline
LocalRiemannFluidStructure::LocalRiemannFluidStructure() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

inline
void LocalRiemannFluidStructure::computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej, 
                            double &weightj, int it)
{

  int dim = 5;
	
  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};


  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi);

  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf); //caution: U_i will not be modified!

  Wstar[0]  = R_i;                     Wstar[dim]    = Wstar[0];
  Wstar[1]  = vti[0]+U_i*nphi[0];      Wstar[dim+1]  = Wstar[1];
  Wstar[2]  = vti[1]+U_i*nphi[1];      Wstar[dim+2]  = Wstar[2];
  Wstar[3]  = vti[2]+U_i*nphi[2];      Wstar[dim+3]  = Wstar[3];
  Wstar[4]  = P_i;                     Wstar[dim+4]  = Wstar[4];

  if(it==1){
    weightj += 1.0;
    for (int k=0; k<5; k++)
      rupdatej[k] += Wstar[k];
  }

}

//------------------------------------------------------------------------------

inline
void LocalRiemannFluidStructure::computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej,
                            double &weightj, int it)
{

  int dim = 5;

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

  if(it==1){
    weightj += 1.0;
    for (int k=0; k<5; k++)
      rupdatej[k] += Wstar[k];
  }

}

//------------------------------------------------------------------------------

inline
void LocalRiemannFluidStructure::eriemannfs(double rho, double u, double p,
                                            double &rhoi, double ui, double &pi,
                                            VarFcn *vf) //Caution: "ui" will not be modified!
{

  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double gamma = vf->getGamma();
  double pref  = vf->getPressureConstant();

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
/*    double temp = 2.0/((gamma+1)*rho*(ui-u)*(ui-u));
    pi = p + 0.5*(1.0+sqrt(8.0*gamma*temp*(p+pref)/(gamma+1.0)+1.0))/temp;
*/
    double temp = ((gamma+1)*rho*(ui-u)*(ui-u))/2.0;
    pi = p + 0.5*temp + sqrt(0.25*temp*temp + 2.0*gamma*temp*(p+pref)/(gamma+1.0));
    temp = (gamma-1.0)/(gamma+1.0);
    double pstarbar = pi + pref;
    double pbar = p + pref;
    rhoi = rho*(pstarbar/pbar+temp)/(temp*pstarbar/pbar+1);
    if (pi<p) {fprintf(stderr,"ERROR: Wrong solution to FS Riemann problem. Aborting.\n"); exit(-1);}
  }
}

//------------------------------------------------------------------------------

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




#endif
