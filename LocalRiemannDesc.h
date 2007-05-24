#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include <LocalRiemann.h>
#include <VarFcn.h>
#include <math.h>

//------------------------------------------------------------------------------

class LocalRiemannGfmpGasGas : public LocalRiemann {

public:
	LocalRiemannGfmpGasGas();
	~LocalRiemannGfmpGasGas();

	void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
			double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it);
};

//------------------------------------------------------------------------------

inline
LocalRiemannGfmpGasGas::LocalRiemannGfmpGasGas() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

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
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class LocalRiemannGfmpTaitTait : public LocalRiemann {

public:
	LocalRiemannGfmpTaitTait();
	~LocalRiemannGfmpTaitTait();

	void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
      double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it);
};

//------------------------------------------------------------------------------

inline
LocalRiemannGfmpTaitTait::LocalRiemannGfmpTaitTait() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class LocalRiemannGfmparGasGas : public LocalRiemann {

public:
	LocalRiemannGfmparGasGas();
	~LocalRiemannGfmparGasGas();

	void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
      double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it);
};

//------------------------------------------------------------------------------

inline
LocalRiemannGfmparGasGas::LocalRiemannGfmparGasGas() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class LocalRiemannGfmparGasTait: public LocalRiemann {

public:
	LocalRiemannGfmparGasTait();
	~LocalRiemannGfmparGasTait();

	void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
      double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it);
};

//------------------------------------------------------------------------------

inline
LocalRiemannGfmparGasTait::LocalRiemannGfmparGasTait() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

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
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class LocalRiemannGfmparTaitTait: public LocalRiemann {

public:
	LocalRiemannGfmparTaitTait();
	~LocalRiemannGfmparTaitTait();

	void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
      double *rupdatei, double *rupdatej, double &weighti, double &weightj, int it);
};

//------------------------------------------------------------------------------

inline
LocalRiemannGfmparTaitTait::LocalRiemannGfmparTaitTait() : LocalRiemann()
{

}

//------------------------------------------------------------------------------

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
//------------------------------------------------------------------------------

#endif
