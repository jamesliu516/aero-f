#include <FluxFcnDescTait.h>

#include <LinkF77.h>

#include <stdlib.h>
#include <stdio.h>



//------------------------------------------------------------------------------
// fortran routines located in f77src folder

extern "C" {
  void F77NAME(roeflux1water)(const double&, const double&, const double&,
                              const double&, const double&, double*,
                              const double&, double*, double*, double*);
  void F77NAME(roeflux5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*,
                         const double&, const double&, const double&, const double&, const int&);
  void F77NAME(roejac5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, 
			 const double&, const double&, const double&, 
			 const double&, const int&);
  void F77NAME(genbcfluxtait)(const int&, const double&, const double&, const double&, 
			 const double&,  double*, const double&, double*, double*, double*);
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
inline
void flux3Dwater(int type, VarFcnBase *vf, double *normal, double normalVel, double *V, double *Ub, double *flux){
 
  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double VV[5];
  double machb = vf->computeMachNumber(Vb);
  double cb = vf->computeSoundSpeed(Vb);
  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2];//-nVel;
  double un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
  double c = vf->computeSoundSpeed(V);
  double rhoun, p, rhoE;
  
  if(un==0.0){                   // SLIP-WALL LIKE, if boundary velocity is in the same plane as the face of the boundary
    VV[0]=vf->getDensity(Vb);
    VV[1]=0.0; VV[2]=0.0; VV[3]=0.0; VV[4]=0.0;

    rhoun=0.0;
    p = vf->getPressure(VV);
    rhoE = 0.0;
  }else{
    if(un<0.0){ 		//INLET, as the normal is going outward.
      if(-un-c>0.0){ 			//SUPERSONIC
	VV[0]=vf->getDensity(Vb);
	VV[1]=Vb[1];
	VV[2]=Vb[2];
	VV[3]=Vb[3];
	VV[4]=vf->computeTemperature(Vb);
      }else{        			//SUBSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
      }
    }
    else{                  //OUTLET
      if(un-c>0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
      }else{                              //SUBSONIC
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
        VV[4] = vf->computeTemperature(V);
      }  
    }
    rhoun = VV[0] * ( VV[1]*n[0] + VV[2]*n[1] + VV[3]*n[2] - nVel );
    p = vf->getPressure(VV);
    rhoE = vf->computeRhoEnergy(VV);
  }
  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*VV[1] + p*n[0]);
  flux[2] = S * (rhoun*VV[2] + p*n[1]);
  flux[3] = S * (rhoun*VV[3] + p*n[2]);
  flux[4] = S * ((rhoE + p) * rhoun/VV[0] + p*nVel); 

}

//------------------------------------------------------------------------------



template<int dim>
inline
void jacflux3Dwater(int type, VarFcnBase *vf, FluxFcnBase::Type localTypeJac, double *normal,
		      double normalVel, double *V, double *Ub, double *jac){

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double VV[dim];
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2]; 
  double cb = vf->computeSoundSpeed(Vb);
  double un = V[1]*n[0]+V[2]*n[1]+V[3]*n[2];
  double c  = vf->computeSoundSpeed(V);

  double _dfdV[dim][dim];
  double *dfdV = reinterpret_cast<double *>(_dfdV);
  for(int k=0; k<dim*dim; ++k)
    dfdV[k] = 0.0;


  if(un == 0.0){                 //SLIP-WALL LIKE
    
    //nothing to do, the flux has only pressure terms
    //which depend on density only, which is the one 
    //at the boundary (not inside the volume)

  }else{  
    if(un < 0.0){ 		//INLET, as the normal is going outward.
      if(-un-c > 0.0){                  //SUPERSONIC
	
	//nothing to do: jacobian is null
	
      }else{                         //SUBSONIC
	
	//derivative wrt pressure
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	double unp = VV[1]*n[0] +VV[2]*n[1] + VV[3]*n[2] - nVel;
	_dfdV[0][0] = S * unp;
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[2][0] = S * (VV[2]*unp + cp2*n[1]);
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
      }
    }else{                  //OUTLET
      if(un-c > 0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	
	_dfdV[0][0] = S * unp;
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	_dfdV[2][0] = S * (VV[2]*unp +cp2*n[1]);
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
	
      }else{//subsonic outlet
	
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	
	
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
      }
    }
  }
  
  
  
  if (localTypeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<dim*dim; ++k)
      jac[k] = dfdV[k]; 
  
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnTaitFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitFDJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(roeflux5waterdissprec)(0, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
          normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(),
          sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
 
}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, bool useLimiter)
{
//void roejacappr3Dwater(int type, double gamma, VarFcn* varFcn, double vfcv, double vfa, 
//                       double vfb, double vfp, FluxFcnBase::Type typeJac, double* normal, 
//                       double normalVel, double* VL, double* VR, double* jacL, 
//                       double* jacR, double irey, double betaRef, double k1, 
//                       double cmach, double shockreducer, double length, int prec, int flag)
  const int dim = 5;
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;
  const int type = 0; //no turbulence

  double dfdUL[dim2], dfdUR[dim2];
  for (int kk = 0; kk<dim2; kk++) dfdUR[kk] = 0.0;
  double n[3] = {normal[0], normal[1], normal[2]};
  F77NAME(roejac5waterdissprec)(type, gamma,
                                vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                                n, normalVel, 
                                VL, VR, dfdUL, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  F77NAME(roejac5waterdissprec)(type, gamma,
                                vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                                n, -normalVel,
                                VR, VL, dfdUR, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
  for (int k=0; k<dim2; k++)
    dfdUR[k] = -dfdUR[k];
  
  int k;

  if (typeJac == FluxFcnBase::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) { 
      jacL[k] = dfdUL[k]; 
      jacR[k] = dfdUR[k]; 
    }
  }
  else {
    vf->postMultiplyBydUdV(VL, dfdUL, jacL);
    vf->postMultiplyBydUdV(VR, dfdUR, jacR);
  }

}

//------------------------------------------------------------------------------

void FluxFcnTaitExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					double *VL, double *VR, double *flux, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitExactJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitExactJacRoeEuler3D::computeJacobians not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  double P = vf->getPrefWater() + vf->getAlphaWater()*pow(V[0], vf->getBetaWater());
  flux[0] = 0.0;
  flux[1] = P * normal[0];
  flux[2] = P * normal[1];
  flux[3] = P * normal[2];
  flux[4] = P * normalVel;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxtait)(0, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{
  fprintf(stderr, "*** Error: FluxFcnTaitInflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnTaitOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux, bool useLimiter)
{
  fprintf(stderr, "*** Error: FluxFcnTaitOutflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, bool useLimiter)
{

  flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacflux3Dwater<5>(0,vf,typeJac,normal,normalVel,V,Ub,jacL);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, bool useLimiter)
{

  flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, bool useLimiter)
{
 
  jacflux3Dwater<5>(0,vf,typeJac,normal,normalVel,V,Ub,jacL);

}

//------------------------------------------------------------------------------

