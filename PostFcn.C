#include <PostFcn.h>

#include <WallFcn.h>
#include <Vector3D.h>

#include <cstdlib>
#include <cstdio>

const double PostFcnEuler::third = 1.0/3.0;
//-----------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	only in the way the constructor creates the PostFcnEuler class
//	and to compute certain values for NS using lame coefficients (lambda and mu)
//------------------------------------------------------------------------------

PostFcn::PostFcn(VarFcn *vf)
{

  varFcn = vf;

}

//------------------------------------------------------------------------------

double PostFcn::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{

  fprintf(stderr, "*** Warning: computeNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcn::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *X, double *dX, double *V, double *dV, double phi)
{

  fprintf(stderr, "*** Warning: computeDerivativeOfNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcn::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					  Vec3D& n, double d2w[3], double* Vwall, 
					  double* Vface[3], double* Vtet[4])
{

  fprintf(stderr, "*** Warning: computeFaceScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

PostFcnEuler::PostFcnEuler(IoData &iod, VarFcn *vf) : PostFcn(vf)
{

  if (iod.eqs.fluidModel.fluid == FluidModelData::GAS){
    mach = iod.ref.mach;
    pinfty = iod.bc.inlet.pressure;

// Included (MB)
    dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

    if (iod.problem.mode == ProblemData::DIMENSIONAL) {
      dimFlag = true;
    }
    else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
      dimFlag = false;
      if (iod.sa.apressFlag == false) {
        if (iod.sa.pressFlag == false) {
          dPin = dpinfty;
        }
        else {
          dPin = 0.0;
        }      
      }
      else {
        dPin = 0.0;
      }
    }

  } else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID){
    mach = iod.ref.mach;
    double P = vf->getPrefWater();
    double a = vf->getAlphaWater();
    double b = vf->getBetaWater();
    pinfty = (P+a*pow(iod.bc.inlet.density, b));
  } else if (iod.eqs.fluidModel.fluid == FluidModelData::JWL){
    mach = iod.ref.mach;
    pinfty = iod.bc.inlet.pressure;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::rstVar(IoData &iod, Communicator *com)
{

  mach = iod.ref.mach;
  pinfty = iod.bc.inlet.pressure;
  dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

  if (iod.problem.mode == ProblemData::DIMENSIONAL) {
    dimFlag = true;
  }
  else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
    dimFlag = false;
    if (iod.sa.apressFlag == false) {
      if (iod.sa.pressFlag == false) {
        dPin = dpinfty;
      }
      else {
        dPin = 0.0;
      }      
    }
    else {
      dPin = 0.0;
    }
  }
  
}

//------------------------------------------------------------------------------

double PostFcnEuler::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{
  double q = 0.0;
  double n[3];

  if (type == DENSITY)
    q = varFcn->getDensity(V, fluidId);
  else if (type == MACH)
    q = varFcn->computeMachNumber(V, fluidId);
  else if (type == WTMACH)
    q = varFcn->computeWtMachNumber(V, fluidId);
  else if (type == SPEED)
    q = sqrt(varFcn->computeU2(V));
  else if (type == WTSPEED)
    q = sqrt(varFcn->computeWtU2(V));
  else if (type == PRESSURE)
    q = varFcn->getPressure(V, fluidId);
  else if (type == DIFFPRESSURE)
    q = varFcn->getPressure(V, fluidId)-pinfty;
  else if (type == TEMPERATURE)
    q = varFcn->computeTemperature(V, fluidId);
  else if (type == TOTPRESSURE)
    q = varFcn->computeTotalPressure(mach, V, fluidId);
  else if (type == NUT_TURB)
    q = varFcn->getTurbulentNuTilde(V, fluidId);
  else if (type == K_TURB)
    q = varFcn->getTurbulentKineticEnergy(V, fluidId);
  else if (type == EPS_TURB)
    q = varFcn->getTurbulentDissipationRate(V,fluidId);
  else if(type == HYDROSTATICPRESSURE)
    q = varFcn->hydrostaticPressure(V[0],X);
  else if(type == HYDRODYNAMICPRESSURE)
    q = varFcn->hydrodynamicPressure(V,X,fluidId);
  else if(type == PRESSURECOEFFICIENT)
    q = varFcn->computePressureCoefficient(V, pinfty, mach, dimFlag,fluidId);
  else if(type == PHILEVEL)
    //q = static_cast<double>(fluidId);
    q = phi[0]/varFcn->getDensity(V, fluidId);
  else if (type == FLUIDID)
    q = static_cast<double>(fluidId);
 // Included (MB)
  else if (type == VELOCITY_NORM)
    q = varFcn->getVelocityNorm(V,fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double q = 0.0;

  if (type == DERIVATIVE_DENSITY)
    q = varFcn->getDensity(dV);
  else if (type == DERIVATIVE_MACH)
    q = varFcn->computeDerivativeOfMachNumber(V, dV, dS[0]);
  else if (type == DERIVATIVE_PRESSURE)
    q = varFcn->getPressure(dV);
  else if (type == DERIVATIVE_TEMPERATURE)
    q = varFcn->computeDerivativeOfTemperature(V, dV);
  else if (type == DERIVATIVE_TOTPRESSURE)
    q = varFcn->computeDerivativeOfTotalPressure(mach, dS[0], V, dV, dS[0]);
  else if (type == DERIVATIVE_NUT_TURB)
    q = varFcn->getTurbulentNuTilde(dV);
  else if (type == DERIVATIVE_VELOCITY_SCALAR)
    q = varFcn->getDerivativeOfVelocityNorm(V, dV);
  else
  {
    // Error message
    fprintf(stderr, "*** Warning: PostFcnEuler::computeDerivativeOfNodeScalarQuantity does not define the type %d\n", type);
  }

  return q;

}

//------------------------------------------------------------------------------

void PostFcnEuler::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                                double *Vwall, double *Vface[3], double *Vtet[4],
				double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro,
				int fid)
{

  double pcg[3], p[3];
  double pcgin;
  int i;

  switch(hydro)
    {
    case 0:
      for(i=0;i<3;i++) pcg[i] = varFcn->getPressure(Vface[i],fid);
      break;
    case 1: // hydrostatic pressure
      for(i=0;i<3;i++) pcg[i] = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
      break;
    case 2: // hydrodynamic pressure
      for(i=0;i<3;i++) pcg[i] = varFcn->hydrodynamicPressure(Vface[i],Xface[i]);
      break;
    default:
      fprintf(stderr,"hydro parameter is not correct. Pressure at the face cannot be computed. hydro = %d\n",hydro);
      exit(-1);
    }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

// ##################
 Vec3D x0, x1, x2, x;
 double temp;
 for(i = 0; i<3; i++)
 {
  x0[i] = Xface[0][i];
  x1[i] = Xface[1][i];
  x2[i] = Xface[2][i];
 }

// for node 0
 i=0; // i represents node
 x = (7.0/18.0)*(x1 + x2 - 2.0*x0);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi0 = (1.0/6.0*temp)*n;

// for node 1
 i=1;
 x = (7.0/18.0)*(x2 + x0 - 2.0*x1);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi1 = (1.0/6.0*temp)*n;

// for node 2
 i=2;
 x = (7.0/18.0)*(x0 + x1 - 2.0*x2);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi2 = (1.0/6.0*temp)*n;

 Fv = 0.0;

}
//------------------------------------------------------------------------------

void PostFcnEuler::computeForceEmbedded(int orderOfAccuracy, double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
					double *Vwall, double *Vface[3], double *Vtet[4],
					double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  Vec3D p;
  double pcgin;
  int i;

  switch(hydro)
    {
    case 0:
      for(i=0;i<3;i++) p[i] = varFcn->getPressure(Vface[i],fid);
      break;
    case 1: // hydrostatic pressure
      for(i=0;i<3;i++) p[i] = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
      break;
    case 2: // hydrodynamic pressure
      for(i=0;i<3;i++) p[i] = varFcn->hydrodynamicPressure(Vface[i],Xface[i]);
      break;
    default:
      fprintf(stderr,"hydro parameter is not correct. Pressure at the face cannot be computed. hydro = %d\n",hydro);
      exit(-1);
    }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p -= pcgin;

 // Computes Int_{T} (N_i * Sum_{j\in T} p_j N_j) dx
 // At first order, N_j = Chi_j (constant per control volume)
 // At second order, N_j = Phi_j (P1 lagrangian basis functions within T)
 switch (orderOfAccuracy)
   {
   case 1: 
     Fi0 = (p[0]/3.0)*n;
     Fi1 = (p[1]/3.0)*n;
     Fi2 = (p[2]/3.0)*n;
     break;
   case 2:
     double c1 = 1.0/6.0, c2 = 1.0/12.0;
     Fi0 = (c1*p[0]+c2*(p[1]+p[2]))*n;
     Fi1 = (c1*p[1]+c2*(p[2]+p[0]))*n;
     Fi2 = (c1*p[2]+c2*(p[0]+p[1]))*n;
     break;
   }     
     
 Fv = 0.0;
}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
					    double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],  
					    double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  double pcg[3], p[3];
  double dPcg[3], dP[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]);
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++) {
        pcg[i]  = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
	dPcg[i] = varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
     }
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i]  = varFcn->hydrodynamicPressure(Vface[i], Xface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]) 
              - varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  if (pin) {
    if (dimFlag)
      dPcgin = (-2.0*(*pin)/mach) * dS[0];
    else
      dPcgin = dPin * dS[0];
  }
  else {
    if (hydro == 0)
      dPcgin = dPinfty;
    else
      dPcgin = 0.0;
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  dP[0] = (dPcg[0] - dPcgin) ;
  dP[1] = (dPcg[1] - dPcgin) ;
  dP[2] = (dPcg[2] - dPcgin) ;

 Vec3D x0, x1, x2, x;
 Vec3D dx0, dx1, dx2, dx;
 double temp, dTemp;
 for(int i = 0; i<3; i++)
 {
  x0[i] = Xface[0][i];
  x1[i] = Xface[1][i];
  x2[i] = Xface[2][i];
  dx0[i] = dXface[0][i];
  dx1[i] = dXface[1][i];
  dx2[i] = dXface[2][i];
 }

// for node 0
 i=0; // i represents node
 x = (7.0/18.0)*(x1 + x2 - 2.0*x0);
 dx = (7.0/18.0)*(dx1 + dx2 - 2.0*dx0);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi0 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

// for node 1
 i=1;
 x = (7.0/18.0)*(x2 + x0 - 2.0*x1);
 dx = (7.0/18.0)*(dx2 + dx0 - 2.0*dx1);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi1 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

// for node 2
 i=2;
 x = (7.0/18.0)*(x0 + x1 - 2.0*x2);
 dx = (7.0/18.0)*(dx0 + dx1 - 2.0*dx2);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi2 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

  dFv = 0.0;

}

//-----------------------------------------------------------------------------------

void PostFcnEuler::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
					   double *Vwall, double *Vface[3], double *Vtet[4],
					   double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro,int fid)
{

  double pcg[3], p[3];
  double pcgin;
  int i;

  if (hydro == 0) {
    for(i=0;i<3;i++)
      pcg[i] = varFcn->getPressure(Vface[i],fid);
  }
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++)
        pcg[i] = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++)
      pcg[i] = varFcn->hydrodynamicPressure(Vface[i],Xface[i]);

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

// ##########################

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 // Computing p_1C, p_0C and p_3C
 i = 0;
 x = x1-xC;
 p_1C = p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xC;
 p_0C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xC;
 p_3C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_1S, p_0S and p_2S
 i = 1;
 x = x1-xS;
 p_1S = p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xS;
 p_0S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x2-xS;
 p_2S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_2P, p_0P and p_3P
 i = 2;
 x = x2-xP;
 p_2P = p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xP;
 p_0P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xP;
 p_3P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0 ;
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0 ;
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0 ;

 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0 ;
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0 ;
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0 ;

 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0 ;

 Fi0 = (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi1 = (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi2 = (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);

 Fv = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::computeDerivativeOfForceTransmitted(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
					    double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],  
					    double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  double pcg[3], p[3];
  double dPcg[3], dP[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]);
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++) {
        pcg[i]  = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
	dPcg[i] = varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
     }
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i]  = varFcn->hydrodynamicPressure(Vface[i], Xface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]) 
              - varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  if (pin) {
    if (dimFlag)
      dPcgin = (-2.0*(*pin)/mach) * dS[0];
    else
      dPcgin = dPin * dS[0];
  }
  else {
    if (hydro == 0)
      dPcgin = dPinfty;
    else
      dPcgin = 0.0;
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  dP[0] = (dPcg[0] - dPcgin) ;
  dP[1] = (dPcg[1] - dPcgin) ;
  dP[2] = (dPcg[2] - dPcgin) ;

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 Vec3D dxC, dxS, dxP, dx1, dx2, dx3, dx0, dx;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double dP_1C, dP_1S, dP_2S, dP_2P, dP_3P, dP_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;
 double dP_C, dP_S, dP_P, dP_0C, dP_0S, dP_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
  dxC[i] = dXface[0][i];
  dxS[i] = dXface[1][i];
  dxP[i] = dXface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;
 dx1 = (dxC+dxS)/2.0;
 dx2 = (dxS+dxP)/2.0;
 dx3 = (dxP+dxC)/2.0;
 dx0 = (dxC+dxS+dxP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 dP_C = dP[0], dP_S = dP[1], dP_P = dP[2];
 // Computing p_1C, p_0C and p_3C, and computing dP_1C, dP_0C and dP_3C
 i = 0;
 x = x1-xC;
 dx = dx1-dxC;
 p_1C = p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_1C = dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xC;
 dx = dx0-dxC;
 p_0C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0C =  dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x3-xC;
 dx = dx3-dxC;
 p_3C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_3C =  dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // Computing p_1S, p_0S and p_2S, and computing dP_1S, dP_0S and dP_2S
 i = 1;
 x = x1-xS;
 dx = dx1-dxS;
 p_1S = p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_1S = dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xS;
 dx = dx0-dxS;
 p_0S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0S =  dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x2-xS;
 dx = dx2-dxS;
 p_2S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_2S =  dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // Computing p_2P, p_0P and p_3P, and computing dP_2P, dP_0P and dP_3P
 i = 2;
 x = x2-xP;
 dx = dx2-dxP;
 p_2P = p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_2P = dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xP;
 dx = dx0-dxP;
 p_0P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0P =  dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x3-xP;
 dx = dx3-dxP;
 p_3P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_3P =  dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 double dPhi_C, dPhi_S, dPhi_P, dPhi_1, dPhi_2, dPhi_3, dPhi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0 ;
 dPhi_C = (2.0*dP_C + dP_0C + dP_3C)/12.0 + (2.0*dP_C + dP_0C + dP_1C)/12.0 ;
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0 ;
 dPhi_S = (2.0*dP_S + dP_0S + dP_1S)/12.0 + (2.0*dP_S + dP_0S + dP_2S)/12.0 ;
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0 ;
 dPhi_P = (2.0*dP_P + dP_0P + dP_2P)/12.0 + (2.0*dP_P + dP_0P + dP_3P)/12.0 ;

 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0 ;
 dPhi_1 = (2.0*dP_1C + dP_0C + dP_C)/12.0 + (2.0*dP_1S + dP_0S + dP_S)/12.0 ;
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0 ;
 dPhi_2 = (2.0*dP_2S + dP_0S + dP_S)/12.0 + (2.0*dP_2P + dP_0P + dP_P)/12.0 ;
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0 ;
 dPhi_3 = (2.0*dP_3P + dP_0P + dP_P)/12.0 + (2.0*dP_3C + dP_0C + dP_C)/12.0 ;

 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0 ;
 dPhi_0 = (2.0*dP_0C + dP_1C + dP_C)/12.0 + (2.0*dP_0C + dP_3C + dP_C)/12.0 + (2.0*dP_0S + dP_1S + dP_S)/12.0 + (2.0*dP_0S + dP_2S + dP_S)/12.0 + (2.0*dP_0P + dP_2P + dP_P)/12.0 + (2.0*dP_0P + dP_3P + dP_P)/12.0 ;

 dFi0 = (dPhi_C + dPhi_1/2.0 + dPhi_3/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*dn);
 dFi1 = (dPhi_S + dPhi_1/2.0 + dPhi_2/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0) * (1.0/6.0*dn);
 dFi2 = (dPhi_P + dPhi_2/2.0 + dPhi_3/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*dn);

 dFv = 0.0;

}

//--------------------------------------------------------------------------------------

double PostFcnEuler::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				      double* Vwall, double* Vface[3], double* Vtet[4])
{

  return 0.0;

}

//--------------------------------------------------------------------------------------

double PostFcnEuler::computeHeatFluxRelatedValues(double dp1dxj[4][3], Vec3D& n, double d2w[3],
                                               double* Vwall, double* Vface[3], double* Vtet[4], bool includeKappa)
{

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				      double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcnEuler::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
					  double d2w[3], double* Vwall, double* Vface[3], 
					  double* Vtet[4], double pin)
{
  double p = third * ( varFcn->getPressure(Vface[0]) + varFcn->getPressure(Vface[1]) +
		       varFcn->getPressure(Vface[2]) ) - pin;
  double W = - ndot * p;

  return W;
}

//------------------------------------------------------------------------------

PostFcnNS::PostFcnNS(IoData &iod, VarFcn *vf) 
  : PostFcnEuler(iod, vf), NavierStokesTerm(iod, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, PostFcn::varFcn, viscoFcn);
  else
    wallFcn = 0;

}

//------------------------------------------------------------------------------

PostFcnNS::~PostFcnNS()
{
  if (wallFcn) delete wallFcn;
}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnNS::rstVar(IoData &iod, Communicator *com)
{

  PostFcnEuler::rstVar(iod, com);

  NavierStokesTerm::rstVarNS(iod, com);
  
  if (wallFcn)
    wallFcn->rstVar(iod, com);

}

//------------------------------------------------------------------------------

double PostFcnNS::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					    Vec3D& n, double d2w[3], double* Vwall, 
					    double* Vface[3], double* Vtet[4])
{

  double q = 0.0;

  if (type == DELTA_PLUS) {
#if defined(HEAT_FLUX)
    q = computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet) / sqrt(n*n);
#else
    if (wallFcn)
      q = wallFcn->computeDeltaPlus(n, d2w, Vwall, Vface);
    else
      fprintf(stderr, "*** Warning: yplus computation not implemented\n");
#endif
  }
  else if (type == SKIN_FRICTION) {
    Vec3D t(1.0, 0.0, 0.0);
    Vec3D F = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
    q = 2.0 * t * F / sqrt(n*n);
  }

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double q = 0.0;

  q = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return q;

}

//------------------------------------------------------------------------------

Vec3D PostFcnNS::computeViscousForce(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				     double* Vwall, double* Vface[3], double* Vtet[4])
{

  Vec3D Fv;

  if (wallFcn && Vwall)
    Fv = wallFcn->computeForce(n, d2w, Vwall, Vface);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    Fv[0] = tij[0][0] * n[0] + tij[0][1] * n[1] + tij[0][2] * n[2];
    Fv[1] = tij[1][0] * n[0] + tij[1][1] * n[1] + tij[1][2] * n[2];
    Fv[2] = tij[2][0] * n[0] + tij[2][1] * n[1] + tij[2][2] * n[2];
  }

  return -1.0 * Fv;

}

//------------------------------------------------------------------------------

Vec3D PostFcnNS::computeViscousForceCVBoundary(Vec3D& n,  double* Vi, double dudxj[3][3])
{

  Vec3D Fv;
  // Could be useful later…
  /*
  if (wallFcn)
    Fv = wallFcn->computeForce(n, d2w, Vwall, Vface);
  else {
  */
  double T;
  computeTemperature(Vi,T);
  
  double mu = ooreynolds_mu * viscoFcn->compute_mu(T);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(T, mu);
  
  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);
  
  Fv[0] = tij[0][0] * n[0] + tij[0][1] * n[1] + tij[0][2] * n[2];
  Fv[1] = tij[1][0] * n[0] + tij[1][1] * n[1] + tij[1][2] * n[2];
  Fv[2] = tij[2][0] * n[0] + tij[2][1] * n[1] + tij[2][2] * n[2];
  
  return -1.0 * Fv;
}

//------------------------------------------------------------------------------

// Included (MB)
Vec3D PostFcnNS::computeDerivativeOfViscousForce(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3],
				     double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  Vec3D dFv;

  if (wallFcn)
    dFv = wallFcn->computeDerivativeOfForce(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double du[4][3], ducg[3];
    computeDerivativeOfVelocity(dVtet, du, ducg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double ddudxj[3][3];
    computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];

    double dooreynolds_lambda = -1.0 / ( reynolds_lambdaNS * reynolds_lambdaNS ) * dRe_lambdadMachNS * dS[0];

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);

    double dmu = dooreynolds_mu * viscoFcn->compute_mu(Tcg) + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dS[0]);

    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double dlambda = dooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu) + ooreynolds_lambda * viscoFcn->compute_lambdaDerivative(mu, dmu, dS[0]);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    double dtij[3][3];
    computeDerivativeOfStressTensor(mu, dmu, lambda, dlambda, dudxj, ddudxj, dtij);

    dFv[0] = dtij[0][0] * n[0] + dtij[0][1] * n[1] + dtij[0][2] * n[2] + tij[0][0] * dn[0] + tij[0][1] * dn[1] + tij[0][2] * dn[2];
    dFv[1] = dtij[1][0] * n[0] + dtij[1][1] * n[1] + dtij[1][2] * n[2] + tij[1][0] * dn[0] + tij[1][1] * dn[1] + tij[1][2] * dn[2];
    dFv[2] = dtij[2][0] * n[0] + dtij[2][1] * n[1] + dtij[2][2] * n[2] + tij[2][0] * dn[0] + tij[2][1] * dn[1] + tij[2][2] * dn[2];

  }

  return -1.0 * dFv;

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
			     double *Vwall, double *Vface[3], double *Vtet[4], 
                    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  PostFcnEuler::computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro,fid);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForceEmbedded(int orderOfAccuracy,double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
				     double *Vwall, double *Vface[3], double *Vtet[4], 
				     double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  PostFcnEuler::computeForceEmbedded(orderOfAccuracy,dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro,fid);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void PostFcnNS::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                            double *Vface[3], double *dVface[3], double *Vtet[4],
                                            double *dVtet[4], double dS[3], double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2,
                                            Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  PostFcnEuler::computeDerivativeOfForce(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx,  ddPdx, hydro);

  dFv = computeDerivativeOfViscousForce(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                             double *Vwall, double *Vface[3], double *Vtet[4],
                    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  PostFcnEuler::computeForceTransmitted(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro,fid);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void PostFcnNS::computeDerivativeOfForceTransmitted(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                            double *Vface[3], double *dVface[3], double *Vtet[4],
                                            double *dVtet[4], double dS[3], double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2,
                                            Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  PostFcnEuler::computeDerivativeOfForceTransmitted(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx,  ddPdx, hydro);

  dFv = computeDerivativeOfViscousForce(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

}

//------------------------------------------------------------------------------

double PostFcnNS::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				   double* Vwall, double* Vface[3], double* Vtet[4])
{

  double hp = 0.0;

  if (wallFcn)
    hp = wallFcn->computeHeatPower(n, d2w, Vwall, Vface);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);
    double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    hp = qj[0]*n[0] + qj[1]*n[1] + qj[2]*n[2]; 
}

  return hp;

}
//------------------------------------------------------------------------------
double PostFcnNS::computeHeatFluxRelatedValues(double dp1dxj[4][3], Vec3D& n, double d2w[3],
                                   double* Vwall, double* Vface[3], double* Vtet[4], bool includeKappa)
{
  double hp = 0.0;

  if (wallFcn)
    hp = wallFcn->computeHeatPower(n, d2w, Vwall, Vface);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);  

    double kappa = -1; //The fact that it is negative balances the minus sign in NavierStokesTerm::computeHeatFluxVector
    if(includeKappa == true)
      {
        kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);
      }
     
       double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    hp = qj[0]*n[0] + qj[1]*n[1] + qj[2]*n[2]; 
    }
  return hp;

}


//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				   double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  double dhp = 0.0;

  if (wallFcn)
    dhp = wallFcn->computeDerivativeOfHeatPower(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double ddTdxj[3];
    computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);
    double kappa = ooreynolds * thermalCondFcn->compute(Tcg);
    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];
    double dkappa = dooreynolds_mu * thermalCondFcn->compute(Tcg) + ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcg, dS[0]);
    double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    double dqj[3];
    computeDerivativeOfHeatFluxVector(kappa, dkappa, dTdxj, ddTdxj, dqj);
    dhp = dqj[0]*n[0] + qj[0]*dn[0] + dqj[1]*n[1] + qj[1]*dn[1] + dqj[2]*n[2] + qj[2]*dn[2]; 
  }

  return dhp;

}

//------------------------------------------------------------------------------

double PostFcnNS::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
				       double d2w[3], double* Vwall, double* Vface[3], 
				       double* Vtet[4], double pin)
{

  double W = PostFcnEuler::computeInterfaceWork(dp1dxj, n, ndot, d2w, Vwall, Vface, Vtet, pin);

  if (wallFcn)
    W += wallFcn->computeInterfaceWork(n, d2w, Vwall, Vface);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
    double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    W += (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * n[0] +
      (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * n[1] +
      (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * n[2];
  }

  return W;

}

//------------------------------------------------------------------------------

PostFcnSA::PostFcnSA(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), SATerm(iod)
{
}

//------------------------------------------------------------------------------

PostFcnSA::~PostFcnSA()
{
}
//------------------------------------------------------------------------------

// Included (MB)
void PostFcnSA::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarSA(iod);

}

//------------------------------------------------------------------------------

double PostFcnSA::computeNodeScalarQuantity(ScalarType type, double *V, double *X,  int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnSA::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
  {
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity
         (type, dS, V, dV, X, dX);
  }

  return dq;

}

// Included (MB)
void PostFcnDES::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarDES(iod);

}

//------------------------------------------------------------------------------
PostFcnDES::PostFcnDES(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), DESTerm(iod)
{

}

//------------------------------------------------------------------------------

                                                                                                
double PostFcnDES::computeNodeScalarQuantity(ScalarType type, double *V, double *X,int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnDES::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------

PostFcnKE::PostFcnKE(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), KEpsilonTerm(iod)
{

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnKE::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarKE(iod);

}

//------------------------------------------------------------------------------

double PostFcnKE::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY)
    q = computeTurbulentViscosity(V);
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnKE::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY)
    dq = computeDerivativeOfTurbulentViscosity(V, dV, dS[0]);
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------
