#include <AgglomeratedFace.h>

#include "NavierStokesTerm.h"

AgglomeratedFace::AgglomeratedFace() : node(0), code(0), normal(0.0) {

}

AgglomeratedFace::AgglomeratedFace(int node, int code) : node(node), code(code), normal(0.0) {

}

AgglomeratedFace::AgglomeratedFace(const AgglomeratedFace& oth) :
  node(oth.node), code(oth.code), normal(oth.normal) {

}

AgglomeratedFace::~AgglomeratedFace() {

}

template<int dim>
void AgglomeratedFace::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U)
{
  int k, j;

  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED)
    for (k=0; k<dim; ++k) {
      U[k] = Uin[node][k];
    }
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k) {
      U[k] = Uout[node][k];
    }
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;
}


template<int dim>
void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
	                                       SVec<double,dim> &V, 
			                       double *Ub, SVec<double,dim> &fluxes) {

  if(fluxFcn[code]){
    double flux[dim];
    fluxFcn[code]->compute(0.0, 0.0,  normal, 0.0, 
  		           V[node], Ub, flux);
    for (int k=0; k<dim; ++k){
      fluxes[ node ][k] += flux[k];
    } 
  }
}

void AgglomeratedFace::setNodeType(int* priority, int* nodeType) {
 
  if (priority[code] > priority[ nodeType[ node ] ])
    nodeType[ node ] = code;
}

AgglomeratedFaceSet::AgglomeratedFaceSet(int size) : numFaces(size) {

  myFaces = new AgglomeratedFace[numFaces];
}

AgglomeratedFaceSet::~AgglomeratedFaceSet() {

  delete [] myFaces;
}

template<int dim, class Scalar, int neq>
void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, 
			 	                       SVec<double,dim> &V, 
 				                       double *Ub, GenMat<Scalar,neq> &A) {

  double jac[neq*neq];

  double normVel= 0.0;

  fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[node], Ub, jac);
  Scalar *Aii = A.getElem_ii(node);
  for (int k=0; k<neq*neq; ++k) 
    Aii[k] += jac[k];

}

template<int dim>
void AgglomeratedFace::
computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm* ns,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,dim> &fluxes) {

  if (code == BC_SYMMETRY) return;

  double ooreynolds_mu = ns->get_ooreynolds_mu(); 
  double area = normal.norm();
  
  Vec3D xhat = normal;
  xhat /= area;

  double Tcg = varFcn->computeTemperature(V[node]);

  double Tg[5];
  varFcn->computeTemperatureGradient(V[node],Tg);

  double mu     = ns->getViscoFcn()->compute_mu(Tcg);
  double lambda = ns->getViscoFcn()->compute_lambda(Tcg,mu);
  double kappa  = ns->getThermalCondFcn()->compute(Tcg);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double uhat = V[node][1]*xhat[0]+V[node][2]*xhat[1]+V[node][3]*xhat[2];

  // Not necessarily p
  double dpdxhat = dX[node][4]*xhat[0]+dY[node][4]*xhat[1]+dZ[node][4]*xhat[2];
  // Not necessarily rho
  double drhodxhat = dX[node][0]*xhat[0]+dY[node][0]*xhat[1]+dZ[node][0]*xhat[2];

  double dudxhat[3] = {dX[node][1]*xhat[0]+dY[node][1]*xhat[1]+dZ[node][1]*xhat[2],
                       dX[node][2]*xhat[0]+dY[node][2]*xhat[1]+dZ[node][2]*xhat[2],
                       dX[node][3]*xhat[0]+dY[node][3]*xhat[1]+dZ[node][3]*xhat[2] };

  double duhatdx = dX[node][1]*xhat[0]+dX[node][2]*xhat[1]+dX[node][3]*xhat[2];
  double duhatdy = dY[node][1]*xhat[0]+dY[node][2]*xhat[1]+dY[node][3]*xhat[2]; 
  double duhatdz = dZ[node][1]*xhat[0]+dZ[node][2]*xhat[1]+dZ[node][3]*xhat[2]; 
  
  double duhatdxhat = duhatdx*xhat[0]+duhatdy*xhat[1]+duhatdz*xhat[2];

  double dTdxhat = drhodxhat*Tg[0]+dudxhat[0]*Tg[1]+dudxhat[1]*Tg[2]+
                   dudxhat[2]*Tg[3]+dpdxhat*Tg[4];
 
  //double velnorm = sqrt(V[node][1]*V[node][1]+V[node][2]*V[node][2]+V[node][3]*V[node][3]);

  double flux[5];
  double f1 = (lambda+2.0*mu)*duhatdxhat;
  flux[0] = 0.0;
  for (int k = 0; k < 3; ++k) {

    flux[k+1] = xhat[k]*f1+mu*(dudxhat[k]-duhatdxhat*xhat[k]);
  }
  
  flux[4] = uhat*f1+mu*(V[node][1]*dudxhat[0]+V[node][2]*dudxhat[1]+V[node][3]*dudxhat[2] -
                        uhat*duhatdxhat)-kappa*dTdxhat ;

  for (int k = 0; k < 5; ++k) {

    fluxes[node][k] += flux[k]*area;
  }
}

template<int dim,class Scalar,int neq>
void AgglomeratedFace::
computeJacobianThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm* ns,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,neq*neq> &JacX,
                                        SVec<double,neq*neq> &JacY,
                                        SVec<double,neq*neq> &JacZ,
                                        Vec<double>& ctrlVol,
                                        GenMat<Scalar,neq>& A) {

  if (code == BC_SYMMETRY) return;

  Scalar* jac = A.getElem_ii(node);
  
  Scalar* jacx = JacX[node];
  Scalar* jacy = JacY[node];
  Scalar* jacz = JacZ[node];

  //memset(jac,0,sizeof(double)*neq*neq);
  memset(jacx,0,sizeof(Scalar)*neq*neq);
  memset(jacy,0,sizeof(Scalar)*neq*neq);
  memset(jacz,0,sizeof(Scalar)*neq*neq);

  double ooreynolds_mu = ns->get_ooreynolds_mu(); 
  double area = normal.norm();
  
  Vec3D xhat = normal;
  xhat /= area;

  double Tcg = varFcn->computeTemperature(V[node]);

  double Trr,Trp,Tpp;

  double Tg[5];
  varFcn->computeTemperatureGradient(V[node],Tg);
  varFcn->computeTemperatureHessian(V[node],Trr,Trp,Tpp);

  double mu     = ns->getViscoFcn()->compute_mu(Tcg);
  double lambda = ns->getViscoFcn()->compute_lambda(Tcg,mu);
  double kappa  = ns->getThermalCondFcn()->compute(Tcg);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double uhat = V[node][1]*xhat[0]+V[node][2]*xhat[1]+V[node][3]*xhat[2];

  // Not necessarily p
  double dpdxhat = dX[node][4]*xhat[0]+dY[node][4]*xhat[1]+dZ[node][4]*xhat[2];
  // Not necessarily rho
  double drhodxhat = dX[node][0]*xhat[0]+dY[node][0]*xhat[1]+dZ[node][0]*xhat[2];

  double dudxhat[3] = {dX[node][1]*xhat[0]+dY[node][1]*xhat[1]+dZ[node][1]*xhat[2],
                       dX[node][2]*xhat[0]+dY[node][2]*xhat[1]+dZ[node][2]*xhat[2],
                       dX[node][3]*xhat[0]+dY[node][3]*xhat[1]+dZ[node][3]*xhat[2] };

  double duhatdx = dX[node][1]*xhat[0]+dX[node][2]*xhat[1]+dX[node][3]*xhat[2];
  double duhatdy = dY[node][1]*xhat[0]+dY[node][2]*xhat[1]+dY[node][3]*xhat[2]; 
  double duhatdz = dZ[node][1]*xhat[0]+dZ[node][2]*xhat[1]+dZ[node][3]*xhat[2]; 
  
  double duhatdxhat = duhatdx*xhat[0]+duhatdy*xhat[1]+duhatdz*xhat[2];

  double dTdxhat = drhodxhat*Tg[0]+dudxhat[0]*Tg[1]+dudxhat[1]*Tg[2]+
                   dudxhat[2]*Tg[3]+dpdxhat*Tg[4];
 
  //double velnorm = sqrt(V[node][1]*V[node][1]+V[node][2]*V[node][2]+V[node][3]*V[node][3]);

  double flux[5];
  double f1 = (lambda+2.0*mu)*duhatdxhat;
  flux[0] = 0.0;
  double xhatsum = xhat[0]+xhat[1]+xhat[2];
  for (int k = 0; k < 3; ++k) {

    for (int l = 0; l < 3; ++l) {
      jacx[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[0];
      jacy[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[1];
      jacz[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[2];
  
      jacx[(k+1)*neq+l+1] += mu*((k==l?xhat[0]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[0]*xhat[k]);
      jacy[(k+1)*neq+l+1] += mu*((k==l?xhat[1]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[1]*xhat[k]);
      jacz[(k+1)*neq+l+1] += mu*((k==l?xhat[2]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[2]*xhat[k]);
    } 
    //flux[k+1] = xhat[k]*f1+mu*(dudxhat[k]-duhatdxhat*xhat[k]);
  }

  for (int l = 0; l < 3; ++l) {
    jacx[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[0])*uhat;
    jacy[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[1])*uhat;
    jacz[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[2])*uhat;

    jac[4*neq+l+1] += xhat[l]*f1*area;

    jac[4*neq+l+1] += mu*(dudxhat[l]-xhat[l]*duhatdxhat)*area;

    jacx[4*neq+l+1] += mu*(V[node][1]*(l==0?xhat[0]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[0]);
    jacy[4*neq+l+1] += mu*(V[node][2]*(l==1?xhat[1]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[1]);
    jacz[4*neq+l+1] += mu*(V[node][3]*(l==2?xhat[2]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[2]);

    jacx[4*neq] = -kappa*Tg[0]*xhat[0];
    jacx[4*neq+4] = -kappa*Tg[4]*xhat[0];
    jacy[4*neq] = -kappa*Tg[0]*xhat[1];
    jacy[4*neq+4] = -kappa*Tg[4]*xhat[1];
    jacz[4*neq] = -kappa*Tg[0]*xhat[2];
    jacz[4*neq+4] = -kappa*Tg[4]*xhat[2];

    jac[4*neq] += -kappa*(drhodxhat*Trr+dpdxhat*Trp)*area;
    jac[4*neq+4] += -kappa*(drhodxhat*Trp+dpdxhat*Tpp)*area;
  } 
   
  //flux[4] = uhat*f1+mu*(V[node][1]*dudxhat[0]+V[node][2]*dudxhat[1]+V[node][3]*dudxhat[2] -
  //                      uhat*duhatdxhat)-kappa*dTdxhat ;

  for (int k = 0; k < neq*neq; ++k) {

    jacx[k] *= area;
    jacy[k] *= area;
    jacz[k] *= area;
  }

  double tmpx[neq*neq],tmpy[neq*neq],tmpz[neq*neq];
  varFcn->postMultiplyBydVdU(V[node], jacx,tmpx);
  varFcn->postMultiplyBydVdU(V[node], jacy,tmpy);
  varFcn->postMultiplyBydVdU(V[node], jacz,tmpz);
  
  double invvol = 1.0/ctrlVol[node];
  for (int k = 0; k < neq*neq; ++k) {

    jac[k] += invvol*(normal[0]*jacx[k]+normal[1]*jacy[k]+normal[2]*jacz[k]);
  }
}

template<int dim>
void AgglomeratedFaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
                                               SVec<double,dim> &V, 
			                       SVec<double,dim>& Ub, SVec<double,dim> &fluxes) {

  for (int i = 0; i < numFaces; ++i) {

    myFaces[i].computeFiniteVolumeTerm(fluxFcn, V, Ub[i], fluxes);
  }
}

template<int dim, class Scalar, int neq>
void AgglomeratedFaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn,
			                               SVec<double,dim> &V, 
				                       SVec<double,dim>& Ub, GenMat<Scalar,neq> &A) {

  for (int i = 0; i < numFaces; ++i) {
 
    myFaces[i].computeJacobianFiniteVolumeTerm(fluxFcn,V, Ub[i], A);
  }
}

template<int dim>
void AgglomeratedFaceSet::
computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm* ns,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,dim> &fluxes) {

  for (int i = 0; i < numFaces; ++i) {
 
    myFaces[i].computeThinLayerViscousFiniteVolumeTerm(ns,varFcn, V,dX,dY, dZ,fluxes);
  }
}


#define INST_HELPER(dim) \
template void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
         SVec<double,dim> &V, double* Ub, SVec<double,dim> &fluxes); \
template void AgglomeratedFaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
         SVec<double,dim> &V, SVec<double,dim>& Ub, SVec<double,dim> &fluxes); \
template void AgglomeratedFace::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U); \
  template void  AgglomeratedFace::computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm*, \
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ, \
			                SVec<double,dim> &fluxes);\
  template void  AgglomeratedFaceSet::computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm*, \
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ, \
			                SVec<double,dim> &fluxes);

#define INST_HELPER2(dim,Scalar,neq) \
template void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
                         SVec<double,dim> &V, double *Ub, GenMat<Scalar,neq> &A); \
template  void AgglomeratedFace:: \
  computeJacobianThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm* ns,\
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ,\
			                SVec<double,neq*neq> &JacX, \
                                        SVec<double,neq*neq> &JacY, \
                                        SVec<double,neq*neq> &JacZ, \
                                        Vec<double>& ctrlVol, \
                                        GenMat<Scalar,neq>& A);\
template void AgglomeratedFaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
                  SVec<double,dim> &V, SVec<double,dim>& Ub, GenMat<Scalar,neq> &A);


INST_HELPER(1);
INST_HELPER(2);
INST_HELPER(3);
INST_HELPER(5);
INST_HELPER(6);
INST_HELPER(7);

INST_HELPER2(1,double,1);
INST_HELPER2(2,double,2);
INST_HELPER2(5,double,5);
INST_HELPER2(6,double,6);
INST_HELPER2(7,double,7);
