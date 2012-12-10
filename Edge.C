#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

#include <Edge.h>
#include <BcDef.h>
#include <FluxFcn.h>
#include <RecFcn.h>
#include <NodalGrad.h>
#include <EdgeGrad.h>
#include <Elem.h>
#include <ExactRiemannSolver.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <LowMachPrec.h>
#include "LevelSet/LevelSetStructure.h"
#include "FluidSelector.h"
#include "DenseMatrixOps.h"
#include "FemEquationTermDesc.h"

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                              SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, Vec<double> &idtv,
                              TimeLowMachPrec &tprec)
{

  double Vmid[dim];
  double Xmid[3];
  double vis = 0.0;
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  double locbeta=0.0;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double S = sqrt(normal[l] * normal[l]);
    double invS = 1.0 / S;

    Vec3D n = invS * normal[l];
    double ndot = invS * normalVel[l];

    for (int k=0; k<dim; ++k)
      Vmid[k] = 0.5 * (V[i][k] + V[j][k]);
    for (int k=0; k<3; k++)
      Xmid[k] = 0.5 *(X[i][k]+X[j][k]);

    Vec3D u = varFcn->getVelocity(Vmid);
    double a = varFcn->computeSoundSpeed(Vmid);

    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(Vmid);
    locbeta = tprec.getBeta(locMach,true);

    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    idti[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
    idti[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;

    /*
    if(fet) vis = fet->computeViscousTimeStep(Xmid, Vmid)*S*S;
    idtv[i] += vis;
    idtv[j] += vis;
    */
  }

}

template<int dim>
void EdgeSet::computeTimeStep2(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                              SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, Vec<double> &idtv,
                              TimeLowMachPrec &tprec, Vec<double>& edgeArea)
{

  double Vmid[dim];
  double Xmid[3];
  double vis = 0.0;
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  double locbeta=0.0;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double S = sqrt(normal[l] * normal[l]);
    double invS = 1.0 / S;

    Vec3D n = invS * normal[l];
    double ndot = invS * normalVel[l];

    for (int k=0; k<dim; ++k)
      Vmid[k] = 0.5 * (V[i][k] + V[j][k]);
    for (int k=0; k<3; k++)
      Xmid[k] = 0.5 *(X[i][k]+X[j][k]);

    Vec3D u = varFcn->getVelocity(Vmid);
    double a = varFcn->computeSoundSpeed(Vmid);

    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(Vmid);
    locbeta = tprec.getBeta(locMach,true);

    double beta2 = locbeta * locbeta;
    double coeff1 = fabs((1.0+beta2)*un);
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    idti[i] += 0.5*(coeff1+coeff2) * edgeArea[l];
    idti[j] += 0.5*(coeff1+coeff2) * edgeArea[l];

    /*
    if(fet) vis = fet->computeViscousTimeStep(Xmid, Vmid)*S*S;
    idtv[i] += vis;
    idtv[j] += vis;
    */
  }

}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void EdgeSet::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                              SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV,
			      Vec<double> &dIdti, Vec<double> &dIdtv, double dMach,
                              TimeLowMachPrec &tprec)
{

  double Vmid[dim], dVmid[dim];
  double Xmid[3], dXmid[3];
  double vis = 0.0;
  double dvis = 0.0;
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<Vec3D>& dNormal = geoState.getdEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();
  Vec<double>& dNormalVel = geoState.getdEdgeNormalVel();
  double locbeta=0.0;
  double dLocbeta=0.0;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double S = sqrt(normal[l] * normal[l]);
    double dS = normal[l] * dNormal[l] / sqrt(normal[l] * normal[l]);
    double invS = 1.0 / S;
    double dInvS = -dS / (S*S);

    Vec3D n = invS * normal[l];
    Vec3D dn = dInvS * normal[l] + invS * dNormal[l];
    double ndot = invS * normalVel[l];
    double dndot = dInvS * normalVel[l] + invS * dNormalVel[l];

    for (int k=0; k<dim; ++k) {
      Vmid[k] = 0.5 * (V[i][k] + V[j][k]);
      dVmid[k] = 0.5 * (dV[i][k] + dV[j][k]);
    }
    for (int k=0; k<3; k++) {
      Xmid[k] = 0.5 *(X[i][k]+X[j][k]);
      dXmid[k] = 0.5 *(dX[i][k]+dX[j][k]);
    }

    Vec3D u = varFcn->getVelocity(Vmid);
    Vec3D du = varFcn->getVelocity(dVmid);
    double a = varFcn->computeSoundSpeed(Vmid);
    double da = varFcn->computeDerivativeOfSoundSpeed(Vmid, dVmid, dMach);

    double un = u * n - ndot;
    double dun = du * n + u * dn - dndot;
    double locMach = varFcn->computeMachNumber(Vmid);
    locbeta = tprec.getBeta(locMach,true);
    double dLocMach = varFcn->computeDerivativeOfMachNumber(Vmid, dVmid, dMach);
    dLocbeta = tprec.getdBeta(locMach,dLocMach,true);

    double beta2 = locbeta * locbeta;
    double dbeta2 = 2.0*locbeta * dLocbeta;
    double coeff1 = (1.0+beta2)*un;
    double dCoeff1 = dbeta2*un + (1.0+beta2)*dun;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);
    double dCoeff2 = (((1.0-beta2)*un) * (-dbeta2*un + (1.0-beta2)*dun) + (2.0*locbeta*a) * (2.0*dLocbeta*a + 2.0*locbeta*da))  / pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    if (min(0.5*(coeff1-coeff2), 0.0) != 0.0) {
      dIdti[i] += 0.5*(dCoeff1-dCoeff2) * S + 0.5*(coeff1-coeff2) * dS;
      dIdti[j] += 0.5*(-dCoeff1-dCoeff2) * S + 0.5*(-coeff1-coeff2) * dS;
    }

    if(fet) {
      vis = fet->computeViscousTimeStep(Xmid, Vmid)*S*S;
      dvis = fet->computeDerivativeOfViscousTimeStep(Xmid, dXmid, Vmid, dVmid, dMach)*S*S + fet->computeViscousTimeStep(Xmid, Vmid)*2.0*S*dS;
    }
    dIdtv[i] += dvis;
    dIdtv[j] += dvis;

  }

}

//------------------------------------------------------------------------------
template<int dim>
void EdgeSet::computeTimeStep(VarFcn *varFcn, GeoState &geoState,
                              SVec<double,dim> &V, Vec<double> &dt,
                              TimeLowMachPrec &tprec, Vec<int> &fluidId, int subnum,
			      Vec<double>* umax)
{
  double Vmid[dim];
  double S, invS, ndot;
  double a, un, mach;
  double locbeta, beta2, coeff1, coeff2;
  int i,j;
  Vec3D n, u;

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();


  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    i = ptr[l][0];
    j = ptr[l][1];

    S = sqrt(normal[l] * normal[l]);
    invS = 1.0/S;
    n = invS * normal[l];
    ndot = invS * normalVel[l];

    if(fluidId[i]==fluidId[j]){
      for (int k=0; k<dim; ++k)
        Vmid[k] = 0.5 * (V[i][k] + V[j][k]);

      u = varFcn->getVelocity(Vmid);
      a = varFcn->computeSoundSpeed(Vmid, fluidId[i]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(Vmid, fluidId[i]);

      locbeta = tprec.getBeta(mach,true);
      beta2 = locbeta*locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
      if (umax) {
        (*umax)[i] += min(un,-1.0e-12)*S;
        (*umax)[j] += min(-un,-1.0e-12)*S;
      }
    }else{
      u = varFcn->getVelocity(V[i]);
      a = varFcn->computeSoundSpeed(V[i], fluidId[i]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[i],fluidId[i]);

      locbeta = tprec.getBeta(mach,true);
      beta2 = locbeta * locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;

      if (umax)
	(*umax)[i] += min(un,-1.0e-12)*S;

      u = varFcn->getVelocity(V[j]);
      a = varFcn->computeSoundSpeed(V[j], fluidId[j]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[j],fluidId[j]);

      locbeta = tprec.getBeta(mach,true);
      beta2 = locbeta * locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
      
      if (umax)
        (*umax)[j] += min(-un,-1.0e-12)*S;
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(int* locToGlobNodeMap, Vec<double> &irey, FluxFcn** fluxFcn,
                                     RecFcn* recFcn, ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     SVec<double,dim>& fluxes, SVec<int,2>& tag, int failsafe, int rshift)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double edgeirey, length;

  int ierr = 0;
	int l;

  for (int iEdge=0; iEdge<numSampledEdges; ++iEdge) {
		l = sampleMesh ? edgesConnectedToSampleNode[iEdge]: iEdge;

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (egrad)
      egrad->compute(l, i, j, elems, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    edgeirey = 0.5*(irey[i]+irey[j]);

    if (!rshift && locToGlobNodeMap)
    // check for negative pressure or density //
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    fluxFcn[BC_INTERNAL]->compute(length, edgeirey, normal[l], normalVel[l], Vi, Vj, flux);

    for (int k=0; k<dim; ++k) {
      fluxes[i][k] += flux[k];
      fluxes[j][k] -= flux[k];
    }

  }

  return(ierr);

}

template<int dim>
int EdgeSet::computeThinLayerViscousFiniteVolumeTerm(int* locToGlobNodeMap,
                                     VarFcn* varFcn,
                                     FemEquationTerm *fet,
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V,
                                     SVec<double,dim>& fluxes)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  double length;

  int ierr = 0;
  int l;

  double Fuhat;
  Vec3D Fvhat;
  double Tcg,Tj,Ti;
  double mu,lambda,kappa;
  double flux[dim];
  flux[0] = 0.0;
  double fluxl[5][3];
  fluxl[0][0] = fluxl[0][1] = fluxl[0][2] = 0.0;

  FemEquationTermNS* ns = dynamic_cast<FemEquationTermNS*>(fet);
  FemEquationTermSA* sa = dynamic_cast<FemEquationTermSA*>(fet);

  NavierStokesTerm* nsterm = NULL; 
  if (ns)
    nsterm = dynamic_cast<NavierStokesTerm*>(ns);
  else if (sa)
    nsterm = dynamic_cast<NavierStokesTerm*>(sa);
  else  {

    fprintf(stderr,"Error - Cannot construct a NavierStokesTerm in "
                        "EdgeSet::computeThinLayerViscousFiniteVolumeTerm");
  }

  double ooreynolds_mu,mutilde,mut,lambdat,kappat;
  if (ns) ooreynolds_mu = ns->get_ooreynolds_mu(); 
  if (sa) ooreynolds_mu = sa->get_ooreynolds_mu(); 
  
  for (int l=0; l<numSampledEdges; ++l) {    

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    // Compute the interal terms
    double area = normal[l].norm();
  
    //Vec3D xhat = normal[l];
    //xhat /= area;
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    Vec3D xhat(dx[0]/length,dx[1]/length,dx[2]/length);

    Ti = varFcn->computeTemperature(V[i]);
    Tj = varFcn->computeTemperature(V[j]);
    Tcg = 0.5*(Ti+Tj);

    Vec3D Vi = varFcn->getVelocity(V[i]);
    Vec3D Vj = varFcn->getVelocity(V[j]);
    
    Vec3D Vcg = 0.5*(Vi+Vj);
    
    double Ui = Vi*xhat, Uj = Vj*xhat;
    Vec3D vi = Vi-Ui*xhat, vj = Vj-Uj*xhat;

    mu     = nsterm->getViscoFcn()->compute_mu(Tcg);
    lambda = nsterm->getViscoFcn()->compute_lambda(Tcg,mu);
    kappa  = nsterm->getThermalCondFcn()->compute(Tcg);
    if (sa) {

      int nodeNum[3] = {i,j,i};
      double* Vl[] = {V[i],V[j],V[i],V[j]};
      sa->computeTurbulentTransportCoefficients(Vl, nodeNum, X, mu,lambda,
                 kappa, mutilde, mut, lambdat, kappat); 
      mu += mut;
      lambda += lambdat;
      kappa += kappat;
    }

    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;
    kappa  *= ooreynolds_mu;

/*    Fuhat = (lambda+2.0*mu)*(Uj-Ui)/length;
    Fvhat = mu*(vj-vi)/length;
    
    for (int k = 0; k < 3; ++k)
      flux[k+1] = xhat[k]*Fuhat + Fvhat[k];
    flux[4] = (0.5*(Ui+Uj)*(Uj-Ui)*(lambda+2.0*mu)+mu*0.5*(vj.normsq()-vi.normsq())+
               kappa*(Tj-Ti))/length;

    for (int k = 0; k < dim; ++k) {
      flux[k] *= area;
      fluxes[i][k] += flux[k];
      fluxes[j][k] -= flux[k];
    }
*/

    double gradu[3][3];
    for (int k = 0; k < 3; ++k) {
      for (int m = 0; m < 3; ++m)
        gradu[m][k] = (Vj[m]-Vi[m])/length*xhat[k];
    }
    double divu = gradu[0][0]+gradu[1][1]+gradu[2][2];
    
    for (int k = 0; k < 3; ++k) {
      for (int m = 0; m < 3; ++m) {

        fluxl[m+1][k] = lambda*divu*(k==m?1.0:0.0)+mu*(gradu[k][m]+gradu[m][k]);
      }
      fluxl[4][k] = Vcg[0]*fluxl[1][k]+Vcg[1]*fluxl[2][k]+Vcg[2]*fluxl[3][k] + kappa*(Tj-Ti)/length*xhat[k];
    }

    for (int k = 0; k < dim; ++k) {

      double ft = fluxl[k][0]*normal[l][0]+fluxl[k][1]*normal[l][1]+fluxl[k][2]*normal[l][2];
      fluxes[i][k] -= ft;
      fluxes[j][k] += ft;
    }
  }

  return 0;
}

template<int dim,class Scalar,int neq>
int EdgeSet::computeJacobianThinLayerViscousFiniteVolumeTerm(int* locToGlobNodeMap,
                                     VarFcn* varFcn,
                                     FemEquationTerm *fet,
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V,
                                     Vec<double>& ctrlVol,
/*                                     SVec<double,3>& faceJacX,
                                     SVec<double,3>& faceJacY,
                                     SVec<double,3>& faceJacZ,
                                     bool* boundaryFlag,
*/
                                     GenMat<Scalar,neq>& A)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  double length;

  int ierr = 0;
  int l;

  double Fuhat;
  Vec3D Fvhat;
  double Tcg,Tj,Ti;
  double mu,lambda,kappa;
  double flux[dim];
  flux[0] = 0.0;
  double jac_ii[dim*dim],jac_ij[dim*dim];
  double tmp[dim*dim];
  double tmp2[dim*dim];
  double tmp3[dim*dim];
  memset(jac_ii,0,sizeof(double)*dim*dim);
  memset(jac_ij,0,sizeof(double)*dim*dim);
  double Tgi[5],Tgj[5];

  FemEquationTermNS* ns = dynamic_cast<FemEquationTermNS*>(fet);
  FemEquationTermSA* sa = dynamic_cast<FemEquationTermSA*>(fet);
  FemEquationTermSAmean* sa_mean = dynamic_cast<FemEquationTermSAmean*>(fet);

  NavierStokesTerm* nsterm = NULL; 
  if (ns)
    nsterm = dynamic_cast<NavierStokesTerm*>(ns);
  else if (sa)
    nsterm = dynamic_cast<NavierStokesTerm*>(sa);
  else if (sa_mean)
    nsterm = dynamic_cast<NavierStokesTerm*>(sa_mean);
  else  {

    fprintf(stderr,"Error - Cannot construct a NavierStokesTerm in "
                        "EdgeSet::computeJacobianThinLayerViscousFiniteVolumeTerm");
  }

  double ooreynolds_mu = nsterm->get_ooreynolds_mu(); 
  double mutilde,mut,lambdat,kappat;
  //if (ns) ns->get_ooreynolds_mu(); 
  //if (sa) sa->get_ooreynolds_mu(); 
  
  for (int l=0; l<numSampledEdges; ++l) {    

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    // Compute the interal terms
    double area = normal[l].norm();
  
    Vec3D xhat = normal[l];
    xhat /= area;
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    Ti = varFcn->computeTemperature(V[i]);
    Tj = varFcn->computeTemperature(V[j]);
    varFcn->computeTemperatureGradient(V[i],Tgi);
    varFcn->computeTemperatureGradient(V[j],Tgj);
    Tcg = 0.5*(Ti+Tj);

    Vec3D Vi = varFcn->getVelocity(V[i]);
    Vec3D Vj = varFcn->getVelocity(V[j]);
    
    Vec3D Vcg = 0.5*(Vi+Vj);
    
    double Ui = Vi*xhat, Uj = Vj*xhat;
    Vec3D vi = Vi-Ui*xhat, vj = Vj-Uj*xhat;
    
    mu     = nsterm->getViscoFcn()->compute_mu(Tcg);
    lambda = nsterm->getViscoFcn()->compute_lambda(Tcg,mu);
    kappa  = nsterm->getThermalCondFcn()->compute(Tcg);

    if (sa) {

      int nodeNum[3] = {i,j,i};
      double* Vl[] = {V[i],V[j],V[i],V[j]};
      sa->computeTurbulentTransportCoefficients(Vl, nodeNum, X, mu,lambda,
                 kappa, mutilde, mut, lambdat, kappat); 
      mu += mut;
      lambda += lambdat;
      kappa += kappat;
    }


    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;
    kappa  *= ooreynolds_mu;
/*
    Fuhat = (lambda+2.0*mu)*(Uj-Ui)/length;
    Fvhat = mu*(vj-vi)/length;   
 
    for (int k = 0; k < 3; ++k) {
      
      for (int ll = 0; ll < 3; ++ll) {
        jac_ii[(k+1)*neq+ll+1] = -xhat[k]*(lambda+2.0*mu)/length*xhat[ll] - mu*((k==ll?1.0:0.0)-xhat[ll]*xhat[k])/length;
        jac_ij[(k+1)*neq+ll+1] = xhat[k]*(lambda+2.0*mu)/length*xhat[ll] + mu*((k==ll?1.0:0.0)-xhat[ll]*xhat[k])/length;
      }
      //flux[k+1] = xhat[k]*Fuhat + Fvhat[k];
    }

    for (int ll = 0; ll < 3; ++ll) {

      jac_ii[4*neq+ll+1] = (-(lambda+2.0*mu)*Ui*xhat[ll]-mu*(vi*(Vec3D((ll==0?1.0:0.0)-xhat[0]*xhat[ll],(ll==1?1.0:0.0)-xhat[1]*xhat[ll],(ll==2?1.0:0.0)-xhat[2]*xhat[ll]))))/length;
      jac_ij[4*neq+ll+1] = ((lambda+2.0*mu)*Uj*xhat[ll]+mu*(vj*(Vec3D((ll==0?1.0:0.0)-xhat[0]*xhat[ll],(ll==1?1.0:0.0)-xhat[1]*xhat[ll],(ll==2?1.0:0.0)-xhat[2]*xhat[ll]))))/length;
    }

    jac_ii[4*neq] = kappa*(Tgi[0])/length;
    jac_ij[4*neq] = -kappa*(Tgj[0])/length;
    jac_ii[4*neq+4] = kappa*(Tgi[4])/length;
    jac_ij[4*neq+4] = -kappa*(Tgj[4])/length;
//    flux[4] = (0.5*(Ui+Uj)*(Uj-Ui)*(lambda+2.0*mu)+mu*0.5*(vj.normsq()-vi.normsq())-
//               kappa*(Tj-Ti))/length;

*/
    double gradu[3][3];
    for (int k = 0; k < 3; ++k) {
      for (int m = 0; m < 3; ++m)
        gradu[m][k] = (Vj[m]-Vi[m])/length*xhat[k];
    }
    double divu = gradu[0][0]+gradu[1][1]+gradu[2][2];
    
    double ndx = normal[l][0]*xhat[0]+normal[l][1]*xhat[1]+normal[l][2]*xhat[2];
    jac_ii[4*dim] = kappa/length*(-Tgi[0])*ndx;
    jac_ij[4*dim] = kappa/length*(Tgj[0])*ndx; 
    
    jac_ii[4*dim+4] = kappa/length*(-Tgi[4])*ndx;
    jac_ij[4*dim+4] = kappa/length*(Tgj[4])*ndx; 

    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        jac_ii[(p+1)*dim+(q+1)] = lambda*(-xhat[q]/length)*normal[l][p]+
                                  mu*(-normal[l][q]/length*xhat[p]-(p==q?1.0:0.0)*ndx);
        jac_ij[(p+1)*dim+(q+1)] = lambda*(xhat[q]/length)*normal[l][p]+
                                  mu*(normal[l][q]/length*xhat[p]-(p==q?1.0:0.0)*ndx);
      }
    }

    for (int p = 0; p < 3; ++p) {

      jac_ii[4*dim+(p+1)] = Vcg[0]*jac_ii[(1)*dim+(p+1)] + 
                           Vcg[1]*jac_ii[(2)*dim+(p+1)] + 
                           Vcg[2]*jac_ii[(3)*dim+(p+1)];
      jac_ij[4*dim+(p+1)] = Vcg[0]*jac_ij[(1)*dim+(p+1)] + 
                           Vcg[1]*jac_ij[(2)*dim+(p+1)] + 
                           Vcg[2]*jac_ij[(3)*dim+(p+1)];

      jac_ii[4*dim+(p+1)] += 0.5*(lambda*divu*normal[l][p]+
                            mu*(gradu[0][p]*normal[l][0]+gradu[1][p]*normal[l][1]+gradu[2][p]*normal[l][2]+
                                gradu[p][0]*normal[l][0]+gradu[p][1]*normal[l][1]+gradu[p][2]*normal[l][2]));
      jac_ij[4*dim+(p+1)] += -0.5*(lambda*divu*normal[l][p]+
                            mu*(gradu[0][p]*normal[l][0]+gradu[1][p]*normal[l][1]+gradu[2][p]*normal[l][2]+
                                gradu[p][0]*normal[l][0]+gradu[p][1]*normal[l][1]+gradu[p][2]*normal[l][2]));
    }

    
 
    Scalar* jaci = A.getElem_ii(i);
    Scalar* jacj = A.getElem_ii(j);
    Scalar* jacij = A.getElem_ij(l);
    Scalar* jacji = A.getElem_ji(l);
    double voli = 1.0 / ctrlVol[i];
    double volj = 1.0 / ctrlVol[j];
    varFcn->postMultiplyBydVdU(V[i], jac_ii,tmp);
    if (masterFlag[l]) {
      for (int k = 0; k < neq; ++k) {
        for (int k2 = 0; k2 < neq; ++k2) {
          jaci[k*neq+k2] -= tmp[k*dim+k2];
        }
      }
    }
    for (int k = 0; k < neq; ++k) {
      for (int k2 = 0; k2 < neq; ++k2) {
        jacji[k*neq+k2] += tmp[k*dim+k2]*volj;
      }
    }
    varFcn->postMultiplyBydVdU(V[j], jac_ij,tmp);
    if (masterFlag[l]) {
      for (int k = 0; k < neq; ++k) {
        for (int k2 = 0; k2 < neq; ++k2) {
          jacj[k*neq+k2] += tmp[k*dim+k2];
        }
      }
    }
    for (int k = 0; k < neq; ++k) {
      for (int k2 = 0; k2 < neq; ++k2) {
        jacij[k*neq+k2] -= tmp[k*dim+k2]*voli;
      }
    }
 /* 
    if (boundaryFlag[i]) {
      varFcn->postMultiplyBydVdU(V[i],faceJacX[i], tmp);
      varFcn->postMultiplyBydVdU(V[i],faceJacY[i], tmp2);
      varFcn->postMultiplyBydVdU(V[i],faceJacZ[i], tmp3);
      if (masterFlag[l]) {
        for (int k = 0; k < neq*neq; ++k) {
     
          jaci[k] += 0.5*tmp[k]*voli*normal[l][0];
          jaci[k] += 0.5*tmp2[k]*voli*normal[l][1];
          jaci[k] += 0.5*tmp3[k]*voli*normal[l][2];
        }
      } 
      varFcn->postMultiplyBydVdU(V[j],faceJacX[i], tmp);
      varFcn->postMultiplyBydVdU(V[j],faceJacY[i], tmp2);
      varFcn->postMultiplyBydVdU(V[j],faceJacZ[i], tmp3);
      for (int k = 0; k < neq*neq; ++k) {
     
        jacij[k] += 0.5*tmp[k]*voli*voli*normal[l][0];
        jacij[k] += 0.5*tmp2[k]*voli*voli*normal[l][1];
        jacij[k] += 0.5*tmp3[k]*voli*voli*normal[l][2]; 
      }
    }
    if (boundaryFlag[j]) {
      varFcn->postMultiplyBydVdU(V[j],faceJacX[j], tmp);
      varFcn->postMultiplyBydVdU(V[j],faceJacY[j], tmp2);
      varFcn->postMultiplyBydVdU(V[j],faceJacZ[j], tmp3);
      if (masterFlag[l]) {
        for (int k = 0; k < neq*neq; ++k) {
     
          jacj[k] -= 0.5*tmp[k]*voli*normal[l][0];
          jacj[k] -= 0.5*tmp2[k]*voli*normal[l][1];
          jacj[k] -= 0.5*tmp3[k]*voli*normal[l][2];
        }
      } 
      varFcn->postMultiplyBydVdU(V[i],faceJacX[j], tmp);
      varFcn->postMultiplyBydVdU(V[i],faceJacY[j], tmp2);
      varFcn->postMultiplyBydVdU(V[i],faceJacZ[j], tmp3);
      for (int k = 0; k < neq*neq; ++k) {
     
        jacji[k] -= 0.5*tmp[k]*volj*volj*normal[l][0];
        jacji[k] -= 0.5*tmp2[k]*volj*volj*normal[l][1];
        jacji[k] -= 0.5*tmp3[k]*volj*volj*normal[l][2]; 
      }
    }
 */
  }
 
  return 0;
}
//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTermRestrict(int* locToGlobNodeMap, Vec<double>
		&irey, FluxFcn** fluxFcn, RecFcn* recFcn, ElemSet& elems, GeoState&
		geoState, SVec<double,3>& X, SVec<double,dim>& V, NodalGrad<dim>& ngrad,
		EdgeGrad<dim>* egrad, SVec<double,dim>& fluxes, SVec<int,2>& tag, int
		failsafe, int rshift) 
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double edgeirey, length;

  int ierr = 0;

	int l;
  for (int iSampledEdge=0; iSampledEdge<edgesConnectedToSampleNode.size(); ++iSampledEdge) {	//TODO: only some edges
		l = edgesConnectedToSampleNode[iSampledEdge];

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (egrad)
      egrad->compute(l, i, j, elems, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    edgeirey = 0.5*(irey[i]+irey[j]);

    if (!rshift)
    // check for negative pressure or density //
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    fluxFcn[BC_INTERNAL]->compute(length, edgeirey, normal[l], normalVel[l], Vi, Vj, flux);

    for (int k=0; k<dim; ++k) {
      fluxes[i][k] += flux[k];
      fluxes[j][k] -= flux[k];
    }

  }

  return(ierr);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void EdgeSet::computeDerivativeOfFiniteVolumeTerm(Vec<double> &irey, Vec<double> &dIrey, FluxFcn** fluxFcn, RecFcn* recFcn,
				      ElemSet& elems, GeoState& geoState, SVec<double,3>& X, SVec<double,3>& dX,
				      SVec<double,dim>& V, SVec<double,dim>& dV, NodalGrad<dim>& ngrad,
				      EdgeGrad<dim>* egrad, double dMach, SVec<double,dim>& dFluxes)
{

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<Vec3D>& dNormal = geoState.getdEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();
  Vec<double>& dNormalVel = geoState.getdEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  SVec<double,dim>& ddVdx = ngrad.getXderivative();
  SVec<double,dim>& ddVdy = ngrad.getYderivative();
  SVec<double,dim>& ddVdz = ngrad.getZderivative();

  double ddVij[dim], dddVij[dim], ddVji[dim], dddVji[dim], Vi[2*dim], dVi[2*dim], Vj[2*dim], dVj[2*dim], flux[dim], dFlux[dim];

//  double ddVijp[dim], ddVjip[dim], ddVijm[dim], ddVjim[dim], Vip[2*dim], Vim[2*dim], Vjp[2*dim], Vjm[2*dim];

  double edgeirey, dedgeirey;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    if (egrad)
      egrad->computeDerivative(l, i, j, elems, X, dX, V, dV, dVdx, dVdy, dVdz, ddVdx, ddVdy, ddVdz, dddVij, dddVji);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      double ddx[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};
      for (int k=0; k<dim; ++k) {
         ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
         ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
         dddVij[k] = ddx[0]*dVdx[i][k] + dx[0]*ddVdx[i][k] + ddx[1]*dVdy[i][k] + dx[1]*ddVdy[i][k] + ddx[2]*dVdz[i][k] + dx[2]*ddVdz[i][k];
         dddVji[k] = ddx[0]*dVdx[j][k] + dx[0]*ddVdx[j][k] + ddx[1]*dVdy[j][k] + dx[1]*ddVdy[j][k] + ddx[2]*dVdz[j][k] + dx[2]*ddVdz[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    recFcn->computeDerivative(V[i], dV[i], ddVij, dddVij, V[j], dV[j], ddVji, dddVji, dVi, dVj);

    edgeirey = 0.5*(irey[i]+irey[j]);

    dedgeirey = 0.5*(dIrey[i]+dIrey[j]);

    int k;
    for (k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
      dVi[k+dim] = dV[i][k];
      dVj[k+dim] = dV[j][k];
    }

    fluxFcn[BC_INTERNAL]->computeDerivative(edgeirey, dedgeirey, normal[l], dNormal[l], normalVel[l], dNormalVel[l], Vi, dVi, Vj, dVj, dMach, flux, dFlux);

    for (int k=0; k<dim; ++k) {
      dFluxes[i][k] += dFlux[k];
      dFluxes[j][k] -= dFlux[k];
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, Vec<int> &fluidId,
                                     FluidSelector &fluidSelector, 
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
				     SVec<double,dimLS>& phi, // needed for higher order computations
                                     NodalGrad<dimLS>& ngradLS,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,dimLS>& dPdx = ngradLS.getX();
  SVec<double,dimLS>& dPdy = ngradLS.getY();
  SVec<double,dimLS>& dPdz = ngradLS.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wi[2*dim], Wj[2*dim];
  double fluxi[dim], fluxj[dim];
  double gradphi[3];
  double gphii[3];
  double gphij[3];
  double Udummy[dim];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  int ierr=0;
  riemann.reset(it);

  programmedBurn = fluidSelector.getProgrammedBurn();


  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    for (int i = 0; i < dim; ++i) {
      fluxi[i] = fluxj[i] = 0.0;
    }

    bool isAtInterface = false;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    if (!higherOrderMF) {
      if (fluidId[i] == fluidId[j])
	recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
      else {
	for (int k = 0; k < dim; ++k) {
	  Vi[k] = V[i][k];
	  Vj[k] = V[j][k];
	}     
      }
    } else {

      if (fluidId[i] == fluidId[j] && !higherOrderMF->isCellCut(i) &&
	  !higherOrderMF->isCellCut(j)) {
	recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
      } else {
	isAtInterface = true;
	if (higherOrderMF->isCellCut(i) &&higherOrderMF->isCellCut(j))
	  continue;
      }
    }

    //if(Phi[i]*Phi[j] > 0.0) recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    //else                    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj, Phi[i], Phi[j]);

    // check for negative pressure or density //

    if (programmedBurn) {

      //std::cout << "Hello" << std::endl;
      varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
      varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);
    }

    if (!rshift && !isAtInterface)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], fluidId[i], fluidId[j]);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (fluidId[i]==fluidId[j] && !isAtInterface) { 	// same fluid
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, fluidId[i]);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
      riemann.resetInterfacialW(l);
    }
    else{			// interface
      //ngradLS returns nodal gradients of primitive phi
      // need fluidSelector to determine which level set to look at knowing which two fluids are considered at this interface
      int lsdim, burnTag;

      if (!(programmedBurn && programmedBurn->isDetonationInterface(fluidId[i],fluidId[j],burnTag)) ) {
	if (!higherOrderMF ||
	    (!higherOrderMF->isCellCut(i) && !higherOrderMF->isCellCut(j)))
	  lsdim = fluidSelector.getLevelSetDim(fluidId[i],fluidId[j],locToGlobNodeMap[i]+1,locToGlobNodeMap[j]+1);
	else if (higherOrderMF->isCellCut(j))
	  lsdim = fluidSelector.getLevelSetDim(fluidId[i],
					       higherOrderMF->getOtherFluidId(j,fluidId[i]),locToGlobNodeMap[i]+1,locToGlobNodeMap[j]+1);	  
	else
	  lsdim = fluidSelector.getLevelSetDim(higherOrderMF->getOtherFluidId(i,fluidId[j]),
					       fluidId[j],
					       locToGlobNodeMap[i]+1,locToGlobNodeMap[j]+1);	  

	gphii[0] = -dPdx[i][lsdim];
	gphii[1] = -dPdy[i][lsdim];
	gphii[2] = -dPdz[i][lsdim];
	gphij[0] = -dPdx[j][lsdim];
	gphij[1] = -dPdy[j][lsdim];
	gphij[2] = -dPdz[j][lsdim];
	for (int k=0; k<3; k++)
	  gradphi[k] = 0.5*(gphii[k]+gphij[k]);
	double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
	for (int k=0; k<3; k++)
	  gradphi[k] /= normgradphi;
      } else {
	
	double xmid[3];
	for (int k=0; k<3; k++)
	  xmid[k] = (X[j][k]+X[i][k])*0.5;

	programmedBurn->getDetonationNormal(burnTag,fluidId[i],fluidId[j], xmid, gradphi);
      }
      if (higherOrderMF) {
	
	assert(fluidId[i] != fluidId[j] || !higherOrderMF->isCellCut(i) ||
	       !higherOrderMF->isCellCut(j));
	
	// There are two cases.  In the first case the surrogate interface is the
	// same for both fluids.  This implies that the edge in question is cut by
	// the material interface.
	double iloc[3];
	// Step 1: Compute the intersection location (where the 0-contour of the level
	// set crosses this edge.
	double s = phi[j][lsdim]/(phi[j][lsdim]-phi[i][lsdim]);
	for (int k=0; k<3; k++)
	  iloc[k] = X[i][k]*s+X[j][k]*(1.0-s);
	if (!higherOrderMF->isCellCut(i) && 
	    !higherOrderMF->isCellCut(j)) {
	  
	  // Step 2: Extrapolate the values from cell i and cell j to the interface
	  if (s < 0.9) {
	    for (int k = 0; k < dim; ++k) {
	      Vi[k] = V[i][k]+dVdx[i][k]*(iloc[0]-X[i][0])+
		dVdy[i][k]*(iloc[1]-X[i][1])+
		dVdz[i][k]*(iloc[2]-X[i][2]);
	    }
	  } else {
	    for (int k = 0; k < dim; ++k) {
	      Vi[k] = V[i][k];
	    }
	  }

	  if (s > 0.1) {
	    for (int k = 0; k < dim; ++k) {
	      Vj[k] = V[j][k]+dVdx[j][k]*(iloc[0]-X[j][0])+
		dVdy[j][k]*(iloc[1]-X[j][1])+
		dVdz[j][k]*(iloc[2]-X[j][2]);
	    }
	  } else {
	    for (int k = 0; k < dim; ++k) {
	      Vj[k] = V[j][k];
	    }
	  }
	
          // Check for negative pressures/densities.
          // If a negative value is detected, drop back to first order extrapolation 
          // (i.e., the Riemann solution)
          if (Vi[0] <= 0.0)
            Vi[0] = V[i][0];
          if (Vi[4] <= 0.0)
            Vi[4] = V[i][4];
          if (Vj[0] <= 0.0)
            Vj[0] = V[j][0];
          if (Vj[4] <= 0.0)
            Vj[4] = V[j][4];
        
    
	  riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
					 Wi,Wj,i,j,l,dx,true);

	  // Step 3: Interpolate/Extrapolate back to the surrogate interface.
	  
	  for (int k = 0; k < dim; ++k) {
	    if (s < 0.9)
	      Vi[k] = (V[i][k]*(0.5-s)+Wi[k]*(0.5))/(1.0-s);

	    if (s > 0.1)
	      Vj[k] = (V[j][k]*(-0.5+s)+Wj[k]*(0.5))/s;
	  }

          // Check for negative pressures/densities.
          // If a negative value is detected, drop back to first order extrapolation 
          // (i.e., the Riemann solution)
          if (Vi[0] <= 0.0)
            Vi[0] = Wi[0];
          if (Vi[4] <= 0.0)
            Vi[4] = Wi[4];
          if (Vj[0] <= 0.0)
            Vj[0] = Wj[0];
          if (Vj[4] <= 0.0)
            Vj[4] = Wj[4];
        
	    
	  if (s < 0.9)
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  Vi, Vi, fluxi, fluidId[i]);
	  else
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  V[i], Wi, fluxi, fluidId[i]);
	  
	  if (s > 0.1)
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  Vj, Vj, fluxj, fluidId[j]);
	  else
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  Wj, V[j], fluxj, fluidId[j]);
	  

	  // Now extrapolate back to compute the riemann update for cells i/j
	  if (it == 1) {
	    SVec<double,dim> &rupdate = riemann.getRiemannUpdate();
	    Vec<double> &weight = riemann.getRiemannWeight();
            double updatei[dim],updatej[dim];
	    for (int k = 0; k < dim; ++k) {
              updatei[k] = (2.0*Vj[k]-V[j][k]);
              updatej[k] = (2.0*Vi[k]-V[i][k]);
            }
            if (updatei[0] <= 0.0 || updatei[4] <= 0.0) {
	      for (int k = 0; k < dim; ++k)
                updatei[k] = Vj[k];
            }
            if (updatej[0] <= 0.0 || updatej[4] <= 0.0) {
	      for (int k = 0; k < dim; ++k) 
                updatej[k] = Vi[k];
            }

	    for (int k = 0; k < dim; ++k) {
	      rupdate[i][k] += updatei[k];
	      rupdate[j][k] += updatej[k];
	    }
	    weight[i] += 1.0;
	    weight[j] += 1.0;	  
	  }

	} else if (higherOrderMF->isCellCut(j)) {
	  //std::cout << "Cut cell " << locToGlobNodeMap[j]+1 << std::endl;

	  // In this case, the cell j is cut be the interface, so its value does not exist
	  if (fabs(s-0.5) < 1.5) {
	    int fidj = higherOrderMF->getOtherFluidId(j, fluidId[i]);
	    higherOrderMF->template computeCutCellExtrapolations<dim>(j,fluidId[i],fidj, iloc, Vi,Vj,X);
	    /*if (locToGlobNodeMap[j]+1 == 6583) {
	      std::cout << "--------------------------" << std::endl;
	      std::cout << "Phi: " << phi[i][lsdim] << " " << phi[j][lsdim] << " " << s << std::endl;
	      std::cout << locToGlobNodeMap[i]+1 << "-" << locToGlobNodeMap[j]+1 << " (" <<  fluidId[i] << " , " << fidj << ") " << std::endl;
	      std::cout << "[" << iloc[0] << " " << iloc[1] << " " << iloc[2] << "] ";
	      std::cout << "[" << X[i][0] << " " << X[i][1] << " " << X[i][2] << "] " << std::endl;
	      std::cout << "[" << X[j][0] << " " << X[j][1] << " " << X[j][2] << "] " << std::endl;
	      for (int k = 0; k < dim; ++k) {
		std::cout << V[i][k] << " " << Vi[k] << "\t" << Vj[k] << std::endl;
	      }
	      higherOrderMF->template printCutCellData<dim>(j);
	      }*/
          
	    riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fidj,gradphi,varFcn,
					   Wi,Wj,i,j,l,dx,true);
	    
	    for (int k = 0; k < dim; ++k) {
	      Vi[k] = (V[i][k]*(0.5-s)+Wi[k]*(0.5))/(1.0-s);
	    }
	    
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  Vi, Vi, fluxi, fluidId[i]);
	    
	    // Now extrapolate back to compute the riemann update for cells i/j
	    if (it == 1 /*&& fluidId[j] != fluidId[i]*/) {
	      SVec<double,dim> &rupdate = riemann.getRiemannUpdate();
	      Vec<double> &weight = riemann.getRiemannWeight();
              double updatej[dim];
   	      for (int k = 0; k < dim; ++k) {
                updatej[k] = (2.0*Vi[k]-V[i][k]);
              }
              if (updatej[0] <= 0.0 || updatej[4] <= 0.0) {
	        for (int k = 0; k < dim; ++k) 
                  updatej[k] = Vi[k];
              }
	      for (int k = 0; k < dim; ++k) {
		rupdate[j][k] += updatej[k];
	      }
	      weight[j] += 1.0;
	    }
	  } else {

	    double Vxt[dim][3],Vjl[dim];
	    higherOrderMF->template getCutCellData<dim>(j, fluidId[i], Vjl, Vxt);
	    for (int k=0; k<dim; ++k) {
	      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
	      ddVji[k] = dx[0]*Vxt[k][0] + dx[1]*Vxt[k][1] + dx[2]*Vxt[k][2];
	    }
	    recFcn->compute(V[i], ddVij, Vjl, ddVji, Vi, Vj);
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, fluxi, fluidId[i]);
	    /*if (it == 1 && fluidId[j] != fluidId[i]) {
	      SVec<double,dim> &rupdate = riemann.getRiemannUpdate();
	      Vec<double> &weight = riemann.getRiemannWeight();
	      for (int k = 0; k < dim; ++k) {
		rupdate[j][k] += (2.0*Vi[k]-V[i][k]);
	      }
	      weight[j] += 1.0;
	    }*/
	  }

	} else if (higherOrderMF->isCellCut(i)) {
	  //std::cout << "Cut cell " << locToGlobNodeMap[i]+1 << std::endl;

	  // In this case, the cell j is cut be the interface, so its value does not exist
	  if (fabs(s-0.5) < 1.5) {
	    int fidi = higherOrderMF->getOtherFluidId(i, fluidId[j]);
	    higherOrderMF->template computeCutCellExtrapolations<dim>(i, fidi,fluidId[j],iloc, Vi,Vj,X); 
	    /*if (locToGlobNodeMap[i]+1 == 6583) {
	      std::cout << "--------------------------" << std::endl;
	      std::cout << "Phi: " << phi[i][lsdim] << " " << phi[j][lsdim] << " " << s << std::endl;
	      std::cout << locToGlobNodeMap[j]+1 << "-" << locToGlobNodeMap[i]+1 << " (" <<  fluidId[j] << " , " << fidi << ") " << std::endl;
	      std::cout << "[" << iloc[0] << " " << iloc[1] << " " << iloc[2] << "] ";
	      std::cout << "[" << X[i][0] << " " << X[i][1] << " " << X[i][2] << "] " << std::endl;
	      std::cout << "[" << X[j][0] << " " << X[j][1] << " " << X[j][2] << "] " << std::endl;	    
	      for (int k = 0; k < dim; ++k) {
		std::cout << V[j][k] << " " << Vj[k] << "\t" << Vi[k] << std::endl;
	      }
	      higherOrderMF->template printCutCellData<dim>(i);
	      }*/
	    
	    riemann.computeRiemannSolution(Vi,Vj,fidi,fluidId[j],gradphi,varFcn,
					   Wi,Wj,i,j,l,dx,true);
	    
	    for (int k = 0; k < dim; ++k) {
	      Vj[k] = (V[j][k]*(-0.5+s)+Wj[k]*(0.5))/s;
	    }
	    
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
					  Vj, Vj, fluxj, fluidId[j]);
	    
	    // Now extrapolate back to compute the riemann update for cells i/j
	    if (it == 1) {
	      SVec<double,dim> &rupdate = riemann.getRiemannUpdate();
	      Vec<double> &weight = riemann.getRiemannWeight();
              double updatei[dim];
   	      for (int k = 0; k < dim; ++k) {
                updatei[k] = (2.0*Vj[k]-V[j][k]);
              }
              if (updatei[0] <= 0.0 || updatei[4] <= 0.0) {
	        for (int k = 0; k < dim; ++k) 
                  updatei[k] = Vj[k];
              }
	      for (int k = 0; k < dim; ++k) {
		rupdate[i][k] += updatei[k];
	      }
	      weight[i] += 1.0;
	    }
	  } else { 

	    double Vxt[dim][3],Vjl[dim];
	    higherOrderMF->template getCutCellData<dim>(i, fluidId[j], Vjl, Vxt);
	    for (int k=0; k<dim; ++k) {
	      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
	      ddVij[k] = dx[0]*Vxt[k][0] + dx[1]*Vxt[k][1] + dx[2]*Vxt[k][2];
	    }
	    recFcn->compute(Vjl, ddVij, V[j], ddVji, Vi, Vj);
	    fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, fluxj, fluidId[i]);	  
	    /*if (it == 1) {
	      SVec<double,dim> &rupdate = riemann.getRiemannUpdate();
	      Vec<double> &weight = riemann.getRiemannWeight();
	      for (int k = 0; k < dim; ++k) {
		rupdate[i][k] += (2.0*Vj[k]-V[j][k]);
	      }
	      weight[i] += 1.0;
	    }*/
	  }
	}
      }	else {

	riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
				       Wi,Wj,i,j,l,dx,false);
	fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
				      Vi, Wi, fluxi, fluidId[i]);
	fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
				      Wj, Vj, fluxj, fluidId[j]);
      }

      for (int k=0; k<dim; k++){
        fluxes[i][k] += fluxi[k];
        fluxes[j][k] -= fluxj[k];
      }
    }
  }
  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji,
                                     LevelSetStructure& LSS, bool linRecAtInterface, Vec<int> &fluidId,
                                     int Nriemann, SVec<double,3>* Nsbar, FluidSelector &fluidSelector,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     NodalGrad<dimLS>& ngradLS,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{
  // ------------------------------------------------
  //  Preparation -- General Info. 
  // ------------------------------------------------
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  int ierr=0;
  int burnTag;
  riemann.reset(it);

  // ------------------------------------------------
  //  Preparation -- Fluid/Structure Part
  // ------------------------------------------------
  int farfieldFluid = 0; //assume that intersector assigns Id=0 for outside fluid.
  double Wstar[2*dim]; //FS Riemann solution
  Vec3D normalDir; //Normal direction used for setting the Riemann problem

  // ------------------------------------------------
  //  Preparation -- Fluid/Fluid Part
  // ------------------------------------------------
  SVec<double,dimLS>& dPdx = ngradLS.getX();
  SVec<double,dimLS>& dPdy = ngradLS.getY();
  SVec<double,dimLS>& dPdz = ngradLS.getZ();
  double Wi[2*dim], Wj[2*dim]; //FF Riemann solution
  double gradphi[3], gphii[3], gphij[3], Udummy[dim];

  // ------------------------------------------------
  //  THE MAIN EDGE LOOP...
  // ------------------------------------------------
  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,l);
    bool iActive, jActive;
    if(LSS.withCracking()) { //assumption: only occluded nodes are inactive in the case of cracking simulations
      iActive = !LSS.isOccluded(0.0,i);
      jActive = !LSS.isOccluded(0.0,j);
    } else {
      iActive = LSS.isActive(0.0,i);
      jActive = LSS.isActive(0.0,j);
    } 

    if(!iActive && !jActive) continue; //this edge is inside a solid body!

    // ------------------------------------------------
    //  Reconstruction without crossing the FS interface.
    // ------------------------------------------------
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    if (!intersect)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.
    else { // linRec at interface using Wstar
      if (!linRecAtInterface) // just set Vi = V[i], Vj = V[j]
        for(int k=0; k<dim; k++) {
          Vi[k] = V[i][k];
          Vj[k] = V[j][k];
        }
      else {//linRec at interface using Wstar
        double Vtemp[2*dim];
        if (Wstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
          for (int k=0; k<dim; k++) Vi[k] = V[i][k];
        } else recFcn->compute(V[i], ddVij, Wstarij[l], ddVji, Vi, Vtemp);
        if (Wstarji[l][0]<1e-8) {
          for (int k=0; k<dim; k++) Vj[k] = V[j][k];
        } else recFcn->compute(Wstarji[l], ddVij, V[j], ddVji, Vtemp, Vj);
      }
    }    
    
    programmedBurn = fluidSelector.getProgrammedBurn();

    if (programmedBurn) {

      //std::cout << "Hello" << std::endl;
      varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
      varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);
    }

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], fluidId[i], fluidId[j]);
    if (ierr) continue;

    if(it>0)
      for(int k=0;k<dim;k++)
        Wstarij[l][k] = Wstarji[l][k] = 0.0; //clean-up for re-fill.

    for (int k=0; k<dim; ++k) { // 1st-order values
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    // --------------------------------------------------------
    //  Compute the flux along this edge.
    //    Step 1. If e(i,j) intersects the structure -> FS flux
    //    Step 2. If Idi!=Idj -> FF flux
    //    Step 3. Otherwise, the usual single-phase flux. 
    // --------------------------------------------------------
    if(intersect) {

      // for node i
      if(iActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
        switch (Nriemann) { // normal should point to this node (i). 
          case 0: //structure normal
            normalDir = (normal[l]*resij.gradPhi>=0.0) ? -1.0*resij.gradPhi : resij.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(LSS.fluidModel(0.0,i)==farfieldFluid)  normalDir =      Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            else                                     normalDir = -1.0*Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }

        riemann.computeFSIRiemannSolution(Vi,resij.normVel,normalDir,varFcn,Wstar,j,fluidId[i]);

        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi, fluidId[i], false);
          for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
        switch (Nriemann) {
          case 0: //structure normal
            normalDir = (normal[l]*resji.gradPhi>=0.0) ? resji.gradPhi : -1.0*resji.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(LSS.fluidModel(0.0,j)==farfieldFluid)  normalDir =      Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            else                                     normalDir = -1.0*Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }

        riemann.computeFSIRiemannSolution(Vj,resji.normVel,normalDir,varFcn,Wstar,i,fluidId[j]);

        if (it>0)
          for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj, fluidId[j], false);
          for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
        }

      }
    }
    else if(fluidId[i]!=fluidId[j]) { //NOTE: It's NOT equivalent with checking Phi_i x Phi_j < 0!
      if(!masterFlag[l]) continue;
    
      // Force constant reconstruction at the interface.
      for (int k = 0; k < dim; ++k) {
        Vi[k] = V[i][k];
        Vj[k] = V[j][k];
      }      

      //ngradLS returns nodal gradients of primitive phi
      // need fluidSelector to determine which level set to look at knowing which two fluids are considered at this interface   
      if (!(programmedBurn && programmedBurn->isDetonationInterface(fluidId[i],fluidId[j],burnTag)) ) {

	int lsdim = fluidSelector.getLevelSetDim(fluidId[i],fluidId[j],locToGlobNodeMap[i]+1,locToGlobNodeMap[j]+1);
	gphii[0] = -dPdx[i][lsdim];
	gphii[1] = -dPdy[i][lsdim];
	gphii[2] = -dPdz[i][lsdim];
	gphij[0] = -dPdx[j][lsdim];
	gphij[1] = -dPdy[j][lsdim];
	gphij[2] = -dPdz[j][lsdim];
	for (int k=0; k<3; k++)
	  gradphi[k] = 0.5*(gphii[k]+gphij[k]);
	double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
	for (int k=0; k<3; k++)
	  gradphi[k] /= normgradphi;
      } else {
	double xmid[3];
	for (int k=0; k<3; k++)
	  xmid[k] = (X[j][k]+X[i][k])*0.5;

	programmedBurn->getDetonationNormal(burnTag,fluidId[i],fluidId[j], xmid, gradphi);
      }

      //if (myrank == 46 && (fluidId[i] == 3 ||fluidId[j] == 3) )
      //	std::cout << "difffluid" << std::endl;

      riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
				     Wi,Wj,i,j,l,dx,false);

      checkReconstructedValues(i, j, Wi, Wj, varFcn, locToGlobNodeMap,
      			       failsafe, tag, Vi, Vj, fluidId[i], fluidId[j]);

      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
				    Vi, Wi, fluxi, fluidId[i]);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
				    Wj, Vj, fluxj, fluidId[j]);

      for (int k=0; k<dim; k++){
        fluxes[i][k] += fluxi[k];
        fluxes[j][k] -= fluxj[k];
      }
    }

    else {
      if(!masterFlag[l]) continue;
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, fluidId[i]);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
      riemann.resetInterfacialW(l);
    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
inline
void extendedLinearExtrapolationToIntersection(ElemSet& elems, int idxTet, int idxFace, 
		double face_r, double face_t, SVec<double,3>& X, SVec<double,dim>& V, double* Wstar, 
		double alpha, double length, int i) {
  int n0_loc = elems[idxTet].faceDef(idxFace, 0);
  int n1_loc = elems[idxTet].faceDef(idxFace, 1);
  int n2_loc = elems[idxTet].faceDef(idxFace, 2);
  int n0 = elems[idxTet].nodeNum(n0_loc);
  int n1 = elems[idxTet].nodeNum(n1_loc);
  int n2 = elems[idxTet].nodeNum(n2_loc);

  double Vface[dim];
  double Xface[3];

  for (int k=0; k<dim; ++k)
    Vface[k] = V[n2][k]+face_r*(V[n0][k]-V[n2][k])+face_t*(V[n1][k]-V[n2][k]);
  for (int k=0; k<3; ++k)
    Xface[k] = X[n2][k]+face_r*(X[n0][k]-X[n2][k])+face_t*(X[n1][k]-X[n2][k]);
  double alpha_f = sqrt((Xface[0]-X[i][0])*(Xface[0]-X[i][0])+
				   		(Xface[1]-X[i][1])*(Xface[1]-X[i][1])+
				   		(Xface[2]-X[i][2])*(Xface[2]-X[i][2]))/length;
  for (int k=0; k<dim; ++k)
    Wstar[k] = Wstar[k]+((0.5-alpha)/(1.0+alpha_f-alpha))*(Vface[k]-Wstar[k]);
  return;
}

//------------------------------------------------------------------------------

inline
bool notAllActive(Elem& elem, int idxFace, LevelSetStructure& LSS) {
  int n0_loc = elem.faceDef(idxFace, 0);
  int n1_loc = elem.faceDef(idxFace, 1);
  int n2_loc = elem.faceDef(idxFace, 2);
  int n0 = elem.nodeNum(n0_loc);
  int n1 = elem.nodeNum(n1_loc);
  int n2 = elem.nodeNum(n2_loc);
  return ((!LSS.isActive(0.0,n0))||(!LSS.isActive(0.0,n1))||(!LSS.isActive(0.0,n2)));
}

//------------------------------------------------------------------------------

//inline
//int computePrimitiveJacEuler3D(char dir, double* V, double** J) {
//  double gamma = 1.4;
//  double rho = V[0];
//  if (rho<=0.0)
//    return 1;
//  double u   = V[1];
//  double v   = V[2];
//  double w   = V[3];
//  double p   = V[4];
//  double nx  = 0.0;
//  double ny  = 0.0;
//  double nz  = 0.0;
//  switch (dir) {
//  	case 'x':
//  		nx = 1.0;	break;
//  	case 'y':
//  		ny = 1.0;	break;
//  	case 'z':
//  		nz = 1.0;	break;
//  }
//  J[0][0] = nx*u+ny*v+nz*w;
//  J[0][1] = rho*nx;	J[0][2] = rho*ny;	J[0][3] = rho*nz;	
//  J[0][4] = 0.0;
//  
//  J[1][0] = 0.0;
//  J[1][1] = nx*u+ny*v+nz*w;
//  J[1][2] = 0.0;	J[1][3] = 0.0;	J[1][4] = nx/rho;
//  
//  J[2][0] = 0.0;	J[2][1] = 0.0;
//  J[2][2] = nx*u+ny*v+nz*w;
//  J[2][3] = 0.0;	J[2][4] = ny/rho;
//  
//  J[3][0] = 0.0;	J[3][1] = 0.0;	J[3][2] = 0.0;
//  J[3][3] = nx*u+ny*v+nz*w;
//  J[3][4] = nz/rho;
//  
//  J[4][0] = 0.0;
//  J[4][1] = gamma*p*nx;	J[4][2] = gamma*p*ny;	J[4][3] = gamma*p*nz;
//  J[4][4] = nx*u+ny*v+nz*w;
//  return 0;
//}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim>& Wstarij,
                                     SVec<double,dim>& Wstarji, LevelSetStructure &LSS, bool linRecAtInterface,
                                     Vec<int> &fluidId, int Nriemann, SVec<double,3> *Nsbar,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{
  int farfieldFluid = 0; 
//  int BuggyNode = 340680; //KEVIN: DEBUG ONLY

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim],Udummy[dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  Vec3D normalDir;

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,l);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);

//    if(locToGlobNodeMap[i]+1==BuggyNode || locToGlobNodeMap[i]+1==BuggyNode) {
//      fprintf(stderr,"--- (%d(%d), %d(%d)): intersect = %d, normal = %e %e %e\n", locToGlobNodeMap[i]+1, iActive, locToGlobNodeMap[j]+1, jActive, intersect, normal[l][0], normal[l][1], normal[l][2]);
//    }
    
    if( !iActive && !jActive ) {
      if(it>0) for(int k=0;k<dim;k++) Wstarij[l][k] = Wstarji[l][k] = 0.0; //clean-up Wstar
      continue;
    }

    // ------------------------------------------------
    //  Reconstruction without crossing the interface.
    // ------------------------------------------------
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    if (iActive && jActive && !intersect)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.

    else { // linRec at interface using Wstar
      if (!linRecAtInterface) // just set Vi = V[i], Vj = V[j]
        for(int k=0; k<dim; k++) {
          Vi[k] = V[i][k];
          Vj[k] = V[j][k];
        }
      else { // linRec at interface using Wstar        
        // Case 1: i and j are both active (but i-j intersects)
        if (iActive&&jActive) {
          double Vtemp[2*dim];
          if (Wstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) Vi[k] = V[i][k];
          } else recFcn->compute(V[i], ddVij, Wstarij[l], ddVji, Vi, Vtemp);
          if (Wstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) Vj[k] = V[j][k];
          } else recFcn->compute(Wstarji[l], ddVij, V[j], ddVji, Vtemp, Vj);
        }
        // Case 2: i is active, j is inactive
        else if (iActive)
          if (Wstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(V[i], ddVij, Wstarij[l], ddVji, Vi, Vj);
        // Case 3: i is inactive, j is active
        else if (jActive)
          if (Wstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(Wstarji[l], ddVij, V[j], ddVji, Vi, Vj);
      }
    } 
      
    varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
    varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);

    if(it>0)
      for(int k=0;k<dim;k++)
        Wstarij[l][k] = Wstarji[l][k] = 0.0;

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], fluidId[i], fluidId[j]); //also checking reconstructed values acrossinterface.

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }


//    if(locToGlobNodeMap[i]+1==BuggyNode || locToGlobNodeMap[i]+1==BuggyNode) {
//      fprintf(stderr,"--- %d: %e %e %e %e %e | %d: %e %e %e %e %e\n", locToGlobNodeMap[i]+1, Vi[0], Vi[1], Vi[2], Vi[3], Vi[4], 
//                      locToGlobNodeMap[j]+1, Vj[0], Vj[1], Vj[2], Vj[3], Vj[4]);
//    }
    // --------------------------------------------------------
    //                   Compute fluxes
    // --------------------------------------------------------
    if (!intersect) {  // same fluid
      if (!masterFlag[l]) continue; //not a master edge

      if(!(iActive && jActive)) {
        fprintf(stderr,"Really odd... (%d(%d),%d(%d): intersect = %d; iSwept = %d, jSwept = %d, iOccluded = %d, jOcculded = %d, model = %d/%d.\n",
                locToGlobNodeMap[i]+1, iActive, locToGlobNodeMap[j]+1, jActive, intersect, LSS.isSwept(0.0,i), LSS.isSwept(0.0,j), 
                LSS.isOccluded(0.0,i), LSS.isOccluded(0.0,j), LSS.fluidModel(0.0,i), LSS.fluidModel(0.0,j));
        continue;
      }

      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, fluidId[i]);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
//      if(locToGlobNodeMap[i]+1==BuggyNode || locToGlobNodeMap[i]+1==BuggyNode)
//        fprintf(stderr,"--- A flux(%d,%d) = %e %e %e %e %e\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, flux[0], flux[1], flux[2], flux[3], flux[4]);
      
    }
    else{// interface

      // for node i
      if(iActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);

        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      resij.gradPhi;
            else                                normalDir = -1.0*resij.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());

        riemann.computeFSIRiemannSolution(Vi,resij.normVel,normalDir,varFcn,Wstar,j,fluidId[i]);

        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi, fluidId[i], false);
          for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
//          if(locToGlobNodeMap[i]+1==BuggyNode || locToGlobNodeMap[i]+1==BuggyNode)
//            fprintf(stderr,"--- B fluxi(%d,%d) = %e %e %e %e %e\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, fluxi[0], fluxi[1], fluxi[2], fluxi[3], fluxi[4]);
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);

        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      resji.gradPhi;
            else                                normalDir = -1.0*resji.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());
  
        riemann.computeFSIRiemannSolution(Vj,resji.normVel,normalDir,varFcn,Wstar,i,fluidId[j]);
  
        if (it>0)
          for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj, fluidId[j], false);
          for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
//          if(locToGlobNodeMap[i]+1==BuggyNode || locToGlobNodeMap[i]+1==BuggyNode)
//            fprintf(stderr,"--- C fluxj(%d,%d) = %e %e %e %e %e\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, fluxj[0], fluxj[1], fluxj[2], fluxj[3], fluxj[4]);
        }
      }

    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn, ElemSet& elems, 
									 GeoState& geoState, SVec<double,3>& X, SVec<double,dim>& V, 
									 SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji, 
									 Vec<int>& countWstarij, Vec<int>& countWstarji, 
									 LevelSetStructure &LSS, bool linRecAtInterface, 
									 Vec<int> &fluidId, int Nriemann, SVec<double,3> *Nsbar, 
									 double dt, double alpha, NodalGrad<dim>& ngrad, 
									 EdgeGrad<dim>* egrad, SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift, 
									 V6NodeData (*v6data)[2])
{
  int farfieldFluid = 0; 

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim],Udummy[dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  Vec3D normalDir;

  int ierr=0;
  riemann.reset(it);

  //double clip_alpha_min = 1e-1;
  //double clip_alpha_max = 1e-1;

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue; //not a master edge
    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,l);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);

    if( !iActive && !jActive ) {
      if(it>0) for(int k=0;k<dim;k++) Wstarij[l][k] = Wstarji[l][k] = 0.0; //clean-up Wstar
      continue;
    }

    // ------------------------------------------------
    //  Reconstruction without crossing the interface.
    // ------------------------------------------------
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    if (iActive && jActive && !intersect)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.
    else { // linRec at interface using Wstar
      if (!linRecAtInterface) // just set Vi = V[i], Vj = V[j]
        for(int k=0; k<dim; k++) {
          Vi[k] = V[i][k];
          Vj[k] = V[j][k];
        }
      else { // linRec at interface using Wstar        
		if (iActive) {
		  LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0,l,true);
		  for (int k=0; k<dim; k++) Vi[k] = V[i][k]+(1.0-resij.alpha)*ddVij[k];
		}
		else
		  for (int k=0; k<dim; k++) Vi[k] = V[i][k];
		if (jActive) {
		  LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0,l,false);
		  for (int k=0; k<dim; k++) Vj[k] = V[j][k]-(1.0-resji.alpha)*ddVji[k];
		}
		else
		  for (int k=0; k<dim; k++) Vj[k] = V[j][k];
      }
    }
      
    varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
    varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], fluidId[i], fluidId[j]); //also checking reconstructed values acrossinterface.

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    // --------------------------------------------------------
    //                   Compute fluxes
    // --------------------------------------------------------
    if (!intersect) {  // same fluid
      if(!(iActive && jActive)) {
        fprintf(stderr,"Really odd... (%d(%d),%d(%d): intersect = %d; iSwept = %d, jSwept = %d, iOccluded = %d, jOcculded = %d, model = %d/%d.\n",
                locToGlobNodeMap[i]+1, iActive, locToGlobNodeMap[j]+1, jActive, intersect, LSS.isSwept(0.0,i), LSS.isSwept(0.0,j), 
                LSS.isOccluded(0.0,i), LSS.isOccluded(0.0,j), LSS.fluidModel(0.0,i), LSS.fluidModel(0.0,j));
        continue;
      }

      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, fluidId[i]);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{// interface
      if(iActive) {// for node i
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      resij.gradPhi;
            else                                normalDir = -1.0*resij.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());
		if (countWstarij[l] == 0)
          riemann.computeFSIRiemannSolution(Vi,resij.normVel,normalDir,varFcn,Wstar,j,fluidId[i]);
		else {
		  for (int k=0; k<dim; k++) Wstar[k] = Wstarij[l][k];
		  double dVdU[dim*dim], Jacx[dim*dim], Jacy[dim*dim], Jacz[dim*dim];
		  double axisnrm[3];
		  for (int k=0; k<dim*dim; k++) dVdU[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacx[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacy[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacz[k] = 0.0;
		  varFcn->getVarFcnBase(fluidId[i])->computedVdU(Wstar,dVdU);
		  axisnrm[0] = 1.0; axisnrm[1] = 0.0; axisnrm[2] = 0.0;
		  varFcn->getVarFcnBase(fluidId[i])->computedFdV(axisnrm,Wstar,Jacx);
		  axisnrm[0] = 0.0; axisnrm[1] = 1.0; axisnrm[2] = 0.0;
		  varFcn->getVarFcnBase(fluidId[i])->computedFdV(axisnrm,Wstar,Jacy);
		  axisnrm[0] = 0.0; axisnrm[1] = 0.0; axisnrm[2] = 1.0;
		  varFcn->getVarFcnBase(fluidId[i])->computedFdV(axisnrm,Wstar,Jacz);
		  for (int k=0; k<dim; k++)
		    for (int k1=0; k1<dim; k1++)
			  for (int k2=0; k2<dim; k2++)
		        Wstar[k] -= dt*(dVdU[k*dim+k1]*Jacx[k1*dim+k2]*dVdx[i][k2]+
		  	  	  		        dVdU[k*dim+k1]*Jacy[k1*dim+k2]*dVdy[i][k2]+
		  	  			        dVdU[k*dim+k1]*Jacz[k1*dim+k2]*dVdz[i][k2]);
		}
		if (it>0) {
		  for (int k=0; k<dim; k++) Wstarij[l][k] = Wstar[k];
		  countWstarij[l] += 1;
		}
		if (v6data==NULL)
		  for (int k=0; k<dim; k++) Wstar[k] = V[i][k]+(0.5/max(1.0-resij.alpha,alpha))*(Wstar[k]-V[i][k]);
		else {
		  int idxTet = v6data[l][0].tet;
		  int idxFace = v6data[l][0].face;
		  double face_r = v6data[l][0].r;
		  double face_t = v6data[l][0].t;
		  if ((idxTet<0)||(idxTet>=elems.size())||notAllActive(elems[idxTet],idxFace,LSS))
			for (int k=0; k<dim; k++) Wstar[k] = V[i][k]+(0.5/max(1.0-resij.alpha,alpha))*(Wstar[k]-V[i][k]);
		  else
			extendedLinearExtrapolationToIntersection<dim>(elems,idxTet,idxFace,face_r,face_t,
			  		X,V,Wstar,resij.alpha,length,i);
		}

        //fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Wstar, fluxi, fluidId[i], false);
		double fluxnrm[3] = {normal[l][0],normal[l][1],normal[l][2]};
		varFcn->getVarFcnBase(fluidId[i])->computeFofV(fluxnrm,Wstar,fluxi);
        for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
      }
      if(jActive){// for node j
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      resji.gradPhi;
            else                                normalDir = -1.0*resji.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());
		if (countWstarji[l] == 0)
          riemann.computeFSIRiemannSolution(Vj,resji.normVel,normalDir,varFcn,Wstar,i,fluidId[j]);
		else {
		  for (int k=0; k<dim; k++) Wstar[k] = Wstarji[l][k];
		  double dVdU[dim*dim], Jacx[dim*dim], Jacy[dim*dim], Jacz[dim*dim];
		  double axisnrm[3];
		  for (int k=0; k<dim*dim; k++) dVdU[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacx[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacy[k] = 0.0;
		  for (int k=0; k<dim*dim; k++) Jacz[k] = 0.0;
		  varFcn->getVarFcnBase(fluidId[j])->computedVdU(Wstar,dVdU);
		  axisnrm[0] = 1.0; axisnrm[1] = 0.0; axisnrm[2] = 0.0;
		  varFcn->getVarFcnBase(fluidId[j])->computedFdV(axisnrm,Wstar,Jacx);
		  axisnrm[0] = 0.0; axisnrm[1] = 1.0; axisnrm[2] = 0.0;
		  varFcn->getVarFcnBase(fluidId[j])->computedFdV(axisnrm,Wstar,Jacy);
		  axisnrm[0] = 0.0; axisnrm[1] = 0.0; axisnrm[2] = 1.0;
		  varFcn->getVarFcnBase(fluidId[j])->computedFdV(axisnrm,Wstar,Jacz);
		  for (int k=0; k<dim; k++)
		    for (int k1=0; k1<dim; k1++)
			  for (int k2=0; k2<dim; k2++)
		        Wstar[k] -= dt*(dVdU[k*dim+k1]*Jacx[k1*dim+k2]*dVdx[j][k2]+
		  	  	  		        dVdU[k*dim+k1]*Jacy[k1*dim+k2]*dVdy[j][k2]+
		  	  			        dVdU[k*dim+k1]*Jacz[k1*dim+k2]*dVdz[j][k2]);
		}
		if (it>0) {
		  for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
		  countWstarji[l] += 1;
		}
  
		if (v6data==NULL)
		  for (int k=0; k<dim; k++) Wstar[k] = V[j][k]+(0.5/max(1.0-resji.alpha,alpha))*(Wstar[k]-V[j][k]);
		else {
		  int idxTet = v6data[l][1].tet;
		  int idxFace = v6data[l][1].face;
		  double face_r = v6data[l][1].r;
		  double face_t = v6data[l][1].t;
		  if ((idxTet<0)||(idxTet>=elems.size())||notAllActive(elems[idxTet],idxFace,LSS))
			for (int k=0; k<dim; k++) Wstar[k] = V[j][k]+(0.5/max(1.0-resji.alpha,alpha))*(Wstar[k]-V[j][k]);
		  else
			extendedLinearExtrapolationToIntersection<dim>(elems,idxTet,idxFace,face_r,face_t,
			  		X,V,Wstar,resji.alpha,length,j);
		}
        //fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Wstar, fluxj, fluidId[j], false);
		double fluxnrm[3] = {normal[l][0],normal[l][1],normal[l][2]};
		varFcn->getVarFcnBase(fluidId[j])->computeFofV(fluxnrm,Wstar,fluxj);
        for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
      }
    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void EdgeSet::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                      ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, Vec<int>& fluidId, NodalGrad<dim>& ngrad,
                                      NodalGrad<dimLS> &ngradLS,
                                      EdgeGrad<dim>* egrad, SVec<double,dimLS>& Phi, SVec<double,dimLS>& PhiF,
                                      LevelSetStructure *LSS)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx  = ngrad.getX();
  SVec<double,dim>& dVdy  = ngrad.getY();
  SVec<double,dim>& dVdz  = ngrad.getZ();
  // in this routine Phi denotes the "conservative phi" ie (rho*phi)
  SVec<double,dimLS>& dPhidx = ngradLS.getX();
  SVec<double,dimLS>& dPhidy = ngradLS.getY();
  SVec<double,dimLS>& dPhidz = ngradLS.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim];
  double ddPij[dimLS], ddPji[dimLS], Pi[2*dimLS], Pj[2*dimLS];
  double Uni, Unj, Uwi, Uwj;
  double Phia;
  double srho1, srho2, srhod;
  double unroe, rroe;
  double uroe;
  int k;

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;
    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS ? LSS->edgeIntersectsStructure(0.0,l) : false;
    int iCovered, jCovered;
    if(LSS && !LSS->withCracking()) {
      iCovered = LSS->fluidModel(0.0,i);
      jCovered = LSS->fluidModel(0.0,j);
    } else
      iCovered = jCovered = 0; //when i(j)Covered = 0, we have valid Phi (and U) on this node.

    if(iCovered && jCovered) continue;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }
    for (k=0; k<dimLS; ++k) {
      ddPij[k] = dx[0]*dPhidx[i][k] + dx[1]*dPhidy[i][k] + dx[2]*dPhidz[i][k];
      ddPji[k] = dx[0]*dPhidx[j][k] + dx[1]*dPhidy[j][k] + dx[2]*dPhidz[j][k];
    }

    if(!iCovered && !jCovered) { //the usual linear reconstruction
      if(intersect) {
        for(k=0; k<dim; k++) {
          Vi[k] = V[i][k];
          Vj[k] = V[j][k];
        }
        for(k=0; k<dimLS; k++) {
          Pi[k] = Phi[i][k];
          Pj[k] = Phi[j][k];    
        }
      }
      else {
        recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
        recFcnLS->compute(Phi[i], ddPij, Phi[j], ddPji, Pi, Pj);
      }
    } else { //const reconstruction
      for(k=0; k<dim; k++) {
        Vi[k] = V[i][k];
        Vj[k] = V[j][k];
      }
      for(k=0; k<dimLS; k++) {
        Pi[k] = Phi[i][k];
        Pj[k] = Phi[j][k];
      }
    }

    Uni = Vi[1]*normal[l][0] + Vi[2]*normal[l][1] + Vi[3]*normal[l][2] - normalVel[l];
    Unj = Vj[1]*normal[l][0] + Vj[2]*normal[l][1] + Vj[3]*normal[l][2] - normalVel[l];
/*    if(intersect) {
      LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
      LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
      Uwi = normal[l]*resij.normVel;
      Uwj = normal[l]*resji.normVel; 
    } else Uwi = Uwj = 0.0;
*/

    // roe flux
    if(1/*!iCovered && !jCovered && !intersect*/)
      for (k=0; k<dimLS; k++){

        // Original
        /*if(fabs(Pj[k]-Pi[k])<1.e-12*fabs(Pi[k]) || Pj[k]==Pi[k])
          uroe     = 0.5*(Unj + Uni);
        else
          uroe     = (Pj[k]*Unj - Pi[k]*Uni)/(Pj[k]-Pi[k]);
        Phia = Uni*Pi[k] + Unj*Pj[k] - fabs(uroe)*(Pj[k]-Pi[k]);
        PhiF[i][k] += 0.5*Phia;
        PhiF[j][k] -= 0.5*Phia;
        */

        double uav = 0.5*(Uni+Unj);
        if (uav > 0.0) {
          Phia = Pi[k]*uav;
        } else if (uav < 0.0) {
          Phia = Pj[k]*uav;
        } 

        PhiF[i][k] += Phia;
        PhiF[j][k] -= Phia;
      }
    else {
      if(!iCovered) {
        double sign = -1.0;
        if (fluidId[i] != 0)
          sign = 1.0;

        for(k=0; k<dimLS; k++)
          PhiF[i][k] += 0.5*Uni* Pi[k];//*/sign;
      }
      if(!jCovered) {
        double sign = -1.0;
        if (fluidId[j] != 0)
          sign = 1.0;

        for(k=0; k<dimLS; k++)
          PhiF[j][k] -= 0.5*Unj* Pj[k];//*/ sign;
      }    
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState,
                                              Vec<double> &irey, SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A)
{

  int k;
  double edgeirey;

  double dfdUi[neq*neq], dfdUj[neq*neq];

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  double length;

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    edgeirey = 0.5*(irey[i]+irey[j]);
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if(fluxFcn){
      fluxFcn[BC_INTERNAL]->computeJacobians(length, edgeirey, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj);

      if (masterFlag[l]) {
        Scalar *Aii = A.getElem_ii(i);
        Scalar *Ajj = A.getElem_ii(j);

        for (k=0; k<neq*neq; ++k) {
          Aii[k] += dfdUi[k];
          Ajj[k] -= dfdUj[k];
        }
      }

      Scalar *Aij = A.getElem_ij(l);
      Scalar *Aji = A.getElem_ji(l);

      if (Aij && Aij) {

        double voli = 1.0 / ctrlVol[i];
        double volj = 1.0 / ctrlVol[j];
        for (k=0; k<neq*neq; ++k) {
          Aij[k] += dfdUj[k] * voli;
          Aji[k] -= dfdUi[k] * volj;
        }
      }

    }

  }
}
//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState,
                                              Vec<double> &irey, SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A,
                                              int *nodeType)
{

  /* in this function, rhs has already the values extrapolated at the inlet nodes
   * if we are in the case of water simulations
   * we are computing the jacobian matrix
   */
  int k,m;
  Scalar *Aii;
  Scalar *Ajj;
  Scalar *Aij;
  Scalar *Aji;
  double edgeirey, length;

  double dfdUi[neq*neq], dfdUj[neq*neq];
  bool atInleti, atInletj;

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    edgeirey = 0.5*(irey[i]+irey[j]);
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if(nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED ||
       nodeType[i] == BC_INLET_MOVING ||nodeType[i] == BC_OUTLET_MOVING)
       atInleti = true;
    else
       atInleti = false;

    if(nodeType[j] == BC_INLET_FIXED || nodeType[j] == BC_OUTLET_FIXED ||
       nodeType[j] == BC_INLET_MOVING ||nodeType[j] == BC_OUTLET_MOVING)
       atInletj = true;
    else
      atInletj = false;

    fluxFcn[BC_INTERNAL]->computeJacobians(length, edgeirey, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj);

/* first case: the two nodes are interior nodes
 * then the routine remains the same as usual
 * ie, fluxes are added to both jacobians (Aii and Ajj) and both crossed jacobians(Aij and Aji)
 * and the rhs need no changes.
 */
    if (!atInleti && !atInletj){

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

        for (k=0; k<neq*neq; ++k) {
          Aii[k] += dfdUi[k];
          Ajj[k] -= dfdUj[k];
        }
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);
      
      if (Aij && Aji) {
	double voli = 1.0 / ctrlVol[i];
	double volj = 1.0 / ctrlVol[j];
	for (k=0; k<neq*neq; ++k) {
	  Aij[k] += dfdUj[k] * voli;
	  Aji[k] -= dfdUi[k] * volj;
	}
      }

     }
/* second case: node i is an interior node, but j is an inlet node
 * the routine remains the same for Aii and Aij (which represents the influence of j on i)
 * the routine changes for Ajj and Aji. Aji is set to 0.0 since there is no influence of i on j
 * (value at j will be extrapolated), and Ajj is set to 1.0 (in the rhs term, the corresponding term
 * will be exactly the value it should take!)
 * the rhs for j need to be changed (this has been done previously in recomputeRHS)
 */
    else if (!atInleti && atInletj){

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

        for (k=0; k<neq*neq; k++)
          Aii[k] += dfdUi[k];

        for (k=0; k<neq; k++)
          Ajj[k*(neq+1)] = 1.0*ctrlVol[j];
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);

      if (Aij && Aji) {
        double voli = 1.0 / ctrlVol[i];
        for (k=0; k<neq*neq; k++) {
          Aij[k] += dfdUj[k] * voli;
          Aji[k] = 0.0;
        }
      }

    }

/* third case: same as case 2, but i and j have inversed roles
 */
    else if (atInleti && !atInletj){

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

        for (k=0; k<neq*neq; k++)
          Ajj[k] -= dfdUj[k];

        for (k=0; k<neq; k++)
          Aii[k*(neq+1)] = 1.0 * ctrlVol[i];
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);

      if (Aij && Aji) {
        double volj = 1.0 / ctrlVol[j];
        for (k=0; k<neq*neq; k++){
          Aij[k] = 0.0;
          Aji[k] -= dfdUi[k] * volj;
        }
      }

   }

/* fourth case: both nodes i and j are inletNodes
 * the routine is different for both of them, as both have no influence on
 * each other, and they are not subject to the influence of any nodes.
 * Their jacobians are set to 1.0 and the crossed jacobians to 0.0
 * the rhs for i and j need to be changed (this has been done previously in
 * recomputeRHS )
 */

    else if (atInleti && atInletj){

      if (masterFlag[l]){
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

        for (k=0; k<neq; k++){
          Aii[k*(neq+1)] = 1.0 *ctrlVol[i];
          Ajj[k*(neq+1)] = 1.0 *ctrlVol[j];
        }
      }

        // both Aij and Aji are kept to 0
        // and no change in the rhs is needed since they are both extrapolated!

    }

  }
}

//----------------------------------------------------------
template<int dim, class Scalar, int neq, int dimLS>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                              FluxFcn **fluxFcn, GeoState &geoState,
                                              NodalGrad<dim> &ngrad, NodalGrad<dimLS> &ngradLS,
                                              SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, FluidSelector &fluidSelector,
                                              Vec<int> &fluidId)
{
  // it is assumed that dim=5, ie no turbulence possible
  int k,m,q;

  double gradphi[3];
  double gphii[3];
  double gphij[3];
  double length;

  double dfdUi[neq*neq], dfdUj[neq*neq];
  double dfdUk[neq*neq], dfdUl[neq*neq];
  double Vi[2*dim], Vj[2*dim];
  double Wi[2*dim], Wj[2*dim];
  double dWidWi[neq*neq], dWjdWj[neq*neq],dWidWj[neq*neq],dWjdWi[neq*neq];
  double dWidUi[neq*neq], dWjdUj[neq*neq],dWidUj[neq*neq],dWjdUi[neq*neq];
  double dUidUi[neq*neq], dUjdUj[neq*neq],dUidUj[neq*neq],dUjdUi[neq*neq];
  double dii[neq*neq],djj[neq*neq];

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,dimLS>& dPdx = ngradLS.getX();
  SVec<double,dimLS>& dPdy = ngradLS.getY();
  SVec<double,dimLS>& dPdz = ngradLS.getZ();

  double jacii[neq*neq], jacij[neq*neq], jacji[neq*neq], jacjj[neq*neq];
  double jaciiu[neq*neq], jaciju[neq*neq], jacjiu[neq*neq], jacjju[neq*neq];

  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

  riemann.reset(0);

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (fluidId[i]==fluidId[j]) {
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, fluidId[i]);
     } else {
      //ngradLS returns nodal gradients of phi
		  //ngradLS returns nodal gradients of phi
      // need fluidSelector to determine which level set to look at knowing which two fluids are considered at this interface
      int lsdim = fluidSelector.getLevelSetDim(fluidId[i],fluidId[j]);
      gphii[0] = -dPdx[i][lsdim];
      gphii[1] = -dPdy[i][lsdim];
      gphii[2] = -dPdz[i][lsdim];
      gphij[0] = -dPdx[j][lsdim];
      gphij[1] = -dPdy[j][lsdim];
      gphij[2] = -dPdz[j][lsdim];
      for (k=0; k<3; k++)
        gradphi[k] = 0.5*(gphii[k]+gphij[k]);
      double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
      for (k=0; k<3; k++)
        gradphi[k] /= normgradphi;

      for(k=0; k<5; k++){
        Vi[k] = V[i][k];
        Vj[k] = V[j][k];
        Vi[k+5] = Vi[k];
        Vj[k+5] = Vj[k];
      }
      riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                     Wi,Wj,i,j,l,dx,false);
      riemann.computeRiemannJacobian(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                     Wi,Wj,i,j,l,dx,dWidWi, dWidWj,dWjdWi, dWjdWj );
      varFcn->postMultiplyBydVdU(Vi, dWidWi, dWidUi,fluidId[i]);
      varFcn->postMultiplyBydVdU(Vi, dWjdWi, dWjdUi,fluidId[i]);
      varFcn->postMultiplyBydVdU(Vj, dWidWj, dWidUj,fluidId[j]);
      varFcn->postMultiplyBydVdU(Vj, dWjdWj, dWjdUj,fluidId[j]);		
      
      varFcn->preMultiplyBydUdV(Wi, dWidUi, dUidUi,fluidId[i]);
      varFcn->preMultiplyBydUdV(Wj, dWjdUi, dUjdUi,fluidId[j]);
      varFcn->preMultiplyBydUdV(Wi, dWidUj, dUidUj,fluidId[i]);
      varFcn->preMultiplyBydUdV(Wj, dWjdUj, dUjdUj,fluidId[j]);
      
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Wi, dfdUi, dfdUk, fluidId[i]);
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wj, Vj, dfdUl, dfdUj, fluidId[j]);
    }
    
    Scalar *Aii;      
    Scalar *Ajj;
    if (masterFlag[l]) {
      Aii = A.getElem_ii(i);
      Ajj = A.getElem_ii(j);
      
      if (fluidId[i]==fluidId[j]) {
	for (k=0; k<neq*neq; ++k) {
	  Aii[k] += (dfdUi[k]);
	  Ajj[k] -= (dfdUj[k]);
	}
      } else {
	DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUi, 0, &dii,0);
	DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUj, 0, &djj,0);
	
	for (k=0; k<neq*neq; ++k) {
	  Aii[k] += (dfdUi[k]+dii[k]);
	  Ajj[k] -= (dfdUj[k]+djj[k]);
	}
      }
    }
    
    Scalar *Aij = A.getElem_ij(l);
    Scalar *Aji = A.getElem_ji(l);
    if (Aij && Aji) {
      double voli = 1.0 / ctrlVol[i];
      double volj = 1.0 / ctrlVol[j];			
      
      if (fluidId[i]==fluidId[j]) {
	for (k=0; k<neq*neq; ++k) {
	  Aij[k] += (dfdUj[k])*voli;
	  Aji[k] -= (dfdUi[k])*volj;
	}
      } else {
	DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUj,0,&dii,0);
	DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUi,0,&djj,0);
	for (k=0; k<neq*neq; ++k) {
	  Aij[k] += (dii[k]) * voli;
	  Aji[k] -= (djj[k]) * volj;
	}
      }
    }

    //riemann.resetInterfacialW(l);
  }
}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq, int dimLS>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                              FluxFcn **fluxFcn, GeoState &geoState,
                                              NodalGrad<dim> &ngrad, NodalGrad<dimLS> &ngradLS,
                                              SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, FluidSelector &fluidSelector,
                                              Vec<int> &fluidId, int *nodeType)
{
  /* in this function, rhs has already the values extrapolated at the inlet nodes
   * if we are in the case of water simulations
   * we are computing the jacobian matrix
   */

  int k,m;
  Scalar *Aii;
  Scalar *Ajj;
  Scalar *Aij;
  Scalar *Aji;

  double dfdUi[neq*neq], dfdUj[neq*neq];
  double dfdUk[neq*neq], dfdUl[neq*neq];
  double Vi[2*dim], Vj[2*dim];
  double Wi[2*dim], Wj[2*dim];
  bool atInleti, atInletj;
  double length;

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,dimLS>& dPdx = ngradLS.getX();
  SVec<double,dimLS>& dPdy = ngradLS.getY();
  SVec<double,dimLS>& dPdz = ngradLS.getZ();

  VarFcn *varFcn;

  double dWidWi[neq*neq], dWjdWj[neq*neq],dWidWj[neq*neq],dWjdWi[neq*neq];
  double dWidUi[neq*neq], dWjdUj[neq*neq],dWidUj[neq*neq],dWjdUi[neq*neq];
  double dUidUi[neq*neq], dUjdUj[neq*neq],dUidUj[neq*neq],dUjdUi[neq*neq];
  double dii[neq*neq],djj[neq*neq];
  double dij[neq*neq],dji[neq*neq];

  assert(false);

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    double gradphi[3];
    double gphii[3];
    double gphij[3];

    if(nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED ||
       nodeType[i] == BC_INLET_MOVING ||nodeType[i] == BC_OUTLET_MOVING)
       atInleti = true;
    else
       atInleti = false;

    if(nodeType[j] == BC_INLET_FIXED || nodeType[j] == BC_OUTLET_FIXED ||
       nodeType[j] == BC_INLET_MOVING ||nodeType[j] == BC_OUTLET_MOVING)
       atInletj = true;
    else
      atInletj = false;

    bool sameFluid = (fluidId[i]==fluidId[j]);

    if (sameFluid) {
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, fluidId[i]);
      riemann.resetInterfacialW(l);
    } else {
      //ngradLS returns nodal gradients of phi
      // need fluidSelector to determine which level set to look at knowing which two fluids are considered at this interface
      int lsdim = fluidSelector.getLevelSetDim(fluidId[i],fluidId[j]);
      gphii[0] = -dPdx[i][lsdim];
      gphii[1] = -dPdy[i][lsdim];
      gphii[2] = -dPdz[i][lsdim];
      gphij[0] = -dPdx[j][lsdim];
      gphij[1] = -dPdy[j][lsdim];
      gphij[2] = -dPdz[j][lsdim];
      for (k=0; k<3; k++)
        gradphi[k] = 0.5*(gphii[k]+gphij[k]);
      double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
      for (k=0; k<3; k++)
        gradphi[k] /= normgradphi;

      for(k=0; k<5; k++){
        Vi[k] = V[i][k];
        Vj[k] = V[j][k];
        Vi[k+5] = Vi[k];
        Vj[k+5] = Vj[k];
      }
      varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
      riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                     Wi,Wj,i,j,l,dx,false);

      riemann.computeRiemannJacobian(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                     Wi,Wj,i,j,l,dx,dWidWi, dWidWj,dWjdWi, dWjdWj );
      varFcn->postMultiplyBydVdU(Vi, dWidWi, dWidUi,fluidId[i]);
      varFcn->postMultiplyBydVdU(Vi, dWjdWi, dWjdUi,fluidId[i]);
      varFcn->postMultiplyBydVdU(Vj, dWidWj, dWidUj,fluidId[j]);
      varFcn->postMultiplyBydVdU(Vj, dWjdWj, dWjdUj,fluidId[j]);		
      
      varFcn->preMultiplyBydUdV(Wi, dWidUi, dUidUi,fluidId[i]);
      varFcn->preMultiplyBydUdV(Wj, dWjdUi, dUjdUi,fluidId[j]);
      varFcn->preMultiplyBydUdV(Wi, dWidUj, dUidUj,fluidId[i]);
      varFcn->preMultiplyBydUdV(Wj, dWjdUj, dUjdUj,fluidId[j]);
      
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Wi, dfdUi, dfdUk, fluidId[i]);
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wj, Vj, dfdUl, dfdUj, fluidId[j]);
      DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUi, 0, &dii,0);
      DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUj, 0, &djj,0);
      DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUj,0,&dij,0);
      DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUi,0,&dji,0);
    }

/* first case: the two nodes are interior nodes
 * then the routine remains the same as usual
 * ie, fluxes are added to both jacobians (Aii and Ajj) and both crossed jacobians(Aij and Aji)
 * and the rhs need no changes.
 */
    if (!atInleti && !atInletj){

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

	if (sameFluid) {
	  for (k=0; k<neq*neq; ++k) {
	    Aii[k] += dfdUi[k];
	    Ajj[k] -= dfdUj[k];
	  }
	} else {
	  for (k=0; k<neq*neq; ++k) {
	    Aii[k] += dfdUi[k]+dii[k];
	    Ajj[k] -= dfdUj[k]+djj[k];
	  }
	}
	  
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);
      
      if (Aij && Aji) {
	double voli = 1.0 / ctrlVol[i];
	double volj = 1.0 / ctrlVol[j];
	if (sameFluid) {
	  for (k=0; k<neq*neq; ++k) {
	    Aij[k] += dfdUj[k] * voli;
	    Aji[k] -= dfdUi[k] * volj;
	  }
	} else {
	  for (k=0; k<neq*neq; ++k) {
	    Aij[k] += (dfdUj[k]+dij[k]) * voli;
	    Aji[k] -= (dfdUi[k]+dji[k]) * volj;
	  }
	}
      }

     }
/* second case: node i is an interior node, but j is an inlet node
 * the routine remains the same for Aii and Aij (which represents the influence of j on i)
 * the routine changes for Ajj and Aji. Aji is set to 0.0 since there is no influence of i on j
 * (value at j will be extrapolated), and Ajj is set to 1.0 (in the rhs term, the corresponding term
 * will be exactly the value it should take!)
 * the rhs for j need to be changed (this has been done previously in recomputeRHS)
 */
    else if (!atInleti && atInletj) {

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

	if (sameFluid) {
	  for (k=0; k<neq*neq; k++)
	    Aii[k] += dfdUi[k];
	} else {
	  for (k=0; k<neq*neq; k++)
	    Aii[k] += dfdUi[k]+dii[k];	  
	}

        for (k=0; k<neq; k++)
          Ajj[k*(neq+1)] = 1.0*ctrlVol[j];
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);

      if (Aij && Aji) {
        double voli = 1.0 / ctrlVol[i];
	if (sameFluid) {
	  for (k=0; k<neq*neq; k++) {
	    Aij[k] += dfdUj[k] * voli;
	    Aji[k] = 0.0;
	  }
	} else {
	  for (k=0; k<neq*neq; k++) {
	    Aij[k] += (dfdUj[k]+dij[k]) * voli;
	    Aji[k] = 0.0;
	  }
	}
      }

    }
/* third case: same as case 2, but i and j have inversed roles
 */
    else if (atInleti && !atInletj){

      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

	if (sameFluid) {
	  for (k=0; k<neq*neq; k++)
	    Ajj[k] -= dfdUj[k];
	} else {
	  for (k=0; k<neq*neq; k++)
	    Ajj[k] -= dfdUj[k]+djj[k];
	}

        for (k=0; k<neq; k++)
          Aii[k*(neq+1)] = 1.0 * ctrlVol[i];
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);

      if (Aij && Aji) {
        double volj = 1.0 / ctrlVol[j];
	if (sameFluid) {
	  for (k=0; k<neq*neq; k++){
	    Aij[k] = 0.0;
	    Aji[k] -= dfdUi[k] * volj;
	  }
	} else {
	  for (k=0; k<neq*neq; k++){
	    Aij[k] = 0.0;
	    Aji[k] -= (dfdUi[k]+dji[k]) * volj;
	  }
	}
      }

   }

/* fourth case: both nodes i and j are inletNodes
 * the routine is different for both of them, as both have no influence on
 * each other, and they are not subject to the influence of any nodes.
 * Their jacobians are set to 1.0 and the crossed jacobians to 0.0
 * the rhs for i and j need to be changed (this has been done previously in
 * recomputeRHS )
 */

    else if (atInleti && atInletj){

      if (masterFlag[l]){
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);

        for (k=0; k<neq; k++){
          Aii[k*(neq+1)] = 1.0 *ctrlVol[i];
          Ajj[k*(neq+1)] = 1.0 *ctrlVol[j];
        }
      }

        // both Aij and Aji are kept to 0
        // and no change in the rhs is needed since they are both extrapolated!

    }
  }
}

template<class Scalar,int dim,int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                             FluxFcn** fluxFcn,
                                             GeoState& geoState, SVec<double,3>& X,
                                             SVec<double,dim>& V, Vec<double>& ctrlVol,
                                             LevelSetStructure &LSS,
                                             Vec<int> &fluidId, int Nriemann, SVec<double,3> *Nsbar,
                                             GenMat<Scalar,neq>& A,Vec<double>& irey) {

  int farfieldFluid = 0; 

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  double Vi[2*dim], Vj[2*dim];
  double Wstar[2*dim];
  double dUdU[dim*dim],dfdUi[dim*dim],dfdUj[dim*dim],dkk[dim*dim],dWdW[dim*dim],dWdU[dim*dim];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;
  Scalar * Aii,*Ajj,*Aij,*Aji;
  Vec3D normalDir;
  int k;
 
  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,l);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);
    
    double edgeirey = 0.5*(irey[i]+irey[j]);
    if( !iActive && !jActive ) {
      continue;
    }

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    for(k=0; k<dim; k++) {
      Vi[k] = V[i][k];
      Vj[k] = V[j][k];
      Vi[k+dim] = Vi[k];
      Vj[k+dim] = Vj[k];
    }
 
    if (!intersect) {  // same fluid

      if(!(iActive && jActive)) {
        fprintf(stderr,"Really odd too... (-(%d),-(%d): intersect = %d; iSwept = %d, jSwept = %d, iOccluded = %d, jOcculded = %d, model = %d/%d.\n",
                iActive, jActive, intersect, LSS.isSwept(0.0,i), LSS.isSwept(0.0,j), 
                LSS.isOccluded(0.0,i), LSS.isOccluded(0.0,j), LSS.fluidModel(0.0,i), LSS.fluidModel(0.0,j));
        continue;
      }

      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Vj, dfdUi, dfdUj, fluidId[i]);
      Aii = A.getElem_ii(i);
      Ajj = A.getElem_ii(j);

      if (masterFlag[l]) {
        for (k=0; k<dim*dim; ++k) {
          Aii[k] += dfdUi[k];
	  Ajj[k] -= dfdUj[k];
        }
      }

      Aij = A.getElem_ij(l);
      Aji = A.getElem_ji(l);
      if (Aij && Aji) {
        double voli = 1.0 / ctrlVol[i];
        double volj = 1.0 / ctrlVol[j];
 
        for (k=0; k<dim*dim; ++k) {
	  Aij[k] += dfdUj[k] * voli;
	  Aji[k] -= dfdUi[k] * volj;
	}
      }
    } else if (masterFlag[l]) {

      if(iActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);

        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      resij.gradPhi;
            else                                normalDir = -1.0*resij.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[i]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());

        riemann.computeFSIRiemannSolution(Vi,resij.normVel,normalDir,varFcn,Wstar,j,fluidId[i]);
        riemann.computeFSIRiemannJacobian(Vi,resij.normVel,normalDir,varFcn,Wstar,j,dWdW,fluidId[i]);
        
        varFcn->postMultiplyBydVdU(Wstar, dWdW, dWdU,fluidId[i]);
        varFcn->preMultiplyBydUdV(Vi, dWdU, dUdU,fluidId[i]);

        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Wstar, dfdUi, dfdUj, fluidId[i],false);
        DenseMatrixOp<double, dim, dim*dim>::applyToDenseMatrix(&dfdUj,0,&dUdU, 0, &dkk,0);
        Aii = A.getElem_ii(i);
        for (k=0; k<dim*dim; ++k) {
          Aii[k] += dfdUi[k]+dkk[k];
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);

        switch (Nriemann) {
          case 0: //structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      resji.gradPhi;
            else                                normalDir = -1.0*resji.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(fluidId[j]==farfieldFluid)       normalDir =      Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            else                                normalDir = -1.0*Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());
  
        riemann.computeFSIRiemannSolution(Vj,resji.normVel,normalDir,varFcn,Wstar,i,fluidId[j]);
        riemann.computeFSIRiemannJacobian(Vj,resji.normVel,normalDir,varFcn,Wstar,i,dWdW,fluidId[j]);
        
        varFcn->postMultiplyBydVdU(Wstar, dWdW, dWdU,fluidId[j]);
        varFcn->preMultiplyBydUdV(Vj, dWdU, dUdU,fluidId[j]);
  
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wstar,Vj, dfdUi, dfdUj, fluidId[j],false); 
        DenseMatrixOp<double, dim, dim*dim>::applyToDenseMatrix(&dfdUi,0,&dUdU, 0, &dkk,0);
        Ajj = A.getElem_ii(j);
        for (k=0; k<dim*dim; ++k) {
          Ajj[k] -= dfdUj[k]+dkk[k];
        }
      }
    }
  }
}

template<class Scalar,int dim, int dimLS,int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn,
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, 
                                     LevelSetStructure& LSS, Vec<int> &fluidId,
                                     int Nriemann, SVec<double,3>* Nsbar, FluidSelector &fluidSelector,
                                     NodalGrad<dimLS>& ngradLS,Vec<double>& ctrlVol,
                                     GenMat<Scalar,neq>& A) {

  // ------------------------------------------------
  //  Preparation -- General Info. 
  // ------------------------------------------------
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  double Vi[2*dim], Vj[2*dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;
  double dfdUi[neq*neq],dfdUj[neq*neq],dkk[neq*neq],dUidUi[neq*neq],dUidUj[neq*neq],dUjdUj[neq*neq],dUjdUi[neq*neq],dii[neq*neq],djj[neq*neq];
  double dWidUi[neq*neq],dWjdUj[neq*neq],dWidUj[neq*neq], dWjdUi[neq*neq]; 
  int ierr=0,k;
  double dWidWi[neq*neq],dWjdWj[neq*neq],dWidWj[neq*neq], dWjdWi[neq*neq],dfdUk[neq*neq],dfdUl[neq*neq]; 
  double dWdW[neq*neq],dWdU[neq*neq], dUdW[neq*neq],dUdU[neq*neq];
  riemann.reset(0);

  // ------------------------------------------------
  //  Preparation -- Fluid/Structure Part
  // ------------------------------------------------
  int farfieldFluid = 0; //assume that intersector assigns Id=0 for outside fluid.
  double Wstar[2*dim]; //FS Riemann solution
  Vec3D normalDir; //Normal direction used for setting the Riemann problem

  // ------------------------------------------------
  //  Preparation -- Fluid/Fluid Part
  // ------------------------------------------------
  SVec<double,dimLS>& dPdx = ngradLS.getX();
  SVec<double,dimLS>& dPdy = ngradLS.getY();
  SVec<double,dimLS>& dPdz = ngradLS.getZ();
  double Wi[2*dim], Wj[2*dim]; //FF Riemann solution
  double gradphi[3], gphii[3], gphij[3];

  // ------------------------------------------------
  //  THE MAIN EDGE LOOP...
  // ------------------------------------------------
  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,l);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);

    if(!iActive && !jActive) continue; //this edge is inside a solid body!

    for(k=0; k<dim; k++){
      Vi[k] = V[i][k];
      Vj[k] = V[j][k];
      Vi[k+dim] = Vi[k];
      Vj[k+dim] = Vj[k];
    }
    // ------------------------------------------------
    //  Reconstruction without crossing the FS interface.
    // ------------------------------------------------
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    // --------------------------------------------------------
    //  Compute the flux along this edge.
    //    Step 1. If e(i,j) intersects the structure -> FS flux
    //    Step 2. If Idi!=Idj -> FF flux
    //    Step 3. Otherwise, the usual single-phase flux. 
    // --------------------------------------------------------
    if(intersect) {

      // for node i
      if(iActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
        switch (Nriemann) { // normal should point to this node (i). TODO: better to compute a projection to determine the sign.
          case 0: //structure normal
            if(LSS.fluidModel(0.0,i)==farfieldFluid)  normalDir =      resij.gradPhi;
            else                                     normalDir = -1.0*resij.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(LSS.fluidModel(0.0,i)==farfieldFluid)  normalDir =      Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            else                                     normalDir = -1.0*Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }

        riemann.computeFSIRiemannSolution(Vi,resij.normVel,normalDir,varFcn,Wstar,j,fluidId[i]);
        riemann.computeFSIRiemannJacobian(Vi,resij.normVel,normalDir,varFcn,Wstar,j,dWdW,fluidId[i]);
        
        varFcn->postMultiplyBydVdU(Wstar, dWdW, dWdU,fluidId[i]);
        varFcn->preMultiplyBydUdV(Vi, dWdU, dUdU,fluidId[i]);
  
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi,Wstar, dfdUi, dfdUj, fluidId[i],false); 
        DenseMatrixOp<double, dim, dim*dim>::applyToDenseMatrix(&dfdUj,0,&dUdU, 0, &dkk,0);
        Scalar* Aii = A.getElem_ii(i);
        for (k=0; k<dim*dim; ++k) {
          Aii[k] += dfdUi[k]+dkk[k];
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
        switch (Nriemann) {
          case 0: //structure normal
            if(LSS.fluidModel(0.0,j)==farfieldFluid)  normalDir =      resji.gradPhi;
            else                                     normalDir = -1.0*resji.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            if(LSS.fluidModel(0.0,j)==farfieldFluid)  normalDir =      Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            else                                     normalDir = -1.0*Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        riemann.computeFSIRiemannSolution(Vj,resji.normVel,normalDir,varFcn,Wstar,i,fluidId[j]);
        riemann.computeFSIRiemannJacobian(Vj,resji.normVel,normalDir,varFcn,Wstar,i,dWdW,fluidId[j]);
        
        varFcn->postMultiplyBydVdU(Wstar, dWdW, dWdU,fluidId[j]);
        varFcn->preMultiplyBydUdV(Vj, dWdU, dUdU,fluidId[j]);
  
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wstar,Vj, dfdUi, dfdUj, fluidId[j],false); 
        DenseMatrixOp<double, dim, dim*dim>::applyToDenseMatrix(&dfdUi,0,&dUdU, 0, &dkk,0);
        Scalar* Ajj = A.getElem_ii(j);
        for (k=0; k<dim*dim; ++k) {
          Ajj[k] -= dfdUj[k]+dkk[k];
        }
      }
    }

    else {
  
      if(fluidId[i]!=fluidId[j]) { //NOTE: It's NOT equivalent with checking Phi_i x Phi_j < 0!
        //ngradLS returns nodal gradients of primitive phi
        // need fluidSelector to determine which level set to look at knowing which two fluids are considered at this interface
        int lsdim = fluidSelector.getLevelSetDim(fluidId[i],fluidId[j],locToGlobNodeMap[i]+1,locToGlobNodeMap[j]+1);
        gphii[0] = -dPdx[i][lsdim];
        gphii[1] = -dPdy[i][lsdim];
        gphii[2] = -dPdz[i][lsdim];
        gphij[0] = -dPdx[j][lsdim];
        gphij[1] = -dPdy[j][lsdim];
        gphij[2] = -dPdz[j][lsdim];
        for (int k=0; k<3; k++)
          gradphi[k] = 0.5*(gphii[k]+gphij[k]);
        double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
        for (int k=0; k<3; k++)
          gradphi[k] /= normgradphi;

        riemann.computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                       Wi,Wj,i,j,l,dx,false);
        riemann.computeRiemannJacobian(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
                                       Wi,Wj,i,j,l,dx,dWidWi, dWidWj,dWjdWi, dWjdWj );
        varFcn->postMultiplyBydVdU(Vi, dWidWi, dWidUi,fluidId[i]);
        varFcn->postMultiplyBydVdU(Vi, dWjdWi, dWjdUi,fluidId[i]);
        varFcn->postMultiplyBydVdU(Vj, dWidWj, dWidUj,fluidId[j]);
        varFcn->postMultiplyBydVdU(Vj, dWjdWj, dWjdUj,fluidId[j]);		
      
        varFcn->preMultiplyBydUdV(Wi, dWidUi, dUidUi,fluidId[i]);
        varFcn->preMultiplyBydUdV(Wj, dWjdUi, dUjdUi,fluidId[j]);
        varFcn->preMultiplyBydUdV(Wi, dWidUj, dUidUj,fluidId[i]);
        varFcn->preMultiplyBydUdV(Wj, dWjdUj, dUjdUj,fluidId[j]);
      
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Wi, dfdUi, dfdUk, fluidId[i]);
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wj, Vj, dfdUl, dfdUj, fluidId[j]);
      } else {
   
        fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, fluidId[i]);
      }

      Scalar *Aii;      
      Scalar *Ajj;
      if (masterFlag[l]) {
        Aii = A.getElem_ii(i);
        Ajj = A.getElem_ii(j);
      
        if (fluidId[i]==fluidId[j]) {
  	  for (k=0; k<neq*neq; ++k) {
	    Aii[k] += (dfdUi[k]);
	    Ajj[k] -= (dfdUj[k]);
	  }
        } else {
  	  DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUi, 0, &dii,0);
	  DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUj, 0, &djj,0);
	
	  for (k=0; k<neq*neq; ++k) {
	    Aii[k] += (dfdUi[k]+dii[k]);
	    Ajj[k] -= (dfdUj[k]+djj[k]);
	  }
        }
      }
    
      Scalar *Aij = A.getElem_ij(l);
      Scalar *Aji = A.getElem_ji(l);
      if (Aij && Aji) {
        double voli = 1.0 / ctrlVol[i];
        double volj = 1.0 / ctrlVol[j];			
      
        if (fluidId[i]==fluidId[j]) {
	  for (k=0; k<neq*neq; ++k) {
	    Aij[k] += (dfdUj[k])*voli;
	    Aji[k] -= (dfdUi[k])*volj;
	  }
        } else {
	  DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUk,0,&dUidUj,0,&dii,0);
	  DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUl,0,&dUjdUi,0,&djj,0);
	  for (k=0; k<neq*neq; ++k) {
	    Aij[k] += (dii[k]) * voli;
	    Aji[k] -= (djj[k]) * volj;
	  }
        }
      }
    }
  }
}

template<class Scalar, int dim, int dimLS>
void EdgeSet::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
						GeoState& geoState, SVec<double,3>& X,
						SVec<double,dim>& V, NodalGrad<dim>& ngrad,
						NodalGrad<dimLS> &ngradLS,
						EdgeGrad<dim>* egrad,
						Vec<double> &ctrlVol, SVec<double,dimLS>& Phi,
						GenMat<Scalar,dimLS> &A,LevelSetStructure* LSS)
{
  int k;
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx  = ngrad.getX();
  SVec<double,dim>& dVdy  = ngrad.getY();
  SVec<double,dim>& dVdz  = ngrad.getZ();
  // in this routine Phi denotes the "conservative phi" ie (rho*phi)
  SVec<double,dimLS>& dPhidx = ngradLS.getX();
  SVec<double,dimLS>& dPhidy = ngradLS.getY();
  SVec<double,dimLS>& dPhidz = ngradLS.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim];
  double ddPij[dimLS], ddPji[dimLS], Pi[2*dimLS], Pj[2*dimLS];
  double Uni, Unj;
  double Phia;
  double srho1, srho2, srhod;
  double unroe, rroe;
  double uroe,uroed[2];
  Scalar df[2];

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    
    int iCovered = LSS ? LSS->fluidModel(0.0,i) : 0;
    int jCovered = LSS ? LSS->fluidModel(0.0,j) : 0;
    
    if (iCovered && jCovered)
      continue;

    for(k=0; k<dim; k++){
      Vi[k] = V[i][k];
      Vj[k] = V[j][k];
      Vi[k+dim] = Vi[k];
      Vj[k+dim] = Vj[k];
      }
    for(k=0; k<dimLS; k++){
      Pi[k] = Phi[i][k];
      Pj[k] = Phi[j][k];
    }

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }
    
    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    
    Uni      = Vi[1]*normal[l][0]  +Vi[2]*normal[l][1]  +Vi[3]*normal[l][2] - normalVel[l];
    Unj      = Vj[1]*normal[l][0]  +Vj[2]*normal[l][1]  +Vj[3]*normal[l][2] - normalVel[l];
    //Roe averaged variables
    if (1/*!iCovered && !jCovered*/) {
      for (k = 0; k < dimLS; ++k) {

        df[0] = df[1] = 0.0;

        double uav = 0.5*(Uni+Unj);
        if (uav > 0)
	  df[0] = uav;
        else if (uav < 0)
	  df[1] = uav;
      
        if (masterFlag[l]) {
	  Scalar* Aii = A.getElem_ii(i);
	  Scalar* Ajj = A.getElem_ii(j);
	  Aii[k*dimLS+k] += df[0];
	  Ajj[k*dimLS+k] += -df[1];
        }
      
        Scalar* Aij = A.getElem_ij(l);
        Scalar* Aji = A.getElem_ji(l);
        if (Aij && Aji) {
  	  double voli = 1.0 / ctrlVol[i];
	  double volj = 1.0 / ctrlVol[j];
	  Aji[k*dimLS+k] += (-df[0])*volj;
	  Aij[k*dimLS+k] += (df[1])*voli;
        }
      }
    } else {
      if(!iCovered) {
        if (masterFlag[l]) {
	  Scalar* Aii = A.getElem_ii(i);
          for(k=0; k<dimLS; k++)
            Aii[k*dimLS+k] += 0.0;//0.5*Uni;
        }
      }
      if(!jCovered) {
        if (masterFlag[l]) {
	  Scalar* Ajj = A.getElem_ii(j);
          for(k=0; k<dimLS; k++)
            Ajj[k*dimLS+k] -= 0.0;//0.5*Unj;
        }
      }    
    }
  }
}

//------------------------------------------------------------------------------
template<int dimLS>
void EdgeSet::TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi)
{

  int tag = 1;
  for(int l=0; l<numEdges; l++){
    int i = ptr[l][0];
    int j = ptr[l][1];

    if(Phi[i][lsdim]*Phi[j][lsdim]<=0.0){
      Tag[i] = tag;
      Tag[j] = tag;
    }
  }

}

//------------------------------------------------------------------------------
/*
template<int dimLS>
void EdgeSet::TagInterfaceNodes(int lsdim, Vec<int> &Tag1, Vec<int> &Tag2, SVec<double,dimLS> &Phi, LevelSetStructure *LSS)
{
  int tag = 1;
  int i,j;
  for(int l=0; l<numEdges; l++) {
    i = ptr[l][0];
    j = ptr[l][1]; 

    if(Phi[i][lsdim]*Phi[j][lsdim]<=0.0){
      if(LSS->edgeIntersectsStructure(0.0, l))
        Tag1[i] = tag;
      else
        Tag2[i] = tag;
    } else { //just for debug. can be removed.
      if(LSS->edgeIntersectsStructure(0.0,l))
        fprintf(stderr,"BUG: (i,j) intersects bug phi[i]=%e, phi[j]=%e.\n", Phi[i][lsdim], Phi[j][lsdim]);
    }
  }
}
*/
//------------------------------------------------------------------------------
