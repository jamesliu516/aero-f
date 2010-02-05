#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

#include <Edge.h>
#include <BcDef.h>
#include <FluxFcnDescPerfectGas.h>
#include <FluxFcnDescWaterCompressible.h>
#include <FluxFcnDescGasInGas.h>
#include <FluxFcnDescLiquidInLiquid.h>
#include <FluxFcnDescGasInLiquid.h>
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
    locbeta = tprec.getBeta(locMach);

    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    idti[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
    idti[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;

    if(fet) vis = fet->computeViscousTimeStep(Xmid, Vmid)*S*S;
    idtv[i] += vis;
    idtv[j] += vis;

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
    locbeta = tprec.getBeta(locMach);
    double dLocMach = varFcn->computeDerivativeOfMachNumber(Vmid, dVmid, dMach);
    dLocbeta = tprec.getdBeta(locMach,dLocMach);

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
                              TimeLowMachPrec &tprec, Vec<double> &Phi, int subnum)
{
  double Phimid, Vmid[dim];
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

    if(Phi[i]*Phi[j]>=0.0){
      for (int k=0; k<dim; ++k)
        Vmid[k] = 0.5 * (V[i][k] + V[j][k]);
      Phimid = 0.5 * (Phi[i] + Phi[j]);

      u = varFcn->getVelocity(Vmid);
      a = varFcn->computeSoundSpeed(Vmid, Phimid);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(Vmid,Phimid);

      locbeta = tprec.getBeta(mach);
      beta2 = locbeta*locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
    }else{
      u = varFcn->getVelocity(V[i]);
      a = varFcn->computeSoundSpeed(V[i], Phi[i]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[i],Phi[i]);

      locbeta = tprec.getBeta(mach);
      beta2 = locbeta * locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;

      u = varFcn->getVelocity(V[j]);
      a = varFcn->computeSoundSpeed(V[j], Phi[j]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[j],Phi[j]);

      locbeta = tprec.getBeta(mach);
      beta2 = locbeta * locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
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

  for (int l=0; l<numEdges; ++l) {

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


template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, Vec<double> &Phi,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     NodalGrad<1>& ngradLS,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<double,dim>* interfaceFlux,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,1>& dPdx = ngradLS.getX();
  SVec<double,1>& dPdy = ngradLS.getY();
  SVec<double,1>& dPdz = ngradLS.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wi[2*dim], Wj[2*dim];
  double fluxi[dim], fluxj[dim];
  double gradphi[3];
  double gphii[3];
  double gphij[3];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    //if(Phi[i]*Phi[j] > 0.0) recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    //else                    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj, Phi[i], Phi[j]);

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], Phi[i], Phi[j]);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (Phi[i]*Phi[j] > 0.0) { 	// same fluid
      if (Phi[i] > 0.0 && Phi[j] > 0.0)  {
        fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);
      }
      if (Phi[i] < 0.0 && Phi[j] < 0.0)  {
        fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);
      }
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{			// interface
      //ngradLS returns nodal gradients of primitive phi
      gphii[0] = dPdx[i][0];
      gphii[1] = dPdy[i][0];
      gphii[2] = dPdz[i][0];
      gphij[0] = dPdx[j][0];
      gphij[1] = dPdy[j][0];
      gphij[2] = dPdz[j][0];
      for (int k=0; k<3; k++)
        gradphi[k] = 0.5*(gphii[k]+gphij[k]);
      double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
      for (int k=0; k<3; k++)
        gradphi[k] /= normgradphi;

      int epsi = 0;
      int epsj = 0;
      riemann.computeRiemannSolution(Vi,Vj,Phi[i],Phi[j],gradphi,varFcn,
                                    epsi,epsj,Wi,Wj,i,j);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Vi, Wi, fluxi, epsi);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Wj, Vj, fluxj, epsj);
      for (int k=0; k<dim; k++){
        fluxes[i][k] += fluxi[k];
        fluxes[j][k] -= fluxj[k];
      }
      // in order to check mass conservation
      if(interfaceFlux)
        for (int k=0; k<dim; k++){
          (*interfaceFlux)[i][k] += fluxi[k];
          (*interfaceFlux)[j][k] -= fluxj[k];
        }

    }
  }
  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim>& Wstarij, 
                                     SVec<double,dim>& Wstarji, LevelSetStructure &LSS,
                                     bool linRecAtInterface, int Nriemann, NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  double gphii[3];
  double gphij[3];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  SVec<double,dim> tempWstarij(Wstarij);
  SVec<double,dim> tempWstarji(Wstarji);
  if (it>0) { //if it>0 (i.e. not called from computeResidualNorm), clear Wstarij and Wstarji for update.
    Wstarij = 0.0;
    Wstarji = 0.0;
  }

  Vec3D normalDir; //normal direction fed to the Riemann solver. Related to the choice of "Nriemann".

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool iIsActive = LSS.isActive(0, i);
    bool jIsActive = LSS.isActive(0, j);
    bool intersect = LSS.edgeIntersectsStructure(0,i,j);

    if( !iIsActive && !jIsActive ) 
      continue;

    // ------------------------------------------------
    //  Reconstruction without crossing the interface.
    // ------------------------------------------------
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    if (iIsActive && jIsActive && !intersect)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.

    else { // now at interface

      if (!linRecAtInterface) // just set Vi = V[i], Vj = V[j]
        for(int k=0; k<dim; k++) {
          Vi[k] = V[i][k];
          Vj[k] = V[j][k];
        }
     
      else { // linRec at interface using Wstar        
        // Case 1: i and j are both active (but i-j intersects)
        if (iIsActive&&jIsActive) {
          double Vtemp[2*dim];
          if (tempWstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) Vi[k] = V[i][k]; 
          } else recFcn->compute(V[i], ddVij, tempWstarij[l], ddVji, Vi, Vtemp);
          if (tempWstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) Vj[k] = V[j][k];
          } else recFcn->compute(tempWstarji[l], ddVij, V[j], ddVji, Vtemp, Vj);
        }
        // Case 2: i is active, j is inactive
        else if (iIsActive)
          if (tempWstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(V[i], ddVij, tempWstarij[l], ddVji, Vi, Vj);
        // Case 3: i is inactive, j is active
        else if (jIsActive)
          if (tempWstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(tempWstarji[l], ddVij, V[j], ddVji, Vi, Vj);
      }
    }

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j]); //also checking reconstructed values across interface.

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }
    

    // --------------------------------------------------------
    //                   Compute fluxes
    // --------------------------------------------------------
    if (!intersect) {  // same fluid
      if (!masterFlag[l]) continue; //not a master edge

      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{// interface
      if (iIsActive) {
        LevelSetResult res = LSS.getLevelSetDataAtEdgeCenter(0.0, i, j);
        normalDir = (Nriemann) ? -1.0/(normal[l].norm())*normal[l] : res.gradPhi;
        riemann.computeFSIRiemannSolution(Vi,res.normVel,normalDir,varFcn,Wstar,j);

        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k]; 
        if (masterFlag[l]) { 
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi);
          for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
        }
      }
      if (jIsActive) {
        LevelSetResult res = LSS.getLevelSetDataAtEdgeCenter(0.0, j,i);
        normalDir = (Nriemann) ? 1.0/(normal[l].norm())*normal[l] : res.gradPhi;
        riemann.computeFSIRiemannSolution(Vj,res.normVel,res.gradPhi,varFcn,Wstar,i);
        
        if (it>0)
          for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj);
          for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
        }
      }
    }
  }
  return ierr;
}

//------------------------------------------------------------------------------
/* old flux solver that can only link to Michel's Intersector.
template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim>& Wstarij, 
                                     SVec<double,dim>& Wstarji, LevelSetStructure &LSS,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  double gphii[3];
  double gphij[3];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  SVec<double,dim> tempWstarij(Wstarij);
  SVec<double,dim> tempWstarji(Wstarji);
  if (it>0) { //if it>0 (i.e. not called from computeResidualNorm), clear Wstarij and Wstarji for update.
    Wstarij = 0.0;
    Wstarji = 0.0;
  }

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool iIsActive = LSS.isActive(0, i);
    bool jIsActive = LSS.isActive(0, j);

    if( !iIsActive && !jIsActive ) 
      continue;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

// Reconstruction without crossing the interface.
    if (iIsActive && jIsActive)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.
    else { // now at interface
      if (iIsActive)
        if (tempWstarij[l][0]<1e-8 && tempWstarij[l][4]<1e-8) {// no riemann sol. (first time-step)
          for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
//          fprintf(stderr,"linRec turned off for Node %d on Edge (%d->%d).\n",
//                         locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
        } else recFcn->compute(V[i], ddVij, tempWstarij[l], ddVji, Vi, Vj);
      if (jIsActive)
        if (tempWstarji[l][0]<1e-8 && tempWstarji[l][4]<1e-8) {
          for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
//          fprintf(stderr,"linRec turned off for Node %d on Edge (%d->%d).\n",
//                         locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
        } else recFcn->compute(tempWstarji[l], ddVij, V[j], ddVji, Vi, Vj);
    }

    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j]); //also checking reconstructed values across interface.

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (!LSS.edgeIntersectsStructure(0, i, j)) {  // same fluid
      if (!masterFlag[l]) continue;
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{// interface
      if (iIsActive) {
        LevelSetResult res = LSS.getLevelSetDataAtEdgeCenter(0.0, i, j);
        riemann.computeFSIRiemannSolution(Vi,res.normVel,res.gradPhi,varFcn,Wstar,j);
        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k]; 
        if (!masterFlag[l]) continue;
        fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi);
        for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
      }
      if (jIsActive) {
        LevelSetResult res = LSS.getLevelSetDataAtEdgeCenter(0.0, j,i);
        riemann.computeFSIRiemannSolution(Vj,res.normVel,res.gradPhi,varFcn,Wstar,i);
        if (it>0)
          for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
        if (!masterFlag[l]) continue;
        fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj);
        for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
      }
    }
  }

  return ierr;

}
*/
//------------------------------------------------------------------------------
// for FSF (fluid-shell-fluid).
template<int dim>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim>& Wstarij,
                                     SVec<double,dim>& Wstarji, LevelSetStructure &LSS, Vec<double> &nodeTag,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  double gphii[3];
  double gphij[3];
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  SVec<double,dim> tempWstarij(Wstarij);
  SVec<double,dim> tempWstarji(Wstarji);
  Wstarij = 0.0;
  Wstarji = 0.0;
  Vec3D vStar(10.0,0.0,0.0);
  Vec3D gradPhi(1.0,0.0,0.0);

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];
    int typei = (nodeTag[i]>0.0) ? 1 : -1;
    int typej = (nodeTag[j]>0.0) ? 1 : -1;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

// Reconstruction without crossing the interface.
    if ((typei*typej)>0)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.
    else { // now at interface
      if (tempWstarij[l][0]<1e-8 && tempWstarij[l][4]<1e-8) {// no riemann sol. (first time-step)
        for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
        fprintf(stderr,"linRec turned off for Node %d on Edge (%d->%d).\n",
                       locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
      } else {double tempVj[dim*2]; recFcn->compute(V[i], ddVij, tempWstarij[l], ddVji, Vi, tempVj);}

      if (tempWstarji[l][0]<1e-8 && tempWstarji[l][4]<1e-8) {
        for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
        fprintf(stderr,"linRec turned off for Node %d on Edge (%d->%d).\n",
                       locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
      } else {double tempVi[dim*2]; recFcn->compute(tempWstarji[l], ddVij, V[j], ddVji, tempVi, Vj);}
    }

/*    if ((typei*typej)>0)
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj); //Vi and Vj are reconstructed states.
    else for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
*/
    // check for negative pressure or density //
    if (!rshift)
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, V[i], V[j], nodeTag[i], nodeTag[j]); //also checking reconstructed values acrossinterface.

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (typei==typej) {  // same fluid
      if (typei<0) fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);
      if (typei>0) fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{                       // interface
//      fprintf(stderr,"edge(%d(%f),%d(%f)). LEFT: %e %e | RIGHT: %e %e\n", i,nodeTag[i],j,nodeTag[j], V[i][4], Vi[4], V[j][4], Vj[4]);
      if (nodeTag[i]<0.0)
        riemann.computeFSIRiemannSolution((int)nodeTag[i],Vi,vStar,gradPhi,varFcn,Wstar,j);
      else
        riemann.computeFSIRiemannSolution((int)nodeTag[i],Vi,vStar,-1.0*gradPhi,varFcn,Wstar,j);
      fprintf(stderr,"Vi = (%e %e %e %e %e), vStar = (%e %e %e %e %e). nodeTag[i]=%e.\n", Vi[0], Vi[1],Vi[2], Vi[3], Vi[4], Wstar[0], Wstar[1], Wstar[2], Wstar[3], Wstar[4], nodeTag[i]);
      for (int k=0; k<dim; k++) Wstarij[l][k] = Wstar[k]; //stores Wstar for later use.


/*      double area = normal[l].norm(); //area of c.v. surface
      area *= normal[l]*gradPhi/normal[l].norm(); //projected to structure surface
      LSS.totalForce[0] += Wstar[4]*gradPhi[0]*area;
      LSS.totalForce[1] += Wstar[4]*gradPhi[1]*area;
      LSS.totalForce[2] += Wstar[4]*gradPhi[2]*area;
*/      LSS.totalForce[0] += Wstar[4]*normal[l][0];
      LSS.totalForce[1] += Wstar[4]*normal[l][1];
      LSS.totalForce[2] += Wstar[4]*normal[l][2];


      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi);
      for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];

      if (nodeTag[j]<0.0)
        riemann.computeFSIRiemannSolution((int)nodeTag[j],Vj,vStar,gradPhi,varFcn,Wstar,i);
      else
        riemann.computeFSIRiemannSolution((int)nodeTag[j],Vj,vStar,-1.0*gradPhi,varFcn,Wstar,i);
      fprintf(stderr,"Vj = (%e %e %e %e %e), vStar = (%e %e %e %e %e). nodeTag[j]=%e.\n", Vj[0], Vj[1],Vj[2], Vj[3], Vj[4], Wstar[0], Wstar[1], Wstar[2], Wstar[3], Wstar[4], nodeTag[j]);

      for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
/*      LSS.totalForce[0] += -Wstar[4]*gradPhi[0]*area;
      LSS.totalForce[1] += -Wstar[4]*gradPhi[1]*area;
      LSS.totalForce[2] += -Wstar[4]*gradPhi[2]*area;
*/      LSS.totalForce[0] += -Wstar[4]*normal[l][0];
      LSS.totalForce[1] += -Wstar[4]*normal[l][1];
      LSS.totalForce[2] += -Wstar[4]*normal[l][2];


      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj);
      for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
    }
  }
  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                      ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, NodalGrad<dim>& ngrad,
                                      NodalGrad<1> &ngradLS,
                                      EdgeGrad<dim>* egrad, SVec<double,1>& Phi, Vec<double>& PhiF)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx  = ngrad.getX();
  SVec<double,dim>& dVdy  = ngrad.getY();
  SVec<double,dim>& dVdz  = ngrad.getZ();
  // in this routine Phi denotes the "conservative phi" ie (rho*phi)
  SVec<double,1>& dPhidx = ngradLS.getX();
  SVec<double,1>& dPhidy = ngradLS.getY();
  SVec<double,1>& dPhidz = ngradLS.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim];
  double ddPij[1], ddPji[1], Pi[2], Pj[2];
  double Uni, Unj;
  double Phia;
  double srho1, srho2, srhod;
  double unroe, rroe;
  double uroe;

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;
    int i = ptr[l][0];
    int j = ptr[l][1];
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }
    ddPij[0] = dx[0]*dPhidx[i][0] + dx[1]*dPhidy[i][0] + dx[2]*dPhidz[i][0];
    ddPji[0] = dx[0]*dPhidx[j][0] + dx[1]*dPhidy[j][0] + dx[2]*dPhidz[j][0];

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    recFcnLS->compute(Phi[i], ddPij, Phi[j], ddPji, Pi, Pj);

    Uni      = Vi[1]*normal[l][0]  +Vi[2]*normal[l][1]  +Vi[3]*normal[l][2] - normalVel[l];
    Unj      = Vj[1]*normal[l][0]  +Vj[2]*normal[l][1]  +Vj[3]*normal[l][2] - normalVel[l];
    //Roe averaged variables
    if(fabs(Pj[0]-Pi[0])<1.e-12*fabs(Pi[0]) || Pj[0]==Pi[0])
      uroe     = 0.5*(Unj + Uni);
    else
      uroe     = (Pj[0]*Unj - Pi[0]*Uni)/(Pj[0]-Pi[0]);

    // roe flux
    Phia = Uni*Pi[0] + Unj*Pj[0] - fabs(uroe)*(Pj[0]-Pi[0]);

    PhiF[i] += 0.5*Phia;
    PhiF[j] -= 0.5*Phia;
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
template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                              FluxFcn **fluxFcn, GeoState &geoState,
                                              NodalGrad<dim> &ngrad, NodalGrad<1> &ngradLS,
                                              SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<double> &Phi)
{
	// it is assumed that dim=5, ie no turbulence possible
  int k;

  double gradphi[3];
  double gphii[3];
  double gphij[3];
  double length;

  double dfdUi[neq*neq], dfdUj[neq*neq];
  double dfdUk[neq*neq], dfdUl[neq*neq];
  double Vi[2*dim], Vj[2*dim];
  double Wi[2*dim], Wj[2*dim];

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,1>& dPdx = ngradLS.getX();
  SVec<double,1>& dPdy = ngradLS.getY();
  SVec<double,1>& dPdz = ngradLS.getZ();

  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

  riemann.reset(0);

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (Phi[i]*Phi[j] > 0.0) {
       if (Phi[i] > 0.0  && Phi[j] > 0.0){
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       }
       if (Phi[i] < 0.0  && Phi[j] < 0.0){
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
       }
    } else {
      //ngradLS returns nodal gradients of phi
      gphii[0] = dPdx[i][0];
      gphii[1] = dPdy[i][0];
      gphii[2] = dPdz[i][0];
      gphij[0] = dPdx[j][0];
      gphij[1] = dPdy[j][0];
      gphij[2] = dPdz[j][0];
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

      int epsi = 0;
      int epsj = 0;
      riemann.computeRiemannSolution(Vi,Vj,Phi[i],Phi[j],gradphi,varFcn,
                                     epsi,epsj,Wi,Wj,i,j);
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Vi, Wi, dfdUi, dfdUk, epsi);
      fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], Wj, Vj, dfdUl, dfdUj, epsj);

    }

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
    if (Aij && Aji) {
      double voli = 1.0 / ctrlVol[i];
      double volj = 1.0 / ctrlVol[j];
      for (k=0; k<neq*neq; ++k) {
        Aij[k] += dfdUj[k] * voli;
        Aji[k] -= dfdUi[k] * volj;
      }
    }

  }
}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                              FluxFcn **fluxFcn, GeoState &geoState,
                                              NodalGrad<dim> &ngrad, NodalGrad<1> &ngradLS,
                                              SVec<double,3> &X,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<double> &Phi,
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

  double dfdUi[neq*neq], dfdUj[neq*neq];
  bool atInleti, atInletj;
  double length;

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,1>& dPdx = ngradLS.getX();
  SVec<double,1>& dPdy = ngradLS.getY();
  SVec<double,1>& dPdz = ngradLS.getZ();

  VarFcn *varFcn;

  double alpha, beta, pref, gam;
  double R_g,U_g,P_g;
  double R_w,U_w,P_w, T_w;
  double P_i,U_i,R_il,R_ir;
  double Vwater[dim];
  double Vgas[dim];
  double dfdUk[neq*neq], dfdUl[neq*neq];

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    double gradphi[3];
    double gphii[3];
    double gphij[3];
    //ngradLS returns nodal gradients of phi
    gphii[0] = dPdx[i][0];
    gphii[1] = dPdy[i][0];
    gphii[2] = dPdz[i][0];
    gphij[0] = dPdx[j][0];
    gphij[1] = dPdy[j][0];
    gphij[2] = dPdz[j][0];
    for (k=0; k<3; k++)
      gradphi[k] = 0.5*(gphii[k]+gphij[k]);
    double normgradphi = sqrt(gradphi[0]*gradphi[0]+gradphi[1]*gradphi[1]+gradphi[2]*gradphi[2]);
    for (k=0; k<3; k++)
      gradphi[k] /= normgradphi;


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

    if (Phi[i]*Phi[j] >   0.0) {
       if (Phi[i] >  0.0  && Phi[j] >  0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       if (Phi[i] <  0.0  && Phi[j] <  0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
    } else {
       varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
       if (varFcn->getType() == VarFcn::GASINLIQUID) {
				 fprintf(stdout, "Big Error here\n");
         //RiemannJacobianGasTait(i,j,V,Phi[i],Phi[j],gradphi,normal[l],normalVel[l],varFcn,fluxFcn,dfdUi,dfdUj);
       }
       else{
        if (Phi[i] >= 0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
        else
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);

        if (Phi[j] >= 0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
        else
         fluxFcn[BC_INTERNAL]->computeJacobians(length, 0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);

       }
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

//------------------------------------------------------------------------------
