#include <math.h>

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
                              TimeLowMachPrec &tprec, Vec<int> &fluidId, int subnum)
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

      locbeta = tprec.getBeta(mach);
      beta2 = locbeta*locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
    }else{
      u = varFcn->getVelocity(V[i]);
      a = varFcn->computeSoundSpeed(V[i], fluidId[i]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[i],fluidId[i]);

      locbeta = tprec.getBeta(mach);
      beta2 = locbeta * locbeta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;

      u = varFcn->getVelocity(V[j]);
      a = varFcn->computeSoundSpeed(V[j], fluidId[j]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[j],fluidId[j]);

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

template<int dim, int dimLS>
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, Vec<int> &fluidId,
                                     FluidSelector &fluidSelector, 
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     NodalGrad<dimLS>& ngradLS,
                                     SVec<double,dim>& fluxes, int it,
                                     SVec<double,dim>* interfaceFlux,
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
                                       failsafe, tag, V[i], V[j], fluidId[i], fluidId[j]);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (fluidId[i]==fluidId[j]) { 	// same fluid
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
                                    Wi,Wj,i,j,l,dx);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Vi, Wi, fluxi, fluidId[i]);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Wj, Vj, fluxj, fluidId[j]);
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
                                     SVec<double,dim>* interfaceFlux,
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
  double gradphi[3], gphii[3], gphij[3];

  // ------------------------------------------------
  //  THE FAMOUS, MOUTH-WATERING EDGE LOOP...
  // ------------------------------------------------
  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,i,j);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);

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
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, i, j);
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

        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi, fluidId[i], false);
          for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, j,i);
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
                                    Wi,Wj,i,j,l,dx);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Vi, Wi, fluxi, fluidId[i]);
      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l],
                                    Wj, Vj, fluxj, fluidId[j]);
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
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, SVec<double,dim>& Wstarij, 
                                     SVec<double,dim>& Wstarji, LevelSetStructure &LSS,
                                     bool linRecAtInterface, int Nriemann, SVec<double,3>* Nsbar,
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
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  Vec3D normalDir; //normal direction fed to the Riemann solver. Related to the choice of "Nriemann".

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool iIsActive = LSS.isActive(0, i);
    bool jIsActive = LSS.isActive(0, j);
    bool intersect = LSS.edgeIntersectsStructure(0,i,j);

    if( !iIsActive && !jIsActive ) {
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
          if (Wstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) Vi[k] = V[i][k]; 
          } else recFcn->compute(V[i], ddVij, Wstarij[l], ddVji, Vi, Vtemp);
          if (Wstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) Vj[k] = V[j][k];
          } else recFcn->compute(Wstarji[l], ddVij, V[j], ddVji, Vtemp, Vj);
        }
        // Case 2: i is active, j is inactive
        else if (iIsActive)
          if (Wstarij[l][0]<1e-8) {// no riemann sol. (first time-step)
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(V[i], ddVij, Wstarij[l], ddVji, Vi, Vj);
        // Case 3: i is inactive, j is active
        else if (jIsActive)
          if (Wstarji[l][0]<1e-8) {
            for (int k=0; k<dim; k++) {Vi[k] = V[i][k]; Vj[k] = V[j][k];}
          } else recFcn->compute(Wstarji[l], ddVij, V[j], ddVji, Vi, Vj);
      }
    }

    if(it>0)
      for(int k=0;k<dim;k++)
        Wstarij[l][k] = Wstarji[l][k] = 0.0;

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
        switch (Nriemann) {
          case 0: //structure normal
            normalDir = res.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = -1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            normalDir = Vec3D((*Nsbar)[i][0], (*Nsbar)[i][1], (*Nsbar)[i][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());

        riemann.computeFSIRiemannSolution(Vi,res.normVel,normalDir,varFcn,Wstar,j);

        if (it>0) //if it>0 (i.e. not called in computeResidualNorm), store Wstarij.
          for (int k=0; k<dim; k++)  Wstarij[l][k] = Wstar[k]; 
        if (masterFlag[l]) { 
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Wstar, fluxi, 0, false);
          for (int k=0; k<dim; k++) fluxes[i][k] += fluxi[k];
        }
      }
      if (jIsActive) {
        LevelSetResult res = LSS.getLevelSetDataAtEdgeCenter(0.0, j,i);
        switch (Nriemann) {
          case 0: //structure normal
            normalDir = res.gradPhi;
            break;
          case 1: //fluid normal
            normalDir = 1.0/(normal[l].norm())*normal[l];
            break;
          case 2: //cell-averaged structure normal
            normalDir = Vec3D((*Nsbar)[j][0], (*Nsbar)[j][1], (*Nsbar)[j][2]);
            break;
          default:
            fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
            exit(-1);
        }
        if(std::abs(1.0-normalDir.norm())>0.1)
          fprintf(stderr,"KW: normalDir.norm = %e. This is too bad...\n", normalDir.norm());

        riemann.computeFSIRiemannSolution(Vj,res.normVel,normalDir,varFcn,Wstar,i);
        
        if (it>0)
          for (int k=0; k<dim; k++) Wstarji[l][k] = Wstar[k];
        if (masterFlag[l]) {
          fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Wstar, Vj, fluxj, 0, false);
          for (int k=0; k<dim; k++)  fluxes[j][k] -= fluxj[k];
        }
      }
    }
  }
  return ierr;
}

//------------------------------------------------------------------------------
// for FSF (fluid-shell-fluid).
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

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Wstar[2*dim];
  double fluxi[dim], fluxj[dim];  for (int i=0; i<dim; i++) fluxi[i] = fluxj[i] = 0.0;
  VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();
  double length;

  Vec3D normalDir;

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    bool intersect = LSS.edgeIntersectsStructure(0,i,j);
    bool iActive = LSS.isActive(0.0,i);
    bool jActive = LSS.isActive(0.0,j);
    
    if(!iActive && !jActive) continue;

    // ------------------------------------------------
    //  Reconstruction without crossing the interface.
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


    // --------------------------------------------------------
    //                   Compute fluxes
    // --------------------------------------------------------
    if (!intersect) {  // same fluid
      if (!masterFlag[l]) continue; //not a master edge

      fluxFcn[BC_INTERNAL]->compute(length, 0.0, normal[l], normalVel[l], Vi, Vj, flux, fluidId[i]);
      for (int k=0; k<dim; ++k) {
        fluxes[i][k] += flux[k];
        fluxes[j][k] -= flux[k];
      }
    }
    else{// interface

      // for node i
      if(iActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, i, j);

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
        }
      }

      // for node j
      if(jActive){
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, j,i);

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
        }
      }
    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void EdgeSet::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                      ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, NodalGrad<dim>& ngrad,
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
  double Uni, Unj;
  double Phia;
  double srho1, srho2, srhod;
  double unroe, rroe;
  double uroe;
  int k;

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;
    int i = ptr[l][0];
    int j = ptr[l][1];
    int iCovered = LSS ? LSS->fluidModel(0.0,i) : 0;
    int jCovered = LSS ? LSS->fluidModel(0.0,j) : 0;
      //when i(j)Covered = 0, we have valid Phi (and U) on this node.

    if(iCovered && jCovered) continue;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }
    for (k=0; k<dimLS; ++k) {
      ddPij[k] = dx[k]*dPhidx[i][k] + dx[1]*dPhidy[i][k] + dx[2]*dPhidz[i][k];
      ddPji[k] = dx[k]*dPhidx[j][k] + dx[1]*dPhidy[j][k] + dx[2]*dPhidz[j][k];
    }

    if(!iCovered && !jCovered) { //the usual linear reconstruction
      recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
      recFcnLS->compute(Phi[i], ddPij, Phi[j], ddPji, Pi, Pj);
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

    // roe flux
    if(!iCovered && !jCovered)
      for (k=0; k<dimLS; k++){
        // Original
        if(fabs(Pj[k]-Pi[k])<1.e-12*fabs(Pi[k]) || Pj[k]==Pi[k])
          uroe     = 0.5*(Unj + Uni);
        else
          uroe     = (Pj[k]*Unj - Pi[k]*Uni)/(Pj[k]-Pi[k]);
        Phia = Uni*Pi[k] + Unj*Pj[k] - fabs(uroe)*(Pj[k]-Pi[k]);
        PhiF[i][k] += 0.5*Phia;
        PhiF[j][k] -= 0.5*Phia;
  
        //Alex    
  /*      if (Uni > 0 && Unj > 0)
  	Phia = (Pi[k]*Uni);
        else if (Uni < 0 && Unj < 0)
          Phia = (Pj[k]*Unj);
        else
          Phia = 0.0;
      
        PhiF[i][k] += Phia;
        PhiF[j][k] -= Phia;*/
      }
    else {
      if(!iCovered) {
        for(k=0; k<dimLS; k++)
          PhiF[i][k] += 0.5*Uni*Pi[k];
      }
      if(!jCovered) {
        for(k=0; k<dimLS; k++)
          PhiF[j][k] -= 0.5*Unj*Pj[k];
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
                                     Wi,Wj,i,j,l,dx);
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
                                     Wi,Wj,i,j,l,dx);

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

extern double implbeta_alex;
template<class Scalar, int dim, int dimLS>
void EdgeSet::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
						GeoState& geoState, SVec<double,3>& X,
						SVec<double,dim>& V, NodalGrad<dim>& ngrad,
						NodalGrad<dimLS> &ngradLS,
						EdgeGrad<dim>* egrad,
						Vec<double> &ctrlVol, SVec<double,dimLS>& Phi,
						GenMat<Scalar,dimLS> &A)
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
    for (k=0; k<dimLS; ++k) {
      ddPij[k] = dx[k]*dPhidx[i][k] + dx[1]*dPhidy[i][k] + dx[2]*dPhidz[i][k];
      ddPji[k] = dx[k]*dPhidx[j][k] + dx[1]*dPhidy[j][k] + dx[2]*dPhidz[j][k];
    }
    
    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    //recFcnLS->compute(Phi[i], ddPij, Phi[j], ddPji, Pi, Pj);
    
    Uni      = Vi[1]*normal[l][0]  +Vi[2]*normal[l][1]  +Vi[3]*normal[l][2] - normalVel[l];
    Unj      = Vj[1]*normal[l][0]  +Vj[2]*normal[l][1]  +Vj[3]*normal[l][2] - normalVel[l];
    //Roe averaged variables
    for (k = 0; k < dimLS; ++k) {

      df[0] = df[1] = 0.0;

      if (Uni > 0 && Unj > 0)
	df[0] = Uni;
      else if (Uni < 0 && Unj < 0)
	df[1] = Unj;
      
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
