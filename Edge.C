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
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                              SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, Vec<double> &idtv,
                              double beta, double k1, double cmach)
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
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    locbeta = fmin(fmax(k1*locMach, beta),cmach);
    
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

template<int dim>
void EdgeSet::computeTimeStep(VarFcn *varFcn, GeoState &geoState,
                              SVec<double,dim> &V, Vec<double> &dt,
                              double beta, double k1, double cmach, Vec<double> &Phi)
{
  double phia;
  double Vmid[dim];
  double locmach=0.0;
                                                                                                        
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
                                                                                                        
  for (int l=0; l<numEdges; ++l) {
                                                                                                        
    if (!masterFlag[l]) continue;
                                                                                                        
    int i = ptr[l][0];
    int j = ptr[l][1];
    
    phia  = 0.5*(Phi[i]  +Phi[j]);                                                                                                    
    double S = sqrt(normal[l] * normal[l]);
    double invS = 1.0 / S;
                                                                                                        
    Vec3D n = invS * normal[l];
    double ndot = invS * normalVel[l];
                                                                                                        
    for (int k=0; k<dim; ++k)
      Vmid[k] = 0.5 * (V[i][k] + V[j][k]);
                                                                                                        
    Vec3D u = varFcn->getVelocity(Vmid);
    double a = varFcn->computeSoundSpeed(Vmid, phia);
    double un = u * n - ndot;
    double mach = varFcn->computeMachNumber(Vmid, phia);
                                                                                                        
    double beta2 = beta * beta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);
                                                                                                        
    dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
    dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
                                                                                                        
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
  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double edgeirey;

  int ierr = 0;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    if (egrad)
      egrad->compute(l, i, j, elems, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    edgeirey = 0.5*(irey[i]+irey[j]);

   if (!rshift) {
     double rho[2], p[2];
     rho[0] = Vi[0];
     p[0]   = Vi[4];
     rho[1] = Vj[0];
     p[1]   = Vj[4];

     ierr += checkReconstruction(rho, p, i, j, locToGlobNodeMap, failsafe, tag);
   }

    if (ierr) continue;

    int k;
    for (k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    
    fluxFcn[BC_INTERNAL]->compute(edgeirey, normal[l], normalVel[l], Vi, Vj, flux);

    for (k=0; k<dim; ++k) {
      fluxes[i][k] += flux[k];
      fluxes[j][k] -= flux[k];
    }

  }

  return(ierr);

}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(int* locToGlobNodeMap, FluxFcn** fluxFcn, RecFcn* recFcn,
                                     ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, Vec<double> &Phi,
                                     NodalGrad<dim>& ngrad,
                                     EdgeGrad<dim>* egrad, SVec<double,dim>& fluxes,
                                     SVec<int,2>& tag, int failsafe, int rshift)
{

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();

  int iside;
  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], flux[dim];
  double Vwater[2*dim], Vgas[2*dim];
  double T_w, P_g, E_w, P_w, T_g, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir, c2_il, c2_ir, alpha, beta, pref, gam;
  double R_r, R_l;
  double uf_g, vf_g, wf_g;
  double uf_w, vf_w, wf_w;
  double V_g, W_g, V_w, W_w;
  double G_v, W_v;
  double unorm,vnorm,wnorm;
  double velnorm;
  double unorm2,vnorm2,wnorm2;
  double velnorm2;
  double eps;
  VarFcn *varFcn;

  eps = 1.e-6;
  int ierr=0;

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];

    if (egrad)
      egrad->compute(l, i, j, elems, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    if (!rshift) {
      double rho[2], p[2];
     rho[0] = Vi[0];
     p[0]   = Vi[4];
     rho[1] = Vj[0];
     p[1]   = Vj[4];

     ierr += checkReconstruction(rho, p, i, j, locToGlobNodeMap, failsafe, tag);
   }
   if (ierr) continue;

    int k;
    for (k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (Phi[i]*Phi[j] >   0.0) {
      if (Phi[i] >  0.0  && Phi[j] >  0.0)  {
        fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);
      }
      if (Phi[i] <  0.0  && Phi[j] <  0.0)  {
        fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);
      }
      for (k=0; k<dim; ++k) {
        if (masterFlag[l]) fluxes[i][k] += flux[k];
        if (masterFlag[l]) fluxes[j][k] -= flux[k];
      }
    }
    else{
      varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
      if (varFcn->getType() == VarFcn::GASINLIQUID) {
        alpha   = varFcn->getAlphaWater();
        beta    = varFcn->getBetaWater();
        pref    = varFcn->getPrefWater();
        gam     = varFcn->getGamma();
        if (Phi[i] >= 0.0)  {
          R_g  = Vj[0];  R_w  = Vi[0];
          velnorm  = max(sqrt(Vj[1]*Vj[1]  +Vj[2]*Vj[2]  +Vj[3]*Vj[3]), eps);
          unorm    = Vj[1]/velnorm; vnorm  = Vj[2]/velnorm; wnorm  = Vj[3]/velnorm;
          velnorm2 = max(sqrt(Vi[1]*Vi[1]  +Vi[2]*Vi[2]  +Vi[3]*Vi[3]), eps);
          unorm2   = Vi[1]/velnorm2; vnorm2  = Vi[2]/velnorm2; wnorm2  = Vi[3]/velnorm2;
          //U_g  = Vj[1];  U_w  = Vi[1];
          U_g  = velnorm;  U_w  = velnorm2;
          P_g  = varFcn->getPressure(Vj, Phi[j]);
          P_w  = varFcn->getPressure(Vi, Phi[i]);
          fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
          for (k=0; k<2*dim-1; k++){
            Vwater[k] = Vj[k];
          }
          Vwater[0]  = R_ir;   Vwater[dim]    = R_ir;
          //Vwater[1]  = U_i;   Vwater[dim+1]  = U_i;
          Vwater[1]  = U_i*unorm;   Vwater[dim+1]  = U_i*unorm;
          Vwater[2]  = U_i*vnorm;   Vwater[dim+2]  = U_i*vnorm;
          Vwater[3]  = U_i*wnorm;   Vwater[dim+3]  = U_i*wnorm;
          Vwater[4]  = P_i;   Vwater[2*dim-1] = P_i;
          T_w  = varFcn->computeTemperature(Vwater, Phi[j]);
          Vwater[4]  = T_w;    Vwater[2*dim-1]= T_w;
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vwater, flux, 1);
        }
        else{
          R_g  = Vi[0];  R_w  = Vj[0];
          velnorm  = max(sqrt(Vj[1]*Vj[1]  +Vj[2]*Vj[2]  +Vj[3]*Vj[3]), eps);
          unorm    = Vj[1]/velnorm; vnorm  = Vj[2]/velnorm; wnorm  = Vj[3]/velnorm;
          velnorm2 = max(sqrt(Vi[1]*Vi[1]  +Vi[2]*Vi[2]  +Vi[3]*Vi[3]), eps);
          unorm2   = Vi[1]/velnorm2; vnorm2  = Vi[2]/velnorm2; wnorm2  = Vi[3]/velnorm2;
          //U_g  = Vi[1];  U_w  = Vj[1];
          U_g  = velnorm2;  U_w  = velnorm;
          P_g  = varFcn->getPressure(Vi, Phi[i]);
          P_w  = varFcn->getPressure(Vj, Phi[j]);
          fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
          for (k=0; k<2*dim-1; k++)
            Vgas[k] = Vj[k];
          Vgas[0]  = R_il;   Vgas[dim]    = R_il;
          //Vgas[1]  = U_i;    Vgas[dim+1]  = U_i;
          Vgas[1]  = U_i*unorm;    Vgas[dim+1]  = U_i*unorm;
          Vgas[2]  = U_i*vnorm;    Vgas[dim+2]  = U_i*vnorm;
          Vgas[3]  = U_i*wnorm;    Vgas[dim+3]  = U_i*wnorm;
          Vgas[4]  = P_i;    Vgas[2*dim-1]= P_i;
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vgas, flux, -1);
        }
        for (k=0; k<dim; ++k)
          if (masterFlag[l]) fluxes[i][k] += flux[k];

        if (Phi[j] >= 0.0)  {
          R_g  = Vi[0];  R_w  = Vj[0];
          velnorm  = max(sqrt(Vj[1]*Vj[1]  +Vj[2]*Vj[2]  +Vj[3]*Vj[3]), eps);
          unorm    = Vj[1]/velnorm; vnorm  = Vj[2]/velnorm; wnorm  = Vj[3]/velnorm;
          velnorm2 = max(sqrt(Vi[1]*Vi[1]  +Vi[2]*Vi[2]  +Vi[3]*Vi[3]), eps);
          unorm2   = Vi[1]/velnorm2; vnorm2  = Vi[2]/velnorm2; wnorm2  = Vi[3]/velnorm2;
          //U_g  = Vi[1];  U_w  = Vj[1];
          U_g  = velnorm2;  U_w  = velnorm;
          P_g  = varFcn->getPressure(Vi, Phi[i]);
          P_w  = varFcn->getPressure(Vj, Phi[j]);
          fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
          for (k=0; k<2*dim-1; k++)
            Vwater[k] = Vi[k];
          Vwater[0]  = R_ir;   Vwater[dim]    = R_ir;
          //Vwater[1]  = U_i;   Vwater[dim+1]  = U_i;
          Vwater[1]  = U_i*unorm2;   Vwater[dim+1]  = U_i*unorm2;
          Vwater[2]  = U_i*vnorm2;   Vwater[dim+2]  = U_i*vnorm2;
          Vwater[3]  = U_i*wnorm2;   Vwater[dim+3]  = U_i*wnorm2;
          Vwater[4]  = P_i;   Vwater[2*dim-1] = P_i;
          T_w  = varFcn->computeTemperature(Vwater, Phi[i]);
          Vwater[4]  = T_w;    Vwater[2*dim-1]= T_w;
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vwater, Vj, flux, 1);
        }
        else   {
          R_g  = Vj[0];  R_w  = Vi[0];
          velnorm  = max(sqrt(Vj[1]*Vj[1]  +Vj[2]*Vj[2]  +Vj[3]*Vj[3]), eps);
          unorm    = Vj[1]/velnorm; vnorm  = Vj[2]/velnorm; wnorm  = Vj[3]/velnorm;
          velnorm2 = max(sqrt(Vi[1]*Vi[1]  +Vi[2]*Vi[2]  +Vi[3]*Vi[3]), eps);
          unorm2   = Vi[1]/velnorm2; vnorm2  = Vi[2]/velnorm2; wnorm2  = Vi[3]/velnorm2;
          U_g  = velnorm;  U_w  = velnorm2;
          //U_g  = Vj[1];  U_w  = Vi[1];
          P_g  = varFcn->getPressure(Vj, Phi[j]);
          P_w  = varFcn->getPressure(Vi, Phi[i]);
          fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
          for (k=0; k<2*dim-1; k++)
            Vgas[k] = Vi[k];
          Vgas[0]  = R_il;   Vgas[dim]    = R_il;
          //Vgas[1]  = U_i;    Vgas[dim+1]  = U_i;
          Vgas[1]  = U_i*unorm2;    Vgas[dim+1]  = U_i*unorm2;
          Vgas[2]  = U_i*vnorm2;    Vgas[dim+2]  = U_i*vnorm2;
          Vgas[3]  = U_i*wnorm2;    Vgas[dim+3]  = U_i*wnorm2;
          Vgas[4]  = P_i;    Vgas[2*dim-1]= P_i;
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vgas, Vj, flux, -1);
        }
        for (k=0; k<dim; ++k)
          if (masterFlag[l]) fluxes[j][k] -= flux[k];
      }
      else {
        if (Phi[i] >= 0.0)  {
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);}
        else{
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);}

        for (k=0; k<dim; ++k)
          if (masterFlag[l]) fluxes[i][k] += flux[k];
        if (Phi[j] >= 0.0)  {
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);}
        else   {
          fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);}

        for (k=0; k<dim; ++k)
          if (masterFlag[l]) fluxes[j][k] -= flux[k];
      }
    }
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::storeGhost(SVec<double,dim> &V, SVec<double,dim> &Vgf, Vec<double> &Phi)
{
  for (int i=0; i<Phi.size(); i++){
    if (Phi[i] < 0.0) {
       Vgf[i][1]  = V[i][1];
       Vgf[i][2]  = V[i][2];
       Vgf[i][3]  = V[i][3];
       for (int j=0; j<Phi.size(); j++){
         if (Phi[j] >=0.0  && Phi[j] < 0.01) {
             //for (int k=0; k<dim; k++) {
             //  Vgf[i][k]  = min(Vgf[i][k], V[j][k]);
             Vgf[i][0]  = min(Vgf[i][0], V[j][0]);
             Vgf[i][dim-1]  = min(Vgf[i][dim-1], V[j][dim-1]);
             //}
         }
       }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
					ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
					SVec<double,dim>& V, NodalGrad<dim>& ngrad,
					NodalGrad<dim> &ngrad1,
					EdgeGrad<dim>* egrad, Vec<double>& Phi, Vec<double>& PhiF,
					SVec<double,dim> &PhiS)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  Vec<double>& normalVel = geoState.getEdgeNormalVel();

  SVec<double,dim>& dVdx  = ngrad.getX();
  SVec<double,dim>& dVdy  = ngrad.getY();
  SVec<double,dim>& dVdz  = ngrad.getZ();
  SVec<double,dim>& dVdx1 = ngrad1.getX();
  SVec<double,dim>& dVdy1 = ngrad1.getY();
  SVec<double,dim>& dVdz1 = ngrad1.getZ();

  double ddVij[dim], ddVji[dim], Vi[2*dim], Vj[2*dim], Vij[dim], Phia, Un, Uni, Unj;
  double ddVij1[dim], ddVji1[dim], Vi1[2*dim], Vj1[2*dim], Vij1[dim];
  double Phia1, Phia2, Phiaa, Phi1, Phi2;
  double rho1, rho2, srho1, srho2, srhod;
  double unroe, rroe;
                                                                                                  
  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;
    int i = ptr[l][0];
    int j = ptr[l][1];
    if (egrad)
      egrad->compute(l, i, j, elems, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    }

    if (egrad)
      egrad->compute(l, i, j, elems, X, PhiS, dVdx1, dVdy1, dVdz1, ddVij1, ddVji1);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddVij1[k] = dx[0]*dVdx1[i][k] + dx[1]*dVdy1[i][k] + dx[2]*dVdz1[i][k];
        ddVji1[k] = dx[0]*dVdx1[j][k] + dx[1]*dVdy1[j][k] + dx[2]*dVdz1[j][k];
      }
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    recFcnLS->compute(PhiS[i], ddVij1, PhiS[j], ddVji1, Vi1, Vj1);

    for (int k=0; k<1; ++k) {
      Uni      = Vi[1]*normal[l][0]  +Vi[2]*normal[l][1]  +Vi[3]*normal[l][2];
      Unj      = Vj[1]*normal[l][0]  +Vj[2]*normal[l][1]  +Vj[3]*normal[l][2];
      Phi1     = Vi1[0];
      Phi2     = Vj1[0];
      Un       = 0.5*(Uni  +Unj);
      Phiaa    = 0.5*(Phi1  +Phi2);
      srho1    = sqrt(Vi[0]);  srho2  = sqrt(Vj[0]); srhod = 1.0/(srho1  +srho2);
      unroe    = (srho1*Uni  +srho2*Unj)*srhod;
      rroe     = srho1*srho2;
      Phia     = Uni*Phi1  + Unj*Phi2  - fabs(rroe*unroe)*(Phi2 -Phi1);
//      Phia     = rroe*Uni*Phi1  + rroe*Unj*Phi2  - fabs(rroe*unroe)*(Phi2 -Phi1);
//      Phia     = (Uni  +fabs(Uni))*Phi1  +(Unj  -fabs(Unj))*Phi2;
//      Phia     = (Un  +fabs(Un))*Phiaa  +(Un  -fabs(Un))*Phiaa;
//    Phia     = Uni*Phi1  + Unj*Phi2  - fabs(Un)*(Phi2 -Phi1);
//    Phia     = Un*Phi1  + Un*Phi2  - fabs(Un)*(Phi2 -Phi1);
//    Phia     = Uni*Phiaa  + Unj*Phiaa - fabs(Un)*(Phi2 -Phi1);
//    Phia     = (srho1*Uni*Phi1  + srho2*Unj*Phi2)/(srho1  +srho2)  - fabs(Un)*(Phi2 -Phi1);
//      Phia     = 0.5*(Uni*Phi1  +Unj*Phi2);
      PhiF[i] += 0.5*Phia;
      PhiF[j] -= 0.5*Phia;
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState, 
					      Vec<double> &irey, 
					      Vec<double> &ctrlVol, SVec<double,dim> &V, 
					      GenMat<Scalar,neq> &A)
{

  int k;
  double edgeirey;

  double dfdUi[neq*neq], dfdUj[neq*neq];

  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();

  for (int l=0; l<numEdges; ++l) {

    int i = ptr[l][0];
    int j = ptr[l][1];
    edgeirey = 0.5*(irey[i]+irey[j]);

    if(fluxFcn){
      fluxFcn[BC_INTERNAL]->computeJacobians(edgeirey, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj);

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
					      Vec<double> &irey, 
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A,
                                              int *nodeType)
{
                                                                                                  
        /* in this function, rhs has already the values extrapolated at the inlet nodes
         * if we are in the case of water simulations
         * we are computing the jacobian matrix
         */
  if(!nodeType){
    fprintf(stderr, "exiting program since this function is only used for extrapolation\n");
    exit(1);
  }
                                                                                                  
  int k,m;
  Scalar *Aii;
  Scalar *Ajj;
  Scalar *Aij;
  Scalar *Aji;
  double edgeirey;
                                                                                                  
  double dfdUi[neq*neq], dfdUj[neq*neq];
  bool atInleti, atInletj;
                                                                                                  
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
                                                                                                  
  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
    edgeirey = 0.5*(irey[i]+irey[j]);

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
                                                                                                  
                                                                                                  
    fluxFcn[BC_INTERNAL]->computeJacobians(edgeirey, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj);
                                                                                                  
                                                                                                  
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

/* //It should be possible to modify the rhs and the jacobian at the same time
 * // see block matrices... but it does not seem to work...bug
 *       double voli = 1.0 / ctrlVol[i];
 *       // both Aij and Aji are kept to 0
 *       // but a change in the rhs is needed
 *       for (k=0; k<neq; k++){
 *         for (m=0; m<neq; m++)
 *           rhs[i][k] -= dfdUj[k*neq+m] * rhs[j][m];
 *         rhs[i][k] *= voli;
 *       }
 */
                                                                                                  
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
                                                                                                  
                                                                                                  
/*        double volj = 1.0 / ctrlVol[j];
 *        // both Aij and Aji are kept to 0
 *        // but a change in the rhs is needed
 *        for (k=0; k<neq; k++){
 *          for (m=0; m<neq; m++)
 *            rhs[j][k] += dfdUi[k*neq+m] * rhs[i][m];
 *          rhs[j][k] *= volj;
 *        }
 */
                                                                                                  
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
    else
      std::cout<<"***Error: nodes are not well defined wrt inletnodes\n";
                                                                                                  
                                                                                                  
                                                                                                  
  }
}

//----------------------------------------------------------                                                                                                      
template<int dim, class Scalar, int neq>
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<double> &Phi)
{
  int k;
  double dfdUi[neq*neq], dfdUj[neq*neq];
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();

  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];

    if (Phi[i]*Phi[j] > 0.0) { 
       if (Phi[i] > 0.0  && Phi[j] > 0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       if (Phi[i] < 0.0  && Phi[j] < 0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
    } else {
       VarFcn *varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
       if (varFcn->getType() != VarFcn::GASINLIQUID) {
         if (Phi[i] >= 0.0)
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
         else
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
                                                                                                 
         if (Phi[j] >= 0.0)
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
         else
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
       }else{
         fprintf(stderr, "*** Error: jacobian not computed for GasInWater simulations with no inlet nodes\n");
         exit(1);
       }
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
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<double> &Phi,
                                              int *nodeType)
{
                                                                                                                                                                                                     
        /* in this function, rhs has already the values extrapolated at the inlet nodes
         * if we are in the case of water simulations
         * we are computing the jacobian matrix
         */
  if(!nodeType){
    fprintf(stderr, "exiting program since this function is only used for extrapolation\n");
    exit(1);
  }
                                                                                                                                                                                                     
  int k,m;
  Scalar *Aii;
  Scalar *Ajj;
  Scalar *Aij;
  Scalar *Aji;
                                                                                                                                                                                                     
  double dfdUi[neq*neq], dfdUj[neq*neq];
  bool atInleti, atInletj;
                                                                                                                                                                                                     
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();

  double Vwater[dim], Vgas[dim];
  double dfdUk[neq*neq], dfdUl[neq*neq];
  double T_w, E_w;
  double R_g, R_w, U_g, U_w, P_g, P_w, P_i, U_i, R_il, R_ir;
  double alpha, beta, pref, gam;
  VarFcn *varFcn;                                                                                                                                                                                               
  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];
                                                                                                                                                                                                     
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
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       if (Phi[i] <  0.0  && Phi[j] <  0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
    } else {
       varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
       if (varFcn->getType() == VarFcn::GASINLIQUID) {
         alpha   = varFcn->getAlphaWater();
         beta    = varFcn->getBetaWater();
         pref    = varFcn->getPrefWater();
         gam     = varFcn->getGamma();

         if (Phi[i] >=  0.0) {
           R_g  = V[j][0];  R_w  = V[i][0];
           U_g  = V[j][1];  U_w  = V[i][1];
           P_g  = varFcn->getPressure(V[j], Phi[j]);
           P_w  = varFcn->getPressure(V[i], Phi[i]);
           fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
           for (k=0; k<dim; k++)
             Vwater[k] = V[j][k];
           Vwater[0]  = R_ir;  
           Vwater[1]  = U_i;   
           Vwater[4]  = P_i;   
           T_w  = varFcn->computeTemperature(Vwater, Phi[j]);
           Vwater[4]  = T_w;   

           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], Vwater, dfdUi, dfdUl, 1);
         }
         else {
           R_g  = V[i][0];  R_w  = V[j][0];
           U_g  = V[i][1];  U_w  = V[j][1];
           P_g  = varFcn->getPressure(V[i], Phi[i]);
           P_w  = varFcn->getPressure(V[j], Phi[j]);
           fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
           for (k=0; k<dim; k++)
             Vgas[k] = V[j][k];
           Vgas[0]  = R_il;   
           Vgas[1]  = U_i;    
           Vgas[4]  = P_i;    
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], Vgas, dfdUi, dfdUl,-1);
         }
         if (Phi[j] >= 0.0) {
           R_g  = V[i][0];  R_w  = V[j][0];
           U_g  = V[i][1];  U_w  = V[j][1];
           P_g  = varFcn->getPressure(V[i], Phi[i]);
           P_w  = varFcn->getPressure(V[j], Phi[j]);
           fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
           for (k=0; k<dim; k++)
             Vwater[k] = V[i][k];
           Vwater[0]  = R_ir;  
           Vwater[1]  = U_i;   
           Vwater[4]  = P_i;   
           T_w  = varFcn->computeTemperature(Vwater, Phi[i]);
           Vwater[4]  = T_w;    
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], Vwater, V[j], dfdUk, dfdUj, 1);
         }
         else {
           R_g  = V[j][0];  R_w  = V[i][0];
           U_g  = V[j][1];  U_w  = V[i][1];
           P_g  = varFcn->getPressure(V[j], Phi[j]);
           P_w  = varFcn->getPressure(V[i], Phi[i]);
           fluxFcn[BC_INTERNAL]->compute_riemann(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
           for (k=0; k<dim; k++)
            Vgas[k] = V[i][k];
           Vgas[0]  = R_il;  
           Vgas[1]  = U_i;  
           Vgas[4]  = P_i;  
           fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], Vgas, V[j], dfdUk, dfdUj,-1);
         }
       }
       else{
        if (Phi[i] >= 0.0)  
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
        else
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);

        if (Phi[j] >= 0.0) 
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
        else  
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);

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
                                                                                                                                                                                                     
/* //It should be possible to modify the rhs and the jacobian at the same time
 * // see block matrices... but it does not seem to work...bug
 *       double voli = 1.0 / ctrlVol[i];
 *       // both Aij and Aji are kept to 0
 *       // but a change in the rhs is needed
 *       for (k=0; k<neq; k++){
 *         for (m=0; m<neq; m++)
 *           rhs[i][k] -= dfdUj[k*neq+m] * rhs[j][m];
 *         rhs[i][k] *= voli;
 *       }
 */
                                                                                                                                                                                                     
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
                                                                                                                                                                                                     
                                                                                                                                                                                                     
/*        double volj = 1.0 / ctrlVol[j];
 *        // both Aij and Aji are kept to 0
 *        // but a change in the rhs is needed
 *        for (k=0; k<neq; k++){
 *          for (m=0; m<neq; m++)
 *            rhs[j][k] += dfdUi[k*neq+m] * rhs[i][m];
 *          rhs[j][k] *= volj;
 *        }
 */
                                                                                                                                                                                                     
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
    else
      std::cout<<"***Error: nodes are not well defined wrt inletnodes\n";
  }
}

//------------------------------------------------------------------------------
