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
#include <ExactRiemannSolver.h>
#include <Tet.h>
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
		assert(S>0.0);
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
                              double beta, double k1, double cmach, Vec<double> &Phi, int subnum)
{
  double Phimid, Vmid[dim];
  double S, invS, ndot;
  double a, un, mach;
  double beta2, coeff1, coeff2;
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

      beta2 = beta * beta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;
      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
    }else{
      u = varFcn->getVelocity(V[i]);
      a = varFcn->computeSoundSpeed(V[i], Phi[i]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[i],Phi[i]);

      beta2 = beta * beta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

      dt[i] += min(0.5*(coeff1-coeff2), 0.0) * S;

      u = varFcn->getVelocity(V[j]);
      a = varFcn->computeSoundSpeed(V[j], Phi[j]);
      un = u * n - ndot;
      mach = varFcn->computeMachNumber(V[j],Phi[j]);

      beta2 = beta * beta;
      coeff1 = (1.0+beta2)*un;
      coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

      dt[j] += min(0.5*(-coeff1-coeff2), 0.0) * S;
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
int EdgeSet::computeFiniteVolumeTerm(int* locToGlobNodeMap, Vec<double> &irey, FluxFcn** fluxFcn,
                                     RecFcn* recFcn, TetSet& tets, GeoState& geoState, SVec<double,3>& X,
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
  double edgeirey;

  int ierr = 0;

  for (int l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

    if (egrad)
      egrad->compute(l, i, j, tets, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
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
int EdgeSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     TetSet& tets, GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, Vec<double> &Phi,
                                     NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                     NodalGrad<1>& ngradLS,
                                     SVec<double,dim>& fluxes, int it,
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

  int ierr=0;
  riemann.reset(it);

  for (int l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;

    int i = ptr[l][0];
    int j = ptr[l][1];

  //fprintf(stdout, "Edge.C5\n");
    //if (egrad)
    //  egrad->compute(l, i, j, tets, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    //else {
    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }
    //}

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    if (!rshift)
    // check for negative pressure or density //
      ierr += checkReconstructedValues(i, j, Vi, Vj, varFcn, locToGlobNodeMap,
                                       failsafe, tag, Phi[i], Phi[j]);

    if (ierr) continue;

    for (int k=0; k<dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (Phi[i]*Phi[j] > 0.0) { 	// same fluid
      if (Phi[i] > 0.0 && Phi[j] > 0.0)  {
        fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, 1);
      }
      if (Phi[i] < 0.0 && Phi[j] < 0.0)  {
        fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l], Vi, Vj, flux, -1);
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
      fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l],
                                    Vi, Wi, fluxi, epsi);
      fluxFcn[BC_INTERNAL]->compute(0.0, normal[l], normalVel[l],
                                    Wj, Vj, fluxj, epsj);
      for (int k=0; k<dim; k++){
        fluxes[i][k] += fluxi[k];
        fluxes[j][k] -= fluxj[k];
      }
    }
  }
  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void EdgeSet::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                      TetSet& tets, GeoState& geoState, SVec<double,3>& X,
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
    //if (egrad)
    //  egrad->compute(l, i, j, tets, X, V, dVdx, dVdy, dVdz, ddVij, ddVji);
    //else {
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
        ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      }
    //}

    //if (egradLS)
    //  egradLS->compute(l, i, j, tets, X, Phi, dPhidx, dPhidy, dPhidz, ddPij, ddPji);
    //else {
      //double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
        ddPij[0] = dx[0]*dPhidx[i][0] + dx[1]*dPhidy[i][0] + dx[2]*dPhidz[i][0];
        ddPji[0] = dx[0]*dPhidx[j][0] + dx[1]*dPhidy[j][0] + dx[2]*dPhidz[j][0];
      //}
    //}

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
    recFcnLS->compute(Phi[i], ddPij, Phi[j], ddPji, Pi, Pj);

    Uni      = Vi[1]*normal[l][0]  +Vi[2]*normal[l][1]  +Vi[3]*normal[l][2];
    Unj      = Vj[1]*normal[l][0]  +Vj[2]*normal[l][1]  +Vj[3]*normal[l][2];
    //Roe averaged variables
    if(fabs(Pj[0]-Pi[0])<1.e-12*fabs(Pi[0]) || Pj[0]==Pi[0])
      uroe     = 0.5*(Unj + Uni);
    else
      uroe     = (Pj[0]*Unj - Pi[0]*Uni)/(Pj[0]-Pi[0]);

    //dynamic mesh inclusion
    uroe -= normalVel[l];

    // roe flux
    Phia = Uni*Pi[0] + Unj*Pj[0] - fabs(uroe)*(Pj[0]-Pi[0]);

    PhiF[i] += 0.5*Phia;
    PhiF[j] -= 0.5*Phia;
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
void EdgeSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, GeoState &geoState,
					      NodalGrad<dim> &ngrad, NodalGrad<1> &ngradLS,
                                              Vec<double> &ctrlVol, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<double> &Phi)
{
  int k;
  double dfdUi[neq*neq], dfdUj[neq*neq];
  Vec<Vec3D> &normal = geoState.getEdgeNormal();
  Vec<double> &normalVel = geoState.getEdgeNormalVel();
  SVec<double,dim>& dVdx = ngrad.getX();
  SVec<double,dim>& dVdy = ngrad.getY();
  SVec<double,dim>& dVdz = ngrad.getZ();
  SVec<double,1>& dPdx = ngradLS.getX();
  SVec<double,1>& dPdy = ngradLS.getY();
  SVec<double,1>& dPdz = ngradLS.getZ();


  for (int l=0; l<numEdges; ++l) {
    int i = ptr[l][0];
    int j = ptr[l][1];


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


    if (Phi[i]*Phi[j] > 0.0) { 
       if (Phi[i] > 0.0  && Phi[j] > 0.0){
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       }
       if (Phi[i] < 0.0  && Phi[j] < 0.0){
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
       }
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
				 fprintf(stdout, "Big Error here\n");
        // RiemannJacobianGasTait(i,j,V,Phi[i],Phi[j],gradphi,normal[l],normalVel[l],varFcn,fluxFcn,dfdUi,dfdUj);
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
					      NodalGrad<dim> &ngrad, NodalGrad<1> &ngradLS,
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
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, 1);
       if (Phi[i] <  0.0  && Phi[j] <  0.0)
         fluxFcn[BC_INTERNAL]->computeJacobians(0.0, normal[l], normalVel[l], V[i], V[j], dfdUi, dfdUj, -1);
    } else {
       varFcn  = fluxFcn[BC_INTERNAL]->getVarFcn();
       if (varFcn->getType() == VarFcn::GASINLIQUID) {
				 fprintf(stdout, "Big Error here\n");
         //RiemannJacobianGasTait(i,j,V,Phi[i],Phi[j],gradphi,normal[l],normalVel[l],varFcn,fluxFcn,dfdUi,dfdUj);
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
