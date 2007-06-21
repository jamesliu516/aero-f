#include <Face.h>

#include <FluxFcnDescPerfectGas.h>
#include <FluxFcnDescWaterCompressible.h>
#include <FluxFcnDescGasInGas.h>
#include <FluxFcnDescLiquidInLiquid.h>
#include <FluxFcnDescGasInLiquid.h>
#include <FemEquationTerm.h>
#include <BcData.h>
#include <BcDef.h>
#include <Tet.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>

#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

//------------------------------------------------------------------------------

template<class NodeMap>
void Face::renumberNodes(NodeMap &nodemap)
{

  for (int j=0; j<3; ++j)
    nodeNum[j] = nodemap[ nodeNum[j] ];

}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U)
{
  int k;

  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = third*(Uin[nodeNum[0]][k] + Uin[nodeNum[1]][k] + Uin[nodeNum[2]][k]);
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = third*(Uout[nodeNum[0]][k] + Uout[nodeNum[1]][k] + Uout[nodeNum[2]][k]);
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;
  

}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignFreeStreamValues(double *Uin, double *Uout, double *U)
{

  int k;

  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uin[k];
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uout[k];
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeFaceBcValue(SVec<double,dim> &Unode, double *Uface)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED)
    for (int k=0; k<dim; ++k)
      Uface[k] = third * (Unode[nodeNum[0]][k] + Unode[nodeNum[1]][k] + Unode[nodeNum[2]][k]);

}

//------------------------------------------------------------------------------

template<int dim1, int dim2>
void Face::computeNodeBcValue(SVec<double,3> &X, double *Uface, SVec<double,dim2> &Unode)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);

    for (int j=0; j<3; ++j) {
      Unode[ nodeNum[j] ][0] += S;
      for (int k=1; k<dim2; ++k) 
	Unode[ nodeNum[j] ][k] += S * Uface[dim1-dim2+k];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeForce(TetSet &tets, PostFcn *postFcn, SVec<double,3> &X, 
			Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, 
			double *nodalForceWeight, int hydro)
{

  Vec3D n;
  computeNormal(X, n);

  double dp1dxj[4][3];
  if (postFcn->doesFaceNeedGradientP1Function())
    tets[elemNum].computeGradientP1Function(X, dp1dxj);

  double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
  double *Vface[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};
  double *Vtet[4] = {V[tets[elemNum][0]], V[tets[elemNum][1]], 
		     V[tets[elemNum][2]], V[tets[elemNum][3]]};
  double *Xface[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};

  postFcn->computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, nodalForceWeight, hydro);

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeNodalForce(TetSet &tets, PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			     double pin, SVec<double,3> &F, double *nodalForceWeight)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
    || code == BC_SLIP_WALL_MOVING) {
    Vec3D Fi0, Fi1, Fi2, Fv;
    computeForce(tets, postFcn, X, d2wall, Vwall, V, &pin, Fi0, Fi1, Fi2, Fv, nodalForceWeight);
    
    Vec3D Ftot[3];
    Ftot[0] = Fi0 + third*Fv; 
    Ftot[1] = Fi1 + third*Fv;
    Ftot[2] = Fi2 + third*Fv;

    for (int j=0; j<3; ++j) {
      F[ nodeNum[j] ][0] += Ftot[j][0];
      F[ nodeNum[j] ][1] += Ftot[j][1];
      F[ nodeNum[j] ][2] += Ftot[j][2];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeNodalHeatPower(TetSet& tets, PostFcn* postFcn, SVec<double,3>& X, 
				 Vec<double>& d2wall, double* Vwall, 
				 SVec<double,dim>& V, Vec<double>& P)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING) {
    Vec3D n;
    computeNormal(X, n);

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      tets[elemNum].computeGradientP1Function(X, dp1dxj);

    double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
    double* Vface[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};
    double* Vtet[4] = {V[tets[elemNum][0]], V[tets[elemNum][1]], 
		       V[tets[elemNum][2]], V[tets[elemNum][3]]};

    double hp = third * postFcn->computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet);

    for (int j=0; j<3; ++j)
      P[ nodeNum[j] ] += hp;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeForceAndMoment(TetSet &tets, PostFcn *postFcn, SVec<double,3> &X, 
				 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
				Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
				double *nodalForceWeight, int hydro)
{

    Vec3D x[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;

    Vec3D fi0,fi1,fi2,fv;
    computeForce(tets, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv,nodalForceWeight, hydro);

    Fi += fi0 + fi1 + fi2;
    Mi += (dx[0] ^ fi0) + (dx[1] ^ fi1) + (dx[2] ^ fi2);

    Fv += fv;
    Mv += dxm ^ fv;

}

//------------------------------------------------------------------------------
/*
template<int dim>
void Face::computeCenterOfForce(TetSet &tets, PostFcn *postFcn, SVec<double,3> &X,
                                 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                                Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv,
                                double *nodalForceWeight, int hydro)  {

    Vec3D x[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;
    Vec3D xm = third*(x[0]+x[1]+x[2]);

    Vec3D fi0,fi1,fi2,fv;
    computeForce(tets, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv,nodalForceWeight, hydro);

    Fi += fi0 + fi1 + fi2;

    Mi += dx[0] ^ fi0 + dx[1] ^ fi1 + dx[2] ^ fi2;

    Fv = fi0 + fi1 + fi2;
    for (int k = 0; k < 3; k++)
      Fv[k] *= xm[k];
    
    Mv += dxm ^ fv;

}
*/
//------------------------------------------------------------------------------

template<int dim>
double Face::computeInterfaceWork(TetSet& tets, PostFcn* postFcn, SVec<double,3>& X, 
				  Vec<double>& d2wall, double ndot, 
				  double* Vwall, SVec<double,dim>& V, double pin)
{

  double W = 0.0;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING) {
    Vec3D n;
    computeNormal(X, n);

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      tets[elemNum].computeGradientP1Function(X, dp1dxj);

    double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
    double* Vface[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};
    double* Vtet[4] = {V[tets[elemNum][0]], V[tets[elemNum][1]],
		       V[tets[elemNum][2]], V[tets[elemNum][3]]};

    W = postFcn->computeInterfaceWork(dp1dxj, n, ndot, d2w, Vwall, Vface, Vtet, pin);
  }

  return W;

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeScalarQuantity(PostFcn::ScalarType type, TetSet& tets, PostFcn *postFcn, 
				 SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
				 SVec<double,dim> &V, SVec<double,2> &Q)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      tets[elemNum].computeGradientP1Function(X, dp1dxj);

    double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
    double *Vface[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};
    double *Vtet[4] = {V[tets[elemNum][0]], V[tets[elemNum][1]],
		       V[tets[elemNum][2]], V[tets[elemNum][3]]};

    double q = postFcn->computeFaceScalarQuantity(type, dp1dxj, n, d2w, Vwall, Vface, Vtet);

    for (int j=0; j<3; ++j) {
      Q[ nodeNum[j] ][0] += S;
      Q[ nodeNum[j] ][1] += S * q;
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(VarFcn *varFcn, Vec3D &normal, double normalVel,
			   SVec<double,dim> &V, Vec<double> &dt, double beta,
                           double k1, double cmach)
{

  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * normalVel;

  for (int l=0; l<3; ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum[l] ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum[l] ]);
    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum[l] ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    beta = fmin(fmax(k1*locMach, beta), cmach);
    
    double beta2 = beta * beta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

    dt[ nodeNum[l] ] += min(0.5*(coeff1-coeff2), 0.0)* third * S;

  }
    
}
//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec3D &normal, double normalVel,
			   SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, Vec<double> &idtv, double beta,
                           double k1, double cmach)
{

  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * normalVel;

  for (int l=0; l<3; ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum[l] ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum[l] ]);
    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum[l] ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    beta = fmin(fmax(k1*locMach, beta), cmach);
    
    double beta2 = beta * beta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

    idti[ nodeNum[l] ] += min(0.5*(coeff1-coeff2), 0.0)* third * S;

    double vis = 0.0;
    if (fet) vis = fet->computeViscousTimeStep(X[nodeNum[l]],V[nodeNum[l]]);
    idtv[ nodeNum[l] ] += vis*S*S;

  }
    
}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal, 
				   double normalVel, SVec<double,dim> &V, 
				   double *Ub, SVec<double,dim> &fluxes)
{
  if(fluxFcn[code]){
    double flux[dim];
    for (int l=0; l<3; ++l) {
      fluxFcn[code]->compute(0.0, normal, normalVel, V[nodeNum[l]], Ub, flux);
      for (int k=0; k<dim; ++k){
        fluxes[ nodeNum[l] ][k] += third * flux[k];
      }
    }
  }
}

//------------------------------------------------------------------------------
                                                                                                      
template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal,
                                   double normalVel, SVec<double,dim> &V,
                                   double *Ub, Vec<double> &Phi, 
                                   SVec<double,dim> &fluxes)
{
  double flux[dim];

  if(fluxFcn[code]){
    for (int l=0; l<3; ++l) {
      if (Phi[nodeNum[l]] >= 0.0) 
        fluxFcn[code]->compute(0.0, normal, normalVel, V[nodeNum[l]], Ub, flux, 1);
      if (Phi[nodeNum[l]] <  0.0) 
        fluxFcn[code]->compute(0.0, normal, normalVel, V[nodeNum[l]], Ub, flux, -1);
      for (int k=0; k<dim; ++k)
        fluxes[ nodeNum[l] ][k] += third * flux[k];
    }
  }
                                                                                                      
}

//------------------------------------------------------------------------------
                                                                                        
                                                                                        
template<int dim>
inline
void Face::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec3D &normal,
                                   double normalVel, SVec<double,dim> &V,
                                   Vec<double> &Phi, Vec<double> &PhiF)
{
  double flux;
  double PhiF1, Uf;                                                                                                           
  for (int l=0; l<3; ++l) {
//    fluxFcn[code]->computeLS(normal, normalVel, V[nodeNum[l]], Phi[nodeNum[l]], flux);
    Uf   = V[nodeNum[l]][1]*normal[0] + V[nodeNum[l]][2]*normal[1] +V[nodeNum[l]][3]*normal[2];
//    Un   = 0.5*(Uf  +fabs(Uf));
    PhiF1  = Uf*Phi[nodeNum[l]];
    for (int k=0; k<1; ++k)
//      PhiF[ nodeNum[l] ] += third * flux;
      PhiF[ nodeNum[l] ] += third * PhiF1;
    if (isnan(Uf)) fprintf(stderr, "in face computeFiniteVolumeTerm for Uf %d \n",nodeNum[l]);
    if (isnan(PhiF1)) fprintf(stderr, "in face computeFiniteVolumeTerm %d \n",nodeNum[l]);
  }

}
//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal, 
					   double normalVel, SVec<double,dim> &V, 
					   double *Ub, GenMat<Scalar,neq> &A)
{

  double jac[neq*neq];
  for (int l=0; l<3; ++l) {
    fluxFcn[code]->computeJacobian(0.0, normal, normalVel, V[nodeNum[l]], Ub, jac);
    Scalar *Aii = A.getElem_ii(nodeNum[l]);
    for (int k=0; k<neq*neq; ++k) 
      Aii[k] += third * jac[k];
  }


}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec3D &normal,
                                           double normalVel, SVec<double,dim> &V,
                                           double *Ub, GenMat<Scalar,neq> &A, int* nodeType)
{

  double jac[neq*neq];
  for (int l=0; l<3; ++l) {
    if(!(code == BC_INLET_MOVING || code == BC_OUTLET_MOVING ||
         code == BC_INLET_FIXED  || code == BC_OUTLET_FIXED)) {

      fluxFcn[code]->computeJacobian(0.0, normal, normalVel, V[nodeNum[l]], Ub, jac);
      Scalar *Aii = A.getElem_ii(nodeNum[l]);
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += third * jac[k];

    }
  }
}

//------------------------------------------------------------------------------
template<int dim, class Scalar>
inline
void Face::computeJacobianFiniteVolumeTermLS(Vec3D &normal,
                                           double normalVel, SVec<double,dim> &V,
                                           GenMat<Scalar,1> &A)
{

  double jac;
  for (int l=0; l<3; ++l) {
//    fluxFcn[code]->computeJacobian(normal, normalVel, V[nodeNum[l]], Ub, jac);
    jac  = V[nodeNum[l]][1]*normal[0]  + V[nodeNum[l]][2]*normal[1]  + V[nodeNum[l]][3]*normal[2];
    jac *= V[nodeNum[l]][0];
    Scalar *Aii = A.getElem_ii(nodeNum[l]);
    for (int k=0; k<1; ++k)
      Aii[k] += third * jac;
  }
}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeGalerkinTerm(TetSet &tets, FemEquationTerm *fet, SVec<double,3> &X, 
			       Vec<double> &d2wall, double *Vwall,
			       SVec<double,dim> &V, SVec<double,dim> &R)
{

  if (!fet->doesFaceTermExist(code)) return;

  Vec3D n;
  computeNormal(X, n);

  if (fet->doesFaceNeedGradientP1Function())
    tets[elemNum].computeFaceGalerkinTerm(fet, nodeNum, code, n, X, d2wall, Vwall, V, R);
  else {
    double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
    double *v[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};
    double r[dim];
    fet->computeSurfaceTerm(code, n, d2w, Vwall, v, r);

    for (int l=0; l<3; ++l)
      for (int k=0; k<dim; ++k)
	R[ nodeNum[l] ][k] -= third * r[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianGalerkinTerm(TetSet &tets, FemEquationTerm *fet, 
				       SVec<double,3> &X, Vec<double> &ctrlVol,
				       Vec<double> &d2wall, double *Vwall, 
				       SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  if (!fet->doesFaceTermExist(code)) return;

  Vec3D n;
  computeNormal(X, n);

  if (fet->doesFaceNeedGradientP1Function())
    tets[elemNum].computeFaceJacobianGalerkinTerm(fet, nodeNum, code, n, X, ctrlVol, d2wall, Vwall, V, A);
  else {
    double d2w[3] = {d2wall[nodeNum[0]], d2wall[nodeNum[1]], d2wall[nodeNum[2]]};
    double *v[3] = {V[nodeNum[0]], V[nodeNum[1]], V[nodeNum[2]]};

    double dRdU[3][neq*neq];
    fet->computeJacobianSurfaceTerm(code, n, d2w, Vwall, v, reinterpret_cast<double *>(dRdU));

    for (int k=0; k<3; ++k) {
      Scalar *Aii = A.getElem_ii(nodeNum[k]);
      for (int m=0; m<neq*neq; ++m)
	Aii[m] -= third * dRdU[k][m];
    }

    for (int l=0; l<3; ++l) {
      int i, j;
      if (nodeNum[ edgeEnd[l][0] ] < nodeNum[ edgeEnd[l][1] ]) {
	i = edgeEnd[l][0];
	j = edgeEnd[l][1];
      } 
      else {
	i = edgeEnd[l][1];
	j = edgeEnd[l][0];
      }

      Scalar *Aij = A.getElem_ij(edgeNum[l]);
      Scalar *Aji = A.getElem_ji(edgeNum[l]);

      if (Aij) {

        double cij = third / ctrlVol[ nodeNum[i] ];
	for (int m=0; m<neq*neq; ++m)
	  Aij[m] -= cij * dRdU[j][m];
      }

      if (Aji) {

	double cji = third / ctrlVol[ nodeNum[j] ];
	for (int m=0; m<neq*neq; ++m)
	  Aji[m] -= cji * dRdU[i][m];
      }
    }
    
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,dim> &V, Vec<double> &dt,
                              double beta, double k1, double cmach)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i].computeTimeStep(varFcn, n[i], ndot[i], V, dt, beta, k1, cmach);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti,Vec<double> &idtv,
                              double beta, double k1, double cmach)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i].computeTimeStep(fet, varFcn, n[i], ndot[i], X, V, idti, idtv, beta, k1, cmach);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V, 
				      SVec<double,dim> &fluxes)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i)  {
    faces[i].computeFiniteVolumeTerm(fluxFcn, n[i], ndot[i], V, Ub[i], fluxes);
  }
}

//------------------------------------------------------------------------------
template<int dim>
void FaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,dim> &V,
                                      Vec<double> &Phi, SVec<double,dim> &fluxes)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                                                                                                       
  for (int i=0; i<numFaces; ++i)  {
    faces[i].computeFiniteVolumeTerm(fluxFcn, n[i], ndot[i], V, Ub[i], Phi, fluxes);
  }

                                                                                                       
}

//------------------------------------------------------------------------------
                                                                                                                                        
template<int dim>
void FaceSet::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,dim> &V,
                                      Vec<double> &Phi, Vec<double> &PhiF)
{
                                                                                                                                        
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                                                                                                                                        
  for (int i=0; i<numFaces; ++i)  {
    faces[i].computeFiniteVolumeTermLS(fluxFcn, n[i], ndot[i], V, Phi, PhiF);
  }
                                                                                                                                        
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
					      GeoState &geoState, SVec<double,dim> &V, 
					      GenMat<Scalar,neq> &A)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i) 
    faces[i].computeJacobianFiniteVolumeTerm(fluxFcn, n[i], ndot[i], V, Ub[i], A);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                              GeoState &geoState, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, int* nodeType)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                                                                                                  
  for (int i=0; i<numFaces; ++i)
    faces[i].computeJacobianFiniteVolumeTerm(fluxFcn, n[i], ndot[i], V, Ub[i], A, nodeType);
                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                                                        
template<int dim, class Scalar>
void FaceSet::computeJacobianFiniteVolumeTermLS(
                                              GeoState &geoState, SVec<double,dim> &V,
                                              GenMat<Scalar,1> &A)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  for (int i=0; i<numFaces; ++i)
    faces[i].computeJacobianFiniteVolumeTermLS(n[i], ndot[i], V, A);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeGalerkinTerm(TetSet &tets, FemEquationTerm *fet, BcData<dim> &bcData, 
				  GeoState &geoState, SVec<double,3> &X, 
				  SVec<double,dim> &V, SVec<double,dim> &R)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i].computeGalerkinTerm(tets, fet, X, d2wall, Vwall[i], V, R);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianGalerkinTerm(TetSet &tets, FemEquationTerm *fet, 
					  BcData<dim> &bcData, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol, 
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i].computeJacobianGalerkinTerm(tets, fet, X, ctrlVol, d2wall, Vwall[i], V, A);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeForceDerivs(TetSet &tets, VarFcn *varFcn, SVec<double,3> &X, 
                              SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
                              SVec<double,3> **localMX)  {

  static double third = 1.0/3.0;
  int j, k;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
   || code == BC_SLIP_WALL_MOVING) {

    Vec3D F[dim];

    Vec3D x[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};

    computeFDerivs(tets, varFcn, X, V, F);

    double ducg[dim];
    for (k = 0; k < dim; ++k) {
      ducg[k] = 0;
      for ( j=0; j<3; ++j)
        ducg[k] += third*(deltaU[ nodeNum[j] ][k]);
    }

    int nStrMode = modalF.len;

    for (int iMode = 0; iMode < nStrMode; iMode++)  {

      SVec<double, 3> *locMX = localMX[iMode];
      for ( j=0; j<3; ++j) {
        for ( k=0; k<dim; ++k) {
          modalF[iMode] += third*(F[k][0]*((*locMX)[nodeNum[j] ][0]))*ducg[k];
          modalF[iMode] += third*(F[k][1]*((*locMX)[nodeNum[j] ][1]))*ducg[k];
          modalF[iMode] += third*(F[k][2]*((*locMX)[nodeNum[j] ][2]))*ducg[k];
        }
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, TetSet &tets, SVec<double,3> &X, 
                   SVec<double,dim> &V, Vec<double> &d2wall, SVec<double, dim> &Vwall, double pInfty, 
		   Vec3D &CFi, Vec3D &CMi, Vec3D &CFv, Vec3D &CMv, double *nodalForceWeight) 
{

  static double third = 1.0/3.0;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
   || code == BC_SLIP_WALL_MOVING) {

    Vec3D x[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};
    Vec3D xcg = third * (x[0] + x[1] + x[2]);
    Vec3D Fi0,Fi1,Fi2,Fv;

    double *vWall = reinterpret_cast<double *>(Vwall.data());
    computeForce(tets, postFcn, X, d2wall, vWall, V, &pInfty, Fi0,Fi1,Fi2,Fv,nodalForceWeight);

    CFi += Fi0 + Fi1 + Fi2;
    CMi += (x[0] - x0) ^ Fi0 + (x[1] - x0) ^ Fi1 + (x[2] - x0) ^ Fi2;

    CFv += Fv;
    CMv += (xcg - x0) ^ Fv;

  }

}

//------------------------------------------------------------------------------


template<int dim>
inline
void Face::computeFDerivs(TetSet &tets, VarFcn *varFcn, SVec<double,3> &X, 
                          SVec<double,dim> &Vgl, Vec3D (*F))  {

  static double third = 1.0/3.0;

  for(int i=0; i<5; ++i)
    F[i] = 0.0;

  Vec3D x[3] = {X[nodeNum[0]], X[nodeNum[1]], X[nodeNum[2]]};

  Vec3D n = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

  Vec3D vel = third * (  varFcn->getVelocity(Vgl[nodeNum[0]]) +
                         varFcn->getVelocity(Vgl[nodeNum[1]]) +
                         varFcn->getVelocity(Vgl[nodeNum[2]]) );

  double dens = third * (  varFcn->getDensity(Vgl[nodeNum[0]]) +
                         varFcn->getDensity(Vgl[nodeNum[1]]) +
                         varFcn->getDensity(Vgl[nodeNum[2]]) );

  double gam1 = varFcn->getGamma1();

  F[0] = (0.5*gam1*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])) * n;
  F[1] = (-gam1*vel[0]) * n;
  F[2] = (-gam1*vel[1]) * n;
  F[3] = (-gam1*vel[2]) * n;
  F[4] = ( gam1) * n;

}

//------------------------------------------------------------------------------
