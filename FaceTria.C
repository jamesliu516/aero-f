#include <FaceTria.h>

#include <FluxFcnDescPerfectGas.h>
#include <FluxFcnDescWaterCompressible.h>
#include <FluxFcnDescGasInGas.h>
#include <FluxFcnDescLiquidInLiquid.h>
#include <FluxFcnDescGasInLiquid.h>
#include <FemEquationTerm.h>
#include <BcData.h>
#include <BcDef.h>
#include <Elem.h>
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

template<int dim>
inline
void FaceTria::computeForce(ElemSet &elems,
			    PostFcn *postFcn, SVec<double,3> &X, 
			    Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, 
			    double *nodalForceWeight, int hydro)
{

  Vec3D n;
  computeNormal(X, n);
  Elem& elem = elems[elemNum];

  double dp1dxj[4][3];
  if (postFcn->doesFaceNeedGradientP1Function())
    elem.computeGradientP1Function(X, dp1dxj);

  double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
  double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
  double *Vtet[4] = {V[elem[0]], V[elem[1]], 
		     V[elem[2]], V[elem[3]]};
  double *Xface[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};

  postFcn->computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, 
			Fi0, Fi1, Fi2, Fv, nodalForceWeight, hydro);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeNodalForce(ElemSet &elems,
				 PostFcn *postFcn, SVec<double,3> &X, 
				 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
				 double pin, SVec<double,3> &F, double *nodalForceWeight)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
    || code == BC_SLIP_WALL_MOVING) {
    Vec3D Fi0, Fi1, Fi2, Fv;
    computeForce(elems, postFcn, X, d2wall, Vwall, V, &pin, Fi0, Fi1, Fi2, Fv, nodalForceWeight );
    
    Vec3D Ftot[3];
    Ftot[0] = Fi0 + third*Fv; 
    Ftot[1] = Fi1 + third*Fv;
    Ftot[2] = Fi2 + third*Fv;

    for (int j=0; j<3; ++j) {
      F[ nodeNum(j) ][0] += Ftot[j][0];
      F[ nodeNum(j) ][1] += Ftot[j][1];
      F[ nodeNum(j) ][2] += Ftot[j][2];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeNodalHeatPower(ElemSet& elems,
				     PostFcn* postFcn, SVec<double,3>& X, 
				     Vec<double>& d2wall, double* Vwall, 
				     SVec<double,dim>& V, Vec<double>& P)
{
  if (code == BC_ISOTHERMAL_WALL_MOVING) {
    Vec3D n;
    computeNormal(X, n);
    Elem& elem = elems[elemNum];

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      elem.computeGradientP1Function(X, dp1dxj);

    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double* Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double* Vtet[4] = {V[elem[0]], V[elem[1]], 
		       V[elem[2]], V[elem[3]]};

    double hp = third * postFcn->computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet);

    for (int j=0; j<3; ++j)
      P[ nodeNum(j) ] += hp;
  }
}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, 
				     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
				     SVec<double,dim> &V, Vec3D &x0, Vec3D &Fi, Vec3D &Mi, 
				     Vec3D &Fv, Vec3D &Mv, double *nodalForceWeight, int hydro)
{
    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;

    Vec3D fi0,fi1,fi2,fv;
    computeForce(elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv,nodalForceWeight, hydro);

    Fi += fi0 + fi1 + fi2;
    Mi += dx[0] ^ fi0 + dx[1] ^ fi1 + dx[2] ^ fi2;

    Fv += fv;
    Mv += dxm ^ fv;
}

//------------------------------------------------------------------------------

template<int dim>
double FaceTria::computeInterfaceWork(ElemSet& elems, PostFcn* postFcn, 
				      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
				      double* Vwall, SVec<double,dim>& V, double pin)
{
  double W = 0.0;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING) {
    Vec3D n;
    computeNormal(X, n);
    Elem& elem = elems[elemNum];

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      elem.computeGradientP1Function(X, dp1dxj);

    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double* Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double* Vtet[4] = {V[elem[0]], V[elem[1]],
		       V[elem[2]], V[elem[3]]};

    W = postFcn->computeInterfaceWork(dp1dxj, n, ndot, d2w, Vwall, Vface, Vtet, pin);
  }

  return W;
}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeScalarQuantity(PostFcn::ScalarType stype, ElemSet& elems,
				     PostFcn *postFcn, SVec<double,3> &X, 
				     Vec<double> &d2wall, double *Vwall, 
				     SVec<double,dim> &V, SVec<double,2> &Q)
{
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);
    Elem& elem = elems[elemNum];

    double dp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function())
      elem.computeGradientP1Function(X, dp1dxj);
    
    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double *Vtet[4] = {V[elem[0]], V[elem[1]],
		       V[elem[2]], V[elem[3]]};

    double q = postFcn->computeFaceScalarQuantity(stype, dp1dxj, n, d2w, Vwall, Vface, Vtet);

    for (int j=0; j<3; ++j) {
      Q[ nodeNum(j) ][0] += S;
      Q[ nodeNum(j) ][1] += S * q;
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
//inline
void FaceTria::computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &d2wall, double *Vwall,
				   SVec<double,dim> &V, SVec<double,dim> &R)
{
  if (!fet->doesFaceTermExist(code)) return;

  Vec3D n;
  computeNormal(X, n);

  if (fet->doesFaceNeedGradientP1Function())
    elems[elemNum].computeFaceGalerkinTerm(fet, nodeNum(), code, n, X, d2wall, Vwall, V, R);
  else {
    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double *v[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double r[dim];
    fet->computeSurfaceTerm(code, n, d2w, Vwall, v, r);

    for (int l=0; l<3; ++l)
      for (int k=0; k<dim; ++k)
	R[ nodeNum(l) ][k] -= third * r[k];
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
//inline
void FaceTria::computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
					   SVec<double,3> &X, Vec<double> &ctrlVol,
					   Vec<double> &d2wall, double *Vwall, 
					   SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{
  if (!fet->doesFaceTermExist(code)) return;

  Vec3D n;
  computeNormal(X, n);

  if (fet->doesFaceNeedGradientP1Function()) {
    elems[elemNum].computeFaceJacobianGalerkinTerm(fet, nodeNum(), code, n, X, ctrlVol, d2wall, Vwall, V, A);
  } else {
    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double *v[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};

    double dRdU[3][neq*neq];
    fet->computeJacobianSurfaceTerm(code, n, d2w, Vwall, v, reinterpret_cast<double *>(dRdU));

    for (int k=0; k<3; ++k) {
      Scalar *Aii = A.getElem_ii(nodeNum(k));
      for (int m=0; m<neq*neq; ++m)
	Aii[m] -= third * dRdU[k][m];
    }

    for (int l=0; l<3; ++l) {
      int i, j;
      if (nodeNum( edgeEnd(l,0) ) < nodeNum( edgeEnd(l,1) )) {
	i = edgeEnd(l,0);
	j = edgeEnd(l,1);
      } 
      else {
	i = edgeEnd(l,1);
	j = edgeEnd(l,0);
      }

      Scalar *Aij = A.getElem_ij(edgeNum(l));
      Scalar *Aji = A.getElem_ji(edgeNum(l));

      if (Aij) {

        double cij = third / ctrlVol[ nodeNum(i) ];
	for (int m=0; m<neq*neq; ++m)
	  Aij[m] -= cij * dRdU[j][m];
      }

      if (Aji) {

	double cji = third / ctrlVol[ nodeNum(j) ];
	for (int m=0; m<neq*neq; ++m)
	  Aji[m] -= cji * dRdU[i][m];
      }
    }
    
  }
}

//------------------------------------------------------------------------------

template<int dim>
//inline
void FaceTria::computeForceDerivs(ElemSet &elems, VarFcn *varFcn, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &deltaU, Vec<double> &modalF, 
				  SVec<double,3> **localMX)  
{
  static double third = 1.0/3.0;
  int j, k;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
   || code == BC_SLIP_WALL_MOVING) {

    Vec3D F[dim];

    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};

    computeFDerivs(elems, varFcn, X, V, F);

    double ducg[dim];
    for (k = 0; k < dim; ++k) {
      ducg[k] = 0;
      for ( j=0; j<3; ++j)
        ducg[k] += third*(deltaU[ nodeNum(j) ][k]);
    }

    int nStrMode = modalF.len;

    for (int iMode = 0; iMode < nStrMode; iMode++)  {

      SVec<double, 3> *locMX = localMX[iMode];
      for ( j=0; j<3; ++j) {
        for ( k=0; k<dim; ++k) {
          modalF[iMode] += third*(F[k][0]*((*locMX)[nodeNum(j) ][0]))*ducg[k];
          modalF[iMode] += third*(F[k][1]*((*locMX)[nodeNum(j) ][1]))*ducg[k];
          modalF[iMode] += third*(F[k][2]*((*locMX)[nodeNum(j) ][2]))*ducg[k];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
//inline
void FaceTria::computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems,
					SVec<double,3> &X, 
					SVec<double,dim> &V, Vec<double> &d2wall, 
					SVec<double, dim> &Vwall, double pInfty, 
					Vec3D &CFi, Vec3D &CMi, Vec3D &CFv, Vec3D &CMv, 
					double *nodalForceWeight) 
{
  static double third = 1.0/3.0;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
   || code == BC_SLIP_WALL_MOVING) {

    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D xcg = third * (x[0] + x[1] + x[2]);
    Vec3D Fi0,Fi1,Fi2,Fv;

    double *vWall = reinterpret_cast<double *>(Vwall.data());
    computeForce(elems, postFcn, X, d2wall, vWall, V, &pInfty, Fi0,Fi1,Fi2,Fv,nodalForceWeight);

    CFi += Fi0 + Fi1 + Fi2;
    CMi += (x[0] - x0) ^ Fi0 + (x[1] - x0) ^ Fi1 + (x[2] - x0) ^ Fi2;

    CFv += Fv;
    CMv += (xcg - x0) ^ Fv;

  }
}

//------------------------------------------------------------------------------

template<int dim>
//inline
void FaceTria::computeFDerivs(ElemSet &elems,
			      VarFcn *varFcn, SVec<double,3> &X, 
			      SVec<double,dim> &Vgl, Vec3D (*F))  
{
  static double third = 1.0/3.0;

  for(int i=0; i<5; ++i)
    F[i] = 0.0;

  Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};

  Vec3D n = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

  Vec3D vel = third * (  varFcn->getVelocity(Vgl[nodeNum(0)]) +
                         varFcn->getVelocity(Vgl[nodeNum(1)]) +
                         varFcn->getVelocity(Vgl[nodeNum(2)]) );

  double dens = third * (  varFcn->getDensity(Vgl[nodeNum(0)]) +
                         varFcn->getDensity(Vgl[nodeNum(1)]) +
                         varFcn->getDensity(Vgl[nodeNum(2)]) );

  double gam1 = varFcn->getGamma1();

  F[0] = (0.5*gam1*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])) * n;
  F[1] = (-gam1*vel[0]) * n;
  F[2] = (-gam1*vel[1]) * n;
  F[3] = (-gam1*vel[2]) * n;
  F[4] = ( gam1) * n;
}

//------------------------------------------------------------------------------
