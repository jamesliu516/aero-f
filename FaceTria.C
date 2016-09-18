#include <FaceTria.h>

#include <FluxFcn.h>
#include <FemEquationTerm.h>
#include <BcData.h>
#include <BcDef.h>
#include <Elem.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <VectorSet.h>
#include <RectangularSparseMatrix.h>

#include <cmath>

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
			    double* gradP[3], int hydro)
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

// Vamshi
   double dPdx[3][3];
   for(int i=0; i<3; i++) // i represents node
   {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
   }

  postFcn->computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, 
			Fi0, Fi1, Fi2, Fv, dPdx, hydro);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FaceTria::computeForce(ExactRiemannSolver<dim>& riemann, 
                            VarFcn *varFcn, Vec<Vec3D> &normals, Vec<double> &normalVel, ElemSet &elems,
                            PostFcn *postFcn, SVec<double,3> &X,
                            Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                            double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv,
                            double* gradP[3], int hydro)
{

  Vec3D n;
  computeNormal(X, n);
  Elem& elem = elems[elemNum];

  double dp1dxj[4][3];
  if (postFcn->doesFaceNeedGradientP1Function())
    elem.computeGradientP1Function(X, dp1dxj);

  double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};

  //compute new Vface and plug into Vtet.
  //need: normalVel, normal, varFcn.
  double Wstar[3][2*dim];
  double Vi[2*dim];
  for (int l=0; l<3; ++l) {
    int k = nodeNum(l);
    for(int iDim=0; iDim<dim; iDim++)
      Vi[iDim] = Vi[iDim+dim] = V[k][iDim];
    Vec3D unitNormal = getNormal(normals, l)/(getNormal(normals,l).norm());
    Vec3D wallVel = getNormalVel(normalVel, l)/(getNormal(normals,l).norm())*unitNormal;
    riemann.computeFSIRiemannSolution(Vi, wallVel, -1.0*unitNormal, varFcn, Wstar[l], 0); 
  }

  double *Vface[3] = {Wstar[0], Wstar[1], Wstar[2]};
  double *Vtet[4] = {V[elem[0]], V[elem[1]],
                     V[elem[2]], V[elem[3]]};
  for (int i=0; i<4; i++)
    for(int j=0; j<3; j++)
      if (elem[i]==nodeNum(j))
        Vtet[i] = Wstar[j];

//  double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
//  double *Vtet[4] = {V[elem[0]], V[elem[1]],
//                     V[elem[2]], V[elem[3]]};


  double *Xface[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};

// Vamshi
   double dPdx[3][3];
   for(int i=0; i<3; i++) // i represents node
   {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
   }

  postFcn->computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin,
                        Fi0, Fi1, Fi2, Fv, dPdx, hydro);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void FaceTria::computeDerivativeOfForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall,
                                    double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
                                    double dS[3], double *pin, Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, 
			            double* gradP[3], double* dGradP[3], int hydro)
{

  Vec3D n;
  Vec3D dn;

  computeNormalAndDerivative(X, dX, n, dn);
  Elem& elem = elems[elemNum];

  double dp1dxj[4][3], ddp1dxj[4][3];
  if (postFcn->doesFaceNeedGradientP1Function()) {
    elem.computeGradientP1Function(X, dp1dxj);
    elem.computeDerivativeOfGradientP1Function(X, dX, ddp1dxj);
  }

  double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
  double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
  double *Vtet[4] = {V[elem[0]], V[elem[1]],
		     V[elem[2]], V[elem[3]]};
  double *Xface[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
  double *dVface[3] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)]};
  double *dVtet[4] = {dV[elem[0]], dV[elem[1]],
		     dV[elem[2]], dV[elem[3]]};
  double *dXface[3] = {dX[nodeNum(0)], dX[nodeNum(1)], dX[nodeNum(2)]};

  double dPdx[3][3], ddPdx[3][3];
  for(int i=0; i<3; i++) // i represents node
  {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
    ddPdx[i][0] = dGradP[0][nodeNum(i)];
    ddPdx[i][1] = dGradP[1][nodeNum(i)];
    ddPdx[i][2] = dGradP[2][nodeNum(i)];
  }

  postFcn->computeDerivativeOfForce(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx, ddPdx, hydro);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FaceTria::computeForceTransmitted(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
				       Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
				       double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv,
				       double* gradP[3], int hydro)
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

// Vamshi
   double dPdx[3][3];
   for(int i=0; i<3; i++) // i represents node
   {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
   }
   postFcn->computeForceTransmitted(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FaceTria::computeDerivativeOfForceTransmitted(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall,
                                    double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
                                    double dS[3], double *pin, Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, 
                                    double* gradP[3], double* dGradP[3], int hydro)
{

  Vec3D n;
  Vec3D dn;

  computeNormalAndDerivative(X, dX, n, dn);
  Elem& elem = elems[elemNum];
  double dndX[3][3][3] = {0};
  compute_dndX(X, dndX);

  double dp1dxj[4][3], ddp1dxj[4][3];
  if (postFcn->doesFaceNeedGradientP1Function()) {
    elem.computeGradientP1Function(X, dp1dxj);
    elem.computeDerivativeOfGradientP1Function(X, dX, ddp1dxj);
  }

  double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
  double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
  double *Vtet[4] = {V[elem[0]], V[elem[1]], V[elem[2]], V[elem[3]]};
  double *Xface[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
  double *dVface[3] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)]};
  double dVfacedV[3][dim] = {0};
  for(int i=0; i<3; ++i) for(int j=0; j<dim; ++j) dVfacedV[i][j] = 1.0;
  double *dVtet[4] = {dV[elem[0]], dV[elem[1]], dV[elem[2]], dV[elem[3]]};
  double *dXface[3] = {dX[nodeNum(0)], dX[nodeNum(1)], dX[nodeNum(2)]};
  double dXfacedX[3][3] = {0};
  for(int i=0; i<3; i++) for(int j=0; j<3; ++j) dXfacedX[i][j] = 1.0;


  double dPdx[3][3], ddPdx[3][3], ddPdxdGradP[3][3] = {0};
  for(int i=0; i<3; i++) // i represents node
  {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
    ddPdx[i][0] = dGradP[0][nodeNum(i)];
    ddPdx[i][1] = dGradP[1][nodeNum(i)];
    ddPdx[i][2] = dGradP[2][nodeNum(i)];
    for(int j=0; j<3; ++j)
      ddPdxdGradP[i][j] = 1.0;
  }

//  Vec3D dFi00(dFi0), dFi11(dFi1), dFi22(dFi2), dFvv(dFv), diff0, diff1, diff2, diffv;
  postFcn->computeDerivativeOfForceTransmitted(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx, ddPdx, hydro);

  double dFi0dn[3] = {0}, dFi0dS[3][3] = {0}, dFi0dVface[3][3][5] = {0}, dFi0ddPdx[3][3][3] = {0}, dFi0dXface[3][3][3] = {0};
  double dFi1dn[3] = {0}, dFi1dS[3][3] = {0}, dFi1dVface[3][3][5] = {0}, dFi1ddPdx[3][3][3] = {0}, dFi1dXface[3][3][3] = {0};
  double dFi2dn[3] = {0}, dFi2dS[3][3] = {0}, dFi2dVface[3][3][5] = {0}, dFi2ddPdx[3][3][3] = {0}, dFi2dXface[3][3][3] = {0};

  postFcn->computeDerivativeOperatorsOfForceTransmitted(Xface, n, Vface, pin, dPdx, hydro, 
                                                        dFi0dn, dFi0dS, dFi0dVface, dFi0ddPdx, dFi0dXface,
                                                        dFi1dn, dFi1dS, dFi1dVface, dFi1ddPdx, dFi1dXface,
                                                        dFi2dn, dFi2dS, dFi2dVface, dFi2ddPdx, dFi2dXface);
/*
 dFi00 = 0.0;
 for(int l=0; l<3; ++l) {
//   dFi00[l] += dFi0dn[l]*dn[l];
   for(int j=0; j<3; ++j) {
     dFi00[l] += dFi0dS[l][j]*dS[j];
     for(int k=0; k<5; ++k)
       dFi00[l] += dFi0dVface[l][j][k]*dVfacedV[j][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi00[l] += dFi0ddPdx[l][j][k]*ddPdxdGradP[j][k]*dGradP[k][nodeNum(j)] + dFi0dXface[l][j][k]*dXfacedX[j][k]*dX[nodeNum(j)][k];
       dFi00[l] += dFi0dn[l]*dndX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

 dFi11 = 0.0;
 for(int l=0; l<3; ++l) {
//   dFi11[l] += dFi1dn[l]*dn[l];
   for(int j=0; j<3; ++j) {
     dFi11[l] += dFi1dS[l][j]*dS[j];
     for(int k=0; k<5; ++k)
       dFi11[l] += dFi1dVface[l][j][k]*dVfacedV[j][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi11[l] += dFi1ddPdx[l][j][k]*ddPdxdGradP[j][k]*dGradP[k][nodeNum(j)] + dFi1dXface[l][j][k]*dXfacedX[j][k]*dX[nodeNum(j)][k];
       dFi11[l] += dFi1dn[l]*dndX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

 dFi22 = 0.0;
 for(int l=0; l<3; ++l) {
//   dFi22[l] += dFi2dn[l]*dn[l];
   for(int j=0; j<3; ++j) {
     dFi22[l] += dFi2dS[l][j]*dS[j];
     for(int k=0; k<5; ++k)
       dFi22[l] += dFi2dVface[l][j][k]*dVfacedV[j][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi22[l] += dFi2ddPdx[l][j][k]*ddPdxdGradP[j][k]*dGradP[k][nodeNum(j)] + dFi2dXface[l][j][k]*dXfacedX[j][k]*dX[nodeNum(j)][k];
       dFi22[l] += dFi2dn[l]*dndX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

  dFvv = 0.0;
  diff0 = dFi0 - dFi00;
  diff1 = dFi1 - dFi11;  
  diff2 = dFi2 - dFi22;
  diffv = dFv - dFvv;
  double diff0norm = diff0.norm();
  double diff1norm = diff1.norm();
  double diff2norm = diff2.norm();
  double diffvnorm = diffv.norm();
  double dFi0norm = dFi0.norm();
  double dFi1norm = dFi1.norm();
  double dFi2norm = dFi2.norm();
  double dFvnorm = dFv.norm();
  
  if(dFi0norm != 0) {
    double reldiff = diff0norm/dFi0norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi0 is %e\n", diff0norm/dFi0norm);
  } else {
    if(diff0norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi0 is %e\n", diff0norm);
  }

  if(dFi1norm != 0) {
    double reldiff = diff1norm/dFi1norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi1 is %e\n", diff1norm/dFi1norm);
  } else {
    if(diff1norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi1 is %e\n", diff1norm);
  }

  if(dFi2norm != 0) {
    double reldiff = diff2norm/dFi2norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi2 is %e\n", diff2norm/dFi2norm);
  } else {
    if(diff2norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi2 is %e\n", diff2norm);
  }

  if(dFvnorm != 0) {
    double reldiff = diffvnorm/dFvnorm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFv is %e\n", diffvnorm/dFvnorm);
  } else {
    if(diffvnorm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFv is %e\n", diffvnorm);
  }
*/
}

//------------------------------------------------------------------------------

template<int dim>
inline
void FaceTria::computeDerivativeOperatorsOfForceTransmitted(PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, 
                                                            double *pin, double* gradP[3], int hydro,
                                                            double dFi0dV[3][3][dim], double dFi0dGradP[3][3][3], double dFi0dX[3][3][3],
                                                            double dFi0dn[3], double dFi0dS[3][3],
                                                            double dFi1dV[3][3][dim], double dFi1dGradP[3][3][3], double dFi1dX[3][3][3],
                                                            double dFi1dn[3], double dFi1dS[3][3],
                                                            double dFi2dV[3][3][dim], double dFi2dGradP[3][3][3], double dFi2dX[3][3][3],
                                                            double dFi2dn[3], double dFi2dS[3][3])
{

  Vec3D n;
  double dndX[3][3][3] = {0};
  computeNormal(X, n);
  compute_dndX(X, dndX);

  double *Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
  double *Xface[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
  double dVfacedV[3][dim] = {0}, dXfacedX[3][3] = {0};
  for(int i=0; i<3; ++i) for(int j=0; j<dim; ++j) dVfacedV[i][j] = 1.0;
  for(int i=0; i<3; i++) for(int j=0; j<3; ++j) dXfacedX[i][j] = 1.0;


  double dPdx[3][3], ddPdx[3][3], ddPdxdGradP[3][3] = {0};
  for(int i=0; i<3; i++) // i represents node
  {
    dPdx[i][0] = gradP[0][nodeNum(i)];
    dPdx[i][1] = gradP[1][nodeNum(i)];
    dPdx[i][2] = gradP[2][nodeNum(i)];
    for(int j=0; j<3; ++j)
      ddPdxdGradP[i][j] = 1.0;
  }

  double dFi0dVface[3][3][5] = {0}, dFi0ddPdx[3][3][3] = {0}, dFi0dXface[3][3][3] = {0};
  double dFi1dVface[3][3][5] = {0}, dFi1ddPdx[3][3][3] = {0}, dFi1dXface[3][3][3] = {0};
  double dFi2dVface[3][3][5] = {0}, dFi2ddPdx[3][3][3] = {0}, dFi2dXface[3][3][3] = {0};

  postFcn->computeDerivativeOperatorsOfForceTransmitted(Xface, n, Vface, pin, dPdx, hydro, 
                                                        dFi0dn, dFi0dS, dFi0dVface, dFi0ddPdx, dFi0dXface,
                                                        dFi1dn, dFi1dS, dFi1dVface, dFi1ddPdx, dFi1dXface,
                                                        dFi2dn, dFi2dS, dFi2dVface, dFi2ddPdx, dFi2dXface);
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     for(int k=0; k<dim; ++k) {
       dFi0dV[j][l][k] = dFi0dVface[l][j][k]*dVfacedV[j][k];
     }
     for(int k=0; k<3; ++k) {
       dFi0dGradP[j][l][k] = dFi0ddPdx[l][j][k]*ddPdxdGradP[j][k];
       dFi0dX[j][l][k] = dFi0dXface[l][j][k]*dXfacedX[j][k] +  dFi0dn[l]*dndX[j][l][k];
     }
   }
 }

 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     for(int k=0; k<dim; ++k) {
       dFi1dV[j][l][k] = dFi1dVface[l][j][k]*dVfacedV[j][k];
     }
     for(int k=0; k<3; ++k) {
       dFi1dGradP[j][l][k] = dFi1ddPdx[l][j][k]*ddPdxdGradP[j][k];
       dFi1dX[j][l][k] = dFi1dXface[l][j][k]*dXfacedX[j][k] + dFi1dn[l]*dndX[j][l][k];
     }
   }
 }

 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     for(int k=0; k<dim; ++k) {
       dFi2dV[j][l][k] = dFi2dVface[l][j][k]*dVfacedV[j][k];
     }
     for(int k=0; k<3; ++k) {
       dFi2dGradP[j][l][k] = dFi2ddPdx[l][j][k]*ddPdxdGradP[j][k];
       dFi2dX[j][l][k] = dFi2dXface[l][j][k]*dXfacedX[j][k] + dFi2dn[l]*dndX[j][l][k];
     }
   }
 }

}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeNodalForce(ElemSet &elems,
				 PostFcn *postFcn, SVec<double,3> &X, 
				 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
				 double pin, SVec<double,3> &F, double* gradP[3])
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
    || code == BC_SLIP_WALL_MOVING || code == BC_POROUS_WALL_MOVING) {
    Vec3D Fi0, Fi1, Fi2, Fv;

    computeForceTransmitted(elems, postFcn, X, d2wall, Vwall, V, &pin, Fi0, Fi1, Fi2, Fv, gradP);
    
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

// Included (MB)
template<int dim>
void FaceTria::computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3], double* dGradP[3])
{

  SVec<double,3> dF2(dF), diff2(dF);

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
    || code == BC_SLIP_WALL_MOVING || code == BC_POROUS_WALL_MOVING) {

    Vec3D dFi0, dFi1, dFi2, dFv;
    Vec3D dFi00(dFi0), dFi11(dFi1), dFi22(dFi2), dFvv(dFv), diff0, diff1, diff2, diffv;
    computeDerivativeOfForceTransmitted(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, 0, dFi0, dFi1, dFi2, dFv, gradP, dGradP);

    double dFi0dn[3] = {0}, dFi0dS[3][3] = {0}, dFi0dV[3][3][dim] = {0}, dFi0dGradP[3][3][3] = {0}, dFi0dX[3][3][3] = {0};
    double dFi1dn[3] = {0}, dFi1dS[3][3] = {0}, dFi1dV[3][3][dim] = {0}, dFi1dGradP[3][3][3] = {0}, dFi1dX[3][3][3] = {0};
    double dFi2dn[3] = {0}, dFi2dS[3][3] = {0}, dFi2dV[3][3][dim] = {0}, dFi2dGradP[3][3][3] = {0}, dFi2dX[3][3][3] = {0};

    computeDerivativeOperatorsOfForceTransmitted(postFcn, X, V, 0, gradP, 0,
                                                 dFi0dV, dFi0dGradP, dFi0dX, dFi0dn, dFi0dS,
                                                 dFi1dV, dFi1dGradP, dFi1dX, dFi1dn, dFi1dS,
                                                 dFi2dV, dFi2dGradP, dFi2dX, dFi2dn, dFi2dS);

 dFi00 = 0.0;
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi00[l] += dFi0dS[l][j]*dS[j];
     for(int k=0; k<dim; ++k)
       dFi00[l] += dFi0dV[j][l][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi00[l] += dFi0dGradP[j][l][k]*dGradP[k][nodeNum(j)] + dFi0dX[j][l][k]*dX[nodeNum(j)][k];
//       dFi00[l] += dFi0dX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

 dFi11 = 0.0;
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi11[l] += dFi1dS[l][j]*dS[j];
     for(int k=0; k<dim; ++k)
       dFi11[l] += dFi1dV[j][l][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi11[l] += dFi1dGradP[j][l][k]*dGradP[k][nodeNum(j)] + dFi1dX[j][l][k]*dX[nodeNum(j)][k];
//       dFi11[l] += dFi1dX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

 dFi22 = 0.0;
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi22[l] += dFi2dS[l][j]*dS[j];
     for(int k=0; k<dim; ++k)
       dFi22[l] += dFi2dV[j][l][k]*dV[nodeNum(j)][k];
     for(int k=0; k<3; ++k) {
       dFi22[l] += dFi2dGradP[j][l][k]*dGradP[k][nodeNum(j)] + dFi2dX[j][l][k]*dX[nodeNum(j)][k];
//       dFi22[l] += dFi2dX[j][l][k]*dX[nodeNum(j)][k];
     }
   }
 }

  dFvv = 0.0;
  diff0 = dFi0 - dFi00;
  diff1 = dFi1 - dFi11;  
  diff2 = dFi2 - dFi22;
  diffv = dFv - dFvv;
  double diff0norm = diff0.norm();
  double diff1norm = diff1.norm();
  double diff2norm = diff2.norm();
  double diffvnorm = diffv.norm();
  double dFi0norm = dFi0.norm();
  double dFi1norm = dFi1.norm();
  double dFi2norm = dFi2.norm();
  double dFvnorm = dFv.norm();
/*  
  if(dFi0norm != 0) {
    double reldiff = diff0norm/dFi0norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi0 is %e\n", diff0norm/dFi0norm);
  } else {
    if(diff0norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi0 is %e\n", diff0norm);
  }

  if(dFi1norm != 0) {
    double reldiff = diff1norm/dFi1norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi1 is %e\n", diff1norm/dFi1norm);
  } else {
    if(diff1norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi1 is %e\n", diff1norm);
  }

  if(dFi2norm != 0) {
    double reldiff = diff2norm/dFi2norm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFi2 is %e\n", diff2norm/dFi2norm);
  } else {
    if(diff2norm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFi2 is %e\n", diff2norm);
  }

  if(dFvnorm != 0) {
    double reldiff = diffvnorm/dFvnorm;
    if(reldiff > 1.0e-12) fprintf(stderr, " ... rel. diff for dFv is %e\n", diffvnorm/dFvnorm);
  } else {
    if(diffvnorm > 1.0e-12) fprintf(stderr, " ... abs. diff for dFv is %e\n", diffvnorm);
  }
*/

   Vec3D dFtot[3];
   dFtot[0] = dFi0 + third*dFv; 
   dFtot[1] = dFi1 + third*dFv;
   dFtot[2] = dFi2 + third*dFv;
/*
   for (int j=0; j<3; ++j) {
     dF[ nodeNum(j) ][0] += dFtot[j][0];
     dF[ nodeNum(j) ][1] += dFtot[j][1];
     dF[ nodeNum(j) ][2] += dFtot[j][2];
   }
*/
   dF[ nodeNum(0) ][0] += dFi00[0] + third*dFvv[0];
   dF[ nodeNum(0) ][1] += dFi00[1] + third*dFvv[1]; 
   dF[ nodeNum(0) ][2] += dFi00[2] + third*dFvv[2];
   dF[ nodeNum(1) ][0] += dFi11[0] + third*dFvv[0];
   dF[ nodeNum(1) ][1] += dFi11[1] + third*dFvv[1];
   dF[ nodeNum(1) ][2] += dFi11[2] + third*dFvv[2];
   dF[ nodeNum(2) ][0] += dFi22[0] + third*dFvv[0]; 
   dF[ nodeNum(2) ][1] += dFi22[1] + third*dFvv[1];
   dF[ nodeNum(2) ][2] += dFi22[2] + third*dFvv[2];

  }
/*
  diff2 = dF2 - dF;
  double diff2norm = diff2.norm();
  double dF2norm = dF2.norm();
  double dFnorm = dF.norm();
  if(dFnorm != 0) fprintf(stderr, " ... rel. diff = %e\n", diff2norm/dFnorm);
  else fprintf(stderr, " ... abs. diff = %e\n", diff2norm);
*/

}

//------------------------------------------------------------------------------
// Included (YC)
template<int dim>
void FaceTria::computeDerivativeOperatorsOfNodalForce(PostFcn *postFcn, SVec<double,3> &X, 
                          SVec<double,dim> &V, double pin, double* gradP[3], 
                          RectangularSparseMat<double,3,3> &dForcedX,
                          RectangularSparseMat<double,3,3> &dForcedGradP,
                          RectangularSparseMat<double,dim,3> &dForcedV,
                          RectangularSparseMat<double,3,3> &dForcedS
                         )
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
    || code == BC_SLIP_WALL_MOVING || code == BC_POROUS_WALL_MOVING) {

    double dFi0dn[3] = {0}, dFi0dS[3][3] = {0}, dFi0dV[3][3][dim] = {0}, dFi0dGradP[3][3][3] = {0}, dFi0dX[3][3][3] = {0};
    double dFi1dn[3] = {0}, dFi1dS[3][3] = {0}, dFi1dV[3][3][dim] = {0}, dFi1dGradP[3][3][3] = {0}, dFi1dX[3][3][3] = {0};
    double dFi2dn[3] = {0}, dFi2dS[3][3] = {0}, dFi2dV[3][3][dim] = {0}, dFi2dGradP[3][3][3] = {0}, dFi2dX[3][3][3] = {0};

    computeDerivativeOperatorsOfForceTransmitted(postFcn, X, V, 0, gradP, 0,
                                                 dFi0dV, dFi0dGradP, dFi0dX, dFi0dn, dFi0dS,
                                                 dFi1dV, dFi1dGradP, dFi1dX, dFi1dn, dFi1dS,
                                                 dFi2dV, dFi2dGradP, dFi2dX, dFi2dn, dFi2dS);
/*
   dF[ nodeNum(0) ][0] += dFi0[0] + third*dFv[0];
   dF[ nodeNum(0) ][1] += dFi0[1] + third*dFv[1]; 
   dF[ nodeNum(0) ][2] += dFi0[2] + third*dFv[2];
   dF[ nodeNum(1) ][0] += dFi1[0] + third*dFv[0];
   dF[ nodeNum(1) ][1] += dFi1[1] + third*dFv[1];
   dF[ nodeNum(1) ][2] += dFi1[2] + third*dFv[2];
   dF[ nodeNum(2) ][0] += dFi2[0] + third*dFv[0]; 
   dF[ nodeNum(2) ][1] += dFi2[1] + third*dFv[1];
   dF[ nodeNum(2) ][2] += dFi2[2] + third*dFv[2];
*/
   dForcedX.addContrib(nodeNum(0), nodeNum(0), dFi0dX[0][0]);
   dForcedX.addContrib(nodeNum(0), nodeNum(1), dFi0dX[1][0]);
   dForcedX.addContrib(nodeNum(0), nodeNum(2), dFi0dX[2][0]);
   dForcedX.addContrib(nodeNum(1), nodeNum(0), dFi1dX[0][0]);
   dForcedX.addContrib(nodeNum(1), nodeNum(1), dFi1dX[1][0]);
   dForcedX.addContrib(nodeNum(1), nodeNum(2), dFi1dX[2][0]);
   dForcedX.addContrib(nodeNum(2), nodeNum(0), dFi2dX[0][0]);
   dForcedX.addContrib(nodeNum(2), nodeNum(1), dFi2dX[1][0]);
   dForcedX.addContrib(nodeNum(2), nodeNum(2), dFi2dX[2][0]);

   dForcedGradP.addContrib(nodeNum(0), nodeNum(0), dFi0dGradP[0][0]);
   dForcedGradP.addContrib(nodeNum(0), nodeNum(1), dFi0dGradP[1][0]);
   dForcedGradP.addContrib(nodeNum(0), nodeNum(2), dFi0dGradP[2][0]);
   dForcedGradP.addContrib(nodeNum(1), nodeNum(0), dFi1dGradP[0][0]);
   dForcedGradP.addContrib(nodeNum(1), nodeNum(1), dFi1dGradP[1][0]);
   dForcedGradP.addContrib(nodeNum(1), nodeNum(2), dFi1dGradP[2][0]);
   dForcedGradP.addContrib(nodeNum(2), nodeNum(0), dFi2dGradP[0][0]);
   dForcedGradP.addContrib(nodeNum(2), nodeNum(1), dFi2dGradP[1][0]);
   dForcedGradP.addContrib(nodeNum(2), nodeNum(2), dFi2dGradP[2][0]);

   dForcedV.addContrib(nodeNum(0), nodeNum(0), dFi0dV[0][0]);
   dForcedV.addContrib(nodeNum(0), nodeNum(1), dFi0dV[1][0]);
   dForcedV.addContrib(nodeNum(0), nodeNum(2), dFi0dV[2][0]);
   dForcedV.addContrib(nodeNum(1), nodeNum(0), dFi1dV[0][0]);
   dForcedV.addContrib(nodeNum(1), nodeNum(1), dFi1dV[1][0]);
   dForcedV.addContrib(nodeNum(1), nodeNum(2), dFi1dV[2][0]);
   dForcedV.addContrib(nodeNum(2), nodeNum(0), dFi2dV[0][0]);
   dForcedV.addContrib(nodeNum(2), nodeNum(1), dFi2dV[1][0]);
   dForcedV.addContrib(nodeNum(2), nodeNum(2), dFi2dV[2][0]);

   dForcedS.addContrib(nodeNum(0), 0, dFi0dS[0]); 
   dForcedS.addContrib(nodeNum(1), 0, dFi1dS[0]); 
   dForcedS.addContrib(nodeNum(2), 0, dFi2dS[0]); 

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
double FaceTria::computeHeatFluxes(ElemSet& elems,
                                     PostFcn* postFcn, SVec<double,3>& X,
                                     Vec<double>& d2wall, double* Vwall,
                                     SVec<double,dim>& V)
{  double hp = 0; 
  if (code == BC_ISOTHERMAL_WALL_MOVING){   
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

    hp = postFcn->computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet);
  }
   return hp;
}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeNodalHeatFluxRelatedValues(ElemSet& elems,
                                     PostFcn* postFcn, SVec<double,3>& X,
                                     Vec<double>& d2wall, double* Vwall,
                                     SVec<double,dim>& V, Vec<double>& P,
                                     Vec<double>& N, bool includeKappa)
{
  if (code == BC_ISOTHERMAL_WALL_MOVING){
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

    double hp = postFcn->computeHeatFluxRelatedValues(dp1dxj, n, d2w, Vwall, Vface, Vtet, includeKappa);
    double norm = sqrt(n*n);
     
        for (int j=0; j<3; ++j)
        {
          if(hp != 0){
            if(N[ nodeNum(j) ] <0){
                 N[ nodeNum(j) ] = 0;
            }         
             P[ nodeNum(j) ] += hp;
             N[ nodeNum(j) ] += norm;
          }
        }
  }
}

//------------------------------------------------------------------------------


// Included (MB)
template<int dim>
void FaceTria::computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, SVec<double,3>& dX, 
				 Vec<double>& d2wall, double* Vwall, double* dVwall, 
				 SVec<double,dim>& V, SVec<double,dim>& dV, double dS[3], Vec<double>& dP)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING) {
    Vec3D n;
    Vec3D dn;

    computeNormalAndDerivative(X, dX, n, dn);
    Elem& elem = elems[elemNum];

    double dp1dxj[4][3], ddp1dxj[4][3];
    if (postFcn->doesFaceNeedGradientP1Function()) {
      elem.computeGradientP1Function(X, dp1dxj);
      elem.computeDerivativeOfGradientP1Function(X, dX, ddp1dxj);
    }

    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double* Vface[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double* Vtet[4] = {V[elem[0]], V[elem[1]], 
		       V[elem[2]], V[elem[3]]};

    double *dVface[3] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)]};
    double *dVtet[4] = {dV[elem[0]], dV[elem[1]],
	  	        dV[elem[2]], dV[elem[3]]};

    double dhp = third * postFcn->computeDerivativeOfHeatPower(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

    for (int j=0; j<3; ++j)
      dP[ nodeNum(j) ] += dhp;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, 
				     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
				     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
				     double* gradP[3], int hydro,
                                     SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                     Vec<double> *genCF)
{
    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;

// Lift, Drag and Forces needed to be computed based on computeForce alone
    Vec3D fi0,fi1,fi2,fv;
    computeForce(elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv, gradP, hydro);
    Fi += fi0 + fi1 + fi2;
    Fv += fv;

// XML Debug
    Vec3D ff = fi0+fi2+fi1;
// Moments need to be computed based on Transmitted forces, which takes care of spatial distribution of force
    computeForceTransmitted(elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv, gradP, hydro);
    Mi += (dx[0] ^ fi0) + (dx[1] ^ fi1) + (dx[2] ^ fi2); // dont remove paranthesis, they are needed for priority reasons
    Mv += dxm ^ fv;

    if(genCF != 0) {
      for(int i = 0; i < genCF->size(); i++) {
        Vec3D d[3] = { (*mX)[i][nodeNum(0)],  (*mX)[i][nodeNum(1)], (*mX)[i][nodeNum(2)] };
        (*genCF)[i] += d[0]*fi0 + d[1]*fi1 + d[2]*fi2;
      }
    }
}

//------------------------------------------------------------------------------

template<int dim>
void FaceTria::computeForceAndMoment(ExactRiemannSolver<dim>& riemann,
                                     VarFcn *varFcn, Vec<Vec3D> &n, Vec<double> &nVel,
                                     ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                                     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                                     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv,
                                     double* gradP[3], int hydro,
                                     SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                     Vec<double> *genCF)
{
    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;

// Lift, Drag and Forces needed to be computed based on computeForce alone
    Vec3D fi0,fi1,fi2,fv;
    computeForce(riemann, varFcn, n, nVel, elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv, gradP, hydro);
    Fi += fi0 + fi1 + fi2;
    Fv += fv;

// XML Debug
    Vec3D ff = fi0+fi2+fi1;
// Moments need to be computed based on Transmitted forces, which takes care of spatial distribution of force
    computeForceTransmitted(elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv, gradP, hydro);
    Mi += (dx[0] ^ fi0) + (dx[1] ^ fi1) + (dx[2] ^ fi2); // dont remove paranthesis, they are needed for priority reasons
    Mv += dxm ^ fv;

    if(genCF != 0) {
      for(int i = 0; i < genCF->size(); i++) {
        Vec3D d[3] = { (*mX)[i][nodeNum(0)],  (*mX)[i][nodeNum(1)], (*mX)[i][nodeNum(2)] };
        (*genCF)[i] += d[0]*fi0 + d[1]*fi1 + d[2]*fi2;
      }
    }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceTria::computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                             SVec<double,3> &X, SVec<double,3> &dX,
                                             Vec<double> &d2wall, double *Vwall, double *dVwall,
                                             SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                             Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
				             double* gradP[3], double* dGradP[3], int hydro)
{

    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D dxds[3] = {dX[nodeNum(0)], dX[nodeNum(1)], dX[nodeNum(2)]};
    Vec3D dx[3] = {x[0]-x0 , x[1]-x0 , x[2]-x0};
    Vec3D ddx[3] = {dxds[0] , dxds[1] , dxds[2]};
    Vec3D dxm = third*(x[0]+x[1]+x[2]) - x0;
    Vec3D ddxm = third*(dxds[0]+dxds[1]+dxds[2]);

    Vec3D dfi0, dfi1, dfi2, dfv;
    computeDerivativeOfForce(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, 0, dfi0, dfi1, dfi2, dfv, gradP, dGradP, hydro);
    dFi += dfi0 + dfi1 + dfi2;
    dFv += dfv;

    Vec3D fi0,fi1,fi2,fv;
    computeForceTransmitted(elems, postFcn, X, d2wall, Vwall, V, 0, fi0,fi1,fi2, fv, gradP, hydro);
    computeDerivativeOfForceTransmitted(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, 0, dfi0, dfi1, dfi2, dfv, gradP, dGradP, hydro);
    dMi += (ddx[0] ^ fi0) + (dx[0] ^ dfi0) + (ddx[1] ^ fi1) + (dx[1] ^ dfi1) + (ddx[2] ^ fi2) + (dx[2] ^ dfi2);
    dMv += (ddxm ^ fv) + (dxm ^ dfv);

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
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED)
  {
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
				   SVec<double,dim> &V, SVec<double,dim> &R,LevelSetStructure *LSS)
{ 
  // In the case of an embedded simulation, check if the face is actually active
  bool isFaceInactive=true;
  if(LSS) 
    {
      for(int i=0;i<3;++i)  isFaceInactive    = (isFaceInactive && !(LSS->isActive(0,nodeNum(i))));
      if(isFaceInactive) return;
    }

  // UH (07/2012) 
  // This test keeps only some faces
  // with a code in {BC_ADIABATIC_WALL_*, BC_ISOTHERMAL_WALL_*}
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

// Included (MB)
template<int dim>
//inline
void FaceTria::computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			       Vec<double> &d2wall, double *Vwall, double *dVwall,
			       SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  // UH (07/2012) 
  // This test keeps only some faces
  // with a code in {BC_ADIABATIC_WALL_*, BC_ISOTHERMAL_WALL_*}
  if (!fet->doesFaceTermExist(code)) return;

  Vec3D n;
  Vec3D dn;

  computeNormalAndDerivative(X, dX, n, dn);


  if (fet->doesFaceNeedGradientP1Function())
    elems[elemNum].computeDerivativeOfFaceGalerkinTerm(fet, nodeNum(), code, n, dn, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  else {
    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double *v[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};
    double *dv[3] = {dV[nodeNum(0)], dV[nodeNum(1)], dV[nodeNum(2)]};

    double dr[dim];
    fet->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, v, dv, dMach, dr);

    for (int l=0; l<3; ++l)
      for (int k=0; k<dim; ++k)
	dR[ nodeNum(l) ][k] -= third * dr[k];
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

  // UH (07/2012) 
  // This test keeps only some faces
  // with a code in {BC_ADIABATIC_WALL_*, BC_ISOTHERMAL_WALL_*}
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

// Included (MB)
template<int dim>
void FaceTria::computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                        double *Vwall, double *dVwall, SVec<double,dim> &V)
{

  // UH (07/2012) 
  // This test keeps only some faces
  // with a code in {BC_ADIABATIC_WALL_*, BC_ISOTHERMAL_WALL_*}
  if (!fet->doesFaceTermExist(code)) return;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED || 
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);

    double d2w[3] = {d2wall[nodeNum(0)], d2wall[nodeNum(1)], d2wall[nodeNum(2)]};
    double *v[3] = {V[nodeNum(0)], V[nodeNum(1)], V[nodeNum(2)]};

    fet->computeBCsJacobianWallValues(code, n, d2w, Vwall, dVwall, v);
    
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
   || code == BC_SLIP_WALL_MOVING || code == BC_POROUS_WALL_MOVING) {

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
					Vec3D &CFi, Vec3D &CMi, Vec3D &CFv, Vec3D &CMv, double* gradP[3],
                                        VecSet< SVec<double,3> > *mX, Vec<double> *genCF) 
{
  static double third = 1.0/3.0;

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ADIABATIC_WALL_MOVING
   || code == BC_SLIP_WALL_MOVING || code == BC_POROUS_WALL_MOVING) {

    Vec3D x[3] = {X[nodeNum(0)], X[nodeNum(1)], X[nodeNum(2)]};
    Vec3D xcg = third * (x[0] + x[1] + x[2]);
    Vec3D Fi0,Fi1,Fi2,Fv;
    double *vWall = reinterpret_cast<double *>(Vwall.data());

//  Forces needed to be computed based on computeForce alone
    computeForce(elems, postFcn, X, d2wall, vWall, V, &pInfty,  Fi0,Fi1,Fi2,Fv, gradP);
    CFi += Fi0 + Fi1 + Fi2;
    CFv += Fv;

// Moments need to be computed based on Transmitted forces, which takes care of spatial distribution of force
    computeForceTransmitted(elems, postFcn, X, d2wall, vWall, V, &pInfty, Fi0,Fi1,Fi2,Fv, gradP);
    CMi += (x[0] - x0) ^ Fi0 + (x[1] - x0) ^ Fi1 + (x[2] - x0) ^ Fi2; // dont remove paranthesis, they are needed for priority reasons
    CMv += (xcg - x0) ^ Fv;

    if(genCF != 0) {
     for(int i = 0; i < genCF->size(); i++) {
       Vec3D d[3] = { (*mX)[i][nodeNum(0)],  (*mX)[i][nodeNum(1)], (*mX)[i][nodeNum(2)] };
       (*genCF)[i] += d[0]*Fi0 + d[1]*Fi1 + d[2]*Fi2;
     }
    }

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
