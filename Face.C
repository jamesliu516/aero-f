#include <Face.h>

#include <FluxFcn.h>
#include <FemEquationTerm.h>
#include <BcData.h>
#include <BcDef.h>
#include <Elem.h>
#include <ExactRiemannSolver.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <LowMachPrec.h>
#include "LevelSet/LevelSetStructure.h"

#include <cmath>

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
  for (int j=0; j<numNodes(); ++j)
    nodeNum(j) = nodemap[ nodeNum(j) ];
}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U)
{
  int k, j;

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");
  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED)
    for (k=0; k<dim; ++k) {
      for (j=0, U[k] = 0.0; j<numNodes(); ++j) 
	U[k] += Uin[nodeNum(j)][k];
      U[k] /= numNodes();
    }
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k) {
      for (j=0, U[k] = 0.0; j<numNodes(); ++j) 
	U[k] += Uout[nodeNum(j)][k];
      U[k] /= numNodes();
    }
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
  int j, k;

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED)
    for (k=0; k<dim; ++k) {
      for (j=0, Uface[k] = 0.0; j<numNodes(); ++j) 
	Uface[k] += Unode[nodeNum(j)][k];
      Uface[k] /= numNodes();
    }
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

    for (int j=0; j<numNodes(); ++j) {
      Unode[ nodeNum(j) ][0] += S;
      for (int k=1; k<dim2; ++k) 
	Unode[ nodeNum(j) ][k] += S * Uface[dim1-dim2+k];
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim1, int dim2>
void Face::computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, double *Uface, double *dUface, SVec<double,dim2> &dUnode)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {

    Vec3D n;
    Vec3D dn;

    computeNormalAndDerivative(X, dX, n, dn);

    double S = sqrt(n*n);
    double dS = 1.0/(2.0*S) * (dn*n + n*dn);

    for (int j=0; j<numNodes(); ++j) {
      dUnode[ nodeNum(j) ][0] += dS;
      for (int k=1; k<dim2; ++k)
	dUnode[ nodeNum(j) ][k] += dS * Uface[dim1-dim2+k] + S * dUface[dim1-dim2+k];
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Face::computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, double *dUfaceSA, SVec<double,dim> &dUnodeSA)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);

    for (int j=0; j<numNodes(); ++j) {
      dNormSA[ nodeNum(j) ][0] += S;
      for (int k=0; k<dim; ++k) 
	dUnodeSA[ nodeNum(j) ][k] += S * dUfaceSA[k];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(VarFcn *varFcn, Vec<Vec3D> &normals, Vec<double> &normalVel,
			   SVec<double,dim> &V, Vec<double> &dt, 
			   TimeLowMachPrec &tprec, Vec<int> &fluidId)
{
  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * getNormalVel(normalVel);

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

  for (int l=0; l<numNodes(); ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double locbeta = tprec.getBeta(locMach);
    
    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    dt[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();
  }   
}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D> &normals, 
			   Vec<double> &normalVel, SVec<double,3> &X, SVec<double,dim> &V, 
			   Vec<double> &idti, Vec<double> &idtv,
			   TimeLowMachPrec &tprec)
{
  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * getNormalVel(normalVel);

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

  for (int l=0; l<numNodes(); ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ]);
    double un = u * n - ndot;

    // Low-Mach Preconditioner
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    double locbeta = tprec.getBeta(locMach);
    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    idti[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();

    // Adam 2010.06.09 Correction

    //  double oldDt = min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();
    //   double oldDt = 15.0*min(0.5*(coeff1-coeff2),-(fabs(u*n)+a))*S/numNodes();
    //double oldDt = 15.0*min(0.5*(coeff1-coeff2),0.0)*S/numNodes();
    // double newDt = 5.0*min(0.5*(coeff1-coeff2),-(fabs(u*n)+a))*S;
 


    //   double veryNewDt =  -(sqrt(u*u)+a)*S/3.0;
    //idti[ nodeNum(l) ] += veryNewDt;

    // idti[ nodeNum(l) ] += 5.0*min(0.5*(coeff1-coeff2),-(fabs(u*n)+a))*S;///numNodes();
    
    // cout<<"beta: "<<locbeta<<" coeff1: "<<coeff1<<" coeff2: "<<coeff2<<" u.n+a: "<<fabs(u*n)+a<<endl;
    // std::cin.get();

    // Previous Version
    // idti[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();

    // double vis = 0.0;
    // if (fet) vis = fet->computeViscousTimeStep(X[nodeNum(l)],V[nodeNum(l)]);
    // idtv[ nodeNum(l) ] += vis*S*S;
  }
    
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Face::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D>  &normals, Vec<Vec3D>  &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
			   SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			   Vec<double> &dIdti, Vec<double> &dIdtv, double dMach,
                           TimeLowMachPrec &tprec)
{

  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  Vec3D dNormal = getdNormal(dNormals);
  double dS = (normal * dNormal) / sqrt(normal * normal);
  double invS = 1.0 / S;
  double dInvS = - dS / (S*S);

  Vec3D n = invS * normal;
  Vec3D dn = dInvS * normal + invS * dNormal;
  double ndot = invS * getNormalVel(normalVel);
  double dndot = dInvS * getNormalVel(normalVel) + invS *getdNormalVel(dNormalVel);

  for (int l=0; l<numNodes(); ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ]);
    Vec3D du = varFcn->getVelocity(dV[ nodeNum(l) ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ]);
    double da = varFcn->computeDerivativeOfSoundSpeed(V[ nodeNum(l) ], dV[ nodeNum(l) ], dMach);
    double un = u * n - ndot;
    double dun = du * n + u * dn - dndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    double locbeta = tprec.getBeta(locMach);
    double dLocMach = varFcn->computeDerivativeOfMachNumber(V[ nodeNum(l) ], dV[ nodeNum(l) ], dMach);
    double dbeta = tprec.getdBeta(locMach,dLocMach);

    double beta2 = locbeta * locbeta;
    double dbeta2 = 2.0 * locbeta * dbeta;
    double coeff1 = (1.0+beta2)*un;
    double dCoeff1 = dbeta2*un + (1.0+beta2)*dun;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);
    double dCoeff2 = (((1.0-beta2)*un)*((-dbeta2*un) + ((1.0-beta2)*dun)) + (2.0*locbeta*a)*(2.0*dbeta*a+2.0*locbeta*da)) / pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    if (min(0.5*(coeff1-coeff2), 0.0) != 0.0)
      dIdti[ nodeNum(l) ] += 0.5*(dCoeff1-dCoeff2)* S/numNodes() + 0.5*(coeff1-coeff2)* dS/numNodes();
      
    double vis = 0.0;
    double dvis = 0.0;
    if (fet) {
      vis = fet->computeViscousTimeStep(X[nodeNum(l)],V[nodeNum(l)]);
      dvis = fet->computeDerivativeOfViscousTimeStep(X[nodeNum(l)],dX[nodeNum(l)],V[nodeNum(l)],dV[nodeNum(l)],dMach);
    }
    dIdtv[ nodeNum(l) ] += dvis*S*S + vis*2.0*S*dS;

  }
    
}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals, 
				   Vec<double> &normalVel, SVec<double,dim> &V, 
				   double *Ub, SVec<double,dim> &fluxes)
{
  if(fluxFcn[code]){
    double flux[dim];
    for (int l=0; l<numNodes(); ++l) {
      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
			     V[nodeNum(l)], Ub, flux);
      for (int k=0; k<dim; ++k){
        fluxes[ nodeNum(l) ][k] += flux[k];
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                   FluxFcn **fluxFcn, Vec<Vec3D> &normals,
                                   Vec<double> &normalVel, SVec<double,dim> &V,
                                   double *Ub, SVec<double,dim> &fluxes)
{
  if(code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED ||
     code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED      ||
     code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
  // FS Riemann based flux calculation.
    double flux[dim], Wstar[2*dim], Vi[2*dim];
    int k;
    Vec3D wallVel, unitNormal;
    VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

    for (int l=0; l<numNodes(); ++l) {
      k = nodeNum(l);
      unitNormal = getNormal(normals,l)/(getNormal(normals,l).norm());
      wallVel = getNormalVel(normalVel, l)/(getNormal(normals,l).norm())*unitNormal;
      for(int iDim=0; iDim<dim; iDim++)
        Vi[iDim] = Vi[iDim+dim] = V[k][iDim];

      riemann.computeFSIRiemannSolution(Vi, wallVel, -1.0*unitNormal, varFcn, Wstar, k/*dummy variable here*/);

//      fluxFcn[BC_INTERNAL]->compute(0.0, 0.0, getNormal(normals,l), getNormalVel(normalVel,l), V[k], Wstar, flux);
      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), Wstar, Ub, flux);

      for (int i=0; i<dim; ++i)
        fluxes[k][i] += flux[i];
    }
    return;
  }

  if(fluxFcn[code]){
    double flux[dim];
    for (int l=0; l<numNodes(); ++l) {
      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l),
                             V[nodeNum(l)], Ub, flux);
      for (int k=0; k<dim; ++k){
        fluxes[ nodeNum(l) ][k] += flux[k];
      }
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void Face::computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<Vec3D> &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
				      SVec<double,dim> &V, double *Ub,
				      double *dUb, SVec<double,dim> &dFluxes)
{

  if(fluxFcn[code]){
    double flux[dim];
    double dFlux[dim];
    for (int l=0; l<numNodes(); ++l) {
      fluxFcn[code]->computeDerivative(0.0, 0.0, getNormal(normals, l), getdNormal(dNormals, l), getNormalVel(normalVel, l), getdNormalVel(dNormalVel, l), V[nodeNum(l)], Ub, dUb, flux, dFlux);
      for (int k=0; k<dim; ++k){
        dFluxes[ nodeNum(l) ][k] += dFlux[k];
      }
    }
  }

}

//------------------------------------------------------------------------------
                                                                                                      
template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				   Vec<double> &normalVel, SVec<double,dim> &V,
				   double *Ub, Vec<int> &fluidId, 
				   SVec<double,dim> &fluxes, LevelSetStructure *LSS)
{
  Vec3D normal = getNormal(normals);
  double flux[dim];
  bool cracking = LSS ? LSS->withCracking() : false;

  if(fluxFcn[code]){
    for (int l=0; l<numNodes(); ++l) {
      if(cracking) {
        if(LSS->isOccluded(0.0,nodeNum(l))) continue;}
      else {
        if(LSS && !LSS->isActive(0.0, nodeNum(l))) continue;}

      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
                             V[nodeNum(l)], Ub, flux, fluidId[nodeNum(l)]);
      for (int k=0; k<dim; ++k)
        fluxes[ nodeNum(l) ][k] += flux[k];
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
inline
void Face::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				     Vec<double> &normalVel, SVec<double,dim> &V,
				     SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiF)
{
  double Uf = 0.0;
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED  ||
      code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED       ||
      code == BC_SYMMETRY
      ) {
    //at wall either U.n = Uwall.n (Euler) or U = Uwall (Navier-Stokes)
    //and thus the flux is 0.0
  }else{
    for (int l=0; l<numNodes(); ++l) {
      Vec3D normal = getNormal(normals, l);
      Uf   = ( V[nodeNum(l)][1]*normal[0] +
               V[nodeNum(l)][2]*normal[1] +
               V[nodeNum(l)][3]*normal[2] ) -
             getNormalVel(normalVel, l);
      for (int k=0; k<dimLS; k++)
        PhiF[ nodeNum(l) ][k] += Uf*Phi[nodeNum(l)][k];
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals, 
					   Vec<double> &normalVel, SVec<double,dim> &V, 
					   double *Ub, GenMat<Scalar,neq> &A)
{

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<neq*neq; ++k) 
      Aii[k] += jac[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
					   Vec<double> &normalVel, SVec<double,dim> &V,
					   double *Ub, GenMat<Scalar,neq> &A, int* nodeType)
{

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    if(!(code == BC_INLET_MOVING || code == BC_OUTLET_MOVING ||
         code == BC_INLET_FIXED  || code == BC_OUTLET_FIXED)) {
      Vec3D normal = getNormal(normals, l);
      double normVel= getNormalVel(normalVel, l);

      fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += jac[k];

    }
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
					   Vec<double> &normalVel, SVec<double,dim> &V, 
					   double *Ub, GenMat<Scalar,neq> &A, Vec<int> &fluidId,
                                           LevelSetStructure* LSS)
{

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    if(LSS && !LSS->isActive(0.0, nodeNum(l))) continue;
    fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac, fluidId[nodeNum(l)]);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<neq*neq; ++k) 
      Aii[k] += jac[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
                                           Vec<double> &normalVel, SVec<double,dim> &V,
                                           double *Ub, GenMat<Scalar,neq> &A, 
                                           Vec<int> &fluidId, int* nodeType)
{

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    if(!(code == BC_INLET_MOVING || code == BC_OUTLET_MOVING ||
         code == BC_INLET_FIXED  || code == BC_OUTLET_FIXED)) {
      Vec3D normal = getNormal(normals, l);
      double normVel= getNormalVel(normalVel, l);

      fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac, fluidId[nodeNum(l)]);
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += jac[k];

    }
  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int dimLS>
inline
void Face::computeJacobianFiniteVolumeTermLS(Vec<Vec3D> &normals,
					     Vec<double> &normalVel, SVec<double,dim> &V,
					     GenMat<Scalar,dimLS> &A)
{

  double Uf = 0.0;
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED  ||
      code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED       ||
      code == BC_SYMMETRY
      ) {
    //at wall either U.n = Uwall.n (Euler) or U = Uwall (Navier-Stokes)
    //and thus the flux is 0.0
  }else{
    for (int l=0; l<numNodes(); ++l) {
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      Vec3D normal = getNormal(normals, l);
      Uf   = ( V[nodeNum(l)][1]*normal[0] +
               V[nodeNum(l)][2]*normal[1] +
               V[nodeNum(l)][3]*normal[2] ) -
             getNormalVel(normalVel, l);
      *Aii /*+= PhiF[ nodeNum(l) ]*/ += Uf;
    }

  }
  /*  double jac;
  for (int l=0; l<numNodes(); ++l) {
    Vec3D normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    jac  = V[nodeNum(l)][1]*normal[0]  + V[nodeNum(l)][2]*normal[1]  + V[nodeNum(l)][3]*normal[2];
    jac *= V[nodeNum(l)][0]; // why did sriram write this? not true!!!
    exit(1);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<dimLS*dimLS; ++k)
      Aii[k] += jac;
      }*/
}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,dim> &V, Vec<double> &dt,
			      TimeLowMachPrec &tprec,
                              Vec<int> &fluidId)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeTimeStep(varFcn, n, ndot, V, dt, tprec, fluidId);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti,
			      Vec<double> &idtv, TimeLowMachPrec &tprec)
{
  
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeTimeStep(fet, varFcn, n, ndot, X, V, idti, idtv, tprec);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			      Vec<double> &dIdti,Vec<double> &dIdtv, double dMach, 
                              TimeLowMachPrec &tprec)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<Vec3D> &dn = geoState.getdFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  Vec<double> &dndot = geoState.getdFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfTimeStep(fet, varFcn, n, dn, ndot, dndot, X, dX, V, dV, dIdti, dIdtv, dMach, tprec);

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

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<numSampledFaces; ++iFace) {
			i = facesConnectedToSampleNode[iFace];
			faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes);
		}
	}
	else {
		for (int i=0; i<numFaces; ++i) 
			faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes);
	}
}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, 
                                      FluxFcn **fluxFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,dim> &V,
                                      SVec<double,dim> &fluxes)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeFiniteVolumeTerm(riemann, fluxFcn, n, ndot, V, Ub[i], fluxes);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V,
				      SVec<double,dim> &dFluxes)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<Vec3D> &dn = geoState.getdFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  Vec<double> &dndot = geoState.getdFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  SVec<double,dim> &dUb = bcData.getdFaceStateVector();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfFiniteVolumeTerm(fluxFcn, n, dn, ndot, dndot, V, Ub[i], dUb[i], dFluxes);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V,
				      Vec<int> &fluidId, SVec<double,dim> &fluxes,
                                      LevelSetStructure *LSS)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                                                                                                         
  for (int i=0; i<numFaces; ++i)  {
    faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluidId, fluxes, LSS);
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void FaceSet::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, BcData<dim> &bcData,
					GeoState &geoState, SVec<double,dim> &V,
					SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiF)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i)  {
    faces[i]->computeFiniteVolumeTermLS(fluxFcn, n, ndot, V, Phi, PhiF);
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
    faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A);
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
    faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, nodeType);
                                                                                                  
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
					      GeoState &geoState, SVec<double,dim> &V, 
					      GenMat<Scalar,neq> &A, Vec<int> &fluidId,
                                              LevelSetStructure* LSS)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId,LSS);

}
//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                              GeoState &geoState, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<int> &fluidId, 
                                              int* nodeType)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                                                                                                  
  for (int i=0; i<numFaces; ++i)
    faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId, nodeType);
                                                                                                  
}


//------------------------------------------------------------------------------

template<int dim, class Scalar, int dimLS>
void FaceSet::computeJacobianFiniteVolumeTermLS(GeoState &geoState, SVec<double,dim> &V,
						GenMat<Scalar,dimLS> &A)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  for (int i=0; i<numFaces; ++i)
    faces[i]->computeJacobianFiniteVolumeTermLS(n, ndot, V, A);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				  BcData<dim> &bcData, GeoState &geoState, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &R, LevelSetStructure *LSS)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<numSampledFaces; ++iFace) {
			i = facesConnectedToSampleNode[iFace];
			faces[i]->computeGalerkinTerm(elems, fet, X, d2wall, Vwall[i], V, R);
		}
	}
	else {
		for (int i=0; i<numFaces; ++i) 
    faces[i]->computeGalerkinTerm(elems, fet, X, d2wall, Vwall[i], V, R, LSS);
	}

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dMach, dR);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
					  BcData<dim> &bcData, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol, 
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall[i], V, A);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, 
					  BcData<dim> &bcData, GeoState &geoState, 
					  SVec<double,3> &X, SVec<double,dim> &V)
{

  SVec<double,dim> &dVwallface = bcData.getdFaceStateVectorSA();

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall[i], dVwallface[i], V);

}

//------------------------------------------------------------------------------
