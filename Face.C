#include <Face.h>

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

template<int dim>
void Face::computeTimeStep(VarFcn *varFcn, Vec<Vec3D> &normals, Vec<double> &normalVel,
			   SVec<double,dim> &V, Vec<double> &dt, double beta,
			   double k1, double cmach)
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
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    beta = fmin(fmax(k1*locMach, beta), cmach);
    
    double beta2 = beta * beta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

    dt[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();
  }   
}
//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D> &normals, 
			   Vec<double> &normalVel, SVec<double,3> &X, SVec<double,dim> &V, 
			   Vec<double> &idti, Vec<double> &idtv, double beta,
			   double k1, double cmach)
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
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    beta = fmin(fmax(k1*locMach, beta), cmach);
    
    double beta2 = beta * beta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*beta*a,2.0),0.5);

    idti[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();

    double vis = 0.0;
    if (fet) vis = fet->computeViscousTimeStep(X[nodeNum(l)],V[nodeNum(l)]);
    idtv[ nodeNum(l) ] += vis*S*S;
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
      fluxFcn[code]->compute(0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
			     V[nodeNum(l)], Ub, flux);
      for (int k=0; k<dim; ++k){
        fluxes[ nodeNum(l) ][k] += flux[k];
      }
    }
  }
}

//------------------------------------------------------------------------------
                                                                                                      
template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				   Vec<double> &normalVel, SVec<double,dim> &V,
				   double *Ub, Vec<double> &Phi, 
				   SVec<double,dim> &fluxes)
{
  Vec3D normal = getNormal(normals);

  double flux[dim];

  if(fluxFcn[code]){
    for (int l=0; l<numNodes(); ++l) {
      if (Phi[nodeNum(l)] >= 0.0) 
        fluxFcn[code]->compute(0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
			       V[nodeNum(l)], Ub, flux, 1);
      if (Phi[nodeNum(l)] <  0.0) 
        fluxFcn[code]->compute(0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
			       V[nodeNum(l)], Ub, flux, -1);
      for (int k=0; k<dim; ++k)
        fluxes[ nodeNum(l) ][k] += flux[k];
    }
  }
                                                                                                      
}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				     Vec<double> &normalsVel, SVec<double,dim> &V,
				     Vec<double> &Phi, Vec<double> &PhiF)
{
  double PhiF1, Uf;

  for (int l=0; l<numNodes(); ++l) {
    Vec3D normal = getNormal(normals, l);    
    
    Uf   = ( V[nodeNum(l)][1]*normal[0] +
	     V[nodeNum(l)][2]*normal[1] +
	     V[nodeNum(l)][3]*normal[2] );
    // Un   = 0.5*(Uf  +fabs(Uf));
    // fluxFcn[code]->computeLS(normal, normalVel, V[nodeNum(l)], Phi[nodeNum(l)], flux);
    
    PhiF1  = Uf*Phi[nodeNum(l)];
    for (int k=0; k<1; ++k)
      // PhiF[ nodeNum(l) ] += flux;
      PhiF[ nodeNum(l) ] += PhiF1;

    if (isnan(Uf)) fprintf(stderr, "in face computeFiniteVolumeTerm for Uf %d \n",nodeNum(l));
    if (isnan(PhiF1)) fprintf(stderr, "in face computeFiniteVolumeTerm %d \n",nodeNum(l));
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

    fluxFcn[code]->computeJacobian(0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
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

      fluxFcn[code]->computeJacobian(0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += jac[k];

    }
  }
}

//------------------------------------------------------------------------------
template<int dim, class Scalar>
inline
void Face::computeJacobianFiniteVolumeTermLS(Vec<Vec3D> &normals,
					     Vec<double> &normalVel, SVec<double,dim> &V,
					     GenMat<Scalar,1> &A)
{

  double jac;
  for (int l=0; l<numNodes(); ++l) {
    Vec3D normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    //    fluxFcn[code]->computeJacobian(normal, normalVel, V[nodeNum(l)], Ub, jac);
    jac  = V[nodeNum(l)][1]*normal[0]  + V[nodeNum(l)][2]*normal[1]  + V[nodeNum(l)][3]*normal[2];
    jac *= V[nodeNum(l)][0];
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<1; ++k)
      Aii[k] += jac;
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
    faces[i]->computeTimeStep(varFcn, n, ndot, V, dt, beta, k1, cmach);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti,
			      Vec<double> &idtv, double beta, double k1, double cmach)
{
  
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeTimeStep(fet, varFcn, n, ndot, X, V, idti, idtv, beta, k1, cmach);

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

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes);

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
    faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], Phi, fluxes);
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

template<int dim, class Scalar>
void FaceSet::computeJacobianFiniteVolumeTermLS(GeoState &geoState, SVec<double,dim> &V,
						GenMat<Scalar,1> &A)
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
				  SVec<double,dim> &R)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeGalerkinTerm(elems, fet, X, d2wall, Vwall[i], V, R);

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

