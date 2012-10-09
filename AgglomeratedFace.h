#pragma once

#include <Vector.h>
#include <FluxFcn.h>
#include <GenMatrix.h>

class AgglomeratedFace {

 public:

  AgglomeratedFace();

  AgglomeratedFace(int node, int code);

  AgglomeratedFace(const AgglomeratedFace&);

  ~AgglomeratedFace();
  
  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
			       SVec<double,dim> &V, 
			       double *Ub, SVec<double,dim> &fluxes);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, 
				       SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A);

  template <int dim>
  void computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm*,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,dim> &fluxes);

  template<int dim,class Scalar,int neq>
  void
  computeJacobianThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm* ns,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,neq*neq> &JacX,
                                        SVec<double,neq*neq> &JacY,
                                        SVec<double,neq*neq> &JacZ,
                                        Vec<double>& ctrlVol,
                                        GenMat<Scalar,neq>& A);


  Vec3D& getNormal() { return normal; }

  int getCode() const { return code; }
  int getNode() const { return node; }

  template <int dim>
  void assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U);

  void setNodeType(int*, int*);

 private:

  int code;

  int node;

  Vec3D normal;
};

class AgglomeratedFaceSet {

 public:

  AgglomeratedFaceSet(int size);

  ~AgglomeratedFaceSet();

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
			       SVec<double,dim> &V, 
			       SVec<double,dim>& Ub, SVec<double,dim> &fluxes);

  template <int dim>
  void computeThinLayerViscousFiniteVolumeTerm(class NavierStokesTerm*,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
			                SVec<double,dim> &fluxes);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn,
				       SVec<double,dim> &V, 
				       SVec<double,dim>& Ub, GenMat<Scalar,neq> &A);

  AgglomeratedFace& operator [] (int i) { return myFaces[i]; }

  int size() const { return numFaces; }

 private:

  AgglomeratedFace* myFaces;

  int numFaces;

};
