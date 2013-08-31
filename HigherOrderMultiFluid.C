/* HigherOrderMultiFluid.C

 */

#include <iostream>

#include "V6NodeData.h"
#include <Elem.h>

//#include <EdgeGrad.h>


inline
HigherOrderMultiFluid::HigherOrderMultiFluid() {
  
  lastPhaseChangeState = NULL;
  limitExtrap = false;
}

inline
HigherOrderMultiFluid::~HigherOrderMultiFluid() {

}

inline 
void HigherOrderMultiFluid::setLimitedExtrapolation() {

  limitExtrap = true;
}

template<int dim>
void HigherOrderMultiFluid::initialize(int numNodes,ElemSet& e, V6NodeData (*v)[2]) {

  SVec<double,dim>* V = new SVec<double,dim>(numNodes);
  *V = 0.0;

  lastPhaseChangeState = static_cast<void*>(V);

  v6data = v;
  
  elems = &e;
  
}

template <int dim>
bool HigherOrderMultiFluid::hasLastPhaseChangeValue(int nodeId) {

  SVec<double,dim>* V = static_cast<SVec<double,dim>*>(lastPhaseChangeState);
  for (int i = 0; i < dim; ++i) {

    if ((*V)[nodeId][i] != 0.0)
      return true;
  }
  
  return false;
}

template <int dim>
const double* HigherOrderMultiFluid::
getLastPhaseChangeValue(int nodeId) {

  SVec<double,dim>* V = 
    static_cast<SVec<double,dim>*>(lastPhaseChangeState);
  
  return (*V)[nodeId];
}

template <int dim>
void HigherOrderMultiFluid::
setLastPhaseChangeValue(int nodeId,const double* v) {

  SVec<double,dim>* V = 
    static_cast<SVec<double,dim>*>(lastPhaseChangeState);

  memcpy((*V)[nodeId], v ,sizeof(double)*dim);
}

template <int dim>
double HigherOrderMultiFluid::
estimateR(int l, int vertex, 
	  int i, SVec<double,dim>& V, 
	  NodalGrad<dim>& dVdx, SVec<double,3>& X,
	  Vec<int>& fluidId) {

  int idxTet = v6data[l][vertex].tet;
  int idxFace = v6data[l][vertex].face;
  double face_r = v6data[l][vertex].r;
  double face_t = v6data[l][vertex].t;

  if ((idxTet<0)||(idxTet>=elems->size()))
    return 0.0;

  int myfid = fluidId[i];

  Elem& elem = (*elems)[idxTet];
  for (int k = 0; k < 4; ++k) {

    if (fluidId[elem.nodeNum(k)] != myfid)
      return 0.0;
  }

  int n0_loc = (*elems)[idxTet].faceDef(idxFace, 0);
  int n1_loc = (*elems)[idxTet].faceDef(idxFace, 1);
  int n2_loc = (*elems)[idxTet].faceDef(idxFace, 2);
  int n0 = (*elems)[idxTet].nodeNum(n0_loc);
  int n1 = (*elems)[idxTet].nodeNum(n1_loc);
  int n2 = (*elems)[idxTet].nodeNum(n2_loc);
  
  double xf[3] = {0,0,0};
  for (int k = 0; k < 3; ++k)
    xf[k] = X[n2][k] + face_r*(X[n0][k] -X[n2][k])+
             face_t*( X[n1][k] - X[n2][k]);

  double vec[3];
  for (int k = 0; k < 3; ++k)
    vec[k] = X[i][k] - xf[k];

  double Vf[dim],Vfg[dim][3];
  for (int k = 0; k < dim; ++k)
    Vf[k] = V[n2][k] + face_r*(V[n0][k]-V[n2][k])+
      face_t*(V[n1][k]-V[n2][k]);

  for (int k = 0; k < dim; ++k) {
    Vfg[k][0] = dVdx.getX()[n2][k] + face_r*(dVdx.getX()[n0][k]-dVdx.getX()[n2][k])+
      face_t*(dVdx.getX()[n1][k]-dVdx.getX()[n2][k]);
    Vfg[k][1] = dVdx.getY()[n2][k] + face_r*(dVdx.getY()[n0][k]-dVdx.getY()[n2][k])+
      face_t*(dVdx.getY()[n1][k]-dVdx.getY()[n2][k]);
    Vfg[k][2] = dVdx.getZ()[n2][k] + face_r*(dVdx.getZ()[n0][k]-dVdx.getZ()[n2][k])+
      face_t*(dVdx.getZ()[n1][k]-dVdx.getZ()[n2][k]);
  }

  double r = 2.0;
  for (int k = 0; k < dim; ++k) {

    if (fabs(V[i][k]-Vf[k] ) > 1.0e-8)
      r = std::min<double>(r, 2.0*(Vfg[k][0]*vec[0]+
			       Vfg[k][1]*vec[1]+
			       Vfg[k][2]*vec[2])/(V[i][k]-Vf[k])-1.0);
    else
      r = 0.0;
  }
  
  return 0.0;
  //  return std::max<double>(r,0.0);
}

template <int dim>
double HigherOrderMultiFluid::
computeAlpha(int nodeId,const double* currentV, 
	     const double* neighborV) {

  if (!hasLastPhaseChangeValue<dim>(nodeId))
    return 0.0;

  const double* vlast = getLastPhaseChangeValue<dim>(nodeId);

  double alpha = 1.0;
  for (int k = 0; k < dim; ++k) {

    //std::cout << currentV[k] << " " << vlast[k]<< " " << 
    //  neighborV[k] << std::endl;
    alpha = std::min<double>(alpha, fabs(currentV[k]-vlast[k])/
			     std::max<double>(1e-8,fabs(neighborV[k]-currentV[k])));
  }

  return alpha;
}

template <int dim>
void HigherOrderMultiFluid::
extrapolateV6(int l, int vertex, 
	      int i, SVec<double,dim>& V, 
	      double* Vsurrogate,const double* W, SVec<double,3>& X,
	      double alpha,double length,
	      Vec<int>& fluidId, double beta) {

  int idxTet = v6data[l][vertex].tet;
  int idxFace = v6data[l][vertex].face;
  double face_r = v6data[l][vertex].r;
  double face_t = v6data[l][vertex].t;

  if ((idxTet<0)||(idxTet>=elems->size())){
    memcpy(Vsurrogate,W, sizeof(double)*dim);
    return;
  }

  int myfid = fluidId[i];

  Elem& elem = (*elems)[idxTet];
  for (int k = 0; k < 4; ++k) {

    if (fluidId[elem.nodeNum(k)] != myfid) {
      memcpy(Vsurrogate,W, sizeof(double)*dim);
      return;
    }
  }

  int n0_loc = (*elems)[idxTet].faceDef(idxFace, 0);
  int n1_loc = (*elems)[idxTet].faceDef(idxFace, 1);
  int n2_loc = (*elems)[idxTet].faceDef(idxFace, 2);
  int n0 = (*elems)[idxTet].nodeNum(n0_loc);
  int n1 = (*elems)[idxTet].nodeNum(n1_loc);
  int n2 = (*elems)[idxTet].nodeNum(n2_loc);
  
  double xf[3] = {0,0,0};
  for (int k = 0; k < 3; ++k)
    xf[k] = X[n2][k] + face_r*(X[n0][k] -X[n2][k])+
             face_t*( X[n1][k] - X[n2][k]);

  double vec[3];
  for (int k = 0; k < 3; ++k)
    vec[k] = X[i][k] - xf[k];

  double Vf[dim],Vfg[dim][3];
  for (int k = 0; k < dim; ++k)
    Vf[k] = V[n2][k] + face_r*(V[n0][k]-V[n2][k])+
      face_t*(V[n1][k]-V[n2][k]);

  double alpha_f = sqrt((xf[0]-X[i][0])*(xf[0]-X[i][0])+
			(xf[1]-X[i][1])*(xf[1]-X[i][1])+
			(xf[2]-X[i][2])*(xf[2]-X[i][2]))/length;

  double frac = 0.5/(1.0+alpha_f-alpha);
  for (int k = 0; k < dim; ++k) {

    Vsurrogate[k] = ((1.0-2.0*alpha)*frac*Vf[k]+(1.0+2.0*alpha_f)*frac*W[k])*beta + 
      (1.0-beta)*W[k];
  }
}
