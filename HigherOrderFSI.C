/* HigherOrderMultiFluid.C */

#include <iostream>
#include "V6NodeData.h"
#include <Elem.h>

//#include <EdgeGrad.h>


inline
HigherOrderFSI::HigherOrderFSI() {
  
  lastPhaseChangeState = NULL;
  limitExtrap = false;
}

inline
HigherOrderFSI::~HigherOrderFSI() {

}

inline 
void HigherOrderFSI::setLimitedExtrapolation() {
  limitExtrap = true;
}

template<int dim>
void HigherOrderFSI::initialize(int numNodes, ElemSet& e, V6NodeData (*v)[2]) {

  SVec<double,dim>* V = new SVec<double,dim>(numNodes);
  *V = 0.0;

  lastPhaseChangeState = static_cast<void*>(V);

  v6data = v;
  
  elems = &e;
  
}

template <int dim>
bool HigherOrderFSI::hasLastPhaseChangeValue(int nodeId) {

  SVec<double,dim>* V = static_cast<SVec<double,dim>*>(lastPhaseChangeState);
  for (int i = 0; i < dim; ++i) {

    if ((*V)[nodeId][i] != 0.0)
      return true;
  }
  
  return false;
}

template <int dim>
const double* HigherOrderFSI::
getLastPhaseChangeValue(int nodeId) {

  SVec<double,dim>* V = 
    static_cast<SVec<double,dim>*>(lastPhaseChangeState);
  
  return (*V)[nodeId];
}

template <int dim>
void HigherOrderFSI::
setLastPhaseChangeValue(int nodeId,const double* v) {

  SVec<double,dim>* V = 
    static_cast<SVec<double,dim>*>(lastPhaseChangeState);

  //if (v[0] != 0.0)
  //  std::cout << "phase change value = " << v[0]  << std::endl;
  memcpy((*V)[nodeId], v ,sizeof(double)*dim);
}

template <int dim>
void HigherOrderFSI::
setLastPhaseChangeValues(SVec<double,dim>& update,Vec<double>& weight) {

  for (int i = 0; i < update.size(); ++i) {

    if (hasLastPhaseChangeValue<dim>(i) || weight[i] <= 0.0)
      continue;

    double v[dim];
    for (int k = 0; k < dim; ++k)
      v[k] = update[i][k] / weight[i];

    setLastPhaseChangeValue<dim>(i, v);
  }
}

//-------------------------------------------------
template <int dim>
void HigherOrderFSI::estimateR(int l, int vertex, 
			       int i, SVec<double,dim>& V, 
			       NodalGrad<dim>& dVdx, SVec<double,3>& X,
			       Vec<int>& fluidId, 
			       double* r) {

  int idxTet    = v6data[l][vertex].tet;
  int idxFace   = v6data[l][vertex].face;
  double face_r = v6data[l][vertex].r;
  double face_t = v6data[l][vertex].t;

  if( (idxTet<0) || (idxTet>=elems->size()) ) {
    memset(r, 0, sizeof(double)*dim);
    return;
  }

  int myfid = fluidId[i];

  Elem& elem = (*elems)[idxTet];

  for (int k = 0; k < 4; ++k) {
    if (fluidId[elem.nodeNum(k)] != myfid) {
      memset(r, 0, sizeof(double)*dim);
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
    xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])+
                       face_t*(X[n1][k] - X[n2][k]);

  double vec[3];
  for (int k = 0; k < 3; ++k)
    vec[k] = X[i][k] - xf[k];

  double Vf[dim],Vfg[dim][3];
  for (int k = 0; k < dim; ++k)
    Vf[k] = V[n2][k] + face_r*(V[n0][k] - V[n2][k])+
                       face_t*(V[n1][k] - V[n2][k]);

  for (int k = 0; k < dim; ++k) {
    Vfg[k][0] = dVdx.getX()[n2][k] + face_r*(dVdx.getX()[n0][k] - dVdx.getX()[n2][k])+
                                     face_t*(dVdx.getX()[n1][k] - dVdx.getX()[n2][k]);

    Vfg[k][1] = dVdx.getY()[n2][k] + face_r*(dVdx.getY()[n0][k] - dVdx.getY()[n2][k])+
                                     face_t*(dVdx.getY()[n1][k] - dVdx.getY()[n2][k]);

    Vfg[k][2] = dVdx.getZ()[n2][k] + face_r*(dVdx.getZ()[n0][k] - dVdx.getZ()[n2][k])+
                                     face_t*(dVdx.getZ()[n1][k] - dVdx.getZ()[n2][k]);
  }
  
  for (int k = 0; k < dim; ++k) {

    if (fabs(V[i][k]-Vf[k] ) > 1.0e-8)
      r[k] = std::min<double>(2.0, (Vfg[k][0]*vec[0]+
				    Vfg[k][1]*vec[1]+
				    Vfg[k][2]*vec[2])/(V[i][k]-Vf[k]) );
    else
      r[k] = 2.0;
    //else
    //  r = 0.0;
    r[k] = std::max<double>(r[k],0.0);
  }
  
  //r[4] =  1.0;
  
}

//--------------------------------------------------


//-------------------------------------------------
template <int dim>
void HigherOrderFSI::estimateRderivative(int l, int vertex, 
					 int i, SVec<double,dim>& V, 
					 NodalGrad<dim>& dVdx, SVec<double,3>& X,
					 Vec<int>& fluidId, 
					 double* dV, double* dV_g) {

  /* // No longer-used
  int idxTet    = v6data[l][vertex].tet;
  int idxFace   = v6data[l][vertex].face;
  double face_r = v6data[l][vertex].r;
  double face_t = v6data[l][vertex].t;

  if( (idxTet<0) || (idxTet>=elems->size()) ) {
    memset(dV,   0, sizeof(double)*dim);
    memset(dV_g, 0, sizeof(double)*dim);
    return;
  }

  int myfid = fluidId[i];

  Elem& elem = (*elems)[idxTet];

  for (int k = 0; k < 4; ++k) {
    if (fluidId[elem.nodeNum(k)] != myfid) {
      memset(dV,   0, sizeof(double)*dim);
      memset(dV_g, 0, sizeof(double)*dim);
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
    xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])+
                       face_t*(X[n1][k] - X[n2][k]);

  double vec[3];
  for (int k = 0; k < 3; ++k)
    vec[k] = X[i][k] - xf[k];

  double Vf[dim],Vfg[dim][3];
  for (int k = 0; k < dim; ++k)
    Vf[k] = V[n2][k] + face_r*(V[n0][k] - V[n2][k])+
                       face_t*(V[n1][k] - V[n2][k]);

  for (int k = 0; k < dim; ++k) {
    Vfg[k][0] = dVdx.getX()[n2][k] + face_r*(dVdx.getX()[n0][k] - dVdx.getX()[n2][k])+
                                     face_t*(dVdx.getX()[n1][k] - dVdx.getX()[n2][k]);

    Vfg[k][1] = dVdx.getY()[n2][k] + face_r*(dVdx.getY()[n0][k] - dVdx.getY()[n2][k])+
                                     face_t*(dVdx.getY()[n1][k] - dVdx.getY()[n2][k]);

    Vfg[k][2] = dVdx.getZ()[n2][k] + face_r*(dVdx.getZ()[n0][k] - dVdx.getZ()[n2][k])+
                                     face_t*(dVdx.getZ()[n1][k] - dVdx.getZ()[n2][k]);
  }
  
  for (int k = 0; k < dim; ++k) {

    dV[k]   = V[i][k] - Vf[k];
    dV_g[k] = Vfg[k][0]*vec[0] + Vfg[k][1]*vec[1] + Vfg[k][2]*vec[2];    

    if (fabs(V[i][k]-Vf[k] ) > 1.0e-8){

      r[k] = (Vfg[k][0]*vec[0] + Vfg[k][1]*vec[1] + Vfg[k][2]*vec[2])/(V[i][k] - Vf[k]);

      if( r[k] > 2.0 ){
  	 r[k] = 2.0;
	delta = 0.0;
      }
      
    } else {

       r[k] = 2.0
      delta = 0.0;

    }

    if(r[k]<0.0) delta = 0.0;

    dV[k]   *= delta;
    dV_g[k] *= delta;

  }  
  */
}

//--------------------------------------------------


template <int dim>
double HigherOrderFSI::computeAlpha(int nodeId,const double* currentV, 
				    const double* neighborV) {

  if (!hasLastPhaseChangeValue<dim>(nodeId)) {

    //std::cout << "alpha = 0!!" << std::endl;
    return 0.0;
  }

  const double* vlast = getLastPhaseChangeValue<dim>(nodeId);

  double alpha = 1.0;
  for (int k = 0; k < dim; ++k) {

    //std::cout << currentV[k] << " " << vlast[k]<< " " << 
    //  neighborV[k] << std::endl;
    alpha = std::min<double>(alpha, fabs(currentV[k]-vlast[k])/
			     std::max<double>(1e-8,fabs(neighborV[k]-currentV[k])));
  }

  //std::cout << "alpha = " << alpha << std::endl;

  return alpha;
}

//----------------------------------------------------

template <int dim>
void HigherOrderFSI::extrapolateV6(int l, int vertex, int i, 
				   SVec<double,dim>& V, 
				   double* Vsurrogate, const double* W, 
				   SVec<double,3>& X,
				   double alpha, double length,
				   Vec<int>& fluidId, double* beta) {

  int idxTet    = v6data[l][vertex].tet;
  int idxFace   = v6data[l][vertex].face;
  double face_r = v6data[l][vertex].r;
  double face_t = v6data[l][vertex].t;

  if ( (idxTet<0) || (idxTet>=elems->size()) ){
    memcpy(Vsurrogate, W, sizeof(double)*dim);
    return;
  }

  int myfid = fluidId[i];

  Elem& elem = (*elems)[idxTet];
  for (int k = 0; k < 4; ++k) {

    if (fluidId[elem.nodeNum(k)] != myfid) {
      memcpy(Vsurrogate, W, sizeof(double)*dim);
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
    xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])+
                       face_t*(X[n1][k] - X[n2][k]);
		       
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

  double frac = 0.5/(1.0 + alpha_f - alpha);
  for (int k = 0; k < dim; ++k) {
    Vsurrogate[k] = ( (1.0-2.0*alpha)*frac*Vf[k] + (1.0+2.0*alpha_f)*frac*W[k] )*beta[k] + 
                    (1.0-beta[k])*W[k];
  }
}

//----------------------------------------------------

template <int dim>
void HigherOrderFSI::RcnExtrap(int l, int vertex, int i, 
			       double length, double alphaij, 
			       SVec<double,dim>& p,
			       double* ddpij,
			       SVec<double,3>& X, Vec<int>& fluidId, double* beta,
			       double* dVsdV, 
			       double *dVpij, double *dVpji) {

  double lc, alpha, alpha_lim = 0.1;

  if(vertex == 0) {    
    lc =  1.0;
    alpha = alphaij;
  }else if(vertex == 1) {    
    lc = -1.0;
    alpha = 1.0 - alphaij;
  }else{
    fprintf(stderr, "Error vertex must be 0 or 1\n");
    exit(-1);
  }

  int k;

  double tmp[dim];
  for (k = 0; k < dim; ++k)
    tmp[k] = p[i][k] + lc*(1.0 - alphaij)*ddpij[k];

  double tmp2[dim], dtmp[dim*dim];

  for (k = 0; k < dim*dim; ++k) dtmp[k] = dVsdV[k];
  DenseMatrixOp<double, dim, dim*dim>::applyToVector(&dtmp, 0, &tmp, 0, &tmp2, 0);
  
  if (v6data == NULL) {

    for (k = 0; k < dim; ++k){
      dVpij[k] = tmp[k];
      dVpji[k] = p[i][k] + (0.5/max(1.0 - alphaij, alpha_lim)) * (tmp2[k] - p[i][k]);  
    }
  
  } else {

    int idxTet    = v6data[l][vertex].tet;
    int idxFace   = v6data[l][vertex].face;
    double face_r = v6data[l][vertex].r;
    double face_t = v6data[l][vertex].t;

    Elem& elem = (*elems)[idxTet];
 
    if ( (idxTet<0) || (idxTet>=elems->size()) ){
      for (k = 0; k < dim; ++k){     
	dVpij[k] = tmp2[k];
	dVpji[k] = tmp2[k];
      }
      return;
    }

    int myfid = fluidId[i];
    for (int ni = 0; ni < 4; ++ni) {
      if (fluidId[elem.nodeNum(ni)] != myfid) {
	for (k = 0; k < dim; ++k){ 
	  dVpij[k] = tmp2[k];
	  dVpji[k] = tmp2[k];
	}
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
    for (k = 0; k < 3; ++k){
      xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])+
	                 face_t*(X[n1][k] - X[n2][k]);
    }
    
    double vec[3];
    for (k = 0; k < 3; ++k)
      vec[k] = X[i][k] - xf[k];
    
    double alpha_f = sqrt((xf[0]-X[i][0])*(xf[0]-X[i][0])+
			  (xf[1]-X[i][1])*(xf[1]-X[i][1])+
			  (xf[2]-X[i][2])*(xf[2]-X[i][2]))/length;

    double frac = 0.5/(1.0 + alpha_f - alpha);

    double eta1 = (1.0 - 2.0*alpha  )*frac;
    double eta2 = (1.0 + 2.0*alpha_f)*frac;

    double pf[dim];
    for (k = 0; k < dim; ++k){
      pf[k] = p[n2][k] + face_r*(p[n0][k]-p[n2][k])+
	                 face_t*(p[n1][k]-p[n2][k]);
    }

    for (k = 0; k < dim; ++k){
      dVpij[k] = beta[k]*(eta1*pf[k] + eta2*tmp2[k]) + (1.0 - beta[k])*tmp2[k];
      dVpji[k] = tmp2[k];
    }

  }


}

//----------------------------------------------------
template <int dim>
void HigherOrderFSI::derivativeofHOFSI(int l, int vertex, int i, 
				       SVec<double,dim>& V, 
                     		       double*  Vf,    double*  Vs,
				       double* dVf_ds, double* dVs_ds,
				       SVec<double,3>& X,
				       double alphaij, double dalphaij_ds,
				       double length, Vec<int>& fluidId, double* beta,
				       double* dVfluid_ds, double* dVstar_ds){

  int k;
  double alpha, dalpha_ds;
  double alpha_lim = 0.1;

  if(vertex == 0) {    
       alpha = alphaij;
   dalpha_ds = dalphaij_ds;
  } else if(vertex == 1) {
        alpha = 1.0 - alphaij;
    dalpha_ds =     -dalphaij_ds;
  } else { 
    fprintf(stderr, "Error vertex must be 0 or 1\n");
    exit(-1);
  }

  if (v6data == NULL) {

    double alpha_b, dalpha_b;

    alpha_b = max(1.0 - alphaij, alpha_lim);

    if( (1.0 - alphaij) > alpha_lim )
      dalpha_b = -dalphaij_ds;
    else
      dalpha_b = 0.0;
    
    for (k=0; k<dim; ++k) { 

      dVfluid_ds[k] = dVf_ds[k];

      dVstar_ds[k] = -(0.5/(alpha_b*alpha_b))*dalpha_b*(Vs[k] - V[i][k])
                   +  (0.5/(alpha_b))*dVs_ds[k];
      
    }
    
  } else {

    for (k = 0; k < dim; ++k)
      dVstar_ds[k] = dVs_ds[k];

    int idxTet    = v6data[l][vertex].tet;
    int idxFace   = v6data[l][vertex].face;
    double face_r = v6data[l][vertex].r;
    double face_t = v6data[l][vertex].t;

    Elem& elem = (*elems)[idxTet];
 
    if ( (idxTet<0) || (idxTet>=elems->size()) ){
      for (k = 0; k < dim; ++k) 
	dVfluid_ds[k] = dVs_ds[k];
      return;
    }

    int myfid = fluidId[i];
    for (int ni = 0; ni < 4; ++ni) {
      if (fluidId[elem.nodeNum(ni)] != myfid) {
	for (k = 0; k < dim; ++k)
	  dVfluid_ds[k] = dVs_ds[k];
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
    for (k = 0; k < 3; ++k){
      xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])+
	                 face_t*(X[n1][k] - X[n2][k]);
    }
    
    double vec[3];
    for (k = 0; k < 3; ++k)
      vec[k] = X[i][k] - xf[k];

    double alpha_f = sqrt((xf[0]-X[i][0])*(xf[0]-X[i][0])+
			  (xf[1]-X[i][1])*(xf[1]-X[i][1])+
			  (xf[2]-X[i][2])*(xf[2]-X[i][2]))/length;

    double frac = 0.5/(1.0 + alpha_f - alpha);

    double dfrac_ds = 0.5*dalpha_ds
                    /( (1.0 + alpha_f - alpha)*(1.0 + alpha_f - alpha) );

    double  eta1 = (1.0 - 2.0*alpha)*frac;
    double deta1 = -2.0*dalpha_ds*frac + (1.0 - 2.0*alpha)*dfrac_ds;

    double  eta2 = (1.0 + 2.0*alpha_f)*frac;
    double deta2 = (1.0 + 2.0*alpha_f)*dfrac_ds;

    double Vf[dim];
    for(k = 0; k < dim; ++k)
      Vf[k] = V[n2][k] + face_r*(V[n0][k]-V[n2][k])+
                         face_t*(V[n1][k]-V[n2][k]);

    for (k = 0; k < dim; ++k) {
      dVfluid_ds[k] = beta[k]*(deta1*Vf[k] + deta2*Vs[k] + eta2*dVs_ds[k]) 
                    + (1.0 - beta[k])*dVs_ds[k];

    }

  }

}

