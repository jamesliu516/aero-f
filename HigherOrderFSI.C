/* HigherOrderMultiFluid.C */

#include <iostream>
#include "V6NodeData.h"
#include <Elem.h>
#include <DenseMatrixOps.h>

inline
HigherOrderFSI::HigherOrderFSI(const IoData & iod)
{
  locationhalfriemannproblem = iod.embed.locationhalfriemannproblem;
  geomTol = 1.0e-10;
  lastPhaseChangeState = NULL;
  SIData_p = NULL;
  SIData_m = NULL;
  FEMData_p = NULL;
  FEMData_m = NULL;
  limitExtrap = false;
}

inline
HigherOrderFSI::~HigherOrderFSI() {
    std::cout <<" HigherOrderFSI desctructor " << std::endl;
	if(SIData_p)  delete[] SIData_p;
    if(SIData_m)  delete[] SIData_m;
    if(FEMData_p) delete[] FEMData_p;
    if(FEMData_m) delete[] FEMData_m;
    //if(lastPhaseChangeState) delete lastPhaseChangeState;
}

inline 
void HigherOrderFSI::setLimitedExtrapolation() 
{
  limitExtrap = true;
}

template<int dim>
void HigherOrderFSI::initialize(int numNodes, ElemSet& e, V6NodeData (*v)[2]) 
{

  SVec<double,dim>* V = new SVec<double,dim>(numNodes);
  *V = 0.0;

  lastPhaseChangeState = static_cast<void*>(V);

  v6data = v;
  
  elems = &e;
  
}

template <int dim>
void HigherOrderFSI::initialize(IoData &iod, Communicator *comm, int numNodes, ElemSet& e,
										  bool linearReconstruction, bool hoViscReconstruction)
{

	SVec<double,dim>* V = new SVec<double,dim>(numNodes);
	*V = 0.0;

	lastPhaseChangeState = static_cast<void*>(V);

	elems = &e;

	HOtreatment = linearReconstruction;

	viscQuadRcn = hoViscReconstruction;

}

template <int dim>
bool HigherOrderFSI::hasLastPhaseChangeValue(int nodeId) 
{

  SVec<double,dim>* V = static_cast<SVec<double,dim>*>(lastPhaseChangeState);

  for(int i = 0; i < dim; ++i) 
  {
    if((*V)[nodeId][i] != 0.0) return true;
  }
  
  return false;
}

template <int dim>
const double* HigherOrderFSI::getLastPhaseChangeValue(int nodeId) 
{

  SVec<double,dim>* V = static_cast<SVec<double,dim>*>(lastPhaseChangeState);
  
  return (*V)[nodeId];
}

template <int dim>
void HigherOrderFSI::setLastPhaseChangeValue(int nodeId, const double* v) 
{

  SVec<double,dim>* V = static_cast<SVec<double,dim>*>(lastPhaseChangeState);

  //if (v[0] != 0.0)
  //  std::cout << "phase change value = " << v[0]  << std::endl;
  memcpy((*V)[nodeId], v ,sizeof(double)*dim);
}

template <int dim>
void HigherOrderFSI::setLastPhaseChangeValues(SVec<double,dim>& update,Vec<double>& weight) 
{

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
double HigherOrderFSI::computeAlpha(int nodeId, const double* currentV, const double* neighborV) 
{

  if(!hasLastPhaseChangeValue<dim>(nodeId)) 
  {
    //std::cout << "alpha = 0!!" << std::endl;
    return 0.0;
  }

  const double* vlast = getLastPhaseChangeValue<dim>(nodeId);

  double alpha = 1.0;

  for(int k = 0; k < dim; ++k) 
  {
    //std::cout << currentV[k] << " " << vlast[k]<< " " << neighborV[k] << std::endl;
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

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
//inline
//void HigherOrderFSI::setSIstencil(V6NodeData *SIstencilData)
//{
//	SIData = SIstencilData;
//}
//
////------------------------------------------------------------------------------
//inline
//void HigherOrderFSI::setFEMstencil(V6NodeData *FEMstencilData_p,
//											  V6NodeData *FEMstencilData_m)
//{
//	FEMData_p = FEMstencilData_p;
//	FEMData_m = FEMstencilData_m;
//}

//------------------------------------------------------------------------------

template <int dim>
void HigherOrderFSI::extrapolateToWall_1(int l, int n, int Fid, VarFcn *varFun, 
													  SVec<double,dim>& V, NodalGrad<dim>& dV, double* V_n, 
													  SVec<double,3>& X, Vec3D &xWall, Vec3D &normWall,Vec3D &Xij,
													  double* V_ext, bool externalSI)
{	

	//	Extrapolate the solution from Xf to xWall using Vf and dVf 

	for(int k=0; k<dim; ++k) 
	{
		V_ext[k]     = V_n[k];
		V_ext[k+dim] = V_n[k];
	}

	if(!HOtreatment) return;

	Vec3D Xf;

	double Vf[dim];

	double dVf[dim][3];

    int idxTet = -1, idxFace =  -1, face_r = 0.0, face_t = 0.0;
    bool pOrNot;
    Vec3D dir1 = Xij - xWall, dir2(X[n][0] - xWall[0], X[n][1] - xWall[1], X[n][2] - xWall[2]);
    double norm = sqrt(dir1*dir1);
    if(norm > geomTol && externalSI){
        pOrNot = (dir1*normWall >= 0);
    }else {
        pOrNot = (dir2 * normWall >= 0);
    }

//    if(!pOrNot){std::cout << "****H interpolateToWall_1 pOrNot, node "<< n << " l " << l << " norm " << norm
//                          << "\n x[n] " << X[n][0] << " " << X[n][1] << " "  << X[n][2]
//                          << "\n dir2 "<< dir2[0] << " " << dir2[1] << " " << dir2[2]
//                          << "\n dir1 "<< dir1[0] << " " << dir1[1] << " " << dir1[2]
//                          << "\n norm "<< normWall[0] << " " << normWall[1] << " " << normWall[2]
//                          << "\n xWall"<< xWall[0] << " " << xWall[1] << " " << xWall[2]<< std::endl;}

    if(pOrNot) {
        idxTet = SIData_p[l].tet;
        idxFace = SIData_p[l].face;
        face_r = SIData_p[l].r;
        face_t = SIData_p[l].t;
    }else{
        idxTet = SIData_m[l].tet;
        idxFace = SIData_m[l].face;
        face_r = SIData_m[l].r;
        face_t = SIData_m[l].t;
    }

	if(idxTet < 0) 
	{
        //std::cout << "*** Warning: edge "<< l << " does not have stencil to compute second order Euler flux based on the closest point " << std::endl;
		return;

/*
		for(int k=0; k<3; ++k) Xf[k] = X[n][k];

		for(int k=0; k<dim; ++k) Vf[k] = V[n][k];
		
		for(int k=0; k<dim; ++k)
		{		
		 	dVf[k][0] = dV.getX()[n][k];
			dVf[k][1] = dV.getY()[n][k];
			dVf[k][2] = dV.getZ()[n][k];
		}
*/
	}
	else
	{
		int n0_loc = (*elems)[idxTet].faceDef(idxFace, 0);
		int n1_loc = (*elems)[idxTet].faceDef(idxFace, 1);
		int n2_loc = (*elems)[idxTet].faceDef(idxFace, 2);

		int n0 = (*elems)[idxTet].nodeNum(n0_loc);
		int n1 = (*elems)[idxTet].nodeNum(n1_loc);
		int n2 = (*elems)[idxTet].nodeNum(n2_loc);

		for(int k=0; k<3; ++k)
			Xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])
				              + face_t*(X[n1][k] - X[n2][k]);
		
		for(int k=0; k<dim; ++k)
			Vf[k] = V[n2][k] + face_r*(V[n0][k] - V[n2][k])
				              + face_t*(V[n1][k] - V[n2][k]);

		for(int k=0; k<dim; ++k) 
		{		
			dVf[k][0] = dV.getX()[n2][k] 
		      		 + face_r*(dV.getX()[n0][k] - dV.getX()[n2][k])
				       + face_t*(dV.getX()[n1][k] - dV.getX()[n2][k]);
			
			dVf[k][1] = dV.getY()[n2][k] 
				       + face_r*(dV.getY()[n0][k] - dV.getY()[n2][k])
				       + face_t*(dV.getY()[n1][k] - dV.getY()[n2][k]);
			
			dVf[k][2] = dV.getZ()[n2][k] 
				       + face_r*(dV.getZ()[n0][k] - dV.getZ()[n2][k])
				       + face_t*(dV.getZ()[n1][k] - dV.getZ()[n2][k]);
		}

	}	

	Vec3D Dir = xWall - Xf;

	double DV[dim];
	for(int k=0; k<dim; ++k) DV[k] = dVf[k][0]*Dir[0] + dVf[k][1]*Dir[1] + dVf[k][2]*Dir[2];

	for(int k=0; k<dim; ++k) {
		V_ext[k]     = Vf[k] + DV[k];
		V_ext[k+dim] = Vf[k] + DV[k];
	}


	if(V_ext[0] <= 0  || V_ext[4] <= 0 || face_r < geomTol || face_t < geomTol || (1 - face_r - face_t) < geomTol)
	{   //Negative pressure, negative density or the stencil has ghost nodes
        for(int k=0; k<dim; ++k)
        {
            V_ext[k]     = V_n[k];
            V_ext[k+dim] = V_n[k];
        }
	}

}

//------------------------------------------------------------------------------

template <int dim>
void HigherOrderFSI::interpolateToSI(int l, int n, int Fid, VarFcn *varFun, 
												 SVec<double,dim>& V, double* Vstar, NodalGrad<dim>& dV,
												 SVec<double,3>& X, Vec3D &xWall, Vec3D &normWall, Vec3D &Xij,
												 double* Vsi, double limiter, bool externalSI)
{

	if(Vstar[0] < 0 || Vstar[4] < 0)
	{
		fprintf(stderr, "*: negative density/pressure at Xw=[%f,%f,%f]\n", xWall[0], xWall[1], xWall[2]);
		exit(-1);
	}

	for(int k=0; k<2*dim; ++k) Vsi[k] = Vstar[k];

	if(!HOtreatment) return;

    int idxTet = -1, idxFace =  -1, face_r = 0.0, face_t = 0.0;
    bool pOrNot;
    Vec3D dir1 = Xij - xWall, dir2(X[n][0] - xWall[0], X[n][1] - xWall[1], X[n][2] - xWall[2]);
    double norm = sqrt(dir1*dir1);
    if(norm > geomTol && externalSI){
       pOrNot = dir1*normWall >= 0;
    }else {
       pOrNot = dir2 * normWall >= 0;
    }
    if(pOrNot) {
        idxTet = SIData_p[l].tet;
        idxFace = SIData_p[l].face;
        face_r = SIData_p[l].r;
        face_t = SIData_p[l].t;
    }else{
        idxTet = SIData_m[l].tet;
        idxFace = SIData_m[l].face;
        face_r = SIData_m[l].r;
        face_t = SIData_m[l].t;
    }



	if(idxTet < 0)	return;

	int n0_loc = (*elems)[idxTet].faceDef(idxFace, 0);
	int n1_loc = (*elems)[idxTet].faceDef(idxFace, 1);
	int n2_loc = (*elems)[idxTet].faceDef(idxFace, 2);

	int n0 = (*elems)[idxTet].nodeNum(n0_loc);
	int n1 = (*elems)[idxTet].nodeNum(n1_loc);
	int n2 = (*elems)[idxTet].nodeNum(n2_loc);

	Vec3D Xf;	
	for(int k=0; k<3; ++k)
		Xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])
			              + face_t*(X[n1][k] - X[n2][k]);

	double Vf[dim]; 		
	for(int k=0; k<dim; ++k)
		Vf[k] = V[n2][k] + face_r*(V[n0][k] - V[n2][k])
		    	           + face_t*(V[n1][k] - V[n2][k]);		

	double dVf[dim][3];
	for(int k=0; k<dim; ++k) 
	{		
		dVf[k][0] = dV.getX()[n2][k] 
		       	 + face_r*(dV.getX()[n0][k] - dV.getX()[n2][k])
			       + face_t*(dV.getX()[n1][k] - dV.getX()[n2][k]);
		
		dVf[k][1] = dV.getY()[n2][k] 
			       + face_r*(dV.getY()[n0][k] - dV.getY()[n2][k])
				    + face_t*(dV.getY()[n1][k] - dV.getY()[n2][k]);
		
		dVf[k][2] = dV.getZ()[n2][k] 
			       + face_r*(dV.getZ()[n0][k] - dV.getZ()[n2][k])
			       + face_t*(dV.getZ()[n1][k] - dV.getZ()[n2][k]);
	}

  	//Vec3D Dir_wSI = xWall - Xij;
	Vec3D Dir_wf  = xWall - Xf;
    Vec3D dx = 2.0*(Xij -Xf);

	//double eta = sqrt(Dir_wSI * Dir_wSI) / sqrt(Dir_wf  * Dir_wf);

    //alpha is XjxWall/XjXi
    double alpha = sqrt(Dir_wf  * Dir_wf)/sqrt(dx * dx);
	// Linear interpolation
	//double Vsi_Test[dim];
    //for(int k=0; k<dim; ++k) Vsi[k] = (1.0 - eta)*Vstar[k] + eta*Vf[k];
    for(int k=0; k<dim; ++k) Vsi[k] = Vf[k] + 0.5/max(alpha, limiter) *(Vstar[k] - Vf[k]);




	 if(Vsi[0] <= 0 || Vsi[4] <= 0)//failsafe
	 {

         for(int k=0; k<dim; ++k)
         {
             Vsi[k]     = Vstar[k];
             Vsi[k+dim] = Vstar[k];
         }
	 }



}

//------------------------------------------------------------------------------

template <int dim>
bool HigherOrderFSI::setFEGhostPoint(int dir, int i, VarFcn *varFun, SVec<double,dim>& U, 
												 NodalGrad<dim>& dV, SVec<double,3>& X, Vec<int> &fluidId,
												 Vec3D &xWall, Vec3D &vWall, Vec3D &nWall, 
												 bool isIsoTherm, double TWall, 
												 FemEquationTerm *fet, double* Vg, int &fId)

{

	// NOTE: U is the vector of the conservative variables
   //       dV is the nodal gradient of the primitive variables
	
	int idxTet, idxFace;
	double face_r, face_t;

	if(dir > 0)
	{
		idxTet  = FEMData_p[i].tet; // stencil tet
		idxFace = FEMData_p[i].face; // stencil face
		face_r  = FEMData_p[i].r; //stencil face coordinate1
		face_t  = FEMData_p[i].t; //stencil face coordinate2
	}
	else
	{
		idxTet  = FEMData_m[i].tet;
		idxFace = FEMData_m[i].face;
		face_r  = FEMData_m[i].r;
		face_t  = FEMData_m[i].t;
	}

	if(idxTet < 0)	return false;

	int n0_loc = (*elems)[idxTet].faceDef(idxFace, 0);
	int n1_loc = (*elems)[idxTet].faceDef(idxFace, 1);
	int n2_loc = (*elems)[idxTet].faceDef(idxFace, 2);

	int n0 = (*elems)[idxTet].nodeNum(n0_loc);
	int n1 = (*elems)[idxTet].nodeNum(n1_loc);
	int n2 = (*elems)[idxTet].nodeNum(n2_loc);

	double V0[dim], V1[dim], V2[dim];
	varFun->conservativeToPrimitive(U[n0], V0, fluidId[n0]);
	varFun->conservativeToPrimitive(U[n1], V1, fluidId[n1]);
	varFun->conservativeToPrimitive(U[n2], V2, fluidId[n2]);
    //stencil fluid id
	fId = (int) (fluidId[n2] + face_r*(fluidId[n0] - fluidId[n2])
					             + face_t*(fluidId[n1] - fluidId[n2]));

	Vec3D Xf;
	for(int k=0; k<3; ++k)//stencil point
		Xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])
                       + face_t*(X[n1][k] - X[n2][k]);

	double* Vf = new double[dim];

	for (int k=0; k<dim; ++k)//stencil primitive variable
		Vf[k] = V2[k] + face_r*(V0[k] - V2[k])
                    + face_t*(V1[k] - V2[k]);

	double** dVf = new double* [dim];
	for(int k=0; k<dim; ++k) dVf[k] = new double [3];

	for(int k=0; k<dim; ++k) //stencil primitive gradient
	{		
		dVf[k][0] = dV.getX()[n2][k] 
			       + face_r*(dV.getX()[n0][k] - dV.getX()[n2][k])
			       + face_t*(dV.getX()[n1][k] - dV.getX()[n2][k]);
		
		dVf[k][1] = dV.getY()[n2][k] 
			       + face_r*(dV.getY()[n0][k] - dV.getY()[n2][k])
			       + face_t*(dV.getY()[n1][k] - dV.getY()[n2][k]);

		dVf[k][2] = dV.getZ()[n2][k] 
			       + face_r*(dV.getZ()[n0][k] - dV.getZ()[n2][k])
			       + face_t*(dV.getZ()[n1][k] - dV.getZ()[n2][k]);
	}

	// Replace the gradient of the pressure with 
	// the gradient of the temperature (5th variable)...
	for(int j=0; j<3; ++j)
	{
		double dV_[dim];
		for(int k=0; k<dim; ++k) dV_[k] = dVf[k][j];

		double dT = varFun->computeDerivativeOfTemperature(Vf, dV_);//todo check

		dVf[4][j] = dT;
	}

	// ...and replace the pressure variable with the temperature (5th variable)
	// Do not change the order of operations for dVf[4] and Vf[4]!!!
	double T;
	T = varFun->computeTemperature(Vf, fluidId[i]);
	Vf[4] = T;

	Vec3D Dir = xWall - Xf;
	double xi = sqrt(Dir*Dir); //wall distance
	if(xi != 0) Dir *= 1.0 / xi; //unit vector from Xf to xWall

	Vec3D d;
	for(int k=0; k<3; ++k) d[k] = X[i][k] - Xf[k];
	double eta = sqrt(d*d); // distance from ghost node to stencil

	double coeff_1, coeff_2, coeff_3;
	 
	//	Density: dummy value (used only to compute P from T
	//	and viceversa; P isn't relevant for viscous flux)	
	Vg[0] = Vf[0];

	Vec3D nW = nWall/nWall.norm(), tgW1, tgW2;
	double dudn, dTdn;
	bool wWF;
	if(fet->withWallFcn()) 
	{
		computeWallVersors(Vf, nW, varFun, tgW1, tgW2);
		
		//fprintf(stdout, "%f,%f,%f ", xWall[0],xWall[1],xWall[2]);
		wWF = computeDuDTwf(Vf, varFun, xi, vWall, TWall, 
								  tgW1, tgW2, fet, dudn, dTdn);
	}

	// Velocity	
	if(fet->withWallFcn() && wWF) /* Wall functions */
		setVelocityWfG(Vf, dVf, varFun, Dir, xi, eta, 
							vWall, nW, tgW1, tgW2, dudn, Vg);
	else
		setVelocityG(Vf, dVf, Dir, xi, eta, vWall, Vg);
	
	// Temperature
	setTemperatureG(Vf, dVf, Dir, xi, eta, isIsoTherm, TWall, Vg);

	// Turbulent stuff
	if(dim > 5) setTurboG(Vf, dVf, Dir, xi, eta, dim, Vg);

	// Replace back the 4th variable with the pressure
	//T = Vg[4];
	//varFun->getV4FromTemperature(Vg, T, fluidId[i]);
	// ...this is done in ghostpoints::reduce

	delete [] Vf;

	for(int k=0; k<dim; ++k) delete [] dVf[k];
	delete [] dVf;

	return true;
}

//------------------------------------------------------------------------------

	/* 
		Quadratic inter/extra-polation of the solution at the ghost node

		V(eta) = coeff_1*eta^2 + coeff_2*eta + coeff_3
		V'(eta) = 2*coeff_1*eta + coeff_2

		eta = 0 @ Xf; eta = 1 @ X_i
		
		For velocity, turbulent quantities and temperature with isothermal wall
		*** V(0) = Vf; V'(0) = dVf(0); V(xi_Wall) = V_wall

		For temperature with adiabatic wall 
		*** T(0) = Tf; V'(0) = dTf(0); T'(xi_wall) = 0
  
      ----- 

		Linear inter/extra-polation of the solution at the ghost node

		V(eta) = coeff_1*eta + coeff_2
		V'(eta) = coeff_1*eta

		eta = 0 @ Xf; eta = 1 @ X_i
		
		For velocity, turbulent quantities and temperature with isothermal wall
		*** V(0) = Vf; V(xi_Wall) = V_wall

		For temperature with adiabatic wall 
		*** T(0) = Tf;  T'(xi_wall) = 0
	*/

inline
void HigherOrderFSI::setVelocityG(double* Vf, double** dVf,
											 Vec3D Dir, double xi, double eta,
	                               Vec3D &vWall, double* Vg)
{

	double coeff_1, coeff_2, coeff_3;
	
	for(int k=1; k<4; ++k)
	{
		double du_dn;
		if(viscQuadRcn)
		{
			coeff_3 = Vf[k];
			coeff_2 = dVf[k][0]*Dir[0] + dVf[k][1]*Dir[1] + dVf[k][2]*Dir[2]; 
  		   coeff_1 = (vWall[k-1] - coeff_2*xi - coeff_3)/(xi*xi);	

			Vg[k] = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		}
		else
		{
			coeff_2 = Vf[k];
			coeff_1 = (vWall[k-1] - coeff_2)/xi;

			Vg[k] = coeff_1*eta + coeff_2;
		}
	}

}

//------------------------------------------------------------------------------

inline
void HigherOrderFSI::setVelocityWfG(double* Vf, double** dVf, VarFcn *vf,
												Vec3D Dir, double xi, double eta,
												Vec3D &vWall, Vec3D &nW, 
												Vec3D &tgW1, Vec3D &tgW2,
												double dudn, double* Vg)
{

	double coeff_1, coeff_2, coeff_3;

	double u1, uw, ug_n, ug_tg2, ug_tg1, dudxi;

	Vec3D uf = vf->getVelocity(Vf);
	
	if(viscQuadRcn)
	{	
		// Normal component 
		u1 =    uf * nW;		
		uw = vWall * nW;
		
		dudxi = 0.0;
		for(int k=0; k<3; ++k) 
			dudxi += (dVf[1][k]*nW[0] + dVf[2][k]*nW[1] + dVf[3][k]*nW[2])*Dir[k];

		coeff_3 = u1;
		coeff_2 = dudxi;
		coeff_1 = (uw - coeff_3 - coeff_2*xi)/(xi*xi);
		
		ug_n = coeff_1*eta*eta + coeff_2*eta + coeff_3;

		// Second tangential component
		u1 =    uf * tgW2;		
		uw = vWall * tgW2;

		dudxi = 0.0;
		for(int k=0; k<3; ++k) 
			dudxi += (dVf[1][k]*tgW2[0] + dVf[2][k]*tgW2[1] + dVf[3][k]*tgW2[2])*Dir[k];

		coeff_3 = u1;
		coeff_2 = dudxi;		
		coeff_1 = (uw - coeff_3 - coeff_2*xi)/(xi*xi);
		
		ug_tg2 = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		
		// First tangential component		
      /*
		u1 = uf * tgW1;
				
		dudxi = 0.0;
		for(int k=0; k<3; ++k) 
			dudxi += (dVf[1][k]*tgW1[0] + dVf[2][k]*tgW1[1] + dVf[3][k]*tgW1[2])*Dir[k];

		coeff_3 = u1;
		coeff_2 = dudxi; 
		coeff_1 = -0.5*(dudn + coeff_2)/xi;	// eta has opposite direction to n_wall

		ug_tg1 = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		*/
		u1 = uf * tgW1;
		
		coeff_2 = u1;
		coeff_1 = -dudn; // eta has opposite direction to n_wall
		
		ug_tg1 = coeff_1*eta + coeff_2;
		
	}
	else
	{
		// Normal component
		u1 =    uf * nW;		
		uw = vWall * nW;

		coeff_2 = u1;
		coeff_1 = (uw - coeff_2)/xi;

		ug_n = coeff_1*eta + coeff_2;
		
		// Second tangential component
		u1 =    uf * tgW2;		
		uw = vWall * tgW2;
		
		coeff_2 = u1;
		coeff_1 = (uw - coeff_2)/xi;

		ug_tg2 = coeff_1*eta + coeff_2;

		// First tangential component
		u1 = uf * tgW1;
		
		coeff_2 = u1;
		coeff_1 = -dudn; // eta has opposite direction to n_wall
		
		ug_tg1 = coeff_1*eta + coeff_2;

		//fprintf(stdout, "%f, %f, %f, %f, %f, %f\n", tgW1[1],xi,eta, u1, dudn, ug_tg1);
	}

	Vec3D ug;
	ug =  ug_n   * nW;
	ug += ug_tg1 * tgW1;
	ug += ug_tg2 * tgW2;

	for(int k=0; k<3; ++k) Vg[k+1] = ug[k];
}

//------------------------------------------------------------------------------

inline
void HigherOrderFSI::setTemperatureG(double *Vf, double** dVf,
												 Vec3D Dir, double xi, double eta,
												 bool isIsoTherm, double TWall, double* Vg)
{

	double coeff_1, coeff_2, coeff_3;
	
	if(isIsoTherm)
	{	
		if(viscQuadRcn)
		{
			coeff_3 = Vf[4];
		   coeff_2 = dVf[4][0]*Dir[0] + dVf[4][1]*Dir[1] + dVf[4][2]*Dir[2]; 
		   coeff_1 = (TWall - coeff_2*xi - coeff_3)/(xi*xi);

		   Vg[4] = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		}
		else
		{
			coeff_2 = Vf[4];
			coeff_1 = (TWall - coeff_2)/xi;

			Vg[4] = coeff_1*eta + coeff_2;
		}
	}
	else /*Adiabatic*/
	{
		if(viscQuadRcn)
		{
			coeff_3 = Vf[4];
		   coeff_2 = dVf[4][0]*Dir[0] + dVf[4][1]*Dir[1] + dVf[4][2]*Dir[2]; 
		   coeff_1 = -0.5*coeff_2/xi;

		   Vg[4] = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		}
		else
		{
			coeff_2 = Vf[4];
			coeff_1 = 0.0;  // Normal gradient of T at wall

			Vg[4] = coeff_1*eta + coeff_2;
		}
	}

}

//------------------------------------------------------------------------------

inline
void HigherOrderFSI::setTurboG(double* Vf, double** dVf,
										 Vec3D Dir, double xi, double eta,
										 int dim, double* Vg)
{
	
	double coeff_1, coeff_2, coeff_3;
	
	for(int k=5; k<dim; ++k)
	{				
		if(viscQuadRcn)
		{
			coeff_3 = Vf[k];
			coeff_2 = dVf[k][0]*Dir[0] + dVf[k][1]*Dir[1] + dVf[k][2]*Dir[2]; 
			coeff_1 = (-coeff_2*xi - coeff_3)/(xi*xi); // Turbulent quantities are zero on the wall

			Vg[k] = coeff_1*eta*eta + coeff_2*eta + coeff_3;
		}
		else
		{
			coeff_2 = Vf[k];
			coeff_1 = -coeff_2/xi;

			Vg[k] = coeff_1*eta + coeff_2;
		}
	}

}

//------------------------------------------------------------------------------
// un is the normal velocity at stencil, ut is the tangential velocity at the stencil
// tgw1 is the same direction of ut
// tgw2=tgw1 X nwall
inline
void HigherOrderFSI::computeWallVersors(double *V1, Vec3D &nW, VarFcn *vf,
                                        Vec3D &tgW1, Vec3D &tgW2)
{

	double norm;
	const double tol = 1.0e-12;


	Vec3D uf = vf->getVelocity(V1);
	//uf += tol;//todo what is this??
	
	Vec3D un = (uf*nW)*nW; //normal velocity
	tgW1 = uf - un;         // tangential velocity direction

	norm = tgW1.norm();
    if(norm > tol)
	    tgW1 *= (1.0/norm);
    else{
    //std::cout << "else in HIGHORDERFSI" << std::endl;
        //tgW1 is any unit vector perpendicular to uW
        Vec3D e1(1.0,0.0,0.0),e2(0.0,1.0,0.0);
        e1 = nW^e1;
        if(e1.norm( ) > tol)
            tgW1 = e1/e1.norm();
        else {
            e2 = nW ^ e2;
            tgW1 = e2/e2.norm();
        }
        //std::cout << "uf " << uf[0] <<" "<<uf[1] <<" " <<uf[2] <<"un " << un[0] <<" "
        //          << un[1]<<" " <<un[2]<<"uW " << nW[0] <<" "<< nW[1]<<" " <<nW[2]<<" tgW1 " << tgW1[0] <<" " << tgW1[1]<<" " << tgW1[2] <<" "<<norm <<    std::endl;

    }

	tgW2 = tgW1 ^ nW;	
	norm = tgW2.norm();
	tgW2 *= (1.0/norm);

}

//------------------------------------------------------------------------------

inline
bool HigherOrderFSI::computeDuDTwf(double *V1, VarFcn *vf, double d2w, 
											  Vec3D &vWall, double TWall, 
											  Vec3D &tgW1, Vec3D &tgW2,
											  FemEquationTerm *fet, 
											  double &dudn, double &dTdn)
{

	Vec3D u1 = vf->getVelocity(V1);
	
	double un  =    u1 * tgW1;//tangential velocity
	double unW = vWall * tgW1;//tangential wall velocity

	double Du = un - unW;
    //dzh
    //std::cout << Du << std::endl;
	
	double rho = vf->getDensity(V1);
	double T   = vf->computeTemperature(V1);

	double DT = T - TWall;
	double utau = fet->computeNormDerivWallFcn(rho, T, Du, DT, d2w, dudn, dTdn);

	//double yp = d2w*utau*rho/1.0e-6;
	double cf_ = 2.0*rho*utau*utau;

	bool withWF = true;
	//if(yp <= 5.0) withWF = false;

	//fprintf(stdout, " %f\n", cf_); 

	return withWF;

}





inline
double HigherOrderFSI::vanAlbada(double a, double b)
{
    //a = Vj -vi
    //b = dVi

    double beta = 0.5, eps = 1.0e-16;

    b = (1 - 2*beta)*a + 2*beta*b;


    double a2 = a * a;
    double b2 = b * b;

    if (a*b > 0.0)
        return (a*(b2+eps) + b*(a2+eps)) / (a2 + b2 + 2.0*eps);
    else if (a*b == 0.0)
        return 0.5 * eps * (a + b) / (a2 + b2 + 2.0*eps);
    else
        return 0.0;

}


//----------------------------------------------------------------------------------------
// Vi(dVi)-------------Ve---(alpha)-------Vghost
// primitive value Vi , Vghost,
// dV = grad Vi*(Xghost - Xi) ~ Vghost - Vi if ij is true,
// dV = grad Vi*(Xi - Xghost) ~ Vi - Vghost if ij is false,
// alpha = |Ve - Vghost|/|Vi -Vghost|
// extrapolation to Ve
// if the limiter is on, use Van Albada average
// if the limiter is off, use pure extrapolation
inline
void HigherOrderFSI::safeExtrapolation(int dim, const double* Vi, const double* Vghost, const double *dV, bool ij, double alpha, double* Ve)
{
    double* dVij = new double[dim];
    for (int k = 0; k < dim; k++)
        dVij[k] = (ij ? dV[k]: -dV[k]);


    if(!limitExtrap || !Vghost) {
        for (int k = 0; k < dim; k++)
            Ve[k] = Vi[k] + (1.0 - alpha) * dVij[k];
    }else{
        for (int k = 0; k < dim; k++) 
            Ve[k] = Vi[k] + (1.0 - alpha) * vanAlbada(Vghost[k] - Vi[k], dVij[k]);
    }

//Fail Safe, if negative pressure or density, revert back to constant extrapolation
    if(Ve[0] <= 0.0 || Ve[4] <= 0.0)
        for (int k = 0; k < dim; k++)
            Ve[k] = Vi[k];


    delete []dVij;

}


//----------------------------------------------------------------------------------------
// Ve = V + beta*dV if no negative pressure or density, and HOtreatment is on
// else Ve = V
inline
void HigherOrderFSI::safeExtrapolation(int dim, const double* V, const double* dV, double beta, double* Ve)
{
    for (int k = 0; k < dim; k++)
        Ve[k] = V[k] + beta * dV[k];

//Fail Safe, if negative pressure or density, revert back to constant extrapolation
    if(Ve[0] <= 0.0 || Ve[4] <= 0.0 || !HOtreatment)
        for (int k = 0; k < dim; k++)
            Ve[k] = V[k];

}
