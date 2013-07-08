/* HigherOrderMultiFluid.C

 */

#include <iostream>

#include "V6NodeData.h"
#include <Elem.h>

//#include <EdgeGrad.h>


inline
HigherOrderMultiFluid::HigherOrderMultiFluid(Vec<CutCellState*>& vs) : cutCells(vs) {
  
  cutCells = static_cast<CutCellState*>(0);
  lastPhaseChangeState = NULL;
}

inline
HigherOrderMultiFluid::~HigherOrderMultiFluid() {

}

/*template<class Scalar, int dim,int dimLS>
HigherOrderMultiFluid<Scalar,dim,dimLS>::HigherOrderMultiFluid() {

}

template<class Scalar, int dim,int dimLS>
HigherOrderMultiFluid<Scalar,dim,dimLS>::~HigherOrderMultiFluid() {

}


template<class Scalar, int dim,int dimLS>
void HigherOrderMultiFluid<Scalar,dim,dimLS>::
setVectors(DistVec<int>* fid, DistSVec<Scalar,dim>* _u,
           DistSVec<Scalar,dimLS>* _ls, DistNodalGrad<Scalar,dim>* _du,
           DistSVec<Scalar,dimLS>* _dls) {

  myFid = fid;
  U = _u;
  levelSet = _ls;
  nodalGrad = _du;
  levelSetGrad = _dls;
  }*/

/*
template<class Scalar, int dim,int dimLS>
void HigherOrderMultiFluid<Scalar,dim,dimLS>::
computeLevelSetIntersection(const double x0[3], const double x1[3],
                            int subId, int nodeId, double iloc[3],
                            int lsid) {

  double dir[3];
  for (int i = 0; i < 3; ++i) 
    dir[i] = x1[i] - x0[i];

  NodalGrad<Scalar,dimLS>& g = *(*levelSetGrad)(iSub);
  double phi = (*levelSet)(iSub)[nodeId][lsid];
  double phigrad[3] = { g.getX()[nodeId][lsid],
                        g.getY()[nodeId][lsid],
                        g.getZ()[nodeId][lsid] };

  // expand phi as phi+gradphi*(x-x1) = 0
  // where x is x0+(x1-x0)*s => phi + gradphi*(x1-x0)*(1-s) = 0
  // => 1-s = -phi/(gradphi*(x1-x0))
  // => s = 1+phi/(gradphi*(x1-x0))

  double dot = 0.0;
  for (int i = 0; i < 3; ++i)
    dot += gradphi[i]*dir[i];
  double s = 1.0+phi/dot;

  for (int i = 0; i < 3; ++i)
    iloc[i] = x0[i] + dir[i]*s;
}

template<class Scalar, int dim,int dimLS>
void HigherOrderMultiFluid<Scalar,dim,dimLS>::
computeExtrapolation(const double x0[3], 
                     const double x1[3],
                     const double iloc[3],
                     std::pair<int,int> node0,
                     std::pair<int,int> node1, 
                     Scalar Uext1[dim],
                     Scalar Uext2[dim]) {

  CutCellState* C = myCutCellState(node1.first)[node1.second];
  if (C) {

    int fid = myFid(node0.first)[node0.second];
    int lsid = 0;
    int othfid = 0;
    if (C->fid1 == fid) {
      othfid = C->fid2;
    } else {
      othfid = C->fid1;
    }
 
    for (int k = 0; k < dim; ++k) {

      Uext1[k] = C->U[fid][k];
      Uext2[k] = C->U[othfid][k];
      
      for (int i = 0; i < 3; ++i) {
        Uext1[k] += C->dU[fid][k][i]*(iloc[i]-x1[i]);
        Uext2[k] += C->dU[othfid][k][i]*(iloc[i]-x1[i]);
      }
    }
    
  } else {

     // node0.first = node1.first
     SVec<Scalar,dim>& dX = (*nodalGrad)(node0.first).getX();
     SVec<Scalar,dim>& dY = (*nodalGrad)(node0.first).getY();
     SVec<Scalar,dim>& dZ = (*nodalGrad)(node0.first).getZ();
     SVec<Scalar,dim>& u = (*U)(node0.first);
     const Scalar* u_v_0 = u[node0.second];
     const Scalar* u_v_1 = u[node1.second];
     for (int k = 0; k < dim; ++k) {
 
       Uext1[k] = u_v_0[k] + dX[node0.second][k]*(iloc[0]-x0[0]) +
                             dY[node0.second][k]*(iloc[1]-x0[1]) +
                             dZ[node0.second][k]*(iloc[2]-x0[2]); 
       Uext2[k] = u_v_1[k] + dX[node1.second][k]*(iloc[0]-x1[0]) +
                             dY[node1.second][k]*(iloc[1]-x1[1]) +
                             dZ[node1.second][k]*(iloc[2]-x1[2]); 
     } 
  }
}

template<class Scalar, int dim,int dimLS>
void HigherOrderMultiFluid<Scalar,dim,dimLS>::
computeInterpolation(const double x0[3],
                     const double iloc[3],
                     const double xmid[3],
                     const Scalar U0[dim],
                     const Scalar Ustar[dim],
                     Scalar Ui[dim]) {

  double dxg[3] = {iloc[0]-x0[0],iloc[1]-x0[1],iloc[2]-x0[2]};
  double dxc[3] = {xmid[0]-x0[0],xmid[1]-x0[0],xmid[2]-x0[2]};
  double lg = sqrt(dxg[0]*dxg[0]+dxg[1]*dxg[1]+dxg[2]*dxg[2]);
  double lc = sqrt(dxc[0]*dxc[0]+dxc[1]*dxc[1]+dxc[2]*dxc[2]);
  double s = lc/lg;
  for (int k = 0; k < dim; ++k) {

    Ui[k] = s*Ustar[k] + (1.0-s)*U0[k];
  }

}
 
*/
template<int dim>
void HigherOrderMultiFluid::
computeCutCellExtrapolations(int cutCellId,int fidi,int fidj, const double iloc[3],
                             double* Vi, double* Vj, SVec<double,3>& X) {


  const double* x0 = X[cutCellId];

  CutCellState* C = cutCells[cutCellId];

  for (int i = 0; i < dim; ++i) {

    Vi[i] = static_cast<CutCellStateData<dim>*>(C->cutCellData)->V[fidi][i];
    Vj[i] = static_cast<CutCellStateData<dim>*>(C->cutCellData)->V[fidj][i];
    for (int l = 0; l < 3; ++l) {
      Vi[i] += static_cast<CutCellStateData<dim>*>(C->cutCellData)->dV[fidi][i][l]*(iloc[l]-x0[l]);
      Vj[i] += static_cast<CutCellStateData<dim>*>(C->cutCellData)->dV[fidj][i][l]*(iloc[l]-x0[l]);
    }
  }
}

template <int dim>
void HigherOrderMultiFluid::
setCutCellFlags(int lsdim, Vec<int>& status) {

  for (int i = 0; i < status.size(); ++i) {
    
    if (status[i]) {
      if (!cutCells[i]) {
	cutCells[i] = new CutCellState;
      
	cutCells[i]->cutCellData = new CutCellStateData<dim>;
      }
      cutCells[i]->fid1 = 0;
      cutCells[i]->fid2 = lsdim+1;
      
      numCutCells++;
    } else if (cutCells[i]) {

      delete (static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData));
      delete (cutCells[i]);
      cutCells[i] = 0;
      
    }
  }
}

template <int dim>
void HigherOrderMultiFluid::
clearCutCellFlags() {

  for (int i = 0; i < cutCells.size(); ++i) {

    if (cutCells[i]) {

      delete (static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData));
      delete (cutCells[i]);
      cutCells[i] = 0;
    }
    
  }      
  numCutCells = 0;
}


template <int dim>
void HigherOrderMultiFluid::printCutCellData(int i) {

  CutCellStateData<dim>* data = static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData);
  std::cout << "Cut cell data: " << std::endl;
  for (int l = 0; l < 2; ++l) {
    std::cout << "fluid " << l << std::endl;
    for (int k = 0; k < dim; ++k) {
      std::cout << data->V[l][k] << " " << data->dV[l][k][0] << " " << data->dV[l][k][1] << " " << data->dV[l][k][2] << std::endl;
    }
  }
}

template<int dim>
void HigherOrderMultiFluid::getCutCellData(int cut,int fid,double V[dim], double x[dim][3]) {

  CutCellStateData<dim>* data = static_cast<CutCellStateData<dim>*>(cutCells[cut]->cutCellData);
  int l = (fid == 0 ? 0 : 1);
  
  for (int k = 0; k < dim; ++k) {
    V[k] = data->V[l][k];
    x[k][0] = data->dV[l][k][0];
    x[k][1] = data->dV[l][k][1];
    x[k][2] = data->dV[l][k][2];
  }
}

inline
int HigherOrderMultiFluid::getNumCutCells() {

  return numCutCells;
}

template <int dim>
void HigherOrderMultiFluid::storeCutCellData(SVec<double,dim>* cutCell[2],
					     NodalGrad<dim,double>* cutGrad[2],
					     Vec<int>* counts[2]) {

  for (int i = 0; i < cutCells.size(); ++i) {

    if (!cutCells[i])
      continue;

    CutCellStateData<dim>* data = static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData);

    for (int l = 0; l < 2; ++l) {
      if (counts[l]->operator[](i) > 0) {
        for (int k = 0; k < dim; ++k) {
	  data->V[l][k] = cutCell[l]->operator[](i)[k] / counts[l]->operator[](i);
	  data->dV[l][k][0] = cutGrad[l]->getX()[i][k] /  counts[l]->operator[](i);
	  data->dV[l][k][1] = cutGrad[l]->getY()[i][k] /  counts[l]->operator[](i);
	  data->dV[l][k][2] = cutGrad[l]->getZ()[i][k] /  counts[l]->operator[](i);
        }
      } else {

        delete (static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData));
        delete (cutCells[i]);
        cutCells[i] = 0;
        break;
      }
    }
  }
}

template<int dim>
void HigherOrderMultiFluid::setCutCellData(SVec<double,dim>& V, Vec<int>& fid) {

  for (int i = 0; i < cutCells.size(); ++i) {

    if (!cutCells[i])
      continue;

    CutCellStateData<dim>* data = static_cast<CutCellStateData<dim>*>(cutCells[i]->cutCellData);

    int fididx = (fid[i] == 0 ? 0 : 1);
    for (int k = 0; k < dim; ++k) {
      V[i][k] = data->V[fididx][k];
    }
  }     
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
    xf[k] = (X[n0_loc][k] + X[n1_loc][k] + X[n2_loc][k])/3.0;

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
      r = std::min<double>(r, (Vfg[k][0]*vec[0]+
			       Vfg[k][1]*vec[1]+
			       Vfg[k][2]*vec[2])/(V[i][k]-Vf[k])-1.0);
    else
      r = 0.0;
  }
  
  return std::max<double>(r,0.0);
}
