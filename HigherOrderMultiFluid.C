/* HigherOrderMultiFluid.C

 */

inline
HigherOrderMultiFluid::HigherOrderMultiFluid() {

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
      Vi[i] += static_cast<CutCellStateData<dim>*>(C->cutCellData)->dV[fidi][i][l]*(iloc[i]-x0[i]);
      Vj[i] += static_cast<CutCellStateData<dim>*>(C->cutCellData)->dV[fidj][i][l]*(iloc[i]-x0[i]);
    }
  }
}

template <int dim>
void HigherOrderMultiFluid::
setCutCellFlags(int lsdim, Vec<int>& status) {

  for (int i = 0; i < cutCells.size(); ++i) {
   
    cutCells[i] = new CutCellState;
    
    cutCells[i]->cutCellData = new CutCellStateData<dim>;
    cutCells[i]->fid1 = 0;
    cutCells[i]->fid2 = lsdim+1;

    numCutCells++;
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

inline
int HigherOrderMultiFluid::getNumCutCells() {

  return numCutCells;
}
