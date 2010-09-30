#include "LevelSet/LevelSetStructure.h"
#include <assert.h>

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::initializeFluidIds(DistSVec<double,dim> &Phin, DistSVec<double,dim> &Phinm1, DistSVec<double,dim> &Phinm2){
  getFluidId(Phin);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) getFluidId(*fluidIdnm1, Phinm1);
  if(fluidIdnm2) getFluidId(*fluidIdnm2, Phinm2);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistSVec<double,dim> &Phi){
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      tag[iNode] = 0;
      for(int i=0; i<dim; i++)
        if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(int &tag, double *phi){
  assert(dim<=numPhases-1);
  tag = 0;
  for(int i=0; i<dim; i++)
    if(phi[i]>0.0) {tag = i+1; return; }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(Vec<int> &tag, SVec<double,dim> &phi){
  assert(dim<=numPhases-1);
  for(int iNode=0; iNode<phi.size(); iNode++){
    tag[iNode] = 0;
    for(int i=0; i<dim; i++)
      if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistVec<int> &Tag, DistSVec<double,dim> &Phi, DistVec<int>* fsId){
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = Tag.subData(iSub);
    int     *fsid      = fsId ? fsId->subData(iSub) : 0;
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      tag[iNode] = 0;
      if(fsid && fsid[iNode]!=0) //isolated by structure. use fsId as fluidId.
        tag[iNode] = fsid[iNode];
      else for(int i=0; i<dim; i++)
          if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFS(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV)
{
  assert(dim<=numPhases-1);
  DistVec<int> &fsId(distLSS->getStatus());
#pragma omp parallel for
  for (int iSub=0; iSub<PhiV.numLocSub(); ++iSub) {
    Vec<int> &subfsId(fsId(iSub));
    Vec<int> &subId((*fluidId)(iSub));
    SVec<double,dim> &subPhiV(PhiV(iSub));
    LevelSetStructure &LSS((*distLSS)(iSub));

    for(int i=0; i<subPhiV.size(); i++) {
      int Id = LSS.fluidModel(0.0,i);
      bool swept = LSS.isSwept(0.0,i);

      if(Id==0) { // not isolated by structure. need to consider level-set
        if(swept) {
          fluidId[i] = 0;
          for(int k=0; k<dim; k++)
            if(subPhiV[i][k]>0.0) {
              fluidId[i] = k+1;
              break;
            }
        } else {/* not swept. nothing to be done :) */}
      } else // isolated by structure. Id determined by intersector
        fluidId[i] = Id;
    }
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFF(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi)
{
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
  DistVec<int> &fsId(distLSS->getStatus());
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    int     *fsid      = fsId.subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      if(fsid[iNode]!=0) {
        if(fsid[iNode]!=tag[iNode]){
          fprintf(stderr,"This is WRONG!\n"); exit(-1);}
        continue;
      }
      tag[iNode] = 0;
      for(int i=0; i<dim; i++)
        if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
    }
  }
}

//------------------------------------------------------------------------------

