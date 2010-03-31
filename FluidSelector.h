#ifndef _FLUID_SELECTOR_H_
#define _FLUID_SELECTOR_H_

#include "DistVector.h"
#include "Vector.h"
#include "IoData.h"

#include <assert.h>
//#define NDEBUG // if commented, assert statements are evaluated

class Domain;

//--------------------------------------------------------------------------
//
// This algorithm class allows to determine in which fluid a node is,
// by determining the fluid identification from the different 
// available level sets.
//
// The convention is as follows. For n different fluids, only
// n-1 (=dim) level-sets are necessary.
// Phi[i] is positive where fluid i+1 (=varFcn[i+1]) is found. 
// fluid 0 (=varFcn[0]) is found where all levelsets are negative.
//--------------------------------------------------------------------------

class FluidSelector {

  int numPhases;

public:
  DistVec<int> *fluidId;
  DistVec<int> *fluidIdn;
  DistVec<int> *fluidIdnm1;
  DistVec<int> *fluidIdnm2;

public:

  FluidSelector(const int nPhases, IoData &ioData, Domain *domain = 0);
  ~FluidSelector();

  template<int dim>
  void initializeFluidIds(DistSVec<double,dim> &Phin, DistSVec<double,dim> &Phinm1, DistSVec<double,dim> &Phinm2){
    getFluidId(Phin);
    *fluidIdn = *fluidId;
    if(fluidIdnm1) getFluidId(*fluidIdnm1, Phinm1);
    if(fluidIdnm2) getFluidId(*fluidIdnm2, Phinm2);
  }
  void update(){
    if(fluidIdnm2) *fluidIdnm2 = *fluidIdnm1;
    if(fluidIdnm1) *fluidIdnm1 = *fluidIdn;
    *fluidIdn = *fluidId;
  }

  void getFluidId(DistVec<double> &Phi){
    int numLocSub = Phi.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double *phi = Phi.subData(iSub);
      int    *tag = fluidId->subData(iSub);
      for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
        tag[iNode] = (phi[iNode]<0.0) ? 0 : 1; 
    }
  }

  template<int dim>
  void getFluidId(DistSVec<double,dim> &Phi){
    assert(dim==numPhases-1);
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



// ---- obtaining fluidId from levelset phi ---- //
// ---- and vice-versa                      ---- //

  /* Node-Level Operators */
  void getFluidId(int &tag, double phi){ tag = (phi<0.0) ? 0 : 1; }

  void getFluidId(int &tag, double *phi){
    tag = 0;
    for(int i=0; i<numPhases-1; i++)
      if(phi[i]>0.0) {tag = i+1; return; }
  }

  template<int dim>
  void getFluidId(int &tag, double *phi){
    assert(dim==numPhases-1);
    tag = 0;
    for(int i=0; i<dim; i++)
      if(phi[i]>0.0) {tag = i+1; return; }
  }

  int getLevelSetDim(int fluidId1, int fluidId2){
    if(fluidId1 == fluidId2){
      fprintf(stdout, "*** Error: getLevelSetDim should not be called when there is no interface.\n");
      exit(1);
    }
    if(fluidId1 < 0 || fluidId2 < 0){
      fprintf(stdout, "*** Error: arguments of getLevelSetDim (%d %d)should be positive\n", fluidId1, fluidId2);
      exit(1);
    }
    if(fluidId1 * fluidId2 != 0){
      fprintf(stdout, "*** Error: it  is assumed that all interfaces are between any fluid and fluid 0\n");
      exit(1);
    }
    int fId = fluidId1 + fluidId2;
    return fId-1;
    
  }


  /* Non-Distributed Operators */
  void getFluidId(Vec<int> &tag, Vec<double> &phi){
    for(int iNode=0; iNode<phi.size(); iNode++)
      tag[iNode] = (phi[iNode]<0.0) ? 0 : 1; 
  }

  template<int dim>
  void getFluidId(Vec<int> &tag, SVec<double,dim> &phi){
    assert(dim==numPhases-1);
    for(int iNode=0; iNode<phi.size(); iNode++){
      tag[iNode] = 0;
      for(int i=0; i<dim; i++)
        if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
    }
  }
    

  /* Distributed Operators */
  void getFluidId(DistVec<int> &Tag, DistVec<double> &Phi){
    int numLocSub = Phi.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double *phi = Phi.subData(iSub);
      int    *tag = Tag.subData(iSub);
      for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
        tag[iNode] = (phi[iNode]<0.0) ? 0 : 1; 
    }
  }

  template<int dim>
  void getFluidId(DistVec<int> &Tag, DistSVec<double,dim> &Phi){
    assert(dim==numPhases-1);
    int numLocSub = Phi.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double (*phi)[dim] = Phi.subData(iSub);
      int     *tag       = Tag.subData(iSub);
      for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
        tag[iNode] = 0;
        for(int i=0; i<dim; i++)
          if(phi[iNode][i]>0.0) { tag[iNode] = i+1; break; }
      }
    }
  }

// Debug
  void printFluidId(){
    int numLocSub = fluidId->numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      int    *tag = fluidId->subData(iSub);
      for (int i=0; i<fluidId->subSize(iSub); i++)
        fprintf(stdout, "fluidId[%d] = %d\n", i, tag[i]);
    }
  }

};

//------------------------------------------------------------------------------

#endif
