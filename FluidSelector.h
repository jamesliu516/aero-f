#ifndef _FLUID_SELECTOR_H_
#define _FLUID_SELECTOR_H_

#include "DistVector.h"
#include "Vector.h"
#include "IoData.h"
#include "Domain.h"

#include <assert.h>

//#define NDEBUG // if commented, assert statements are evaluated

//--------------------------------------------------------------------------
//
// This algorithm class allows to determine in which fluid a node is,
// by determining the fluid identification from the different 
// available level sets.
//
//--------------------------------------------------------------------------

class FluidSelector {

  int numPhases;

public:
  DistVec<int> fluidId;
  DistVec<int> fluidIdn;
  DistVec<int> *fluidIdnm1;
  DistVec<int> *fluidIdnm2;

public:
  
  FluidSelector(const int nPhases, IoData &ioData, Domain *domain) :
    fluidId(domain->getNodeDistInfo()), fluidIdn(domain->getNodeDistInfo())
  { 
    numPhases = nPhases; 
    fluidIdnm1 = 0;
    fluidIdnm2 = 0;
    if(ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF){
      fluidIdnm1 = new DistVec<int>(domain->getNodeDistInfo());
      *fluidIdnm1 = 0;
    }
    else if(ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF){
      fluidIdnm1 = new DistVec<int>(domain->getNodeDistInfo());
      *fluidIdnm1 = 0;
      fluidIdnm2 = new DistVec<int>(domain->getNodeDistInfo());
      *fluidIdnm2 = 0;
    }
  }
  ~FluidSelector() { delete fluidIdnm1; delete fluidIdnm2; }


  void initializeFluidIds(DistSVec<double,1> &Phin, DistSVec<double,1> &Phinm1, DistSVec<double,1> &Phinm2){
    getFluidId(Phin);
    fluidIdn = fluidId;
    if(fluidIdnm1) getFluidId(*fluidIdnm1, Phinm1);
    if(fluidIdnm2) getFluidId(*fluidIdnm2, Phinm2);
  }
  void update(){
    if(fluidIdnm2) *fluidIdnm2 = *fluidIdnm1;
    if(fluidIdnm1) *fluidIdnm1 = fluidIdn;
    fluidIdn = fluidId;
  }

  void getFluidId(DistVec<double> &Phi){
    int numLocSub = Phi.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double *phi = Phi.subData(iSub);
      int    *tag = fluidId.subData(iSub);
      for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
        tag[iNode] = (phi[iNode]>=0.0) ? 0 : 1; 
    }
  }

  template<int dim>
  void getFluidId(DistSVec<double,dim> &Phi){
    assert(dim==numPhases-1);
    int numLocSub = Phi.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double (*phi)[dim] = Phi.subData(iSub);
      int     *tag       = fluidId.subData(iSub);
      for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
        tag[iNode] = numPhases-1;
        for(int i=0; i<dim; i++)
          if(phi[iNode][i]>=0.0) { tag[iNode] = i; break; }
      }
    }
  }



// ---- obtaining fluidId from levelset phi ---- //

  /* Node-Level Operators */
  void getFluidId(int &tag, double phi){ tag = (phi>=0.0) ? 0 : 1; }

  void getFluidId(int &tag, double *phi){
    for(int i=0; i<numPhases-1; i++)
      if(phi[i]>=0.0) {tag = i; return; }
    tag = numPhases-1;
  }

  template<int dim>
  void getFluidId(int &tag, double *phi){
    assert(dim==numPhases-1);
    for(int i=0; i<dim; i++)
      if(phi[i]>=0.0) {tag = i; return; }
    tag = dim;
  }


  /* Non-Distributed Operators */
  void getFluidId(Vec<int> &tag, Vec<double> &phi){
    for(int iNode=0; iNode<phi.size(); iNode++)
      tag[iNode] = (phi[iNode]>=0.0) ? 0 : 1; 
  }

  template<int dim>
  void getFluidId(Vec<int> &tag, SVec<double,dim> &phi){
    assert(dim==numPhases-1);
    for(int iNode=0; iNode<phi.size(); iNode++){
      tag[iNode] = numPhases-1;
      for(int i=0; i<dim; i++)
        if(phi[iNode][i]>=0.0) { tag[iNode] = i; break; }
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
        tag[iNode] = (phi[iNode]>=0.0) ? 0 : 1; 
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
        tag[iNode] = numPhases-1;
        for(int i=0; i<dim; i++)
          if(phi[iNode][i]>=0.0) { tag[iNode] = i; break; }
      }
    }
  }

// Debug
  void printFluidId(){
    int numLocSub = fluidId.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      int    *tag = fluidId.subData(iSub);
      for (int i=0; i<fluidId.subSize(iSub); i++)
        fprintf(stdout, "fluidId[%d] = %d\n", i, fluidId[i]);
    }
  }

};

//------------------------------------------------------------------------------

#endif
