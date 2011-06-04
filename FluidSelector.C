#include "FluidSelector.h"

#include "Domain.h"


FluidSelector::FluidSelector(const int nPhases, IoData &ioData, Domain *domain)
{ 
  numPhases = nPhases; 
  fluidId  = 0;
  fluidIdn = 0;
  if(domain){
    fluidId  = new DistVec<int>(domain->getNodeDistInfo());
    fluidIdn = new DistVec<int>(domain->getNodeDistInfo());
  }
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
FluidSelector::~FluidSelector() { delete fluidId; delete fluidIdn; delete fluidIdnm1; delete fluidIdnm2; }
