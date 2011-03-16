/* ProgrammedBurn.C

*/

#include "ProgrammedBurn.h"


template <int dim>
void ProgrammedBurn::setCJInitialState(ProgrammedBurn::Burn& B,VarFcn* vf,DistSVec<double,dim>& U,DistVec<int>& fid,DistVec<int>& fidn) {
  
  if (B.x0subdom >= 0 && B.x0id >= 0) {
    double V0[dim];
    assert(dim == 5);
    V0[1] = V0[2] = V0[3] = 0.0;
    V0[0] = B.pgData->cjDensity;
    V0[4] = B.pgData->cjPressure;

    std::cout << "Setting burn values = " << B.pgData->cjDensity << " " << B.pgData->cjPressure << std::endl;
    double* UU = U.subData(B.x0subdom)[B.x0id];
    vf->primitiveToConservative(V0,U.subData(B.x0subdom)[B.x0id], B.pgData->burnedEOS);
    fidn.subData(B.x0subdom)[B.x0id] = fid.subData(B.x0subdom)[B.x0id] = B.pgData->burnedEOS;
  }
    
}

template <int dim>
void ProgrammedBurn::setCurrentTime(double t,VarFcn* vf,DistSVec<double,dim>& U,DistVec<int>& fid,DistVec<int>& fidn) {

  for (int j = 0; j < myBurns.size(); ++j) {
    Burn& B = myBurns[j];
    if (B.pgData->ignitionTime <= t && 
	!B.ignited) {
      B.ignited = true;
      setCJInitialState<dim>(B, vf, U, fid,fidn);
    }
  }
}
