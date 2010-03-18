#ifndef _DIST_EXACT_RIEMANN_H_
#define _DIST_EXACT_RIEMANN_H_


class IoData;
class Domain;
class Communicator;
class SparseGridCluster;


template<class Scalar, int dim> class DistSVec;
template<int dim> class ExactRiemannSolver;
//------------------------------------------------------------------------------

template<int dim>
class DistExactRiemannSolver {

  MultiFluidData::TypePhaseChange phaseChangeType_;
  MultiFluidData::Method algorithmType_;
  int numLocSub;
  Domain *domain;

  DistSVec<double,dim> *riemannupdate; //node based
  DistVec<double> *weight;             //node based
  // riemannupdate is used only when using GFMPAR
  // it stores the future value at each cell if
  // that cell changes phases between two time steps
  DistSVec<double,dim-2> *interfacialWi; //edge based
  DistSVec<double,dim-2> *interfacialWj; //edge based

  ExactRiemannSolver<dim> **subExactRiemannSolver;

  SparseGridCluster *tabulationC;

public:
  DistExactRiemannSolver(IoData &iod, Domain *dom, VarFcn *vf);
  ~DistExactRiemannSolver();

  DistSVec<double,dim> *getRiemannUpdate() const { return riemannupdate; }
  DistVec<double> *getRiemannWeight() const { return weight; }

  void updatePhaseChange(DistSVec<double,dim> &V, DistVec<int> &fluidId, DistVec<int> &fluidIdn); 
  void storePreviousPrimitive(DistSVec<double,dim> &V, DistVec<int> &fluidId, DistSVec<double,3> &X);
  template<int dimLS>
  void avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin);

  ExactRiemannSolver<dim> &operator() (int i) 
    const { return *subExactRiemannSolver[i]; }

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistExactRiemannSolver.C>
#endif

#endif

