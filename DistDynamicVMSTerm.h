#ifndef _DIST_DYNAMIC_VMS_TERM_H_
#define _DIST_DYNAMIC_VMS_TERM_H_

#include <DistMacroCell.h>
#include <DistVector.h>
#include <Domain.h>
#include <DynamicVMSTerm.h>
#include <IoData.h>

class DistGeoState;
class FemEquationTerm;
class Domain;
class VarFcn;
template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<class Scalar, int dim> class DistSVec;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

//--------------------------------------------------------------------------

class Communicator;

template<int dim>
class DistDynamicVMSTerm {


 private:

  VarFcn *varFcn;
  int          scopeWidth, scopeDepth1, scopeDepth2;
  int          lastConfig;
  int          numLocSub;
  int          it0;
  int          lastIt;
  int          method;

  Domain            *domain;
  Communicator      *com;
  DistMacroCellSet  *macroCells;
  bool              **masterFlag;

  DistNodalGrad<dim, double> *ngrad;
  DistEdgeGrad<dim> *egrad;
  DynamicVMSTerm    *dvmst;

  DistSVec<double,dim> **VBar;
  DistSVec<double,1>   **volRatio;

  DistSVec<double,dim> *M;
  DistSVec<double,dim> *MBar;
  DistSVec<double,dim> *dWdt;
  DistSVec<double,dim> *dWBardt;
  DistSVec<double,dim> *RBar;
  DistSVec<double,dim> *S;

 public:

  DistDynamicVMSTerm( VarFcn *, IoData &, Domain *);
  ~DistDynamicVMSTerm();

  void compute(FluxFcn**, RecFcn*, FemEquationTerm *, int, DistVec<double> &, 
               DistBcData<dim> &, DistGeoState &, DistTimeState<dim> *, 
               DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, 
               DistSVec<double,dim> &, int, int);

  void computeCs(FluxFcn**, RecFcn*, FemEquationTerm *,int, DistVec<double> &,
                 DistBcData<dim> &, DistGeoState &, DistTimeState<dim> *,
                 DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                 DistSVec<double,dim> &, DistVec<double> *, int, int);

  void printDelRatios(DistVec<double> &);

};

//------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistDynamicVMSTerm.C>
#endif

#endif

