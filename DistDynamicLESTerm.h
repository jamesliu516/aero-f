#ifndef _DIST_DYNAMIC_LES_TERM_H_
#define _DIST_DYNAMIC_LES_TERM_H_

#include <DistVector.h>
#include <Domain.h>
#include <IoData.h>

//------------------------------------------------------------------------

template<int dim>
class DistDynamicLESTerm {

 private:

  int               numLocSub;
  Domain            *domain;
  double            gam, R;

  DistSVec<double,dim> *VCap;
  DistSVec<double,16> *Mom_Test;
  DistSVec<double,6> *Eng_Test;

 public:

  DistDynamicLESTerm(IoData &, Domain *);
  ~DistDynamicLESTerm();

  void computeTestFilterValues(DistSVec<double,2> &, DistVec<double> &, DistSVec<double,3> &, 
                               DistSVec<double,dim> &);

};

//------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistDynamicLESTerm.C>
#endif

#endif
