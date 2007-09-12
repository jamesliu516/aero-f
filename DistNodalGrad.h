#ifndef _DIST_NODAL_GRAD_H_
#define _DIST_NODAL_GRAD_H_

#include <IoData.h>

class RecFcn;
class Domain;

#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim, class Scalar>
class DistNodalGrad {

  SchemeData::Gradient typeGradient;

  int failSafeNewton;

  int lastConfig;

  int numLocSub;

  DistVec<bool> *tag;
  DistVec<bool> *backuptag;

  DistSVec<Scalar,dim> *Vmin;
  DistSVec<Scalar,dim> *Vmax;
  DistSVec<Scalar,dim> *phi;

  DistSVec<Scalar,3>* sensor;
  DistVec<Scalar>* sigma;

  DistSVec<double,6> *R;
  
  DistSVec<double,3> *wii;
  DistSVec<double,3> *wij;
  DistSVec<double,3> *wji;

  DistSVec<Scalar,dim> *ddx;
  DistSVec<Scalar,dim> *ddy;
  DistSVec<Scalar,dim> *ddz;

  Domain *domain;

  NodalGrad<dim, Scalar> **subNodalGrad;

public:

  DistNodalGrad(IoData &, Domain *);
  DistNodalGrad(IoData &, Domain *, int);
  ~DistNodalGrad();

  NodalGrad<dim, Scalar> &operator() (int i) const { return *subNodalGrad[i]; }

  void computeWeights(DistSVec<double,3> &);

  template<class Scalar2>
  void compute(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &);

  template<class Scalar2>
  void compute(int, DistSVec<double,3> &, DistVec<double> &,
               DistVec<double> &, DistSVec<Scalar2,dim> &);

  //template<class Scalar2>
  //void computeLS(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &);

  void compute(int config, DistSVec<double,3> &X, DistSVec<double,dim> &Psi);

  template<class Scalar2> 
  void computeT(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &, 
       DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar2>
  void limit(RecFcn *, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &);

  void fix(DistSVec<bool,2>&);
  void fix(DistSVec<int,2>&);

  void resetTag();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistNodalGrad.C>
#endif

#endif
