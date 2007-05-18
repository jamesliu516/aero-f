#ifndef _DIST_VOLUMIC_FORCE_TERM_H_
#define _DIST_VOLUMIC_FORCE_TERM_H_

template<int dim> class VolumicForceTerm;

//------------------------------------------------------------------------------

template<int dim>
class DistVolumicForceTerm {

  int lastConfig;
  int it0;
  int lastIt;
  int numLocSub;
  SubDomain** subDomain;
  VolumicForceTerm<dim>** subVolumicForceTerm;
  Communicator *com;


public:

  DistVolumicForceTerm(IoData&, Domain*);
  ~DistVolumicForceTerm();

  VolumicForceTerm<dim>& operator() (int i) const { return *subVolumicForceTerm[i]; }

  void computeVolumeTerm();
  void computeJacobianVolumeTerm();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistVolumicForceTerm.C>
#endif
                                                                                           
                                                                                           
#endif
