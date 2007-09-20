#ifndef _VOLUMIC_FORCE_TERM_H_
#define _VOLUMIC_FORCE_TERM_H_


#include <stdlib.h>
#include <stdio.h>

class IoData;

//------------------------------------------------------------------------------

class VolumicForceTerm {

  double gravity;
  double dir[3];
  double alpha, beta;

public:

  VolumicForceTerm(IoData &);
  ~VolumicForceTerm() {}

  void computeVolumeTerm(double ctrlVol, double *V, double *R);
  void computeJacobianVolumeTerm(int dim, double ctrlVol, double *V, double *jac);

// Included (MB)
  void computeDerivativeOfVolumeTerm(double ctrlVol, double dCtrlVol, double *V, double *dV, double *dR);

};
  
//------------------------------------------------------------------------------

#endif
