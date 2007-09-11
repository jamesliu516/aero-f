#include <VolumicForceTerm.h>

#include<IoData.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
using std::min;
#endif

//------------------------------------------------------------------------------

VolumicForceTerm::VolumicForceTerm(IoData& iod)
{

  gravity = iod.bc.hydro.gravity;
  alpha = iod.bc.hydro.alpha;
  beta  = iod.bc.hydro.beta;
  dir[0] = cos(alpha)*cos(beta);
  dir[1] = cos(alpha)*sin(beta);
  dir[2] = sin(alpha);

}

//------------------------------------------------------------------------------

void VolumicForceTerm::computeVolumeTerm(double ctrlVol, double *V, double *flux)
{

  double momentum = V[0]*gravity*ctrlVol;
  double energy   = momentum*(V[1]*dir[0]+V[2]*dir[1]+V[3]*dir[2]);


  flux[0] = 0.0;
  flux[1] = -momentum*dir[0];
  flux[2] = -momentum*dir[1];
  flux[3] = -momentum*dir[2];
  flux[4] = -energy;

}

//------------------------------------------------------------------------------

void VolumicForceTerm::computeJacobianVolumeTerm(int dim, double ctrlVol,
						 double *V, double *jac)
{

  for (int i=0; i<dim*dim; i++)
    jac[i] = 0.0;
  
  jac[dim]  = gravity*dir[0];
  jac[2*dim] = gravity*dir[1];
  jac[3*dim] = gravity*dir[2];
  jac[4*dim+1] = jac[dim];
  jac[4*dim+2] = jac[2*dim];
  jac[4*dim+3] = jac[3*dim];

}

//------------------------------------------------------------------------------

