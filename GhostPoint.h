#ifndef _GHOST_POINT_H_
#define _GHOST_POINT_H_
  
//#include<VarFcn.h>
#include<Vector.h>
#include<iostream>
using std::cout;
using std::endl;

template<class Scalar> class Vec;

template<int dim>
class GhostPoint {
  Vec<double> Vg; // Sum of the primitive States at the ghost-point. 
  int ng; // Number of neighbours in the fluid. State at GP is then equal to Vg/ng.
  // After all GP have been populated, Vg /= ng and ng=1.
 public:
  ~GhostPoint() {};
 GhostPoint() :
  Vg(dim)
  {
    ng = 0;
    Vg = 0.0;
  }
  void addNeighbour(Vec<double> &Vi,double distanceRate, int tag=0)
  {
    // Ui is the state at the neighbour. 
    // distanceRate is the rate of the distances from the GP and the neighbour to the interface = dg/di\
    
    ng++;

    // We want the velocity to be zero at the interface and we obtain the 
    // state at the GP by linear interpolation.
    Vg[0]   += Vi[0];
    for(int i=1;i<dim-1;++i) Vg[i] -= distanceRate*Vi[i];
    Vg[dim-1] += Vi[dim-1];
  }
  double* getPrimitiveState()
  {
    //    fprintf(stderr,"State: %f %f %f %f %f",Vg.v[0],Vg.v[1],Vg.v[2],Vg.v[3],Vg.v[4]);
    return Vg.v;    
  }
  void reduce()
  {
    Vg /= (double) ng;
    ng = 1;
  }
};

#endif
