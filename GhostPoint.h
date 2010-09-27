#ifndef _GHOST_POINT_H_
#define _GHOST_POINT_H_

//#include<VarFcn.h>
#include <Vector.h>
#include <Vector3D.h>
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;

template<class Scalar> class Vec;

template<int dim>
class GhostPoint {
  Vec<double> Vg; // Sum of the primitive States at the ghost-point. 
  int ng; // Number of neighbours in the fluid. State at GP is then equal to Vg/ng.
  // After all GP have been populated, Vg /= ng and ng=1.
  int ghostTag; // We store here the tag of the surrounding nodes. All the tags of the neighbours 
  // should be the same. In the case of a complex multiphase flow simulation with Fluid/Structure 
  // Interaction, this might be no longer true. To be done...
 public:
  ~GhostPoint() {};
 GhostPoint() :
  Vg(dim)
  {
    ng = 0;
    Vg = 0.0;
    ghostTag = -2; // Inactive nodes tag
  }
  void addNeighbour(Vec<double> &Vi,double distanceRate, Vec3D interfaceVelocity, int tag)
  {
    // Ui is the state at the neighbour. 
    // distanceRate is the rate of the distances from the GP and the neighbour to the interface = dg/di\
    
    ng++;

    // We want the velocity to be zero at the interface and we obtain the 
    // state at the GP by linear interpolation.
    Vg[0]   += Vi[0];
    for(int i=1;i<4;++i) Vg[i] += interfaceVelocity[i-1] - distanceRate*(Vi[i]-interfaceVelocity[i-1]);
    Vg[4]   += Vi[4];
    if(dim == 6) // Turbulent Viscosity
      {
	//	Vg[5] -= distanceRate*Vi[5];
	Vg[5] = 0.0;
      }

    // Tag check
    if(ghostTag < 0)
      {
	ghostTag = tag;
      }
    else if(ghostTag != tag)
      {
	fprintf(stderr,"We have a ghost node here with two active neighbours having different tags\n");
	fprintf(stderr,"ghostTag: %i, neighbourTag: %i",ghostTag,tag);
	exit(-1);
      }
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
