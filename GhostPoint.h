#ifndef _GHOST_POINT_H_
#define _GHOST_POINT_H_

//#include<VarFcn.h>
#include <Vector.h>
#include <Vector3D.h>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

template<class Scalar> class Vec;

template<int dim>
class GhostPoint {
 public:
  int ng_last;
  Vec<double> Vg; // Sum of the primitive States at the ghost-point. 
  Vec<double> dVg; // Sum of the second order correction term. 
  int ng; // Number of neighbours in the fluid. State at GP is then equal to Vg/ng.
  // After all GP have been populated, Vg /= ng and ng=1.
  int ghostTag; // We store here the tag of the surrounding nodes. All the tags of the neighbours 
  // should be the same. In the case of a complex multiphase flow simulation with Fluid/Structure 
  // Interaction, this might be no longer true. To be done...
  ~GhostPoint() {};
 GhostPoint() :
  Vg(dim), dVg(dim)
  {
    ng = 0;
    Vg = 0.0;
    dVg = 0.0;
    ghostTag = -2; // Inactive nodes tag
  }
  GhostPoint<dim> & operator=(const GhostPoint<dim> &GP)
    {
      Vg = GP.Vg;
      dVg = GP.dVg;
      ng = GP.ng;
      ng_last = ng;
      ghostTag = GP.ghostTag;
      return *this;
    }
  GhostPoint<dim> & operator+=(const GhostPoint<dim> &GP)
    {
      if(ghostTag<0) ghostTag = GP.ghostTag;
      else if(ghostTag != GP.ghostTag) 
	{
	  fprintf(stderr,"The two ghost States refer to different Fluids\n");
	  fprintf(stderr,"ghostTag: %i, GP.ghostTag: %i",ghostTag,GP.ghostTag);
	  exit(-1);
	}
      Vg += GP.Vg;
      dVg += GP.dVg;
      ng += GP.ng;
      ng_last = ng;
      return *this;
    }
  void addNeighbour(Vec<double> &Vi,Vec<double> &dVi,double distanceRate, Vec3D interfaceVelocity, int tag)
  {
    // Ui is the state at the neighbour. 
    // distanceRate is the rate of the distances from the GP and the neighbour to the interface = dg/di\
    
    ng++;
    ng_last = ng;

    // We want the velocity to be zero at the interface and we obtain the 
    // state at the GP by linear interpolation.
    Vg[0]   += Vi[0];
    dVg[0]  += dVi[0];
    for(int i=1;i<4;++i) Vg[i]  += 2.0*interfaceVelocity[i-1] - Vi[i];
    for(int i=1;i<4;++i) dVg[i] += (1.0 - distanceRate)*(Vi[i]-interfaceVelocity[i-1]);
    Vg[4]   += Vi[4];
    dVg[4]  += dVi[4];
    if(dim == 6) // One Equation Turbulent Model
      {
	//	Vg[5] -= distanceRate*Vi[5];
	Vg[5] = 0.0;
	dVg[5] = 0.0;
      }
    else if(dim == 7) // Two Equations Turbulent Model
      {
	Vg[5] = 0.0;
	Vg[6] = 0.0;
	dVg[5] = 0.0;
	dVg[6] = 0.0;
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
    dVg /= (double) ng;
    Vg = Vg + dVg;
    ng = 1;
  }
  
  int lastCount() { return ng_last; }
};

#endif
