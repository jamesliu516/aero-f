#ifndef _GHOST_POINT_H_
#define _GHOST_POINT_H_

#include<VarFcn.h>
#include <Vector.h>
#include <Vector3D.h>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

class VarFcn;

template<class Scalar> class Vec;

template<int dim>
class GhostPoint{
 protected:
  VarFcn *varFcn;

 public:
  double* Vg;	// Sum of weighted states (rho,u,v,w,T) at the ghost-point. 
  		// After population, it is set to the state at the ghost-point.
  double *V;    // Stores the final primitive states
  double* Ws; 	// Sum of the weights 
  int ng; 	// Number of neighbours in the fluid.
  		// After all GP have been populated, ng=0.
  int ghostTag; // We store here the tag of the surrounding nodes. All the tags of the neighbours 
  		// should be the same. In the case of a complex multiphase flow simulation with Fluid/Structure 
  		// Interaction, this might be no longer true. To be done...
//  ~GhostPoint();

//=============================================================================

  GhostPoint(VarFcn *vf) : varFcn(vf) {
    ng = 0;
    Vg = new double[dim];
    V  = new double[dim];
    Ws = new double[dim];
    for(int i=0;i<dim;++i) {
      Vg[i] = 0.0;
      V[i]  = 0.0;
      Ws[i] = 0.0;
    }
    ghostTag = -2; // Inactive nodes tag
  }
//=============================================================================

  GhostPoint<dim> & operator=(const GhostPoint<dim> &GP) {
    varFcn = GP.varFcn;
    Vg = GP.Vg;
    V  = GP.V;
    Ws = GP.Ws;
    ng = GP.ng;
    ghostTag = GP.ghostTag;
    return *this;
  }
//=============================================================================

  GhostPoint<dim> & operator+=(const GhostPoint<dim> &GP) {
    if(ghostTag<0) ghostTag = GP.ghostTag;
    else if(ghostTag != GP.ghostTag) 
      {
        fprintf(stderr,"The two ghost States refer to different Fluids\n");
        fprintf(stderr,"ghostTag: %i, GP.ghostTag: %i",ghostTag,GP.ghostTag);
        exit(-1);
      }
    Vg += GP.Vg;
    Ws += GP.Ws;
    ng += GP.ng;
    return *this;
  }
//=============================================================================

  void addNeighbour(double *Vi, double *Wi, int tag) {

// We want to satisfy interface condition in least squares manner 
    for(int i=0;i<dim;++i) {
      Vg[i] += Wi[i]*Vi[i];
      Ws[i] += Wi[i];
    }

//    if(dim == 6) { // One Equation Turbulent Model
//      Vg[5] = 0.0;
//      Ws[5] = 1.0;
//    }
//    else if(dim == 7) {// Two Equations Turbulent Model
//      Vg[5] = 0.0;
//      Vg[6] = 0.0;
//      Ws[5] = 1.0;
//      Ws[6] = 1.0;
//    }

    // Tag check
    if(ghostTag < 0) {
      ghostTag = tag;
    }
    else if(ghostTag != tag) {
      fprintf(stderr,"We have a ghost node here with two active neighbours having different tags\n");
      fprintf(stderr,"ghostTag: %i, neighbourTag: %i",ghostTag,tag);
      exit(-1);
    }

    ng++;

  }
//=============================================================================

  double* getState()
  {
    return Vg;    
  }
//=============================================================================

  double* getPrimitiveState()
  {
    return V;    
  }
//=============================================================================

  void reduce()
  {
    for(int i=0;i<dim;++i) {
      Vg[i] /= Ws[i];
      Ws[i] = 1.0;
    }
    ng = 0;

// populate primitive state vector
    for (int i=0; i<dim; ++i) V[i] = Vg[i];
    varFcn->getV4FromTemperature(V,Vg[4],ghostTag);
  }
//=============================================================================

  ~GhostPoint() {
    delete [] Vg;
    delete [] V;
    delete [] Ws;
  }
//=============================================================================
};

#endif
