//
//       Class LevelSet   created by Francois Courty (C.A.S. Boulder)
//       March 2004
//
///////////////////////////////////////////////////////////////////////

#ifndef _LEVEL_SET_H_
#define _LEVEL_SET_H_

#include <IoData.h>
#include <DistInfo.h>
#include <DistVector.h>

#include <Domain.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Vector.h>
#include <Communicator.h>

class LevelSet {

  Communicator *com;

 public:

  double spheres[10][4];
  int nspheres;

  DistVec<double> Phi;
  DistVec<double> Phi1;
  DistVec<double> Phi2;
  DistVec<double> PhiRet;

  LevelSet(DistSVec<double,3> & X, IoData &iod, Domain *domain);
  ~LevelSet(); 

  DistVec<double> &getPhi();

  template<int dim>
  void conservativeToPrimitive(DistSVec<double,dim> &U);

  template<int dim>
  void primitiveToConservative(DistSVec<double,dim> &U);
};

#ifdef TEMPLATE_FIX
#include <LevelSet.C>
#endif

#endif


