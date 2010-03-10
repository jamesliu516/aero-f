/***********************************************************************************
*                                                                                  *
* Different explicit time-integration schemes are available:                       *
*  - Forward Euler where both fluid states and level set are advanced together     *
*  - RungeKutta2   where both fluid states and level set are advanced together     *
*                                                          (cf LeTallec's comment) *
*  - RungeKutta2   where fluid states and level set are advanced in a staggered    *
*                  fashion (original idea of Pr Farhat to increase stability,      *
*                  but then accuracy is lost)                                      *
*  - RungeKutta4   similar to staggered RungeKutta2                                *
*                                                                                  *
*                                                                                  *
*                                                                                  *
*                                                                                  *
***********************************************************************************/

#ifndef _EXPLICIT_LEVELSET_TS_DESC_H_
#define _EXPLICIT_LEVELSET_TS_DESC_H_

#include <LevelSetTsDesc.h>

#include <IoData.h>
#include <Domain.h>

struct DistInfo;

class GeoSource;
class LevelSet; 
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------

template<int dim>
class ExplicitLevelSetTsDesc : public LevelSetTsDesc<dim> {

 private:
  ExplicitData::Type timeType;

  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;

  DistSVec<double,1> Phi0;
  DistVec<int> fluidId0;
  DistSVec<double,1> p1;
  DistSVec<double,1> p2;
  DistSVec<double,1> p3;
  DistSVec<double,1> p4;

// mesh motion modification for RK2
// otherwise equal to U and Phi respectively
  DistSVec<double,dim> ratioTimesU;
  DistSVec<double,1> ratioTimesPhi;

 public:
  ExplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitLevelSetTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &U);

 private:
  void solveNLSystemOneBlock(DistSVec<double,dim> &U);
  void solveNLSystemTwoBlocks(DistSVec<double,dim> &U);

// for solving the total system in one block (U and Phi at the same time)
  void solveNLAllFE(DistSVec<double,dim> &U);
  void solveNLAllRK2(DistSVec<double,dim> &U);
  void solveNLAllRK2bis(DistSVec<double,dim> &U);

// for solving the total system in two blocks (U first, Phi second)
  void solveNLEuler(DistSVec<double,dim> & U);
  void solveNLEulerRK2(DistSVec<double,dim> &U);
  void solveNLEulerRK4(DistSVec<double,dim> &U);
  void solveNLLevelSet(DistSVec<double,dim> &U);
  void solveNLLevelSetRK2(DistSVec<double,dim> &U);
  void solveNLLevelSetRK4(DistSVec<double,dim> &U);


  void computeRKUpdate(DistSVec<double,dim>& Ulocal, 
                       DistSVec<double,dim>& dU, int it);
  void computeRKUpdateLS(DistSVec<double,1>& Philocal, DistVec<int> &localFluidId, DistSVec<double,1>& dPhi, 
                         DistSVec<double,dim>& U);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitLevelSetTsDesc.C>
#endif

#endif
