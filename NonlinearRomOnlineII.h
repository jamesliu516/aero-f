#ifndef _NONLINEAR_ROM_ONLINE_II_H_
#define _NONLINEAR_ROM_ONLINE_II_H_

#include <NonlinearRom.h>

template <int dim>
class NonlinearRomOnlineII : public NonlinearRom<dim> {

  protected:


  public:

  NonlinearRomOnlineII(Communicator *, IoData &, Domain &);
  ~NonlinearRomOnlineII();

  void closestCenter(DistSVec<double, dim> &, int*);
  void updateBasis(int, DistSVec<double, dim> &);
  void appendNonStateDataToBasis(int, char *); 
  void readClusterOnlineQuantities(int);
  void readDistanceCalcInfo();

  //void appendVectorToBasis(DistSVec<double, dim>*, int numVec = 0);
 
};

#include "NonlinearRomOnlineII.C"
#endif
