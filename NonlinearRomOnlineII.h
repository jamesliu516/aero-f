#ifndef _NONLINEAR_ROM_ONLINE_II_H_
#define _NONLINEAR_ROM_ONLINE_II_H_

#include <NonlinearRom.h>

template <int dim>
class NonlinearRomOnlineII : public NonlinearRom<dim> {

  protected:


  public:

  NonlinearRomOnlineII(Communicator *, IoData &, Domain &);
  ~NonlinearRomOnlineII();

  void updateBasis(int, DistSVec<double, dim> &);
  void appendNonStateDataToBasis(int, char *, bool relProjError = false); 
  void readClusteredOnlineQuantities(int);

  void readClosestCenterInfoModelII();

  //void appendVectorToBasis(DistSVec<double, dim>*, int numVec = 0);
 
};

#include "NonlinearRomOnlineII.C"
#endif
