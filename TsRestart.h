#ifndef _TS_RESTART_H_
#define _TS_RESTART_H_

class IoData;
class RefVal;
class DistGeoState;
class LevelSet;

template<int dim> class DistTimeState;

//------------------------------------------------------------------------------

class TsRestart {

  RefVal *refVal;

  int index;

public:

  int iteration;
  double etime;
  double residual;
  double energy[2];

  char *solutions[3];
  char *positions[3];
  char *levelsets[3];
  char *data[3];

  int frequency;

public:

  TsRestart(IoData &, RefVal *);
  TsRestart();

  template<int dim>
  void writeToDisk(int, bool, int, double, double, 
		   DistTimeState<dim> &, DistGeoState &, 
		   LevelSet *LS);

// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TsRestart.C>
#endif

#endif
