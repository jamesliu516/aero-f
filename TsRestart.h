#ifndef _TS_RESTART_H_
#define _TS_RESTART_H_

#include<Vector.h>
#include<Vector3D.h>

class IoData;
class RefVal;
class DistGeoState;

template<int dimLS> class LevelSet;
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

  char *structPos;

  int frequency;

public:

  TsRestart(IoData &, RefVal *);
  TsRestart();

  template<int dim, int dimLS>
  void writeToDisk(int, bool, int, double, double, 
		   DistTimeState<dim> &, DistGeoState &, LevelSet<dimLS> *levelSet = 0);

  /** Function to write the structure positions to disk. Used for the embedded-method only. */
  void writeStructPosToDisk(int, bool, Vec<Vec3D>&);
 
// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TsRestart.C>
#endif

#endif
