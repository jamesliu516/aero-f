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
  double frequency_dt, prtout;

  bool deleteCharStar;

private:
  bool toWrite(int it, bool lastIt, double t);

public:

  TsRestart(IoData &, RefVal *);
  TsRestart();

  void updatePrtout(double t);

  template<int dim, int dimLS>
  void writeToDisk(int, bool, int, double, double, 
		   DistTimeState<dim> &, DistGeoState &, LevelSet<dimLS> *levelSet = 0);

  /** Function to write the structure positions to disk. Used for the embedded-method only. */
  void writeStructPosToDisk(int, bool, Vec<Vec3D>&);
 
// Included (MB)
  void rstVar(IoData &);
  ~TsRestart()
    {
      int last=deleteCharStar ? 3 : 1;
      for(int i=0;i<last;++i)
	{
	  delete[] data[i];
	  delete[] solutions[i];
	  delete[] positions[i];
	  delete[] levelsets[i];
	}
      delete[] structPos;
    }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TsRestart.C>
#endif

#endif
