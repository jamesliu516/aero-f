#include <cstdio>

#include <TsRestart.h>

#include <IoData.h>
#include <RefVal.h>
#include <DistTimeState.h>
#include <DistGeoState.h>
#include <MeshMotionHandler.h>
#include <LevelSet.h>

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void TsRestart::writeToDisk(int cpuNum, bool lastIt, int it, double t, double dt,
			    DistTimeState<dim> &timeState, DistGeoState &geoState,
			    LevelSet<dimLS> *levelSet)
{

  iteration = it;
  etime = t;
  double dt_nm1 = timeState.getData().dt_nm1;
  double dt_nm2 = timeState.getData().dt_nm2;

  if (toWrite(iteration, lastIt, etime)) {

    if (lastIt) 
      index = 0;

    if (solutions[index][0] != 0)
      timeState.writeToDisk(solutions[index]);
    if (positions[index][0] != 0)
      geoState.writeToDisk(positions[index]);
    if (levelsets[index][0] != 0 && levelSet)
      levelSet->writeToDisk(levelsets[index]);
  
    if (cpuNum == 0 && data[index][0] != 0) {
      FILE *fp = fopen(data[index], "w");
      if (!fp) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", data[index]);
	exit(1);
      }
      
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
	fprintf(fp, "// Restart file (values are non dimensional)\n\n");
      else
	fprintf(fp, "// Restart file (values are dimensional)\n\n");
      fprintf(fp, "under RestartParameters {\n");
      fprintf(fp, "  Iteration = %d;\n", iteration);
      fprintf(fp, "  Time = %e;\n", etime * refVal->time);
      fprintf(fp, "  TimeStep1 = %e;\n", dt_nm1 * refVal->time);
      fprintf(fp, "  TimeStep2 = %e;\n", dt_nm2 * refVal->time);
      fprintf(fp, "  Residual = %e;\n", residual);
      fprintf(fp, "  Energy = %e;\n", energy[0] * refVal->energy);
      fprintf(fp, "}\n");
      fclose(fp);
    }

    if (index == 1)
      index = 2;
    else if (index == 2)
      index = 1;

  }

}

//------------------------------------------------------------------------------


