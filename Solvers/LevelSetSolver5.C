#include "Solvers.h"
#include "LevelSetSolver.h"

template <>
void
LevelSetSolver<5>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startLevelSetSolver<5>(ioData, geoSource, domain);
}
