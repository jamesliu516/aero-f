#include "Solvers.h"
#include "StructLevelSetSolver.h"

template <>
void
LevelSetSolver<5>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startStructLevelSetSolver<5>(ioData, geoSource, domain);
}
