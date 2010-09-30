#include "Solvers.h"
#include "NavierStokesMultiPhysicsEmbedded.h"

template <>
void
NavierStokesMultiPhysicsEmbedded<5,1>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesMultiPhysicsEmbedded<5,1>(ioData, geoSource, domain);
}
