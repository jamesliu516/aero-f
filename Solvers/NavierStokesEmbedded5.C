#include "Solvers.h"
#include "NavierStokesEmbedded.h"

template <>
void
NavierStokesEmbedded<5>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbedded<5>(ioData, geoSource, domain);
}
