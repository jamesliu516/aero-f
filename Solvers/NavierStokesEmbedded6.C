#include "Solvers.h"
#include "NavierStokesEmbedded.h"

template <>
void
NavierStokesEmbedded<6>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbedded<6>(ioData, geoSource, domain);
}
