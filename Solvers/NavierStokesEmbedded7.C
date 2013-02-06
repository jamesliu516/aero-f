#include "Solvers.h"
#include "NavierStokesEmbedded.h"

template <>
void
NavierStokesEmbedded<7>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbedded<7>(ioData, geoSource, domain);
}
