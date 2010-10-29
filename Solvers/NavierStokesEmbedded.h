#ifndef _NAVIER_STOKES_EMBEDDED_H_
#define _NAVIER_STOKES_EMBEDDED_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitEmbeddedTsDesc.h>
template<int dim>
void startNavierStokesEmbedded(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.ts.type == TsData::IMPLICIT) {
    com->fprintf(stderr, "***Error: wrong time integrator for EulerStructGhostFluid method\n");   
  }
  else{
    ExplicitEmbeddedTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitEmbeddedTsDesc<dim> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }

}

#endif
