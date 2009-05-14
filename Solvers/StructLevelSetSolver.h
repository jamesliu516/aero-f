#ifndef _LEVELSETSOLVER_H_
#define _LEVELSETSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitStructLevelSetTsDesc.h>
template<int dim>
void startStructLevelSetSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.ts.type == TsData::IMPLICIT) {
    com->fprintf(stderr, "***Error: wrong time integrator for EulerStructGhostFluid method\n");   
  }
  else{
    ExplicitStructLevelSetTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitStructLevelSetTsDesc<dim> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }

}

#endif
