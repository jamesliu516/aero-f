#ifndef _NAVIERSTOKESCOUPLEDSOLVER_H_
#define _NAVIERSTOKESCOUPLEDSOLVER_H_

#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitTsDesc.h>
#include <ImplicitCoupledTsDesc.h>
// Included (MB)
#include <FluidSensitivityAnalysisHandler.h>

template<int dim>
void startNavierStokesCoupledSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

// Modified (MB)
    if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {
      FluidSensitivityAnalysisHandler<dim> fsah(ioData, geoSource, &domain);
      TsSolver<FluidSensitivityAnalysisHandler<dim> > tsSolver(&fsah);
      tsSolver.fsaSolve(ioData);
    }
    else if (ioData.ts.type == TsData::IMPLICIT) {
      ImplicitCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitCoupledTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    }
    else if (ioData.ts.type == TsData::EXPLICIT) {
      ExplicitTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ExplicitTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    }
    else
      com->fprintf(stderr, "*** Error: wrong time-integrator\n");
      
}

#endif
