#ifndef _NAVIERSTOKESCOUPLEDSOLVER_H_
#define _NAVIERSTOKESCOUPLEDSOLVER_H_

#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitTsDesc.h>
#include <ImplicitCoupledTsDesc.h>
#include <ImplicitPGTsDesc.h>
//#include <ImplicitGalerkinTsDesc.h>
#include <ImplicitGnatTsDesc.h>
#include <ImplicitRomPostproTsDesc.h>
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
		else if ((ioData.problem.alltype == ProblemData::_NONLINEAR_ROM_) || (ioData.problem.alltype == ProblemData::_ACC_NONLINEAR_ROM_) ) {
			if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == 0) {
				ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);
				TsSolver<ImplicitPGTsDesc<dim> > tsSolver(&tsDesc);
				tsSolver.solve(ioData);
			}
			else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == 1) {
					ImplicitGnatTsDesc<dim> tsDesc(ioData, geoSource, &domain);
					TsSolver<ImplicitGnatTsDesc<dim> > tsSolver(&tsDesc);
					tsSolver.solve(ioData);
			}
			/*else if (ioData.rom.projection == 1 && ioData.rom.systemApproximation == 0) {
					ImplicitGalerkinTsDesc<dim> tsDesc(ioData, geoSource, &domain);
					TsSolver<ImplicitGalerkinTsDesc<dim> > tsSolver(&tsDesc);
					tsSolver.solve(ioData);
			}*/
			else
				com->fprintf(stderr, "*** Error: this type of nonlinear ROM simulation is not currently supported\n");
		}
    else if (ioData.problem.alltype == ProblemData::_NONLINEAR_ROM_POST_) {
				ImplicitRomPostproTsDesc <dim> tsDesc(ioData, geoSource, &domain);
				TsSolver<ImplicitRomPostproTsDesc<dim> > tsSolver(&tsDesc);
				tsSolver.solve(ioData);
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
