#ifndef _NAVIERSTOKESCOUPLEDSOLVER_H_
#define _NAVIERSTOKESCOUPLEDSOLVER_H_

#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitTsDesc.h>
#include <ImplicitCoupledTsDesc.h>
#include <ImplicitPGTsDesc.h>
#include <ImplicitGalerkinTsDesc.h>
#include <ImplicitBroydenTsDesc.h>
#include <ImplicitGnatTsDesc.h>
#include <ImplicitCollLSTsDesc.h>
#include <ImplicitCollGalTsDesc.h>
#include <ImplicitRomPostproTsDesc.h>
#include <ImplicitProjErrorTsDesc.h>
#include <MultiGridSolver.h>
#include <MultiGridCoupledTsDesc.h>
#include <FluidShapeOptimizationHandler.h>  // YC
// Included (MB)
#include <FluidSensitivityAnalysisHandler.h>


template<int dim>
void startNavierStokesCoupledSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);





  if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_)
  {
// Modified (MB)
      FluidSensitivityAnalysisHandler<dim> fsah(ioData, geoSource, &domain);
      TsSolver<FluidSensitivityAnalysisHandler<dim> > tsSolver(&fsah);
      tsSolver.fsaSolve(ioData);
  }
	else if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_) { // YC
      FluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
	}
	else if (ioData.problem.alltype == ProblemData::_NONLINEAR_ROM_ ) {
		if (ioData.rom.projection == 0 && ioData.rom.systemApproximation == 0) {
			ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitPGTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 0 && ioData.rom.systemApproximation == 1) {
			ImplicitGnatTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitGnatTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 0 && ioData.rom.systemApproximation == 2) {
			ImplicitCollLSTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitCollLSTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 0 && ioData.rom.systemApproximation == 3) {
			ImplicitBroydenTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitBroydenTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 1 && ioData.rom.systemApproximation == 0) {
			ImplicitGalerkinTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitGalerkinTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 1 && ioData.rom.systemApproximation == 2) {
			ImplicitCollGalTsDesc<dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitCollGalTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else if (ioData.rom.projection == 2) {
			ImplicitProjErrorTsDesc <dim> tsDesc(ioData, geoSource, &domain);
			TsSolver<ImplicitProjErrorTsDesc<dim> > tsSolver(&tsDesc);
			tsSolver.solve(ioData);
		}
		else
			com->fprintf(stderr, "*** Error: incompatible nonlinear ROM\n");
	}
  else if (ioData.problem.alltype == ProblemData::_NONLINEAR_ROM_POST_) {
		ImplicitRomPostproTsDesc <dim> tsDesc(ioData, geoSource, &domain);
		TsSolver<ImplicitRomPostproTsDesc<dim> > tsSolver(&tsDesc);
		tsSolver.solve(ioData);
	}
  else if (ioData.ts.type == TsData::IMPLICIT) {
    if (ioData.problem.solutionMethod == ProblemData::TIMESTEPPING) {
      ImplicitCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitCoupledTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    } else {
      MultiGridCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      MultiGridSolver<MultiGridCoupledTsDesc<dim> > mgSolver(&tsDesc);
      mgSolver.solve(ioData);
    }
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
