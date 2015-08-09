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
#include <ImplicitCollocationTsDesc.h>
#include <ImplicitMetricTsDesc.h>
#include <ImplicitRomPostproTsDesc.h>
#include <MultiGridSolver.h>
#include <MultiGridCoupledTsDesc.h>
#include <FluidShapeOptimizationHandler.h>
#include <FluidRomShapeOptimizationHandler.h>  // MZ

template<int dim>
void startNavierStokesCoupledSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_) { // YC
      FluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
  }
  else if (ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_) { // YC
      FluidShapeOptimizationHandler<dim> fsisoh(ioData, geoSource, &domain);
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsisoh);
      tsSolver.fsisoSolve(ioData);
  }
  else if (ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_) { // MZ
      FluidRomShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidRomShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
  }
  else if ((ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_) || 
           (ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_) ||
           (ioData.problem.alltype == ProblemData::_ACC_UNSTEADY_NONLINEAR_ROM_) ||
           (ioData.problem.alltype == ProblemData::_FORCED_NONLINEAR_ROM_)) {
    if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) { 
        ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitPGTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::GNAT) {
        ImplicitGnatTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitGnatTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION) {
        ImplicitCollocationTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitCollocationTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
        ImplicitMetricTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitMetricTsDesc<dim> > tsSolver(&tsDesc);
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
