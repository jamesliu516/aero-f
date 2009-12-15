#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitTsDesc.h>
#include <ImplicitCoupledTsDesc.h>
#include <ImplicitSegTsDesc.h>
#include <ImplicitLevelSetTsDesc.h>
#include <ExplicitLevelSetTsDesc.h>
#include <SparseGridGeneratorDesc.h>

// Included (MB)
#include <FluidSensitivityAnalysisHandler.h>

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void startNavierStokesSegSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  domain.createVecPat(dim);
  domain.createRhsPat(dim, ioData);


  ImplicitSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);

  TsSolver<ImplicitSegTsDesc<dim,neq1,neq2> > tsSolver(&tsDesc);

  tsSolver.solve(ioData);

}

//------------------------------------------------------------------------------

template<int dim>
void startLevelSetSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.ts.type == TsData::IMPLICIT) {
    ImplicitLevelSetTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    TsSolver<ImplicitLevelSetTsDesc<dim> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }
  else{
    ExplicitLevelSetTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitLevelSetTsDesc<dim> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }

}

//-----------------------------------------------------------------------------

void startNavierStokesSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator* com = domain.getCommunicator();

  if (ioData.eqs.numPhase == 1){
    if (ioData.eqs.type == EquationsData::EULER)
      startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
    else if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
      if (ioData.eqs.tc.type == TurbulenceClosureData::NONE ||
	  ioData.eqs.tc.type == TurbulenceClosureData::LES)
	startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
      else if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
	if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	    ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
	  if (ioData.ts.type == TsData::IMPLICIT &&
              ioData.ts.implicit.coupling == ImplicitData::WEAK)
	    startNavierStokesSegSolver<6,5,1>(ioData, geoSource, domain);
	  else
	    startNavierStokesCoupledSolver<6>(ioData, geoSource, domain);
	}
	else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
	  if (ioData.ts.implicit.coupling == ImplicitData::WEAK)
	    startNavierStokesSegSolver<7,5,2>(ioData, geoSource, domain);
	  else
	    startNavierStokesCoupledSolver<7>(ioData, geoSource, domain);
	}
	else {
	  com->fprintf(stderr, "*** Error: wrong turbulence model type\n");
	  exit(1);
	}
      }
      else {
	com->fprintf(stderr, "*** Error: wrong turbulence closure type\n");
	exit(1);
      }
    }
    else {
      com->fprintf(stderr, "*** Error: wrong equation type\n");
      exit(1);
    }
  }
  else if (ioData.eqs.numPhase == 2){
    startLevelSetSolver<5>(ioData, geoSource, domain);
  }else
    com->fprintf(stderr, "*** Error: wrong number of phases\n");

}

//------------------------------------------------------------------------------

void startSparseGridGeneration(IoData &ioData, Domain &domain)
{

  Communicator* com = domain.getCommunicator();

  fprintf(stdout, "*** Warning: Generating a sparse grid\n");
  SparseGridGeneratorDesc sgDesc(ioData, com);
  sgDesc.tabulate(ioData);

}

//------------------------------------------------------------------------------
