#include <IoData.h>
#include <Domain.h>
#include "Solvers/Solvers.h"
#include <SparseGridGeneratorDesc.h>

//-----------------------------------------------------------------------------

void startNavierStokesSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();
  if (ioData.strucIntersect.intersectorName != 0) {
    com->fprintf(stderr, "*** NOTE: Running an Embedded Fluid-Structure simulation\n");
    StructLevelSetSolver<5>::solve(ioData, geoSource, domain);
  } else if (ioData.eqs.numPhase == 1){
    if (ioData.eqs.type == EquationsData::EULER) 
      NavierStokesCoupledSolver<5>::solve(ioData, geoSource, domain);
    else if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
      if (ioData.eqs.tc.type == TurbulenceClosureData::NONE ||
	  ioData.eqs.tc.type == TurbulenceClosureData::LES)
	//startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
	NavierStokesCoupledSolver<5>::solve(ioData, geoSource, domain);
      else if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
	if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	    ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
	  if (ioData.ts.type == TsData::IMPLICIT &&
              ioData.ts.implicit.coupling == ImplicitData::WEAK)
	   // startNavierStokesSegSolver<6,5,1>(ioData, geoSource, domain);
	    NavierStokesSegSolver<6,5,1>::solve(ioData, geoSource, domain);
	  else
	    NavierStokesCoupledSolver<6>::solve(ioData, geoSource, domain);
	}
	else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
	  if (ioData.ts.implicit.coupling == ImplicitData::WEAK)
	    NavierStokesSegSolver<7,5,2>::solve(ioData, geoSource, domain);
	  else
	    NavierStokesCoupledSolver<7>::solve(ioData, geoSource, domain);
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
    com->fprintf(stdout, "*** Warning: number of phases is %d\n", ioData.eqs.numPhase);
    LevelSetSolver<5,1>::solve(ioData, geoSource, domain);
  }
  else if (ioData.eqs.numPhase == 3){
    com->fprintf(stdout, "*** Warning: number of phases is %d\n", ioData.eqs.numPhase);
    LevelSetSolver<5,2>::solve(ioData, geoSource, domain);
  }

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
