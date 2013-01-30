#ifndef _NAVIERSTOKESSEGSOLVER_H_
#define _NAVIERSTOKESSEGSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ImplicitSegTsDesc.h>
#include <MultiGridSegTsDesc.h>
#include <MultiGridSolver.h>

template<int dim, int neq1, int neq2>
void startNavierStokesSegSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  domain.createVecPat(dim);
  domain.createRhsPat(dim, ioData);

  if (ioData.problem.solutionMethod == ProblemData::TIMESTEPPING) {
    ImplicitSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);

    TsSolver<ImplicitSegTsDesc<dim,neq1,neq2> > tsSolver(&tsDesc);

    tsSolver.solve(ioData);
  } else {
    
    MultiGridSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);
    MultiGridSolver<MultiGridSegTsDesc<dim,neq1,neq2> > mgSolver(&tsDesc);
    mgSolver.solve(ioData);

  }
}

#endif
