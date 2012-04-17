#ifndef _NAVIERSTOKESSEGSOLVER_H_
#define _NAVIERSTOKESSEGSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ImplicitSegTsDesc.h>

template<int dim, int neq1, int neq2>
void startNavierStokesSegSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  domain.createVecPat(dim);
  domain.createRhsPat(dim, ioData);


  ImplicitSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);

  TsSolver<ImplicitSegTsDesc<dim,neq1,neq2> > tsSolver(&tsDesc);

  tsSolver.solve(ioData);

}

#endif
