#include <DistVolumicForceTerm.h>

#include <stdio.h>

//------------------------------------------------------------------------------

template<int dim>
DistVolumicForceTerm<dim>::DistVolumicForceTerm(IoData& iod, Domain* domain)
{

  numLocSub  = domain->getNumLocSub();
  subDomain  = domain->getSubDomain();
  it0        = iod.restart.iteration;
  lastIt     = it0;
  lastConfig = -1;

  subVolumicForceTerm = new VolumicForceTerm<dim>*[numLocSub];
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subVolumicForceTerm[iSub] = new VolumicForceTerm<dim>(iod);

}
//------------------------------------------------------------------------------
template<int dim>
DistVolumicForceTerm<dim>::~DistVolumicForceTerm()
{

  if (subVolumicForceTerm) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      if (subVolumicForceTerm[iSub])
        delete subVolumicForceTerm[iSub];
    delete [] subVolumicForceTerm;
  }
}

//------------------------------------------------------------------------------

template<int dim>
void DistVolumicForceTerm<dim>::compute(int config, DistSVec<double,3> &X)
{

  if ((config != lastConfig) || (lastIt == it0)) {
    for (int iSub = 0; iSub < numLocSub; ++iSub)

    if (lastConfig != config)
      lastConfig = config;
    if (lastIt == it0)
      lastIt = -1;
  }
}

//------------------------------------------------------------------------------
