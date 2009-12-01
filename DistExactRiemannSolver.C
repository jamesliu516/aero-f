#include <DistExactRiemannSolver.h>


#include <ExactRiemannSolver.h>
#include <IoData.h>
#include <Domain.h>

#include <stdio.h>

//------------------------------------------------------------------------------

template<int dim>
DistExactRiemannSolver<dim>::DistExactRiemannSolver(IoData &ioData, Domain *dom,
                                                    VarFcn *vf)
{

  updatePhase = false;
  if(ioData.mf.method          == MultiFluidData::GHOSTFLUID_WITH_RIEMANN &&
     ioData.mf.typePhaseChange == MultiFluidData::RIEMANN_SOLUTION)
    updatePhase = true;

  numLocSub = dom->getNumLocSub();
  firstpass= true;

  riemannupdate = new DistSVec<double,dim>(dom->getNodeDistInfo());
  weight        = new DistVec<double>(dom->getNodeDistInfo());

  subExactRiemannSolver = new ExactRiemannSolver<dim>*[numLocSub];
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subExactRiemannSolver[iSub] =
      new ExactRiemannSolver<dim>(ioData, (*riemannupdate)(iSub),
                         	(*weight)(iSub), vf);

}
//------------------------------------------------------------------------------
template<int dim>
DistExactRiemannSolver<dim>::~DistExactRiemannSolver()
{

  if (riemannupdate) delete riemannupdate;
  if (weight) delete weight;

  if (subExactRiemannSolver) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      if (subExactRiemannSolver[iSub]) 
        delete subExactRiemannSolver[iSub];

    delete [] subExactRiemannSolver;
  }

}
//------------------------------------------------------------------------------
