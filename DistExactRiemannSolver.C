#include <DistExactRiemannSolver.h>


#include <ExactRiemannSolver.h>
#include <IoData.h>
#include <Domain.h>
#include <Communicator.h>
#include "SparseGridCluster.h"

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

  if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){

    double *refIn  = new double[2];
    double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = pow(ioData.ref.rv.density,-ioData.eqs.fluidModel2.jwlModel.omega)*ioData.ref.rv.velocity*ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.velocity;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(1,refIn, refOut, "SparseGridClusterRiemannInvariant", dom->getCommunicator()->cpuNum());

  }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){

    double *refIn = new double[5]; double *refOut = new double[2];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = ioData.ref.rv.pressure;
    refIn[2] = ioData.ref.rv.density;
    refIn[3] = ioData.ref.rv.pressure;
    refIn[4] = ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.density;
    refOut[1] = ioData.ref.rv.density;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(62,refIn, refOut, "SparseGridClusterRiemannSolution", dom->getCommunicator()->cpuNum());
  }else{
    tabulationC = 0;
  }

  subExactRiemannSolver = new ExactRiemannSolver<dim>*[numLocSub];
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subExactRiemannSolver[iSub] =
      new ExactRiemannSolver<dim>(ioData, (*riemannupdate)(iSub),
                         	(*weight)(iSub), vf, tabulationC);

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

  delete tabulationC;

}
//------------------------------------------------------------------------------
