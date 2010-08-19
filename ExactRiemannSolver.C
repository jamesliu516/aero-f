#include <ExactRiemannSolver.h>
#include <IoData.h>
#include <Vector.h>
#include <Vector3D.h>
#include <LocalRiemannDesc.h>

#include <math.h>


//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::ExactRiemannSolver(IoData &iod, SVec<double,dim> &_rupdate,
                                            Vec<double> &_weight, SVec<double,dim-2> &_interfacialWi,
                                            SVec<double,dim-2> &_interfacialWj, VarFcn *vf,
                                            SparseGridCluster *sgCluster):
                                            rupdate(_rupdate), weight(_weight),
                                            interfacialWi(_interfacialWi),
                                            interfacialWj(_interfacialWj)
{

  iteration = -1;
  lriemann = 0;
  numLriemann = 0;
  fsiRiemann = 0;

// FSI Riemann problem
  if(iod.problem.framework==ProblemData::EMBEDDED) 
    fsiRiemann = new LocalRiemannFluidStructure<dim>(); //NOTE(KW): The following lines will still be 
                                                   //  executed. Currently they are never used but in
                                                   //  future if we have both FS and FF, they are needed.

// Multiphase Riemann problem
// Assumption: there are only numPhase-1 interfaces (between fluid 0 and each one of the other fluids)
  if(iod.eqs.numPhase > 1){
    numLriemann = iod.eqs.numPhase-1;
    lriemann = new LocalRiemann*[numLriemann];
    for(int iPhase=0; iPhase<numLriemann; iPhase++){
      map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(iPhase+1);
      if(it == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", iPhase+1);
        exit(1);
      }
      if(iod.mf.method == MultiFluidData::GHOSTFLUID_FOR_POOR){
        if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
           it->second->fluid == FluidModelData::GAS)
          lriemann[iPhase] = new LocalRiemannGfmpGasGas(vf,0,iPhase+1);
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
                it->second->fluid == FluidModelData::LIQUID)
          lriemann[iPhase] = new LocalRiemannGfmpTaitTait(vf,0,iPhase+1);
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
                it->second->fluid == FluidModelData::JWL)
          lriemann[iPhase] = new LocalRiemannGfmpJWLJWL(vf,0,iPhase+1);
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
                it->second->fluid == FluidModelData::JWL)
          lriemann[iPhase] = new LocalRiemannGfmpGasJWL(vf,0,iPhase+1);
        else{
          fprintf(stdout, "*** Error: no gfmp possible for that simulation\n");
          exit(1);
        }
      }else if(iod.mf.method == MultiFluidData::GHOSTFLUID_WITH_RIEMANN){
        if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
           it->second->fluid == FluidModelData::GAS){
          lriemann[iPhase] = new LocalRiemannGfmparGasGas(vf,0,iPhase+1, iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparGasGas\n", iPhase);
        }
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
                it->second->fluid == FluidModelData::GAS){
          lriemann[iPhase] = new LocalRiemannGfmparGasTait(vf,iPhase+1,0, iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparGasTait\n", iPhase);
        }
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
                it->second->fluid == FluidModelData::LIQUID){
          lriemann[iPhase] = new LocalRiemannGfmparTaitTait(vf,0,iPhase+1, iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparTaitTait\n", iPhase);
        }
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
                it->second->fluid == FluidModelData::LIQUID){
          lriemann[iPhase] = new LocalRiemannGfmparGasTait(vf,0,iPhase+1, iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparGasTait\n", iPhase);
        }
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
                it->second->fluid == FluidModelData::JWL){
          lriemann[iPhase] = new LocalRiemannGfmparGasJWL(vf,0,iPhase+1,sgCluster,iod.mf.riemannComputation,iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparGasJwl\n", iPhase);
        }
        else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
                it->second->fluid == FluidModelData::JWL){
          lriemann[iPhase] = new LocalRiemannGfmparJWLJWL(vf,0,iPhase+1, iod.mf.typePhaseChange);
          //fprintf(stdout, "Debug: created %d - LocalRiemannGfmparJwlJwl\n", iPhase);
        }
        else{
          fprintf(stdout, "*** Error: no gfmpar possible for that simulation\n");
          exit(1);
        }
      } 
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::~ExactRiemannSolver() 
{

  for(int iLriemann=0; iLriemann<numLriemann; iLriemann++)
    delete lriemann[iLriemann];
  delete [] lriemann;
  delete fsiRiemann;

}

//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::updatePhaseChange(SVec<double,dim> &V, Vec<int> &fluidId,
                                                Vec<int> &fluidIdn)
{

  for(int i=0; i<V.size(); i++){
    lriemann[0]->updatePhaseChange(V[i],fluidId[i],fluidIdn[i],rupdate[i],weight[i]);
    // lriemann[0] can be used for all interfaces, because  this routine does
    // not need to consider which interface has traversed that node (this
    // was done previously when computing the riemann problem)
  }

}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeRiemannSolution(double *Vi, double *Vj,
      int IDi, int IDj, double *nphi, VarFcn *vf,
      double *Wi, double *Wj, int i, int j, int edgeNum,
      double dx[3])
{

  //fprintf(stdout, "Debug: calling computeRiemannSolution with IDi = %d - IDj = %d for LocalRiemann[%d]\n", IDi, IDj, IDi+IDj-1);
  lriemann[IDi+IDj-1]->computeRiemannSolution(Vi,Vj,IDi,IDj,nphi,interfacialWi[edgeNum],interfacialWj[edgeNum],
          Wi,Wj,rupdate[i],rupdate[j],weight[i],weight[j],
          dx,iteration);

}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeFSIRiemannSolution(double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej, int Id)

{
  fsiRiemann->computeRiemannSolution(Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration, Id);
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeFSIRiemannSolution(int tag, double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej)

{
  // Adam 2010.08.18
  // This function doesn't seem to be used anymore.
  // To be removed in a couple of months
  fprintf(stderr,"Oh Sorry ! Please uncomment the function (ExactRiemannSolver.C:159). I thought it wasn't needed anymore\n");
  exit(-1);
  /*
  fsiRiemann->computeRiemannSolution(tag, Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration);
  */
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::reset(int it)
{

  iteration = it;
  if(iteration==1){
    rupdate = 0.0;
    weight  = 0.0;
  }

}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::resetInterfacialW(int edgeNum)
{

  //dim-2 since interfacialW store only rho, u, p
  for(int idim=0; idim<dim-2; idim++){
    interfacialWi[edgeNum][idim] = 0.0;
    interfacialWj[edgeNum][idim] = 0.0;
  }

}
//------------------------------------------------------------------------------
