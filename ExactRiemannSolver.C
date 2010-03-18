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
  fsiRiemann = 0;

// FSI Riemann problem
  if(1) // FIX
  //if(iod.strucIntersect.intersectorName != 0) // FIX
    fsiRiemann = new LocalRiemannFluidStructure();

// Multiphase Riemann problem
  if(iod.mf.method == MultiFluidData::GHOSTFLUID_FOR_POOR){
    if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
       iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmpGasGas(vf);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmpTaitTait(vf);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmpJWLJWL(vf);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmpGasJWL(vf);
    else{
      //fprintf(stdout, "*** Error: no gfmp possible for that simulation\n");
      //exit(1);
    }
  }else if(iod.mf.method == MultiFluidData::GHOSTFLUID_WITH_RIEMANN){
    if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
       iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmparGasGas(vf, iod.mf.typePhaseChange);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmparGasTait(vf, iod.mf.typePhaseChange);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmparTaitTait(vf, iod.mf.typePhaseChange);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmparGasTait(vf, iod.mf.typePhaseChange);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmparGasJWL(vf,sgCluster,iod.mf.riemannComputation,iod.mf.typePhaseChange);
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmparJWLJWL(vf, iod.mf.typePhaseChange);
    else{
      fprintf(stdout, "*** Error: no gfmpar possible for that simulation\n");
      exit(1);
    }
  } 
		

}

//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::~ExactRiemannSolver() 
{

  if(lriemann) delete lriemann;
  if(fsiRiemann) delete fsiRiemann;

}

//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::updatePhaseChange(SVec<double,dim> &V, Vec<int> &fluidId,
                                                Vec<int> &fluidIdn)
{

  for(int i=0; i<V.size(); i++){
    lriemann->updatePhaseChange(V[i],fluidId[i],fluidIdn[i],rupdate[i],weight[i]);
  }

}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeRiemannSolution(double *Vi, double *Vj,
      int IDi, int IDj, double *nphi, VarFcn *vf,
      double *Wi, double *Wj, int i, int j, int edgeNum,
      double dx[3])
{

  lriemann->computeRiemannSolution(Vi,Vj,IDi,IDj,nphi,interfacialWi[edgeNum],interfacialWj[edgeNum],
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
  fsiRiemann->computeRiemannSolution(tag, Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration);
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
