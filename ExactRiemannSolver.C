#include <ExactRiemannSolver.h>
#include <IoData.h>
#include <Vector.h>
#include <Vector3D.h>
#include <LocalRiemannDesc.h>

#include <math.h>


//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::ExactRiemannSolver(IoData &iod, SVec<double,dim> &_rupdate,
                                            Vec<double> &_weight):
                                            rupdate(_rupdate), weight(_weight)
{

  iteration = -1;
  lriemann = 0;
  if(iod.mf.method == MultiFluidData::GHOSTFLUID_FOR_POOR){
    if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
       iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmpGasGas();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmpTaitTait();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmpJWLJWL();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmpGasJWL();
    else{
      fprintf(stdout, "*** Error: no gfmp possible for that simulation\n");
      exit(1);
    }
  }else if(iod.mf.method == MultiFluidData::GHOSTFLUID_WITH_RIEMANN){
    if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
       iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmparGasGas();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::GAS)
      lriemann = new LocalRiemannGfmparGasTait();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::LIQUID &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmparTaitTait();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::LIQUID)
      lriemann = new LocalRiemannGfmparGasTait();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::GAS &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmparGasJWL();
    else if(iod.eqs.fluidModel.fluid  == FluidModelData::JWL &&
            iod.eqs.fluidModel2.fluid == FluidModelData::JWL)
      lriemann = new LocalRiemannGfmparJWLJWL();
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

}

//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeRiemannSolution(double *Vi, double *Vj,
      double Phii, double Phij, double *nphi, VarFcn *vf,
      int &epsi, int &epsj, double *Wi, double *Wj, int i, int j)
{

  lriemann->computeRiemannSolution(Vi,Vj,Phii,Phij,nphi,vf,
          epsi,epsj,Wi,Wj,rupdate[i],rupdate[j],weight[i],weight[j], iteration);

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
