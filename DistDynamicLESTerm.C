#include <DistDynamicLESTerm.h>

#include <DistVector.h>
#include <Domain.h>

#include <stdio.h>

//------------------------------------------------------------------------

template<int dim>
DistDynamicLESTerm<dim>::DistDynamicLESTerm(IoData &iod, Domain *dom) : domain(dom)

{

  numLocSub  = domain->getNumLocSub();

  gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  R = iod.eqs.fluidModel.gasModel.idealGasConstant;

  VCap = new DistSVec<double,dim>(domain->getNodeDistInfo()); 
  Mom_Test = new DistSVec<double,16>(domain->getNodeDistInfo());
  Eng_Test = new DistSVec<double,6>(domain->getNodeDistInfo());

}

//------------------------------------------------------------------------

template<int dim>
DistDynamicLESTerm<dim>::~DistDynamicLESTerm()

{

  if (VCap) delete VCap;
  if (Mom_Test) delete Mom_Test;
  if (Eng_Test) delete Eng_Test;

}

//------------------------------------------------------------------------

template<int dim>
void DistDynamicLESTerm<dim>::computeTestFilterValues(DistSVec<double,2> &Cs, DistVec<double> &VolSum,
			DistSVec<double,3> &X, DistSVec<double,dim> &V)
				

{

  *VCap = 0.0;      // contains test filtered values of the primitive variables //
  *Mom_Test = 0.0;  // contains test filtered values for computing cs from momentum equation //
  *Eng_Test = 0.0;  // contains test filtered values for computing pt from energy equation //

  domain->computeTestFilterValues(*VCap, *Mom_Test, *Eng_Test, Cs, VolSum, X, V, gam, R);

}

//------------------------------------------------------------------------
