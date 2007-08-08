#include <Elem.h>
#include <Face.h>

#include <FemEquationTerm.h>
#include <MacroCell.h>
#include <VMSLESTerm.h>
#include <DynamicVMSTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <GenMatrix.h>
#include <math.h>
#include <GeoState.h>


//------------------------------------------------------------------------------
//--------------functions in ElemSet class
//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &R)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeGalerkinTerm(fet, X, d2wall, V, R);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeMBarAndM(DynamicVMSTerm *dvmst,
			      SVec<double,dim> **VBar,
			      SVec<double,1> **volRatio,
			      SVec<double,3> &X,
			      SVec<double,dim> &V,
			      SVec<double,dim> &MBar,
			      SVec<double,dim> &M)
{

  for (int i=0; i<numElems; ++i)
   elems[i]->computeMBarAndM(dvmst, VBar, volRatio, X, V, MBar, M);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeDynamicVMSTerm(DynamicVMSTerm *dvmst,
				    SVec<double,dim> **VBar,
				    SVec<double,3> &X,
				    SVec<double,dim> &V, SVec<double,dim> &S,
				    Vec<double> &CsDelSq, Vec<double> &PrT,
				    Vec<double> *Cs, Vec<double> &Delta)
{

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicVMSTerm(dvmst, VBar, X, V, S, CsDelSq, PrT, Cs, Delta);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeVMSLESTerm(VMSLESTerm *vmst,
				SVec<double,dim> &VBar,
				SVec<double,3> &X,
				SVec<double,dim> &V,
				SVec<double,dim> &Sigma)
                                                                                                                          
{
                                                                                                                          
  for (int i=0; i<numElems; ++i)
    elems[i]->computeVMSLESTerm(vmst, VBar, X, V, Sigma);
                                                                                                                          
}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
					SVec<double,dim> &V, SVec<double,dim> &R)

{
  for (int i=0; i<numElems; ++i)
    elems[i]->computeSmagorinskyLESTerm(smag, X, V, R);
}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
				SVec<double,dim> &V, SVec<double,dim> &R)

{
  for (int i=0; i<numElems; ++i)
    elems[i]->computeWaleLESTerm(wale, X, V, R);
}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, Vec<double> &VolSum,
				    SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R)

{

 for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicLESTerm(dles, Cs, VolSum, X, V, R);

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void ElemSet::computeJacobianGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol,
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeJacobianGalerkinTerm(fet, X, ctrlVol, d2wall, V, A);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeTestFilterAvgs(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
				    SVec<double,6> &Eng_Test, SVec<double,3> &X, SVec<double,dim> &V, double gam,
				    double R)

{

 for (int i=0; i<numElems; ++i)
   elems[i]->computeP1Avg(VCap, Mom_Test, Eng_Test, X, V, gam, R);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeCsValues(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
			      SVec<double,6> &Eng_Test, SVec<double,2> &Cs, 
			      Vec<double> &VolSum, SVec<double,3> &X, double gam, double R)

{

 for (int i=0; i<numElems; ++i)
   elems[i]->computeCsValues(VCap, Mom_Test, Eng_Test, Cs, VolSum, X, gam, R);

}

//------------------------------------------------------------------------------

