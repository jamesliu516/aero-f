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
#include "LevelSet/LevelSetStructure.h"


//------------------------------------------------------------------------------
//--------------functions in ElemSet class
//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeTimeStep(FemEquationTerm *fet, SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idtv)
{
  int *nodeNumber,numberOfNodes;
  double Vmid[dim],oneOnDofs;
  double Xmid[3];
  double localDt;
  double nGrad[4][3];
  int    nodeID;
  for (int i=0; i<numElems; ++i)
    {
      nodeNumber    = elems[i]->nodeNum();
      numberOfNodes = elems[i]->numNodes();

      oneOnDofs = 1.0/((double) numberOfNodes);

      for (int k=0; k<dim; ++k) 
	{
	  Vmid[k] = 0.0;
	  for(int j=0;j<numberOfNodes;++j)
	    { 
	      nodeID   = nodeNumber[j];
	      Vmid[k] += V[nodeID][k];	  
	    }
	  Vmid[k] *= oneOnDofs;
	}
      // localDt_i = mu/rho*||n_i||^2/|T|
      // This line gives: 9.0*mu/rho*|T|
      localDt = 9.0*fet->computeViscousTimeStep(Xmid, Vmid)*elems[i]->computeVolume(X);
      // gradPhi = n_i/(3.0*|T|)
      elems[i]->computeGradientP1Function(X,nGrad);
      for(int j=0;j<numberOfNodes;++j)
	{
	  nodeID = nodeNumber[j];
	  // And we got what we need
	  idtv[nodeID]  += localDt*(nGrad[j][0]*nGrad[j][0] + nGrad[j][1]*nGrad[j][1] + nGrad[j][2]*nGrad[j][2]);
	}
    }
}
//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &R,
				  Vec<GhostPoint<dim>*> *ghostPoints,LevelSetStructure *LSS)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  for (int i=0; i<numElems; ++i) {
    elems[i]->computeGalerkinTerm(fet, X, d2wall, V, R, ghostPoints,LSS);
  }
}

//------------------------------------------------------------------------------

// Included
template<int dim>
void ElemSet::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, GeoState &geoState,
				 SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV,
				 double dMach, SVec<double,dim> &dR)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDerivativeOfGalerkinTerm(fet, X, dX, d2wall, V, dV, dMach, dR);

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
void ElemSet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, 
                                    SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R)

{

 for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicLESTerm(dles, Cs, X, V, R);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void ElemSet::computeJacobianGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol,
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                          Vec<GhostPoint<dim>*>* ghostPoints,LevelSetStructure *LSS)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeJacobianGalerkinTerm(fet, X, ctrlVol, d2wall, V, A, ghostPoints,LSS);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeTestFilterAvgs(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                   SVec<double,6> &Sij_Test, Vec<double> &modS_Test, 
                                   SVec<double,8> &Eng_Test, SVec<double,3> &X, SVec<double,dim> &V, 
                                   double gam, double R)
{

 for (int i=0; i<numElems; ++i)
   elems[i]->computeP1Avg(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R);

}

//------------------------------------------------------------------------------
// Level Set Reinitialization functions

template<int dimLS>
void ElemSet::computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                       SVec<double,dimLS> &ddz,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//------------------------------------------------------------------------------
template<int dimLS>
void ElemSet::recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                       SVec<double,dimLS> &ddz,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//------------------------------------------------------------------------------
template<int dimLS>
void ElemSet::computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                       SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);

}
// End of Level Set Reinitialization functions
//-------------------------------------------------------------------------------

