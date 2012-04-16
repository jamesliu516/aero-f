#include <cstdio>
#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif

#include <FluxFcn.h>
#include <RecFcn.h>
#include <MacroCell.h>
#include <VMSLESTerm.h>
#include <DynamicVMSTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <FemEquationTerm.h>
#include <VolumicForceTerm.h>
#include <PostFcn.h>
#include <BcFcn.h>
#include <BcDef.h>
#include <NodalGrad.h>
#include <EdgeGrad.h>
#include <Extrapolation.h>
#include <ExactRiemannSolver.h>
#include <BcData.h>
#include <GeoState.h>
#include <Vector.h>
#include <MvpMatrix.h>
#include <SparseMatrix.h>
#include <DenseMatrixOps.h>
#include <Connectivity.h>
#include <MemoryPool.h>
#include <Communicator.h>
#include <BinFileHandler.h>
#include <VectorSet.h>
#include <LinkF77.h>
#include <LowMachPrec.h>
#include "FluidSelector.h"
#include <GhostPoint.h>
#include <DenseMatrixOps.h>
#include <limits>
#include <PolygonReconstructionData.h> 

extern "C" {
  void F77NAME(mvp5d)(const int &, const int &, int *, int *, int (*)[2],
		      double (*)[25], double (*)[5], double (*)[5]);
  void F77NAME(torsionspring)(double (*)[3], int [4], double (*)[12], double &invCoef);
  void F77NAME(ballvertex)(double (*)[3], double(*)[3], int [4], double (*)[12], double &invCoef);
};

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &dt,
                                Vec<double> &idti, Vec<double> &idtv,
                                TimeLowMachPrec &tprec)
{
  dt = 0.0;
  idti = 0.0;
  idtv = 0.0;
  edges.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec);
  // Adam 2010.06.01: compute time step by tetrahedra, as the viscous term is computed
  if(fet) elems.computeTimeStep(fet,X,V,idtv);
  faces.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV,
                                Vec<double> &dIdti, Vec<double> &dIdtv, double dMach,
                                TimeLowMachPrec &tprec)
{

  dIdti = 0.0;
  dIdtv = 0.0;

  edges.computeDerivativeOfTimeStep(fet, varFcn, geoState, X,  dX, V, dV, dIdti, dIdtv, dMach, tprec);
  faces.computeDerivativeOfTimeStep(fet, varFcn, geoState, X,  dX, V, dV, dIdti, dIdtv, dMach, tprec);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,dim> &V, Vec<double> &dt,
				Vec<double> &idti, Vec<double> &idtv,
                                TimeLowMachPrec &tprec,
				Vec<int> &fluidId, Vec<double>* umax)
{

  dt = 0.0;

  edges.computeTimeStep(varFcn, geoState, V, dt, tprec, fluidId, globSubNum,umax);
  faces.computeTimeStep(varFcn, geoState, V, dt, tprec, fluidId);

}

//------------------------------------------------------------------------------

inline
void computeLocalWeightsLeastSquares(double dx[3], double *R, double *W)
{

  if(R[0]*R[3]*R[5] == 0.0) fprintf(stderr, "Going to divide by 0 %e %e %e\n",
         R[0], R[3], R[5]);
  double or11 = 1.0 / R[0];
  double or22 = 1.0 / R[3];
  double or33 = 1.0 / R[5];

  double r12or11 = R[1] * or11;
  double r23or22 = R[4] * or22;

  double psi = (R[1]*R[4] - R[2]*R[3]) * or11* or22;

  double alpha1 = dx[0] * or11 * or11;
  double alpha2 = (dx[1] - r12or11*dx[0]) * or22 * or22;
  double alpha3 = (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * or33;

  W[0] = alpha1 - r12or11*alpha2 + psi*alpha3;
  W[1] = alpha2 - r23or22*alpha3;
  W[2] = alpha3;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void computeDerivativeOfLocalWeightsLeastSquares(double dx[3], double ddx[3], double *R, double *dR, double *W, double *dW)
{

  if(R[0]*R[3]*R[5] == 0.0) fprintf(stderr, "Going to divide by 0 %e %e %e\n",
         R[0], R[3], R[5]);
  double or11 = 1.0 / R[0];
  double dor11 = -1.0 / ( R[0]*R[0] )*dR[0];
  double or22 = 1.0 / R[3];
  double dor22 = -1.0 / ( R[3]*R[3] )*dR[3];
  double or33 = 1.0 / R[5];
  double dor33 = -1.0 / ( R[5]*R[5] )*dR[5];

  double r12or11 = R[1] * or11;
  double dr12or11 = dR[1] * or11 + R[1] * dor11;
  double r23or22 = R[4] * or22;
  double dr23or22 = dR[4] * or22 + R[4] * dor22;

  double psi = (R[1]*R[4] - R[2]*R[3]) * or11 * or22;
  double dpsi = (dR[1]*R[4] + R[1]*dR[4] - dR[2]*R[3] - R[2]*dR[3]) * or11 * or22 + (R[1]*R[4] - R[2]*R[3]) * dor11 * or22 + (R[1]*R[4] - R[2]*R[3]) * or11 * dor22;

  double alpha1 = dx[0] * or11 * or11;
  double dalpha1 = ddx[0] * or11 * or11 + dx[0] * dor11 * or11 + dx[0] * or11 * dor11;
  double alpha2 = (dx[1] - r12or11*dx[0]) * or22 * or22;
  double dalpha2 = (ddx[1] - dr12or11*dx[0] - r12or11*ddx[0]) * or22 * or22 + (dx[1] - r12or11*dx[0]) * dor22 * or22 + (dx[1] - r12or11*dx[0]) * or22 * dor22;
  double alpha3 = (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * or33;
  double dalpha3 = (ddx[2] - dr23or22*dx[1] - r23or22*ddx[1] + dpsi*dx[0] + psi*ddx[0]) * or33 * or33 + (dx[2] - r23or22*dx[1] + psi*dx[0]) * dor33 * or33 + (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * dor33;

  W[0] = alpha1 - r12or11*alpha2 + psi*alpha3;
  W[1] = alpha2 - r23or22*alpha3;
  W[2] = alpha3;

  dW[0] = dalpha1 - dr12or11*alpha2 - r12or11*dalpha2 + dpsi*alpha3  + psi*dalpha3;
  dW[1] = dalpha2 - dr23or22*alpha3 - r23or22*dalpha3;
  dW[2] = dalpha3;
}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsLeastSquares(SVec<double,3> &X, SVec<double,6> &R,
	        SVec<Scalar,dim> &var, SVec<Scalar,dim> &ddx,
	        SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)  {	//KTC: CANNOT MODIFY

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	if (sampleMesh) {

		for (int iEdge=0; iEdge<edges.getNumTwoLayersEdges(); ++iEdge) {

			int l = edges.edgesTwoLayersSampleNode[iEdge];

			if (!edgeFlag[l])
				continue;

			int i = edgePtr[l][0];
			int j = edgePtr[l][1];

			double Wi[3], Wj[3];
			Scalar deltaVar;

			double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
			computeLocalWeightsLeastSquares(dx, R[i], Wi);

			dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
			computeLocalWeightsLeastSquares(dx, R[j], Wj);

			for (int k=0; k<dim; ++k) {
				deltaVar = var[j][k] - var[i][k];

				ddx[i][k] += Wi[0] * deltaVar;
				ddy[i][k] += Wi[1] * deltaVar;
				ddz[i][k] += Wi[2] * deltaVar;
				ddx[j][k] -= Wj[0] * deltaVar;
				ddy[j][k] -= Wj[1] * deltaVar;
				ddz[j][k] -= Wj[2] * deltaVar;
			}
		}

	}

	else {
		for (int l=0; l<edges.size(); ++l) {

			if (!edgeFlag[l])
				continue;

			int i = edgePtr[l][0];
			int j = edgePtr[l][1];

			double Wi[3], Wj[3];
			Scalar deltaVar;

			double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
			computeLocalWeightsLeastSquares(dx, R[i], Wi);

			dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
			computeLocalWeightsLeastSquares(dx, R[j], Wj);

			for (int k=0; k<dim; ++k) {
				deltaVar = var[j][k] - var[i][k];

				ddx[i][k] += Wi[0] * deltaVar;
				ddy[i][k] += Wi[1] * deltaVar;
				ddz[i][k] += Wi[2] * deltaVar;
				ddx[j][k] -= Wj[0] * deltaVar;
				ddy[j][k] -= Wj[1] * deltaVar;
				ddz[j][k] -= Wj[2] * deltaVar;
			}
		}
	}

}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow and FSI)
template<int dim, class Scalar>
void SubDomain::computeGradientsLeastSquares(SVec<double,3> &X,
                const Vec<int> &fluidId, SVec<double,6> &R,
                SVec<Scalar,dim> &var, SVec<Scalar,dim> &ddx,
                SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz,
                bool linRecFSI, LevelSetStructure *LSS)  {

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    if(fluidId[i]!=fluidId[j] || (LSS && LSS->edgeIntersectsStructure(0.0,l))) continue;

    if (higherOrderMF && (higherOrderMF->isCellCut(i) || 
                          higherOrderMF->isCellCut(j)))
      continue;

    double Wi[3], Wj[3];
    Scalar deltaVar;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    if(R[i][0]>0.0 && fabs(R[i][0]*R[i][3]*R[i][5]) > 1.0e-10) // should be positive for a well posed least square problem
      computeLocalWeightsLeastSquares(dx, R[i], Wi);
    else{ // gradient is set to 0.0
      Wi[0] = 0.0;
      Wi[1] = 0.0;
      Wi[2] = 0.0;
    }

    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    if(R[j][0]>0.0 && fabs(R[j][0]*R[j][3]*R[j][5]) > 1.0e-10) // should be positive for a well posed least square problem
      computeLocalWeightsLeastSquares(dx, R[j], Wj);
    else{ // gradient is set to 0.0
      Wj[0] = 0.0;
      Wj[1] = 0.0;
      Wj[2] = 0.0;
    }

    for (int k=0; k<dim; ++k) {
      deltaVar = var[j][k] - var[i][k];

      /*if (fabs(deltaVar) > 1000)
        std::cout << "Error: deltaVar is huge: " << deltaVar << std::endl;
 
      if (fabs(Wi[1]+Wi[0]+Wi[2]) > 1000) {
        std::cout << "Error: Wi[0] is huge: " << Wi[0] << " " << Wi[1] << " " << Wi[2] << std::endl;
        std::cout << "R[5] = " << R[i][5] << std::endl;
      }*/

      ddx[i][k] += Wi[0] * deltaVar;
      ddy[i][k] += Wi[1] * deltaVar;
      ddz[i][k] += Wi[2] * deltaVar;
      ddx[j][k] -= Wj[0] * deltaVar;
      ddy[j][k] -= Wj[1] * deltaVar;
      ddz[j][k] -= Wj[2] * deltaVar;
    }
  }

//KW: set gradients = 0 for cells near interface.
  if(!linRecFSI)
    for (int l=0; l<edges.size(); ++l) {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];
      if (fluidId[i]!=fluidId[j] || (LSS && LSS->edgeIntersectsStructure(0.0,l))) {
        for (int k=0; k<dim; ++k)
          ddx[i][k] = ddy[i][k] = ddz[i][k] = ddx[j][k] = ddy[j][k] = ddz[j][k] = 0.0;
      }
    }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOfGradientsLeastSquares(SVec<double,3> &X, SVec<double,3> &dX,
                         SVec<double,6> &R, SVec<double,6> &dR,
					     SVec<Scalar,dim> &var, SVec<Scalar,dim> &dvar, SVec<Scalar,dim> &dddx,
                         SVec<Scalar,dim> &dddy, SVec<Scalar,dim> &dddz)
{

  dddx = (Scalar) 0.0;
  dddy = (Scalar) 0.0;
  dddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double Wi[3], Wj[3], deltaVar;
    double dWi[3], dWj[3], dDeltaVar;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    double ddx[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};
    computeDerivativeOfLocalWeightsLeastSquares(dx, ddx, R[i], dR[i], Wi, dWi);

    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    ddx[0] = -ddx[0]; ddx[1] = -ddx[1]; ddx[2] = -ddx[2];
    computeDerivativeOfLocalWeightsLeastSquares(dx, ddx, R[j], dR[j], Wj, dWj);

    for (int k=0; k<dim; ++k) {
      deltaVar = var[j][k] - var[i][k];
      dDeltaVar = dvar[j][k] - dvar[i][k];

      dddx[i][k] += (dWi[0] * deltaVar + Wi[0] * dDeltaVar);
      dddy[i][k] += (dWi[1] * deltaVar + Wi[1] * dDeltaVar);
      dddz[i][k] += (dWi[2] * deltaVar + Wi[2] * dDeltaVar);
      dddx[j][k] -= (dWj[0] * deltaVar + Wj[0] * dDeltaVar);
      dddy[j][k] -= (dWj[1] * deltaVar + Wj[1] * dDeltaVar);
      dddz[j][k] -= (dWj[2] * deltaVar + Wj[2] * dDeltaVar);
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsGalerkin(Vec<double> &ctrlVol, SVec<double,3> &wii,
		SVec<double,3> &wij, SVec<double,3> &wji, SVec<Scalar,dim> &var,
                SVec<Scalar,dim> &ddx, SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)
{

  int i, j, k, l;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0];
      ddy[i][k] = var[i][k] * wii[i][1];
      ddz[i][k] = var[i][k] * wii[i][2];
    }
  }

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {

    //if (!edgeFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wij[l][0] * var[j][k];
      ddx[j][k] += wji[l][0] * var[i][k];

      ddy[i][k] += wij[l][1] * var[j][k];
      ddy[j][k] += wji[l][1] * var[i][k];

      ddz[i][k] += wij[l][2] * var[j][k];
      ddz[j][k] += wji[l][2] * var[i][k];

    }

  }

  for (i=0; i<var.size(); ++i)  {

    double coef = 1.0 / (4.0*ctrlVol[i]);

    for (k=0; k<dim; ++k) {
      ddx[i][k] *= coef;
      ddy[i][k] *= coef;
      ddz[i][k] *= coef;
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOfGradientsGalerkin(Vec<double> &ctrlVol, Vec<double> &dCtrlVol,
					 SVec<double,3> &wii, SVec<double,3> &wij, SVec<double,3> &wji,
					 SVec<double,3> &dwii, SVec<double,3> &dwij, SVec<double,3> &dwji,
					 SVec<Scalar,dim> &var, SVec<Scalar,dim> &dvar, SVec<Scalar,dim> &ddx,
                     SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz, SVec<Scalar,dim> &dddx,
                     SVec<Scalar,dim> &dddy, SVec<Scalar,dim> &dddz)
{

  int i, j, k, l;

  ddx=0.0;
  ddy=0.0;
  ddz=0.0;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0];
      ddy[i][k] = var[i][k] * wii[i][1];
      ddz[i][k] = var[i][k] * wii[i][2];

      dddx[i][k] = var[i][k] * dwii[i][0] + dvar[i][k] * wii[i][0];
      dddy[i][k] = var[i][k] * dwii[i][1] + dvar[i][k] * wii[i][1];
      dddz[i][k] = var[i][k] * dwii[i][2] + dvar[i][k] * wii[i][2];
    }
  }

  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {

    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wij[l][0] * var[j][k];
      ddx[j][k] += wji[l][0] * var[i][k];

      ddy[i][k] += wij[l][1] * var[j][k];
      ddy[j][k] += wji[l][1] * var[i][k];

      ddz[i][k] += wij[l][2] * var[j][k];
      ddz[j][k] += wji[l][2] * var[i][k];

      dddx[i][k] += dwij[l][0] * var[j][k] + wij[l][0] * dvar[j][k];
      dddx[j][k] += dwji[l][0] * var[i][k] + wji[l][0] * dvar[i][k];

      dddy[i][k] += dwij[l][1] * var[j][k] + wij[l][1] * dvar[j][k];
      dddy[j][k] += dwji[l][1] * var[i][k] + wji[l][1] * dvar[i][k];

      dddz[i][k] += dwij[l][2] * var[j][k] + wij[l][2] * dvar[j][k];
      dddz[j][k] += dwji[l][2] * var[i][k] + wji[l][2] * dvar[i][k];

    }

  }

  for (i=0; i<var.size(); ++i)  {

    double coef = 1.0 / (4.0*ctrlVol[i]);

    double dcoef = -1.0 / (4.0*ctrlVol[i]*ctrlVol[i]) * dCtrlVol[i];

    for (k=0; k<dim; ++k) {
      dddx[i][k] *= coef;
      dddx[i][k] += ddx[i][k]*dcoef;
      dddy[i][k] *= coef;
      dddy[i][k] += ddy[i][k]*dcoef;
      dddz[i][k] *= coef;
      dddz[i][k] += ddz[i][k]*dcoef;
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsGalerkinT(Vec<double> &ctrlVol,
                SVec<double,3> &wii, SVec<double,3> &wij, SVec<double,3> &wji,
                SVec<Scalar,dim> &var, SVec<Scalar,dim> &var1,
                SVec<Scalar,dim> &var2,SVec<Scalar,dim> &ddx,
                SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)  {

  int i, j, k, l;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0] / (4.0*ctrlVol[i]);
      ddy[i][k] = var1[i][k] * wii[i][1] / (4.0*ctrlVol[i]);
      ddz[i][k] = var2[i][k] * wii[i][2] / (4.0*ctrlVol[i]);
    }
  }

  //bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wji[l][0] * var[j][k] / (4.0*ctrlVol[j]);
      ddx[j][k] += wij[l][0] * var[i][k] / (4.0*ctrlVol[i]);

      ddy[i][k] += wji[l][1] * var1[j][k] / (4.0*ctrlVol[j]);
      ddy[j][k] += wij[l][1] * var1[i][k] / (4.0*ctrlVol[i]);

      ddz[i][k] += wji[l][2] * var2[j][k] / (4.0*ctrlVol[j]);
      ddz[j][k] += wij[l][2] * var2[i][k] / (4.0*ctrlVol[i]);


    }
  }

}


//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMinMaxStencilValues(SVec<double,dim> &V, SVec<double,dim> &Vmin,
					   SVec<double,dim> &Vmax)
{

  Vmin = V;
  Vmax = V;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    for (int k=0; k<dim; ++k) {
      Vmin[i][k] = min(Vmin[i][k], V[j][k]);
      Vmax[i][k] = max(Vmax[i][k], V[j][k]);
      Vmin[j][k] = min(Vmin[j][k], V[i][k]);
      Vmax[j][k] = max(Vmax[j][k], V[i][k]);
    }

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfMinMaxStencilValues(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &Vmin, SVec<double,dim> &dVmin,
					   SVec<double,dim> &Vmax, SVec<double,dim> &dVmax)
{

  Vmin = V;
  Vmax = V;
  dVmin = dV;
  dVmax = dV;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    for (int k=0; k<dim; ++k) {

      if (min(Vmin[i][k], V[j][k]) == V[j][k]) {
        Vmin[i][k] = V[j][k];
        dVmin[i][k] = dV[j][k];
      }

      if (max(Vmax[i][k], V[j][k]) == V[j][k]) {
        Vmax[i][k] = V[j][k];
        dVmax[i][k] = dV[j][k];
      }

      if (min(Vmin[j][k], V[i][k]) == V[i][k]) {
        Vmin[j][k] = V[i][k];
        dVmin[j][k] = dV[i][k];
      }

      if (max(Vmax[j][k], V[i][k]) == V[i][k]) {
        Vmax[j][k] = V[i][k];
        dVmax[j][k] = dV[i][k];
      }

    }

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMultiDimLimiter(RecLimiter *recFcn, SVec<double,3> &X,
				       Vec<double> &ctrlVol, SVec<double,dim> &V,
				       SVec<double,dim> &dVdx, SVec<double,dim> &dVdy,
				       SVec<double,dim> &dVdz, SVec<double,dim> &Vmin,
				       SVec<double,dim> &Vmax, SVec<double,dim> &phi)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim];

  phi = 1.0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    recFcn->computeLimiter(Vmax[i], Vmin[i], V[i], Vi, ctrlVol[i],
			   Vmax[j], Vmin[j], V[j], Vj, ctrlVol[j], phi[i], phi[j]);

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfMultiDimLimiter(RecLimiter *recFcn, SVec<double,3> &X, SVec<double,3> &dX,
				       Vec<double> &ctrlVol, Vec<double> &dCtrlVol, SVec<double,dim> &V, SVec<double,dim> &dV,
				       SVec<double,dim> &dVdx, SVec<double,dim> &dVdy, SVec<double,dim> &dVdz,
				       SVec<double,dim> &ddVdx, SVec<double,dim> &ddVdy, SVec<double,dim> &ddVdz,
                       SVec<double,dim> &Vmin, SVec<double,dim> &dVmin, SVec<double,dim> &Vmax,
                       SVec<double,dim> &dVmax, SVec<double,dim> &phi, SVec<double,dim> &dphi)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim];

  double dddVij[dim], dddVji[dim], dVi[dim], dVj[dim];

  phi = 1.0;
  dphi = 0.0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    double ddx[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      dddVij[k] = ddx[0]*dVdx[i][k] + dx[0]*ddVdx[i][k] + ddx[1]*dVdy[i][k] + dx[1]*ddVdy[i][k] + ddx[2]*dVdz[i][k] + dx[2]*ddVdz[i][k];
      dddVji[k] = ddx[0]*dVdx[j][k] + dx[0]*ddVdx[j][k] + ddx[1]*dVdy[j][k] + dx[1]*ddVdy[j][k] + ddx[2]*dVdz[j][k] + dx[2]*ddVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    recFcn->computeDerivative(V[i], dV[i], ddVij, dddVij, V[j], dV[j], ddVji, dddVji, dVi, dVj);

    recFcn->computeDerivativeOfLimiter(Vmax[i], dVmax[i], Vmin[i], dVmin[i], V[i], dV[i], Vi, dVi, ctrlVol[i], dCtrlVol[i],
			   Vmax[j], dVmax[j], Vmin[j], dVmin[j], V[j], dV[j], Vj, dVj, ctrlVol[j], dCtrlVol[j], phi[i], dphi[i], phi[j], dphi[j]);

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computePressureSensor(SVec<double,3>& X, SVec<double,dim>& V,
				      SVec<double,dim>& dVdx, SVec<double,dim>& dVdy,
				      SVec<double,dim>& dVdz, SVec<double,3>& sensor)
{

  sensor = 0.0;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dpi = dx[0]*dVdx[i][4] + dx[1]*dVdy[i][4] + dx[2]*dVdz[i][4];
    double dpj = dx[0]*dVdx[j][4] + dx[1]*dVdy[j][4] + dx[2]*dVdz[j][4];

    double s0 = fabs(dpj - dpi);
    double s1 = fabs(dpj) + fabs(dpi);
    double s2 = fabs(V[i][4]) + fabs(V[j][4]);

    sensor[i][0] += s0;
    sensor[i][1] += s1;
    sensor[i][2] += s2;
    sensor[j][0] += s0;
    sensor[j][1] += s1;
    sensor[j][2] += s2;
  }

}

//------------------------------------------------------------------------------

/*
@INPROCEEDINGS{dervieux-85,
  author = "Dervieux, A.",
  title = "Steady {E}uler simulations using unstructured meshes",
  booktitle = "Proceedings of the VKI Lectures Series 1985-04",
  series = "16th Computational Fluid Dynamics",
  month = mar,
  year = 1985,
  address = "von Karman Institute, Brussels, Belgium",
}
@INBOOK{barth-92,
  author = "Barth, T. J.",
  title = "Aspects of Unstructured Grids and Finite-Volume Solvers
           for the {E}uler and {N}avier-{S}tokes Equations",
  series = "Special Course on Unstructured Grid Methods for Advection
            Dominated Flows",
  publisher = "AGARD R-787",
  month = may,
  year = 1992,
  chapter = "6",
}
@ARTICLE{vanleer-79,
  author = "van Leer, B.",
  title = "Towards the Ultimate Conservative Difference Scheme. {V}.
           {A} Second-Order Sequel to {G}odunov's Method",
  journal = jcp,
  year = 1979,
  volume = 32,
  pages = "101--136",
}
*/
template<int dim>
int SubDomain::computeFiniteVolumeTerm(Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       SVec<double,dim>& fluxes, SVec<int,2>& tag,
                                       int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, V, ngrad, egrad, fluxes, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluxes);

  return(ierr);

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, 
                                       Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
				       SVec<double,3>& X, SVec<double,dim>& V,
				       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
				       SVec<double,dim>& fluxes, SVec<int,2>& tag,
                                       int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, V, ngrad, egrad, fluxes, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(riemann, fluxFcn, bcData, geoState, V, fluxes);

  return(ierr);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfFiniteVolumeTerm(Vec<double> &irey, Vec<double> &dIrey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,3>& X, SVec<double,3>& dX, SVec<double,dim>& V, SVec<double,dim>& dV,
                                        NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad, double dMach,
                                        SVec<double,dim>& dFluxes)
{

  edges.computeDerivativeOfFiniteVolumeTerm(irey, dIrey, fluxFcn, recFcn, elems, geoState, X, dX, V, dV, ngrad, egrad, dMach, dFluxes);

  faces.computeDerivativeOfFiniteVolumeTerm(fluxFcn, bcData, geoState, V, dFluxes);

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       Vec<int> &fluidId, FluidSelector &fluidSelector,
                                       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
				       SVec<double,dimLS>& phi,
                                       NodalGrad<dimLS>& ngradLS,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift)
{
  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, fluidId, fluidSelector,
                                           ngrad, egrad, phi,ngradLS, fluxes, it,
                                           tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji, 
                                       LevelSetStructure& LSS, bool linRecAtInterface, 
                                       Vec<int> &fluidId, int Nriemann, SVec<double,3>* Nsbar, 
                                       FluidSelector &fluidSelector,
                                       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       NodalGrad<dimLS>& ngradLS,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift)
{
  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, Wstarij, Wstarji, LSS, linRecAtInterface,
                                           fluidId, Nriemann, Nsbar, fluidSelector,
                                           ngrad, egrad, ngradLS, fluxes, it,
                                           tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes, &LSS);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji,
                                       LevelSetStructure &LSS, bool linRecAtInterface, Vec<int> &fluidId,
                                       int Nriemann, SVec<double,3>* Nsbar, NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, Wstarij, Wstarji, LSS, 
                                           linRecAtInterface, fluidId, Nriemann, Nsbar, ngrad, egrad, fluxes, it,
                                           tag, failsafe, rshift);
  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes, &LSS);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void SubDomain::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,3>& X, SVec<double,dim>& V,
                                        NodalGrad<dim>& ngrad, NodalGrad<dimLS> &ngradLS,
                                        EdgeGrad<dim>* egrad,
                                        SVec<double,dimLS>& Phi, SVec<double,dimLS> &PhiF,
                                        LevelSetStructure* LSS)
{
  edges.computeFiniteVolumeTermLS(fluxFcn, recFcn, recFcnLS, elems, geoState, X, V, ngrad, ngradLS,
                                  egrad, Phi, PhiF, LSS);

  faces.computeFiniteVolumeTermLS(fluxFcn, bcData, geoState, V, Phi, PhiF);
    //Note: LSS is not needed because we assume that farfield nodes cannot be covered by structure.
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeBar_Step1(Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                            BcData<dim>& bcData, GeoState& geoState,
                                            SVec<double,3>& X, SVec<double,dim>& VBar,
                                            NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                            SVec<double,dim>& sigma, SVec<int,2> &tag,
                                            int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, VBar, ngrad, egrad, sigma, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, VBar, sigma);

  return(ierr);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFiniteVolumeBar_Step2(MacroCellSet **macroCells,
                                             SVec<double,1> &volRatio,
                                             SVec<double,dim>& sigma,
                                             SVec<double,dim>& fluxes,
                                             int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         int idx = nToMN[i][j];
         for (int k=0; k<dim; ++k)
               fluxes[idx][k] += volRatio[idx][0] * sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState, Vec<double> &irey,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                CommPattern<double>* flag)
{
  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A);

    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A, nodeType);

    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                                FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState, Vec<double> &irey,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                CommPattern<double>* flag)
{
  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A);

    faces.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, bcData, geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A, nodeType);

    fprintf(stderr,"WARNING: Exact Riemann solver at the wall has not been implemented in this case!\n");
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeRHS(VarFcn* vf, SVec<double,dim>& V, SVec<double,dim>& rhs,
                             Extrapolation<dim>* xpol, BcData<dim>& bcData,
                             GeoState& geoState, SVec<double,3> &X)
{
  inletNodes.recomputeRHS(vf, xpol, elems, V, bcData, geoState, rhs, X, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeRHS(VarFcn* vf, SVec<double,dim>& V, Vec<int> &fluidId,
                            SVec<double,dim>& rhs, Extrapolation<dim>* xpol,
                            BcData<dim>& bcData, GeoState& geoState, SVec<double,3> &X)
{
  inletNodes.recomputeRHS(vf, xpol, elems, V, fluidId, bcData, geoState, rhs, X, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeResidual(SVec<double,dim> &F, SVec<double,dim> &Finlet)
{
  Finlet = 0.0;
  inletNodes.recomputeResidual(F,Finlet);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeRealFluidResidual(SVec<double, dim> &F, SVec<double,dim> &Freal, LevelSetStructure &lss)
{
  Freal = 0.0;
  for (int iNode=0; iNode<numNodes(); iNode++)
    if (lss.isActive(0,iNode))
        for (int j=0; j<dim; j++)  Freal[iNode][j] = F[iNode][j];
}




template<class Scalar, int dim>
void SubDomain::checkRHS(Scalar (*rhs)[dim])
{
  int node;
  for (int i = 0; i<inletNodes.size(); i++){
    node = inletNodes[i].getNodeNum();
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                                FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState,
                                                NodalGrad<dim> &ngrad, NodalGrad<dimLS> &ngradLS,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                FluidSelector &fluidSelector,
                                                Vec<int> &fluidId, CommPattern<double>* flag)
{
  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, geoState, ngrad, ngradLS, X, ctrlVol, V, A, fluidSelector, fluidId);
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId);
  }else{
    edges.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, geoState, ngrad, ngradLS, X, ctrlVol, V, A, fluidSelector, fluidId, nodeType);
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}
//-------------------------------------------------------------------------------
template<class Scalar,int dim,int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       LevelSetStructure &LSS,Vec<int> &fluidId,
                                       int Nriemann, SVec<double,3>* Nsbar,
                                       GenMat<Scalar,neq>& A,Vec<double>& irey)
{

  edges.computeJacobianFiniteVolumeTerm(riemann,fluxFcn,
                                        geoState, X, V, ctrlVol, LSS, 
                                        fluidId, Nriemann, Nsbar, A,irey);

  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId, &LSS); 

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
  
}

template<int dim, class Scalar, int neq, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       NodalGrad<dimLS> &ngradLS,
                                       LevelSetStructure &LSS,Vec<int> &fluidId,
                                       int Nriemann, SVec<double,3>* Nsbar,
                                       FluidSelector &fluidSelector,
                                       GenMat<Scalar,neq>& A) {

  edges.computeJacobianFiniteVolumeTerm(riemann,locToGlobNodeMap,
                                        fluxFcn, geoState, X, V, LSS,fluidId,Nriemann, Nsbar,
                                        fluidSelector, ngradLS,ctrlVol, A);
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId,&LSS);
  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}

//-------------------------------------------------------------------------------

template<int dim, class Scalar, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
						  GeoState &geoState,SVec<double,3>& X,SVec<double,dim> &V,
						  NodalGrad<dim>& ngrad, NodalGrad<dimLS> &ngradLS,
						  EdgeGrad<dim>* egrad,Vec<double> &ctrlVol,SVec<double,dimLS>& Phi,
						  GenMat<Scalar,dimLS> &A, LevelSetStructure* LSS, CommPattern<double> * flag)
{

  if (!flag){
    edges.computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS,geoState,X,V,ngrad ,ngradLS,
					    egrad, ctrlVol , Phi, A,LSS);
    faces.computeJacobianFiniteVolumeTermLS(geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS,geoState,X,V,ngrad ,ngradLS,
					    egrad,ctrlVol , Phi, A,LSS);
    faces.computeJacobianFiniteVolumeTermLS(geoState, V, A );
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<dimLS*dimLS; ++k)
      Aii[k] *= voli;
  }
}
//-------------------------------------------------------------------------------

template<class Scalar, int neq>
void SubDomain::finishJacobianGalerkinTerm(Vec<double> &ctrlVol, GenMat<Scalar,neq> &A)  {

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVolumicForceTerm(VolumicForceTerm *volForce, Vec<double> &ctrlVol,
                               SVec<double,dim> &V, SVec<double,dim> &fluxes)
{

  double r[5];
  for (int i=0; i<ctrlVol.size(); i++){
    volForce->computeVolumeTerm(ctrlVol[i], V[i], r);
    for (int k=0; k<5; k++)
      fluxes[i][k] += r[k];
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfVolumicForceTerm(VolumicForceTerm *volForce, Vec<double> &ctrlVol, Vec<double> &dCtrlVol,
                               SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &dFluxes)
{

  double dr[5];
  for (int i=0; i<ctrlVol.size(); i++){
    volForce->computeDerivativeOfVolumeTerm(ctrlVol[i], dCtrlVol[i], V[i], dV[i], dr);
    for (int k=0; k<5; k++)
      dFluxes[i][k] += dr[k];
  }

}

//------------------------------------------------------------------------------
/*
@TECHREPORT{fezoui-lanteri-larrouturou-olivier-89,
  author = "Fezoui, F. and Lanteri, S. and Larrouturou, B. and Olivier, C.",
  title = "R\'esolution num\'erique des \'equations de {N}avier-{S}tokes pour un fluide
  compressible en maillage triangulaire",
  institution = "INRIA, France",
  number = 1033,
  year = 1989,
}
@UNPUBLISHED{barth-91a,
  author = "Barth, T. J.",
  title = "Numerical Aspects of Computing High {R}eynolds Number Flows
           on Unstructured Meshes",
  month = jan,
  year = 1991,
  note = "{AIAA} paper 91-0721",
}
*/
template<int dim>
void SubDomain::computeGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
				    GeoState &geoState, SVec<double,3> &X,
				    SVec<double,dim> &V, SVec<double,dim> &R,
				    Vec<GhostPoint<dim>*> *ghostPoints,LevelSetStructure *LSS)
{

  elems.computeGalerkinTerm(fet, geoState, X, V, R,ghostPoints,LSS);

  faces.computeGalerkinTerm(elems, fet, bcData, geoState, X, V, R,LSS);

}


//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
				    GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				    SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  elems.computeDerivativeOfGalerkinTerm(fet, geoState, X, dX, V, dV, dMach, dR);

  faces.computeDerivativeOfGalerkinTerm(elems, fet, bcData, geoState, X, dX, V, dV, dMach, dR);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
					  SVec<double,dim> &V, SVec<double,dim> &R)
{

  elems.computeSmagorinskyLESTerm(smag, X, V, R);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
			           SVec<double,dim> &V, SVec<double,dim> &R)
{

  elems.computeWaleLESTerm(wale, X, V, R);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTestFilterAvgs(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                      SVec<double,6> &Sij_Test, Vec<double> &modS_Test,
                                      SVec<double,8> &Eng_Test, SVec<double,3> &X,
				      SVec<double,dim> &V, double gam, double R)
{

  elems.computeTestFilterAvgs(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R);

}

//------------------------------------------------------------------------------
// Computes Cs and Pt values in dynamic LES procedure
//
template<int dim>
void SubDomain::computeCsValues(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                SVec<double,6> &Sij_Test, Vec<double> &modS_Test,
                                SVec<double,8> &Eng_Test, SVec<double,2> &Cs,
				Vec<int> &Ni, SVec<double,3> &X, double gam, double R)
{

 for (int i=0; i<nodes.size(); ++i) {
   double vc[5];
   double r_u[3];
   double r_u_u[6];
   double r_s_p[6];
   double ratdelta = pow(Ni[i],(2.0/3.0)); // should precompute this also

   double r_e;
   double r_e_plus_p;
   double r_s_dtdxj[3];
   double dtdxj[3];

   double num, denom;
   double sqrt2S2;
   double Pij[3][3], Bij[3][3], Lij[3][3];
   double Zi[3], Li[3];
   double oogam1 = 1.0/(gam - 1.0);
   double Cp = R*gam*oogam1; // specific heat at constant pressure

   // Compute Smagorinsky Coefficient at each node
   // ----------------------------------------------
   // Ref: Large Eddy Simulation of Bluff-Body flow on Unstructured Grids
   // International Journal of Numerical Methods in Fluids 2002
   // Vol : 40, pgs:1431-1460
   // Authors; Camarri, Salvetti, Koobus, Dervieux
   // ---------------------------------------------------------------------

    r_u[0] = Mom_Test[i][1]/Mom_Test[i][0];
    r_u[1] = Mom_Test[i][2]/Mom_Test[i][0];
    r_u[2] = Mom_Test[i][3]/Mom_Test[i][0];

    r_u_u[0] = Mom_Test[i][4]/Mom_Test[i][0];
    r_u_u[1] = Mom_Test[i][5]/Mom_Test[i][0];
    r_u_u[2] = Mom_Test[i][6]/Mom_Test[i][0];
    r_u_u[3] = Mom_Test[i][7]/Mom_Test[i][0];
    r_u_u[4] = Mom_Test[i][8]/Mom_Test[i][0];
    r_u_u[5] = Mom_Test[i][9]/Mom_Test[i][0];

    r_s_p[0] = Mom_Test[i][10]/Mom_Test[i][0];
    r_s_p[1] = Mom_Test[i][11]/Mom_Test[i][0];
    r_s_p[2] = Mom_Test[i][12]/Mom_Test[i][0];
    r_s_p[3] = Mom_Test[i][13]/Mom_Test[i][0];
    r_s_p[4] = Mom_Test[i][14]/Mom_Test[i][0];
    r_s_p[5] = Mom_Test[i][15]/Mom_Test[i][0];

    vc[0] = VCap[i][0]/Mom_Test[i][0];
    vc[1] = VCap[i][1]/Mom_Test[i][0];
    vc[2] = VCap[i][2]/Mom_Test[i][0];
    vc[3] = VCap[i][3]/Mom_Test[i][0];
    vc[4] = VCap[i][4]/Mom_Test[i][0];

    Pij[0][0] = Sij_Test[i][0]/Mom_Test[i][0];
    Pij[1][1] = Sij_Test[i][1]/Mom_Test[i][0];
    Pij[2][2] = Sij_Test[i][2]/Mom_Test[i][0];
    Pij[0][1] = Sij_Test[i][3]/Mom_Test[i][0];
    Pij[0][2] = Sij_Test[i][4]/Mom_Test[i][0];
    Pij[1][2] = Sij_Test[i][5]/Mom_Test[i][0];
    Pij[1][0] = Pij[0][1];
    Pij[2][0] = Pij[0][2];
    Pij[2][1] = Pij[1][2];

    sqrt2S2 = modS_Test[i]/Mom_Test[i][0];

    r_e = Eng_Test[i][0]/Mom_Test[i][0];

    r_e_plus_p = Eng_Test[i][1]/Mom_Test[i][0];

    r_s_dtdxj[0] = Eng_Test[i][2]/Mom_Test[i][0];
    r_s_dtdxj[1] = Eng_Test[i][3]/Mom_Test[i][0];
    r_s_dtdxj[2] = Eng_Test[i][4]/Mom_Test[i][0];

    dtdxj[0] = Eng_Test[i][5]/Mom_Test[i][0];
    dtdxj[1] = Eng_Test[i][6]/Mom_Test[i][0];
    dtdxj[2] = Eng_Test[i][7]/Mom_Test[i][0];

    // computing Smagorinsky constant //

    computeLij(Lij, r_u, r_u_u, vc);
    computeBij(Bij, r_s_p, sqrt2S2, Pij, ratdelta, vc);


    num = 0.0;
    denom = 0.0;

    // least squares procedure

    for (int j=0; j<3; ++j){
      for (int k=0; k<3; ++k){
         num += Bij[j][k]*Lij[j][k];
         denom += Bij[j][k]*Bij[j][k];
       }
    }

    if(denom < 1.0e-7) denom = 1.0e-7;
    Cs[i][0] = (num/denom);

    // computing turbulent Prandtl number //

    computeZi(Zi, ratdelta, sqrt2S2, dtdxj, r_s_dtdxj, vc, Cp);
    computeLi(Li, r_e, r_e_plus_p, r_u, vc);

    num = 0.0;
    denom = 0.0;

    // least squares procedure

    for (int j=0; j<3; ++j){
       num += Li[j]*Zi[j];
       denom += Zi[j]*Zi[j];
    }

    if(denom < 1.0e-7) denom = 1.0e-7;
    Cs[i][1] = (num/denom);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                      SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R)
{

  elems.computeDynamicLESTerm(dles, Cs, X, V, R);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVMSLES_Step1(VMSLESTerm *vmst,
                                    SVec<double,dim> &VBar,
                                    SVec<double,3> &X,
                                    SVec<double,dim> &V,
                                    SVec<double,dim> &Sigma)
{
  elems.computeVMSLESTerm(vmst, VBar,X, V, Sigma);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVMSLES_Step2(SVec<double,1> &volRatio,
                                    MacroCellSet *macroCells,
                                    SVec<double,dim> &Sigma,
                                    SVec<double,dim> &R,
                                    int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells->containing(i) != -1) {
      // excluding nodes on SD boundary
      // that are assigned to a MC in another SD
      for (int j=0; j<nToMN.num(i); ++j) {
	int idx = nToMN[i][j];
	if (i == idx) {
	  for (int k=0; k<dim; ++k)
	    R[idx][k] += (1.0 - volRatio[idx][0]) * Sigma[i][k];
	}
	else {
	  for (int k=0; k<dim; ++k)
	    R[idx][k] += -1.0 * volRatio[idx][0] * Sigma[i][k];
	}
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeGalerkinBar_Step1(FemEquationTerm *fet,
                                         BcData<dim> &bcData,
                                         GeoState &geoState,
                                         SVec<double,3> &X,
                                         SVec<double,dim> &VBar,
                                         SVec<double,dim> &Sigma)
{

  elems.computeGalerkinTerm(fet, geoState, X, VBar, Sigma);

  faces.computeGalerkinTerm(elems, fet, bcData, geoState, X, VBar, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeGalerkinBar_Step2(MacroCellSet **macroCells,
                                         SVec<double,1> &volRatio,
                                         SVec<double,dim> &Sigma,
                                         SVec<double,dim> &RBar,
                                         int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         int idx = nToMN[i][j];
         for (int k=0; k<dim; ++k)
               RBar[idx][k] += volRatio[idx][0] * Sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMBarAndM_Step1(DynamicVMSTerm *dvmst,
                                      SVec<double,dim> **VBar,
                                      SVec<double,1> **volRatio,
                                      SVec<double,3> &X,
                                      SVec<double,dim> &V,
                                      SVec<double,dim> &SigmaBar,
                                      SVec<double,dim> &Sigma)
{

  elems.computeMBarAndM(dvmst, VBar, volRatio, X, V, SigmaBar, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMBarAndM_Step2(MacroCellSet **macroCells,
                                      SVec<double,1> **volRatio,
                                      SVec<double,dim> &MBar,
                                      SVec<double,dim> &M,
                                      SVec<double,dim> &SigmaBar,
                                      SVec<double,dim> &Sigma,
                                      int scopeDepth1,
                                      int scopeDepth2)
{

  Connectivity &nToMN1 = *nodesToMCNodes[scopeDepth1-1];
  Connectivity &nToMN2 = *nodesToMCNodes[scopeDepth2-1];

  for (int i=0; i<nToMN1.csize(); ++i) {
    if (macroCells[scopeDepth1-1]->containing(i) != -1) {
       for (int j=0; j<nToMN1.num(i); ++j) {
         for (int k=0; k<dim; ++k)
            MBar[ nToMN1[i][j] ][k] += (*volRatio[0])[ nToMN1[i][j] ][0] * SigmaBar[i][k];
         if (i == nToMN1[i][j]) {
            for (int k=0; k<dim; ++k)
               M[ nToMN1[i][j] ][k] += (1.0 - (*volRatio[0])[ nToMN1[i][j] ][0]) * Sigma[i][k];
         }
         else {
            for (int k=0; k<dim; ++k)
               M[ nToMN1[i][j] ][k] += -1.0 * (*volRatio[0])[ nToMN1[i][j] ][0] * Sigma[i][k];
         }
       }
    }
  }

  for (int i=0; i<nToMN2.csize(); ++i) {
    if (macroCells[scopeDepth2-1]->containing(i) != -1) {
       for (int j=0; j<nToMN2.num(i); ++j) {
         for (int k=0; k<dim; ++k)
           MBar[ nToMN2[i][j] ][k] += -1.0 * (*volRatio[1])[ nToMN2[i][j] ][0] * SigmaBar[i][k];
       }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicVMSTerm_Step1(DynamicVMSTerm *dvmst,
                                            SVec<double,dim> **VBar,
                                            SVec<double,3> &X,
                                            SVec<double,dim> &V,
                                            SVec<double,dim> &Sigma,
                                            Vec<double> &CsDelSq,
                                            Vec<double> &PrT,
                                            Vec<double> *Cs,
                                            Vec<double> &Delta)
{

  elems.computeDynamicVMSTerm(dvmst, VBar, X, V, Sigma, CsDelSq, PrT, Cs, Delta);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicVMSTerm_Step2(MacroCellSet **macroCells,
                                            SVec<double,1> **volRatio,
                                            SVec<double,dim> &Sigma,
                                            SVec<double,dim> &R,
                                            int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         if (i == nToMN[i][j]) {
            for (int k=0; k<dim; ++k)
               R[ nToMN[i][j] ][k] += (1.0 - (*volRatio[0])[ nToMN[i][j] ][0]) * Sigma[i][k];
         }
         else {
            for (int k=0; k<dim; ++k)
               R[ nToMN[i][j] ][k] += -1.0 * (*volRatio[0])[ nToMN[i][j] ][0] * Sigma[i][k];
         }
       }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computedWBar_dt(MacroCellSet **macroCells,
                                SVec<double,1> **volRatio,
                                SVec<double,dim> &Sigma,
                                SVec<double,dim> &dWBardt,
                                int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         for (int k=0; k<dim; ++k)
               dWBardt[ nToMN[i][j] ][k] += (*volRatio[0])[ nToMN[i][j] ][0] * Sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeCsDeltaSq(SVec<double,dim> &R,
                                 SVec<double,dim> &RBar,
                                 SVec<double,dim> &M,
                                 SVec<double,dim> &MBar,
                                 SVec<double,dim> &dWdt,
                                 SVec<double,dim> &dWBardt,
                                 Vec<double> &CsDeltaSq,
                                 Vec<double> &PrT,
                                 int method)
{

 for (int i=0; i<nodes.size(); ++i) {

    switch (method) {

    case (0):     // Variational Germano Identity
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += (M[i][j] - MBar[i][j])*((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]));
        denom += (M[i][j] - MBar[i][j])*(M[i][j] - MBar[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      num = CsDeltaSq[i] * (M[i][4] - MBar[i][4]);
      denom = (RBar[i][4] + dWBardt[i][4]) - (R[i][4] + dWdt[i][4]);
      if (fabs(denom) < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    case(1):     // Full Least Squares
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += -(MBar[i][j] * (RBar[i][j] + dWBardt[i][j]))
               -(M[i][j] * (R[i][j] + dWdt[i][j]));
        denom += (MBar[i][j] * MBar[i][j]) + (M[i][j] * M[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      num = -CsDeltaSq[i] * ((R[i][4] + dWdt[i][4])*M[i][4]
                           + (RBar[i][4] + dWBardt[i][4])*MBar[i][4]);
      denom = pow((RBar[i][4] + dWBardt[i][4]),2.0) + pow((R[i][4] + dWdt[i][4]),2.0);
      if (fabs(denom) < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    case(2):   // Special Clipping Procedure
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += (M[i][j] - MBar[i][j])*((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]));
        denom += (M[i][j] - MBar[i][j])*(M[i][j] - MBar[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      if (CsDeltaSq[i] < 0.0) {
        num = 0.0;
        denom = 0.0;
        for(int j=1; j<4; ++j) {
          num += pow((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]) , 2.0);
          denom += pow((M[i][j]) - MBar[i][j], 2.0);
        }
        if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
        else CsDeltaSq[i] = sqrt(num / denom);
      }

      num = CsDeltaSq[i] * (M[i][4] - MBar[i][4]);
      denom = (RBar[i][4] + dWBardt[i][4]) - (R[i][4] + dWdt[i][4]);
      if (denom < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    default:
      fprintf(stderr,"Error :: Method to Solve the Residual Equation is Not Correct...Aborting !!\n");
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuSmag(SmagorinskyLESTerm *smag, SVec<double,3> &X,
                                  SVec<double,dim> &V, Vec<double> &mutOmu)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double mut = smag->computeMutOMu(vol, dp1dxj, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//--------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuVMS(VMSLESTerm *vmst, SVec<double,dim> &VBar, SVec<double,3> &X,
                                 SVec<double,dim> &V, Vec<double> &mutOmu)
{

   for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
     double dp1dxj[4][3];
     double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
     double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
     double *vbar[4] = {VBar[elems[tetNum][0]], VBar[elems[tetNum][1]],
                        VBar[elems[tetNum][2]], VBar[elems[tetNum][3]]};

     double mut = vmst->computeMutOMu(vol, dp1dxj, vbar, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuDynamicVMS(DynamicVMSTerm *dvmst, SVec<double,dim> &VBar, SVec<double,3> &X,
                                        SVec<double,dim> &V, Vec<double> &Cs, Vec<double> &mutOmu)
{

   for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
     double dp1dxj[4][3];
     double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
     double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
     double *vbar[4] = {VBar[elems[tetNum][0]], VBar[elems[tetNum][1]],
                        VBar[elems[tetNum][2]], VBar[elems[tetNum][3]]};
     double cs[4] = {Cs[elems[tetNum][0]], Cs[elems[tetNum][1]],
                     Cs[elems[tetNum][2]], Cs[elems[tetNum][3]]};

     double mut = dvmst->computeMutOMu(vol, dp1dxj, vbar, v, cs, X, elems[tetNum]);
     for (int i=0; i<4; ++i)
       mutOmu[elems[tetNum][i]] += mut * vol;

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuWale(WaleLESTerm *wale, SVec<double,3> &X,
                                  SVec<double,dim> &V, Vec<double> &mutOmu)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double mut = wale->computeMutOMu(vol, dp1dxj, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuDynamicLES(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                  SVec<double,3> &X, SVec<double,dim> &V,
				  Vec<double> &mutOmu)

{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double cs[4] = {Cs[elems[tetNum][0]][0], Cs[elems[tetNum][1]][0],
                    Cs[elems[tetNum][2]][0], Cs[elems[tetNum][3]][0]};
    double pt[4] = {Cs[elems[tetNum][0]][1], Cs[elems[tetNum][1]][1],
                    Cs[elems[tetNum][2]][1], Cs[elems[tetNum][3]][1]};

    double mut = dles->computeMutOMu(vol, dp1dxj, v, cs, pt, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//--------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
					    GeoState &geoState, SVec<double,3> &X,
					    Vec<double> &ctrlVol, SVec<double,dim> &V,
					    GenMat<Scalar,neq> &A,
					    Vec<GhostPoint<dim>*>* ghostPoints,LevelSetStructure *LSS)
{

  elems.computeJacobianGalerkinTerm(fet, geoState, X, ctrlVol, V, A, ghostPoints,LSS);

  faces.computeJacobianGalerkinTerm(elems, fet, bcData, geoState, X, ctrlVol, V, A);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeBCsJacobianWallValues(FemEquationTerm *fet, BcData<dim> &bcData,
					    GeoState &geoState, SVec<double,3> &X,
					    SVec<double,dim> &V)
{

  faces.computeBCsJacobianWallValues(elems, fet, bcData, geoState, X, V);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianVolumicForceTerm(VolumicForceTerm *volForce,
                                           Vec<double> &ctrlVol, SVec<double,dim> &V,
                                           GenMat<Scalar,neq> &A)
{
  /* computation of the jacobian part due to gravity
   * There is only a block diagonal term to compute for each node.
   */
  if (neq>=5){
    Scalar *Aii;
    double jac[neq*neq];
    for (int i=0; i<ctrlVol.size(); i++){
      Aii = A.getElem_ii(i);
      volForce->computeJacobianVolumeTerm(neq, ctrlVol[i], V[i], jac);
      for (int m=0; m<neq*neq; m++)
        Aii[m] += jac[m];

    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::getExtrapolationValue(Extrapolation<dim>* xpol,SVec<double,dim> &V, SVec<double,dim> &Ubc,
				      VarFcn *vf, BcData<dim>& bcData, GeoState& geoState, SVec<double,3>& X)
{
        inletNodes.getExtrapolationValue(xpol, V, Ubc, vf, bcData, geoState, elems, locToGlobNodeMap, X);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyExtrapolationToSolutionVector(Extrapolation<dim>* xpol,SVec<double,dim> &U,
						   SVec<double,dim> &Ubc)
{
        inletNodes.applyExtrapolationToSolutionVector(xpol, U, Ubc, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyBCsToSolutionVector(BcFcn *bcFcn, BcData<dim> &bcData,
                                         SVec<double,dim> &U, LevelSetStructure *LSS)
{
  SVec<double,dim> &Vwall = bcData.getNodeStateVector();

  // In the case of an Embedded Simulation, we want inactive nodes to stay at initial state. Thus no BC should be apply 
  // to them. 
  bool isActive = true; 

  for (int i=0; i<nodes.size(); ++i) {
    if(LSS) isActive = LSS->isActive(0.0,i);
    if (nodeType[i] != BC_INTERNAL && isActive)
      bcFcn->applyToSolutionVector(nodeType[i], Vwall[i], U[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyBCsToResidual(BcFcn *bcFcn, BcData<dim> &bcData,
				   SVec<double,dim> &U, SVec<double,dim> &F, LevelSetStructure *LSS)
{
  SVec<double,dim> &Vwall = bcData.getNodeStateVector();

  // In the case of an Embedded Simulation, we want inactive nodes to stay at initial state. Thus no BC should be apply 
  // to them. 
  bool isActive = true; 

	if (sampleMesh) {
		int i;
		for (int iNode=0; iNode<numSampledNodes; ++iNode) {
			i = locSampleNodes[iNode];
			if (nodeType[i] != BC_INTERNAL && isActive)
				bcFcn->applyToResidualTerm(nodeType[i], Vwall[i], U[i], F[i]);
		}
	}
	else {
		for (int i=0; i<nodes.size(); ++i) {
			if (nodeType[i] != BC_INTERNAL)
				bcFcn->applyToResidualTerm(nodeType[i], Vwall[i], U[i], F[i]);
		}
	}
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::applyBCsToDerivativeOfResidual(BcFcn *bcFcn, BcData<dim> &bcData,
				   SVec<double,dim> &U, SVec<double,dim> &dU, SVec<double,dim> &dF)
{

  SVec<double,dim> &Vwall = bcData.getNodeStateVector();
  SVec<double,dim> &dVwall = bcData.getdNodeStateVector();

  for (int i=0; i<nodes.size(); ++i)
    if (nodeType[i] != BC_INTERNAL)
      bcFcn->applyToDerivativeOfResidualTerm(nodeType[i], Vwall[i], dVwall[i], U[i], dU[i], dF[i]);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToH2Jacobian(BcFcn *bcFcn, BcData<dim> &bcs,
				   SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();

  int k;
  for (int l=0; l<edges.size(); ++l) {

    if (bcMap.find(l) != bcMap.end())  {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // the Aij is the off-diagonal term for this equation
        if (Aij && Aij_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aij is the diag term for the ith eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aji);

      }
      if (nodeType[j] != BC_INTERNAL) {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]+numBcNodes[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // Aij is the diag term for the jth eq.
        if (Aij && Aij_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]+numBcNodes[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aji is the off-diag term for the jth eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);

      }
    }
  }
  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aii);
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::applyBCsToH2Jacobian(BcFcn *bcFcn, BcData<dim> &bcs,
				   SVec<double,dim> &U, GenMat<Scalar,dim> &A)
{

  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();

  int k;
  for (int l=0; l<edges.size(); ++l) {

    if (bcMap.find(l) != bcMap.end())  {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // the Aij is the off-diagonal term for this equation
        if (Aij && Aij_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aij is the diag term for the ith eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aji);

      }
      if (nodeType[j] != BC_INTERNAL) {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]+numBcNodes[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // Aij is the diag term for the jth eq.
        if (Aij && Aij_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]+numBcNodes[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aji is the off-diag term for the jth eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);

      }
    }
  }
  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aii);
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::applyBCsToProduct(BcFcn *bcFcn, BcData<dim> &bcs, SVec<double,dim> &U, SVec<Scalar,dim> &Prod)
{

  for (int i=0; i<nodes.size(); ++i)
    if (nodeType[i] != BC_INTERNAL)
        bcFcn->applyToProductTerm(nodeType[i], Prod[i]);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToJacobian(BcFcn *bcFcn, BcData<dim> &bcs,
                                   SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getElem_ij(l);
        if (Aij)
          bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
    }

    if (nodeType[j] != BC_INTERNAL) {
      Scalar *Aji = A.getElem_ji(l);
      if (Aji)
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);
    }
  }

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToDiagonalTerm(nodeType[i], Vwall[i], U[i], Aii);
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToJacobianWallValues(BcFcn *bcFcn, BcData<dim> &bcs,
                                   SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{

  SVec<double,dim> &Vwall = bcs.getNodeStateVector();
  SVec<double,dim> &dVwall = bcs.getdNodeStateVectorSA();

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToDiagonalTerm(nodeType[i], Vwall[i], dVwall[i], U[i], Aii);
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
SparseMat<Scalar,dim> *SubDomain::createMaskJacobian(int *ndType, MemoryPool *mp)
{

  nodeToNodeMaskJacobian = createElemBasedConnectivity();

  int *ia = (*nodeToNodeMaskJacobian).ptr();
  int *ja = (*nodeToNodeMaskJacobian)[0];
  int n = nodes.size();
  int nnz = ia[n];

  Scalar (*a)[dim*dim] = 0;

  if (mp)
    a = reinterpret_cast<Scalar (*)[dim*dim]>(mp->request(nnz * dim*dim * sizeof(Scalar)));

  SparseMat<Scalar,dim> *A = new SparseMat<Scalar,dim>(n, nnz, ia, ja, a, 0, ndType);

  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
MvpMat<Scalar,dim> *SubDomain::createMaskMatVecProd(bool nsFlag)
{

// Original
//  int numOffDiagEntries = 0;

// Included (MB*)
  if ((nsFlag) && (!numOffDiagEntries))  {
// Original
//  if (nsFlag)  {
    numBcNodes = new int[edges.size()];
    int (*edgePtr)[2] = edges.getPtr();

    for (int iEdge = 0; iEdge < edges.size(); iEdge++)  {
      int i = edgePtr[iEdge][0];
      int j = edgePtr[iEdge][1];

      numBcNodes[iEdge] = 0;

      if (nodeType[i] != BC_INTERNAL)   {
        bcMap[iEdge] = numOffDiagEntries;
        numOffDiagEntries++;
      }

     if (nodeType[j] != BC_INTERNAL)  {
       if (bcMap.find(iEdge) == bcMap.end())
         bcMap[iEdge] = numOffDiagEntries;
       else
         numBcNodes[iEdge] = 1;
       numOffDiagEntries++;
     }
    }
  }

// Included (MB*)
  if (nsFlag) {
    MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), numOffDiagEntries);

    return A;
  }
  else {
    MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), 0);

    return A;
  }
// Original
//  MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), numOffDiagEntries);
//
//  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DiagMat<Scalar,dim> *SubDomain::createMaskDiagonal(typename DiagMat<Scalar,dim>::Type type,
						   int *ndType)
{

  DiagMat<Scalar,dim> *A = new DiagMat<Scalar,dim>(type, nodes.size(), ndType);

  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
SparseMat<Scalar,dim> *SubDomain::createMaskILU(int fill, int renum, int *ndType)
{

  nodeToNodeMaskILU = createEdgeBasedConnectivity();

  compStruct *nodeRenum = createRenumbering(nodeToNodeMaskILU, renum, 0);

  int *ia = (*nodeToNodeMaskILU).ptr();
  int *ja = (*nodeToNodeMaskILU)[0];
  int n = nodes.size();
  int nnz = ia[n];

  SparseMat<Scalar,dim> *A = new SparseMat<Scalar,dim>(n, nnz, ia, ja, 0, nodeRenum, ndType);

  A->symbolicILU(fill);

  A->createPointers(edges);

  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::computeH1(FluxFcn **fluxFcn, BcData<dim> &bcData,
                          GeoState &geoState, Vec<double> &ctrlVol,
                          SVec<double,dim> &V, GenMat<Scalar,dim> &A)
{

  int k;

  double dfdUi[dim*dim], dfdUj[dim*dim];

  Scalar *Aii, *Ajj, *Aij, *Aji;

  // contribution of the edges

  Vec<Vec3D> &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    fluxFcn[0]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l],
                                 V[i], V[j], dfdUi, dfdUj);

    Aii = A.getElem_ii(i);
    Ajj = A.getElem_ii(j);
    Aij = A.getElem_ij(l);
    Aji = A.getElem_ji(l);

    if (edgeFlag[l])
      for (k=0; k<dim*dim; ++k) { Aii[k] += dfdUi[k]; Ajj[k] -= dfdUj[k]; }

    if (Aij && Aji) {

      double voli = 1.0 / ctrlVol[i];
      double volj = 1.0 / ctrlVol[j];
      for (k=0; k<dim*dim; ++k) {
        Aij[k] += dfdUj[k] * voli;
        Aji[k] -= dfdUi[k] * volj;
      }
    }

  }
  // contribution of the boundary faces
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

  for (int i=0; i<ctrlVol.size(); ++i) {

    double voli = 1.0 / ctrlVol[i];

    Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] *= voli;

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn, BcData<dim> &bcData,
			  GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V,
			  NodalGrad<dim> &ngrad, GenMat<Scalar,neq> &A)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim], dfdVi[dim*dim], dfdVj[dim*dim];

  Scalar *Aij, *Aji;

  // contribution of the edges

  Vec<Vec3D> &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    int k;
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vi, Vj, dfdVi, dfdVj);

    Aij = A.getElem_ij(l);
    Aji = A.getElem_ji(l);

    if (Aij && Aji)  {
      for (k=0; k<neq*neq; ++k) {
        Aij[k] += dfdVj[k];
        Aji[k] += dfdVi[k];
      }
    }

  }

  // contribution of the boundary faces
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::precomputeRec(RecFcn *recFcn, SVec<double,3> &X,
			      SVec<double,dim> &V, NodalGrad<dim> &ngrad,
			      SVec<Scalar,dim> &aij, SVec<Scalar,dim> &aji,
			      SVec<Scalar,dim> &bij, SVec<Scalar,dim> &bji)
{

  double ddVij[dim], ddVji[dim];

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->precompute(V[i], ddVij, V[j], ddVji, aij[l], aji[l], bij[l], bji[l]);
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::computeMatVecProdH1(bool *nodeFlag, GenMat<Scalar,dim> &A,
				    SVec<double,dim> &p, SVec<double,dim> &prod)
{

  int i, j, l;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  Scalar (*a)[dim*dim] = A.data();

  prod = 0.0;

#pragma ivdep
  for (i=0; i<numNodes; ++i)
    if (nodeFlag[i])
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

#pragma ivdep
  for (l=0; l<numEdges; ++l) {

    if (edgeFlag[l]) {

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l, p.v, j, prod.v, i);
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p.v, i, prod.v, j);

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::computeMatVecProdH1(bool *nodeFlag, GenMat<Scalar,dim> &A,
				    SVec<double,dim> &p, SVec<double,dim> &prod,
                                    SVec<double,dim> &ghostP, SVec<double,dim>& ghostProd)
{

  int i, j, l;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  Scalar (*a)[dim*dim] = A.data();

  prod = 0.0;

#pragma ivdep
  for (i=0; i<numNodes; ++i)
    if (nodeFlag[i])
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

#pragma ivdep
  for (l=0; l<numEdges; ++l) {

    if (edgeFlag[l]) {

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l, p.v, j, prod.v, i);
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p.v, i, prod.v, j);

    }

  }

  typename GenMat<Scalar,dim>::AuxilliaryIterator* myItr = A.begin_realNodes();
  if (myItr) {
    do { 
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, ghostP.v, myItr->col , prod.v, myItr->row);
    } while (A.next(myItr));
    A.free(myItr);
  }
  
  myItr = A.begin_ghostNodes();
  if (myItr) {
    do { 
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, p.v, myItr->col , ghostProd.v, myItr->row);
    } while (A.next(myItr));
    A.free(myItr);
  }

  myItr = A.begin_ghostGhostNodes();
  if (myItr) {
    do { 
     DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, ghostP.v, myItr->col , ghostProd.v, myItr->row);
    } while (A.next(myItr));
    A.free(myItr);
  }

}
//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2(RecFcn *recFcn, SVec<double,3> &X,
	        Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
	        SVec<double,dim> &aij, SVec<double,dim> &aji,
	        SVec<double,dim> &bij, SVec<double,dim> &bji,
	        SVec<Scalar2,dim> &p, NodalGrad<dim, Scalar2> &dpdxj,
	        SVec<Scalar2,dim> &prod) {

  int i, j, l;

  Scalar2 ddpij[dim], ddpji[dim], pij[1][dim], pji[1][dim];
  Scalar2 tmp[1][dim], tmpi[1][dim], tmpj[1][dim];

  Scalar1 (*a)[dim*dim] = A.data();

  SVec<Scalar2,dim> &dpdx = dpdxj.getX();
  SVec<Scalar2,dim> &dpdy = dpdxj.getY();
  SVec<Scalar2,dim> &dpdz = dpdxj.getZ();

  prod = (Scalar2) 0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  bool *masterFlag = edges.getMasterFlag();

  if (bcMap.size() > 0)  {
    int index;
    for (l=0; l<numEdges; ++l) {
      if (!masterFlag[l]) continue;
      i = edgePtr[l][0];
      j = edgePtr[l][1];

      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddpij[k] = dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k];
        ddpji[k] = dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k];

      }

      // result of the reconstructed-limited states are in pij, pji
      recFcn->template compute<Scalar2, dim>(p[i], ddpij, p[j], ddpji, aij[l], aji[l],
         				     bij[l], bji[l], pij[0], pji[0]);

      // A is applied to reconstructed-limited states and stored in tmpi, tmpj
      // address of a is shifted by the number of diagonal entries (numnodes)

      if (bcMap.find(l) != bcMap.end())  {
        if (nodeType[i] != BC_INTERNAL)
          index = numNodes+2*numEdges+2*bcMap[l];
        else
          index = numNodes + 2*l;

        DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, index, pji, 0, prod.v, i);
        DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, index+1, pij, 0, prod.v, i);

        if (nodeType[j] != BC_INTERNAL)
          index = numNodes+2*numEdges+2*bcMap[l]+2*numBcNodes[l];
        else
          index = numNodes+2*l;

        DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndSubToVector(a, index, pji, 0, prod.v, j);
        DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndSubToVector(a, index+1, pij, 0, prod.v, j);

      }
      else  {
        DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l, pji, 0, tmpi, 0);
        DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l + 1, pij, 0, tmpj, 0);

        VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);
        VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i);
        VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j);
      }
    }
  }
  else  {

    for (l=0; l<numEdges; ++l) {
      if (!masterFlag[l]) continue;
      i = edgePtr[l][0];
      j = edgePtr[l][1];

      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      for (int k=0; k<dim; ++k) {
        ddpij[k] = dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k];
        ddpji[k] = dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k];

      }

      // result of the reconstructed-limited states are in pij, pji
      recFcn->template compute<Scalar2, dim>(p[i], ddpij, p[j], ddpji, aij[l], aji[l],
                                             bij[l], bji[l], pij[0], pji[0]);

      DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l, pji, 0, tmpi, 0);
      DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l + 1, pij, 0, tmpj, 0);

      VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);
      VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i); 
      VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j); 
    }
  }

  // contribution from diagonal entries of A
  for (i=0; i<numNodes; ++i) {

    DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

    double voli = 1.0 / ctrlVol[i];
    for (int k=0; k<dim; ++k) prod[i][k] *= voli;

  }
}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2T(RecFcn *recFcn, SVec<double,3> &X,
                Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                SVec<double,dim> &aij, SVec<double,dim> &aji,
                SVec<double,dim> &bij, SVec<double,dim> &bji,
                SVec<Scalar2,dim> &p, SVec<Scalar2,dim> &zu,
                SVec<Scalar2,dim> &zgx, SVec<Scalar2,dim> &zgy,
                SVec<Scalar2,dim> &zgz) {
  int i, j, l, k;

  // Same notations as in the Fortran code
  Scalar2 zu_is1[1][dim], zu_is2[1][dim], zgx_is1[1][dim], zgx_is2[1][dim];
  Scalar2 zgy_is1[1][dim], zgy_is2[1][dim], zgz_is1[1][dim], zgz_is2[1][dim];
  //Scalar2 tmp[1][dim];
  Scalar2 tmpi[1][dim], tmpj[1][dim];
  Scalar2 tmp1[1][dim];
  Scalar1 (*a)[dim*dim] = A.data();


  zu = (Scalar2)0.0;
  zgx = (Scalar2)0.0;
  zgy = (Scalar2)0.0;
  zgz = (Scalar2)0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  bool *masterFlag = edges.getMasterFlag();
  for (l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

    for (k=0; k<dim; ++k)
      tmp1[0][k] = p[i][k] - p[j][k];

    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransToVector(a, numNodes + 2*l, tmp1, 0, tmpi, 0);
    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, tmp1, 0, tmpj, 0);
    // These routines do not exist anymore
    //denseMatrixTransTimesVector(a, numNodes + 2*l, tmp1, 0, tmpi, 0);
    //denseMatrixTransTimesVector(a, numNodes + 2*l + 1, tmp1, 0, tmpj, 0);

    recFcn->template computeT<Scalar2, dim> (dx, tmpi[0], tmpj[0], aij[l], aji[l], bij[l], bji[l], i,
                j, zu_is1[0], zu_is2[0], zgx_is1[0], zgx_is2[0], zgy_is1[0], zgy_is2[0], zgz_is1[0], zgz_is2[0]);


    VectorOp<Scalar2,dim>::add(zu_is1, 0, zu.v, i);
    VectorOp<Scalar2,dim>::add(zu_is2, 0, zu.v, j);
    VectorOp<Scalar2,dim>::add(zgx_is1, 0, zgx.v, i);
    VectorOp<Scalar2,dim>::add(zgx_is2, 0, zgx.v, j);
    VectorOp<Scalar2,dim>::add(zgy_is1, 0, zgy.v, i);
    VectorOp<Scalar2,dim>::add(zgy_is2, 0, zgy.v, j);
    VectorOp<Scalar2,dim>::add(zgz_is1, 0, zgz.v, i);
    VectorOp<Scalar2,dim>::add(zgz_is2, 0, zgz.v, j);


/*
    // These routines do not exist anymore
    addVector(p1, 0, prod2.v, i);
    addVector(p2, 0, prod2.v, j);
    addVector(p3, 0, prod.v, i);
    addVector(p4, 0, prod.v, j);
    addVector(p5, 0, prod3.v, i);
    addVector(p6, 0, prod3.v, j);
    addVector(p7, 0, prod4.v, i);
    addVector(p8, 0, prod4.v, j);
*/
  }
  for (i=0; i<numNodes; ++i)
    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransAndAddToVector(a, i, p.v, i, zu.v, i);

    // This routine does not exist anymore
    //addDenseMatrixTransTimesVector(a, i, p.v, i, prod2.v, i);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2Tb(RecFcn *recFcn, SVec<double,3> &X,
                Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                NodalGrad<dim, Scalar2> &dpdxj, SVec<Scalar2,dim> &p,
                SVec<Scalar2,dim> &prod, SVec<Scalar2,dim> &zu)  {

  int i, l;
  Scalar2 p1[1][dim];

  SVec<Scalar2,dim> &dpdx = dpdxj.getX();
  SVec<Scalar2,dim> &dpdy = dpdxj.getY();
  SVec<Scalar2,dim> &dpdz = dpdxj.getZ();

  prod = (Scalar2) 0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<numNodes; ++l) {

    i = l;
    recFcn->computeTb(zu, dpdx, dpdy, dpdz, i, p1[0]);
    VectorOp<Scalar2,dim>::add(p1, 0, prod.v, i);
    //addVector(p1, 0, prod.v, i);
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenNodes(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], sharedNodes->num(iSub)*dim);

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenInletNodes(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], sharedInletNodes->num(iSub)*dim);

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenEdges(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], numSharedEdges[iSub]*dim);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      for (int j = 0; j < dim; ++j)
	buffer[iNode][j] = w[ (*sharedNodes)[iSub][iNode] ][j];
    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)  {	// KTC: COULD DO ONLY FOR ONE LAYER OF NODES AWAY
	w[ (*sharedNodes)[iSub][iNode] ][j] += buffer[iNode][j];
      }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndGhostStates(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{

  double *v;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  v = ghostPoints[ (*sharedNodes)[iSub][iNode] ]->getPrimitiveState();
	  for (int j = 0; j < dim; ++j) buffer[iNode][j] = v[j];
	}
      else 
	{
	  // It can happened that all the edges containing the ghostPoint are not master in this SubDomain.
	  // In which case we do not want to send garbage.
	  for (int j = 0; j < dim; ++j) buffer[iNode][j] = 0.0;
	}
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndGhostWeights(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    double (*buffer)[1] = reinterpret_cast<double (*)[1]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  buffer[iNode][0] = double(ghostPoints[ (*sharedNodes)[iSub][iNode] ]->ng);
	}
      else 
	{
	  buffer[iNode][0] = 0.0;
	}
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvGhostStates(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{
  GhostPoint<dim> *gp;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = sp.recData(rcvChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {

      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  gp = ghostPoints[ (*sharedNodes)[iSub][iNode] ];
	  for (int j = 0; j < dim; ++j)  
	    {
	      gp->Vg[j] += buffer[iNode][j];
	    }
	}
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvGhostWeights(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = sp.recData(rcvChannel[iSub]);
    double (*buffer)[1] = reinterpret_cast<double (*)[1]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      if( buffer[iNode][0] < 1e-14) continue; // if the weight is zero, the whole buffered state is gonna be zero.
      if(!ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  ghostPoints[ (*sharedNodes)[iSub][iNode] ] = new GhostPoint<dim>;
	}  
      ghostPoints[ (*sharedNodes)[iSub][iNode] ]->ng += int(buffer[iNode][0]);    
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class OpType>
void SubDomain::operateRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim], const OpType &oper)
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)  {
        w[ (*sharedNodes)[iSub][iNode] ][j] = OpType::apply(w[ (*sharedNodes)[iSub][iNode]][j], buffer[iNode][j]);
      }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndInletData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        buffer[iNode][j] = w[ (*sharedInletNodes)[iSub][iNode] ][j];
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvInletData(CommPattern<Scalar> &sp, Scalar (*w)[dim], bool ForExtrapolation)
{
  if(!ForExtrapolation){
    /* the values in w are accessed via the inlet node number
     * in the inletNode number chart.
     * This communication of values is used for
     * computation of normals at inlet nodes
     * and computation of number of faces surrounding
     * a given inlet node.
     */

    for (int iSub = 0; iSub < numNeighb; ++iSub) {
      SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
      Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode)
        for (int j = 0; j < dim; ++j)
          w[ (*sharedInletNodes)[iSub][iNode] ][j] += buffer[iNode][j];

    }
  }else{
    /* the values in w are accessed via the inlet node number
     * in the inletNode number chart.
     * This communication of values is used for
     * computation of extrapolation values
     * in the whole domain.
     */

    for (int iSub = 0; iSub < numNeighb; ++iSub) {
      SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
      Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
      int inletnode;
      for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode){
        inletnode = (*sharedInletNodes)[iSub][iNode] ;
        if (w[inletnode][0] == 0.0){
          for (int j = 0; j < dim; ++j)
            w[inletnode][j] += buffer[iNode][j];
        }
        else{
          //if(buffer[iNode][0]!=0.0)
          //      fprintf(stderr, "two shared inlet nodes have their own extrapolation and their densities are %f vs %f\n", w[ (*sharedInletNodes)[iSub][iNode] ][0], buffer[iNode][0]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndInletRhsData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    int inletnode, node;

    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode) {
      inletnode = (*sharedInletNodes)[iSub][iNode] ;
      node = inletNodes[inletnode].getNodeNum();
      //fprintf(stdout, "node in two domains %d for subd %d\n", node, locSubNum);
      for (int j = 0; j < dim; ++j)
        buffer[iNode][j] = w[ node ][j];
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvInletRhsData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
    /* the values in w are accessed via the inlet node number
     * in the subDomain number chart.
     * This communication of values is used for
     * the recomputation of the rhs term, which
     * contains the flux (computeFVT = computeResidual)
     */

        // we want to pass the data of a shared inlet node from one subdomain to a
        // neighbouring subdomain. Only one of the data has a physical value, while
        // the other ones are set to 0.0, except in one case where the intersection
        // between the normal of the inlet node and the opposite face of the tetrahedron
        // gives one of the edges of the tetrahedra.
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    int inletnode, node;
    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode){
      inletnode = (*sharedInletNodes)[iSub][iNode] ;
      node = inletNodes[inletnode].getNodeNum();
      if (w[node][0] == 0.0){
        for (int j = 0; j < dim; ++j)
          w[node][j] += buffer[iNode][j];
      }
      else{
        //if(buffer[iNode][0]!=0.0){
          // fprintf(stderr, "two shared inlet nodes (inlet node %d and node %d) have their own extrapolation\n     and their densitiesRHS are %.14e vs %.14e\n", inletnode, locToGlobNodeMap[node]+1, w[ node ][0], buffer[iNode][0]);
        //}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::minRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        if (buffer[iNode][j] < w[ (*sharedNodes)[iSub][iNode] ][j])
          w[ (*sharedNodes)[iSub][iNode] ][j] = buffer[iNode][j];

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::maxRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        if (buffer[iNode][j] > w[ (*sharedNodes)[iSub][iNode] ][j])
          w[ (*sharedNodes)[iSub][iNode] ][j] = buffer[iNode][j];
  }

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim1, int dim2>
void SubDomain::TagPsiExchangeData(CommPattern<Scalar1> &splevel, Scalar1 (*level)[dim1],
                                   CommPattern<Scalar2> &sppsi, Scalar2 (*psi)[dim2])
{
  /* it is assumed that dim1 = 1, dim2 = 1 */

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar1> sInfolevel = splevel.recData(rcvChannel[iSub]);
    Scalar1 (*blevel)[dim1] = reinterpret_cast<Scalar1 (*)[dim1]>(sInfolevel.data);
    SubRecInfo<Scalar2> sInfopsi = sppsi.recData(rcvChannel[iSub]);
    Scalar2 (*bpsi)[dim2] = reinterpret_cast<Scalar2 (*)[dim2]>(sInfopsi.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      if (level[ (*sharedNodes)[iSub][iNode] ][0] == 0 &&
          blevel[iNode][0] > 0){
        level[ (*sharedNodes)[iSub][iNode] ][0] = blevel[iNode][0];
        psi[ (*sharedNodes)[iSub][iNode] ][0]   = bpsi[iNode][0];
      }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {

      Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);

      for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {

      Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);

      for (int j=0; j<dim*dim; ++j) a[j] += buffer[iNode][j];

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndDiagInletBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{
  int node;
  int type;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      node = (*sharedNodes)[iSub][iNode];
      type = nodeType[node];
      if (!(type == BC_INLET_MOVING || type == BC_OUTLET_MOVING ||
            type == BC_INLET_FIXED  || type == BC_OUTLET_FIXED) ){

        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
  }


}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvDiagInletBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{
  int node;
  int type;

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      node = (*sharedNodes)[iSub][iNode];
      type = nodeType[node];
      if (!(type == BC_INLET_MOVING || type == BC_OUTLET_MOVING ||
            type == BC_INLET_FIXED  || type == BC_OUTLET_FIXED) ){

        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j=0; j<dim*dim; ++j) a[j] += buffer[iNode][j];
      }else{
        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        if (a[0] == 0.0)
          for (int k=0; k<dim*dim; k++) a[k] += buffer[iNode][k];
      }

    }

  }

}
//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndEdgeData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
      for (int k=0; k<dim; ++k)
	buffer[iEdge][k] = w[ sharedEdges[iSub][iEdge].edgeNum ][k];

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvEdgeData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
      for (int k=0; k<dim; ++k)
	w[ sharedEdges[iSub][iEdge].edgeNum ][k] += buffer[iEdge][k];

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.getElem_ij(edgeNum);
	aji = A.getElem_ji(edgeNum);
      }
      else {
	aij = A.getElem_ji(edgeNum);
	aji = A.getElem_ij(edgeNum);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = aij[k];
	  buffer[iEdge][1][k] = aji[k];
	}
      }
      else {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = 0.0;
	  buffer[iEdge][1][k] = 0.0;
	}
      }

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.getElem_ij(edgeNum);
	aji = A.getElem_ji(edgeNum);
      }
      else {
	aij = A.getElem_ji(edgeNum);
	aji = A.getElem_ij(edgeNum);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  aij[k] += buffer[iEdge][0][k];
	  aji[k] += buffer[iEdge][1][k];
	}
      }

    }

  }

}

template<class Scalar, int dim>
void SubDomain::sndGhostOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);
    int (*ptr)[2] = edges.getPtr();

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      int i = ptr[edgeNum][0],j = ptr[edgeNum][1];

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.queryRealNodeElem_ij(i,j);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(i,j);
	aji = A.queryRealNodeElem_ij(j,i);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(j,i);
      }
      else {
	aji = A.queryRealNodeElem_ij(i,j);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(i,j);
	aij = A.queryRealNodeElem_ij(j,i);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(j,i);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = aij[k];
	  buffer[iEdge][1][k] = aji[k];
	}
      }
      else {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = 0.0;
	  buffer[iEdge][1][k] = 0.0;
	}
      }

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvGhostOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    int (*ptr)[2] = edges.getPtr();
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      
      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      int i = ptr[edgeNum][0],j = ptr[edgeNum][1];

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.queryRealNodeElem_ij(i,j);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(i,j);
	aji = A.queryRealNodeElem_ij(j,i);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(j,i);
      }
      else {
	aji = A.queryRealNodeElem_ij(i,j);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(i,j);
	aij = A.queryRealNodeElem_ij(j,i);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(j,i);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  aij[k] += buffer[iEdge][0][k];
	  aji[k] += buffer[iEdge][1][k];
	}
      }
    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
double SubDomain::readTagFromFile(const char *prefix, int no, int *neq, int *nsol)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

  BinFileHandler file(name, "rb");

  int info[3];
  file.read(info, 3);

  if (info[0] != numClusNodes) {
    fprintf(stderr, "*** Error: mismatch in size for \'%s\' (%d vs %d)\n",
	    name, info[0], numClusNodes);
    exit(1);
  }

  *neq = info[1];
  *nsol = info[2];
  double tag = 0.0;

  if (no < *nsol) {
    file.seek(3*sizeof(int) + no*(sizeof(double) + numClusNodes*(*neq)*sizeof(Scalar)));
    file.read(&tag, 1);
  }

  return tag;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::openFileForWriting(const char *prefix, int no)
{

  if (clusSubNum != 0) return;

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

#ifdef SYNCHRO_WRITE
  const char *flag = "ws+";
  if (no == 0) flag = "ws";
#else
  const char *flag = "w+";
  if (no == 0) flag = "w";
#endif

  BinFileHandler file(name, flag);

  int info[3] = {numClusNodes, dim, no + 1};
  file.write(info, 3);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::writeTagToFile(const char *prefix, int no, double tag)
{

  if (clusSubNum != 0) return;

  BinFileHandler::OffType unit = dim * sizeof(Scalar);

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

#ifdef SYNCHRO_WRITE
  BinFileHandler file(name, "ws+");
#else
  BinFileHandler file(name, "w+");
#endif

  file.seek(3*sizeof(int) + no*(sizeof(double) + numClusNodes*unit));
  file.write(&tag, 1);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::readVectorFromFile(const char *prefix, int no, int neq,
				   SVec<Scalar,dim> &U, Scalar* scale)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

  BinFileHandler file(name, "rb");

  BinFileHandler::OffType unit = neq * sizeof(Scalar);

  BinFileHandler::OffType pos = 3*sizeof(int) +
    no*(sizeof(double) + numClusNodes*unit) + sizeof(double);

  Scalar *data;

  if (neq == dim)
    data = reinterpret_cast<Scalar *>(U.data());
  else
    data = new Scalar[U.size()*neq];

  int i, count = 0;
  for (i=0; i<numNodeRanges; ++i) {
    file.seek(pos + nodeRanges[i][1]*unit);
    file.read(data + count*neq, nodeRanges[i][0]*neq);
    count += nodeRanges[i][0];
  }

  if (scale) {
    for (i=0; i<U.size()*neq; ++i)
      data[i] *= *scale;
  }

  if (neq != dim) {
    int minsize = min(neq, dim);
    for (i=0; i<U.size(); ++i) {
      for (int j=0; j<minsize; ++j)
	U[i][j] = data[neq*i + j];
    }
    delete [] data;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::writeVectorToFile(const char *prefix, int no,
				  SVec<Scalar,dim> &U, Scalar* scale)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);
#ifdef SYNCHRO_WRITE
  BinFileHandler file(name, "ws+");
#else
  BinFileHandler file(name, "w+");
#endif

  BinFileHandler::OffType unit = dim * sizeof(Scalar);
  BinFileHandler::OffType pos = 3*sizeof(int) +
    no*(sizeof(double) + numClusNodes*unit) + sizeof(double);

  Scalar (*data)[dim];
  if (scale) {
    data = new Scalar[U.size()][dim];
    Scalar* v = reinterpret_cast<Scalar*>(data);
    Scalar* u = reinterpret_cast<Scalar*>(U.data());
    for (int i=0; i<U.size()*dim; ++i)
      v[i] = (*scale) * u[i];
  }
  else
    data = U.data();

  int count = 0;
  for (int i=0; i<numNodeRanges; ++i) {
    if (nodeRanges[i][2]) {
      file.seek(pos + nodeRanges[i][1]*unit);
      file.write(reinterpret_cast<Scalar *>(data + count), nodeRanges[i][0]*dim);
    }
    count += nodeRanges[i][0];
  }

  if (scale)
    delete [] data;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout,
					SVec<double,dim> &U, SVec<double,dim> &Uinlet)
{
  int node;
  for (int j=0; j<inletNodes.size(); j++){
    node = inletNodes[j].getNodeNum();
    inletNodes[j].template assignFreeStreamValues<dim>(nodeType[node],
				      Uin[node], Uout[node], Uinlet[j]);
  }

  for (int i=0; i<faces.size(); ++i) 
    faces[i].template assignFreeStreamValues2<dim>(Uin, Uout, U[i]);
  
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::assignFreeStreamValues(double *Uin, double *Uout, SVec<double,dim> &U, SVec<double,dim> &Uinlet)
{

  for (int j=0; j<inletNodes.size(); j++)
    inletNodes[j].template assignFreeStreamValues<dim>(nodeType[inletNodes[j].getNodeNum()], Uin, Uout, Uinlet[j]);

  for (int i=0; i<faces.size(); ++i)
    faces[i].template assignFreeStreamValues<dim>(Uin, Uout, U[i]);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setNodeBcValue(double* Vin, SVec<double,dim>& Unode)
{

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] == BC_INLET_MOVING || nodeType[i] == BC_INLET_FIXED) {
      Unode[i][0] = Vin[0];
      Unode[i][1] = Vin[1];
      Unode[i][2] = Vin[2];
      Unode[i][3] = Vin[3];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFaceBcValue(SVec<double,dim> &Unode, SVec<double,dim> &Uface)
{

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeFaceBcValue(Unode, Uface[i]);

}

//------------------------------------------------------------------------------
// compute node values of the last (dim2-1) face values of Uface

template<int dim1, int dim2>
void SubDomain::computeNodeBcValue(SVec<double,3> &X, SVec<double,dim1> &Uface,
				   SVec<double,dim2> &Unode)
{

  Unode = 0.0;

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<faces.getNumSampledFaces(); ++iFace) {
			i = faces.facesConnectedToSampleNode[iFace];
			faces[i].template computeNodeBcValue<dim1,dim2>(X, Uface[i], Unode);
		}

	}
	else {
		for (int i=0; i<faces.size(); ++i)
			faces[i].template computeNodeBcValue<dim1,dim2>(X, Uface[i], Unode);
	}

}

//------------------------------------------------------------------------------
// compute node values of the last (dim2-1) face values of Uface

// Included (MB)
template<int dim1, int dim2>
void SubDomain::computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim1> &Uface, SVec<double,dim1> &dUface,
				   SVec<double,dim2> &dUnode)
{

  dUnode = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].template computeDerivativeOfNodeBcValue<dim1,dim2>(X, dX, Uface[i], dUface[i], dUnode);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, SVec<double,dim> &dUfaceSA, SVec<double,dim> &dUnodeSA)
{

  dUnodeSA = 0.0;
  dNormSA = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodeBCsWallValues(X, dNormSA, dUfaceSA[i], dUnodeSA);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodalForce(PostFcn *postFcn, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X,
				  SVec<double,dim> &V, Vec<double> &Pin,
				  SVec<double,3> &F)
{

  F = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalForce(elems, postFcn, X, d2wall, Vwall[i], V, Pin[i], F, gradP);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfNodalForce(PostFcn *postFcn, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				  SVec<double,dim> &V, SVec<double,dim> &dV, Vec<double> &Pin,
				  double dS[3], SVec<double,3> &dF)
{

  dF = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();

// Remark: Maybe it is necessary to transform conservative in primitive variables.

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, Pin[i], dS, dF, gradP, dGradP);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodalHeatPower(PostFcn* postFcn, BcData<dim>& bcData,
				      GeoState& geoState, SVec<double,3>& X,
				      SVec<double,dim>& V, Vec<double>& P)
{

  P = 0.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalHeatPower(elems, postFcn, X, d2wall, Vwall[i], V, P);

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::computeNodalHeatFluxRelatedValues(PostFcn* postFcn, BcData<dim>& bcData,
                                      GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, Vec<double>& P, Vec<double>& N, bool includeKappa)
{

  P = 0.0;
  N = -1.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalHeatFluxRelatedValues(elems, postFcn, X, d2wall, Vwall[i], V, P, N, includeKappa);

}

//------------------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfNodalHeatPower(PostFcn* postFcn, BcData<dim>& bcData,
				      GeoState& geoState, SVec<double,3>& X, SVec<double,3>& dX,
				      SVec<double,dim>& V, SVec<double,dim>& dV, double dS[3], Vec<double>& dP)
{

  dP = 0.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();
  SVec<double,dim>& dVwall = bcData.getdFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dS, dP);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceAndMoment(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,3> &X,
				      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
				      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
                                      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);
    }
  }

}

//------------------------------------------------------------------------------
// KW: FS Riemann based force calculation.
template<int dim>
void SubDomain::computeForceAndMoment(ExactRiemannSolver<dim> &riemann, VarFcn *varFcn, 
                                      map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,3> &X,
                                      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
                                      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
                                      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(riemann, varFcn, n, ndot, elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);
    }
  }

}

/*
template<int dim>
void SubDomain::computeLiftSurfaces(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,3> &X,
				      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
				      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
                                      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);
    }
  }

}
*/
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::computeHeatFluxes(map<int,int> & surfOutMapHF, PostFcn* postFcn, BcData<dim>& bcData,
                                      GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, double* HF)
{
  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i){
    int idx;
    map<int,int>::iterator it = surfOutMapHF.find(faces[i].getSurfaceID());
    if(it != surfOutMapHF.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING)  
        idx = 0;
      else
        idx = -1;
    }
    if(idx >= 0)  {
   double hp = faces[i].computeHeatFluxes(elems, postFcn, X, d2wall, Vwall[i], V);
    HF[idx] += hp;
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfForceAndMoment(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
				                  GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
			                          SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
				                  Vec3D &x0, Vec3D *dFi, Vec3D *dMi, Vec3D *dFv, Vec3D *dMv, int hydro)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
	 faces[i].getCode() == BC_SLIP_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)
      faces[i].computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dS, x0, dFi[idx], dMi[idx], dFv[idx], dMv[idx], gradP, dGradP, hydro);

  }

}

//------------------------------------------------------------------------------

template<int dim>
double SubDomain::computeInterfaceWork(PostFcn* postFcn, BcData<dim>& bcData,
				       GeoState& geoState, SVec<double,3>& X,
				       SVec<double,dim>& V, Vec<double>& Pin)
{

  Vec<double>& ndot = geoState.getFaceNormalVel();
  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  double E = 0.0;
  for (int i=0; i<faces.size(); ++i)
    E += faces[i].computeInterfaceWork(elems, postFcn, X, d2wall, ndot[i], Vwall[i], V, Pin[i]);

  return E;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFaceScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
					  BcData<dim> &bcData, GeoState &geoState,
					  SVec<double,3> &X, SVec<double,dim> &V,
					  SVec<double,2> &Q)
{

  Q = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall[i], V, Q);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
					  SVec<double,dim> &V, SVec<double,3> &X,
					  Vec<double> &Q)
{
  double phi = 1.0;
  int fluidId = 0;
  for (int i=0; i<Q.size(); ++i)
    Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i],fluidId,&phi);
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void
SubDomain::computeDerivativeOfNodeScalarQuantity(PostFcn::ScalarDerivativeType type, PostFcn *postFcn, double dS[3], SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &dQ)
{

  for (int i=0; i<dQ.size(); ++i)
    dQ[i] = postFcn->computeDerivativeOfNodeScalarQuantity(type, dS, V[i], dV[i], X[i], dX[i], 1.0);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeXP(PostFcn *postFcn, SVec<double,dim> &V, SVec<double,3> &X, Vec<double> &Q, int dir)
{
  double phi = 1.0;
  int fluidId = 0;
  for (int i=0; i<Q.size(); ++i) {
    if (nodeType[i] == BC_ADIABATIC_WALL_MOVING  || BC_ISOTHERMAL_WALL_MOVING)  {
      Q[i] = postFcn->computeNodeScalarQuantity(PostFcn::DIFFPRESSURE, V[i], X[i], &phi, fluidId);
      Q[i] *= X[i][dir];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
inline void SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
                                         SVec<double,dim> &V, SVec<double,3> &X,
                                          Vec<double> &Q, Vec<int> &fluidId,SVec<double,dimLS>* phi) 
{
  
  if (phi) {
    for (int i=0; i<Q.size(); ++i)
      Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],(*phi)[i]);
  } else { 
    for (int i=0; i<Q.size(); ++i)
      Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],NULL);
  }
}

template<int dim, int dimLS>
inline double SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
						   SVec<double,dim> &V, SVec<double,3> &X,
						   Vec<int> &fluidId,int i,SVec<double,dimLS>* phi) 
{
  
  if (phi) {
    return postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],(*phi)[i]);
  } else { 
    return postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],NULL);
  }
}

//------------------------------------------------------------------------------

template<class S1, class S2>
void SubDomain::computeStiffAndForce(DefoMeshMotionData::Element typeElement,
				     SVec<double,3>& X, SVec<double,3>& F, GenMat<S1,3>& K, GenMat<S2,3>* P,
                                     double volStiff, int* ndType=0)
{
  const int MaxSize = (3*Elem::MaxNumNd);
  double kEl[MaxSize*MaxSize];
  double fEl[MaxSize];

  F = 0.0;
  K = 0.0;
  if (P) *P = 0.0;

  int i;
  for (i=0; i<elems.size(); i++)  {
    int j;
    double *fEl_loc;

    // Compute stiffness depending on type of structural analogy
    switch (typeElement) {

    case DefoMeshMotionData::LINEAR_FE : {
      elems[i].computeStiffAndForceLIN(kEl, X, nodes);
      break;
    }

    case DefoMeshMotionData::NON_LINEAR_FE : {
      elems[i].computeStiffAndForce(fEl, kEl, X, nodes, volStiff);
      for (j=0, fEl_loc = fEl;
	   j<elems[i].numNodes();
	   j++, fEl_loc+=3) {
	F[ elems[i][j] ][0] -= fEl_loc[0];
	F[ elems[i][j] ][1] -= fEl_loc[1];
	F[ elems[i][j] ][2] -= fEl_loc[2];
      }
      break;
    }

    case DefoMeshMotionData::TORSIONAL_SPRINGS : {
      elems[i].computeStiffTorsionSpring(kEl, X, volStiff);
      break;
    }

    case DefoMeshMotionData::BALL_VERTEX : {
      elems[i].computeStiffBallVertex(kEl, X, nodes, volStiff);
      break;
    }
      
    case DefoMeshMotionData::NL_BALL_VERTEX : {
      elems[i].computeStiffAndForceBallVertex(fEl, kEl, X, nodes, volStiff);
      for (j=0, fEl_loc = fEl; j<elems[i].numNodes(); j++, fEl_loc+=3) {
	F[ elems[i][j] ][0] -= fEl_loc[0];
	F[ elems[i][j] ][1] -= fEl_loc[1];
	F[ elems[i][j] ][2] -= fEl_loc[2];
      }
      break;
    }

    }

    // Add contribution to global matrix
    K.addContrib(elems[i].numNodes(), elems[i], kEl);
    if (P)
      P->addContrib(elems[i].numNodes(), elems[i], kEl);

  }

  if(ndType){
    for (i = 0; i < nodes.size(); ++i) {
      if (ndType[i] != BC_INTERNAL ) {

        F[i][0] = 0.0;
        F[i][1] = 0.0;
        F[i][2] = 0.0;
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, SVec<double,dim> &U)
{

  int ierr = 0;
  int numclipping = 0;
  double V[dim];
  double rho,p;

  if(varFcn->doVerification())
    for(int i=0; i<U.size(); i++)
      numclipping += varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1,U[i],V);
  else {
    for (int i=0; i<U.size(); ++i) {
      varFcn->conservativeToPrimitive(U[i], V);
      rho = varFcn->getDensity(V);
      p = varFcn->checkPressure(V);
      if (rho <= 0.0) {
        fprintf(stderr, "*** Error: negative density (%e) for node %d\n",
      	        rho, locToGlobNodeMap[i] + 1);
        ++ierr;
      }
      if (p <= 0.0) {
        fprintf(stderr, "*** Error: negative pressure (%e) for node %d (rho = %e)\n",
                p, locToGlobNodeMap[i] + 1, rho);
        ++ierr;
      }
    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, SVec<double,dim> &U, Vec<int> &fluidId)
{
  int ierr = 0;
  int numclipping = 0;
  double V[dim];
  double rho,p;

  if(varFcn->doVerification())
    for(int i=0; i<U.size(); i++)
      numclipping += varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1, U[i], V, fluidId[i]);
  else {
    for (int i=0; i<U.size(); ++i) {
      varFcn->conservativeToPrimitive(U[i], V, fluidId[i]);
      rho = varFcn->getDensity(V, fluidId[i]);
      p = varFcn->checkPressure(V, fluidId[i]);
      if (rho <= 0.0) {
        fprintf(stderr, "*** Error: negative density (%e) for node %d\n",
                rho, locToGlobNodeMap[i] + 1);
        ++ierr;
      }
      if (p <= 0.0) {
        fprintf(stderr, "*** Error: negative pressure (%e) for node %d (rho = %e)\n",
                p, locToGlobNodeMap[i] + 1, rho);
        ++ierr;
      }
    }
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, Vec<double> &ctrlVol, SVec<double,dim> &U,
                             Vec<int> &fluidId, Vec<int> &fluidIdn)
{
  int ierr = 0;
  int numclipping= 0;
  double V[dim];
  double rho, p;

  double *conservation = new double[5];
  for(int k=0; k<5; k++) conservation[k]=0.0;
  for(int i=0; i<U.size(); i++)
    for(int k=0; k<5; k++)
      conservation[k] += ctrlVol[i]*U[i][k];
  //fprintf(stdout, "conservation = %e %e %e %e %e\n", conservation[0],conservation[1],conservation[2],conservation[3],conservation[4]);

  if (!(varFcn->doVerification())){
    for (int i=0; i<U.size(); ++i) {

      if (!(U[i][0] > 0.0)) {
        fprintf(stderr, "*** Error: negative density (%e) for node %d (%d - %d)\n",
              U[i][0], locToGlobNodeMap[i] + 1, fluidId[i], fluidIdn[i]);
        ++ierr;
      }

      varFcn->conservativeToPrimitive(U[i], V, fluidId[i]);
      p = varFcn->checkPressure(V, fluidId[i]);
      if (p < 0.0) {
        fprintf(stderr, "*** Error: negative pressure (%e) for node %d (%d - %d)\n",
              p, locToGlobNodeMap[i] + 1,fluidId[i], fluidIdn[i]);
       ++ierr;
      }
    }
  }
  else{
    for (int i=0; i<U.size(); ++i) {

      if (!(U[i][0] > 0.0)) {
        fprintf(stderr, "*** Error: negative density (%e) for node %d with fluidID=%d (previously %d)\n",
              U[i][0], locToGlobNodeMap[i] + 1, fluidId[i], fluidIdn[i]);
        ++ierr;
      }

      numclipping += varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1, U[i], V, fluidId[i]);
    }
    //if (numclipping > 0) fprintf(stdout, "*** Warning: %d pressure clippings in subDomain %d\n", numclipping, globSubNum);
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::restrictionOnPhi(SVec<double,dim> &initial, Vec<int> &fluidId,
                                 SVec<double,dim> &restriction, int fluidIdTarget){

  int idim;
  restriction = 0.0;
  for (int i=0; i<nodes.size(); i++)
    if(fluidId[i] == fluidIdTarget)
      for(idim=0; idim<dim; idim++) restriction[i][idim] = initial[i][idim];

}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
int SubDomain::fixSolution(VarFcn *varFcn, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<int>* fluidId,int verboseFlag)
{

  int ierr = 0;

  for (int i=0; i<U.size(); ++i) {
    double V[dim];
    double Un[dim];

    for (int j=0; j<dim; ++j)
      Un[j] = U[i][j] + dU[i][j];

      int id = 0;
      if (fluidId)
	id = (*fluidId)[i];
      
      varFcn->conservativeToPrimitive(Un, V,id);
      double rho = varFcn->getDensity(V,id);
      double p = varFcn->checkPressure(V,id);
      
      if (rho <= 0.0) {
      if (verboseFlag == 4)
        fprintf(stderr, "*** Warning: negative density (%e) was fixed for node %d\n", rho, locToGlobNodeMap[i] + 1);

      for (int j=0; j<dim; ++j)
        dU[i][j] = 0.0;

      ++ierr;
    }
    if (p <= 0.0) {
      if (verboseFlag == 4)
        fprintf(stderr, "*** Warning: negative pressure (%e) was fixed for node %d (rho = %e)\n", p, locToGlobNodeMap[i] + 1, rho);

      for (int j=0; j<dim; ++j)
        dU[i][j] = 0.0;

      ++ierr;
    }

  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int neq>
int SubDomain::clipSolution(TsData::Clipping ctype, BcsWallData::Integration wtype,
			    VarFcn* varFcn, double* Uin, bool* flag, SVec<double,dim>& U,
			    int* cmin, int* pmin, double* vmin)
{

  int ierr = 0;

  double V[dim];
  varFcn->conservativeToPrimitive(U[0], V);

  int k;
  for (k=0; k<neq; ++k) {
    cmin[k] = 0;
    pmin[k] = locToGlobNodeMap[0] + 1;
    vmin[k] = V[dim-neq+k];
  }

  for (int i=0; i<U.size(); ++i) {
    varFcn->conservativeToPrimitive(U[i], V);
    double rho = varFcn->getDensity(V);
    double p = varFcn->getPressure(V);

    if (rho <= 0.0) {
      fprintf(stderr, "*** Error: negative density (%e) for node %d\n",
	      rho, locToGlobNodeMap[i] + 1);
      ++ierr;
    }
    if (p <= 0.0) {
      fprintf(stderr, "*** Error: negative pressure (%e) for node %d\n",
	      p, locToGlobNodeMap[i] + 1);
      ++ierr;
    }

    if ((wtype == BcsWallData::WALL_FUNCTION) ||
	(wtype == BcsWallData::FULL &&
	 nodeType[i] != BC_ISOTHERMAL_WALL_MOVING &&
	 nodeType[i] != BC_ISOTHERMAL_WALL_FIXED &&
	 nodeType[i] != BC_ADIABATIC_WALL_MOVING &&
	 nodeType[i] != BC_ADIABATIC_WALL_FIXED)) {
      for (k=0; k<neq; ++k) {
	if (V[dim-neq+k] < 0.0) {
	  if (flag[i]) {
	    cmin[k]++;
	    if (V[dim-neq+k] < vmin[k]) {
	      pmin[k] = locToGlobNodeMap[i] + 1;
	      vmin[k] = V[dim-neq+k];
	    }
	  }
	  if (ctype == TsData::ABS_VALUE)
	    U[i][dim-neq+k] = fabs(U[i][dim-neq+k]);
	  else if (ctype == TsData::FREESTREAM)
	    U[i][dim-neq+k] = Uin[dim-neq+k];
	}
      }
    }
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkFailSafe(VarFcn* varFcn, SVec<double,dim>& U,
                  SVec<bool,2>& tag, Vec<int> *fluidId)
{

  for (int i=0; i<U.size(); ++i) {
    tag[i][0] = false;
    tag[i][1] = false;
    double V[dim];
    double rho, p;
    if(!fluidId){
      varFcn->conservativeToPrimitive(U[i], V);
      rho = varFcn->getDensity(V);
      p = varFcn->getPressure(V);
    }else{
      varFcn->conservativeToPrimitive(U[i], V, (*fluidId)[i]);
      rho = varFcn->getDensity(V,(*fluidId)[i]);
      p = varFcn->getPressure(V, (*fluidId)[i]);
    }
    if (rho <= 0.0 || p <= 0.0)
      tag[i][0] = true;
  }

  int (*edgePtr)[2] = edges.getPtr();
  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];
    tag[i][1] = tag[i][1] || tag[j][0];
    tag[j][1] = tag[j][1] || tag[i][0];
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkGradientsSetUp(SVec<double,3> &X, SVec<double,dim> &V)
{

  double a = 3.0, b = -2.0, c = 1.0, d = -10.0;

  for (int i=0; i<X.size(); ++i)
    for (int k=0; k<dim; ++k)
      V[i][k] = a*X[i][0] + b*X[i][1] + c*X[i][2] + d;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkGradientsWrite(SVec<double,3> &X, NodalGrad<dim> &ngrad)
{

  double a = 3.0, b = -2.0, c = 1.0;

  char fileName[MAXLINE];
  sprintf(fileName, "gradients.%d", globSubNum+1);

  FILE *fp = fopen(fileName, "w");

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  int j = 2;
  for (int i=0; i<dVdx.size(); ++i)
    fprintf(fp, "%d %e %e %e\n", locToGlobNodeMap[i]+1,
	    dVdx[i][j]-a, dVdy[i][j]-b, dVdz[i][j]-c);

  fclose(fp);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkMatVecProd(SVec<double,dim> &prod, const char *msg)
{

  Vec<int> nodeCodes(prod.size());
  nodeCodes = 0;

  int i, j;
  for (i=0; i<faces.size(); ++i)
    for (j=0; j<faces[i].numNodes(); ++j)
      nodeCodes[ faces[i][j] ] += 1;

  char fname1[MAXLINE], fname2[MAXLINE];
  sprintf(fname1, "%s.in.%d", msg, globSubNum+1);
  sprintf(fname2, "%s.bound.%d", msg, globSubNum+1);

  FILE *fp1 = fopen(fname1, "w");
  FILE *fp2 = fopen(fname2, "w");

 /*
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    fprintf(fp1, "neighbor %d:\n", neighb[iSub]);
    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      int glLeft  = sharedEdges[iSub][iEdge].glLeft;
      int glRight = sharedEdges[iSub][iEdge].glRight;
      int sign = sharedEdges[iSub][iEdge].sign;
      int num = sharedEdges[iSub][iEdge].edgeNum;
      int flag = edgeFlag[num];
      int i = locToGlobNodeMap[ edgePtr[num][0] ];
      int j = locToGlobNodeMap[ edgePtr[num][1] ];
      fprintf(fp1, "   %d (%d) %d (%d) %d %d\n", glLeft, i, glRight, j, sign, flag);
    }
  }
 */

  for (i=0; i<prod.size(); ++i) {
    FILE* fp;
    if (nodeCodes[i] == 0)
      fp = fp1;
    else
      fp = fp2;

    fprintf(fp, "%d", locToGlobNodeMap[i]+1);
    for (j=0; j<dim; ++j)
      fprintf(fp, " %e", prod[i][j]);
    fprintf(fp, "\n");
  }

  fclose(fp1);
  fclose(fp2);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceDerivs(VarFcn *varFcn, SVec<double,3> &X,
                                   SVec<double,dim> &V, SVec<double,dim> &deltaU,
                                   Vec<double> &modalF, SVec<double, 3> **locMX)  {

  modalF = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, locMX);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, GeoState &geoState,
                                         BcData<dim> &bcData, SVec<double,3> &X, SVec<double,dim> &V,
					 double pInfty, Vec3D &CFi, Vec3D &CMi, Vec3D &CFv, Vec3D &CMv,
                                         VecSet< SVec<double,3> > *mX , Vec<double> *genCF)
{

  CFi = 0.0;
  CMi = 0.0;
  CFv = 0.0;
  CMv = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty,
                                      CFi, CMi, CFv, CMv, gradP, mX, genCF);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::zeroInternalVals(SVec<double, dim> &v)  {

  for (int i = 0; i < nodes.size(); ++i)
    if (nodeType[i] == 0)
      for (int j = 0; j < dim; j++)
        v[i][j] = 0.0;
}

//------------------------------------------------------------------------------

// HB
template<int dim>
void SubDomain::zeroMeshMotionBCDofs(SVec<double,dim> &x, int* DofType)
{
  int (*dofType)[dim] = reinterpret_cast<int (*)[dim]>(DofType);
  for(int i=0;i<nodes.size(); i++)
    for(int l=0; l<dim; l++)
      if(dofType[i][l]!=BC_FREE) x[i][l] = 0.0;
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setupUVolumesInitialConditions(const int volid, double UU[dim],
                                               SVec<double,dim> &U){

  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == volid)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++)
        for (int idim = 0; idim<dim; idim++)
          U[nodeNums[iNode]][idim] = UU[idim];
    }
  }

}

//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsSphere(FluidModelData &fm,
                       SphereData &ic, SVec<double,3> &X, SVec<double,dim> &U){

  double dist = 0.0;
  double x = ic.cen_x;
  double y = ic.cen_y;
  double z = ic.cen_z;
  double r = ic.radius;

  for (int i=0; i<U.size(); i++){
    dist = (X[i][0] - x)*(X[i][0] - x) + (X[i][1] - y)*(X[i][1] - y) + (X[i][2] - z)*(X[i][2] - z);
    if(sqrt(dist) < r) //it is inside the sphere
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
  }

}
*/
//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
                       PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U){


  double scalar = 0.0;
  double x = ip.cen_x;
  double y = ip.cen_y;
  double z = ip.cen_z;
  double nx = ip.nx;
  double ny = ip.ny;
  double nz = ip.nz;

  for (int i=0; i<U.size(); i++){
    scalar = nx*(X[i][0] - x)+ny*(X[i][1] - y)+nz*(X[i][2] - z);
    if(scalar > 0.0) //node is on the same side indicated by vector
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
  }

}
*/
//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
                       PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U, Vec<int> &nodeTag){

  double scalar = 0.0;
  double x = ip.cen_x;
  double y = ip.cen_y;
  double z = ip.cen_z;
  double nx = ip.nx;
  double ny = ip.ny;
  double nz = ip.nz;

  for (int i=0; i<U.size(); i++){
    scalar = nx*(X[i][0] - x)+ny*(X[i][1] - y)+nz*(X[i][2] - z);
    if(scalar > 0.0) {//node is on the same side indicated by vector
      nodeTag[i] = -1;
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
    }
  }

}
*/
//------------------------------------------------------------------------------
// TODO: should distinguish master nodes and non-master nodes
template<int dim>
void SubDomain::computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights,
                      Vec<double> &Weights, LevelSetStructure &LSS, SVec<double,3> &X, Vec<double> &init, Vec<double> &next_init)
{
  bool* masterEdge = edges.getMasterFlag();
  const Connectivity &nToN = *getNodeToNode();
  for(int currentNode=0;currentNode<numNodes();++currentNode)
    if(init[currentNode]<1.0 && LSS.isActive(0.0,currentNode)){
      for(int j=0;j<nToN.num(currentNode);++j){
        int neighborNode=nToN[currentNode][j];
        if(currentNode == neighborNode || init[neighborNode]<1.0) continue;
        int l = edges.findOnly(currentNode,neighborNode);
        if(LSS.edgeIntersectsStructure(0.0,l) || !masterEdge[l]) continue;
        else if(Weights[currentNode] < 1e-6){
          Weights[currentNode]=1.0;
          next_init[currentNode]=1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] = V[neighborNode][i];
	//  fprintf(stderr,"Node %d Giving Node %d [%e %e %e %e %e]\n",locToGlobNodeMap[neighborNode]+1,locToGlobNodeMap[currentNode]+1,
	//		  V[neighborNode][0],V[neighborNode][1],V[neighborNode][2],V[neighborNode][3],V[neighborNode][4]);
        } else {
          Weights[currentNode] += 1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] += V[neighborNode][i];
        //  fprintf(stderr,"Node %d Giving Node %d [%e %e %e %e %e]\n",locToGlobNodeMap[neighborNode]+1,locToGlobNodeMap[currentNode]+1,
	//		  V[neighborNode][0],V[neighborNode][1],V[neighborNode][2],V[neighborNode][3],V[neighborNode][4]);
        }
      }
    }
}

//------------------------------------------------------------------------------
// who: active and "swept" nodes
// what: phi and U
// how: pull from neighbours which are visible, not swept, and have the same fluid Id.
// TODO: more efficent in an "edge loop".
template<int dim, int dimLS>
void SubDomain::computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights,
                                                SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights, 
                                                Vec<double> &Weights, LevelSetStructure &LSS, SVec<double,3> &X,
                                                Vec<double> &init, Vec<double> &next_init, Vec<int> &fluidId)
{
  bool* masterEdge = edges.getMasterFlag();
  const Connectivity &nToN = *getNodeToNode();
  for(int currentNode=0;currentNode<numNodes();++currentNode)
    if(init[currentNode]<1.0 && LSS.isActive(0.0,currentNode)){
      int myId = fluidId[currentNode]; 
      for(int j=0;j<nToN.num(currentNode);++j){
        int neighborNode=nToN[currentNode][j];
        int yourId = fluidId[neighborNode];
        if(currentNode==neighborNode || init[neighborNode]<1.0 || myId!=yourId) continue;
        int l = edges.findOnly(currentNode,neighborNode);
        if(LSS.edgeIntersectsStructure(0.0,l) || !masterEdge[l]) continue;
        else if(Weights[currentNode] < 1e-6){
          Weights[currentNode]=1.0;
          next_init[currentNode]=1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] = V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] = Phi[neighborNode][i];
        //   fprintf(stderr,"Node %d (%d) Giving Node %d (%d) [%e %e %e %e %e]\n",
	// 		  locToGlobNodeMap[neighborNode]+1,fluidId[neighborNode],locToGlobNodeMap[currentNode]+1,fluidId[currentNode],
	// 		  V[neighborNode][0],V[neighborNode][1],V[neighborNode][2],V[neighborNode][3],V[neighborNode][4]);
        } else {
          Weights[currentNode] += 1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] += V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] += Phi[neighborNode][i];
          // fprintf(stderr,"Node %d (%d) Giving Node %d (%d) [%e %e %e %e %e]\n",
	  //       	  locToGlobNodeMap[neighborNode]+1,fluidId[neighborNode],locToGlobNodeMap[currentNode]+1,fluidId[currentNode],
	  //       	  V[neighborNode][0],V[neighborNode][1],V[neighborNode][2],V[neighborNode][3],V[neighborNode][4]);
        }
      }
}

//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::extrapolatePhiV(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV)
{
  int (*edgePtr)[2] = edges.getPtr();
  bool *masterFlag = edges.getMasterFlag();
  
  for(int l=0; l<edges.size(); l++) {
    if(!masterFlag[l])
      continue;
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    if(!LSS.withCracking()) {
      int Idi = LSS.fluidModel(0.0,i), Idj = LSS.fluidModel(0.0,j);
      if(Idi!=0 || Idj!=0) //meaning one of them is isolated by structure
        continue;
    }

    bool iSwept = LSS.isSwept(0.0,i), jSwept = LSS.isSwept(0.0,j);
 
    // pull data from j to i?
    if(iSwept && !jSwept)
      for(int k=0; k<dimLS; k++)
        PhiV[i][k] += PhiV[j][k];

    // pull data from i to j?
    if(jSwept && !iSwept)
      for(int k=0; k<dimLS; k++)
        PhiV[j][k] += PhiV[i][k];
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::populateGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints,SVec<double,dim> &U,VarFcn *varFcn,LevelSetStructure &LSS,Vec<int> &tag)
{

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++) {
    i = edgePtr[l][0];
    j = edgePtr[l][1];
    if(LSS.edgeIntersectsStructure(0.0,l)) { // at interface
      int tagI = tag[i];
      int tagJ = tag[j];
      bool iIsActive = LSS.isActive(0.0,i);
      bool jIsActive = LSS.isActive(0.0,j);

      if(iIsActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
        Vec<double> Vi(dim);
        varFcn->conservativeToPrimitive(U[i],Vi.v,tagI);
        if(!ghostPoints[j]) // GP has not been created
        {ghostPoints[j]=new GhostPoint<dim>;}
        // If the edge is not a master edge, do nothing. Some other CPU is gonna do the job
        if(edgeFlag[l]) {ghostPoints[j]->addNeighbour(Vi,1.0,resij.normVel,tagI);}
      }
      if(jIsActive) {
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
        Vec<double> Vj(dim);
        varFcn->conservativeToPrimitive(U[j],Vj.v,tagJ);
        if(!ghostPoints[i]) // GP has not been created
        {ghostPoints[i]=new GhostPoint<dim>;}
        // If the edge is not a master edge, do nothing. Some other CPU is gonna do the job
        if(edgeFlag[l]) {ghostPoints[i]->addNeighbour(Vj,1.0,resji.normVel,tagJ);}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::reduceGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints)
{
  for (int i=0; i<nodes.size(); i++)
    {
      if(ghostPoints[i]) // GP has already been created
	{
	  ghostPoints[i]->reduce();
	}
    } 
}

template<int dim,int neq>
void SubDomain::populateGhostJacobian(Vec<GhostPoint<dim>*> &ghostPoints,SVec<double,dim> &U,VarFcn *varFcn,LevelSetStructure &LSS,Vec<int> &tag,
                                      GenMat<double,neq>& A)
{

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  double B[neq*neq],tmp[neq*neq];
  double dUdV[neq*neq],dVdU[neq*neq];
  memset(tmp,0,sizeof(double)*neq*neq);
  memset(dUdV,0,sizeof(double)*neq*neq);
  memset(dVdU,0,sizeof(double)*neq*neq);
  memset(B,0,sizeof(double)*neq*neq);
  B[0] = B[5*neq-1] = 1.0;
  for (k = 1; k < 4; ++k) {
    B[k*neq+k] = -1.0;
  }
  for (int l=0; l<edges.size(); l++) {
    i = edgePtr[l][0];
    j = edgePtr[l][1];
    if(LSS.edgeIntersectsStructure(0.0,l)) { //at interface
	  int tagI = tag[i];
	  int tagJ = tag[j];
      bool iIsActive = LSS.isActive(0.0,i);
      bool jIsActive = LSS.isActive(0.0,j);
	  
	  if(iIsActive) {
        double (*Aij)[neq*neq] = reinterpret_cast< double (*)[neq*neq]>(A.getGhostNodeElem_ij(i,j));
        Vec<double> Vi(dim);
        varFcn->conservativeToPrimitive(U[i],Vi.v,tagI);
        varFcn->getVarFcnBase()->computedVdU(Vi.v,dVdU);
        varFcn->getVarFcnBase()->computedUdV(ghostPoints[j]->getPrimitiveState(), dUdV);
        DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&dUdV, 0, &B, 0, &tmp, 0);
        DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &dVdU, 0, Aij, 0);
      }
      if(jIsActive) {
        double (*Aji)[neq*neq] = reinterpret_cast< double (*)[neq*neq]>(A.getGhostNodeElem_ij(j,i));
        Vec<double> Vj(dim);
        varFcn->conservativeToPrimitive(U[j],Vj.v,tagJ);
        varFcn->getVarFcnBase()->computedVdU(Vj.v,dVdU);
        varFcn->getVarFcnBase()->computedUdV(ghostPoints[i]->getPrimitiveState(), dUdV);
        DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&dUdV, 0, &B, 0, &tmp, 0);
        DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &dVdU, 0, Aji, 0);
	  }
	}
  }

  for (int l = 0; l < ghostPoints.size(); ++l) {
    if (ghostPoints[l]) {
      double *Aii = A.getGhostGhostElem_ij(i,i);
      memset(Aii,0,sizeof(double)*neq*neq);
      for (k = 0; k < neq; ++k)
        Aii[k*neq+k] = ghostPoints[l]->lastCount();
    }
  }
}

//--------------------------------------------------------------------------
//TODO: should distinguish master edges and non-master edges.
template<int dim>
void SubDomain::computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights, 
                      LevelSetStructure &LSS, SVec<double,3> &X)
{

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    if(LSS.edgeIntersectsStructure(0.0,l)){ //at interface
      if (LSS.isActive(0.0,i)) {// add Wstarij on node j.
        if(Weights[j]<1.e-6) {
          Weights[j] = 1.0;
          if(Wstarij[l][0]>1.0e-8) //use Wstarij.
            for(k=0; k<dim; k++)
              VWeights[j][k] = Wstarij[l][k];
          else { 
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[j][k] = V[i][k];
          }
        }else{
          Weights[j] += 1.0;
          if(Wstarij[l][0]>1.0e-8)//use Wstarij
            for(k=0; k<dim; k++)
              VWeights[j][k] += Wstarij[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[j][k] += V[i][k];
          }
        }
      }
 
      if (LSS.isActive(0.0,j)) {// add Wstarji on node i
        if(Weights[i]<1.e-6) {
          Weights[i] = 1.0;
          if(Wstarji[l][0]>1.0e-8) //use Wstarji.
            for(k=0; k<dim; k++)
              VWeights[i][k] = Wstarji[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[i][k] = V[j][k];
          }
        }else{
          Weights[i] += 1.0;
          if(Wstarji[l][0]>1.0e-8)//use Wstarji.
            for(k=0; k<dim; k++)
              VWeights[i][k] += Wstarji[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[i][k] += V[j][k];
          }
        }
      }
 
    }

  }
}

//--------------------------------------------------------------------------
//TODO: should distinguish master edges and non-master edges.
template<int dim, int dimLS>
void SubDomain::computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights,
                      SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights,  
                      LevelSetStructure &LSS, SVec<double,3> &X, Vec<int> &fluidId0, Vec<int> &fluidId)
{
  int i, j, k;
  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    for (int enode = 0; enode<2; enode++) { //loop thru the two nodes on this edge
      SVec<double,dim> &Wstar = (enode==0) ? Wstarji : Wstarij;
      i = edgePtr[l][enode];
      j = edgePtr[l][(int)(!enode)];
      // for V and Phi
      if(LSS.isSwept(0.0,i) && !LSS.isOccluded(0.0,i)) { // phase change occurred && need an update
        if(Wstar[l][0]>1.0e-8 && fluidId0[j]==fluidId[i]) { //use Wstar 
          if(Weights[i]<1.0e-6) { // first touch of node i
            Weights[i] = 1.0;
            for(k=0; k<dim; k++) VWeights[i][k] = Wstar[l][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] = Phi[j][k];
          } else {
            Weights[i] ++;
            for(k=0; k<dim; k++) VWeights[i][k] += Wstar[l][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] += Phi[j][k];
          }
        } else if(!LSS.isSwept(0.0,j) && fluidId[i]==fluidId[j] && !LSS.edgeIntersectsStructure(0.0,l)) { // use V[j]
          if(Weights[i]<1.0e-6) { // first touch of node i
            Weights[i] = 1.0;
            for(k=0; k<dim; k++) VWeights[i][k] = V[j][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] = Phi[j][k];
          } else {
            Weights[i] ++;
            for(k=0; k<dim; k++) VWeights[i][k] += V[j][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] += Phi[j][k];
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
template<int dimLS>
void SubDomain::checkNodePhaseChange(SVec<double,dimLS> &PhiProduct)
{

  for(int i=0; i<nodes.size(); i++)
    for(int idim=0; idim<dimLS; idim++)
      if(PhiProduct[i][idim]<0.0)
        fprintf(stdout, "***Error: node %d (%d) has changed phase during reinitialization\n", i, locToGlobNodeMap[i]+1);

}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::storePreviousPrimitive(SVec<double,dim> &V, Vec<int> &fluidId, 
                                    SVec<double,3> &X, SVec<double,dim> &Vupdate, 
                                    Vec<double> &weight){

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    i = edgePtr[l][0];
    j = edgePtr[l][1];

// selective averaged extrapolation (based on direction of the flow)

    if(fluidId[i] != fluidId[j]){ //at interface
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      double normdx2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      double normUi2 = V[i][1]*V[i][1]+V[i][2]*V[i][2]+V[i][3]*V[i][3];
      double normUj2 = V[j][1]*V[j][1]+V[j][2]*V[j][2]+V[j][3]*V[j][3];
      double udotdx = 0.0;

      if(normdx2*normUj2 > 0.0)
        udotdx = -(dx[0]*V[j][1]+dx[1]*V[j][2]+dx[2]*V[j][3])/sqrt(normdx2*normUj2);
      if(udotdx > 0.0){
        weight[i] += udotdx;
        for(k=0; k<dim; k++)
          Vupdate[i][k] += udotdx*V[j][k];
      }

      if(normdx2*normUi2 > 0.0)
        udotdx = (dx[0]*V[i][1]+dx[1]*V[i][2]+dx[2]*V[i][3])/sqrt(normdx2*normUi2);
      if(udotdx > 0.0){
        weight[j] += udotdx;
        for(k=0; k<dim; k++)
          Vupdate[j][k] += udotdx*V[i][k];
      }
    }

  }

}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U){
// only the pressure in the volumes that have an id=0 (outside part of a cylinder phi>0)
// are updated. It is assumed that the input and output states are uniform!

  bool found = false;
  double Uc[dim];

  for(int i=0; i<nodes.size(); i++) {
    if(!found) {//fill Uc[dim]
      double V[dim];
      vf->conservativeToPrimitive(U[i],V);
      vf->setPressure(p,V); //for Tait it actually modifies density.
      vf->primitiveToConservative(V,Uc);
      found = true;
    }
    for(int k=0; k<dim; k++)
      U[i][k] = Uc[k];
  }

/*  double V[dim];
  vf->conservativeToPrimitive(U[0],V);
  vf->setPressure(p,V);
  double rhoe = vf->computeRhoEnergy(V);


  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == 0)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++)
        U[nodeNums[iNode]][4] = rhoe;
    }
  }
*/
}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U, Vec<int> &fluidId){
// only the pressure with fluidId = 0 (outside part)
// are updated. It is assumed that the input and output states are uniform!

  bool found = false;
  double Uc[dim];

  for(int i=0; i<nodes.size(); i++) {
    if(fluidId[i]!=0)
      continue;
    if(!found) {//fill Uc[dim]
      double V[dim];
      vf->conservativeToPrimitive(U[i],V,fluidId[i]);
      vf->setPressure(p,V,fluidId[i]); //for Tait it actually modifies density.
      vf->primitiveToConservative(V,Uc,fluidId[i]);
      found = true;
    }
    for(int k=0; k<dim; k++)
      U[i][k] = Uc[k];
  }
/*
  bool found = false;
  double rhoe = -1.0; 

// only the pressure with fluidId = 0 (outside part)
// are updated. It is assumed that the input and output states are uniform!
  for(int i=0; i<nodes.size(); i++) {
    if(fluidId[i]!=0) 
      continue;
    if(!found){ //rhoe is not computed yet.
      double V[dim];
      vf->conservativeToPrimitive(U[i],V,fluidId[i]);
      vf->setPressure(p,V,fluidId[i]);
      rhoe = vf->computeRhoEnergy(V,fluidId[i]);
      found = true;
    }    
    U[i][4] = rhoe;
  }
*/
}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::printVariable(SVec<double,dim> &V, VarFcn *vf)
{
  inletNodes.printVariable(V, sharedInletNodes, vf);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::printInletVariable(SVec<double,dim> &V)
{
  inletNodes.printInletVariable(V, sharedInletNodes);
}

//-----------------------------------------------------------------------------

template<int dim>
void SubDomain::printAllVariable(Vec<int> &X, SVec<double,dim> &U, int numSub, int it)
{

  int glob;
  for (int i=0; i<nodes.size(); i++){
    glob = locToGlobNodeMap[i]+1;
    fprintf(stdout, "Tag[%d,%d] = %d\n", i, glob, X[i]);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkExtrapolationValue(SVec<double,dim> &V, VarFcn *vf,
                                        BcData<dim>& bcData, GeoState& geoState)
{
  inletNodes.checkExtrapolationValue(V, sharedInletNodes, nodeType, vf, bcData, geoState);
}

//------------------------------------------------------------------------------

template<class Scalar, int neq>
void SubDomain::printAllMatrix(GenMat<Scalar,neq> &A, int it)
{
    fprintf(stdout, "subDomain %d:\n", locSubNum);
    int glob;
    for (int i=0; i<nodes.size(); i++){
      glob = locToGlobNodeMap[i]+1;
      Scalar *Aii = A.getElem_ii(i);
      fprintf(stdout, "%d %d ", i, glob);
      for (int j = 0; j<neq; j++){
        if(j==0)
          fprintf(stdout, "Aii[%d] = ", i);
        else
          fprintf(stdout, "          ");
        for(int k=0; k<neq; k++){
          fprintf(stdout, "%.6e ", Aii[j*neq+k]);
          fprintf(stdout, "\n");
        }
      }
      fprintf(stdout, "\n");
    }

}
//--------------------------------------------------------------------------
template<int dim>
void SubDomain::padeReconstruction(SVec<double, dim> **dataCoarse, SVec<double, dim> **data, int *stepParam, double *freqCoarse, double deltaFreqFine, int nStrMode,int L, int M, int nPoints)
{


  int i, j;
  int numFreqCoarse = stepParam[1];
  int nSnapsCoarse = stepParam[3];
  int nSteps = stepParam[0];
  int numPadeDeriv = int(ceil((L+M+1.0)/((double)nPoints))-1.0);
  int size = L+M+1;


  bcomp *compMat = new bcomp[numFreqCoarse*(numPadeDeriv+1)];

  double tempMat[dim*nSnapsCoarse];
  bcomp *padeMat = new bcomp[size*size];
  bcomp *padeVec = new bcomp[size];

  bcomp *snaps = new bcomp[nSteps+1];
  int midFreq[1] ;
  double deltaFreqCoarse[numFreqCoarse];
  buildDeltaFreq(deltaFreqCoarse, numFreqCoarse, freqCoarse, midFreq);
  for (int iNode = 0; iNode < (*(data[0])).len; iNode++) {


    extractElementsRelativeToANode(dataCoarse, tempMat, iNode, nSnapsCoarse);
    for (int iDim = 0; iDim < dim; iDim++) {
      for (int iStrMode = 0; iStrMode < nStrMode; iStrMode++) {
        extractElementsRelativeToAComponentAndAMode(tempMat,compMat,iDim,iStrMode,numPadeDeriv,numFreqCoarse,nStrMode,nSnapsCoarse, freqCoarse[0]);
        multiPade(compMat, stepParam, deltaFreqCoarse, padeMat, padeVec, L, M, nPoints, deltaFreqFine, freqCoarse[midFreq[0]], snaps, freqCoarse);
        snapshotsConstruction(data,snaps,nSteps,iDim,iStrMode,iNode,nStrMode,freqCoarse[0]);

      }
    }
  }
}

//--------------------------------------------------------------------------
template<int dim>
void SubDomain::extractElementsRelativeToANode(SVec<double, dim> **dataCoarse, double *tempMat, int iNode, int nSnapsCoarse)
{

  for (int jSnapsCoarse = 0; jSnapsCoarse < nSnapsCoarse; jSnapsCoarse++) {

    SVec<double, dim> *X = dataCoarse[jSnapsCoarse];
    for (int iDim = 0; iDim < dim; iDim++) {
      *(tempMat + iDim*nSnapsCoarse + jSnapsCoarse)  = (*X)[iNode][iDim];  //tempMat[iDim][jSnapsCoarse]

    }

  }
}

//--------------------------------------------------------------------------
template<int dim>
void SubDomain::snapshotsConstruction(SVec<double, dim> **data, bcomp* snaps, int nSteps, int iDim, int iStrMode, int iNode, int nStrMode, double freq1)
{

  SVec<double,dim> *Y = 0;
  SVec<double,dim> *Z = 0;
  if (freq1 == 0.0) {

    for (int jFineFreq=0; jFineFreq<nSteps+1; jFineFreq++) {
      if (jFineFreq==0) {
        Y = data[iStrMode];
        (*Y)[iNode][iDim] = real(snaps[0]);
      }
      else {
        Y = data[nStrMode+2*nStrMode*(jFineFreq-1)+2*iStrMode];
        (*Y)[iNode][iDim] = real(snaps[jFineFreq]);
        Z = data[nStrMode+2*nStrMode*(jFineFreq-1)+2*iStrMode+1];
        (*Z)[iNode][iDim] = imag(snaps[jFineFreq]);
      }
    }
  }
  else {
    for (int jFineFreq=0; jFineFreq<nSteps+1; jFineFreq++) {
      Y = data[2*nStrMode*jFineFreq+2*iStrMode];
      (*Y)[iNode][iDim] = real(snaps[jFineFreq]);
      Z = data[2*nStrMode*jFineFreq+2*iStrMode+1];
      (*Z)[iNode][iDim] = imag(snaps[jFineFreq]);
    }
  }
  Y = 0;
  Z = 0;

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::hardyInterpolationLogMap(SVec<double, dim> ***dataCoarse, SVec<double, dim> **dataInterp, int nData, int numPod, int iDataMin, FullM &B, FullM &b)
{

  double tempMat[dim*nData];
  double *hardyCoefs = new double[nData];
  //for each vector
  for (int iPod = 0; iPod < numPod; ++iPod) {
    //for each node
    for (int iNode = 0; iNode < (*(dataInterp[0])).len; ++iNode) {
      extractElementsRelativeToANodeAndAVector(dataCoarse,tempMat,iNode,nData,iDataMin,iPod);
      // for each component
      for (int iDim = 0; iDim < dim; ++iDim) {
        //hardyCoefs = B*data(component)
        for (int iData = 0; iData < nData; ++iData) {
          hardyCoefs[iData] = 0.0;
          for (int jData = 0; jData < nData; ++jData)
            hardyCoefs[iData] += B[iData][jData]*tempMat[iDim*nData+jData];
        }
        //compute the interpolated analog quantity
        (*(dataInterp[iPod]))[iNode][iDim] = 0.0;
        for (int iData = 0; iData < nData; ++iData)
          (*(dataInterp[iPod]))[iNode][iDim] += b[iData][0]*hardyCoefs[iData];
      }
    }
  }
  delete [] hardyCoefs;

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::extractElementsRelativeToANodeAndAVector(SVec<double, dim> ***dataCoarse, double *tempMat, int iNode, int nData, int jDataMin, int iPod)
{


  for (int jData = 0; jData < nData; ++jData) {
    if (jData !=jDataMin) {
      SVec<double, dim> *X = dataCoarse[jData][iPod];
      for (int iDim = 0; iDim < dim; ++iDim)
        *(tempMat + iDim*nData + jData)  = (*X)[iNode][iDim];  //tempMat[iDim][jData]
    }
    else {
      for (int iDim = 0; iDim < dim; ++iDim)
        *(tempMat + iDim*nData + jData)  = 0.0;
    }

  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::getGradP(NodalGrad<dim>& ngrad)
{

  // gradents stored for future use in computeForce and computeForceTransmitted
  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  for (int i=0;i<nodes.size();i++) {
    gradP[0][i] = dVdx[i][4];
    gradP[1][i] = dVdy[i][4];
    gradP[2][i] = dVdz[i][4];
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::getDerivativeOfGradP(NodalGrad<dim>& ngrad)
{

  SVec<double,dim> &ddVdx = ngrad.getXderivative();
  SVec<double,dim> &ddVdy = ngrad.getYderivative();
  SVec<double,dim> &ddVdz = ngrad.getZderivative();

  for (int i=0;i<nodes.size();i++) {
    dGradP[0][i] = ddVdx[i][4];
    dGradP[1][i] = ddVdy[i][4];
    dGradP[2][i] = ddVdz[i][4];
  }

}

//-----------------------------------------------------------------------------

template<int dim>
void SubDomain::computePrdtWCtrlVolRatio(SVec<double,dim> &ratioTimesU, SVec<double,dim> &U, Vec<double> &ctrlVol, GeoState &geoState)
{
   Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();

   for (int i=0; i<nodes.size(); ++i) {
     double ratio = ctrlVol_n[i]/ctrlVol[i];
     for (int j=0; j<dim; ++j) {
       ratioTimesU[i][j] = ratio * U[i][j];
     }
   }

}

//------------------------------------------------------------------------------
//             LEVEL SET SOLUTION AND REINITIALIZATION                       ---
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin)
{
  if(!NodeToNode)
     NodeToNode = createEdgeBasedConnectivity();
  fprintf(stderr, "***Error: this routine is not valid for multiple subdomains\n");
  exit(1);

  for(int i=0; i<nodes.size(); i++){
    for(int j=0; j<dimLS; j++){
      if(Phi[i][j]*Phin[i][j]<0.0){
        // check if node i HAD a neighbour with a different levelset sign
        bool diffNeigh = false;
        for(int iNeigh=0; iNeigh<NodeToNode->num(i); iNeigh++)
          if(Phin[i][j]*Phin[(*NodeToNode)[i][iNeigh]][j]<0.0){
            diffNeigh = true;
            break;
          }
        if(!diffNeigh) {Phi[i][j] = Phin[i][j]; fprintf(stdout, "node %d has levelset %d clipped to avoid phase creation\n", locToGlobNodeMap[i]+1, j);}
      }
    }
  }

}
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin, Vec<double> &weight, LevelSetStructure *LSS)
{

  for(int i=0; i<nodes.size(); i++){
    int fModel;//if fModel>0 (isolated), Phi is not used at all
    if(LSS && !LSS->withCracking()) 
      fModel = LSS->fluidModel(0.0, i); 
    else
      fModel = 0;
    bool swept = LSS ? LSS->isSwept(0.0, i) : 0; // if swept, Phin is not reliable!
    for(int j=0; j<dimLS; j++){
      if(Phi[i][j]*Phin[i][j]<0.0 && fModel==0 && !swept){
        // check if node i HAD a neighbour with a different levelset sign
        if(weight[i] <= 0.0){
          //fprintf(stdout, "node %d (loc %d in %d) has weight = %f and has levelset %d"
          //                " moving from %e to %e\n", locToGlobNodeMap[i]+1,i,
          //               globSubNum,weight[i],j,Phin[i][j],Phi[i][j]);
          Phi[i][j] = Phin[i][j];
        }
      }
    }
  }

}
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi, int level)
{
  if(!NodeToNode)
     NodeToNode = createEdgeBasedConnectivity();

  if(level==0){
  // tag nodes that are closest to interface by looking at phi[i]*phi[j]
    Tag = 0;
    edges.TagInterfaceNodes(lsdim,Tag,Phi);

  }else{
  // tag nodes that are neighbours of already tagged nodes.
    int nNeighs,nei,k;
    for(int i=0; i<nodes.size(); i++){

      if(Tag[i]==level){

        nNeighs = NodeToNode->num(i);
        for(k=0;k<nNeighs;k++){
          nei = (*NodeToNode)[i][k];
          if(nei==i) continue;
          if(Tag[nei]==0) Tag[nei] = level+1;
        }

      }
    }
  }


}
//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::TagInterfaceNodes(int lsdim, SVec<bool,2> &Tag, SVec<double,dimLS> &Phi, LevelSetStructure *LSS)
{
  Tag = false;

  int i,j;
  int (*ptr)[2] = edges.getPtr(); 
  for(int l=0; l<edges.size(); l++) {
    i = ptr[l][0];
    j = ptr[l][1];

    if(Phi[i][lsdim]*Phi[j][lsdim]<=0.0) {
      if(LSS->edgeIntersectsStructure(0.0,l))
        Tag[i][0] = Tag[j][0] = true;
      else
        Tag[i][1] = Tag[j][1] = true;
    } else {
      if(LSS->edgeIntersectsStructure(0.0,l)) {
        Tag[i][0] = Tag[j][0] = true;
//        fprintf(stderr,"BUG: Sub %d: (%d,%d) intersects but phi[i]=%e, phi[j]=%e.\n", globSubNum, 
//                locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, Phi[i][lsdim], Phi[j][lsdim]);
//        fprintf(stderr,"  %d: occluded(%d), swept(%d).\n", locToGlobNodeMap[i]+1, LSS->isOccluded(0,i),LSS->isSwept(0,i));
//        fprintf(stderr,"  %d: occluded(%d), swept(%d).\n", locToGlobNodeMap[j]+1, LSS->isOccluded(0,j),LSS->isSwept(0,j));
      }
    }
  }
}

//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::printPhi(SVec<double, 3> &X, SVec<double,dimLS> &Phi, int it)
{
  fprintf(stdout, "\nPhi - subDomain %d: \n", locSubNum);
  int glob, sh;
  /*for (int iSub = 0; iSub < numNeighb; ++iSub) {
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      sh = (*sharedNodes)[iSub][iNode] ;
      glob = locToGlobNodeMap[sh]+1;
      fprintf(stderr, "%d %.14e %.14e %.14e %.14e %.14e\n", glob, V[ sh ][0],
                                                                  V[ sh ][1],
                                                                  V[ sh ][2],
                                                                  V[ sh ][3],
                                                                  V[ sh ][4]);
    }
  }*/

  for (int i=0; i<nodes.size(); i++){
    glob = locToGlobNodeMap[i]+1;
    fprintf(stdout, " Phi : %d %d   ", i,glob);
    for(int j=0; j<dimLS; j++)
      fprintf(stdout, "%e ", Phi[i][j]);
    fprintf(stdout, "\n");
  }

}
//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::setupPhiVolumesInitialConditions(const int volid, 
                    const int fluidId, SVec<double,dimLS> &Phi){
  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == volid)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++){
        for (int iDim = 0; iDim < dimLS; iDim++){
          if(iDim+1 == fluidId)
            Phi[nodeNums[iNode]][iDim] = 1.0;
        }
      }
    }
  }

}

//--------------------------------------------------------------------------
// for mesh motion (with RK2 time-integration)
template<int dimLS>
void SubDomain::computePrdtPhiCtrlVolRatio(SVec<double,dimLS> &ratioTimesPhi,
           SVec<double,dimLS> &Phi, Vec<double> &ctrlVol, GeoState &geoState)
{
   Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();

   for (int i=0; i<nodes.size(); ++i) {
     double ratio = ctrlVol_n[i]/ctrlVol[i];
     for (int j=0; j<dimLS; j++)
       ratioTimesPhi[i][j] = ratio * Phi[i][j];
   }

}

//-----------------------------------------------------------------------------
/*
template<int dimLS>
void SubDomain::FinishReinitialization(Vec<int> &Tag, SVec<double,dimLS> &Psi, int level)
{

  if(!NodeToNode)
    NodeToNode = createEdgeBasedConnectivity();

  int nNeighs,nei,k;
  for (int i=0; i<nodes.size(); i++){
    if (Tag[i]==level){

      nNeighs = NodeToNode->num(i);
      for (k=0; k<nNeighs; k++){
        nei = (*NodeToNode)[i][k];
        if(Tag[nei]==0){
          Tag[nei] = level+1;
          Psi[nei][0] = Psi[i][0];
        }else if(Tag[nei]==level+1){
          if( (Psi[i][0] > 0.0 && Psi[i][0] > Psi[nei][0]) ||
              (Psi[i][0] < 0.0 && Psi[i][0] < Psi[nei][0])  )
            Psi[nei][0] = Psi[i][0];
        }

      }
    }
  }

}
*/
//-----------------------------------------------------------------------------
template<int dimLS>
void SubDomain::copyCloseNodes(int lsdim, int level, Vec<int> &Tag,SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{
  for(int i=0; i<nodes.size(); i++) {
    if(Tag[i]==level)
      Psi[i][0] = fabs(Phi[i][lsdim]);
    //DEBUG
//    if(locToGlobNodeMap[i]+1==115697)
//      fprintf(stderr,"Sub %d, level = %d, Tag[%d] = %d. Phi = %e, Psi = %e.\n", globSubNum, level, locToGlobNodeMap[i]+1,
//              Tag[i], Phi[i][0], Psi[i][0]);
  }

}
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       NodalGrad<dimLS> &grad,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{
  for(int i=0; i<nodes.size(); i++){
    if(Tag[i]==1) Tag[i]=-1;
  }
  SVec<double,dimLS>& ddx  = grad.getX();
  SVec<double,dimLS>& ddy  = grad.getY();
  SVec<double,dimLS>& ddz  = grad.getZ();
  elems.computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
}
//-------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                      NodalGrad<dimLS> &grad, SVec<double,dimLS> &Phi,
                                      SVec<double,1> &Psi)
{

  SVec<double,dimLS>& ddx  = grad.getX();
  SVec<double,dimLS>& ddy  = grad.getY();
  SVec<double,dimLS>& ddz  = grad.getZ();
  elems.recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//-------------------------------------------------------------------------------
template<int dimLS>
double SubDomain::computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                       SVec<double,3> &X, SVec<double,1> &Psi,SVec<double,dimLS> &Phi)
{

  if(level==2)
    for(int i=0; i<nodes.size(); i++)
      if(Tag[i]==-1) Tag[i]=1;
  elems.computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);
  double res = 0.0;
  for(int i=0; i<nodes.size(); i++)
    if(Tag[i]==level)
      res += Psi[i][0]*Psi[i][0];

  return res;
}
//-------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::getSignedDistance(int lsdim, SVec<double,1> &Psi, SVec<double,dimLS> &Phi)
{
  for(int i=0; i<nodes.size(); i++){
    if(Phi[i][lsdim]<0.0)
      Psi[i][0] = -Psi[i][0];
    if(Phi[i][lsdim]<0.0 && Psi[i][0]>0.0)
      fprintf(stdout, "globnode %d (%d) has changed phase %e %e\n", locToGlobNodeMap[i]+1,i,Phi[i][lsdim],Psi[i][0]);
  }
}


//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeRecSurfBasedForceLoad(int forceApp, int order, SVec<double,3> &X,
                                             double (*Fs)[3], int sizeFs, LevelSetStructure &LSS, double pInfty, 
                                             SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji, SVec<double,dim> &V, 
                                             Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, VarFcn* vf, Vec<int>* fid)
{
  if (forceApp!=3) 
    {fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); exit(-1);}

  // ---------------
  //   Preparation
  // ---------------
  int T[4];       //nodes in a tet.

  double d2w[3]; // not used, but required by postFcn->computeForces(...)..
  int fid_local[4];
  double dPdx[3][3];  // not used, but required by postFcn->computeForces(...)..
  double dp1dxj[4][3]; // Gradient of the P1 basis functions
  d2w[0] = d2w[1] = d2w[2] = 0.0;
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) dPdx[i][j] = 0.0;
  for(int i=0;i<4;++i) for(int j=0;j<3;++j) dp1dxj[i][j] = 0.0;

  // -------------------------------------------------------------------------
  //   Loop through elements. Reconstruct the interface and compute the force 
  // -------------------------------------------------------------------------
  for (int iElem=0; iElem<elems.size(); iElem++) {
    map<int,int> global2local;
    for (int i=0; i<4; i++) {
      T[i] = elems[iElem][i];
      global2local[elems[iElem][i]] = i;
    }

    if(ghostPoints) //dp1dxj is needed only for the computation of viscous force.
      elems[iElem].computeGradientP1Function(X, dp1dxj);
    PolygonReconstructionData polygons[4];

    // find polygons.
    int numberOfPolygons = getPolygons(elems[iElem], LSS, polygons);
    assert(numberOfPolygons<=4);
    int index,vertex1,vertex2;

    // loop through each polygon
    for(int i=0; i < numberOfPolygons; ++i){

      PolygonReconstructionData& polygon=polygons[i];
      int nEdges = polygon.numberOfEdges;

      // get intersection information.
      LevelSetResult lsRes[nEdges];
      if(!nEdges) continue; //KW: why need this?
      for (int k=0; k<nEdges; ++k) {
        lsRes[k] = LSS.getLevelSetDataAtEdgeCenter(0,polygon.edge[k],polygon.edgeWithVertex[k][0]<polygon.edgeWithVertex[k][1]);
        if (lsRes[k].alpha<0) {fprintf(stderr,"Unable to get intersection results at edge center! Abort...\n"); exit(-1);}
      }

      // determine if we are applying a force from the "real fluid" or "ghost" side. From the "ghost" side, only a pressure
      // force (using pInfty) is applied, even for viscous flows.
      bool applyRealForce = true;
      for(int k=0; k<nEdges; k++)
        if(!LSS.isActive(0,polygon.edgeWithVertex[k][0])) {
          applyRealForce = false;
          break;
        }

      // get the correct states (mix of real and ghost) for this tet. This info is needed for viscous simulation only.
      double *v[4];
      if(ghostPoints && applyRealForce) {
        for(int k=0; k<4; k++)
          v[k] = V[T[k]];

        GhostPoint<dim> *gp;
        for(int k=0; k<nEdges; k++) {
          int vertex2 = polygon.edgeWithVertex[k][1];
          gp = ghostPoints->operator[](vertex2);
          v[global2local[vertex2]] = gp->getPrimitiveState(); 
        }
      }

      // get information at intersection points (3 for triangles, 4 for quadrangles).
      Vec3D Xinter[nEdges];
      double Vinter[nEdges][dim];
      Vec3D start_vertex;
      for(int k=0; k<nEdges; ++k) {
        vertex1 = polygon.edgeWithVertex[k][0];
        vertex2 = polygon.edgeWithVertex[k][1];
        if(k==0) {start_vertex[0]=X[vertex1][0]; start_vertex[1]=X[vertex1][1]; start_vertex[2]=X[vertex1][2];}
        int l = edges.find(polygon.edgeWithVertex[k][0], polygon.edgeWithVertex[k][1]);
        double alpha = lsRes[k].alpha;
        for(int m=0; m<3; m++)
          Xinter[k][m] = alpha*X[vertex1][m] + (1.0-alpha)*X[vertex2][m];

        fid_local[k] = fid?(*fid)[vertex1]:0;
        if(!applyRealForce)
          Vinter[k][4] = pInfty; //this is enough even for Tait (look at postFcn->computeEmbeddedForce(...))
        else {
          if(ghostPoints)  //viscous
            for(int m=0; m<dim; m++)
              Vinter[k][m] = alpha*v[global2local[vertex1]][m] + (1.0-alpha)*v[global2local[vertex2]][m];
          else { // inviscid
            if(vertex1<vertex2) // apply Wstarij
              for(int m=0; m<dim; m++)
                Vinter[k][m] = Wstarij[l][m];
            else //apply Wstarji
              for(int m=0; m<dim; m++)
                Vinter[k][m] = Wstarji[l][m];
          }
        }
      } 

      // compute force...
      double dist13,dist02;
      double *Xface[3], *Vface[3];
      int fid_face[3];
      double oneThird = 1.0/3.0;
      Vec3D nf;
      Vec3D fi0,fi1,fi2,fv,FNodal; // forces storage

      switch(nEdges){
        case 3: //got a triangle
          nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
          if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
          for(int m=0; m<3; m++) {
            Xface[m] = &Xinter[m][0];
            Vface[m] = &Vinter[m][0];
            fid_face[m] = fid_local[m]; 
          }
          postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v,&pInfty,
                                        fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);
          if(ghostPoints && applyRealForce) {
            fv *= oneThird;
            fi0+=fv;
            fi1+=fv;
            fi2+=fv;
          }
          sendLocalForce(fi0,lsRes[0],Fs);
          sendLocalForce(fi1,lsRes[1],Fs);
          sendLocalForce(fi2,lsRes[2],Fs);
          break;
        case 4: //got a quadrangle. cut it into two triangles.
          dist02 = (Xinter[2]-Xinter[0]).norm();
          dist13 = (Xinter[3]-Xinter[1]).norm();
          if(dist02<dist13){ // connect 0,2.
            //for triangle 012
            nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[1][0];  Vface[1] = &Vinter[1][0];  fid_face[1] = fid_local[1];
            Xface[2] = &Xinter[2][0];  Vface[2] = &Vinter[2][0];  fid_face[2] = fid_local[2];
            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v,&pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);
            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[1],Fs);
            sendLocalForce(fi2,lsRes[2],Fs);

            //for triangle 023
            nf = 0.5*(Xinter[2]-Xinter[0])^(Xinter[3]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[2][0];  Vface[1] = &Vinter[2][0];  fid_face[1] = fid_local[2];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];
            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v,&pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);
            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[2],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);

          }else{ // connect 1,3.
            //for triangle 123
            nf = 0.5*(Xinter[2]-Xinter[1])^(Xinter[3]-Xinter[1]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[1][0];  Vface[0] = &Vinter[1][0];  fid_face[0] = fid_local[1]; 
            Xface[1] = &Xinter[2][0];  Vface[1] = &Vinter[2][0];  fid_face[1] = fid_local[2];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];
            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v,&pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);
            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[1],Fs);
            sendLocalForce(fi1,lsRes[2],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);

            //for triangle 013
            nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[3]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[1][0];  Vface[1] = &Vinter[1][0];  fid_face[1] = fid_local[1];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];
            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v,&pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);
            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[1],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);
          }

          break;
#if 0
        default: fprintf(stderr,"ANOTHER PROBLEM!!!\n");break;
#endif
        default: break;
      }
    }
  }
}
//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeCVBasedForceLoad(int forceApp, int orderOfAccuracy, GeoState& geoState,
                                        SVec<double,3> &X, double (*Fs)[3], int sizeFs,
                                        LevelSetStructure &LSS, double pInfty, 
                                        SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji,
                                        SVec<double,dim> &V, Vec<GhostPoint<dim>*> *ghostPoints,
                                        PostFcn *postFcn, NodalGrad<dim, double> &ngrad, VarFcn *vf, Vec<int>* fid)
{
  if (forceApp!=1) 
    {fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); exit(-1);}

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  bool* masterFlag = edges.getMasterFlag();
  int (*ptr)[2];
  int i,j;
  ptr = edges.getPtr();
  bool iActive, jActive, intersect;
  // for the viscous force
  Vec3D vectorIJ, gradP;
  SVec<double,dim> gradX = ngrad.getX();
  SVec<double,dim> gradY = ngrad.getY();
  SVec<double,dim> gradZ = ngrad.getZ();
  double dudxj[3][3];

  for (int l=0; l<edges.size(); l++) {
    if (!masterFlag[l]) continue;
    i = ptr[l][0];
    j = ptr[l][1];
    if(LSS.withCracking()) {
      iActive = !LSS.isOccluded(0,i);
      jActive = !LSS.isOccluded(0,j);
    } else {
      iActive = LSS.isActive(0,i);
      jActive = LSS.isActive(0,j);
    }
    intersect = LSS.edgeIntersectsStructure(0,l);

    if (!iActive && !jActive) continue; //both inside structure
    if (!intersect) continue; 
    // now (i,j) must intersect the structure.

    double *v[2] = {V[i],V[j]};  

    if (iActive) {
      Vec3D flocal(0.0,0.0,0.0); 
      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
      if(ghostPoints) {// Viscous Simulation
        // Replace the state of node j by its corresponding ghost state
        /* This is useless. We use ngrad to compute the velocity gradient. 
        GhostPoint<dim> *gp;
        gp = ghostPoints->operator[](j);
        v[1] = gp->getPrimitiveState();
        */
        for(int m=0;m<3;++m) {
          vectorIJ[m] = X[j][m] - X[i][m];
        }
        gradP[0] = gradX[i][4];
        gradP[1] = gradY[i][4];
        gradP[2] = gradZ[i][4];
        double pp = vf->getPressure(v[0], fid?(*fid)[i]:0);
        flocal = (pp - pInfty + 0.5*(gradP*vectorIJ))*normal[l];
        // Compute Velocity Gradient
        // ***************** Method 1 **********************
        dudxj[0][0] = gradX[i][1]; dudxj[0][1] = gradY[i][1]; dudxj[0][2] = gradZ[i][1];
        dudxj[1][0] = gradX[i][2]; dudxj[1][1] = gradY[i][2]; dudxj[1][2] = gradZ[i][2];
        dudxj[2][0] = gradX[i][3]; dudxj[2][1] = gradY[i][3]; dudxj[2][2] = gradZ[i][3];
        // ***************** Method 2 **********************
        /*
        double norm = vectorIJ.norm();
        double deltaU = (v[0][1]-v[1][1])/(norm*norm);
        dudxj[0][0] = deltaU*vectorIJ[0];dudxj[0][1] = deltaU*vectorIJ[1];dudxj[0][2] = deltaU*vectorIJ[2];
        double deltaU = (v[0][2]-v[1][2])/(norm*norm);
        dudxj[1][0] = deltaU*vectorIJ[0];dudxj[1][1] = deltaU*vectorIJ[1];dudxj[1][2] = deltaU*vectorIJ[2];
        double deltaU = (v[0][3]-v[1][3])/(norm*norm);
        dudxj[2][0] = deltaU*vectorIJ[0];dudxj[2][1] = deltaU*vectorIJ[1];dudxj[2][2] = deltaU*vectorIJ[2];
        */
        flocal += postFcn->computeViscousForceCVBoundary(normal[l],v[0],dudxj);
      }
      else {
        double pp = vf->getPressure(Wstarij[l],fid?(*fid)[i]:0);
        flocal = (pp-pInfty)*normal[l];
      }
      sendLocalForce(flocal, lsRes, Fs);
    }
    if (jActive) {
      Vec3D flocal(0.0,0.0,0.0);
      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
      if(ghostPoints) {// Viscous Simulation
        // Replace the state of node i by its corresponding ghost state
        /* This is useless. We use ngrad to compute the velocity gradient.
        GhostPoint<dim> *gp;
        gp = ghostPoints->operator[](i);
        v[0] = gp->getPrimitiveState();
        */
        for(int m=0;m<3;++m) {
          vectorIJ[m] = X[j][m] - X[i][m];
        }
        gradP[0] = gradX[j][4];
        gradP[1] = gradY[j][4];
        gradP[2] = gradZ[j][4];
        double pp = vf->getPressure(v[1], fid?(*fid)[j]:0);
        // Minus, cause the normal points toward j
        flocal = -(pp - pInfty - 0.5*(gradP*vectorIJ))*normal[l];
        // Compute Velocity Gradient
        // ***************** Method 1 **********************
        dudxj[0][0] = gradX[j][1]; dudxj[0][1] = gradY[j][1]; dudxj[0][2] = gradZ[j][1];
        dudxj[1][0] = gradX[j][2]; dudxj[1][1] = gradY[j][2]; dudxj[1][2] = gradZ[j][2];
        dudxj[2][0] = gradX[j][3]; dudxj[2][1] = gradY[j][3]; dudxj[2][2] = gradZ[j][3];
        // ***************** Method 2 **********************
        /*
        double norm = vectorIJ.norm();
        double deltaU = -(v[1][1]-v[0][1])/(norm*norm);
        dudxj[0][0] = deltaU*vectorIJ[0];dudxj[0][1] = deltaU*vectorIJ[1];dudxj[0][2] = deltaU*vectorIJ[2];
        double deltaU = -(v[1][2]-v[0][2])/(norm*norm);
        dudxj[1][0] = deltaU*vectorIJ[0];dudxj[1][1] = deltaU*vectorIJ[1];dudxj[1][2] = deltaU*vectorIJ[2];
        double deltaU = -(v[1][3]-v[0][3])/(norm*norm);
        dudxj[2][0] = deltaU*vectorIJ[0];dudxj[2][1] = deltaU*vectorIJ[1];dudxj[2][2] = deltaU*vectorIJ[2];
        */
        // Minus, cause the normal points toward j
        flocal -= postFcn->computeViscousForceCVBoundary(normal[l],v[1],dudxj);
      }
      else {
        double pp = vf->getPressure(Wstarji[l],fid?(*fid)[j]:0);
        flocal = -(pp-pInfty)*normal[l];
      }
      sendLocalForce(flocal, lsRes, Fs);
    }
  }
}


//------------------------------------------------------------------------------

template<int dim>
void SubDomain::blur(SVec<double,dim> &U,SVec<double,dim> &U0, Vec<double>& weight)
{
  const Connectivity &nToN = *getNodeToNode(); 
  for(int currentNode=0;currentNode<numNodes();++currentNode) {
        
    for (int k = 0; k < dim; ++k) {
      U0[currentNode][k] = 0.0;
    }

    for(int j=0;j<nToN.num(currentNode);++j){
      int neighborNode=nToN[currentNode][j];
      for (int k = 0; k < dim; ++k) {
	U0[currentNode][k] += U[neighborNode][k];
      }
      
    }
    weight[currentNode] = (double)nToN.num(currentNode);
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::updateFluidIdFS2(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, SVec<bool,3> &poll, 
                                 Vec<int> &fluidId, bool *masterFlag)
{
  const Connectivity &Node2Node = *getNodeToNode();
  // ------- Determine status for grid-points swept by FS interface -------
  // Rule No.1: If this grid point is "occluded", set its status to "numPhases".
  // Rule No.2: If its "visible && !occluded && !swept" neighbors have the same status, use this one. 
  // Rule No.3: Otherwise, consider the sign of "PhiV". (PhiV should have been "blurred".)
  
  for(int i=0; i<PhiV.size(); i++) {
    bool swept = LSS.isSwept(0.0,i);
    bool occluded = LSS.isOccluded(0.0,i);

    //DEBUG
/*    int myNode = 300363;
    if(locToGlobNodeMap[i]+1==myNode){
      fprintf(stderr,"Node %d(%d), Sub %d. master = %d, swept = %d, occluded = %d, id = %d, phi = %e. Poll(%d,%d,%d)\n", myNode, i, globSubNum, masterFlag[i], swept, occluded, fluidId[i], PhiV[i][0], poll[i][0], poll[i][1], poll[i][2]);
      for(int j=0; j<Node2Node.num(i); j++) { 
        if(Node2Node[i][j]==i) continue;
        fprintf(stderr,"  Nei(%d,%d) on Sub %d--> GlobId(%d), occluded(%d), swept(%d), intersect(%d), id(%d), phi(%e).\n", myNode,i,globSubNum,
                          locToGlobNodeMap[Node2Node[i][j]], LSS.isOccluded(0.0,Node2Node[i][j]), LSS.isSwept(0.0,Node2Node[i][j]),
                          LSS.edgeIntersectsStructure(0.0,i,Node2Node[i][j]), fluidId[Node2Node[i][j]], PhiV[Node2Node[i][j]][0]); 
      }

    }
*/
    if(!swept) {//nothing to be done
      if(!poll[i][fluidId[i]]) fprintf(stderr,"TOO BAD!\n");
      continue;
    }

    if(occluded) { // Rule No.1
      fluidId[i] = LSS.numOfFluids();
      if(!poll[i][LSS.numOfFluids()]) fprintf(stderr,"TOO BAD TOO!\n");
      continue;
    }

    int count = (int)poll[i][0] + (int)poll[i][1] + (int)poll[i][2];
    switch (count) {
      case 0: //no info
        fprintf(stderr,"WARNING: More than one layer of nodes are swept in one step (near Node %d).\n", locToGlobNodeMap[i]+1);
        break;
      case 1: // Rule No.2
        for(int j=0; j<3; j++)
          if(poll[i][j]) fluidId[i] = j;
        break;
    } 

    if(count==1) //already applied Rule No.2
      continue;

    // Rule No.3
    bool done = false;
    for(int k=0; k<dimLS; k++)
      if(PhiV[i][k]>0.0) {
        fluidId[i] = k+1;
        done = true;
        break;
      }
    if(!done)
      fluidId[i] = 0;
  }
}

//------------------------------------------------------------------------------
/*
template<int dimLS> 
void SubDomain::updateFluidIdFS2(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, Vec<int> &fluidId, bool *masterFlag)
{
  // ------- Determine status for grid-points swept by FS interface -------
  // Rule No.1: If this grid point is "occluded", set its status to "numPhases".
  // Rule No.2: If its "visible && !occluded && !swept" neighbors have the same status, use this one. 
  // Rule No.3: Otherwise, consider the sign of "PhiV". (PhiV should have been "blurred".)

  const Connectivity &Node2Node = *getNodeToNode();
  
  for(int i=0; i<PhiV.size(); i++) {
    bool swept = LSS.isSwept(0.0,i);
    bool occluded = LSS.isOccluded(0.0,i);

    //DEBUG
    int myNode = 284121;
    if(locToGlobNodeMap[i]+1==myNode){
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      fprintf(stderr,"Node %d(%d), CPU = %d. master = %d, swept = %d, occluded = %d, id = %d, phi = %e.\n", myNode, i, myrank, masterFlag[i], swept, occluded, fluidId[i], PhiV[i][0]);
      for(int j=0; j<Node2Node.num(i); j++) { 
        if(Node2Node[i][j]==i) continue;
        fprintf(stderr,"  Nei(%d,%d) on CPU %d--> GlobId(%d), occluded(%d), swept(%d), intersect(%d), id(%d), phi(%e).\n", myNode,i,myrank,
                          locToGlobNodeMap[Node2Node[i][j]], LSS.isOccluded(0.0,Node2Node[i][j]), LSS.isSwept(0.0,Node2Node[i][j]),
                          LSS.edgeIntersectsStructure(0.0,i,Node2Node[i][j]), fluidId[Node2Node[i][j]], PhiV[Node2Node[i][j]][0]); 
      }

    }


    if(!swept) // nothing to be done.
      continue;
    if(!masterFlag[i]) {//will get Id from another subdomain (TODO: May need something more sophisticated... (KW).)
      fluidId[i] = 0;
      continue;
    }

    // Rule No.1
    if(occluded) {
      fluidId[i] = LSS.numOfFluids();
      continue;
    }

    // Rule No.2
    int myId = -1;
    int iNei, count = 0;
    bool consistent = false;
    for(int j=0; j<Node2Node.num(i); j++) {
      iNei = Node2Node[i][j];
      if(i==iNei)
        continue;
      if(LSS.isOccluded(0.0,iNei) || LSS.isSwept(0.0,iNei) || LSS.edgeIntersectsStructure(0.0,i,iNei))
        continue;
      count++;
      if(myId==-1) {
        myId = fluidId[iNei];
        consistent = true;
      } else if(myId!=fluidId[iNei]) {
        consistent = false;
        break;
      }
    }
    if(count==0)
      fprintf(stderr,"WARNING: More than one layer of nodes are swept in one step (near Node %d).\n", locToGlobNodeMap[i]+1);

    if(consistent) {
      fluidId[i] = myId;
      continue;
    }

    // Rule No.3
    bool done = false;
    for(int k=0; k<dimLS; k++)
      if(PhiV[i][k]>0.0) {
        fluidId[i] = k+1;
        done = true;
      }
    if(!done)
      fluidId[i] = 0;
  }
}
*/

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void SubDomain::debugMultiPhysics(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, Vec<int> &fluidId, SVec<double,dim> &U)
{
  for(int i=0; i<nodes.size(); i++) {
    switch (fluidId[i]) {
      case 0:
        if(PhiV[i][0]>=0.0)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][0]);
        if(LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][0]);
        break;
      case 1:
        if(PhiV[i][0]<=0.0)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][0]);
        if(LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][0]);
        break;
      case 2:
        if(fabs(PhiV[i][0])>1.0e-8)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][0]);
        if(!LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but NOT occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][0]);
        break;
      default:
        fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d!\n",  globSubNum, locToGlobNodeMap[i]+1, fluidId[i]);
    }
  }

/*
  int debugNode1 = 202747, debugNode2 = 202765;
  for(int i=0; i<nodes.size(); i++)
    if(locToGlobNodeMap[i]+1==debugNode1 || locToGlobNodeMap[i]+1==debugNode2) {
      fprintf(stderr,"%d: Node %d: Id = %d, PhiV = %e, occluded(%d), swept(%d), d2wall(%e), U(%e,%e,%e,%e,%e).\n",
              globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][0], LSS.isOccluded(0,i), LSS.isSwept(0,i), LSS.distToInterface(0,i),
              U[i][0], U[i][1], U[i][2], U[i][3], U[i][4]);
      for(int j=0; j<NodeToNode->num(i); j++) {
        int yId = (*NodeToNode)[i][j];
        if(yId==i) continue;
        fprintf(stderr,"%d:    Nei(%d)=%d: Id = %d, PhiV = %e, occluded(%d), swept(%d), d2wall(%e), X(%d).\n", globSubNum,
                locToGlobNodeMap[i]+1, locToGlobNodeMap[yId]+1, fluidId[yId], PhiV[yId][0], LSS.isOccluded(0,yId), LSS.isSwept(0,yId),
                LSS.distToInterface(0,yId), LSS.edgeIntersectsStructure(0,edges.findOnly(i,yId)));
        if((locToGlobNodeMap[i]+1==debugNode1 && locToGlobNodeMap[yId]+1==debugNode2) ||
           (locToGlobNodeMap[i]+1==debugNode2 && locToGlobNodeMap[yId]+1==debugNode1)) {
          LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0,edges.findOnly(i,yId),i<yId);
          fprintf(stderr,"** (%d,%d): alpha(%e), Tr(%d,%d,%d), Normal(%e,%e,%e), xi(%e,%e,%e).\n",locToGlobNodeMap[i]+1,locToGlobNodeMap[yId]+1,
                  resij.alpha, resij.trNodes[0], resij.trNodes[1], resij.trNodes[2], resij.gradPhi[0], resij.gradPhi[1], resij.gradPhi[2],
                  resij.xi[0], resij.xi[1], resij.xi[2]);
        }
      }
    }
*/
}

//------------------------------------------------------------------------------

template<int dim, class Obj>
void SubDomain::integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
				  int npt) 
{
  elems.integrateFunction(obj,X,V,F,npt);
}

template<int dim> 
void SubDomain::interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, 
                                    const std::vector<Vec3D>& locs, double (*sol)[dim],
                                    int* status,int* last,int* nid,
                                    LevelSetStructure* LSS, Vec<GhostPoint<dim>*>* ghostPoints,
                                    VarFcn *varFcn) {

  elems.interpolateSolution(X,U,locs,sol,status,last,LSS,ghostPoints,varFcn);
  for (int i = 0; i < locs.size(); ++i) {
    int eid = last[i],nn;
    Elem& E = elems[eid];
    double mindist = std::numeric_limits<double>::max(),dst;
    for (int j = 0; j < E.numNodes(); ++j) {
      nn = E.nodeNum(j);
      dst = sqrt((X[nn][0]-locs[i][0])*(X[nn][0]-locs[i][0])+
                 (X[nn][1]-locs[i][1])*(X[nn][1]-locs[i][1])+
                 (X[nn][2]-locs[i][2])*(X[nn][2]-locs[i][2]));
      if (dst < mindist) {
        mindist = dst;
        nid[i] = nn;
      }
    }
  }
}

template<int dim>
void SubDomain::interpolatePhiSolution(SVec<double,3>& X, SVec<double,dim>& U,
                                    const std::vector<Vec3D>& locs, double (*sol)[dim],
                                    int* status,int* last,int* nid) {


  elems.interpolateSolution(X,U,locs,sol,status,last);
  for (int i = 0; i < locs.size(); ++i) {
    int eid = last[i],nn;
    Elem& E = elems[eid];
    double mindist = std::numeric_limits<double>::max(),dst;
    for (int j = 0; j < E.numNodes(); ++j) {
      nn = E.nodeNum(j);
      dst = sqrt((X[nn][0]-locs[i][0])*(X[nn][0]-locs[i][0])+
                 (X[nn][1]-locs[i][1])*(X[nn][1]-locs[i][1])+
                 (X[nn][2]-locs[i][2])*(X[nn][2]-locs[i][2]));
      if (dst < mindist) {
        mindist = dst;
        nid[i] = nn;
      }
    }
  }
}


//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::findCutCells(int lsdim,
			     SVec<double,dimLS>& phi,Vec<int>& cutStatus,
			     SVec<double,3>& X) {  

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    int i = edgePtr[l][0],j = edgePtr[l][1];
    double* phii = phi[i],*phij = phi[j];

    if (phii[lsdim]*phij[lsdim] > 0.0)
      continue;
    
    double s = -phii[lsdim]/(phij[lsdim]-phii[lsdim]);
    if (s < 0.1)
      cutStatus[i] = 1;
    else if (s > 0.9)
    cutStatus[j] = 1;
    //if (cutStatus[i] || cutStatus[j])
    //  std::cout << i << " " << j << " " << 
    //	phii[lsdim] << " " << phij[lsdim] << " (" << X[i][0] << " " << X[j][0] << ")" << std::endl;
      
  }

}

template <int dim>
void SubDomain::setCutCellFlags(int lsdim, Vec<int>& cutStatus) {

  higherOrderMF->template setCutCellFlags<dim>(lsdim,cutStatus);
}

inline
int SubDomain::getNumCutCells() {

  return higherOrderMF->getNumCutCells();
}


template <int dim>
void SubDomain::collectCutCellData(SVec<double,dim>* cutCell[2],
				   NodalGrad<dim,double>* cutGrad[2],
				   Vec<int>* counts[2],
				   NodalGrad<dim,double>& grad, Vec<int>& fluidId,
				   SVec<double,dim>& V,SVec<double,3>& X) {

  int (*edgePtr)[2] = edges.getPtr();
  bool* masterFlag = edges.getMasterFlag();

  bool cut;
  int fid,cutnode,numset[2] = {0,0};
  for (int l=0; l<edges.size(); ++l) {

    cut = false;
    int i = edgePtr[l][0],j = edgePtr[l][1];
    double Vext[dim],newgrad[dim][3];
    if (!masterFlag[l])
      continue;

    if (higherOrderMF->isCellCut(i)) {

      if (higherOrderMF->isCellCut(j)) {
	// nothing to do, this edge is between cut cells
	continue;
      }
      
      fid = fluidId[j];
      cutnode = i;
      for (int k = 0; k < dim; ++k) {
	Vext[k] = V[j][k] + grad.getX()[j][k]*(X[i][0]-X[j][0])+ 
	  grad.getY()[j][k]*(X[i][1]-X[j][1])+
	  grad.getZ()[j][k]*(X[i][2]-X[j][2]);
	newgrad[k][0] = grad.getX()[j][k];
	newgrad[k][1] = grad.getY()[j][k];
	newgrad[k][2] = grad.getZ()[j][k];
      }
    } else if (higherOrderMF->isCellCut(j)) {

      fid = fluidId[i];
      cutnode = j;
      for (int k = 0; k < dim; ++k) {
	Vext[k] = V[i][k] + grad.getX()[i][k]*(X[j][0]-X[i][0])+ 
	  grad.getY()[i][k]*(X[j][1]-X[i][1])+
	  grad.getZ()[i][k]*(X[j][2]-X[i][2]);
	newgrad[k][0] = grad.getX()[i][k];
	newgrad[k][1] = grad.getY()[i][k];
	newgrad[k][2] = grad.getZ()[i][k];
      }
      
    } else 
      continue;

    int fididx = (fid == 0? 0 : 1);
    ++(counts[fididx]->operator[](cutnode));
    //std::cout << i << " " << j << " " << fididx << " " << l << " " <<
    //  X[i][0] << " " << X[j][0] << std::endl;
    for (int k = 0; k < dim; ++k) {
      cutCell[fididx]->operator[](cutnode)[k] += Vext[k];
      cutGrad[fididx]->getX()[cutnode][k] += newgrad[k][0];
      cutGrad[fididx]->getY()[cutnode][k] += newgrad[k][1];
      cutGrad[fididx]->getZ()[cutnode][k] += newgrad[k][2];
    } 
  }
  
}
/*
template <int dim>
void SubDomain::collectCutCellData2() {

  bool* masterFlag = edges.getMasterFlag();

  int* count[2] = {new int[getNumCutCells()],new int[getNumCutCells()]};
  int* nodes[2] = {new int[getNumCutCells()],new int[getNumCutCells()]};
  bool cut;
  int fid,cutnode,numset[2] = {0,0};
  for (int l=0; l<edges.size(); ++l) {

    cut = false;
    int i = edgePtr[l][0],j = edgePtr[l][1];
    double Vext[dim],newgrad[dim][3];
    if (!masterFlag[l])
      continue;

    if (higherOrderMF->isCellCut(i)) {

      if (higherOrderMF->isCellCut(j)) {
	// nothing to do, this edge is between cut cells
	continue;
      }
      
      fid = fluidId[i];
      cutnode = i;
      for (int k = 0; k < dim; ++k) {
	Vext[k] = V[j][k] + grad.getX()[j][k]*(X[i][0]-X[j][0])+ 
	  grad.getY()[j][k]*(X[i][1]-X[j][1])+
	  grad.getZ()[j][k]*(X[i][2]-X[j][2]);
	newgrad[k][0] = grad.getX()[j][k];
	newgrad[k][1] = grad.getY()[j][k];
	newgrad[k][2] = grad.getZ()[j][k];
      }
    } else if (higherOrderMF->isCellCut(j)) {

      fid = fluidId[j];
      cutnode = j;
      for (int k = 0; k < dim; ++k) {
	Vext[k] = V[i][k] + grad.getX()[i][k]*(X[j][0]-X[i][0])+ 
	  grad.getY()[i][k]*(X[j][1]-X[i][1])+
	  grad.getZ()[i][k]*(X[j][2]-X[i][2]);
	newgrad[k][0] = grad.getX()[i][k];
	newgrad[k][1] = grad.getY()[i][k];
	newgrad[k][2] = grad.getZ()[i][k];
      }
      
    } else 
      continue;

    int setloc = -1;
    int idx = (fid == 0 ? 0 : 1);
    for (int m = 0; m < numset[idx]; ++m) {
      
      if (nodes[idx][m] == cutnode) {
	setloc = m;
	break;
      }
    }
    
    if (setloc == -1) {
      nodes[idx][numset] = cutnode;
      setloc = numset;
      ++numset;
    }
    
    count[idx][setloc]++;

    for (int k = 0; k < dim; ++k) {

      
    }
    
    
  }
}
*/

template <int dim>
void SubDomain::storeCutCellData(SVec<double,dim>* cutCell[2],
				 NodalGrad<dim,double>* cutGrad[2],
				 Vec<int>* counts[2], Vec<int>& fluidId) {
  
  higherOrderMF->template storeCutCellData<dim>(cutCell,cutGrad, counts);

}

template<int dim>
void SubDomain::setCutCellData(SVec<double,dim>& V, Vec<int>& fid) {

  higherOrderMF->template setCutCellData<dim>(V, fid);
}

// Functions to compute the error (that is, the difference between two state vectors)
template <int dim>
void SubDomain::computeL1Error(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, double error[dim]) {

  for (int k = 0; k < dim; ++k)
    error[k] = 0.0;

  for(int i=0; i<nodes.size(); i++) {

    if (nodeFlag[i]) {
      
      for (int k = 0; k < dim; ++k) {
	
	error[k] += fabs(U[i][k]-Uexact[i][k]);
      }
    }
  }
}

template <int dim>
void SubDomain::computeLInfError(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, double error[dim]) {

  for (int k = 0; k < dim; ++k)
    error[k] = 0.0;

  for(int i=0; i<nodes.size(); i++) {

    if (nodeFlag[i]) {
      
      for (int k = 0; k < dim; ++k) {
	
	error[k] = max(error[k],fabs(U[i][k]-Uexact[i][k]));
      }
    }
  }
}
  
