#include <BcFcn.h>
#include <FluxFcn.h>
#include <RecFcnDesc.h>
#include <DistNodalGrad.h>
#include <DistEdgeGrad.h>
#include <DistExtrapolation.h>
#include <DistMacroCell.h>
#include <DistBcData.h>
#include <DistExactRiemannSolver.h>
#include <SubDomain.h>
#include <DistGeoState.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <Communicator.h>
#include <PostFcn.h>
#include <LowMachPrec.h>
#include <LevelSet/FluidTypeCriterion.h>
#include <GeoState.h>
#include <NodalGrad.h>
#include <FluidSelector.h>

#include <cstdio>
#include <cmath>
#include <unistd.h>

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTimeStep(double cfl, double viscous, FemEquationTerm *fet, VarFcn *varFcn, DistGeoState &geoState,
			     DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
			     DistVec<double> &dt, DistVec<double> &idti, DistVec<double> &idtv,
			     DistVec<double> &irey, TimeLowMachPrec &tprec, SpatialLowMachPrec &sprec)
{

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), X(iSub), V(iSub), dt(iSub), idti(iSub), idtv(iSub), tprec);
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(idti.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(idti.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(idtv.subData(iSub)));

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(idtv.subData(iSub)));

// Included (MB)
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*idtimev) = idtv.subData(iSub);
    double (*idtimei) = idti.subData(iSub);
    double (*dtime) = dt.subData(iSub);
    double (*ireynolds) = irey.subData(iSub);
    double (*volume) = ctrlVol.subData(iSub);
    for (int i = 0; i < ctrlVol.subSize(iSub); ++i) {
      //   idtimev[i] = idtimev[i] / volume[i];
      dtime[i] = cfl *volume[i]/(-1.0*idtimei[i] + viscous*idtimev[i]);
      ireynolds[i] = -sprec.getViscousRatio()*idtimev[i] / idtimei[i];
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfInvReynolds(FemEquationTerm *fet, VarFcn *varFcn, DistGeoState &geoState,
			     DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol,
			     DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim> &dV,
			     DistVec<double> &idti, DistVec<double> &dIdti, DistVec<double> &idtv, DistVec<double> &dIdtv,
			     DistVec<double> &dIrey, double dMach, TimeLowMachPrec&tprec, SpatialLowMachPrec &sprec)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeDerivativeOfTimeStep(fet, varFcn, geoState(iSub), X(iSub), dX(iSub), V(iSub), dV(iSub), dIdti(iSub), dIdtv(iSub), dMach, tprec);
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dIdti.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dIdti.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dIdtv.subData(iSub)));

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dIdtv.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*idtimev) = idtv.subData(iSub);
    double (*dIdtimev) = dIdtv.subData(iSub);
    double (*idtimei) = idti.subData(iSub);
    double (*dIdtimei) = dIdti.subData(iSub);
    double (*dIreynolds) = dIrey.subData(iSub);
    double (*volume) = ctrlVol.subData(iSub);
    double (*dVolume) = dCtrlVol.subData(iSub);
    for (int i = 0; i < ctrlVol.subSize(iSub); ++i) {
      dIdtimev[i] = (dIdtimev[i]*volume[i] - (idtimev[i]*volume[i])*dVolume[i]) / (volume[i]*volume[i]);
      dIreynolds[i] = -sprec.getViscousRatio()*(dIdtimev[i]*idtimei[i] - idtimev[i]*dIdtimei[i]) / (idtimei[i]*idtimei[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTimeStep(double cfl, double viscous, FemEquationTerm *fet, VarFcn *varFcn, DistGeoState &geoState,
                             DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
                             DistVec<double> &dt, DistVec<double> &idti, DistVec<double> &idtv,
			     TimeLowMachPrec &tprec, DistVec<int> &fluidId, DistVec<double>* umax)
{

  int iSub;

  if (umax)
    *umax = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    if (umax)
      subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), V(iSub), dt(iSub), idti(iSub), idtv(iSub), tprec, fluidId(iSub),&(*umax)(iSub));
    else
      subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), V(iSub), dt(iSub), idti(iSub), idtv(iSub), tprec, fluidId(iSub));
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dt.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dt.subData(iSub)));

  if (umax) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {  
      subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(umax->subData(iSub)));
    }

    volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(umax->subData(iSub)));

  }

  dt = -cfl * ctrlVol / dt;

  if (umax)
    *umax = -0.333f * ctrlVol / (*umax);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsLeastSquares(DistSVec<double,3> &X,
					  DistSVec<double,6> &R,
					  DistSVec<Scalar,dim> &var,
					  DistSVec<Scalar,dim> &ddx,
					  DistSVec<Scalar,dim> &ddy,
					  DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsLeastSquares(X(iSub), R(iSub), var(iSub),
						  ddx(iSub), ddy(iSub), ddz(iSub));

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void Domain::computeDerivativeOfGradientsLeastSquares(DistSVec<double,3> &X, DistSVec<double,3> &dX,
					  DistSVec<double,6> &R, DistSVec<double,6> &dR,
					  DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &dvar, DistSVec<Scalar,dim> &dddx,
					  DistSVec<Scalar,dim> &dddy, DistSVec<Scalar,dim> &dddz)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGradientsLeastSquares(X(iSub), dX(iSub), R(iSub), dR(iSub), var(iSub), dvar(iSub),
						  dddx(iSub), dddy(iSub), dddz(iSub));

  timer->addNodalGradTime(t0);

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow and FSI)
template<int dim, class Scalar>
void Domain::computeGradientsLeastSquares(DistSVec<double,3> &X,
                                          DistVec<int> &fluidId,
                                          DistSVec<double,6> &R,
                                          DistSVec<Scalar,dim> &var,
                                          DistSVec<Scalar,dim> &ddx,
                                          DistSVec<Scalar,dim> &ddy,
                                          DistSVec<Scalar,dim> &ddz,
                                          bool linFSI, DistLevelSetStructure *distLSS)
{

  if(distLSS) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, &((*distLSS)(iSub)));
  } else {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, 0);
  }

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsGalerkin(DistVec<double> &ctrlVol, DistSVec<double,3> &wii,
				      DistSVec<double,3> &wij, DistSVec<double,3> &wji,
				      DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &ddx,
				      DistSVec<Scalar,dim> &ddy, DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsGalerkin(ctrlVol(iSub), wii(iSub), wij(iSub), wji(iSub),
					      var(iSub), ddx(iSub), ddy(iSub), ddz(iSub));

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void Domain::computeDerivativeOfGradientsGalerkin(DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
                      DistSVec<double,3> &wii, DistSVec<double,3> &wij,
                      DistSVec<double,3> &wji, DistSVec<double,3> &dwii,
                      DistSVec<double,3> &dwij, DistSVec<double,3> &dwji,
                      DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &dvar, DistSVec<Scalar,dim> &dddx,
                      DistSVec<Scalar,dim> &dddy, DistSVec<Scalar,dim> &dddz)
{

  DistSVec<Scalar,dim> ddx(getNodeDistInfo());
  DistSVec<Scalar,dim> ddy(getNodeDistInfo());
  DistSVec<Scalar,dim> ddz(getNodeDistInfo());

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGradientsGalerkin(ctrlVol(iSub), dCtrlVol(iSub), wii(iSub), wij(iSub), wji(iSub),
                                       dwii(iSub), dwij(iSub), dwji(iSub), var(iSub), dvar(iSub), ddx(iSub), ddy(iSub), ddz(iSub), dddx(iSub), dddy(iSub), dddz(iSub));

  timer->addNodalGradTime(t0);

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsGalerkinT(DistVec<double> &ctrlVol,
                DistSVec<double,3> &wii, DistSVec<double,3> &wij,
                DistSVec<double,3> &wji, DistSVec<Scalar,dim> &var,
                DistSVec<Scalar,dim> &var1, DistSVec<Scalar,dim> &var2,
                DistSVec<Scalar,dim> &ddx, DistSVec<Scalar,dim> &ddy,
                DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsGalerkinT(ctrlVol(iSub), wii(iSub),
                wij(iSub), wji(iSub), var(iSub), var1(iSub), var2(iSub),
                ddx(iSub), ddy(iSub), ddz(iSub));


  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *recFcn, DistSVec<double,3> &X,
				    DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
				    DistSVec<double,dim> &dVdx, DistSVec<double,dim> &dVdy,
				    DistSVec<double,dim> &dVdz, DistSVec<double,dim> &Vmin,
				    DistSVec<double,dim> &Vmax, DistSVec<double,dim> &phi)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMinMaxStencilValues(V(iSub), Vmin(iSub), Vmax(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, Vmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmax.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*vecPat, Vmax.subData(iSub));
    subDomain[iSub]->computeMultiDimLimiter(recFcn, X(iSub), ctrlVol(iSub),
					    V(iSub), dVdx(iSub), dVdy(iSub), dVdz(iSub),
					    Vmin(iSub), Vmax(iSub), phi(iSub));
    subDomain[iSub]->sndData(*vecPat, phi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, phi.subData(iSub));

    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    for (int i=0; i<phi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	locdVdx[i][k] *= locphi[i][k];
	locdVdy[i][k] *= locphi[i][k];
	locdVdz[i][k] *= locphi[i][k];
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfMultiDimLimiter(RecFcnLtdMultiDim<dim> *recFcn, DistSVec<double,3> &X, DistSVec<double,3> &dX,
				    DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim> &dV,
				    DistSVec<double,dim> &dVdx, DistSVec<double,dim> &dVdy, DistSVec<double,dim> &dVdz,
				    DistSVec<double,dim> &ddVdx, DistSVec<double,dim> &ddVdy, DistSVec<double,dim> &ddVdz,
                    DistSVec<double,dim> &Vmin, DistSVec<double,dim> &dVmin, DistSVec<double,dim> &Vmax,
                    DistSVec<double,dim> &dVmax, DistSVec<double,dim> &phi, DistSVec<double,dim> &dphi)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeDerivativeOfMinMaxStencilValues(V(iSub), dV(iSub), Vmin(iSub), dVmin(iSub), Vmax(iSub), dVmax(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, Vmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmax.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*vecPat, Vmax.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dVmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, dVmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dVmax.subData(iSub));
  }

  vecPat->exchange();


#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*vecPat, dVmax.subData(iSub));
    subDomain[iSub]->computeDerivativeOfMultiDimLimiter(recFcn, X(iSub), dX(iSub), ctrlVol(iSub), dCtrlVol(iSub),
					    V(iSub), dV(iSub), dVdx(iSub), dVdy(iSub), dVdz(iSub), ddVdx(iSub), ddVdy(iSub), ddVdz(iSub),
					    Vmin(iSub), dVmin(iSub), Vmax(iSub), dVmax(iSub), phi(iSub), dphi(iSub));
    subDomain[iSub]->sndData(*vecPat, phi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, phi.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dphi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, dphi.subData(iSub));

    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    double (*dlocphi)[dim] = dphi.subData(iSub);
    double (*dlocdVdx)[dim] = ddVdx.subData(iSub);
    double (*dlocdVdy)[dim] = ddVdy.subData(iSub);
    double (*dlocdVdz)[dim] = ddVdz.subData(iSub);

    for (int i=0; i<dphi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	dlocdVdx[i][k] *= locphi[i][k];
	dlocdVdy[i][k] *= locphi[i][k];
	dlocdVdz[i][k] *= locphi[i][k];
	dlocdVdx[i][k] += locdVdx[i][k] * dlocphi[i][k];
	dlocdVdy[i][k] += locdVdy[i][k] * dlocphi[i][k];
	dlocdVdz[i][k] += locdVdz[i][k] * dlocphi[i][k];
      }
    }
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    for (int i=0; i<phi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	locdVdx[i][k] *= locphi[i][k];
	locdVdy[i][k] *= locphi[i][k];
	locdVdz[i][k] *= locphi[i][k];
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computePressureSensor(double threshold, DistSVec<double,3>& X,
				   DistSVec<double,dim>& V, DistSVec<double,dim>& dVdx,
				   DistSVec<double,dim>& dVdy, DistSVec<double,dim>& dVdz,
				   DistSVec<double,3>& sensor, DistVec<double>& sigma)
{

  const int nsmooth = 2;
  const double gamma = 0.01;
  const double eps = 0.5;
  const double omeps = 1.0 - eps;

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computePressureSensor(X(iSub), V(iSub), dVdx(iSub), dVdy(iSub),
					   dVdz(iSub), sensor(iSub));
    subDomain[iSub]->sndData(*vec3DPat, sensor.subData(iSub));
  }
  vec3DPat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec3DPat, sensor.subData(iSub));
    double (*s)[3] = sensor.subData(iSub);
    double* sig = sigma.subData(iSub);
    for (int i=0; i<sigma.subSize(iSub); ++i) {
      double alpha = gamma * s[i][2] / (s[i][0] + s[i][1] + s[i][2]);
      sig[i] = s[i][0] / (s[i][1] + alpha * s[i][2]);
    }
  }

  for (int k=0; k<nsmooth; ++k) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeSmoothedSensor(X(iSub), sigma(iSub), sensor(iSub));
      subDomain[iSub]->sndData(*vec3DPat, sensor.subData(iSub));
    }
    vec3DPat->exchange();
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*vec3DPat, sensor.subData(iSub));
      double (*s)[3] = sensor.subData(iSub);
      double* sig = sigma.subData(iSub);
      for (int i=0; i<sigma.subSize(iSub); ++i)
	sig[i] = eps * sig[i] + omeps * s[i][0] / s[i][1];
    }
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);
    double* sig = sigma.subData(iSub);
    for (int i=0; i<sigma.subSize(iSub); ++i) {
      if (sig[i] >= threshold) {
	sig[i] = 1.0;
	for (int k=0; k<dim; ++k) {
	  locdVdx[i][k] = 0.0;
	  locdVdy[i][k] = 0.0;
	  locdVdz[i][k] = 0.0;
	}
      }
      else
	sig[i] = 0.0;
    }
  }

}
//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol, DistVec<double>& irey,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  //KW&AM: TODO: should add RR as a member.
  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeTerm(irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                     geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                     legrad, (*RR)(iSub), (*tag)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeTerm(irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                 geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                 legrad, (*RR)(iSub), (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeTerm(DistExactRiemannSolver<dim>& riemann,
                                     DistVec<double> &ctrlVol, DistVec<double>& irey,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub), irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                     geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                     legrad, (*RR)(iSub), (*tag)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, V); // bug: this is for one-phase flow!
      //ngrad.compute(geoState.getConfig(), X, ctrlVol, fluidId, V); //where is fluidId?
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub), irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                 geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                 legrad, (*RR)(iSub), (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfFiniteVolumeTerm(DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
				     DistVec<double>& irey, DistVec<double>& dIrey,
				     FluxFcn** fluxFcn, RecFcn* recFcn,
				     DistBcData<dim>& bcData, DistGeoState& geoState,
				     DistSVec<double,3>& X, DistSVec<double,3>& dX, DistSVec<double,dim>& V, DistSVec<double,dim>& dV,
				     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, double dMach,
				     DistSVec<double,dim>& dF)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    subDomain[iSub]->computeDerivativeOfFiniteVolumeTerm(irey(iSub), dIrey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
					     X(iSub), dX(iSub), V(iSub), dV(iSub), ngrad(iSub), legrad, dMach, dF(iSub));
    subDomain[iSub]->sndData(*vecPat, dF.subData(iSub));

  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, dF.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     FluidSelector &fluidSelector,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistNodalGrad<dimLS>& ngradLS,
                                     DistSVec<double,dim>& R, int it,
                                     int failsafe, int rshift)
{
  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual
  DistVec<int> &FluidId(*(fluidSelector.fluidId));

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    Vec<int> &fluidId = FluidId(iSub);
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                             fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                             X(iSub), V(iSub), fluidId,
                                             fluidSelector, ngrad(iSub),
                                             legrad,  ngradLS(iSub), (*RR)(iSub), it,
                                             (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, FluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        Vec<int> &fluidId = FluidId(iSub);
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                     fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                     X(iSub), V(iSub), fluidId,
                                     fluidSelector, ngrad(iSub),
                                     legrad, ngradLS(iSub), (*RR)(iSub), it,
                                     (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

// subdomain communication for riemann update values (cf ExactRiemannSolver.h)
  if(it == 1){
    DistSVec<double,dim> *rupdate = riemann.getRiemannUpdate();
    DistVec<double> *weight= riemann.getRiemannWeight();
    assemble(vecPat,*rupdate);
    assemble(volPat,*weight);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol, DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn, DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V, DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji,
                                     DistLevelSetStructure *distLSS, bool linRecAtInterface, FluidSelector &fluidSelector, int Nriemann,
                                     DistSVec<double,3> *Nsbar, DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistNodalGrad<dimLS>& ngradLS, DistSVec<double,dim>& R, int it, int failsafe, int rshift)
{
 double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual
  DistVec<int> &FluidId(*(fluidSelector.fluidId));

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    SVec<double,3>* nsbar = (Nsbar) ? &((*Nsbar)(iSub)) : 0;
    Vec<int> &fluidId = FluidId(iSub);
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                             fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                             X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), (*distLSS)(iSub), 
                                             linRecAtInterface, fluidId, Nriemann, nsbar,
                                             fluidSelector, ngrad(iSub),
                                             legrad,  ngradLS(iSub), (*RR)(iSub), it,
                                             (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  // failsafe
  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, FluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        SVec<double,3>* nsbar = (Nsbar) ? &((*Nsbar)(iSub)) : 0;
        Vec<int> &fluidId = FluidId(iSub);
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                             fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                             X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), (*distLSS)(iSub), 
                                             linRecAtInterface, fluidId, Nriemann, nsbar,
                                             fluidSelector, ngrad(iSub),
                                             legrad,  ngradLS(iSub), (*RR)(iSub), it,
                                             (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;
  if (RR) delete(RR); // delete temp residual

// subdomain communication for riemann update values (cf ExactRiemannSolver.h)
  if(it == 1){
    DistSVec<double,dim> *rupdate = riemann.getRiemannUpdate();
    DistVec<double> *weight= riemann.getRiemannWeight();
    assemble(vecPat,*rupdate);
    assemble(volPat,*weight);
  }

  timer->addFiniteVolumeTermTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji,
                                     DistLevelSetStructure *LSS, bool linRecAtInterface, DistVec<int> &fluidId, 
                                     int Nriemann, DistSVec<double,3> *Nsbar, DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, int it,
                                     int failsafe, int rshift)
{
  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                             fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                             X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), (*LSS)(iSub),
                                             linRecAtInterface, fluidId(iSub), Nriemann, (Nsbar) ? &((*Nsbar)(iSub)) : 0,  ngrad(iSub),
                                             legrad, (*RR)(iSub), it,
                                             (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, fluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
                                     fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                     X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), (*LSS)(iSub),
                                     linRecAtInterface, fluidId(iSub), Nriemann, (Nsbar) ? &((*Nsbar)(iSub)) : 0, ngrad(iSub),
                                     legrad, (*RR)(iSub), it,
                                     (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistNodalGrad<dim>& ngrad, DistNodalGrad<dimLS>& ngradLS,
                                     DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dimLS>& Phi, DistSVec<double,dimLS> &PhiF,
                                     DistLevelSetStructure *distLSS)
{
  double t0 = timer->getTime();
  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    LevelSetStructure* LSS = (distLSS) ? &((*distLSS)(iSub)) : 0;
    subDomain[iSub]->computeFiniteVolumeTermLS(fluxFcn, recFcn, recFcnLS, bcData(iSub),
                                               geoState(iSub),
                                               X(iSub), V(iSub), ngrad(iSub), ngradLS(iSub),
                                               legrad, Phi(iSub),PhiF(iSub), LSS);
    subDomain[iSub]->sndData(*phiVecPat, PhiF.subData(iSub));
  }

  phiVecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*phiVecPat, PhiF.subData(iSub));

  timer->addLSFiniteVolumeTermTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeBarTerm(DistVec<double> &ctrlVol,
                                        DistVec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                        DistBcData<dim>& bcData, DistGeoState& geoState,
                                        DistSVec<double,3>& X, DistMacroCellSet *macroCells,
                                        DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                                        DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                        DistSVec<double,dim>& RBar, int scopeDepth1, int scopeDepth2,
                                        int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tagBar) {
     tagBar = new DistSVec<int,2>(getNodeDistInfo());
     *tagBar = 0;
  }

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;
  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeBar_Step1(irey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                          X(iSub), VBar(iSub), ngrad(iSub), legrad, (*Sigma)(iSub),
                                                          (*tagBar)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *Sigma = 0.0; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tagBar).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tagBar).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tagBar)(iSub));

      ngrad.fix(*tagBar);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, VBar);
      ngrad.limit(recFcn, X, ctrlVol, VBar);

      if (egrad) egrad->fix(*tagBar);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeBar_Step1(irey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                      X(iSub), VBar(iSub), ngrad(iSub), legrad, (*Sigma)(iSub),
                                                      (*tagBar)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tagBar = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*Sigma).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*Sigma).subData(iSub));


#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];

    for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

     subDomain[iSub]->computeFiniteVolumeBar_Step2(macCells, volRatio(iSub),
                                   (*Sigma)(iSub), RBar(iSub), scopeDepth1);
     subDomain[iSub]->sndData(*vecPat, RBar.subData(iSub));
     delete [] macCells;
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, RBar.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

  delete (Sigma);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             DistGeoState &geoState, DistVec<double> &irey,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol,
                                             DistSVec<double,dim> &V, DistMat<Scalar,neq> &A)
{
  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);


  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq, int dimLS>
void Domain::computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim> &riemann,
                                             FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             DistGeoState &geoState,
                                             DistNodalGrad<dim> &ngrad, DistNodalGrad<dimLS> &ngradLS,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol,
                                             DistSVec<double,dim> &V, DistMat<Scalar,neq> &A,
                                             FluidSelector &fluidSelector)
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
    fprintf(stdout, "with inletRhsPat\n");
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     ngrad(iSub), ngradLS(iSub),
                                                     X(iSub),
                                                     ctrlVol(iSub), V(iSub), A(iSub),
                                                     fluidSelector, 
                                                     (*(fluidSelector.fluidId))(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     ngrad(iSub), ngradLS(iSub),
                                                     X(iSub),
                                                     ctrlVol(iSub), V(iSub), A(iSub),
                                                     fluidSelector, 
                                                     (*(fluidSelector.fluidId))(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);
  }
}

template<class Scalar,int dim,int neq>
void Domain::computeJacobianFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                             DistExactRiemannSolver<dim> &riemann,
                                             FluxFcn** fluxFcn,
                                             DistBcData<dim>& bcData, DistGeoState& geoState,
                                             DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                             DistLevelSetStructure *LSS, DistVec<int> &fluidId, 
                                             int Nriemann, DistSVec<double,3> *Nsbar,
                                             DistMat<Scalar,neq>& A,DistVec<double>& irey) 
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();
 
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      SVec<double,3>* nsbar = (Nsbar) ? &((*Nsbar)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     X(iSub), V(iSub),ctrlVol(iSub),
                                                     (*LSS)(iSub),
                                                     fluidId(iSub),Nriemann,
                                                     nsbar,A(iSub),irey(iSub));
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
}
  
template<int dim, class Scalar, int neq, int dimLS>
void Domain::computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim>& riemann,
                                             FluxFcn** fluxFcn, 
                                             DistBcData<dim>& bcData, DistGeoState& geoState,
                                             DistSVec<double,3>& X, DistSVec<double,dim>& V,DistVec<double>& ctrlVol,
                                             DistNodalGrad<dimLS> &ngradLS,
                                             DistLevelSetStructure *LSS,
                                             int Nriemann, DistSVec<double,3>* Nsbar,
                                             FluidSelector &fluidSelector,
                                             DistMat<Scalar,neq>& A) {

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();
 
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      SVec<double,3>* nsbar = (Nsbar) ? &((*Nsbar)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     X(iSub), V(iSub),ctrlVol(iSub),
                                                     ngradLS(iSub),(*LSS)(iSub),(*(fluidSelector.fluidId))(iSub),
                                                     Nriemann,nsbar, fluidSelector,
                                                     A(iSub));
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int dimLS>
void Domain::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   DistGeoState &geoState,DistSVec<double,3>& X,DistSVec<double,dim> &V,
					   DistNodalGrad<dim>& ngrad,DistNodalGrad<dimLS> &ngradLS,
					   DistEdgeGrad<dim>* egrad,
					   DistVec<double> &ctrlVol,DistSVec<double,dimLS>& Phi,
					   DistMat<Scalar,dimLS> &A,DistLevelSetStructure* distLSS)
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

//	std::cout << "Num Subdomains = " << numLocSub << std::endl;

  if(inletRhsPat){
    fprintf(stdout, "with inletRhsPat\n");
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
      LevelSetStructure* lss = (distLSS) ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS, 
							 geoState(iSub),
							 X(iSub),V(iSub),ngrad(iSub),
							 ngradLS(iSub),
							 legrad,
							 ctrlVol(iSub),
							 Phi(iSub), 
							 A(iSub),lss,
							 inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addLSFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
      LevelSetStructure* lss = (distLSS) ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS, 
							 geoState(iSub),
							 X(iSub),V(iSub),ngrad(iSub),
							 ngradLS(iSub),
							 legrad,
							 ctrlVol(iSub),
							 Phi(iSub), 
							 A(iSub),lss,
							 inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addLSFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);
  }
}

template<int dim>
void Domain::recomputeRHS(VarFcn* vf, DistSVec<double,dim> &V, DistSVec<double,dim> &rhs,
                          DistExtrapolation<dim>* xpol, DistBcData<dim>& bcData,
                          DistGeoState& geoState, DistSVec<double,3> &X)
{
  int iSub;

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->recomputeRHS(vf, V(iSub), rhs(iSub), lxpol,
   				 bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletRhsData(*inletRhsPat, rhs.subData(iSub));
  }

  inletRhsPat->exchange();

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++)
    subDomain[iSub]->addRcvInletRhsData(*inletRhsPat, rhs.subData(iSub));


}
//------------------------------------------------------------------------------

template<int dim>
void Domain::recomputeRHS(VarFcn* vf, DistSVec<double,dim> &V, DistVec<int> &fluidId,
                          DistSVec<double,dim> &rhs, DistExtrapolation<dim>* xpol,
                          DistBcData<dim>& bcData, DistGeoState& geoState, DistSVec<double,3> &X)
{
  int iSub;
#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->recomputeRHS(vf, V(iSub), fluidId(iSub), rhs(iSub), lxpol,
                                  bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletRhsData(*inletRhsPat, rhs.subData(iSub));
  }

  inletRhsPat->exchange();

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++)
    subDomain[iSub]->addRcvInletRhsData(*inletRhsPat, rhs.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
double Domain::recomputeResidual(DistSVec<double,dim> &F, DistSVec<double,dim> &Finlet)
{

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; iSub++)
    subDomain[iSub]->recomputeResidual(F(iSub), Finlet(iSub));

  return Finlet*Finlet;

}

//------------------------------------------------------------------------------

template<int dim>
double Domain::computeRealFluidResidual(DistSVec<double, dim> &F, DistSVec<double,dim> &Freal,
                                        DistLevelSetStructure &dlss)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computeRealFluidResidual(F(iSub), Freal(iSub), dlss(iSub));

  return Freal*Freal;
}



//------------------------------------------------------------------------------

template<class Scalar, int neq>
void Domain::finishJacobianGalerkinTerm(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A)  {

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = A.getDiagMatPat();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->finishJacobianGalerkinTerm(ctrlVol(iSub), A(iSub));
    subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
  }

  double t = timer->addFiniteVolumeJacTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));

  com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &R,
				 DistVec<GhostPoint<dim>*> *ghostPoints,DistLevelSetStructure *LSS)
{

  double t0 = timer->getTime();

//#pragma omp parallel for
  if(ghostPoints)
    {
      if(!LSS) 
	{
	  cout<<"LSS has to be provided in the case of a viscous simulation\n";
	  exit(1);
	}
      //Vec<GhostPoint<dim>*> *gp;
#pragma omp parallel for
      for (int iSub = 0; iSub < numLocSub; ++iSub)
	{
	  //gp     = ghostPoints->operator[](iSub);
	  subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					       X(iSub), V(iSub), R(iSub), ghostPoints->operator[](iSub),
                                               &(LSS->operator()(iSub)));
	}
    }
  else
    {
#pragma omp parallel for
      for (int iSub = 0; iSub < numLocSub; ++iSub)
	subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					     X(iSub), V(iSub), R(iSub));
    }
  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X, DistSVec<double,3> &dX,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &dV, double dMach, DistSVec<double,dim> &dR)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					 X(iSub), dX(iSub), V(iSub), dV(iSub), dMach, dR(iSub));

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeOnlyGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &R)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					 X(iSub), V(iSub), R(iSub));

    subDomain[iSub]->sndData(*vecPat, R.subData(iSub));
  }

  timer->addFiniteElementTermTime(t0);

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, R.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVolumicForceTerm(VolumicForceTerm *volForce, DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &V, DistSVec<double,dim> &R)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeVolumicForceTerm(volForce, ctrlVol(iSub), V(iSub), R(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfVolumicForceTerm(VolumicForceTerm *volForce, DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
                               DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &dR)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeDerivativeOfVolumicForceTerm(volForce, ctrlVol(iSub), dCtrlVol(iSub), V(iSub), dV(iSub), dR(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeGalerkinBarTerm(bool doInitialTasks,
                                    FemEquationTerm *fet, DistBcData<dim> &bcData,
                                    DistGeoState &geoState, DistSVec<double,3> &X,
                                    DistMacroCellSet *macroCells,
                                    DistSVec<double,dim> &VBar,
                                    DistSVec<double,1> &volRatio,
                                    DistSVec<double,dim> &RBar,
                                    int scopeDepth1, int scopeDepth2)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeGalerkinBar_Step1(fet, bcData(iSub), geoState(iSub), X(iSub), VBar(iSub), (*Sigma)(iSub));
    subDomain[iSub]->sndData(*vecPat, (*Sigma).subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*Sigma).subData(iSub));


#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub) {
    MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];
    for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

    subDomain[iSub]->computeGalerkinBar_Step2(macCells, volRatio(iSub), (*Sigma)(iSub), RBar(iSub), scopeDepth1);
    delete [] macCells;
  }

  delete (Sigma);

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computePointWiseSourceTerm(DistGeoState &geoState, DistVec<double> &ctrlVol,
					DistNodalGrad<dim> &ngrad, DistSVec<double,dim> &VV,
					DistSVec<double,dim> &RR)
{

  const double sixth = 1.0/6.;
  const double cb1 = 0.1355;
  const double cb2 = 0.622;
  const double cw2 = 0.3;
  const double cw3 = 2.0;
  const double cv1 = 7.1;
  const double cv2 = 5.0;
  const double sigma = 2.0/3.0;
  const double vkcst = 0.41;
  const double reynolds = 2.91e6;

  const double cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  const double opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  const double cv1_pow3 = cv1*cv1*cv1;
  const double oocv2 = 1.0 / cv2;
  double oosigma = 1.0 / sigma;
  const double oovkcst2 = 1.0 / (vkcst*vkcst);
  double cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;
  const double ooreynolds = 1.0/ reynolds;

  cw1 /= reynolds;
  oosigma /= reynolds;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    Vec<double>& d2w = geoState(iSub).getDistanceToWall();
    double* cvol = ctrlVol.subData(iSub);
    SVec<double,dim>& dVdx = ngrad(iSub).getX();
    SVec<double,dim>& dVdy = ngrad(iSub).getY();
    SVec<double,dim>& dVdz = ngrad(iSub).getZ();
    double (*V)[dim] = VV.subData(iSub);
    double (*R)[dim] = RR.subData(iSub);
    for (int i=0; i<RR.subSize(iSub); ++i) {
      if (d2w[i] < 1.e-10) continue;
      double mutilde = V[i][0]*V[i][5];
      double mul = 1.0;
      double chi = max(mutilde/mul, 0.001);
      double chi3 = chi*chi*chi;
      double fv1 = chi3 / (chi3 + cv1_pow3);
      double fv2 = 1.0 + oocv2*chi;
      fv2 = 1.0 / (fv2*fv2*fv2);
      double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
      double d2wall = d2w[i];
      double ood2wall2 = 1.0 / (d2wall * d2wall);
      double rho = V[i][0];
      double oorho = 1.0 / rho;
      double zz = ooreynolds * oovkcst2 * mutilde * oorho * ood2wall2;
      double s12 = dVdy[i][1] - dVdx[i][2];
      double s23 = dVdz[i][2] - dVdy[i][3];
      double s31 = dVdx[i][3] - dVdz[i][1];
      double s = sqrt(s12*s12 + s23*s23 + s31*s31);
      double Stilde = s*fv3 + zz*fv2;
      double rr = min(zz/Stilde, 2.0);
      double rr2 = rr*rr;
      double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
      double gg2 = gg*gg;
      double fw = opcw3_pow * gg * pow(1.0/(gg2*gg2*gg2 + cw3_pow6), sixth);

      double S = cb1 * Stilde * mutilde;
      //S += oosigma * cb2 * rho * (dVdx[i][5]*dVdx[i][5] + dVdy[i][5]*dVdy[i][5] + dVdz[i][5]*dVdz[i][5]);
      S -= cw1 * fw * oorho * mutilde*mutilde * ood2wall2;
      R[i][5] -= cvol[i] * S;
    }
  }

}

//------------------------------------------------------------------------------
//----------------- All LES Models start  here

template<int dim>
void Domain::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, DistSVec<double,3> &X,
				       DistSVec<double,dim> &V, DistSVec<double,dim> &R)

{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeSmagorinskyLESTerm(smag, X(iSub), V(iSub), R(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDynamicLESTerm(DynamicLESTerm *dles, DistSVec<double,2> &Cs,
                                   DistSVec<double,3> &X, DistSVec<double,dim> &V,
				   DistSVec<double,dim> &R)
{

#pragma omp parallel for
   for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->computeDynamicLESTerm(dles, Cs(iSub), X(iSub), V(iSub), R(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVMSLESTerm(VMSLESTerm *vmst, DistMacroCellSet *macroCells,
                               bool doInitialTasks, DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistSVec<double,dim> &R, int scopeDepth)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;

  // Compute the large-scale component of V (VBar) and //
  // the volume ratios (Vi/VIi) with the macro-cells   //

  macroCells->computeVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth);

  // Exchange VBar and volRatio values //

  assemble(vecPat, VBar);
  assemble(vecPat, volRatio);

  // Compute the contribution to R from the tetrahedra contained within each subdomain //
#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeVMSLES_Step1(vmst, VBar(iSub),X(iSub), V(iSub), (*Sigma)(iSub));

  assemble(vecPat,*Sigma);

  // Add the subgrid scale viscosity to the residual vector R //

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub) {
    MacroCellSet* macCells;
    macCells = macroCells->obtainMacroCell(iSub, scopeDepth-1);
    subDomain[iSub]->computeVMSLES_Step2(volRatio(iSub),macCells,(*Sigma)(iSub), R(iSub), scopeDepth);
  }

  delete (Sigma);

  double t = timer->addVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDynamicVMSTerm(DynamicVMSTerm *dvmst, DistMacroCellSet *macroCells,
                                   bool doInitialTasks, DistVec<double> &ctrlVol, DistSVec<double,dim> **VBar,
                                   DistSVec<double,1> **volRatio, DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                   DistSVec<double,dim> &S, DistSVec<double,dim> &R, DistSVec<double,dim> &RBar,
                                   DistSVec<double,dim> &dWdt, DistSVec<double,dim> &dWBardt,
                                   DistSVec<double,dim> &M, DistSVec<double,dim> &MBar,
                                   int scopeDepth1, int scopeDepth2, int method, DistVec<double> *Cs)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  DistSVec<double,dim>* SigmaBar = new DistSVec<double,dim>(getNodeDistInfo());

  *Sigma = 0.0;
  *SigmaBar = 0.0;

  // compute volume averaged filter width for each node (used in post processing) //

  if (doInitialTasks) {
    Delta = new DistVec<double>(getNodeDistInfo());
    CsDelSq = new DistVec<double>(getNodeDistInfo());
    PrT = new DistVec<double>(getNodeDistInfo());
    WCsDelSq = new DistVec<double>(getNodeDistInfo());
    WPrT = new DistVec<double>(getNodeDistInfo());

   *Delta = 0.0;

    #pragma omp parallel for
      for (int iSub=0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->computeFilterWidth(X(iSub),(*Delta)(iSub));

    applySmoothing(ctrlVol, *Delta);
  }


  *WCsDelSq = 0.0; *WPrT = 0.0;  // variables that store the smoothed CsDelSq and PrT values


  // compute MBar and M //

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,dim> **vBar = new SVec<double,dim>*[2];
      SVec<double,1> **volR = new SVec<double,1>*[2];

      for (int i = 0; i<2; ++i) {
        vBar[i]  =  &(*VBar[i])(iSub);
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      subDomain[iSub]->computeMBarAndM_Step1(dvmst, vBar, volR, X(iSub), V(iSub), (*SigmaBar)(iSub), (*Sigma)(iSub));

      delete [] vBar;
      delete [] volR;
    }


    assemble(vecPat, (*SigmaBar));
    assemble(vecPat, (*Sigma));

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];

      for (int i = 0; i<2; ++i) {
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      subDomain[iSub]->computeMBarAndM_Step2(macCells, volR, MBar(iSub), M(iSub), (*SigmaBar)(iSub),
                                             (*Sigma)(iSub), scopeDepth1, scopeDepth2);

      delete [] volR;
      delete [] macCells;
    }


    assemble(vecPat, MBar);
    assemble(vecPat, M);

  // compute (Cs*Delta)^2 and PrT //

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeCsDeltaSq(R(iSub), RBar(iSub), M(iSub), MBar(iSub),
                                        dWdt(iSub), dWBardt(iSub), (*CsDelSq)(iSub),
                                        (*PrT)(iSub), method);
   }

  // local smoothing of (Cs*Delta)^2 and PrT //

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub) {
       subDomain[iSub]->computeLocalAvg(X(iSub), (*CsDelSq)(iSub), (*WCsDelSq)(iSub));
       subDomain[iSub]->computeLocalAvg(X(iSub), (*PrT)(iSub), (*WPrT)(iSub));
   }

  applySmoothing(ctrlVol, *WCsDelSq);
  applySmoothing(ctrlVol, *WPrT);


  *Sigma = 0.0;

  // computing the reynold stress flux "S" //

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,dim> **vBar = new SVec<double,dim>*[2];
      Vec<double>* cs;

      for (int i = 0; i<2; ++i) {
        vBar[i]  =  &(*VBar[i])(iSub);
      }

      if (Cs) cs = &(*Cs)(iSub);
      else cs = 0;

      subDomain[iSub]->computeDynamicVMSTerm_Step1(dvmst, vBar, X(iSub), V(iSub), (*Sigma)(iSub),
                                                   (*WCsDelSq)(iSub), (*WPrT)(iSub), cs, (*Delta)(iSub));

      delete [] vBar;
    }
    assemble(vecPat, *Sigma);

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];
      Vec<double>* cs;

      for (int i = 0; i<2; ++i) {
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      if (Cs) cs = &(*Cs)(iSub);
      else cs = 0;

      subDomain[iSub]->computeDynamicVMSTerm_Step2(macCells, volR, (*Sigma)(iSub), S(iSub), scopeDepth1);

      delete [] volR;
      delete [] macCells;
    }

    assemble(vecPat, S);

  delete (Sigma);
  delete (SigmaBar);

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWaleLESTerm(WaleLESTerm *wale, DistSVec<double,3> &X,
				DistSVec<double,dim> &V, DistSVec<double,dim> &R)

{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeWaleLESTerm(wale, X(iSub), V(iSub), R(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVBar(DistMacroCellSet *macroCells, bool doInitialTasks, DistGeoState &geoState,
                         DistSVec<double,dim> &VBar, DistSVec<double,dim> &V, int scopeDepth, int n)

{

  double t0 = timer->getTime();

  macroCells->computeVBar(doInitialTasks, geoState, V, VBar, scopeDepth, n);

  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, VBar.subData(iSub));

    vecPat->exchange();

  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvData(*vecPat, VBar.subData(iSub));

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeBarTerm(DistMacroCellSet *macroCells, bool doInitialTasks, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> **VBar, DistSVec<double,1> **volRatio,
                            DistSVec<double,3> &X, DistSVec<double,dim> &V, int scopeDepth1, int scopeDepth2)

{

  double t0 = timer->getTime();

  // Computing the large-scale components VBar and volRatio //

  macroCells->computeDVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth1, scopeDepth2);

   // Exchange of information between subdomains for VBar and volRatio //

   for (int i =0; i < 2; ++i) {

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*VBar[i]).subData(iSub));

    vecPat->exchange();

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvData(*vecPat, (*VBar[i]).subData(iSub));

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->sndData(*vecPat, (*volRatio[i]).subData(iSub));

    vecPat->exchange();

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->addRcvData(*vecPat, (*volRatio[i]).subData(iSub));

   }

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTestFilterValues(DistVec<double> &ctrlVol,
                                     DistSVec<double,dim> &VCap,
                                     DistSVec<double,16> &Mom_Test,
                                     DistSVec<double,6> &Sij_Test,
                                     DistVec<double> &modS_Test,
                                     DistSVec<double,8> &Eng_Test,
                                     DistSVec<double,2> &Cs,
				     DistVec<int> &Ni,
                                     DistBcData<dim> &bcData,
				     DistSVec<double,3> &X, DistSVec<double,dim> &V,
				     double gam, double R)
{

// computing test filtered values for all the required flow variables //

#pragma omp parallel for
   for (int iSub = 0; iSub < numLocSub; ++iSub){
     subDomain[iSub]->computeTestFilterAvgs(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub), modS_Test(iSub),
                                            Eng_Test(iSub), X(iSub), V(iSub), gam, R);
     subDomain[iSub]->sndData(*vecPat, VCap.subData(iSub));
   }

   // START OF ALL EXCHANGES

  vecPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vecPat, VCap.subData(iSub));
    subDomain[iSub]->sndData(*momPat, Mom_Test.subData(iSub));
  }

  momPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*momPat, Mom_Test.subData(iSub));
    subDomain[iSub]->sndData(*engPat, Eng_Test.subData(iSub));
  }

  engPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*engPat, Eng_Test.subData(iSub));
    subDomain[iSub]->sndData(*weightPat, Sij_Test.subData(iSub));
  }

  weightPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*weightPat, Sij_Test.subData(iSub));
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(modS_Test.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(modS_Test.subData(iSub)));

   // END OF ALL EXCHANGES


// computing unsmoothed Cs and Pt values //

#pragma omp parallel for
   for(int iSub = 0; iSub < numLocSub; ++iSub) {
     subDomain[iSub]->computeCsValues(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub),
                                      modS_Test(iSub), Eng_Test(iSub), Cs(iSub),
				      Ni(iSub), X(iSub), gam, R);
   }


// smoothing Cs and Pt values //

   DistSVec<double,2>* W = new DistSVec<double,2>(getNodeDistInfo());
   *W = 0.0;

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeLocalAvg(X(iSub), Cs(iSub), (*W)(iSub));

   applySmoothing(ctrlVol, *W);
   Cs = *W;

   delete (W);

}

//------------------------------------------------------------------------------
// Included (MB)
template<int dim>
void Domain::computeDerivativeOfSmagorinskyLESTerm(SmagorinskyLESTerm *smag, DistSVec<double,3> &X,
				       DistSVec<double,dim> &V, DistSVec<double,dim> &R)

{

  com->fprintf(stderr, "***** Domain::computeDerivativeOfSmagorinskyLESTerm is not implemented!\n");
  exit(1);

}

//----------------------- All LES Models End Here
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------Start of routines that compute MutOMu values

template<int dim>
void Domain::computeMutOMuSmag(SmagorinskyLESTerm *smag, DistVec<double> &ctrlVol,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuSmag(smag, X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuVMS(VMSLESTerm *vms, DistMacroCellSet *macroCells, DistVec<double> &ctrlVol,
                              bool doInitialTasks, DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                              DistSVec<double,3> &X, DistSVec<double,dim> &V, int scopeDepth,
                              DistVec<double> &mutOmu)
{

  macroCells->computeVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth);
  assemble(vecPat, VBar);

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuVMS(vms, VBar(iSub),X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuDynamicVMS(DynamicVMSTerm *dvms, DistVec<double> &ctrlVol,
                                     DistSVec<double,dim> &VBar, DistSVec<double,3> &X,
                                     DistSVec<double,dim> &V, DistVec<double> &Cs, DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuDynamicVMS(dvms, VBar(iSub),X(iSub), V(iSub), Cs(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuWale(WaleLESTerm *wale, DistVec<double> &ctrlVol,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuWale(wale, X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuDynamicLES(DynamicLESTerm *dles, DistVec<double> &ctrlVol,
                                     DistSVec<double,2> &Cs, DistSVec<double,3> &X,
                                     DistSVec<double,dim> &V, DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuDynamicLES(dles, Cs(iSub), X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}


//--------End of routines that compute MutOMu values
//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void Domain::computeJacobianGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
					 DistGeoState &geoState, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
					 DistMat<Scalar,neq> &A,
                                         DistVec<GhostPoint<dim>*> *ghostPoints,DistLevelSetStructure *distLSS)
{

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = A.getOffDiagMatPat();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    Vec<GhostPoint<dim>*>* gp = (ghostPoints? &(*ghostPoints)(iSub) :0);
    LevelSetStructure *LSS = distLSS ? &(distLSS->operator()(iSub)) : 0;
    subDomain[iSub]->computeJacobianGalerkinTerm(fet, bcData(iSub), geoState(iSub), X(iSub),
						 ctrlVol(iSub), V(iSub), A(iSub),gp,LSS);
    subDomain[iSub]->sndOffDiagBlocks(*matPat, A(iSub));
  }

  double t = timer->addFiniteElementJacTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvOffDiagBlocks(*matPat, A(iSub));
  
  if (ghostPoints) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->sndGhostOffDiagBlocks(*matPat, A(iSub));
    }
    
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvGhostOffDiagBlocks(*matPat, A(iSub));
  }

  com->printf(6, "FE Jacobian matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeBCsJacobianWallValues(FemEquationTerm *fet, DistBcData<dim> &bcData,
					  DistGeoState &geoState, DistSVec<double,3> &X,
					  DistSVec<double,dim> &V)
{

  int iSub;

  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeBCsJacobianWallValues(fet, bcData(iSub), geoState(iSub), X(iSub), V(iSub));
  }

  double t = timer->addFiniteElementJacTime(t0);

  com->printf(6, "FE wall BC computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeJacobianVolumicForceTerm(VolumicForceTerm *volForce,
                               DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &V, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeJacobianVolumicForceTerm(volForce,
                                           ctrlVol(iSub), V(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::getExtrapolationValue(DistExtrapolation<dim> *xpol, DistSVec<double,dim> &V,
                                      DistSVec<double,dim> &Ubc, VarFcn* vf,
                                 DistBcData<dim>& bcData, DistGeoState& geoState, DistSVec<double,3>& X)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->getExtrapolationValue(lxpol, V(iSub), Ubc(iSub), vf, bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletData(*inletRhsPat, Ubc.subData(iSub));
  }
 inletRhsPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
   subDomain[iSub]->addRcvInletData(*inletRhsPat, Ubc.subData(iSub), true);
  }


}
//------------------------------------------------------------------------------

template<int dim>
void Domain::applyExtrapolationToSolutionVector(DistExtrapolation<dim> *xpol, DistSVec<double,dim> &U,
                                      DistSVec<double,dim> &Ubc)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->applyExtrapolationToSolutionVector(lxpol, U(iSub), Ubc(iSub));
  }

}
//------------------------------------------------------------------------------

template<int dim>
void Domain::applyBCsToSolutionVector(BcFcn *bcFcn, DistBcData<dim> &bcData,
                                      DistSVec<double,dim> &U, DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    { 
      LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->applyBCsToSolutionVector(bcFcn, bcData(iSub), U(iSub), LSS);
    }
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::applyBCsToResidual(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistSVec<double,dim> &F, DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    {
      LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->applyBCsToResidual(bcFcn, bcData(iSub), U(iSub), F(iSub), LSS);
    }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::applyBCsToDerivativeOfResidual(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,dim> &dF)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToDerivativeOfResidual(bcFcn, bcData(iSub), U(iSub), dU(iSub), dF(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::applyBCsToJacobian(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToJacobian(bcFcn, bcData(iSub), U(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::applyBCsToH2Jacobian(BcFcn *bcFcn, DistBcData<dim> &bcData,
	                          DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToH2Jacobian(bcFcn, bcData(iSub), U(iSub), A(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void Domain::applyBCsToJacobianWallValues(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToJacobianWallValues(bcFcn, bcData(iSub), U(iSub), A(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar2>
void Domain::applyBCsToProduct(BcFcn *bcFcn, DistBcData<dim> &bcData, DistSVec<double,dim> &U, DistSVec<Scalar2,dim> &Prod)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToProduct(bcFcn, bcData(iSub), U(iSub), Prod(iSub));

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::computeH1(FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                       DistGeoState &geoState, DistVec<double> &ctrlVol,
                       DistSVec<double,dim> &V, DistMat<Scalar,dim> &H1)
{

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = H1.getDiagMatPat();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeH1(fluxFcn, bcData(iSub), geoState(iSub),
                               ctrlVol(iSub), V(iSub), H1(iSub));
    subDomain[iSub]->sndDiagBlocks(*matPat, H1(iSub));
  }

  //double t = timer->addH1SetupTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvDiagBlocks(*matPat, H1(iSub));

  //com->printf("Time for computing H1 matrix: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn,
		       DistBcData<dim> &bcData, DistGeoState &geoState,
		       DistSVec<double,3> &X, DistSVec<double,dim> &V,
		       DistNodalGrad<dim, double> &ngrad, DistMat<Scalar,neq> &H2,
		       DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
		       DistSVec<double,dim> &bij, DistSVec<double,dim> &bji)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeH2(fluxFcn, recFcn, bcData(iSub), geoState(iSub),
			       X(iSub), V(iSub), ngrad(iSub), H2(iSub));
    subDomain[iSub]->precomputeRec(recFcn, X(iSub), V(iSub), ngrad(iSub),
				   aij(iSub), aji(iSub), bij(iSub), bji(iSub));
  }

  double t = timer->addH2SetupTime(t0);

  com->printf(6, "H2 matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2(RecFcn *recFcn, DistSVec<double,3> &X,
	     DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
	     DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
	     DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
	     DistSVec<Scalar2,dim> &p, DistNodalGrad<dim, Scalar2> &dpdxj,
	     DistSVec<Scalar2,dim> &prod)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2(recFcn, X(iSub), ctrlVol(iSub), H2(iSub),
		     aij(iSub), aji(iSub), bij(iSub), bji(iSub),
		     p(iSub), dpdxj(iSub), prod(iSub));
    subDomain[iSub]->sndData(*vPat, prod.subData(iSub));
  }

  vPat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));

}
//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2T(RecFcn *recFcn, DistSVec<double,3> &X,
                DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
                DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
                DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod,
                DistSVec<Scalar2,dim> &prod2, DistSVec<Scalar2,dim> &prod3,
                DistSVec<Scalar2,dim> &prod4)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2T(recFcn, X(iSub), ctrlVol(iSub),
        H2(iSub), aij(iSub), aji(iSub), bij(iSub), bji(iSub), p(iSub),
        prod(iSub), prod2(iSub), prod3(iSub), prod4(iSub));

    subDomain[iSub]->sndData(*vPat, prod.subData(iSub));
  }

  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));

 assemble(vPat, prod2);
 assemble(vPat, prod3);
 assemble(vPat, prod4);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2Tb(RecFcn *recFcn, DistSVec<double,3> &X,
                DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                DistNodalGrad<dim, Scalar2> &dpdxj, DistSVec<Scalar2,dim> &p,
                DistSVec<Scalar2,dim> &prod, DistSVec<Scalar2,dim> &prod2)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2Tb(recFcn, X(iSub), ctrlVol(iSub),
        H2(iSub), dpdxj(iSub), p(iSub), prod(iSub), prod2(iSub) );
  }
// No assemble????

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::assemble(CommPattern<Scalar> *commPat, DistSVec<Scalar,dim> &W)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*commPat, W.subData(iSub));

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*commPat, W.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, class OpType >
void Domain::assemble(CommPattern<Scalar> *commPat, DistSVec<Scalar,dim> &W, const OpType &oper)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*commPat, W.subData(iSub));

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->operateRcvData(*commPat, W.subData(iSub), oper);

}

//------------------------------------------------------------------------------

template<class Scalar>
void Domain::assemble(CommPattern<Scalar> *commPat, DistVec<Scalar> &W)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

  commPat->exchange();
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->addRcvData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

}
//------------------------------------------------------------------------------

template<class Scalar, class OpType >
void Domain::assemble(CommPattern<Scalar> *commPat, DistVec<Scalar> &W, const OpType& oper)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->operateRcvData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)),oper);
    }

}

//------------------------------------------------------------------------------
template<int dim>
void Domain::assembleGhostPoints(DistVec<GhostPoint<dim>*> &ghostPoints)
{
  int iSub;
  // Adam 2010.10.27
  // Caution, the order of the calls matters, because a ghost point can lie on a domain boundary, 
  // in which case we may want to create its state after during the exchange. The Ghost weight is 
  // going to be used as a parameter.
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndGhostWeights(*volPat, ghostPoints(iSub));
      subDomain[iSub]->sndGhostStates(*vecPat, ghostPoints(iSub));
    }

  volPat->exchange();
  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->rcvGhostWeights(*volPat, ghostPoints(iSub));
      subDomain[iSub]->rcvGhostStates(*vecPat, ghostPoints(iSub));
    }
}
//------------------------------------------------------------------------------

template<class Scalar, int dim>
bool Domain::readVectorFromFile(const char *prefix, int step, double *tag,
				DistSVec<Scalar,dim> &U, Scalar* scale)
{

  int neq, numSteps;
  double t = subDomain[0]->template readTagFromFile<Scalar,dim>(prefix, step, &neq, &numSteps);
  if (tag) *tag = t;

  if (neq != dim)
    com->printf(1, "*** Warning: mismatch in dim for \'%s\' (%d vs %d)\n", prefix, neq, dim);

  if (step >= numSteps)
    return false;

  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->readVectorFromFile(prefix, step, neq, U(iSub), scale);

  timer->addBinaryReadTime(t0);

  com->printf(3, "Read solution %d from \'%s\'\n", step, prefix);

  return true;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::writeVectorToFile(const char *prefix, int step, double tag,
			       DistSVec<Scalar,dim> &U, Scalar* scale)
{

  int iSub;

  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template openFileForWriting<Scalar,dim>(prefix, step);

  if (step == 0)
    com->barrier();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->writeVectorToFile(prefix, step, U(iSub), scale);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template writeTagToFile<Scalar,dim>(prefix, step, tag);

#ifndef SYNCHRO_WRITE
  sync();
#endif

  timer->addBinaryWriteTime(t0);

  com->printf(4, "Wrote solution %d to \'%s\'\n", step, prefix);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::scaleSolution(DistSVec<Scalar,dim> &data, RefVal* refVal)  {

  int iSub;

  double scale[dim];

  scale[0] = refVal->density;
  scale[1] = refVal->density*refVal->velocity;
  scale[2] = refVal->density*refVal->velocity;
  scale[3] = refVal->density*refVal->velocity;
  scale[4] = refVal->energy;
  
  if (dim == 6) {
    scale[5] = refVal->density*refVal->nutilde;
  }

  if (dim == 7) {
    scale[5] = refVal->density*refVal->kenergy;
    scale[6] = refVal->density*refVal->epsilon;
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)  {

    SVec<Scalar, dim> &subData = data(iSub);

    for (int i = 0; i < data.size(); ++i)
      for (int j = 0; i < dim; ++i)
        subData[i][j] *= scale[j];
  }

}

//------------------------------------------------------------------------------

template<class S1, class S2>
void Domain::computeStiffAndForce(DefoMeshMotionData::Element type, DistSVec<double,3>& X,
                                  DistSVec<double,3>& F, DistMat<S1,3>& K,
                                  DistMat<S2,3>* P, double volStiff, int** ndType=0)  {

  double t0 = timer->getTime();

  CommPattern<S2>* matPat = 0;
  if (P) matPat = P->getDiagMatPat();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    GenMat<S2,3>* p = (P) ? &((*P)(iSub)) : 0;
    int* subNdType = (ndType) ? ndType[iSub] : 0;
    subDomain[iSub]->computeStiffAndForce(type, X(iSub), F(iSub), K(iSub), p, volStiff, subNdType);
    subDomain[iSub]->sndData(*vec3DPat, F.subData(iSub));
    if (P) subDomain[iSub]->sndDiagBlocks(*matPat, (*P)(iSub));
  }

  double t = timer->addMeshAssemblyTime(t0);

  vec3DPat->exchange();
  if (P) matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec3DPat, F.subData(iSub));
    if (P) subDomain[iSub]->addRcvDiagBlocks(*matPat, (*P)(iSub));
  }

  F *= -1.0;

  com->printf(6, "K matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistSVec<double,dim> &U)
{

  int ierr = 0;

#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub));

  com->globalSum(1, &ierr);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistSVec<double,dim> &U, DistVec<int> &fluidId)
{

  int ierr = 0;

#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub), fluidId(iSub));

  com->globalSum(1, &ierr);

  return ierr;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::fixSolution(VarFcn *varFcn, DistSVec<double,dim> &U, DistSVec<double,dim> &dU,DistVec<int>* fluidId)
{

  int verboseFlag = com->getMaxVerbose();

//#pragma omp parallel for reduction(+: ierr)
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    if (fluidId) {
      subDomain[iSub]->fixSolution(varFcn, U(iSub), dU(iSub), &((*fluidId)(iSub)), verboseFlag);
    } else {
  
      subDomain[iSub]->fixSolution(varFcn,U(iSub),dU(iSub),NULL,verboseFlag);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistVec<double> &ctrlVol,
                          DistSVec<double,dim> &U, DistVec<int> &fluidId,
                          DistVec<int> &fluidIdn)
{

  int ierr = 0;

#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, ctrlVol(iSub), U(iSub), fluidId(iSub), fluidIdn(iSub));
  com->globalSum(1, &ierr);
  return ierr;

}
//------------------------------------------------------------------------------

template<int dim, int neq>
int Domain::clipSolution(TsData::Clipping ctype, BcsWallData::Integration wtype,
			 VarFcn* varFcn, double* Uin, DistSVec<double,dim>& U)
{

  const DistInfo& distInfo = U.info();

  int size = neq * distInfo.numGlobSub;
  int sizeint = 2 * size;
  int sizedouble = size;

  int* allint = reinterpret_cast<int *>(alloca(sizeof(int) * (sizeint + 1)));
  double* alldouble = reinterpret_cast<double *>(alloca(sizeof(double) * sizedouble));

  int k;
  for (k=0; k<sizeint; ++k)
    allint[k] = 0;
  for (k=0; k<sizedouble; ++k)
    alldouble[k] = 0.0;

  int (*allcmin)[neq] = reinterpret_cast<int (*)[neq]>(allint);
  int (*allpmin)[neq] = reinterpret_cast<int (*)[neq]>(allint + 1 * size);
  double (*allvmin)[neq] = reinterpret_cast<double (*)[neq]>(alldouble);

  int iSub;
  int ierr = 0;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    int gSub = distInfo.locSubToGlobSub[iSub];
    ierr += subDomain[iSub]->template
      clipSolution<dim,neq>(ctype, wtype, varFcn, Uin, U.getMasterFlag(iSub), U(iSub),
			    allcmin[gSub], allpmin[gSub], allvmin[gSub]);
  }

  allint[sizeint] = ierr;
  com->globalSum(sizeint + 1, allint);
  ierr = allint[sizeint];
  com->globalSum(sizedouble, alldouble);

  int cmin[neq];
  int pmin[neq];
  double vmin[neq];

  for (k=0; k<neq; ++k) {
    cmin[k] = 0;
    pmin[k] = allpmin[0][k];
    vmin[k] = allvmin[0][k];
  }

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
    for (k=0; k<neq; ++k) {
      cmin[k] += allcmin[iSub][k];
      if (allvmin[iSub][k] < vmin[k]) {
	pmin[k] = allpmin[iSub][k];
	vmin[k] = allvmin[iSub][k];
      }
    }
  }

  for (k=0; k<neq; ++k) {
    if (cmin[k] > 0) {
      if (ctype == TsData::NONE)
	com->printf(1, "*** Warning: %d negative %s value%s (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
      else if (ctype == TsData::ABS_VALUE)
	com->printf(1, "*** Warning: %d %s value%s clipped with abs (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
      else if (ctype == TsData::FREESTREAM)
	com->printf(1, "*** Warning: %d %s value%s clipped at freestream (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
    }
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkFailSafe(VarFcn* varFcn, DistSVec<double,dim>& U,
               DistSVec<bool,2>& tag, DistVec<int> *fluidId)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    if(!fluidId)
      subDomain[iSub]->checkFailSafe(varFcn, U(iSub), tag(iSub));
    else
      subDomain[iSub]->checkFailSafe(varFcn, U(iSub), tag(iSub), &(*fluidId)(iSub));
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkGradients(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
			    DistSVec<double,dim> &V, DistNodalGrad<dim> &ngrad)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkGradientsSetUp(X(iSub), V(iSub));

  ngrad.compute(0, ctrlVol, V);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkGradientsWrite(X(iSub), ngrad(iSub));

  com->barrier();

  exit(1);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkMatVecProd(DistSVec<double,dim> &prod, const char *msg)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkMatVecProd(prod(iSub), msg);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeForceDerivs(VarFcn *varFcn, DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                DistSVec<double,dim> &deltaU, Vec<double> &modalF,
                                VecSet< DistSVec<double,3> > &mX)  {

  modalF = 0.0;

  int nStrModes = modalF.len;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    Vec<double> locModalF(nStrModes);

    SVec<double, 3> **locModes = new SVec<double, 3> *[nStrModes];
    for (int iMode = 0; iMode < nStrModes; iMode++)
      locModes[iMode] = &mX[iMode](iSub);

    subDomain[iSub]->computeForceDerivs(varFcn, X(iSub), V(iSub), deltaU(iSub), locModalF, locModes);

    int iMode;
    for(iMode = 0; iMode < nStrModes; ++iMode) {
      #pragma omp critical
      modalF[iMode] += locModalF[iMode];
    }
  }

  com->globalSum(nStrModes, modalF.v);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::zeroInternalVals(DistSVec<double, dim> &v)  {

  int iSub;
//#pragma omp parallel for reduction(+: ierr)
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->zeroInternalVals(v(iSub));
}

//------------------------------------------------------------------------------
template<int dim>
void Domain::printVariable(DistSVec<double,dim>&V, VarFcn *vf)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printVariable(V(iSub), vf);
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::printInletVariable(DistSVec<double,dim>&V)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printInletVariable(V(iSub));
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::printAllVariable(DistVec<int> &X, DistSVec<double,dim>&V, int it)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printAllVariable(X(iSub), V(iSub), numLocSub, it);
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::checkExtrapolationValue(DistSVec<double,dim>&U, VarFcn* vf,
                                     DistBcData<dim>& bcData, DistGeoState& geoState)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->checkExtrapolationValue(U(iSub),  vf, bcData(iSub), geoState(iSub));
  com->barrier();
}

//------------------------------------------------------------------------------
template<class Scalar, int neq>
void Domain::printAllMatrix(DistMat<Scalar, neq> &A, int it)
{

  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printAllMatrix(A(iSub), it);
  com->barrier();
  exit(1);
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::assemble_dWdt(DistSVec<double, dim> &dWdt, DistSVec<double, dim> &Sigma)

{

  assemble(vecPat, dWdt);
  assemble(vecPat, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computedWBar_dt(DistSVec<double, dim> &dWBardt, DistSVec<double, dim> &Sigma,
                             DistMacroCellSet *macroCells, DistSVec<double,1> **volRatio,
                             int scopeDepth)
{

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth];

      for (int i = 0; i<2; ++i)
        volR[i]  =  &(*volRatio[i])(iSub);

      for (int i = 0; i<scopeDepth; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      subDomain[iSub]->computedWBar_dt(macCells, volR, Sigma(iSub), dWBardt(iSub), scopeDepth);

      delete [] volR;
      delete [] macCells;
    }

   assemble(vecPat, dWBardt);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->computeWeightsForEmbeddedStruct(V(iSub),VWeights(iSub),Weights(iSub),
                                                     (*distLSS)(iSub),X(iSub));
  
  assemble(vecPat, VWeights);
  assemble(volPat, Weights);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                             DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                                             DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights, 
                                             DistLevelSetStructure *distLSS, DistVec<int> *fluidId)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeWeightsForEmbeddedStruct(V(iSub),VWeights(iSub),Phi(iSub), PhiWeights(iSub),
                                                     Weights(iSub), (*distLSS)(iSub),X(iSub), (*fluidId)(iSub));

  assemble(vecPat, VWeights);
  assemble(phiVecPat, PhiWeights);
  assemble(volPat, Weights);
}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV)
{
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->extrapolatePhiV((*distLSS)(iSub), PhiV(iSub));

  assemble(phiVecPat, PhiV);
  // Note: here PhiV is not a distance function. It will be used only as an indicator of the
  //       fluidId on newly uncovered nodes.
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints, DistSVec<double,dim> &U, VarFcn *varFcn,DistLevelSetStructure *distLSS,DistVec<int> &tag)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->populateGhostPoints((*ghostPoints)(iSub),U(iSub),varFcn,(*distLSS)(iSub),tag(iSub));
    }

  // Adam 2010.10.26
  assembleGhostPoints(*ghostPoints);

  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->reduceGhostPoints((*ghostPoints)(iSub));
    }
}

template<int dim,int neq>
void Domain::populateGhostJacobian(DistVec<GhostPoint<dim>*> &ghostPoints,DistSVec<double,dim> &U,VarFcn *varFcn,DistLevelSetStructure &LSS,DistVec<int> &tag, DistMat<double,neq>& A) {

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->populateGhostJacobian(ghostPoints(iSub), U(iSub), varFcn, LSS(iSub), tag(iSub), A(iSub));
 
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeRiemannWeightsForEmbeddedStruct(V(iSub), Wstarij(iSub), Wstarji(iSub), 
                                                     VWeights(iSub),Weights(iSub), (*distLSS)(iSub),X(iSub));
  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
}

//-------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
               DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
               DistLevelSetStructure *distLSS, DistVec<int> *fluidId0, DistVec<int> *fluidId)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeRiemannWeightsForEmbeddedStruct(V(iSub), Wstarij(iSub), Wstarji(iSub),
                                                     VWeights(iSub),Weights(iSub), Phi(iSub), PhiWeights(iSub),
                                                     (*distLSS)(iSub),X(iSub), (*fluidId0)(iSub), (*fluidId)(iSub));
  assemble(vecPat, VWeights);
  assemble(phiVecPat, PhiWeights);
  assemble(volPat, Weights);
}

//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::checkNodePhaseChange(DistSVec<double,dimLS> &PhiProduct)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkNodePhaseChange(PhiProduct(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::storePreviousPrimitive(DistSVec<double,dim> &V, DistVec<int> &fluidId, 
                                    DistSVec<double,3> &X, DistSVec<double,dim> &Vupdate, 
                                    DistVec<double> &weight){
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->storePreviousPrimitive(V(iSub), fluidId(iSub), X(iSub), Vupdate(iSub), weight(iSub));

  assemble(vecPat, Vupdate);
  assemble(volPat, weight);

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::IncreasePressure(double p, VarFcn *vf,  DistSVec<double,dim> &U){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->IncreasePressure(p,vf,U(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::IncreasePressure(double p, VarFcn *vf,  DistSVec<double,dim> &U, DistVec<int> &fluidId){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->IncreasePressure(p,vf,U(iSub),fluidId(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::padeReconstruction(VecSet<DistSVec<double, dim> >&snapsCoarse, VecSet<DistSVec<double, dim> >&snaps, int *stepParam, double *freqCoarse, double deltaFreq, int nStrMode, int L, int M, int nPoints)
{
  int nSteps = stepParam[0];
  int nSnaps = stepParam[2];
  int nStepsCoarse = stepParam[1];
  int nSnapsCoarse = stepParam[3];


  int iSub;

  for (iSub = 0; iSub < numLocSub; ++iSub)  {

     SVec<double, dim> **locVecSetCoarse = new SVec<double, dim> *[nSnapsCoarse];
     SVec<double, dim> **locVecSet = new SVec<double, dim> *[nSnaps];
      for (int iSnapsCoarse = 0; iSnapsCoarse < nSnapsCoarse; iSnapsCoarse++){
        locVecSetCoarse[iSnapsCoarse] = &snapsCoarse[iSnapsCoarse](iSub);

      }
      for (int iSnaps = 0; iSnaps < nSnaps; iSnaps++) {

        locVecSet[iSnaps] = &snaps[iSnaps](iSub);


      }


      subDomain[iSub]->padeReconstruction(locVecSetCoarse, locVecSet, stepParam, freqCoarse, deltaFreq, nStrMode, L, M, nPoints);

  }

}

//------------------------------------------------------------------------------
template<int dim>
void Domain::hardyInterpolationLogMap(VecSet<DistSVec<double, dim> >**logMap, VecSet<DistSVec<double, dim> >&logMapInterp, int nData, int numPod, int iDataMin, FullM &B, FullM &b)
{

  SVec<double, dim> ***locVecSet = new SVec<double, dim> **[nData];
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {
    for (int iData=0; iData < nData; ++iData) {
      if (iData != iDataMin) {
        locVecSet[iData] = new SVec<double, dim> *[numPod];
        for (int iPod = 0; iPod < numPod; ++iPod)
          (locVecSet[iData])[iPod] = &((*(logMap[iData]))[iPod])(iSub);
      }
    }
    SVec<double, dim> **locLogMapInterp = new SVec<double, dim> *[numPod];
    for (int iPod = 0; iPod < numPod; ++iPod)
      locLogMapInterp[iPod] = &logMapInterp[iPod](iSub);
    subDomain[iSub]->hardyInterpolationLogMap(locVecSet,locLogMapInterp,nData,numPod,iDataMin,B,b);
    for (int iPod = 0; iPod < numPod; ++iPod)
      locLogMapInterp[iPod] = 0;
    delete[] locLogMapInterp;
  }
  for (int iData=0; iData < nData; ++iData) {
    if (iData != iDataMin){
      for (int iPod=0; iPod<numPod;++iPod)
        (locVecSet[iData])[iPod] = 0;
      delete [] locVecSet[iData];
    }
  }
  delete [] locVecSet;
}
//------------------------------------------------------------------------------
// Included (MB)
template<int dim>
void Domain::getGradP(DistNodalGrad<dim>& ngrad)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getGradP(ngrad(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::getDerivativeOfGradP(DistNodalGrad<dim>& ngrad)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getDerivativeOfGradP(ngrad(iSub));

}

//-----------------------------------------------------------------------------

template<int dim>
void Domain::computeCVBasedForceLoad(int forceApp, int orderOfAccuracy, DistGeoState& geoState,
                                     DistSVec<double,3> &X, double (*Fs)[3], int sizeFs,
                                     DistLevelSetStructure *distLSS,
                                     DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, double pInfty,
				     VarFcn* vf, DistVec<int>* fid)
{
  //PJSA double subFs[numLocSub][sizeFs][3];
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  DistVec<double> pstarij(Wstarij.info());
  DistVec<double> pstarji(Wstarji.info()); //extract p from Wstar
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    int (*ptr)[2] = subDomain[iSub]->getEdges().getPtr();
    for (int i=0; i<Wstarij(iSub).size(); i++) {
      pstarij(iSub)[i] = vf->getPressure(Wstarij(iSub)[i],fid?(*fid)(iSub)[ptr[i][0]]:0);
      pstarji(iSub)[i] = vf->getPressure(Wstarji(iSub)[i],fid?(*fid)(iSub)[ptr[i][1]]:0);
    }
  }
 
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    subDomain[iSub]->computeCVBasedForceLoad(forceApp, orderOfAccuracy, geoState(iSub), X(iSub), subFs[iSub],
                                             sizeFs, (*distLSS)(iSub), pstarij(iSub), pstarji(iSub), pInfty);
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-----------------------------------------------------------------------------

template<int dim>
void Domain::computeCVBasedForceLoadViscous(int forceApp, int orderOfAccuracy, DistGeoState& geoState,
					    DistSVec<double,3> &X, double (*Fs)[3], int sizeFs,
					    DistLevelSetStructure *distLSS,double pInfty,
					    DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
					    DistSVec<double,dim> &V, DistVec<GhostPoint<dim>*> *ghostPoints, 
					    PostFcn *postFcn, DistNodalGrad<dim, double> *ngrad)
{
  //PJSA double subFs[numLocSub][sizeFs][3];
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  Vec<GhostPoint<dim>*> *gp=0;
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    if(ghostPoints) gp = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeCVBasedForceLoadViscous(forceApp, orderOfAccuracy, geoState(iSub), X(iSub), subFs[iSub],
						    sizeFs, (*distLSS)(iSub), pInfty, Wstarij(iSub), Wstarji(iSub),
						    V(iSub),gp,postFcn,(*ngrad)(iSub));
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::computeRecSurfBasedForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X,
					  double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS,
					  DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, double pInfty,
					  VarFcn* vf, DistVec<int>* fid)
{
  //PJSA double subFs[numLocSub][sizeFs][3];
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  DistVec<double> pstarij(Wstarij.info());
  DistVec<double> pstarji(Wstarji.info()); //extract p from Wstar
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    int (*ptr )[2] = subDomain[iSub]->getEdges().getPtr();
    if (fid) {
      for (int i=0; i<Wstarij(iSub).size(); i++) {
	pstarij(iSub)[i] = vf->getPressure(Wstarij(iSub)[i],(*fid)(iSub)[ptr[i][0]]);//Wstarij(iSub)[i][4];
	pstarji(iSub)[i] = vf->getPressure(Wstarji(iSub)[i],(*fid)(iSub)[ptr[i][1]]);//Wstarji(iSub)[i][4];
      }
    } else {
      for (int i=0; i<Wstarij(iSub).size(); i++) {
	pstarij(iSub)[i] = vf->getPressure(Wstarij(iSub)[i],0);//Wstarij(iSub)[i][4];
	pstarji(iSub)[i] = vf->getPressure(Wstarji(iSub)[i],0);//Wstarji(iSub)[i][4];
      }

    }
  }

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    subDomain[iSub]->computeRecSurfBasedForceLoad(forceApp, orderOfAccuracy, X(iSub), subFs[iSub], sizeFs,
                                                  (*distLSS)(iSub), pstarij(iSub), pstarji(iSub), pInfty);
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::computeRecSurfBasedForceLoadViscous(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
						 double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS,
						 double pInfty, DistSVec<double,dim> &V, 
						 DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn)
{
  //PJSA double subFs[numLocSub][sizeFs][3];
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  Vec<GhostPoint<dim>*> *gp;
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    gp     = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeRecSurfBasedForceLoadViscous(forceApp, orderOfAccuracy, X(iSub), subFs[iSub], sizeFs,
							 (*distLSS)(iSub), pInfty, V(iSub), gp, postFcn);
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------
// This one should replace the two previous ones.
template<int dim>
void Domain::computeRecSurfBasedForceLoadNew(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
					     double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS, double pInfty, 
					     DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
					     DistSVec<double,dim> &V, 
					     DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn,DistVec<int>* fid)
{
  //PJSA double subFs[numLocSub][sizeFs][3];
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  Vec<GhostPoint<dim>*> *gp=0;
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    if(ghostPoints) gp = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeRecSurfBasedForceLoadNew(forceApp, orderOfAccuracy, X(iSub), subFs[iSub], sizeFs,
						     (*distLSS)(iSub), pInfty, 
						     Wstarij(iSub), Wstarji(iSub), V(iSub), gp, postFcn,fid?&((*fid)(iSub)):0);
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::computePrdtWCtrlVolRatio(DistSVec<double,dim> &ratioTimesU, DistSVec<double,dim> &U, DistVec<double> &ctrlVol, DistGeoState &geoState) {
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computePrdtWCtrlVolRatio(ratioTimesU(iSub), U(iSub), ctrlVol(iSub), geoState(iSub));

}

//-------------------------------------------------------------------------------
//---------------------- LEVEL SET (PHI) ----------------------------------------
//-------------------------------------------------------------------------------


template<int dimLS>
void Domain::avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->avoidNewPhaseCreation(Phi(iSub), Phin(iSub));

}

//------------------------------------------------------------------------------
template<int dimLS>
void Domain::avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin, DistVec<double> &weight, DistLevelSetStructure *distLSS){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->avoidNewPhaseCreation(Phi(iSub), Phin(iSub),weight(iSub), distLSS ? &((*distLSS)(iSub)) : 0);
}
//------------------------------------------------------------------------------

template<int dimLS>
void Domain::setupPhiVolumesInitialConditions(const int volid, const int fluidId, DistSVec<double,dimLS> &Phi){

  // It is assumed that the initialization using volumes is only
  // called to distinguish nodes that are separated by a material
  // interface (structure). Thus one node cannot be at
  // the boundary of two fluids. A fluid node then gets its
  // id from the element id and there cannot be any problem
  // for parallelization.
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->setupPhiVolumesInitialConditions(volid, fluidId, Phi(iSub));

}

//------------------------------------------------------------------------------
template<int dimLS>
void Domain::TagInterfaceNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,dimLS> &Phi,
                               int level)
{

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->TagInterfaceNodes(lsdim, Tag(iSub),Phi(iSub),level);
    subDomain[iSub]->sndData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
  }

  levelPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->maxRcvData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));



}
//------------------------------------------------------------------------------
/*template<int dimLS>
void Domain::FinishReinitialization(DistVec<int> &Tag, DistSVec<double,dimLS> &Psi,
                                    int level)
{

	int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->FinishReinitialization(Tag(iSub),Psi(iSub),level);
    subDomain[iSub]->sndData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  levelPat->exchange();
  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->TagPsiExchangeData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)),
                                        *volPat, Psi.subData(iSub));

}*/
//------------------------------------------------------------------------------

template<int dimLS>
void Domain::printPhi(DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi, int it)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printPhi(X(iSub), Phi(iSub), numLocSub);
  com->barrier();
}

//-------------------------------------------------------------------------------

template<int dimLS>
void Domain::computePrdtPhiCtrlVolRatio(DistSVec<double,dimLS> &ratioTimesPhi, DistSVec<double,dimLS> &Phi, DistVec<double> &ctrlVol, DistGeoState &geoState) {

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computePrdtPhiCtrlVolRatio(ratioTimesPhi(iSub), Phi(iSub), ctrlVol(iSub), geoState(iSub));

}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::restrictionOnPhi(DistSVec<double,dim> &initial, DistVec<int> &fluidId,
                              DistSVec<double,dim> &restriction, int fluidIdTarget){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->restrictionOnPhi(initial(iSub),fluidId(iSub),
                                      restriction(iSub), fluidIdTarget);

}

//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::getSignedDistance(int lsdim, DistSVec<double,1> &Psi, DistSVec<double,dimLS> &Phi)
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->getSignedDistance(lsdim, Psi(iSub),Phi(iSub));
}
//------------------------------------------------------------------------------
template<int dimLS>
void Domain::computeDistanceCloseNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,3> &X,
                                       DistNodalGrad<dimLS> &lsgrad,
                                       DistSVec<double,dimLS> &Phi,DistSVec<double,1> &Psi,
                                       MultiFluidData::CopyCloseNodes copy)
{
  if(copy==MultiFluidData::FALSE){
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->computeDistanceCloseNodes(lsdim, Tag(iSub), X(iSub), lsgrad(iSub), Phi(iSub), Psi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));
    subDomain[iSub]->recomputeDistanceCloseNodes(lsdim, Tag(iSub), X(iSub), lsgrad(iSub), Phi(iSub), Psi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));
  }else{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->copyCloseNodes(lsdim, 1,Tag(iSub),Phi(iSub),Psi(iSub));
  }

}
//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::computeDistanceLevelNodes(int lsdim, DistVec<int> &Tag, int level,
                                       DistSVec<double,3> &X,DistSVec<double,1> &Psi,
                                       double &_res, DistSVec<double,dimLS> &Phi,
                                       MultiFluidData::CopyCloseNodes copy)
{
  double res(_res);

  if(copy==MultiFluidData::TRUE && level==2){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->copyCloseNodes(lsdim, 2,Tag(iSub),Phi(iSub),Psi(iSub));
    return;
  }

#pragma omp parallel for reduction(+: res)
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    res+=subDomain[iSub]->computeDistanceLevelNodes(lsdim, Tag(iSub), level, X(iSub), Psi(iSub),Phi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));

  com->globalSum(1, &res);
  _res = sqrt(res);
}
//------------------------------------------------------------------------------

template<int dim>
void Domain::setupUVolumesInitialConditions(const int volid, double UU[dim],
                                            DistSVec<double,dim> &U)
{

  // It is assumed that the initialization using volumes is only
  // called to distinguish nodes that are separated by a material
  // interface (structure). Thus one node cannot be at
  // the boundary of two fluids. A fluid node then gets its
  // id from the element id and there cannot be any problem
  // for parallelization.
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->setupUVolumesInitialConditions(volid, UU, U(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::blur(DistSVec<double,dim> &U, DistSVec<double,dim> &U0)
{
  int iSub;

  DistVec<double> loc_weight(getNodeDistInfo());

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->blur(U(iSub),U0(iSub),loc_weight(iSub));
  
  assemble(vecPat, U0);
  assemble(volPat,loc_weight);
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < U.subSize(iSub); ++i) {
      for (int k = 0; k < dim; ++k)
	U0(iSub)[i][k] = 0.5*U0(iSub)[i][k]/loc_weight(iSub)[i] + 0.5*U(iSub)[i][k];
    }

  }

}

template<int dim, class Obj>
void Domain::integrateFunction(Obj* obj,DistSVec<double,3> &X,DistSVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
			       int npt) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->integrateFunction(obj, X(iSub), V(iSub), F, npt);
    subDomain[iSub]->sndData(*volPat, V(iSub) );
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->addRcvData(*volPat, V(iSub));
  }
}
