#include <math.h>

#include <MatVecProd.h>
#include <IoData.h>

#include <BcDef.h>
#include <RecFcnDesc.h>
#include <FluxFcn.h>
#include <DistTimeState.h>
#include <DistGeoState.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MemoryPool.h>

// Included (MB)
#include <FemEquationTermDesc.h>

//------------------------------------------------------------------------------

template<int dim, int neq>
// Included (MB)
MatVecProdFD<dim, neq>::MatVecProdFD(ImplicitData &data, DistTimeState<dim> *ts,
				DistGeoState *gs, SpaceOperator<dim> *spo, 
				Domain *domain, IoData &ioData, bool fdsa) 
  : geoState(gs), Qeps(domain->getNodeDistInfo()), Feps(domain->getNodeDistInfo())
                , Qepstmp(domain->getNodeDistInfo()), Fepstmp(domain->getNodeDistInfo())
                , Q(domain->getNodeDistInfo()), F(domain->getNodeDistInfo())
                , Ftmp(domain->getNodeDistInfo()), iod(&ioData)
{

  com = domain->getCommunicator();

  if (ts)
    timeState = new DistTimeState<dim>(*ts, false, ioData);
  else
    timeState = 0;

  spaceOp = new SpaceOperator<dim>(*spo, false);

  recFcnCon = 0;
  Rn = 0;
  Phi = 0;
  fluidId = 0;
  Riemann = 0;

  if (data.mvp == ImplicitData::H1FD) {
    recFcnCon = new RecFcnConstant<dim>;
    spaceOp->setRecFcn(recFcnCon);
    if (data.type == ImplicitData::CRANK_NICOLSON) {
      Rn = new DistSVec<double,dim>(domain->getNodeDistInfo());
      timeState->setResidual(Rn);
    }
  }

// Included (MB)
  if (fdsa)
    fdOrder = ioData.sa.mvpfdOrdersa; 
  else
    fdOrder = ioData.ts.implicit.fdOrder; 

}

//------------------------------------------------------------------------------

template<int dim, int neq>
MatVecProdFD<dim, neq>::~MatVecProdFD()
{ 

  if (spaceOp) delete spaceOp;
  if (timeState) delete timeState;
  recFcnCon = 0; // deleted by spaceOp
  Rn = 0; // deleted by timeState

}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon) {
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

  }
  else  {
    Feps = f;
  }

  Qeps.strip(Q);
  Feps.strip(F);
  
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon) {
    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  else  {
    Feps = f;
  }

  Qeps.strip(Q);
  Feps.strip(F);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon) {
    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  else  {
    Feps = f;
  }

  Qeps.strip(Q);
  Feps.strip(F);

}

//------------------------------------------------------------------------------
template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv,
                                 DistSVec<double,dim> &q, DistVec<double> &phi,
                                 DistVec<int> &fluidId_, DistExactRiemannSolver<dim> *riemann,
                                 DistSVec<double,dim> &f)
{
  
  X = &x;
  ctrlVol = &cv;
  Qeps = q;
  Phi = &phi;
  fluidId = &fluidId_;
  Riemann = riemann;
                                                                                                                 
  if (recFcnCon) {
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *fluidId, Feps, riemann, 0);
                                                                                                                 
    if (timeState)
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);
                                                                                                                 
    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  else  {
    Feps = f;
  }

  Qeps.strip(Q);
  Feps.strip(F);
  
}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::apply(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  double eps = computeEpsilon(Q, p);

// Included (MB)
  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  if(Phi)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *fluidId, Feps, Riemann, 0);
  else
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);

  if (timeState)
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);

  if (fdOrder == 1) {
    
    prod = (1.0/eps) * (Fepstmp - F);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;
    
    Qepstmp.pad(Qeps);

    if(Phi)
      spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *fluidId, Feps, Riemann, 0);
    else
      spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod = (0.5/eps) * (Fepstmp - Ftmp);

  }

// Original
/*

  //#define MVP_CHECK_ONE_EQ 5
#if MVP_CHECK_ONE_EQ
  int numLocSub = p.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*qeps)[dim] = Qeps.subData(iSub);
    double (*q)[dim] = Q->subData(iSub);
    double (*lp)[dim] = p.subData(iSub);
    for (int i=0; i<p.subSize(iSub); ++i) {
      for (int j=0; j<dim; ++j)
	qeps[i][j] = q[i][j];
      qeps[i][MVP_CHECK_ONE_EQ] = q[i][MVP_CHECK_ONE_EQ] + eps * lp[i][MVP_CHECK_ONE_EQ];
    }
  }
#else
  Qepstmp = Q + eps * p;
  Qepstmp.pad(Qeps);
#endif

  if(Phi)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, Feps, Riemann, 0);
  else
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);


  if (timeState)
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  
  prod = (1.0/eps) * (Fepstmp - F);
*/

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyInviscid(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  double eps = computeEpsilon(Q, p);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

  if (timeState)
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

  Q.pad(Qeps);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  
  if (fdOrder == 1) {

    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);
    
    prod += (1.0/eps) * (Fepstmp - Ftmp);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;

    Qepstmp.pad(Qeps);

    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);

    Q.pad(Qeps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (1.0/eps) * (Fepstmp - Ftmp);

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyViscous(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  double eps = computeEpsilon(Q, p);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  
  if (fdOrder == 1) {

    Q.pad(Qeps);

    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (1.0/eps) * (Fepstmp - Ftmp);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;

    Qepstmp.pad(Qeps);

    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (1.0/eps) * (Fepstmp - Ftmp);

  }

}

//------------------------------------------------------------------------------

template<int dim, int neq>
double MatVecProdFD<dim, neq>::computeEpsilon(DistSVec<double,neq> &U, DistSVec<double,neq> &p)
{

  int iSub, size = 0;
  double eps0 = 1.e-6;

  const DistInfo &distInfo = U.info();

  double *alleps = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) alleps[iSub] = 0.0;

#pragma omp parallel for reduction(+: size)
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    int locsize = 0;
    double loceps = 0.0;

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {

      if (distInfo.masterFlag[locOffset+i]) {
	for (int j=0; j<dim; ++j) {
	  ++locsize;
	  loceps += eps0*fabs(U[locOffset+i][j]) + eps0;
	}
      }

    }

    size += locsize;
    alleps[distInfo.locSubToGlobSub[iSub]] = loceps;

  }

  this->com->globalSum(1, &size);
  this->com->globalSum(distInfo.numGlobSub, alleps);

  double eps = 0.0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) eps += alleps[iSub];

  double norm = sqrt(p*p);

  if (norm > 1.e-14) eps /= double(size) * norm;
  else eps = eps0;
 
  return eps;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::rstSpaceOp(IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, bool typeAlloc, SpaceOperator<dim> *spofd)
{

  spaceOp->rstFluxFcn(ioData);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::MatVecProdH1(DistTimeState<dim> *ts, SpaceOperator<dim> *spo,
					   Domain *domain) : 
  DistMat<Scalar,neq>(domain), timeState(ts), spaceOp(spo)
{

#ifdef _OPENMP 
  this->numLocSub = DistMat<Scalar,neq>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,neq>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,neq>();

    size += double(A[iSub]->numNonZeroBlocks()*neq*neq*sizeof(Scalar)) / (1024.*1024.);
 
  }

  this->com->globalSum(1, &size);
  
  this->com->printf(2, "Memory required for matvec with H1 (dim=%d): %3.2f MB\n", neq, size);

}
//------------------------------------------------------------------------------
                                                                                                  
template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::MatVecProdH1(DistTimeState<dim> *ts, SpaceOperator<dim> *spo,
                                           Domain *domain, IoData &ioData) :
  DistMat<Scalar,neq>(domain), timeState(ts), spaceOp(spo)
{

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,neq>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,neq>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,neq>();
    size += double(A[iSub]->numNonZeroBlocks()*neq*neq*sizeof(Scalar)) / (1024.*1024.);
  }

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory required for matvec with H1 (dim=%d): %3.2f MB\n", neq, size);
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::~MatVecProdH1()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
DistMat<Scalar,neq> &MatVecProdH1<dim,Scalar,neq>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::exportMemory(MemoryPool *mp)
{

  if (!mp) return;

  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    mp->set(A[iSub]->numNonZeroBlocks() * neq*neq * sizeof(Scalar), 
	    reinterpret_cast<void *>(A[iSub]->data()));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluate(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
					    DistSVec<double,dim> &Q, DistSVec<double,dim> &F)
{

  spaceOp->computeJacobian(X, ctrlVol, Q, *this, timeState);

  if (timeState)
    timeState->addToJacobian(ctrlVol, *this, Q);

  spaceOp->applyBCsToJacobian(Q, *this);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluate(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                            DistSVec<double,dim> &Q, DistVec<double> &Phi, DistVec<int> &fluidId,
                                            DistExactRiemannSolver<dim> *riemann,
                                            DistSVec<double,dim> &F)
{

  spaceOp->computeJacobian(X, ctrlVol, Q, *this, fluidId, riemann);

  if (timeState)
    timeState->addToJacobian(ctrlVol, *this, Q);

  spaceOp->applyBCsToJacobian(Q, *this);

}

//------------------------------------------------------------------------------
                                                                                                      
template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluateViscous(int it, DistSVec<double,3> &X,
                                   DistVec<double> &cv)  {

  spaceOp->computeViscousJacobian(X, cv, *this);

}

//------------------------------------------------------------------------------
// note: this can be done in another way (but less efficient) !!
// (1) compute off-diag products (2) assemble (3) compute diag products (->redundancy)

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::apply(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.getMasterFlag(iSub), *A[iSub],
					 p(iSub), prod(iSub));
    this->subDomain[iSub]->sndData(*this->vecPat, prod.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.subData(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::rstSpaceOp(IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, bool typeAlloc, SpaceOperator<dim> *spofd)
{

  spaceOp = spo;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
// Included (MB)
MatVecProdH2<Scalar,dim>::MatVecProdH2(IoData &ioData, VarFcn *varFcn, DistTimeState<dim> *ts,
				       SpaceOperator<dim> *spo, Domain *domain, DistGeoState *gs) :
  DistMat<Scalar,dim>(domain), timeState(ts),
  aij(domain->getEdgeDistInfo()), aji(domain->getEdgeDistInfo()), 
  bij(domain->getEdgeDistInfo()), bji(domain->getEdgeDistInfo())
{

#ifdef _OPENMP 
  this->numLocSub = DistMat<Scalar,dim>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,dim>*[this->numLocSub];

  double size = 0.0;
  double coefsize = double(4*aij.size()*dim*sizeof(double)) / (1024.*1024.);

  // allocate for viscous flux jacobian term
  bool nsFlag = false;
  spaceOp = new SpaceOperator<dim>(*spo, false);

// Included (MB*)
  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.bc.wall.integration != BcsWallData::WALL_FUNCTION))
    nsFlag = true;

// Original
/*
  if (ioData.eqs.type == EquationsData::NAVIER_STOKES)  {
    R = new MatVecProdH1<dim, Scalar ,dim>(ts, spo, domain);
    nsFlag = true;
  }
  else
    R = 0;
*/

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,dim>(nsFlag);

    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(Scalar)) / (1024.*1024.);
 
  }

  this->com->globalSum(1, &size);
  this->com->globalSum(1, &coefsize);
  
  this->com->printf(2, "Memory required for matvec with H2: ");
  this->com->printf(2, "%3.2f+%3.2f=%3.2f MB\n", size, coefsize, size+coefsize);
  

// Included (MB)
  fluxFcn = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1]; 
  fluxFcn -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 12)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");
  fluxFcn[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  spaceOp->setFluxFcn(fluxFcn); //TODO: should avoid doing this!

// Included (MB)
  if (ioData.eqs.type == EquationsData::NAVIER_STOKES)  {
    viscJacContrib = ioData.sa.viscJacContrib;
    vProd = new DistSVec<double,dim>(domain->getNodeDistInfo());

    if (viscJacContrib == 1) {
      R = new MatVecProdH1<dim, Scalar ,dim>(ts, spo, domain);
      RFD = 0;
    }
    else if (viscJacContrib == 2) {
      R = 0;
      if (gs)
        RFD = new MatVecProdFD<dim,dim>(ioData.ts.implicit, ts, gs, spo, domain, ioData);
      else
        RFD = 0;
    }
    else {
      vProd = 0;
      R = 0;
      RFD = 0;
    }
  }
  else {
    vProd = 0;
    R = 0;
    RFD = 0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
MatVecProdH2<Scalar,dim>::~MatVecProdH2()
{ 

  if (spaceOp) delete spaceOp;
  fluxFcn = 0; // deleted by spaceOperator.

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------
template<class Scalar, int dim>
DistMat<Scalar,dim> &MatVecProdH2<Scalar,dim>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

// Included (MB)
  evaluateInviscid(it, x, cv, q, f);
  evaluateViscous(it, x, cv, q, f);

// Original
/*
  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState)
    timeState->addToH2(*ctrlVol, *Q, *this);
  


  // compute viscous flux jacobian
  if (R)  {
    spaceOp->applyBCsToH2Jacobian(*Q, *this);
    R->evaluateViscous(it, *X, *ctrlVol);
    spaceOp->applyBCsToJacobian(*Q, *R);
  }
*/

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState)
    timeState->addToH2(*ctrlVol, *Q, *this);

  spaceOp->applyBCsToH2Jacobian(*Q, *this);

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  // compute viscous flux jacobian
  if (R)  {
    R->evaluateViscous(it, *X, *ctrlVol);
    spaceOp->applyBCsToJacobian(*Q, *R);
  }

  if (RFD) {
    F = &f;
    RFD->evaluateViscous(it, *X, *ctrlVol, *Q, *F);
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evaluate(int it, DistSVec<double,3> &x, 
                               DistVec<double> &cv, DistSVec<double,dim> &q, 
                               DistSVec<double,dim> &F, Scalar shift)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState) {
     switch (it)  {
       //case for the construction of the POD
       case 0:
         timeState->addToH2(*ctrlVol, *Q, *this, shift, 1.0);
         break;   
       case 1:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), 1.0);
         break;
       case 2:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), -1.0);
         break;
       case 3:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(3.0), 2.0);
         break;
       case 4:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(4.0), -2.0);
         break;
       case 5:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), -2.0);
         break;
       case 6:
         timeState->addToH2(*ctrlVol, *Q, *this, (Scalar(4.0)*shift+Scalar(2.0))/(shift*(shift+Scalar(1.0))), 2.0);
         break;
       case 7:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(8.0/3.0), 2.0);
       break;

    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evaluate2(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
                                         DistSVec<double,dim> &q, DistSVec<double,dim> &F)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState) {
    timeState->addToH2Minus(*ctrlVol, *Q, *this);
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::evalH(int it, DistSVec<double,3> &x,
                               DistVec<double> &cv, DistSVec<double,dim> &q)  {

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::apply(DistSVec<double,dim> &p, DistSVec<double,dim> &prod)
{

// Included (MB)
  applyInviscid(p,prod);
  applyViscous(p,prod);

// Original
/*
  spaceOp->applyH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

  if (R)  {
    DistSVec<double, dim> vProd(p);
    vProd = 0.0;
    R->apply(p, vProd);
    prod += vProd;
  }
*/

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::applyInviscid(DistSVec<double,dim> &p, DistSVec<double,dim> &prod)
{

  spaceOp->applyH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::applyViscous(DistSVec<double,dim> &p, DistSVec<double,dim> &prod)
{

  if (R)  {
    *vProd = 0.0;
    R->apply(p, *vProd);
    prod += *vProd;
  }

  if (RFD)  {
    *vProd = 0.0;
    RFD->applyViscous(p, *vProd);
    prod += *vProd;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::apply(DistSVec<bcomp,dim> &p,
                DistSVec<bcomp,dim> &prod)
{

  spaceOp->applyH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::applyT(DistSVec<double,dim> &p,
        DistSVec<double,dim> &prod)
{

  spaceOp->applyH2T(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::applyT(DistSVec<bcomp,dim> &p,
        DistSVec<bcomp,dim> &prod)
{

  spaceOp->applyH2T(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar, int dim>
void MatVecProdH2<Scalar,dim>::rstSpaceOp(IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, bool typeAlloc, SpaceOperator<dim> *spofd)
{

  fluxFcn = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1]; 
  fluxFcn -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 12)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");

  fluxFcn[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn,FluxFcnBase::PRIMITIVE);

  spaceOp->setFluxFcn(fluxFcn);

}

//------------------------------------------------------------------------------

template<class Scalar,int dim, int neq>
MatVecProdLS<Scalar,dim,neq>::MatVecProdLS(IoData &ioData, VarFcn *varFcn, 
                                           DistTimeState<dim> *ts, DistGeoState *gs,
                                           SpaceOperator<dim> *spo, Domain *domain, 
                                           LevelSet *ls) :
  DistMat<Scalar,neq>(domain), timeState(ts), geoState(gs), LS(ls),
  aij(domain->getEdgeDistInfo()), aji(domain->getEdgeDistInfo()),
  bij(domain->getEdgeDistInfo()), bji(domain->getEdgeDistInfo()), 
  Qeps(domain->getNodeDistInfo()), Feps(domain->getNodeDistInfo()),
  QepsV(domain->getNodeDistInfo()), QV(domain->getNodeDistInfo())
{

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,neq>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,neq>*[this->numLocSub];

  double size = 0.0;
  double coefsize = double(4*aij.size()*sizeof(double)) / (1024.*1024.);

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,neq>();
    size += double(A[iSub]->numNonZeroBlocks()*sizeof(Scalar)) / (1024.*1024.);
  }

  this->com->globalSum(1, &size);
  this->com->globalSum(1, &coefsize);

  spaceOp = new SpaceOperator<dim>(*spo, false);
  timeState = new DistTimeState<dim>(*ts, false, ioData);
}

//------------------------------------------------------------------------------
                                                                                                                      
template<class Scalar,int dim, int neq>
MatVecProdLS<Scalar,dim,neq>::~MatVecProdLS()
{
                                                                                                                      
  if (spaceOp) delete spaceOp;
  if (timeState) delete timeState;                                                                                                                      
  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub)
      if (A[iSub]) delete A[iSub];
                                                                                                                      
    delete [] A;
  }
                                                                                                                      
}
//------------------------------------------------------------------------------
// note: this can be done in another way (but less efficient) !!
// (1) compute off-diag products (2) assemble (3) compute diag products (->redundancy)

template<class Scalar,int dim, int neq>
void MatVecProdLS<Scalar,dim,neq>::applyLS(DistVec<double> &p, DistVec<double> &prod)
{
  double eps = computeEpsilon(*Q, p);

  Qeps = (*Q) + eps * p;

  //spaceOp->computeResidualLS(*X, *ctrlVol, Qeps, *U, Feps);
  fprintf(stderr, "*** Error: need to properly call function ::applyLS in MatVecProdLS\n");
  exit(1);

  timeState->add_dAW_dtLS(-1, *geoState, *ctrlVol, Qeps, *Qn, *Qnm1, *Qnm2, Feps);

  prod = (1.0/eps) * (Feps - (*F));

}

//------------------------------------------------------------------------------

template<class Scalar,int dim, int neq>
void MatVecProdLS<Scalar,dim,neq>::evaluateLS(int it, DistSVec<double,3> &x, DistVec<double> &cv,
                                                      DistVec<double> &q,  DistVec<double> &qn,
                                                      DistVec<double> &qnm1, DistVec<double> &qnm2,
						      DistSVec<double,dim> &u, DistVec<double> &f,
                                                      DistVec<int> &fluidId)
{
                                                                                                                    
  X = &x;
  ctrlVol = &cv;
  Q = &q;
  Qn= &qn;
  Qnm1= &qnm1;
  Qnm2= &qnm2;
  U = &u;
  F = &f;                                                                                                                    
  //spaceOp->computeResidualLS(*X, *ctrlVol, *Q, *U, *F);
  fprintf(stderr, "*** Error: need to properly call function ::evaluateLS in MatVecProdLS\n");
  exit(1);
                                                                                                             
  timeState->add_dAW_dtLS(it, *geoState, *ctrlVol, *Q, *Qn, *Qnm1, *Qnm2, *F);
                                                                                                             
                                                                                                                    
}
//------------------------------------------------------------------------------
template<class Scalar,int dim, int neq>
DistMat<Scalar,neq> &MatVecProdLS<Scalar,dim,neq>::operator= (const Scalar x)
{
                                                                                                                      
#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;
                                                                                                                      
  return *this;

}
//------------------------------------------------------------------------------
template<class Scalar,int dim, int neq>
double MatVecProdLS<Scalar,dim,neq>::computeEpsilon(DistVec<double> &u, DistVec<double> &p)                                                                                                             
{
                                                                                                             
  int iSub, size = 0;
  double eps0 = 1.e-6;
                                                                                                             
  const DistInfo &distInfo = u.info();
                                                                                                             
  double *alleps = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));
                                                                                                             
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) alleps[iSub] = 0.0;
                                                                                                             
#pragma omp parallel for reduction(+: size)
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
                                                                                                             
    int locOffset = distInfo.subOffset[iSub];
                                                                                                             
    int locsize = 0;
    double loceps = 0.0;
                                                                                                             
    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
                                                                                                             
      if (distInfo.masterFlag[locOffset+i]) {
        for (int j=0; j<1; ++j) {
          ++locsize;
          loceps += eps0*fabs(u[locOffset+i]) + eps0;
        }
      }
                                                                                                             
    }
                                                                                                             
    size += locsize;
    alleps[distInfo.locSubToGlobSub[iSub]] = loceps;
                                                                                                             
  }
                                                                                                             
  this->com->globalSum(1, &size);
  this->com->globalSum(distInfo.numGlobSub, alleps);
                                                                                                             
  double eps = 0.0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) eps += alleps[iSub];
                                                                                                             
  double norm = sqrt(p*p);
                                                                                                             
  if (norm > 1.e-14) eps /= double(size) * norm;
  else eps = eps0;
                                                                                                             
  return eps;
                                                                                                             
}

//------------------------------------------------------------------------------
