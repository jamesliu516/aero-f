#include <DistTimeState.h>

#include <IoData.h>
#include <VarFcn.h>
#include <TimeState.h>
#include <Domain.h>
#include <DistGeoState.h>
#include <DistMatrix.h>

#include <math.h>
//------------------------------------------------------------------------------

template<int dim>
DistTimeState<dim>::DistTimeState(IoData &ioData, SpaceOperator<dim> *spo, VarFcn *vf,
				  Domain *dom, DistSVec<double,dim> *v) 
  : varFcn(vf), domain(dom)
{
  locAlloc = true;

  if (v) V = v->alias();
  else V = new DistSVec<double,dim>(domain->getNodeDistInfo());

  numLocSub = domain->getNumLocSub();

  data = new TimeData(ioData);  

  dt  = new DistVec<double>(dom->getNodeDistInfo());
  idti = new DistVec<double>(dom->getNodeDistInfo());
  idtv = new DistVec<double>(dom->getNodeDistInfo());
  irey  = new DistVec<double>(dom->getNodeDistInfo());
  viscousCst = ioData.ts.viscousCst;
  Un  = new DistSVec<double,dim>(domain->getNodeDistInfo());

  if (data->use_nm1)
    Unm1 = new DistSVec<double,dim>(domain->getNodeDistInfo());
  else
    Unm1 = Un->alias();

  if (data->use_nm2)
    Unm2 = new DistSVec<double,dim>(domain->getNodeDistInfo());
  else
    Unm2 = Unm1->alias();

  if (ioData.eqs.tc.les.type == LESModelData::DYNAMICVMS) {
    QBar = new DistSVec<double,dim>(domain->getNodeDistInfo());
    VnBar = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Vn = new DistSVec<double,dim>(domain->getNodeDistInfo());
    UnBar = new DistSVec<double,dim>(domain->getNodeDistInfo());
    if (data->use_nm1)
      Unm1Bar = new DistSVec<double,dim>(domain->getNodeDistInfo());
    else
      Unm1Bar = UnBar->alias();
    if (data->use_nm2)
      Unm2Bar = new DistSVec<double,dim>(domain->getNodeDistInfo());
    else
      Unm2Bar = Unm1Bar->alias();
  }
  else {
    QBar = 0;
    VnBar = 0;
    Vn = 0;
    UnBar = 0;
    Unm1Bar = 0;
    Unm2Bar = 0;
  }
                                                                                                                          
  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON)
    Rn = new DistSVec<double,dim>(domain->getNodeDistInfo());
  else
    Rn = Un->alias();

  gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

  fet = spo->getFemEquationTerm();

  //preconditioner setup
  tprec.setup(ioData);
  sprec.setup(ioData);

  subTimeState = new TimeState<dim>*[numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) 
    subTimeState[iSub] = 0;

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {
    dIdti = new DistVec<double>(dom->getNodeDistInfo());
    dIdtv = new DistVec<double>(dom->getNodeDistInfo());
    dIrey = new DistVec<double>(dom->getNodeDistInfo());
    *dIdti = 0.0; 
    *dIdtv = 0.0; 
    *dIrey = 0.0;
  }
  else {
    dIdti = 0; 
    dIdtv = 0; 
    dIrey = 0;
  }

}

//------------------------------------------------------------------------------

template<int dim>
DistTimeState<dim>::DistTimeState(const DistTimeState<dim> &ts, bool typeAlloc, IoData &ioData) 
{

  locAlloc = typeAlloc;

  gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

  //preconditioner setup
  tprec.setup(ioData);
  sprec.setup(ioData);

  varFcn = ts.varFcn;
  fet = ts.fet;
  V = ts.V;

  numLocSub = ts.numLocSub;

  data = ts.data;

  dt  = ts.dt;
  idti = ts.idti;
  idtv = ts.idtv;
  irey = ts.irey;
  viscousCst = ts.viscousCst;
  Un  = ts.Un;
  Unm1 = ts.Unm1;
  Unm2 = ts.Unm2;
  Rn = ts.Rn;

  QBar = ts.QBar;
  VnBar = ts.VnBar;
  Vn = ts.Vn;
  UnBar = ts.UnBar;
  Unm1Bar = ts.Unm1Bar;
  Unm2Bar = ts.Unm2Bar;

  domain = ts.domain;

  subTimeState = ts.subTimeState;

}

//------------------------------------------------------------------------------

template<int dim>
DistTimeState<dim>::~DistTimeState()
{

  if (!locAlloc) return;

  if (data) delete data;
  if (V) delete V;
  if (Un) delete Un;
  if (Unm1) delete Unm1;
  if (Unm2) delete Unm2;
  if (Rn) delete Rn;
  if (QBar) delete QBar;
  if (Vn) delete Vn;
  if (VnBar) delete VnBar;
  if (UnBar) delete UnBar;
  if (Unm1Bar) delete Unm1Bar;
  if (Unm2Bar) delete Unm2Bar;
                                                                                                                          
  if (subTimeState) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      if (subTimeState[iSub]) 
	delete subTimeState[iSub]; 
    
    delete [] subTimeState;
  }

}

//------------------------------------------------------------------------------
/*
template<int dim>
void DistTimeState<dim>::setup(const char *name, double *Ucst, 
                               DistSVec<double,3> &X, DistSVec<double,dim> &U)
{
  Un->set(Ucst);

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1)
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }

  U = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;
                                                                                                               
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub),
					      (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));
}
*/
//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::setup(const char *name, DistSVec<double,dim> &Ufar,
			       DistSVec<double,3> &X, DistSVec<double,dim> &U)
{
  *Un = Ufar;

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1)
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }

  U = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub),
                                              (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));
}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setup(const char *name, DistSVec<double,3> &X,
                               DistSVec<double,dim> &Ufar,
                               DistSVec<double,dim> &U, IoData &iod)
{

  *Un = Ufar;

  setupUVolumesInitialConditions(iod);
  setupUMultiFluidInitialConditions(iod,X);

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1)
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }

  U  = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub),
                                              (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));

}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setupUVolumesInitialConditions(IoData &iod)
{

  // loop on all Volumes to setup U0
  if(!iod.volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator it;
    for (it=iod.volumes.volumeMap.dataMap.begin(); it!=iod.volumes.volumeMap.dataMap.end();it++)
      if(it->second->type==VolumeData::FLUID)
        //each volume (it->first) is setup using Input variables 'volumeInitialConditions'
        //                                 and equation of state 'fluidModel'
        domain->setupUVolumesInitialConditions(it->first, it->second->fluidModel, it->second->volumeInitialConditions, *Un);
  }

}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X)
{
  // These initial conditions are setup using MultiFluid.FluidModel2
  // which can be a gas or a liquid or a jwl gas.
  // Note that Un was already initialized either using the far-field
  // or using the definition of volumes. Therefore the only initialization
  // left to do is for the spheres and other geometric shapes.

  if(iod.mf.initialConditions.nspheres>0)
    for(int i=0; i<iod.mf.initialConditions.nspheres; i++)
      domain->setupUMultiFluidInitialConditionsSphere(iod.mf.fluidModel2, *(iod.mf.initialConditions.sphere[i]),X,*Un);

  if(iod.mf.initialConditions.nplanes>0)
    domain->setupUMultiFluidInitialConditionsPlane(iod.mf.fluidModel2, iod.mf.initialConditions.p1, X, *Un);

}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setup(const char *name, DistSVec<double,3> &X,
                               DistSVec<double,dim> &Ufar,
                               DistSVec<double,dim> &U, DistVec<int> &nodeTag, IoData &iod)
{

  *Un = Ufar;
  nodeTag = 1;

  if (!iod.volumes.volumeMap.dataMap.empty()) {
    Communicator *com = domain->getCommunicator(); //for screen output only.
    com->fprintf(stderr,"WARNING: volume map exists but will be ignored!\n");
  }
  setupUMultiFluidInitialConditions(iod,X,nodeTag);

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1)
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }

  U  = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub),
                                              (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));

}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X, 
                                                           DistVec<int> &nodeTag)
{
  // These initial conditions are setup using MultiFluid.FluidModel2
  // which can be a gas or a liquid or a jwl gas.
  // Note that Un was already initialized using the far-field. Therefore the only initialization
  // left to do is for the spheres and other geometric shapes.

  if (iod.mf.initialConditions.nspheres>0) { 
    Communicator *com = domain->getCommunicator(); //for screen output only.
    com->fprintf(stderr,"WARNING: sphere domains exist but will be ignored!\n");
  }

  if(iod.mf.initialConditions.nplanes>0)
    domain->setupUMultiFluidInitialConditionsPlane(iod.mf.fluidModel2, iod.mf.initialConditions.p1, 
                                                   X, *Un, nodeTag);

}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setup(const char *name, DistSVec<double,dim> &Ufar,
                               double *Ub, DistSVec<double,3> &X,
                               DistVec<double> &Phi,
                               DistSVec<double,dim> &U, IoData &iod)
{
  *Un = Ufar;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub){
    double (*x)[3] = X.subData(iSub);
    double (*u)[dim] = Un->subData(iSub);
    double *phi = Phi.subData(iSub);
    for (int i=0; i<X.subSize(iSub); i++){
      if (phi[i]<0.0)
        for (int j=0; j<dim; j++)
          u[i][j] = Ub[j];
      phi[i] *= u[i][0];
    }
  }

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1)
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }
                                                                                                                      
  U  = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;
                                                                                                                      
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub),
					      (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistTimeState<dim>::rstVar(IoData & ioData) 
{

  //mach = ioData.bc.inlet.mach;

}

//------------------------------------------------------------------------------

template<int dim>
double DistTimeState<dim>::computeTimeStep(double cfl, double* dtLeft, int* numSubCycles,
					   DistGeoState &geoState, DistSVec<double,3> &X,
					   DistVec<double> &ctrlVol, DistSVec<double,dim> &U)
{

  varFcn->conservativeToPrimitive(U, *V);

  domain->computeTimeStep(cfl, viscousCst, fet, varFcn, geoState, X, ctrlVol, *V, *dt, *idti, *idtv, *irey, tprec, sprec);

  double dt_glob;
  if (data->dt_imposed > 0.0) 
    dt_glob = data->dt_imposed;
  else 
    dt_glob = dt->min();

  if (data->typeStartup == ImplicitData::MODIFIED && 
      ((data->typeIntegrator == ImplicitData::THREE_POINT_BDF && !data->exist_nm1) ||
       (data->typeIntegrator == ImplicitData::FOUR_POINT_BDF && (!data->exist_nm2 || !data->exist_nm1)))) {
    if (*dtLeft != 0.0 && dt_glob > *dtLeft)
      dt_glob = *dtLeft / 1000.0;
    else
      dt_glob /= 1000.0;
  }

  if (*dtLeft != 0.0) {
    *numSubCycles = int(*dtLeft / dt_glob);
    if (*numSubCycles == 0 || (*numSubCycles)*dt_glob != *dtLeft) ++(*numSubCycles);
    dt_glob = *dtLeft / double(*numSubCycles);
    *dtLeft -= dt_glob;
  }
  data->computeCoefficients(*dt, dt_glob);

  return dt_glob;

}

//------------------------------------------------------------------------------
                                                                                                         
template<int dim>
double DistTimeState<dim>::computeTimeStep(double cfl, double* dtLeft, int* numSubCycles,
                                           DistGeoState &geoState, DistVec<double> &ctrlVol,
                                           DistSVec<double,dim> &U, DistVec<int> &fluidId)
{
  varFcn->conservativeToPrimitive(U, *V, &fluidId);

  domain->computeTimeStep(cfl, viscousCst, fet, varFcn, geoState, ctrlVol, *V, *dt, *idti, *idtv, tprec, fluidId);
                                                                                                         
  double dt_glob;
  if (data->dt_imposed > 0.0)
    dt_glob = data->dt_imposed;
  else
    dt_glob = dt->min();
                                                                                                         
  if (data->typeStartup == ImplicitData::MODIFIED &&
      ((data->typeIntegrator == ImplicitData::THREE_POINT_BDF && !data->exist_nm1) ||
       (data->typeIntegrator == ImplicitData::FOUR_POINT_BDF && (!data->exist_nm2 || !data->exist_nm1)))) {
    if (*dtLeft != 0.0 && dt_glob > *dtLeft)
      dt_glob = *dtLeft / 1000.0;
    else
      dt_glob /= 1000.0;
  }
  if (*dtLeft != 0.0) {
    *numSubCycles = int(*dtLeft / dt_glob);
    if (*numSubCycles == 0 || (*numSubCycles)*dt_glob != *dtLeft) ++(*numSubCycles);
    dt_glob = *dtLeft / double(*numSubCycles);
    *dtLeft -= dt_glob;
  }
                                                                                                         
  data->computeCoefficients(*dt, dt_glob);
                                                                                                         
  return dt_glob;
                                                                                                         
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::computeCoefficients(double dt_glob)
{

  data->computeCoefficients(*dt, dt_glob);

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::add_dAW_dt(int it, DistGeoState &geoState, 
					    DistVec<double> &ctrlVol,
					    DistSVec<double,dim> &Q, 
					    DistSVec<double,dim> &R)
{

  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->add_dAW_dt(Q.getMasterFlag(iSub), geoState(iSub), 
				   ctrlVol(iSub), Q(iSub), R(iSub));

}
                                                                                                                      
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::add_dAW_dtLS(int it, DistGeoState &geoState,
                                            DistVec<double> &ctrlVol,
                                            DistVec<double> &Q,
					    DistVec<double> &Qn,
					    DistVec<double> &Qnm1,
					    DistVec<double> &Qnm2,
                                            DistVec<double> &R)
{
                                                                                                                      
  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;
                                                                                                                      
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->add_dAW_dtLS(Q.getMasterFlag(iSub), geoState(iSub),
                                   ctrlVol(iSub), Q(iSub), Qn(iSub), Qnm1(iSub),
				   Qnm2(iSub), R(iSub));
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobian(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  if(tprec.timePreconditioner()){
    if(varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
      addToJacobianGasPrec(ctrlVol, A, U);
    else if(varFcn->getType() == VarFcnBase::TAIT)
      addToJacobianLiquidPrec(ctrlVol, A, U);
    else{
      fprintf(stdout, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
      exit(1);
    }
  }else{
     addToJacobianNoPrec(ctrlVol, A, U);
  }

}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianNoPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianNoPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianNoPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, empty);
  }
}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianGasPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianGasPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, gam, pstiff, tprec, (*irey)(iSub), nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianGasPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, gam, pstiff, tprec, (*irey)(iSub), empty);
  }
}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianLiquidPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianLiquidPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, tprec, (*irey)(iSub), nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianLiquidPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, tprec, (*irey)(iSub), empty);
  }
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH1(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH1(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH1(DistVec<double> &ctrlVol,
                DistMat<Scalar,neq> &A, Scalar shift)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH1(V->getMasterFlag(iSub), ctrlVol(iSub),
                A(iSub), shift);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
				 DistMat<Scalar,dim> &A)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub), 
				(*V)(iSub), A(iSub)); 

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
				 DistMat<Scalar,neq> &A)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub), 
				(*V)(iSub), A(iSub)); 

}

//------------------------------------------------------------------------------
                                                                                                                      
template<int dim>
template<class Scalar>
void DistTimeState<dim>::addToH2LS(DistVec<double> &ctrlVol, DistVec<double> &U,
                                 DistMat<Scalar,1> &A)
{
                                                                                                                      
#ifdef DOUBLE_CHECK
  //varFcn->conservativeToPrimitive(U, *V);
#endif
                                                                                                                      
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2LS(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                (*V)(iSub), A(iSub));
                                                                                                                      
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol,
                DistSVec<double,dim> &U, DistMat<Scalar,dim> &A, Scalar shift)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                (*V)(iSub), A(iSub), shift);

}
//------------------------------------------------------------------------------


template<int dim>
template<class Scalar>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol,
                DistSVec<double,dim> &U, DistMat<Scalar,dim> &A, Scalar coefVol, double coefA)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                (*V)(iSub), A(iSub), coefVol, coefA);

}
//------------------------------------------------------------------------------




template<int dim>
template<class Scalar>
void DistTimeState<dim>::addToH2Minus(DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
                                      DistMat<Scalar,dim> &A)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2Minus(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                     (*V)(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::multiplyByTimeStep(DistSVec<double,dim>& dU)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double (*du)[dim] = dU.subData(iSub);
    double* _dt = dt->subData(iSub);
    for (int i=0; i<dU.subSize(iSub); ++i)
      for (int j=0; j<dim; ++j)
	du[i][j] *= _dt[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::multiplyByTimeStep(DistVec<double>& dPhi)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double *dphi = dPhi.subData(iSub);
    double* _dt = dt->subData(iSub);
    for (int i=0; i<dPhi.subSize(iSub); ++i)
	dphi[i] *= _dt[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::multiplyByPreconditioner(DistSVec<double,dim>& U0, DistSVec<double,dim>& dU)
{
  if (tprec.timePreconditioner()){
    if (varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
      multiplyByPreconditionerPerfectGas(U0,dU);
    else if (varFcn->getType() == VarFcnBase::TAIT)
      multiplyByPreconditionerLiquid(U0,dU);
    else{
      fprintf(stdout, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
      exit(1);
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::multiplyByPreconditionerPerfectGas(DistSVec<double,dim>& U0, DistSVec<double,dim>& dU)
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double* _irey = irey->subData(iSub);
    double (*u0)[dim] = U0.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    for (int i=0; i<dU.subSize(iSub); ++i) {
        double ro = u0[i][0];
        double invRho = 1.0/ro;
        double u  = u0[i][1] * invRho;
        double v  = u0[i][2] * invRho;
        double w  = u0[i][3] * invRho;
        double u2 = u*u;
        double v2 = v*v;
        double w2 = w*w;
        double q2 = u2 + v2 + w2;
        double gam1 = gam - 1.0;
        double p  = gam1 * (u0[i][4] - 0.5 * ro * q2)-gam*pstiff;
        double c2 = gam*(p+pstiff)/ro;

        double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
        double locbeta = tprec.getBeta(locMach,_irey[i]);
        //double locbeta = fmax(k1*locMach, beta);
        //locbeta = fmin((1.0+sqrt(_irey[i]))*locbeta,cmach);

        double beta2 =  locbeta * locbeta;
        double qhat2 = (q2 * gam1)/2.0;
                                                                                
        double nu = qhat2/c2;
        double mu = beta2 - 1.0;
        double H = (c2/gam1) + 0.5*q2;
        
 	// initialize output
	double temp[dim];
        for (int j=0; j<dim; j++)
          temp[j] = 0.0;

        // Euler or Navier-Stokes preconditioning
        double P[5][5] = {  {1.0 + nu*mu,  -u*mu*gam1/c2,      -v*mu*gam1/c2,      -w*mu*gam1/c2,       mu*gam1/c2   },
                            {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                            {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                            {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
                            {qhat2*H*mu/c2,  -u*mu*(1+nu),      -v*mu*(1+nu),       -w*mu*(1+nu), 1.0 + mu*gam1*H/c2 } };

       for (int j=0; j<5; ++j)
          for (int k=0; k<5; ++k)
             temp[j] += P[j][k]*du[i][k];

	//turbulence preconditioning
        if(dim==6){
          double t1 = u0[i][5] * invRho;
          double mup = mu*t1*gam1/c2;
          double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
          for (int k=0; k<6; k++)
            temp[5] += Pt[k]*du[i][k];

        }else if(dim==7){
          double t1 = u0[i][5] * invRho;
	  double t2 = u0[i][6] * invRho;
          double mup1 = mu*t1*gam1/c2;
          double mup2 = mu*t2*gam1/c2;
          double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
			      {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
          for (int k=0; k<7; k++){
            temp[5] += Pt[0][k]*du[i][k];
            temp[6] += Pt[1][k]*du[i][k];
          }
	}

	//get output
	for(int j=0; j<dim; ++j)
	   du[i][j] = temp[j];
    
    }
    
  }
}
                                                                                
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::multiplyByPreconditionerLiquid(DistSVec<double,dim> &U, DistSVec<double,dim> &dU)
{
//ARL : turbulence preconditioner never tested...
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double* _irey = irey->subData(iSub);
      double (*u0)[dim] = U.subData(iSub);
      double (*du)[dim] = dU.subData(iSub);
      for (int i=0; i<dU.subSize(iSub); ++i) {
        double V[dim];
        varFcn->conservativeToPrimitive(u0[i],V);
        double pressure = varFcn->getPressure(V);
        double e = varFcn->computeRhoEnergy(V)/V[0];
        double locMach = varFcn->computeMachNumber(V); //local Preconditioning (ARL)
        double locbeta = tprec.getBeta(locMach,_irey[i]);
        //double locbeta = fmax(k1*locMach, beta);
        //locbeta = fmin((1.0+sqrt(_irey[i]))*locbeta, cmach);
        double beta2 = (locbeta*locbeta);
        double beta2m1 = beta2 - 1.0;
        
	// initialize output        
        double temp[dim];
        for (int j=0; j<dim; j++)
          temp[j] = 0.0;

	// preconditioning matrix
        double P[dim][dim];
 	// first initialize to identity
        for (int j=0; j<dim; j++)
          for (int k=0; k<dim; k++)
	     P[j][k] = 0.0;
        for (int j=0; j<dim; j++)
           P[j][j] = 1.0;
	// second get first column terms
        for (int j=0; j<dim; j++)
	  P[j][0] = beta2m1*V[j];
        P[0][0] = beta2;
        P[4][0] = beta2m1*e;
        
	/* The preconditioning matrix is:
         * P[dim][dim] = { { beta2,           0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[1], 1.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[2], 0.0, 1.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[3], 0.0, 0.0, 1.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*e,    0.0, 0.0, 0.0, 1.0 , 0.0, 0.0 },
	 *		   {(beta2-1.0)*V[5], 0.0, 0.0, 0.0, 0.0 , 1.0, 0.0 },
	 *		   {(beta2-1.0)*V[6], 0.0, 0.0, 0.0, 0.0 , 0.0, 1.0 } };
	 * Take the first 5-by-5 matrix to get the Euler preconditioner
	 *      the first 6-by-6 matrix to get the "SA"  preconditioner
	 *      the whole 7-by-7 matrix to get the "k-e" preconditioner
         */
        for (int j=0; j<dim; ++j)
          for (int k=0; k<j+1; ++k)
            temp[j] += P[j][k]*du[i][k];

        for(int j=0; j<dim; ++j)
          du[i][j] = temp[j];
      }                                                                                                                                        
    }
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::update(DistSVec<double,dim> &Q)
{

  data->update();

  if (data->use_nm2 && data->exist_nm1) {
    *Unm2 = *Unm1;
    data->exist_nm2 = true;
  }
  if (data->use_nm1) {
    *Unm1 = *Un;
    data->exist_nm1 = true;
  }
  *Un = Q;
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::update(DistSVec<double,dim> &Q, DistVec<int> &fluidId,
                                DistVec<int> *fluidIdnm1,
                                //DistSVec<double,dim> *Vgf, DistVec<double> *Vgfweight,
                                DistExactRiemannSolver<dim> *riemann)
{
  data->update();

  if (data->use_nm2 && data->exist_nm1) {
    fprintf(stdout, "4pt-BDF has not been studied for 2-phase flow\n");
    exit(1);
  }
  if (data->use_nm1) {
    fprintf(stderr, "*** Error: 3pt-BDF for multiphase must be reviewed\n");
    exit(1);
    varFcn->conservativeToPrimitive(*Un, *V, fluidIdnm1);
    //TODO: to be fixed for 3pt-BDF!!!
    //varFcn->updatePhaseChange(*V, *Unm1, fluidId, fluidIdnm1, Vgf, Vgfweight, riemann);
    //riemann->updatePhaseChange(*V, *Unm1, fluidId, fluidIdnm1);
    data->exist_nm1 = true;
  }
  *Un = Q;

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::computeBar(bool doInitialTasks, DistMacroCellSet *macroCells,
                                    DistGeoState &geoState, int scopeDepth)
{
                                                                                                                          
  if (doInitialTasks) {
    if (data->use_nm1){
      varFcn->conservativeToPrimitive(*Unm1, *Vn);
      domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 2);
      varFcn->primitiveToConservative(*VnBar, *Unm1Bar);
    }
    if (data->use_nm2){
      varFcn->conservativeToPrimitive(*Unm2, *Vn);
      domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 3);
      varFcn->primitiveToConservative(*VnBar, *Unm2Bar);
    }
  }
  else {
    if (data->use_nm2) *Unm2Bar = *Unm1Bar;
    if (data->use_nm1) *Unm1Bar = *UnBar;
  }
                                                                                                                          
  varFcn->conservativeToPrimitive(*Un, *Vn);
  domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 1);
  varFcn->primitiveToConservative(*VnBar, *UnBar);
                                                                                                                          
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::writeToDisk(char *name)
{

  domain->writeVectorToFile(name, 0, 0.0, *Un);
  if (data->use_nm1) 
    domain->writeVectorToFile(name, 1, 0.0, *Unm1);
  if (data->use_nm2) 
    domain->writeVectorToFile(name, 2, 0.0, *Unm2);

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::get_dW_dt(bool doInitialTasks,
                                   DistGeoState &geoState,
                                   DistVec<double> &ctrlVol,
                                   DistSVec<double,dim> &Q,
                                   DistSVec<double,dim> &VBar,
                                   DistSVec<double,dim> &dWdt,
                                   DistSVec<double,dim> &dWBardt,
                                   DistMacroCellSet *macroCells,
                                   DistSVec<double,1> **volRatio,
                                   int scopeDepth)
{
                                                                                                                          
  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(domain->getNodeDistInfo());
                                                                                                                          
  *Sigma = 0.0;
                                                                                                                          
  varFcn->primitiveToConservative(VBar, *QBar);
                                                                                                                          
  // compute dWdt //
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->get_dW_dt(Q.getMasterFlag(iSub), geoState(iSub),
                                   ctrlVol(iSub), Q(iSub), dWdt(iSub));
                                                                                                                          
  // compute dWBardt //
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->get_dWBar_dt((*QBar).getMasterFlag(iSub), geoState(iSub),
                                      ctrlVol(iSub), (*QBar)(iSub), (*UnBar)(iSub),
                                      (*Unm1Bar)(iSub), (*Unm2Bar)(iSub), (*Sigma)(iSub));
                                                                                                                          
  domain->assemble_dWdt(dWdt, *Sigma);
                                                                                                                          
  domain->computedWBar_dt(dWBardt, *Sigma, macroCells, volRatio, scopeDepth);
                                                                                                                          
                                                                                                                          
  delete (Sigma);
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
DistVec<double>* DistTimeState<dim>::getDerivativeOfInvReynolds(DistGeoState &geoState, DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                                          DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim> &dV, double dMach)
{

//Remark: Error mesage for pointers
  if (dIdti == 0) {
    fprintf(stderr, "*** Error: Varible dIdti does not exist!\n");
    exit(1);
  }
  if (dIdtv == 0) {
    fprintf(stderr, "*** Error: Varible dIdtv does not exist!\n");
    exit(1);
  }
  if (dIrey == 0) {
    fprintf(stderr, "*** Error: Varible dIrey does not exist!\n");
    exit(1);
  }

  //domain->computeDerivativeOfInvReynolds(fet, varFcn, geoState, X, dX, ctrlVol, dCtrlVol, V, dV, *idti, *dIdti, *idtv, *dIdtv, *dIrey, dMach, betav, beta, dbeta, k1, cmach);
  domain->computeDerivativeOfInvReynolds(fet, varFcn, geoState, X, dX, ctrlVol, dCtrlVol, V, dV, *idti, *dIdti, *idtv, *dIdtv, *dIrey, dMach, tprec, sprec);

  return dIrey;

}

//------------------------------------------------------------------------------
