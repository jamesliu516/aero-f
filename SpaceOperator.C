#include <SpaceOperator.h>

#include <BcFcn.h>
#include <BcDef.h>
#include <FluxFcn.h>
#include <RecFcnDesc.h>
#include <FemEquationTermDesc.h>
#include <VolumicForceTerm.h>
#include <DistVMSLESTerm.h>
#include <DistDynamicVMSTerm.h>
#include <DynamicLESTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DistDynamicLESTerm.h>
#include <DistNodalGrad.h>
#include <DistEdgeGrad.h>
#include <DistExactRiemannSolver.h>
#include <DistExtrapolation.h>
#include <DistBcData.h>
#include <DistMacroCell.h>
#include <Domain.h>
#include <DistGeoState.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <DistTimeState.h>
#include <Communicator.h>
#include <Timer.h>
#include <FluidSelector.h>


//------------------------------------------------------------------------------

template<int dim>
SpaceOperator<dim>::SpaceOperator(IoData &ioData, VarFcn *vf, DistBcData<dim> *bc,
				  DistGeoState *gs, Domain *dom, DistSVec<double,dim> *v)
  : varFcn(vf), bcData(bc), geoState(gs), domain(dom)
// Included (MB)
, iod(&ioData)
{

  locAlloc = true;
  rshift = dim;
  failsafe = ioData.ts.implicit.newton.failsafe;

  timer = domain->getTimer();
  com = domain->getCommunicator();

  if (v)
    V = v->alias();
  else
    V = new DistSVec<double,dim>(domain->getNodeDistInfo());

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_ || ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_) {
    dU = new DistSVec<double,dim>(domain->getNodeDistInfo());
    dV = new DistSVec<double,dim>(domain->getNodeDistInfo());
    dRm = new DistSVec<double,dim>(domain->getNodeDistInfo());

    *dU = 0.0;
    *dV = 0.0;
    *dRm = 0.0;
  }
  else {
    dU = 0;
    dV = 0;
    dRm = 0;
  }

  bcFcn = createBcFcn(ioData);

  fluxFcn = createFluxFcn(ioData);

  recFcn = createRecFcn(ioData);

  ngrad  = new DistNodalGrad<dim, double>(ioData, domain);

  egrad = 0;
  if (ioData.schemes.ns.dissipation == SchemeData::SIXTH_ORDER ||
      ioData.schemes.ns.gradient == SchemeData::NON_NODAL)
    egrad = new DistEdgeGrad<dim>(ioData, domain);
  xpol = 0;
  if (ioData.schemes.bc.type == BoundarySchemeData::CONSTANT_EXTRAPOLATION ||
      ioData.schemes.bc.type == BoundarySchemeData::LINEAR_EXTRAPOLATION)
    xpol = new DistExtrapolation<dim>(ioData, domain, varFcn);


  smag = 0;
  wale = 0;
  dles = 0;
  vms = 0;
  dvms = 0;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::LES) {
    if (ioData.eqs.tc.les.type == LESModelData::SMAGORINSKY)
      smag = new SmagorinskyLESTerm(ioData, varFcn);
    else if (ioData.eqs.tc.les.type == LESModelData::WALE)
       wale = new WaleLESTerm(ioData, varFcn);
    else if (ioData.eqs.tc.les.type == LESModelData::DYNAMIC){
      dles = new DistDynamicLESTerm<dim>(varFcn, ioData, domain);
    }
    else if (ioData.eqs.tc.les.type == LESModelData::VMS)
      vms = new DistVMSLESTerm<dim>(varFcn, ioData, domain);
    else if (ioData.eqs.tc.les.type == LESModelData::DYNAMICVMS){
      dvms = new DistDynamicVMSTerm<dim>(varFcn, ioData, domain);
    }
  }

  fet = createFemEquationTerm(ioData);
  volForce = createVolumicForceTerm(ioData);

  T = 0;
  if (fet && ioData.embed.viscousinterfaceorder == EmbeddedFramework::SECOND)
    T = new DistVec<double>(domain->getNodeDistInfo());

  if (ioData.problem.type[ProblemData::LINEARIZED])  {
    use_modal = true;
    if (ioData.linearizedData.domain == LinearizedData::FREQUENCY)  {
      use_complex = true;
      compNodalGrad = new DistNodalGrad<dim, bcomp>(ioData, domain);
    }
    else  {
      use_complex = false;
      compNodalGrad = 0;
    }
  }
  else {
    compNodalGrad = 0;
    use_modal = false;

// Included (MB)
    if (ioData.sa.comp3d == SensitivityAnalysis::OFF_COMPATIBLE3D)
      use_modal = true;
  }

  if (use_modal || ioData.ts.form == TsData::DESCRIPTOR)
    descriptorCase = DESCRIPTOR;
  else if (ioData.ts.form == TsData::HYBRID)
    descriptorCase = HYBRID;
  else
    descriptorCase = NONDESCRIPTOR;  


  if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT)
    order = 1;
  else
    order = 2;
}

//------------------------------------------------------------------------------

template<int dim>
SpaceOperator<dim>::SpaceOperator(const SpaceOperator<dim> &spo, bool typeAlloc)
{

  locAlloc = typeAlloc;

  varFcn = spo.varFcn;
  bcData = spo.bcData;
  geoState = spo.geoState;
  V = spo.V;

  bcFcn = spo.bcFcn;
  fluxFcn = spo.fluxFcn;
  recFcn = spo.recFcn;
  ngrad = spo.ngrad;
  compNodalGrad = spo.compNodalGrad;
  egrad = spo.egrad;
  xpol = spo.xpol;
  vms = spo.vms;
  smag = spo.smag;
  wale = spo.wale;
  dles = spo.dles;
  dvms = spo.dvms;
  fet = spo.fet;
  volForce = spo.volForce;

  failsafe = spo.failsafe;
  rshift = spo.rshift;

  domain = spo.domain;

  timer = spo.timer;
  com = spo.com;

  use_modal = spo.use_modal;
  descriptorCase = spo.descriptorCase;
// Included (MB)
  iod = spo.iod;

}

//------------------------------------------------------------------------------

template<int dim>
SpaceOperator<dim>::~SpaceOperator()
{

  if (locAlloc) {
    if (V) delete V;
    if (T) delete T;
    if (bcFcn) delete bcFcn;
    if (fluxFcn) {
        fluxFcn += BC_MIN_CODE;
	for(int i=0;i<BC_MAX_CODE - BC_MIN_CODE + 1;++i)
	  {
	    delete fluxFcn[i];
	  }
        delete [] fluxFcn;
    }
    if (recFcn) delete recFcn;
    if (ngrad) delete ngrad;
    if (compNodalGrad) delete compNodalGrad;
    if (egrad) delete egrad;
    if (xpol) delete xpol;
    if (vms) delete vms;
    if (smag) delete smag;
    if (wale) delete wale;
    if (dles) delete dles;
    if (dvms) delete dvms;
    if (fet) delete fet;
    if (volForce) delete volForce;
  }

}

//------------------------------------------------------------------------------

template<int dim>
BcFcn *SpaceOperator<dim>::createBcFcn(IoData &ioData)
{

  BcFcn *bf = 0;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
    if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	  ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
	bf = new BcFcnSA(ioData);
      else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
	bf = new BcFcnKE(ioData);
    }
    else {
      if (ioData.bc.wall.integration != BcsWallData::WALL_FUNCTION)
	bf = new BcFcnNS;
    }
  }

  return bf;

}

//------------------------------------------------------------------------------

template<int dim>
FluxFcn **SpaceOperator<dim>::createFluxFcn(IoData &ioData)
{

  FluxFcn **ff = 0;

  if (ioData.schemes.ns.advectiveOperator == SchemeData::FINITE_VOLUME)
    rshift = 0;

  double gamma = ioData.schemes.ns.gamma;

  ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
  ff -= BC_MIN_CODE;
  
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 12)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");
  ff[BC_SYMMETRY] = new FluxFcn(rshift,BC_SYMMETRY,ioData,varFcn); 
  ff[BC_OUTLET_MOVING] = new FluxFcn(rshift,BC_OUTLET_MOVING,ioData,varFcn);
  ff[BC_OUTLET_FIXED] = new FluxFcn(rshift,BC_OUTLET_FIXED,ioData,varFcn); 
  ff[BC_INLET_MOVING] = new FluxFcn(rshift,BC_INLET_MOVING,ioData,varFcn);
  ff[BC_INLET_FIXED] = new FluxFcn(rshift,BC_INLET_FIXED,ioData,varFcn);
  ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(rshift,BC_ADIABATIC_WALL_MOVING,ioData,varFcn);
  ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(rshift,BC_ADIABATIC_WALL_FIXED,ioData,varFcn);
  ff[BC_SLIP_WALL_MOVING] = new FluxFcn(rshift,BC_SLIP_WALL_MOVING,ioData,varFcn);
  ff[BC_SLIP_WALL_FIXED] = new FluxFcn(rshift,BC_SLIP_WALL_FIXED,ioData,varFcn);
  ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(rshift,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn);
  ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(rshift,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn);
  ff[BC_INTERNAL] = new FluxFcn(rshift,BC_INTERNAL,ioData,varFcn);

  return ff;

}

//------------------------------------------------------------------------------
// note: only one value for beta is allowed for k-e

template<int dim>
RecFcn *SpaceOperator<dim>::createRecFcn(IoData &ioData)
{

  RecFcn *rf = 0;

  double beta = ioData.schemes.ns.beta;
  double eps = ioData.schemes.ns.eps;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
      if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) {
	if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	  rf = new RecFcnConstant<6>;
      }
      else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
	if (ioData.schemes.ns.limiter == SchemeData::NONE) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLinearConstant<6>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnLinearVanAlbada<6>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnVanAlbadaConstant<6>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnVanAlbada<6>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLtdLinearConstant<6>(beta, eps);
	}
      }
    }
    else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
      if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) {
	if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	  rf = new RecFcnConstant<7>;
      }
      else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
	if (ioData.schemes.ns.limiter == SchemeData::NONE) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLinearConstant<7>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnLinearVanAlbada<7>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnVanAlbadaConstant<7>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnVanAlbada<7>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLtdLinearConstant<7>(beta, eps);
	}
      }
    }
  } else {
    if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT)
      rf = new RecFcnConstant<dim>;
    else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
      if (ioData.schemes.ns.limiter == SchemeData::NONE)
	rf = new RecFcnLinear<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA)
	rf = new RecFcnVanAlbada<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::EXTENDEDVANALBADA)
	rf = new RecFcnExtendedVanAlbada<dim>(beta, eps,ioData.eqs.fluidModelMap.dataMap[0]->pmin ,
                                              ioData.eqs.fluidModelMap.dataMap[0]->rhomin,
                                              ioData.schemes.ns.xip,
                                              ioData.schemes.ns.xirho );
      else if (ioData.schemes.ns.limiter == SchemeData::BARTH)
	rf = new RecFcnBarth<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VENKAT)
	rf = new RecFcnVenkat<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR)
	rf = new RecFcnLtdLinear<dim>(beta, eps);
    }
  }

  if (!rf) {
    com->fprintf(stderr, "*** Error: no valid choice for the reconstruction\n");
    exit(1);
  }

  return rf;

}
//------------------------------------------------------------------------------

template<int dim>
FemEquationTerm *SpaceOperator<dim>::createFemEquationTerm(IoData &ioData)
{

  FemEquationTerm *fem = 0;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
    if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)
	fem = new FemEquationTermSA(ioData, varFcn);
      else if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
	fem = new FemEquationTermDES(ioData, varFcn);
      else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
	fem = new FemEquationTermKE(ioData, varFcn);
    }
    else
      fem = new FemEquationTermNS(ioData, varFcn);
  }

  return fem;

}

//------------------------------------------------------------------------------

template<int dim>
VolumicForceTerm *SpaceOperator<dim>::createVolumicForceTerm(IoData &ioData)
{

  VolumicForceTerm *vft = 0;

  if(varFcn->gravity_value()>0.0)
    vft = new VolumicForceTerm(ioData);

  return vft;

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::setBcFcn(BcFcn *bf)
{

  bcFcn = bf;

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::setFluxFcn(FluxFcn **ff)
{

  fluxFcn = ff;

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::setRecFcn(RecFcn *rf)
{

  recFcn = rf;

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::setFemEquationTerm(FemEquationTerm *fem)
{

  fet = fem;

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::fix(DistSVec<bool,2>& tag)
{

  ngrad->fix(tag);

  if (egrad) egrad->fix(tag);

}


//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::resetTag()
{

  ngrad->resetTag();

  if (egrad) egrad->resetTag();

}

//------------------------------------------------------------------------------

// Modified (MB)
template<int dim>
void SpaceOperator<dim>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistSVec<double,dim> &R,
                                         DistTimeState<dim> *timeState, bool compatF3D)
{
  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    timer->addNodalGradTime(t0);
  }

  if (egrad) {
    egrad->compute(geoState->getConfig(), X);
	}

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms) {
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);
	}

  if (smag) {
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);
	}

  if (wale) {
     domain->computeWaleLESTerm(wale, X, *V, R);
	}

  if (dles) {
    dles->compute(ctrlVol, *bcData, X, *V, R);
	}

  DistVec<double> *irey;
  if(timeState)
    irey = timeState->getInvReynolds();
  else{
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  if (fet) {
    domain->computeGalerkinTerm(fet, *bcData, *geoState, X, *V, R);
    bcData->computeNodeValue(X);
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(ctrlVol, *irey, fluxFcn, recFcn, *bcData, *geoState,
                                  X, *V, *ngrad, egrad, R, failsafe, rshift);

// Included
  domain->getGradP(*ngrad);

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if(dvms)
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState,
                  timeState, X, U, *V, R, failsafe, rshift);

// Modified (MB)
  if (compatF3D) {
    if (descriptorCase != DESCRIPTOR)  {
      int numLocSub = R.numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        switch (descriptorCase) {
          case HYBRID:{
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invsqcv = 1.0 / sqrt(cv[i]);          
              for (int j=0; j<dim; ++j)
                r[i][j] *= invsqcv;
            }
            break; }
          case NONDESCRIPTOR: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invcv = 1.0 / cv[i];
              for (int j=0; j<dim; ++j)
                r[i][j] *= invcv;
            }
            break; }
        }
      }
    }

  }

  // Delete the pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}

template<int dim>
void SpaceOperator<dim>::computeResidualRestrict(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
		DistSVec<double,dim> &U,
		DistSVec<double,dim> &R,
		DistTimeState<dim> *timeState,
		RestrictionMapping<dim> & restrictionMapping, bool compatF3D)
{
	std::vector<std::vector<int> > sampledLocNodes =
		restrictionMapping.getRestrictedToOriginLocNode() ;

	R = 0.0;
	varFcn->conservativeToPrimitive(U, *V);

	// need everywhere
  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    timer->addNodalGradTime(t0); 
	}

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms)
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);

  if (wale)
     domain->computeWaleLESTerm(wale, X, *V, R);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R);

  DistVec<double> *irey;
  if(timeState)
    irey = timeState->getInvReynolds();
  else{
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  if (fet) {	// KTC: made efficient
    domain->computeGalerkinTerm(fet, *bcData, *geoState, X, *V, R);
    bcData->computeNodeValue(X);
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);	// MORE EFF??

	// KTC: made efficient
	domain->computeFiniteVolumeTerm(ctrlVol, *irey, fluxFcn, recFcn,
			*bcData, *geoState, X, *V, *ngrad, egrad, R, failsafe, rshift);

// Included
	// KTC: not needed
  //domain->getGradP(*ngrad);	 // not needed

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if(dvms)
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState,
                  timeState, X, U, *V, R, failsafe, rshift);

// Modified (MB)
  if (compatF3D) {
    if (descriptorCase != DESCRIPTOR)  {
      int i;
      int numLocSub = R.numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        switch (descriptorCase) {
          case HYBRID:{
            for (int iSampledNode=0; iSampledNode<sampledLocNodes[iSub].size(); ++iSampledNode) {
	      i = sampledLocNodes[iSub][iSampledNode];
              double invsqcv = 1.0 / sqrt(cv[i]);
              for (int j=0; j<dim; ++j) 
                r[i][j] *= invsqcv;
            }
            break; }
          case NONDESCRIPTOR: {
            for (int iSampledNode=0; iSampledNode<sampledLocNodes[iSub
].size(); ++iSampledNode) {
              i = sampledLocNodes[iSub][iSampledNode];
              double invcv = 1.0 / cv[i];
              for (int j=0; j<dim; ++j)
                r[i][j] *= invcv;
            }
            break; }
        }
      }
    }
  }

  // Delete the pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}
//------------------------------------------------------------------------------

// Modified (MB)
template<int dim>
void SpaceOperator<dim>::computeResidual(DistExactRiemannSolver<dim> *riemann,
                                         DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistSVec<double,dim> &R,
                                         DistTimeState<dim> *timeState, bool compatF3D)
{
  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    timer->addNodalGradTime(t0);
  }

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms)
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);

  if (wale)
     domain->computeWaleLESTerm(wale, X, *V, R);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R);

  DistVec<double> *irey;
  if(timeState)
    irey = timeState->getInvReynolds();
  else{
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  if (fet) {
    domain->computeGalerkinTerm(fet, *bcData, *geoState, X, *V, R);
    bcData->computeNodeValue(X);
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(*riemann, ctrlVol, *irey, fluxFcn, recFcn, *bcData, *geoState,
                                  X, *V, *ngrad, egrad, R, failsafe, rshift);

// Included
  domain->getGradP(*ngrad);

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if(dvms)
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState,
                  timeState, X, U, *V, R, failsafe, rshift);

// Modified (MB)
  if (compatF3D) {
    if (descriptorCase != DESCRIPTOR) {
      int numLocSub = R.numLocSub();
      int iSub;
#pragma omp parallel for
      for (iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        switch (descriptorCase) {
          case HYBRID:{
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invsqcv = 1.0 / sqrt(cv[i]);
              for (int j=0; j<dim; ++j) 
                r[i][j] *= invsqcv;
            }
            break; }
          case NONDESCRIPTOR: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invcv = 1.0 / cv[i];
              for (int j=0; j<dim; ++j)
                r[i][j] *= invcv;
            }
            break; } 
        }
      }
    }
  }

  // Delete pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::computeDerivativeOfResidual
(
  DistSVec<double,3> &X, DistSVec<double,3> &dX
  , DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol
  , DistSVec<double,dim> &U
  , double dMach
  , DistSVec<double,dim> &R, DistSVec<double,dim> &dR
  , DistTimeState<dim> *timeState
)
{

  dR = 0.0;

  varFcn->conservativeToPrimitive(U, *V);

//Remark: Error mesage for pointers
  if (dV == 0) {
    fprintf(stderr, "*** Error: Variable dV does not exist!\n");
    exit(1);
  }

  *dV = 0.0;

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  {
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    ngrad->computeDerivative(geoState->getConfigSA(), X, dX, ctrlVol, dCtrlVol, *V, *dV);
  }

  if (egrad) {
    egrad->compute(geoState->getConfig(), X);
    egrad->computeDerivative(geoState->getConfig(), X, dX);
  }

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
    xpol->computeDerivative(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms) {
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);
    vms->computeDerivative(geoState->getConfig(), ctrlVol, X, *V, R);
  }

  if (smag) {
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);
    domain->computeDerivativeOfSmagorinskyLESTerm(smag, X, *V, R);
  }

  if (dles){
    com->fprintf(stderr, "***** The equivalent derivatives of the functions dles->computeTestFilterValues and dles->computeTestFilterValues are not implemented!\n");
    exit(1);
  }

  DistVec<double> *irey;
  DistVec<double> *direy;
  if(timeState) 
  {
    irey = timeState->getInvReynolds();
    direy = timeState->getDerivativeOfInvReynolds(*geoState, X, dX, ctrlVol, dCtrlVol, *V, *dV, dMach);
  }
  else 
  {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    direy = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
    *direy = 0.0;
  }

  if (fet) 
  {
    domain->computeDerivativeOfGalerkinTerm(fet, *bcData, *geoState, X, dX, *V, *dV, dMach, dR);
    bcData->computeNodeValue(X);
    bcData->computeDerivativeOfNodeValue(X, dX);
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    ngrad->limitDerivative(recFcn, X, dX, ctrlVol, dCtrlVol, *V, *dV);
//    ngrad->limit(recFcn, X, ctrlVol, *V);
  }

  domain->computeDerivativeOfFiniteVolumeTerm
  (
    ctrlVol, dCtrlVol, *irey, *direy, fluxFcn, recFcn, *bcData, *geoState, 
    X, dX, *V, *dV, *ngrad, egrad, dMach, dR
  );

  domain->getGradP(*ngrad);
  domain->getDerivativeOfGradP(*ngrad);

  if (volForce) 
  {
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);
    domain->computeDerivativeOfVolumicForceTerm(volForce, ctrlVol, dCtrlVol, *V, *dV, dR);
  }

  if(dvms) {
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState, timeState, X, U, *V, R, failsafe, rshift);
    dvms->computeDerivative(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState, timeState, X, U, *V, R, failsafe, rshift);
  }

  if (descriptorCase != DESCRIPTOR)  {
    int numLocSub = dR.numLocSub();
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double *dcv = dCtrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      double (*dr)[dim] = dR.subData(iSub);
      double (*drm)[dim] = (*dRm).subData(iSub);
      switch (descriptorCase) {
        case HYBRID:{
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invsqcv = 1.0 / sqrt(cv[i]);
            double dSqrtInvcv = ( (-0.5) / pow(cv[i], 1.5) ) * dcv[i];
            for (int j=0; j<dim; ++j)
              dr[i][j] = ( ( dr[i][j] * invsqcv ) + ( r[i][j] * dSqrtInvcv ) );
          }
          break; }
        case NONDESCRIPTOR: {
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invcv = 1.0 / cv[i];
            double dInvcv = ( (-1.0) / ( cv[i] * cv[i] ) ) * dcv[i];
            for (int j=0; j<dim; ++j)
              dr[i][j] = ( ( dr[i][j] * invcv ) + ( r[i][j] * dInvcv ) );
          }
          break; }
      }
    }
  }

  // Delete pointers for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
    if (direy)
      delete direy;
  }
  irey = 0;
  direy = 0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::computeInviscidResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &U, DistSVec<double,dim> &R, DistTimeState<dim> *timeState, bool compatF3D)
{

  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms)
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R);

  DistVec<double> *irey;
  if(timeState)
    irey = timeState->getInvReynolds();
  else {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(ctrlVol, *irey, fluxFcn, recFcn, *bcData, *geoState,
                                  X, *V, *ngrad, egrad, R, failsafe, rshift);

  domain->getGradP(*ngrad);

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if(dvms)
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState,
                  timeState, X, U, *V, R, failsafe, rshift);

  if (compatF3D) {
    if (descriptorCase != DESCRIPTOR) {
      int numLocSub = R.numLocSub();
      int iSub;
#pragma omp parallel for
      for (iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        switch (descriptorCase) {
          case HYBRID: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invsqcv = 1.0 / sqrt(cv[i]);
              for (int j=0; j<dim; ++j) 
                r[i][j] *= invsqcv;
            }
            break; }
          case NONDESCRIPTOR: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invcv = 1.0 / cv[i];
              for (int j=0; j<dim; ++j) 
                r[i][j] *= invcv;
            }
            break; }
        }
      }
    }
  }

  // Delete pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }

  irey = 0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::computeViscousResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &U, DistSVec<double,dim> &R, DistTimeState<dim> *timeState, bool compatF3D)
{

  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol){
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);
  }

  if (vms)
    vms->compute(geoState->getConfig(), ctrlVol, X, *V, R);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R);

  if (fet) {
    domain->computeOnlyGalerkinTerm(fet, *bcData, *geoState, X, *V, R);
    bcData->computeNodeValue(X);
  }

  //new source term: need dVdXj (warning for jac if limited rec -> recompute gradients)
  //domain->computePointWiseSourceTerm(*geoState, ctrlVol, *ngrad, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if(dvms)
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState,
                  timeState, X, U, *V, R, failsafe, rshift);

  if (compatF3D) {
    if (descriptorCase != DESCRIPTOR)  {
      int numLocSub = R.numLocSub();
      int iSub;
#pragma omp parallel for
      for (iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        switch (descriptorCase) {
          case HYBRID: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invsqcv = 1.0 / sqrt(cv[i]);
              for (int j=0; j<dim; ++j) 
                r[i][j] *= invsqcv;
            }
            break; } 
          case NONDESCRIPTOR: {
            for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
              double invcv = 1.0 / cv[i];
              for (int j=0; j<dim; ++j)
                r[i][j] *= invcv;
            }
            break; }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistSVec<double,dim> &Wstarij,
                                         DistSVec<double,dim> &Wstarji, DistLevelSetStructure *distLSS,
                                         bool linRecAtInterface, bool viscSecOrder, DistVec<int> &fluidId, 
                                         DistSVec<double,dim> &R, DistExactRiemannSolver<dim> *riemann, 
                                         int Nriemann, DistSVec<double,3> *Nsbar, int it,
                                         DistVec<GhostPoint<dim>*> *ghostPoints)
{
  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V, &fluidId);  

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0){
    double t0 = timer->getTime();
    // compute gradient of V using Phi:
    // for node with Phi, gradient of V is computed using V-values of neighbours
    // that have the same Phi-sign
    bool linFSI = linRecAtInterface || viscSecOrder;
    ngrad->compute(geoState->getConfig(), X, ctrlVol, 
                   fluidId, *V, linFSI, distLSS);
    timer->addNodalGradTime(t0);
  }

  if (fet) {
    if (viscSecOrder) {
      varFcn->computeTemperature(*V,*T,&fluidId);
      ngrad->computeTemperatureGradient(geoState->getConfig(), X, ctrlVol, 
                     fluidId, *T, distLSS);
    }
    this->populateGhostPoints(ghostPoints,X,U,varFcn,distLSS,viscSecOrder,fluidId);
  }

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol) //boundary condition using xpol = extrapolation
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R, ghostPoints, distLSS);

  if (wale)
     domain->computeWaleLESTerm(wale, X, *V, R, ghostPoints, distLSS);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R, ghostPoints, distLSS);

  if (fet) {
      domain->computeGalerkinTerm(fet,*bcData,*geoState,X,*V,R,ghostPoints,distLSS);
      bcData->computeNodeValue(X);
  }

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(ctrlVol, *riemann, fluxFcn, recFcn, *bcData,
                                  *geoState, X, *V, Wstarij, Wstarji, distLSS, linRecAtInterface, fluidId, Nriemann,
                                  Nsbar, *ngrad, egrad, R, it, failsafe,rshift);
  if (descriptorCase != DESCRIPTOR)  {
    int numLocSub = R.numLocSub();
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      switch (descriptorCase) {
        case HYBRID: {
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invsqcv = 1.0 / sqrt(cv[i]);
            for (int j=0; j<dim; ++j) 
              r[i][j] *= invsqcv;
          }
          break; }
        case NONDESCRIPTOR: {
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invcv = 1.0 / cv[i];
            for (int j=0; j<dim; ++j) 
              r[i][j] *= invcv;
          }
          break; }

      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistSVec<double,dim> &Wstarij,
                                         DistSVec<double,dim> &Wstarji, 
										 DistVec<int> &countWstarij, DistVec<int> &countWstarji,
										 DistLevelSetStructure *distLSS, bool linRecAtInterface, bool viscSecOrder, 
										 DistVec<int> &fluidId, DistSVec<double,dim> &R, 
										 DistExactRiemannSolver<dim> *riemann, int Nriemann, 
										 DistSVec<double,3> *Nsbar, double dt, double alpha, 
										 int it, DistVec<GhostPoint<dim>*> *ghostPoints)
{
  R = 0.0;
  varFcn->conservativeToPrimitive(U, *V, &fluidId);  

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0){
    double t0 = timer->getTime();
    // compute gradient of V using Phi:
    // for node with Phi, gradient of V is computed using V-values of neighbours
    // that have the same Phi-sign
        bool linFSI = linRecAtInterface || viscSecOrder;
	ngrad->compute(geoState->getConfig(), X, ctrlVol, fluidId, *V, 
			  	   Wstarij, Wstarji, countWstarij, countWstarji, 
			  	   linFSI, distLSS);
    timer->addNodalGradTime(t0);
  }

  if (fet) {
    if (viscSecOrder) {
      varFcn->computeTemperature(*V,*T,&fluidId);
      ngrad->computeTemperatureGradient(geoState->getConfig(), X, ctrlVol, 
                     fluidId, *T, distLSS);
    }
    this->populateGhostPoints(ghostPoints,X,U,varFcn,distLSS,viscSecOrder,fluidId);
  }

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol) //boundary condition using xpol = extrapolation
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);

  if (smag)
    domain->computeSmagorinskyLESTerm(smag, X, *V, R, ghostPoints, distLSS);

  if (wale)
     domain->computeWaleLESTerm(wale, X, *V, R, ghostPoints, distLSS);

  if (dles)
    dles->compute(ctrlVol, *bcData, X, *V, R, ghostPoints, distLSS);

  if (fet) {
      domain->computeGalerkinTerm(fet,*bcData,*geoState,X,*V,R,ghostPoints,distLSS);
      bcData->computeNodeValue(X);
  }

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(ctrlVol, *riemann, fluxFcn, recFcn, *bcData,
                                  *geoState, X, *V, Wstarij, Wstarji, countWstarij, countWstarji,
								  distLSS, linRecAtInterface, fluidId, Nriemann, Nsbar, dt, alpha,
								  *ngrad, egrad, R, it, failsafe,rshift);
  if (descriptorCase != DESCRIPTOR)  {
    int numLocSub = R.numLocSub();
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      switch (descriptorCase) {
        case HYBRID: {
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invsqcv = 1.0 / sqrt(cv[i]);
            for (int j=0; j<dim; ++j) 
              r[i][j] *= invsqcv;
          }
          break; }
        case NONDESCRIPTOR: {
          for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
            double invcv = 1.0 / cv[i];
            for (int j=0; j<dim; ++j) 
              r[i][j] *= invcv;
          }
          break; }

      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computePostOpDVMS(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                           DistSVec<double,dim> &U, DistVec<double> *Cs,
                                           DistTimeState<dim> *timeState)
{

  DistSVec<double,dim> *R;
  R = new DistSVec<double,dim>(domain->getNodeDistInfo());
  *R = 0.0;
  DistVec<double> *irey = timeState->getInvReynolds();

  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  domain->computeGalerkinTerm(fet, *bcData, *geoState, X, *V, *R);
  bcData->computeNodeValue(X);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  domain->computeFiniteVolumeTerm(ctrlVol, *irey, fluxFcn, recFcn, *bcData, *geoState, X, *V,
                                  *ngrad, egrad, *R, failsafe, rshift);
  dvms->computeCs(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol,
                  *bcData, *geoState, timeState, X, U, *V, *R, Cs, failsafe, rshift);

  if(R) delete R;
  irey = 0;

}

//------------------------------------------------------------------------------

template<int dim>
double SpaceOperator<dim>::recomputeResidual(DistSVec<double,dim> &F, DistSVec<double,dim> &Finlet)
{
 if(xpol){
   return domain->recomputeResidual(F, Finlet);
 }
 return 0.0;

}

//-----------------------------------------------------------------------------

template<int dim>
double SpaceOperator<dim>::computeRealFluidResidual(DistSVec<double, dim> &F, DistSVec<double,dim> &Freal,
                                                    DistLevelSetStructure &dLSS)
{ return domain->computeRealFluidResidual(F, Freal, dLSS); }

//------------------------------------------------------------------------------

// TODO(jontg): Remove X
// TODO(jontg): Fix computeRiemannWeightsForEmbeddedStruct
template<int dim>
void SpaceOperator<dim>::
updateSweptNodes(DistSVec<double,3> &X, int &phaseChangeChoice, int &phaseChangeAlg, 
                 DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                 DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                 DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                 DistLevelSetStructure *distLSS, double *vfar, DistVec<int> *fluidId)
{
  int iSub, numLocSub = this->domain->getNumLocSub();
  DistVec<int> init(domain->getNodeDistInfo()),next_init(domain->getNodeDistInfo());
  SubDomain **subD = this->domain->getSubDomain();

  varFcn->conservativeToPrimitive(U, V, fluidId);
  Weights = 0.0; VWeights = 0.0;

#pragma omp parallel for
  for(iSub=0;iSub<numLocSub;++iSub)
    for(int i=0;i<init(iSub).size();++i)
        init(iSub)[i] = ((*distLSS)(iSub).isSwept(0.0,i) || !(*distLSS)(iSub).isActive(0.0,i) ? 0 : 1);
  next_init = init;

  int iter=0, maxIter=1000000;
  int finished = 0;
  while(finished == 0 && iter < maxIter){++iter;finished = 1;
    switch(phaseChangeChoice){
    case 0: 
	  switch (phaseChangeAlg) {
		case 0:
		    domain->computeWeightsForEmbeddedStruct(X, V, Weights, VWeights, init, next_init, distLSS);
		  break;
		case 1:
		    domain->computeWeightsLeastSquaresForEmbeddedStruct(X, V, Weights, VWeights, init, next_init, distLSS);
//		    domain->computeWeightsForEmbeddedStruct(X, V, Weights, VWeights, init, next_init, distLSS);
		  break;
	  }
      break;
    case 1: com->fprintf(stderr," *** computeReimannWeights temporarily broken / using regular extrapolation ***\n");
      domain->computeWeightsForEmbeddedStruct(X, V, Weights, VWeights, init, next_init, distLSS);
      break;
    }

#pragma omp parallel for
    for(iSub=0;iSub<numLocSub;++iSub) {
	  int* locToGlobNodeMap = subD[iSub]->getNodeMap();
      for(int i=0;i<init(iSub).size();++i) {
//		if ((*distLSS)(iSub).isSwept(0.0,i)&&(*distLSS)(iSub).isActive(0.0,i))
//		  fprintf(stdout,"XY: Node %d needs phase change update:\n",locToGlobNodeMap[i]+1);
        if(init(iSub)[i]<1 && next_init(iSub)[i]>0) {
          const double one_over_weight=(double)1.0/Weights(iSub)[i];
		  switch (phaseChangeAlg) {
			case 0:
          	  for (int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight;
			  break;
			case 1:
			  if (Weights(iSub)[i]>0.0)
          	    for (int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight;
			  else if (Weights(iSub)[i]<0.0)
				for (int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d];
			  else
				fprintf(stderr,"Error: LS phase change update failed at node %d.\n", locToGlobNodeMap[i]+1);
			  break;
		  }
        }
	  }
    }

#pragma omp parallel for
    for(iSub=0;iSub<numLocSub;++iSub)
      for(int i=0;i<init(iSub).size();++i)
        if(init(iSub)[i]<1 && next_init(iSub)[i]>0) finished = 0;
    Weights = 0.0; VWeights = 0.0;
    init = next_init;
    com->globalOp(1,&finished,MPI_PROD);
  }
  if(iter == maxIter) {
    fprintf(stderr,"ERROR: abnormal termination in SpaceOperator::updateSweptNodes. Abort...\n");
    exit(-1);
  }

#pragma omp parallel for
  for(iSub=0;iSub<numLocSub;++iSub){
    for(int i=0;i<init(iSub).size();++i)
      if(init(iSub)[i] < 1 || !(*distLSS)(iSub).isActive(0.0,i))
        for(int d=0; d<dim; d++) V(iSub)[i][d] = vfar[d];
  }
  varFcn->primitiveToConservative(V, U, fluidId);
}

//------------------------------------------------------------------------------

template<int dim> 
void SpaceOperator<dim>::populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints, DistSVec<double,3> &X, DistSVec<double,dim> &U, VarFcn *varFcn,DistLevelSetStructure *distLSS,bool linFSI,DistVec<int> &tag)
{
ghostPoints->deletePointers();
domain->populateGhostPoints(ghostPoints,X,U,ngrad,varFcn,distLSS,linFSI,tag);
}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeCellAveragedStructNormal(DistSVec<double, 3> &Nsbar, DistLevelSetStructure *distLSS)
{
  DistVec<double> weights(domain->getNodeDistInfo()); //currently it's an integer, but might become "double" in future
  Nsbar   = 0.0;
  weights = 0.0;
  domain->computeCellAveragedStructNormal(Nsbar, weights, distLSS);
}

//-----------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
					 DistSVec<double,dim> &U, DistMat<Scalar,neq> &A,
					 DistTimeState<dim> *timeState)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

  A = 0.0;
  DistVec<double> *irey;
  if(timeState) {
    irey = timeState->getInvReynolds();
  }
  else {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  switch (descriptorCase) {
    case DESCRIPTOR: {
      DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
      unitCtrlVol = 1.0;
      if (fet)
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, unitCtrlVol, *V, A);

      domain->computeJacobianFiniteVolumeTerm(fluxFcn, *bcData, *geoState, *irey, X, unitCtrlVol, *V, A);

      if (volForce)
        domain->computeJacobianVolumicForceTerm(volForce, unitCtrlVol, *V, A);
      break; }
    case HYBRID: {
      DistVec<double> sqrtCtrlVol(domain->getNodeDistInfo());
      sqrtCtrlVol.pow(ctrlVol,0.5);
      if (fet)
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, sqrtCtrlVol, *V, A);

      domain->computeJacobianFiniteVolumeTerm(fluxFcn, *bcData, *geoState, *irey, X, sqrtCtrlVol, *V, A);

      if (volForce)
        domain->computeJacobianVolumicForceTerm(volForce, sqrtCtrlVol, *V, A);
  
      break; }
    case NONDESCRIPTOR: {
  {
      if (fet)
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);
  }
      domain->computeJacobianFiniteVolumeTerm(fluxFcn, *bcData, *geoState, *irey, X, ctrlVol, *V, A);

      if (volForce)
        domain->computeJacobianVolumicForceTerm(volForce, ctrlVol, *V, A);
      break; }
  }

  // Delete pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}
//-----------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeJacobian(DistExactRiemannSolver<dim> *riemann, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistMat<Scalar,neq> &A,
                                         DistTimeState<dim> *timeState)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

  A = 0.0;
  DistVec<double> *irey;
  if(timeState) {
    irey = timeState->getInvReynolds();
  }
  else {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  switch (descriptorCase) {
    case DESCRIPTOR: {
      DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
      unitCtrlVol = 1.0;
      if (fet)
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, unitCtrlVol, *V, A);

      domain->computeJacobianFiniteVolumeTerm(*riemann, fluxFcn, *bcData, *geoState, *irey, X, unitCtrlVol, *V, A);

      if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, unitCtrlVol, *V, A);
      break; }
    case HYBRID: {
      DistVec<double> sqrtCtrlVol(domain->getNodeDistInfo());
      sqrtCtrlVol.pow(ctrlVol,0.5);
      if (fet)
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, sqrtCtrlVol, *V, A);

      domain->computeJacobianFiniteVolumeTerm(*riemann, fluxFcn, *bcData, *geoState, *irey, X, sqrtCtrlVol, *V, A);

      if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, sqrtCtrlVol, *V, A);
      break; }
    case NONDESCRIPTOR: {

    if (fet)
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);

    domain->computeJacobianFiniteVolumeTerm(*riemann, fluxFcn, *bcData, *geoState, *irey, X, ctrlVol, *V, A);

    if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, ctrlVol, *V, A);
    break; }

  }

  // Delete pointer for consistency
  if (timeState == 0)
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar,int neq>
void SpaceOperator<dim>::computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U,
                                         DistLevelSetStructure *distLSS,
                                         DistVec<int> &fluidId, 
                                         DistExactRiemannSolver<dim> *riemann, 
                                         int Nriemann, DistSVec<double,3> *Nsbar,
                                         DistVec<GhostPoint<dim>*> *ghostPoints,
                                         DistMat<Scalar,neq>& A,
                                         DistTimeState<dim>* timeState)
{
  A = 0.0;
  varFcn->conservativeToPrimitive(U, *V, &fluidId);  
  
  DistVec<double> *irey;
  if(timeState) {
    irey = timeState->getInvReynolds();
  }
  else {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
  }

  if (fet) {
    domain->computeJacobianGalerkinTerm(fet,*bcData,*geoState,X,ctrlVol, *V,A,ghostPoints,distLSS);
  }
  //if (fet)
  //  domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);
  
  domain->computeJacobianFiniteVolumeTerm(ctrlVol, *riemann, fluxFcn, *bcData, *geoState,
                                          X, *V, distLSS, fluidId, Nriemann, Nsbar, A,*irey);

  if (volForce)
    domain->computeJacobianVolumicForceTerm(volForce, ctrlVol, *V, A);
  
  // Delete pointer for consistency
  if (timeState == 0) 
  {
    if (irey)
      delete irey;
  }
  irey = 0;
}
  

//------------------------------------------------------------------------------
template<int dim>
void SpaceOperator<dim>::recomputeRHS(DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistSVec<double,dim> &rhs)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

  if(xpol)
    domain->recomputeRHS(varFcn, *V, rhs, xpol, *bcData, *geoState, X);

}

//------------------------------------------------------------------------------
template<int dim>
void SpaceOperator<dim>::recomputeRHS(DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> &fluidId, DistSVec<double,dim> &rhs)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V, &fluidId);
#endif

  if(xpol)
    domain->recomputeRHS(varFcn, *V, fluidId, rhs, xpol, *bcData, *geoState, X);

}


//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeViscousJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistMat<Scalar,neq> &A)  {

  A = 0.0;
  if (fet)  {
    switch (descriptorCase) {
      case DESCRIPTOR: {
        DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
        unitCtrlVol = 1.0;
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, unitCtrlVol, *V, A);
        domain->finishJacobianGalerkinTerm(unitCtrlVol, A);
        break; }
      case HYBRID: {
        DistVec<double> sqrtCtrlVol(domain->getNodeDistInfo());
        sqrtCtrlVol.pow(ctrlVol,0.5);;
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, sqrtCtrlVol, *V, A);
        domain->finishJacobianGalerkinTerm(sqrtCtrlVol, A);
        break; }
      case NONDESCRIPTOR: {
        domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);
        domain->finishJacobianGalerkinTerm(ctrlVol, A);
        break; } 
   }

// Included (MB*)
    if ((iod->eqs.type == EquationsData::NAVIER_STOKES) && (iod->eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
      if ((iod->bc.wall.integration == BcsWallData::WALL_FUNCTION) && (iod->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
        domain->computeBCsJacobianWallValues(fet, *bcData, *geoState, X, *V);
        bcData->computeNodeWallValues(X);
      }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::getExtrapolationValue(DistSVec<double,dim> &U,
                                               DistSVec<double,dim> &Ubc, DistSVec<double,3>& X)
{

  if(xpol){
    DistSVec<double,dim>* VV = new DistSVec<double,dim>(domain->getNodeDistInfo());
    varFcn->conservativeToPrimitive(U, *VV); //assumption : only one phase at far-field boundary
    domain->getExtrapolationValue(xpol, *VV, Ubc, varFcn, *bcData, *geoState, X);
    delete VV;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::applyExtrapolationToSolutionVector(DistSVec<double,dim> &U,
                                                            DistSVec<double,dim> &Ubc)
{
  if(xpol){
     domain->applyExtrapolationToSolutionVector(xpol, U, Ubc);
  }
}
//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::applyBCsToSolutionVector(DistSVec<double,dim> &U, DistLevelSetStructure *distLSS)
{

  if (bcFcn)
    domain->applyBCsToSolutionVector(bcFcn, *bcData, U, distLSS);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::applyBCsToResidual(DistSVec<double,dim> &U, DistSVec<double,dim> &R, DistLevelSetStructure *distLSS)
{

  if (bcFcn)
    domain->applyBCsToResidual(bcFcn, *bcData, U, R, distLSS);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::applyBCsToDerivativeOfResidual(DistSVec<double,dim> &U, DistSVec<double,dim> &dR)
{

//Remark: Error mesage for pointers
  if (dU == 0) {
    fprintf(stderr, "*** Error: Variable dU does not exist!\n");
    exit(1);
  }

  *dU = 0.0;

  if (bcFcn)
    domain->applyBCsToDerivativeOfResidual(bcFcn, *bcData, U, *dU, dR);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::applyBCsToJacobian(DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

  if (bcFcn)
    domain->applyBCsToJacobian(bcFcn, *bcData, U, A);

// Included (MB*)
  if (bcFcn)
    if ((iod->eqs.type == EquationsData::NAVIER_STOKES) && (iod->eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
      if ((iod->bc.wall.integration == BcsWallData::WALL_FUNCTION) && (iod->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS))
        domain->applyBCsToJacobianWallValues(bcFcn, *bcData, U, A);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::applyBCsToH2Jacobian(DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

  if (bcFcn)
    domain->applyBCsToH2Jacobian(bcFcn, *bcData, U, A);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
template<class Scalar>
void SpaceOperator<dim>::applyBCsToH2Jacobian(DistSVec<double,dim> &U, DistMat<Scalar,dim> &A)
{

  if (bcFcn)
    domain->applyBCsToH2Jacobian(bcFcn, *bcData, U, A);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void SpaceOperator<dim>::computeH1(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                   DistSVec<double,dim> &U, DistMat<Scalar,dim> &H1)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

  H1 = 0.0;
  
  switch (descriptorCase) {
    case DESCRIPTOR: {
      DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
      unitCtrlVol = 1.0;
      domain->computeH1(fluxFcn, *bcData, *geoState, unitCtrlVol, *V, H1);
      break; }
    case HYBRID: {
      DistVec<double> sqrtCtrlVol(domain->getNodeDistInfo());
      sqrtCtrlVol.pow(ctrlVol,0.5);
      domain->computeH1(fluxFcn, *bcData, *geoState, sqrtCtrlVol, *V, H1);
      break; }
    case NONDESCRIPTOR: {
      domain->computeH1(fluxFcn, *bcData, *geoState, ctrlVol, *V, H1);
      break; }
  } 
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeH2(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
				   DistSVec<double,dim> &U, DistMat<Scalar,neq> &H2,
				   DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
				   DistSVec<double,dim> &bij, DistSVec<double,dim> &bji)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    ngrad->limit(recFcn, X, ctrlVol, *V);
  }
#endif

  H2 = 0.0;

  domain->computeH2(fluxFcn, recFcn, *bcData, *geoState, X, *V, *ngrad,
		    H2, aij, aji, bij, bji);

}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar1, class Scalar2>
void SpaceOperator<dim>::applyH2(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                         DistSVec<double,dim> &U, DistMat<Scalar1,dim> &H2,
                         DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
                         DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
                         DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod)
{

  int numLocSub = p.numLocSub();

  DistSVec<Scalar2, dim> V2(domain->getNodeDistInfo());
  DistNodalGrad<dim, Scalar2> *distNodalGrad = getDistNodalGrad(p);
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locU)[dim] = U.subData(iSub);
    double (*locV)[dim] = V->subData(iSub);
    Scalar2 (*locp)[dim] = p.subData(iSub);
    Scalar2 (*locV2)[dim] = V2.subData(iSub);
    for (int i=0; i<p.subSize(iSub); ++i) {
      varFcn->conservativeToPrimitive(locU[i], locV[i]);
      varFcn->multiplyBydVdU(locV[i], locp[i], locV2[i]);
    }
  }

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    distNodalGrad->compute(geoState->getConfig(), X, ctrlVol, V2);
    distNodalGrad->limit(recFcn, X, ctrlVol, V2);
  }

  switch (descriptorCase) {
    case DESCRIPTOR: {
      DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
      unitCtrlVol = 1.0;
      domain->computeMatVecProdH2(recFcn, X, unitCtrlVol, H2, aij, aji, bij, bji, V2, *distNodalGrad, prod);
      break; }
    case HYBRID: {
      DistVec<double> sqrtCtrlVol(domain->getNodeDistInfo());
      sqrtCtrlVol.pow(ctrlVol,0.5);
      domain->computeMatVecProdH2(recFcn, X, sqrtCtrlVol, H2, aij, aji, bij, bji, V2, *distNodalGrad, prod);
      break; }
    case NONDESCRIPTOR: {
      domain->computeMatVecProdH2(recFcn, X, ctrlVol, H2, aij, aji, bij, bji, V2, *distNodalGrad, prod);
      break; }
  }

// Included (MB*)
  if (bcFcn)
    if ((iod->eqs.type == EquationsData::NAVIER_STOKES) && (iod->eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
      if ((iod->bc.wall.integration == BcsWallData::WALL_FUNCTION) && (iod->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS))
        domain->applyBCsToProduct(bcFcn, *bcData, U, prod);

}

//------------------------------------------------------------------------------

// V is destroyed and replaced by dVdU*p
// dVdxj is destroyed and replaced by the derivatives of V = dVdU*p

template<int dim>
template<class Scalar1, class Scalar2>
void SpaceOperator<dim>::applyH2T(DistSVec<double,3> &X,
                         DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
                         DistMat<Scalar1,dim> &H2, DistSVec<double,dim> &aij,
                         DistSVec<double,dim> &aji, DistSVec<double,dim> &bij,
                         DistSVec<double,dim> &bji, DistSVec<Scalar2,dim> &p,
                         DistSVec<Scalar2,dim> &prod)  {

  DistSVec<Scalar2,dim> prod2(domain->getEdgeDistInfo());
  DistSVec<Scalar2,dim> prod3(domain->getEdgeDistInfo());
  DistSVec<Scalar2,dim> prod4(domain->getEdgeDistInfo());

  varFcn->conservativeToPrimitive(U, *V);

  int numLocSub = p.numLocSub();
  DistNodalGrad<dim, Scalar2> *distNodalGrad = getDistNodalGrad(p);

  domain->computeMatVecProdH2T(recFcn, X, ctrlVol, H2, aij, aji, bij, bji, p,  prod2, prod, prod3, prod4);
  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    distNodalGrad->computeT(geoState->getConfig(), X, ctrlVol, prod, prod3, prod4);
  domain->computeMatVecProdH2Tb(recFcn, X, ctrlVol, H2, *distNodalGrad, p, prod, prod2);
  

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    double (*locV)[dim] = V->subData(iSub);
    Scalar2 (*locp)[dim] = prod.subData(iSub);

    // prod corresponds to zu in the Fortran code
    for (int i=0; i<p.subSize(iSub); ++i)
      varFcn->multiplyBydVdUT(locV[i], locp[i], locp[i]);

  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::rstFluxFcn(IoData &ioData)
{

  FluxFcn **ff = createFluxFcn(ioData);

  setFluxFcn(ff);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::printAllMatrix(DistMat<Scalar,neq> &A, int it)
{

  domain->printAllMatrix(A, it);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::printAllVariable(DistSVec<double,3> &X, DistSVec<double,dim> &U, int it){

  //DistSVec<double,dim> *V = new DistSVec<double,dim>(domain->getNodeDistInfo());
  //varFcn->conservativeToPrimitive(U,*V);
  //domain->printAllVariable(X,U,it);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::printVariable(DistSVec<double,dim> &U){

  domain->printVariable(U, varFcn);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeGradP(DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &U)
{

  varFcn->conservativeToPrimitive(U, *V);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  {
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    timer->addNodalGradTime(t0);
    ngrad->limit(recFcn, X, ctrlVol, *V);
  }

  domain->getGradP(*ngrad);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeDerivativeOfGradP
(
  DistSVec<double,3> &X, DistSVec<double,3> &dX,
  DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
  DistSVec<double,dim> &U, DistSVec<double,dim> &dU
)
{

  varFcn->conservativeToPrimitive(U, *V);
  varFcn->conservativeToPrimitiveDerivative(U, dU, *V, *dV);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  {
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    ngrad->computeDerivative(geoState->getConfigSA(), X, dX, ctrlVol, dCtrlVol, *V, *dV);
    ngrad->limitDerivative(recFcn, X, dX, ctrlVol, dCtrlVol, *V, *dV);
  }

  domain->getDerivativeOfGradP(*ngrad);

}

//------------------------------------------------------------------------------

// UH (08/10) The following function is never called.
//template<int dim>
//void SpaceOperator<dim>::computeDerivativeOfGradP
//(
//  DistSVec<double,3> &X, DistSVec<double,3> &dX,
//  DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
//  DistSVec<double,dim> &U
//)
//{
//
//  varFcn->conservativeToPrimitive(U, *V);
//  *dV = 0.0;
//
//  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  {
//    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
//    ngrad->computeDerivative(geoState->getConfigSA(), X, dX, ctrlVol, dCtrlVol, *V, *dV);
//    ngrad->limitDerivative(recFcn, X, dX, ctrlVol, dCtrlVol, *V, *dV);
//  }
//
//  domain->getDerivativeOfGradP(*ngrad);
//
//}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
					  DistVec<double> &ctrlVol, double (*Fs)[3], int sizeFs, 
					  DistLevelSetStructure *distLSS,
					  DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
					  DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn,
					  DistVec<int>* fid)
{
  double pinternal = iod->aero.pressure;
  
  switch (forceApp)
    {

    case 1: // Control Volume Boundaries
      if(orderOfAccuracy>1) ngrad->compute(geoState->getConfig(), X, ctrlVol,distLSS->getStatus(),*V,true,distLSS);
      domain->computeCVBasedForceLoad(forceApp, orderOfAccuracy, *geoState, X, Fs,
                                      sizeFs, distLSS, pinternal, Wstarij, Wstarji,
                                      *V,ghostPoints,postFcn,ngrad,varFcn,fid);
      break;

    case 2: // Embedded Surface
      domain->computeEmbSurfBasedForceLoad(*iod,forceApp,orderOfAccuracy,X,Fs,sizeFs,
					   distLSS,pinternal,Wstarij,Wstarji,*V,ghostPoints,postFcn,ngrad,varFcn,fid);
      break;

    case 3: // Reconstructed Surface
      domain->computeRecSurfBasedForceLoad(forceApp,orderOfAccuracy,X,Fs,sizeFs,
					   distLSS,pinternal,Wstarij,Wstarji,*V,ghostPoints,postFcn,varFcn,fid);
      break;
    default:
      fprintf(stderr,"ERROR: force approach not specified correctly (%d)! Abort...\n", forceApp); 
      exit(-1);
    }
}

//------------------------------------------------------------------------------
//             MULTIPHASE SPACE OPERATOR
//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhaseSpaceOperator<dim,dimLS>::MultiPhaseSpaceOperator(IoData &ioData, VarFcn *vf, DistBcData<dim> *bc,
				  DistGeoState *gs, Domain *dom, DistSVec<double,dim> *v)
  : SpaceOperator<dim>(ioData, vf, bc, gs, dom, v)
{

  recFcnLS = createRecFcnLS(ioData);
  ngradLS = new DistNodalGrad<dimLS, double>(ioData, this->domain, 1);

  if (ioData.mf.interfaceTreatment == 
      MultiFluidData::SECONDORDER) {

    for (int l = 0; l < 2; ++l) {

      higherOrderData.vals[l] = new DistSVec<double,dim>(dom->getNodeDistInfo());
      higherOrderData.cutgrads[l] = new DistNodalGrad<dim,double>(ioData, this->domain, 1);
      higherOrderData.counts[l] = new DistVec<int>(dom->getNodeDistInfo());
    }
  }

  if (ioData.mf.riemannNormal == MultiFluidData::MESH) {
    for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) {
      dom->getSubDomain()[iSub]->getEdges().setMultiFluidRiemannNormal(EdgeSet::MF_RIEMANN_NORMAL_MESH);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhaseSpaceOperator<dim,dimLS>::MultiPhaseSpaceOperator(const MultiPhaseSpaceOperator<dim,dimLS> &spo, bool typeAlloc)
{
  this->locAlloc = typeAlloc;

  this->varFcn = spo.varFcn;
  this->bcData = spo.bcData;
  this->geoState = spo.geoState;
  this->V = spo.V;

  this->bcFcn = spo.bcFcn;
  this->fluxFcn = spo.fluxFcn;
  this->recFcn = spo.recFcn;
  this->ngrad = spo.ngrad;
  this->compNodalGrad = spo.compNodalGrad;
  this->egrad = spo.egrad;
  this->xpol = spo.xpol;
  this->vms = spo.vms;
  this->smag = spo.smag;
  this->wale = spo.wale;
  this->dles = spo.dles;
  this->dvms = spo.dvms;
  this->fet = spo.fet;
  this->volForce = spo.volForce;

  this->failsafe = spo.failsafe;
  this->rshift = spo.rshift;

  this->domain = spo.domain;

  this->timer = spo.timer;
  this->com = spo.com;

  this->use_modal = spo.use_modal;
  this->descriptorCase = spo.descriptorCase;
// Included (MB)
  this->iod = spo.iod;

  recFcnLS = spo.recFcnLS;
  ngradLS  = spo.ngradLS;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhaseSpaceOperator<dim,dimLS>::~MultiPhaseSpaceOperator()
{

  if (this->locAlloc) {
    delete recFcnLS;
    delete ngradLS;
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
RecFcn *MultiPhaseSpaceOperator<dim,dimLS>::createRecFcnLS(IoData &ioData)
{
  RecFcn *rf = 0;

  double beta = ioData.schemes.ls.beta;
  double eps = ioData.schemes.ls.eps;

  if (ioData.schemes.ls.reconstruction == SchemeData::CONSTANT)
    rf = new RecFcnConstant<dimLS>;
  else if (ioData.schemes.ls.reconstruction == SchemeData::LINEAR) {
    if (ioData.schemes.ls.limiter == SchemeData::NONE)
      rf = new RecFcnLinear<dimLS>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VANALBADA)
      rf = new RecFcnVanAlbada<dimLS>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::BARTH)
      rf = new RecFcnBarth<dimLS>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VENKAT)
      rf = new RecFcnVenkat<dimLS>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::P_SENSOR)
      rf = new RecFcnLtdLinear<dimLS>(beta, eps);
  }

  return rf;
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistSVec<double,dimLS> &Phi,
                                         FluidSelector &fluidSelector, DistSVec<double,dim> &R,
                                         DistExactRiemannSolver<dim> *riemann, 
                                         DistTimeState<dim> * timeState,
                                         int it)
{

  R = 0.0;
  this->varFcn->conservativeToPrimitive(U, *(this->V), fluidSelector.fluidId);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0){
    double t0 = this->timer->getTime();
    this->ngrad->compute(this->geoState->getConfig(), X, ctrlVol, *fluidSelector.fluidId, *(this->V));
    this->timer->addNodalGradTime(t0);
  }
  
/*  this->domain->getCommunicator()->fprintf(stderr,"%lf\n",this->V->norm());  
  int numLocSub = X.numLocSub();
  double gradsize = 0.0;
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      NodalGrad<dim,double>& grad = this->ngrad->operator()(iSub);
      for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
        if (fabs(grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3]) > 1000.0) {
          std::cout << "Error: grad size is giant: " << grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3] <<
                   " " << this->domain->getSubDomain()[iSub]->getNodeMap()[i] + 1;
          if (this->domain->getSubDomain()[iSub]->higherOrderMF->isCellCut(i))
            std::cout << " [cell is cut]";
          std::cout << std::endl;
        }
        gradsize += grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3];
      }
    }

    this->domain->getCommunicator()->globalSum(1,&gradsize);
    this->domain->getCommunicator()->fprintf(stderr,"gradsize = %lf\n",gradsize);
 */

  if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0){
    double t0 = this->timer->getTime();
    ngradLS->compute(this->geoState->getConfig(), X, ctrlVol, Phi);
    this->timer->addNodalGradTime(t0);
  }

  if (this->egrad)
    this->egrad->compute(this->geoState->getConfig(), X);

  if (this->xpol)
    this->xpol->compute(this->geoState->getConfig(),this->geoState->getInletNodeNorm(), X);

  if (this->fet) {
    this->domain->computeGalerkinTerm(this->fet, *(this->bcData), *(this->geoState), X, *(this->V), R);
    if (!this->dvms) this->bcData->computeNodeValue(X);
  }

  if (this->volForce)
    this->domain->computeVolumicForceTerm(this->volForce, ctrlVol, *(this->V), R);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0)
    this->ngrad->limit(this->recFcn, X, ctrlVol, *(this->V),timeState->getFirstOrderNodeSet());
  //if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0)  
  //  ngradLS->limit(recFcnLS, X, ctrlVol, PhiS);

  this->domain->computeFiniteVolumeTerm(ctrlVol, *riemann, this->fluxFcn, this->recFcn, *(this->bcData),
                                  *(this->geoState), X, *(this->V), fluidSelector, *(this->ngrad), this->egrad,
					Phi, *ngradLS, R, it, this->failsafe,this->rshift);

  if (this->descriptorCase != this->DESCRIPTOR)  {
    int numLocSub = R.numLocSub();
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      if (this->descriptorCase == this->HYBRID) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invsqcv = 1.0 / sqrt(cv[i]);
          for (int j=0; j<dim; ++j) 
            r[i][j] *= invsqcv;
        }
      }
      else if (this->descriptorCase == this->NONDESCRIPTOR) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invcv = 1.0 / cv[i];
          for (int j=0; j<dim; ++j)
            r[i][j] *= invcv;
        }
      } 
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::computeResidualLS(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
							   DistSVec<double,dimLS> &Phi, DistVec<int> &fluidId, DistSVec<double,dim> &U,
							   DistSVec<double,dimLS> &PhiF, DistLevelSetStructure *distLSS,
							   bool linRecAtFSInterface,int method)
{
  PhiF = 0.0;

  this->varFcn->conservativeToPrimitive(U, *(this->V), &fluidId);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0){
    double t0 = this->timer->getTime();
    this->ngrad->compute(this->geoState->getConfig(), X, ctrlVol, fluidId, *(this->V), true, distLSS);
    this->timer->addLSNodalWeightsAndGradTime(t0);
  }

  if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0) {
    if(distLSS && !((*distLSS)(0).withCracking()))
      ngradLS->compute(this->geoState->getConfig(), X, ctrlVol, distLSS->getStatus(), Phi, linRecAtFSInterface);
    else
      ngradLS->compute(this->geoState->getConfig(), X, ctrlVol, Phi);
  }

  if (this->egrad)
    this->egrad->compute(this->geoState->getConfig(), X);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0)
    this->ngrad->limit(this->recFcn, X, ctrlVol, *(this->V));

  if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0)
    ngradLS->limit(recFcnLS, X, ctrlVol, Phi);

  this->domain->computeFiniteVolumeTermLS(this->fluxFcn, this->recFcn, recFcnLS, *(this->bcData), *(this->geoState), X, *(this->V),fluidId,
                                    *(this->ngrad), *ngradLS, this->egrad, Phi, PhiF, distLSS);
//  this->domain->getCommunicator()->fprintf(stderr,"PhiF res: %lf\n", PhiF.norm());

  if (this->descriptorCase != this->DESCRIPTOR)  {
    int numLocSub = PhiF.numLocSub();
    int iSub;
    double gradsize = 0.0;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dimLS] = PhiF.subData(iSub);
      NodalGrad<dim,double>& grad = this->ngrad->operator()(iSub);
      if (this->descriptorCase == this->HYBRID) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double sqcv = sqrt(cv[i]);
          double invsqcv = 1.0 / sqcv;
          for (int idim=0; idim<dimLS; idim++) 
            r[i][idim] *= invsqcv;
          if (method == 1) { // DJA: applicable to both cases?
            for (int idim=0; idim<dimLS; idim++)
              r[i][idim] -= sqcv*(grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3])*Phi(iSub)[i][idim];
          }
        }
      }
      else if (this->descriptorCase == this->NONDESCRIPTOR) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invcv = 1.0 / cv[i];
          for (int idim=0; idim<dimLS; idim++) 
            r[i][idim] *= invcv;
          if (method == 1) { 
            for (int idim=0; idim<dimLS; idim++)
              r[i][idim] -= (grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3])*Phi(iSub)[i][idim];
          gradsize += grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3];
          }
        }
      }

//    this->domain->getCommunicator()->globalSum(1,&gradsize);
//    this->domain->getCommunicator()->fprintf(stderr,"gradsize = %lf\n",gradsize);
    }
  }
/*  this->domain->getCommunicator()->fprintf(stderr,"PhiF res: %lf\n", PhiF.norm());
  this->domain->getCommunicator()->fprintf(stderr,"V res: %lf\n", this->V->norm());
  this->domain->getCommunicator()->fprintf(stderr,"Phi res: %lf\n", Phi.norm());
*/
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &U, 
                     DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                     DistLevelSetStructure *distLSS, bool linRecAtInterface, bool viscSecOrder, DistExactRiemannSolver<dim> *riemann, int Nriemann, DistSVec<double,3> *Nsbar,
                     DistSVec<double,dimLS> &PhiV, FluidSelector &fluidSelector, DistSVec<double,dim> &R, int it,
                     DistVec<GhostPoint<dim>*> *ghostPoints)
{
  R = 0.0;
  int numLocSub = R.numLocSub();
  this->varFcn->conservativeToPrimitive(U, *(this->V), fluidSelector.fluidId);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0){
    bool linFSI = linRecAtInterface || viscSecOrder;
    double t0 = this->timer->getTime();
    this->ngrad->compute(this->geoState->getConfig(), X, ctrlVol, *fluidSelector.fluidId, *(this->V), linFSI, distLSS);
    this->timer->addNodalGradTime(t0);
  }

  if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0){
    double t0 = this->timer->getTime();
    ngradLS->compute(this->geoState->getConfig(), X, ctrlVol, distLSS->getStatus(), PhiV, 1/* need to compute grad near FS Interface*/);
    //fluid Id (status) is used only to distinguish different materials. By using distLSS->getStatus(), we allow crossing FF interface but not FS interface.
    //  One can alternatively try to use fluidSelector.fluidId, which avoids crossing both FF and FS interfaces.
    this->timer->addNodalGradTime(t0);
/*
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
     
      SVec<double,dimLS>* grads[3] = {&(ngradLS->getX()(iSub)),
                                      &(ngradLS->getY()(iSub)),
                                      &(ngradLS->getZ()(iSub))}; 
      for (int i = 0; i < normals[0]->subSize(iSub); ++i) {
        for (int k = 0; k < dimLS; ++k) {
          double mag = grads[0]->v[i][k]*grads[0]->v[i][k] + 
                       grads[1]->v[i][k]*grads[1]->v[i][k] +
                       grads[2]->v[i][k]*grads[2]->v[i][k];
          mag = std::max<double>(1.0e-10,sqrt(mag));
          (*normals[0])(iSub)[i][k] =  grads[0]->v[i][k]/mag;
          (*normals[1])(iSub)[i][k] =  grads[1]->v[i][k]/mag;
          (*normals[2])(iSub)[i][k] =  grads[2]->v[i][k]/mag;
        } 
      }
    }

    ngradLS_second[0]->compute(this->geoState->getConfig(), X, ctrlVol, distLSS->getStatus(), 
                               *normals[0]);
    ngradLS_second[1]->compute(this->geoState->getConfig(), X, ctrlVol, distLSS->getStatus(), 
                               *normals[1]);
    ngradLS_second[2]->compute(this->geoState->getConfig(), X, ctrlVol, distLSS->getStatus(), 
                               *normals[2]);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      for (int i = 0; i < normals[0]->subSize(iSub); ++i) {

        for (int k = 0; k < dimLS; ++k) {
          (*curvature)(iSub)[i][k] = ngradLS_second[0]->getX()(iSub)[i][k] +
                                     ngradLS_second[1]->getY()(iSub)[i][k] +
                                     ngradLS_second[2]->getZ()(iSub)[i][k];
        }
      }
    }
  */  
  }

  if (this->fet) {
    if (viscSecOrder) {
      this->varFcn->computeTemperature(*(this->V),*(this->T),fluidSelector.fluidId);
      this->ngrad->computeTemperatureGradient(this->geoState->getConfig(), X, ctrlVol, 
                     *fluidSelector.fluidId, *(this->T), distLSS);
    }
    this->populateGhostPoints(ghostPoints,X,U,this->varFcn,distLSS,viscSecOrder,*fluidSelector.fluidId);
  }

  if (this->egrad)
    this->egrad->compute(this->geoState->getConfig(), X);

  if (this->xpol)
    this->xpol->compute(this->geoState->getConfig(),this->geoState->getInletNodeNorm(), X);

  if (this->fet) {
    this->domain->computeGalerkinTerm(this->fet,*(this->bcData),*(this->geoState),X,*(this->V),R,ghostPoints,distLSS);
    this->bcData->computeNodeValue(X);
  }

  if (this->volForce)
    this->domain->computeVolumicForceTerm(this->volForce, ctrlVol, *(this->V), R);

  if (dynamic_cast<RecFcnConstant<dim> *>(this->recFcn) == 0)
    this->ngrad->limit(this->recFcn, X, ctrlVol, *(this->V));

//  if (dynamic_cast<RecFcnConstant<dimLS> *>(recFcnLS) == 0)
//    ngradLS->limit(recFcnLS, X, ctrlVol, PhiV);

  //Now compute the FV fluxes!
  this->domain->computeFiniteVolumeTerm(ctrlVol, *riemann, this->fluxFcn, this->recFcn, *(this->bcData),
                                  *(this->geoState), X, *(this->V), Wstarij, Wstarji, distLSS, linRecAtInterface, fluidSelector, 
                                  Nriemann, Nsbar, *(this->ngrad), this->egrad,
                                  *ngradLS, R, it, this->failsafe,this->rshift);

  if (this->descriptorCase != this->DESCRIPTOR)  {
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      if (this->descriptorCase == this->HYBRID) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invsqcv = 1.0 / sqrt(cv[i]);
          for (int j=0; j<dim; ++j) 
            r[i][j] *= invsqcv;
        }
      }
      else if (this->descriptorCase == this->NONDESCRIPTOR) {
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invcv = 1.0 / cv[i];
          for (int j=0; j<dim; ++j) 
            r[i][j] *= invcv;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
template<class Scalar, int neq>
void MultiPhaseSpaceOperator<dim,dimLS>::computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistMat<Scalar,neq> &A,
                                         FluidSelector &fluidSelector, DistExactRiemannSolver<dim> *riemann,
                                         DistTimeState<dim>* timeState)
{

  //fprintf(stdout, "going through computeJacobian for two-phase flows\n");
#ifdef DOUBLE_CHECK
  this->varFcn->conservativeToPrimitive(U, *(this->V), fluidSelector.fluidId);
#endif

  A = 0.0;
  
  DistVec<double> *irey;
  if(timeState) {
    irey = timeState->getInvReynolds();
  }
  else {
    irey = new DistVec<double>(this->domain->getNodeDistInfo());
    *irey = 0.0;
  }


  if (this->descriptorCase == this->DESCRIPTOR || this->descriptorCase == this->HYBRID)  {
    fprintf(stderr, "**Error: no modal or descriptor/hybrid form for multiphase flows.. Exiting\n");
    exit(1);
  }
  else  {
    if (this->fet)
      this->domain->computeJacobianGalerkinTerm(this->fet, *(this->bcData), *(this->geoState), X, ctrlVol, *(this->V), A);
    this->domain->computeJacobianFiniteVolumeTerm(*riemann, this->fluxFcn, *(this->bcData), *(this->geoState), *(this->ngrad), *ngradLS, X, ctrlVol, *(this->V), A, fluidSelector);
    if (this->volForce)
      this->domain->computeJacobianVolumicForceTerm(this->volForce, ctrlVol, *(this->V), A);
  }
  
  // Delete pointer for consistency
  if (timeState == 0) 
  {
    if (irey)
      delete irey;
  }
  irey = 0;
}

template<int dim,int dimLS>  
template<class Scalar, int neq>
void MultiPhaseSpaceOperator<dim,dimLS>::computeJacobian(DistExactRiemannSolver<dim>* riemann,
                                                         DistSVec<double,3>& X, DistSVec<double,dim>& U,DistVec<double>& ctrlVol,
                                                         DistLevelSetStructure *distLSS,
                                                         int Nriemann, DistSVec<double,3>* Nsbar,
                                                         FluidSelector &fluidSelector,
                                                         DistMat<Scalar,neq>& A,DistTimeState<dim>* timeState) {

#ifdef DOUBLE_CHECK
  this->varFcn->conservativeToPrimitive(U, *(this->V), fluidSelector.fluidId);
#endif

  A = 0.0;
  
  DistVec<double> *irey;
  if(timeState) {
    irey = timeState->getInvReynolds();
  }
  else {
    irey = new DistVec<double>(this->domain->getNodeDistInfo());
    *irey = 0.0;
  }


  if (this->descriptorCase == this->DESCRIPTOR || this->descriptorCase == this->HYBRID)  {
    fprintf(stderr, "**Error: no modal or descriptor form for multiphase flows.. Exiting\n");
    exit(1);
  }
  else  {
    if (this->fet)
      this->domain->computeJacobianGalerkinTerm(this->fet, *(this->bcData), *(this->geoState), X, ctrlVol, *(this->V), A);
    this->domain->computeJacobianFiniteVolumeTerm(*riemann, this->fluxFcn, *(this->bcData), *(this->geoState),X,*(this->V), ctrlVol,
                                                  *ngradLS, distLSS, Nriemann, Nsbar, fluidSelector,A);
    if (this->volForce)
      this->domain->computeJacobianVolumicForceTerm(this->volForce, ctrlVol, *(this->V), A);
  }
  
  // Delete pointer for consistency
  if (timeState == 0) 
  {
    if (irey)
      delete irey;
  }
  irey = 0;

}


//------------------------------------------------------------------------------

template <int dim,int dimLS>
template<class Scalar>
void MultiPhaseSpaceOperator<dim,dimLS>::computeJacobianLS(DistSVec<double,3> &X,DistSVec<double,dim> &V, DistVec<double> &ctrlVol,
							   DistSVec<double,dimLS> &Phi,DistMat<Scalar,dimLS> &A,DistVec<int> &fluidId,DistLevelSetStructure* distLSS,
							   int method)
{
  A = 0.0;
  this->domain->computeJacobianFiniteVolumeTermLS(this->recFcn, recFcnLS,*(this->geoState),X,V,*(this->ngrad), *ngradLS,this->egrad,ctrlVol, Phi,A,distLSS);

  if (this->descriptorCase == this->NONDESCRIPTOR)  {
    int numLocSub = Phi.numLocSub();
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      NodalGrad<dim,double>& grad = this->ngrad->operator()(iSub);
      for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
	if (method == 1) {
	  for (int idim=0; idim<dimLS; idim++)
	    A(iSub).getElem_ii(i)[idim] -= (grad.getX()[i][1]+grad.getY()[i][2]+grad.getZ()[i][3]);
	}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV)
{
  DistVec<int> &fsId(distLSS->getStatus());

// set PhiV = 0 for 1) isolated regions, and 2) nodes under phase-change, because they will not be used.
// this might be redundant. but it's safer to do it.
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    LevelSetStructure& LSS((*distLSS)(iSub));
    SVec<double,dimLS>& subPhiV(PhiV(iSub));
    Vec<int>& subId(fsId(iSub));

    for(int i=0; i<subPhiV.size(); i++)
      if(subId[i]!=0 || LSS.isSwept(0.0,i))
        for(int k=0; k<dimLS; k++)
          subPhiV[i][k] = 0.0;
  }

  this->domain->extrapolatePhiV(distLSS,PhiV);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::extrapolatePhiV2(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV)
{
  this->domain->extrapolatePhiV(distLSS,PhiV);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::updateSweptNodes(DistSVec<double,3> &X, int &phaseChangeChoice,
                           DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                           DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                           DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
                           DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                           DistLevelSetStructure *distLSS, double *vfar, bool updateWithCracking,
                           DistVec<int> *fluidId0, DistVec<int> *fluidId)
{
  int iSub, numLocSub = this->domain->getNumLocSub();
  DistVec<int> init(this->domain->getNodeDistInfo()),next_init(this->domain->getNodeDistInfo());
  SubDomain **subD = this->domain->getSubDomain();

  this->varFcn->conservativeToPrimitive(U, V, fluidId);
  Weights = 0.0; VWeights = 0.0; PhiWeights = 0.0;

#pragma omp parallel for
  for(iSub=0;iSub<numLocSub;++iSub)
    for(int i=0;i<init(iSub).size();++i)
      init(iSub)[i] = (*distLSS)(iSub).isSwept(0.0,i) || (*fluidId)(iSub)[i]==(*distLSS)(iSub).numOfFluids()  ? 0 : 1;

  next_init = init;

  int iter=0, maxIter=1000000;
  int finished = 0;
  while(finished == 0 && iter < maxIter){++iter;finished = 1;
    switch(phaseChangeChoice){
    case 0: this->domain->computeWeightsForEmbeddedStruct(X, V, Weights, VWeights, Phi, PhiWeights,
                                                          init, next_init, distLSS, fluidId);
      break;
    case 1: this->com->fprintf(stderr," *** computeReimannWeights temporarily broken / using regular extrapolation ***\n");
      this->domain->computeWeightsForEmbeddedStruct(X, V, Weights, VWeights, Phi, PhiWeights,
                                                    init, next_init, distLSS, fluidId);
      break;
    }

    if(!updateWithCracking){
#pragma omp parallel for
      for(iSub=0;iSub<numLocSub;++iSub) {
        for(int i=0;i<init(iSub).size();++i)
          if(init(iSub)[i]<1.0 && next_init(iSub)[i]>0.0 || (init(iSub)[i]<1.0 && (*fluidId)(iSub)[i]==(*distLSS)(iSub).numOfFluids())) {
            if(!(*distLSS)(iSub).isActive(0.0,i)) {
              if(Weights(iSub)[i]>=0.1/*i.e. at least 1*/) {
                const double one_over_weight=(double)1.0/Weights(iSub)[i];
                V(iSub)[i][0] = vfar[0];
                for(int d=1; d<4; d++) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight; //pull the velocity from neighbors
                for(int d=4; d<dim; d++) V(iSub)[i][d] = vfar[d];
                for(int d=0;d<dimLS;++d) Phi(iSub)[i][d] = PhiWeights(iSub)[i][d]*one_over_weight;
              } else {
                for(int d=0; d<dim; d++) V(iSub)[i][d] = vfar[d];
                for(int d=0;d<dimLS;++d) Phi(iSub)[i][d] = -1.0; //not really needed.
              }
            } else {
              const double one_over_weight=(double)1.0/Weights(iSub)[i];
              for(int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight;
              for(int d=0;d<dimLS;++d) Phi(iSub)[i][d] = PhiWeights(iSub)[i][d]*one_over_weight;
            }
          }
      }
    } else {
#pragma omp parallel for
      for(iSub=0;iSub<numLocSub;++iSub) {

        int* locToGlobNodeMap = subD[iSub]->getNodeMap();

        for(int i=0;i<init(iSub).size();++i) {
          if(!((*distLSS)(iSub).isSwept(0.0,i)))
            continue;
          if((init(iSub)[i]<1.0 && next_init(iSub)[i]>0.0) || (init(iSub)[i]<1.0 && (*fluidId)(iSub)[i]==(*distLSS)(iSub).numOfFluids())) {
            if((*fluidId)(iSub)[i]==(*distLSS)(iSub).numOfFluids()) {
              if(!(*distLSS)(iSub).isOccluded(0.0,i)) {fprintf(stderr,"BUG!\n");exit(-1);} //just debug
              if(Weights(iSub)[i]>=0.1/*i.e. at least 1*/) {
                const double one_over_weight=(double)1.0/Weights(iSub)[i];
                V(iSub)[i][0] = vfar[0];
                for(int d=1; d<4; d++) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight; //pull the velocity from neighbors
                for(int d=4; d<dim; d++) V(iSub)[i][d] = vfar[d];
                for(int d=0;d<dimLS;++d) Phi(iSub)[i][d] = PhiWeights(iSub)[i][d]*one_over_weight;
              } else { // might want to print a warning message. I don't see why this can happen
                for(int d=0; d<dim; d++) V(iSub)[i][d] = vfar[d];
                for(int d=0;d<dimLS;++d) Phi(iSub)[i][d] = 0.0; //set phi = 0 because this node is at interface
              }
            } else {
              const double one_over_weight=(double)1.0/Weights(iSub)[i];
              double phiS = (*distLSS)(iSub).distToInterface(0.0,i); //this is the UNSIGNED distance
              if(phiS<0.0) {
                Phi(iSub)[i][dimLS-1] = std::fabs(PhiWeights(iSub)[i][dimLS-1]*one_over_weight);
                for(int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight;
                if((*fluidId)(iSub)[i]==0) Phi(iSub)[i][dimLS-1] *= -1.0;
                for(int d=0;d<dimLS-1;++d) Phi(iSub)[i][d] = PhiWeights(iSub)[i][d]*one_over_weight;
              } else {
                for(int d=0;d<dim;++d) V(iSub)[i][d] = VWeights(iSub)[i][d]*one_over_weight;
                Phi(iSub)[i][dimLS-1] = phiS;
                if((*fluidId)(iSub)[i]==0) Phi(iSub)[i][dimLS-1] *= -1.0;
                for(int d=0;d<dimLS-1;++d) Phi(iSub)[i][d] = PhiWeights(iSub)[i][d]*one_over_weight;
              }

//              if(Phi(iSub)[i][0]<0) {
//                fprintf(stderr,"Warning: got a swept node (%d) far from the interface! PhiS = %e. PhiF = %e. Id0 = %d, Id = %d.\n", 
//                        locToGlobNodeMap[i]+1, Phi(iSub)[i][0], PhiWeights(iSub)[i][0]*one_over_weight, 
//                        (*fluidId0)(iSub)[i], (*fluidId)(iSub)[i]);
//                Phi(iSub)[i][0] = std::fabs(PhiWeights(iSub)[i][0]*one_over_weight);
//              }
            }
          }
        }
      }
    }

#pragma omp parallel for
    for(iSub=0;iSub<numLocSub;++iSub)
      for(int i=0;i<init(iSub).size();++i)
        if(init(iSub)[i]<1.0 && next_init(iSub)[i]>0.0) finished = 0;
    Weights = 0.0; VWeights = 0.0; PhiWeights = 0.0;
    init = next_init;
    this->com->globalOp(1,&finished,MPI_PROD);
  }
  if(iter == maxIter) {
    fprintf(stderr,"ERROR: abnormal termination in SpaceOperator::updateSweptNodes. Abort...\n");
    exit(-1);
  }

  this->varFcn->primitiveToConservative(V, U, fluidId);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::resetFirstLayerLevelSetFS(DistSVec<double,dimLS> &PhiV, DistLevelSetStructure *distLSS,
                                           DistVec<int> &fluidId, DistSVec<bool,2> &Tag)
{
  /* -------------------------------------------------------------
      Reset PhiV for the first layer of nodes near the interface
      Rule #1. If it is occluded, set PhiV to 0. (temporarily removed rule 1, seems to be causing issues)
      Rule #2. If it is near the FS interface but not FF interface, set PhiV to its dist to wall.
      Rule #3. If it is near both FS and FF interface, set PhiV to min(PhiV, d2wall).
      Rule #4. Otherwise, do nothing.
     -------------------------------------------------------------- */

  this->domain->TagInterfaceNodes(dimLS-1,Tag,PhiV,distLSS);

  SubDomain **subD = this->domain->getSubDomain();
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    int* locToGlobNodeMap = subD[iSub]->getNodeMap();
    LevelSetStructure& LSS((*distLSS)(iSub));
    SVec<bool,2> &tag(Tag(iSub));
    Vec<int> &id(fluidId(iSub));
    SVec<double,dimLS> &phiv(PhiV(iSub));

    for(int i=0; i<tag.size(); i++) {
      // Rule #1.
      if(LSS.isOccluded(0.0,i)) {
        //phiv[i][dimLS-1] = 0.0;
        if(id[i]!=LSS.numOfFluids()) {//just for debug
          fprintf(stderr,"BUG: Node %d is occluded but its status is %d. numOfFluids = %d.\n", locToGlobNodeMap[i]+1, id[i], LSS.numOfFluids());
          exit(-1);
        }
        continue;
      }
      // Rule #2.
      if(tag[i][0] && !tag[i][1]) {
        double dist = LSS.distToInterface(0.0,i);
        if(dist<0) {//just for debug
          fprintf(stderr,"BUG: Node %d is near FS interface but its wall distance (%e) is invalid.\n", locToGlobNodeMap[i]+1, dist);
          exit(-1);
        }
        phiv[i][dimLS-1] = (phiv[i][dimLS-1]>0.0) ? dist : -1.0*dist;
        continue;
      }
      // Rule #3.
      if(tag[i][0] && tag[i][1]) {
        double dist = LSS.distToInterface(0.0,i);
        if(dist<0) {//just for debug
          fprintf(stderr,"BUG: Node %d is near FS interface but its wall distance (%e) is invalid.\n", locToGlobNodeMap[i]+1, dist);
          exit(-1);
        }
        dist = min(dist, fabs(phiv[i][dimLS-1]));
        phiv[i][dimLS-1] = (phiv[i][dimLS-1]>0.0) ? dist : -1.0*dist;
        continue;
      }
    }
  }
}

//-----------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhaseSpaceOperator<dim,dimLS>::findCutCells(DistSVec<double,dimLS>& phi,
						      DistVec<int>& status,
						      DistVec<int>& fluidId,
						      DistSVec<double,dim> &U,
						      DistSVec<double,3> &X) {

  SubDomain **subD = this->domain->getSubDomain();

  int iSub;

  this->varFcn->conservativeToPrimitive(U, *(this->V), &fluidId);
   
  for (int k = 0; k < dimLS; ++k) {
      
    status = 0;

#pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->findCutCells(k, phi(iSub), status,X(iSub));
      subD[iSub]->sndData(*this->domain->getLevelPat(), reinterpret_cast<int (*)[1]>(status.subData(iSub)));
    }

    this->domain->getLevelPat()->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->maxRcvData(*this->domain->getLevelPat(), reinterpret_cast<int (*)[1]>(status.subData(iSub)));
      subD[iSub]->template setCutCellFlags<dim>(k,status(iSub));
    }
  }

  /*int numCutCells = 0, cutCellsLoc,tmp;
  int* locoffset = new int[numLocSub];
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    
    locoffset[iSub] = numCutCells;
    numCutCells += subD[iSub]->getNumCutCells();
  }

  MPI_Allreduce(&numCutCells,&tmp, 1, MPI_INT,
                MPI_SUM, this->domain->com->comm);

  numCutCells = tmp;

  fprintf(stdout,"# of cut cells = %i\n",numCutCells);
  
  MPI_Scan(&numCutCells,&cutCellsLoc,1, MPI_INT, 
           MPI_SUM, this->domain->com->comm);
  
  double* grads = new double[numCutCells*3*dim*2];
  double* vals = new double[numCutCells*2*dim];
  int* gl_nodes = new int[numCutCells];
  
  memset(grads,0,sizeof(double)*numCutCells*3*dim*2);
  memset(vals,0,sizeof(double)*numCutCells*2*dim);
  memset(gl_nodes,0,sizeof(int)*numCutCells);
  */

  //std::cout << "Hello!!" << std::endl;
  //exit(-1);
#pragma omp parallel for
  for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
  
    SVec<double,dim>* v[2] = { &higherOrderData.vals[0]->operator()(iSub),
			       &higherOrderData.vals[1]->operator()(iSub)};
    
    NodalGrad<dim,double>* g[2] = { &higherOrderData.cutgrads[0]->operator()(iSub),
				    &higherOrderData.cutgrads[1]->operator()(iSub)};
    
    Vec<int>* c[2] = { &higherOrderData.counts[0]->operator()(iSub),
		       &higherOrderData.counts[1]->operator()(iSub)};

    for (int i = 0; i < 2; ++i) {
      *v[i] = 0.0;
      g[i]->getX() = 0.0;
      g[i]->getY() = 0.0;
      g[i]->getZ() = 0.0;
      *c[i] = 0;
    }
    
    /*subD[iSub]->collectCutCellData(grads + (locOffset[iSub]+cutCellsLoc)*3*2*dim,
					vals + (locOffset[iSub]+cutCellsLoc)*2*dim,
					gl_nodes + (locOffset[iSub]+cutCellsLoc));
    */
    
    subD[iSub]->collectCutCellData(v,g,c,
					(*(this->ngrad))(iSub), fluidId(iSub),
					(*(this->V))(iSub),X(iSub)) ;

  }

  for (int i = 0; i < 2; ++i) {

#pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->sndData(*this->domain->getLevelPat(), reinterpret_cast<int (*)[1]>(higherOrderData.counts[i]->subData(iSub)));
    }
    
    this->domain->getLevelPat()->exchange();
    
#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->addRcvData(*this->domain->getLevelPat(), reinterpret_cast<int (*)[1]>(higherOrderData.counts[i]->subData(iSub)));
    }
  }

  for (int i = 0; i < 2; ++i) {

#pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->sndData(*this->domain->getVecPat(), reinterpret_cast<double (*)[dim]>(higherOrderData.vals[i]->subData(iSub)));
    }
    
    this->domain->getVecPat()->exchange();
    
#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->addRcvData(*this->domain->getVecPat(), reinterpret_cast<double (*)[dim]>(higherOrderData.vals[i]->subData(iSub)));
    }
  }

  for (int i = 0; i < 2; ++i) {

#pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->sndData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getX().v);
    }
    
    this->domain->getVecPat()->exchange();
    
#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->addRcvData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getX().v);
    }

 #pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->sndData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getY().v);
    }
    
    this->domain->getVecPat()->exchange();
    
#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->addRcvData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getY().v);
    }

   #pragma omp parallel for
    for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
      
      subD[iSub]->sndData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getZ().v);
    }
    
    this->domain->getVecPat()->exchange();
    
#pragma omp parallel for
    for (iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {
      subD[iSub]->addRcvData(*this->domain->getVecPat(), higherOrderData.cutgrads[i]->operator()(iSub).getZ().v);
    }
  }

#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {

    SVec<double,dim>* v[2] = { &higherOrderData.vals[0]->operator()(iSub),
			       &higherOrderData.vals[1]->operator()(iSub)};
    
    NodalGrad<dim,double>* g[2] = { &higherOrderData.cutgrads[0]->operator()(iSub),
				    &higherOrderData.cutgrads[1]->operator()(iSub)};
    
    Vec<int>* c[2] = { &higherOrderData.counts[0]->operator()(iSub),
		       &higherOrderData.counts[1]->operator()(iSub)};
    
    subD[iSub]->storeCutCellData(v,g,c,
				 fluidId(iSub)
				 /*(*(this->V))(iSub)*/);
  } 

}

