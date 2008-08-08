#include <SpaceOperator.h>

#include <BcFcn.h>
#include <BcDef.h>
#include <FluxFcnDescWaterCompressible.h>
#include <FluxFcnDescPerfectGas.h>
#include <FluxFcnDescJWL.h>
#include <FluxFcnDescGasInGas.h>
#include <FluxFcnDescJWLInGas.h>
#include <FluxFcnDescLiquidInLiquid.h>
#include <FluxFcnDescGasInLiquid.h>
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
  if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {
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

// Operators for Level-Set Equation
  recFcnLS = createRecFcnLS(ioData);
  ngradLS = new DistNodalGrad<1, double>(ioData, domain, 1);

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
  recFcnLS = spo.recFcnLS;
  ngrad = spo.ngrad;
  ngradLS = spo.ngradLS;
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

// Included (MB)
  iod = spo.iod;

}

//------------------------------------------------------------------------------

template<int dim>
SpaceOperator<dim>::~SpaceOperator()
{

  if (locAlloc) {
    if (V) delete V;
    if (bcFcn) delete bcFcn;
    if (fluxFcn) {
        fluxFcn += BC_MIN_CODE;
        delete [] fluxFcn;
    }
    if (recFcn) delete recFcn;
    if (recFcnLS) delete recFcnLS;
    if (ngrad) delete ngrad;
    if (ngradLS) delete ngradLS;
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
	bf = new BcFcnKE;
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
  int prec;

  if (ioData.problem.prec == ProblemData::PRECONDITIONED) prec = 1;
  else prec = 0;

  double betaRef = ioData.prec.mach;
  double K1 = ioData.prec.k;
  double cmach = ioData.prec.cmach;
  double shockreducer = ioData.prec.shockreducer;

  //for GAS
  if (ioData.eqs.numPhase == 1){
    if (ioData.eqs.fluidModel.fluid == FluidModelData::GAS){
      if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
          ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
        if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS || 
            ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
          ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
          ff -= BC_MIN_CODE;
          if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
            ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasOutflowSA3D(ioData);
            ff[BC_OUTLET_FIXED] = new FluxFcnPerfectGasOutflowSA3D(ioData);
          }
          else {
            ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasInternalOutflowSA3D(ioData);
            ff[BC_OUTLET_FIXED] = new FluxFcnPerfectGasInternalOutflowSA3D(ioData);
          }
          if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
            ff[BC_INLET_MOVING] = new FluxFcnPerfectGasOutflowSA3D(ioData);
            ff[BC_INLET_FIXED] = new FluxFcnPerfectGasOutflowSA3D(ioData);
          }
          else {
            ff[BC_INLET_MOVING] = new FluxFcnPerfectGasInternalInflowSA3D(ioData);
            ff[BC_INLET_FIXED] = new FluxFcnPerfectGasInternalInflowSA3D(ioData);
          }
          ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_SLIP_WALL_MOVING] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_SLIP_WALL_FIXED] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_SYMMETRY] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnPerfectGasWallSA3D(ioData);
          ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnPerfectGasWallSA3D(ioData);

          if (ioData.schemes.ns.flux == SchemeData::ROE) {
            if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacRoeSA3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacRoeSA3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasExactJacRoeSA3D(gamma, ioData);
          }
          else if (ioData.schemes.ns.flux == SchemeData::HLLE) {
            if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacHLLESA3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE) {
              ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacHLLESA3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            }
            else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT) {
              fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
              exit(1);
            }
          }
        }
        else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
          ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
          ff -= BC_MIN_CODE;
          ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasOutflowKE3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnPerfectGasOutflowKE3D(ioData);
          ff[BC_INLET_MOVING] = new FluxFcnPerfectGasOutflowKE3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnPerfectGasOutflowKE3D(ioData);
          ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_SLIP_WALL_MOVING] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_SLIP_WALL_FIXED] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_SYMMETRY] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnPerfectGasWallKE3D(ioData);
          ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnPerfectGasWallKE3D(ioData);


          if (ioData.schemes.ns.flux == SchemeData::ROE) {
            if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacRoeKE3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacRoeKE3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasExactJacRoeKE3D(gamma, ioData);
          }
          else if (ioData.schemes.ns.flux == SchemeData::HLLE) {
            if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
              ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacHLLEKE3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE) {
              ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacHLLEKE3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
            }
            else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT) {
              fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
              exit(1);
            }
          }
        }
      } else {
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        if( ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL ){
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasOutflowEuler3D(ioData);
            ff[BC_OUTLET_FIXED] = new FluxFcnPerfectGasOutflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasGhidagliaEuler3D(ioData);
            ff[BC_OUTLET_FIXED]  = new FluxFcnPerfectGasGhidagliaEuler3D(ioData);
	  } else {
            ff[BC_OUTLET_MOVING] = 0;
            ff[BC_OUTLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL ){
          ff[BC_OUTLET_MOVING] = new FluxFcnPerfectGasInternalOutflowEuler3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnPerfectGasInternalOutflowEuler3D(ioData);
        }
        if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_INLET_MOVING] = new FluxFcnPerfectGasInflowEuler3D(ioData);
            ff[BC_INLET_FIXED] = new FluxFcnPerfectGasInflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_INLET_MOVING] = new FluxFcnPerfectGasGhidagliaEuler3D(ioData);
            ff[BC_INLET_FIXED]  = new FluxFcnPerfectGasGhidagliaEuler3D(ioData);
          }else{
            ff[BC_INLET_MOVING] = 0;
            ff[BC_INLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL ) {
          ff[BC_INLET_MOVING] = new FluxFcnPerfectGasInternalInflowEuler3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnPerfectGasInternalInflowEuler3D(ioData);
        }
        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnPerfectGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnPerfectGasWallEuler3D(ioData);

        if (ioData.schemes.ns.flux == SchemeData::VANLEER)
          ff[BC_INTERNAL] = new FluxFcnPerfectGasVanLeerEuler3D(ioData);
        else if (ioData.schemes.ns.flux == SchemeData::ROE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
            ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
            ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
            ff[BC_INTERNAL] = new FluxFcnPerfectGasExactJacRoeEuler3D(gamma, ioData);
        }
	else if (ioData.schemes.ns.flux == SchemeData::HLLE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
            ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacHLLEEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE) {
            ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacHLLEEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT) {
            fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!"); 
            exit(1);
          }
        }
        else if (ioData.schemes.ns.flux == SchemeData::HLLC) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
            ff[BC_INTERNAL] = new FluxFcnPerfectGasFDJacHLLCEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE) {
            ff[BC_INTERNAL] = new FluxFcnPerfectGasApprJacHLLCEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT) {
            fprintf(stderr,"Error... HLLC with Exact Jacobian not Implemented.. Aborting !!");
            exit(1);
          }
        }
      }
    }//for LIQUID (turbulent parts not yet implemented)
    else if (ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID){
      ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
      ff -= BC_MIN_CODE;
      if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
        if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
          ff[BC_INLET_MOVING] = new FluxFcnWaterCompressibleGhidagliaEuler3D(ioData);
          ff[BC_INLET_FIXED]  = new FluxFcnWaterCompressibleGhidagliaEuler3D(ioData);
        }else{
          ff[BC_INLET_MOVING]  = 0;
          ff[BC_INLET_FIXED]   = 0;
        }
      } else if (ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL) {
        ff[BC_INLET_MOVING] = new FluxFcnWaterCompressibleInternalInflowEuler3D(ioData);
        ff[BC_INLET_FIXED] = new FluxFcnWaterCompressibleInternalInflowEuler3D(ioData);
      }
      if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
        if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
          ff[BC_OUTLET_MOVING] = new FluxFcnWaterCompressibleGhidagliaEuler3D(ioData);
          ff[BC_OUTLET_FIXED]  = new FluxFcnWaterCompressibleGhidagliaEuler3D(ioData);
        }else{
          ff[BC_OUTLET_MOVING]  = 0;
          ff[BC_OUTLET_FIXED]   = 0;
        }
      }else if (ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL) {
        ff[BC_OUTLET_MOVING] = new FluxFcnWaterCompressibleInternalOutflowEuler3D(ioData);
        ff[BC_OUTLET_FIXED] = new FluxFcnWaterCompressibleInternalOutflowEuler3D(ioData);
      }
      ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_SLIP_WALL_MOVING] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_SLIP_WALL_FIXED] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_SYMMETRY] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnWaterCompressibleWallEuler3D(ioData);
      if (ioData.schemes.ns.flux == SchemeData::ROE) {
        if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
          ff[BC_INTERNAL] = new FluxFcnWaterCompressibleFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
        else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
          ff[BC_INTERNAL] = new FluxFcnWaterCompressibleApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
        else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
          ff[BC_INTERNAL] = new FluxFcnWaterCompressibleExactJacRoeEuler3D(gamma, ioData);
      }
    }//for JWL EOS
    else if (ioData.eqs.fluidModel.fluid == FluidModelData::JWL){
      ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
      ff -= BC_MIN_CODE;
      ff[BC_INLET_MOVING] = new FluxFcnJWLGhidagliaEuler3D(ioData);
      ff[BC_INLET_FIXED]  = new FluxFcnJWLGhidagliaEuler3D(ioData);
      ff[BC_OUTLET_MOVING] = new FluxFcnJWLGhidagliaEuler3D(ioData);
      ff[BC_OUTLET_FIXED]  = new FluxFcnJWLGhidagliaEuler3D(ioData);

      ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_SLIP_WALL_MOVING] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_SLIP_WALL_FIXED] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_SYMMETRY] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnJWLWallEuler3D(ioData);
      ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnJWLWallEuler3D(ioData);
      if (ioData.schemes.ns.flux == SchemeData::ROE) {
        if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
          ff[BC_INTERNAL] = new FluxFcnJWLFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
        else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
          ff[BC_INTERNAL] = new FluxFcnJWLApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
        else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
          ff[BC_INTERNAL] = new FluxFcnJWLExactJacRoeEuler3D(gamma, ioData);
      }
    }
  }
  else if (ioData.eqs.numPhase == 2){
    if (ioData.eqs.fluidModel.fluid == FluidModelData::GAS){
      if (ioData.eqs.fluidModel2.fluid == FluidModelData::GAS){
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_OUTLET_MOVING] = new FluxFcnGasInGasOutflowEuler3D(ioData);
            ff[BC_OUTLET_FIXED] = new FluxFcnGasInGasOutflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_OUTLET_MOVING] = new FluxFcnGasInGasGhidagliaEuler3D(ioData);
            ff[BC_OUTLET_FIXED]  = new FluxFcnGasInGasGhidagliaEuler3D(ioData);
          }else{
            ff[BC_OUTLET_MOVING] = 0;
            ff[BC_OUTLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL ){
          ff[BC_OUTLET_MOVING] = new FluxFcnGasInGasInternalOutflowEuler3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnGasInGasInternalOutflowEuler3D(ioData);
        }
        if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_INLET_MOVING] = new FluxFcnGasInGasInflowEuler3D(ioData);
            ff[BC_INLET_FIXED] = new FluxFcnGasInGasInflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_INLET_MOVING] = new FluxFcnGasInGasGhidagliaEuler3D(ioData);
            ff[BC_INLET_FIXED]  = new FluxFcnGasInGasGhidagliaEuler3D(ioData);
          }else{
            ff[BC_INLET_MOVING] = 0;
            ff[BC_INLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL ) {
          ff[BC_INLET_MOVING] = new FluxFcnGasInGasInternalInflowEuler3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnGasInGasInternalInflowEuler3D(ioData);
        }
        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnGasInGasWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnGasInGasWallEuler3D(ioData);
                                                                                                  
        if (ioData.schemes.ns.flux == SchemeData::ROE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE){
            ff[BC_INTERNAL] = new FluxFcnGasInGasFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE){
            ff[BC_INTERNAL] = new FluxFcnGasInGasApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT){
            ff[BC_INTERNAL] = new FluxFcnGasInGasExactJacRoeEuler3D(gamma, ioData);
          }
        }
      }
      else if (ioData.eqs.fluidModel2.fluid == FluidModelData::LIQUID){
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_OUTLET_MOVING] = new FluxFcnGasInLiquidOutflowEuler3D(ioData);
            ff[BC_OUTLET_FIXED] = new FluxFcnGasInLiquidOutflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_OUTLET_MOVING] = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
            ff[BC_OUTLET_FIXED]  = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_OUTLET_MOVING] = 0;
            ff[BC_OUTLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL ){
          ff[BC_OUTLET_MOVING] = new FluxFcnGasInLiquidInternalOutflowEuler3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnGasInLiquidInternalOutflowEuler3D(ioData);
        }
        if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
          if (ioData.schemes.bc.type == BoundarySchemeData::STEGER_WARMING){
            ff[BC_INLET_MOVING] = new FluxFcnGasInLiquidInflowEuler3D(ioData);
            ff[BC_INLET_FIXED] = new FluxFcnGasInLiquidInflowEuler3D(ioData);
          }else if (ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_INLET_MOVING] = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
            ff[BC_INLET_FIXED]  = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_INLET_MOVING] = 0;
            ff[BC_INLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL ) {
          ff[BC_INLET_MOVING] = new FluxFcnGasInLiquidInternalInflowEuler3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnGasInLiquidInternalInflowEuler3D(ioData);
        }
        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnGasInLiquidWallEuler3D(ioData);
                                                                                                  
        if (ioData.schemes.ns.flux == SchemeData::ROE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE){
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE){
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          }
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT){
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidExactJacRoeEuler3D(gamma, ioData);
          }
        }
      }
      else if (ioData.eqs.fluidModel2.fluid == FluidModelData::JWL){
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        ff[BC_OUTLET_MOVING] = new FluxFcnJWLInGasGhidagliaEuler3D(ioData);
        ff[BC_OUTLET_FIXED]  = new FluxFcnJWLInGasGhidagliaEuler3D(ioData);
        ff[BC_INLET_MOVING] = new FluxFcnJWLInGasGhidagliaEuler3D(ioData);
        ff[BC_INLET_FIXED]  = new FluxFcnJWLInGasGhidagliaEuler3D(ioData);

        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnJWLInGasWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnJWLInGasWallEuler3D(ioData);
                                                                                                  
        if (ioData.schemes.ns.flux == SchemeData::ROE) {
          ff[BC_INTERNAL] = new FluxFcnJWLInGasApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
        }
      }
    }
    else if (ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID){
      if (ioData.eqs.fluidModel2.fluid == FluidModelData::LIQUID){
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
          if(ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_OUTLET_MOVING] = new FluxFcnLiquidInLiquidGhidagliaEuler3D(ioData);
            ff[BC_OUTLET_FIXED]  = new FluxFcnLiquidInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_OUTLET_MOVING] = 0;
            ff[BC_OUTLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL ){
          ff[BC_OUTLET_MOVING] = new FluxFcnLiquidInLiquidInternalOutflowEuler3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnLiquidInLiquidInternalOutflowEuler3D(ioData);
        }
        if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
          if(ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_INLET_MOVING] = new FluxFcnLiquidInLiquidGhidagliaEuler3D(ioData);
            ff[BC_INLET_FIXED]  = new FluxFcnLiquidInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_INLET_MOVING] = 0;
            ff[BC_INLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL ) {
          ff[BC_INLET_MOVING] = new FluxFcnLiquidInLiquidInternalInflowEuler3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnLiquidInLiquidInternalInflowEuler3D(ioData);
        }
        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnLiquidInLiquidWallEuler3D(ioData);
                                                                                                                                                                                                     
        if (ioData.schemes.ns.flux == SchemeData::ROE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
            ff[BC_INTERNAL] = new FluxFcnLiquidInLiquidFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
            ff[BC_INTERNAL] = new FluxFcnLiquidInLiquidApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
            ff[BC_INTERNAL] = new FluxFcnLiquidInLiquidExactJacRoeEuler3D(gamma, ioData);
        }
      }
      else if (ioData.eqs.fluidModel2.fluid == FluidModelData::GAS){
        ff = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
        ff -= BC_MIN_CODE;
        if (ioData.bc.outlet.type == BcsFreeStreamData::EXTERNAL) {
          if(ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_OUTLET_MOVING] = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
            ff[BC_OUTLET_FIXED]  = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_OUTLET_MOVING] = 0;
            ff[BC_OUTLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.outlet.type == BcsFreeStreamData::INTERNAL ){
          ff[BC_OUTLET_MOVING] = new FluxFcnGasInLiquidInternalOutflowEuler3D(ioData);
          ff[BC_OUTLET_FIXED] = new FluxFcnGasInLiquidInternalOutflowEuler3D(ioData);
        }
        if (ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL) {
          if(ioData.schemes.bc.type == BoundarySchemeData::GHIDAGLIA){
            ff[BC_INLET_MOVING] = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
            ff[BC_INLET_FIXED]  = new FluxFcnGasInLiquidGhidagliaEuler3D(ioData);
          }else{
            ff[BC_INLET_MOVING] = 0;
            ff[BC_INLET_FIXED]  = 0;
          }
        }
        else if(ioData.bc.inlet.type == BcsFreeStreamData::INTERNAL ) {
          ff[BC_INLET_MOVING] = new FluxFcnGasInLiquidInternalInflowEuler3D(ioData);
          ff[BC_INLET_FIXED] = new FluxFcnGasInLiquidInternalInflowEuler3D(ioData);
        }
        ff[BC_ADIABATIC_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ADIABATIC_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SLIP_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcnGasInLiquidWallEuler3D(ioData);
        ff[BC_SYMMETRY] = new FluxFcnGasInLiquidWallEuler3D(ioData);
                                                                                                                                                                                                     
        if (ioData.schemes.ns.flux == SchemeData::ROE) {
          if (ioData.ts.implicit.jacobian == ImplicitData::FINITE_DIFFERENCE)
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidFDJacRoeEuler3D(gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::APPROXIMATE)
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidApprJacRoeEuler3D(rshift, gamma, betaRef, K1, cmach, shockreducer, prec, ioData);
          else if (ioData.ts.implicit.jacobian == ImplicitData::EXACT)
            ff[BC_INTERNAL] = new FluxFcnGasInLiquidExactJacRoeEuler3D(gamma, ioData);
        }
      }
    }
  }

  if (!ff) {
    com->fprintf(stderr, "*** Error: no valid choice for the flux function\n");
    exit(1);
  }

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
// note: only one value for beta is allowed for k-e
                                                                                                      
template<int dim>
RecFcn *SpaceOperator<dim>::createRecFcnLS(IoData &ioData)
{
  RecFcn *rf = 0;
  
  double beta = ioData.schemes.ls.beta;
  double eps = ioData.schemes.ls.eps;

  if (ioData.schemes.ls.reconstruction == SchemeData::CONSTANT)
    rf = new RecFcnConstant<1>;
  else if (ioData.schemes.ls.reconstruction == SchemeData::LINEAR) {
    if (ioData.schemes.ls.limiter == SchemeData::NONE)
      rf = new RecFcnLinear<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VANALBADA)
      rf = new RecFcnVanAlbada<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::BARTH)
      rf = new RecFcnBarth<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VENKAT)
      rf = new RecFcnVenkat<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::P_SENSOR)
      rf = new RecFcnLtdLinear<1>(beta, eps);
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

  if(ioData.bc.hydro.type  == BcsHydroData::GRAVITY)
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

template<int dim>
void SpaceOperator<dim>::storeGhost(DistSVec<double,dim> &V, DistVec<double> &Phi,
                       DistSVec<double,dim> &Vgf, DistVec<double> &weight)
{
  Vgf = -1.0;
  weight = 1.0;
  domain->storeGhost(V,Vgf,Phi);
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
    if (use_modal == false)  {
      int numLocSub = R.numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invcv = 1.0 / cv[i];
          for (int j=0; j<dim; ++j)
            r[i][j] *= invcv;
        }
      }
    }
  }
  irey = 0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::computeDerivativeOfResidual(DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &U, double dMach, DistSVec<double,dim> &R, DistSVec<double,dim> &dR, DistTimeState<dim> *timeState)
{

  dR = 0.0;

  varFcn->conservativeToPrimitive(U, *V);

//Remark: Error mesage for pointers
  if (dV == 0) {
    fprintf(stderr, "*** Error: Varible dV does not exist!\n");
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
  if(timeState) {
    irey = timeState->getInvReynolds();
    direy = timeState->getDerivativeOfInvReynolds(*geoState, X, dX, ctrlVol, dCtrlVol, *V, *dV, dMach);
  }
  else {
    irey = new DistVec<double>(domain->getNodeDistInfo());
    direy = new DistVec<double>(domain->getNodeDistInfo());
    *irey = 0.0;
    *direy = 0.0;
  }

  if (fet) {
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

  domain->computeDerivativeOfFiniteVolumeTerm(ctrlVol, dCtrlVol, *irey, *direy, fluxFcn, recFcn, *bcData, *geoState, X, dX, *V, *dV, *ngrad, egrad, dMach, dR);
  
  domain->getGradP(*ngrad);
  domain->getDerivativeOfGradP(*ngrad);

  if (volForce) {
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);
    domain->computeDerivativeOfVolumicForceTerm(volForce, ctrlVol, dCtrlVol, *V, *dV, dR);
  }

  if(dvms) {
    dvms->compute(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState, timeState, X, U, *V, R, failsafe, rshift);
    dvms->computeDerivative(fluxFcn, recFcn, fet, geoState->getConfig(), ctrlVol, *bcData, *geoState, timeState, X, U, *V, R, failsafe, rshift);
  }

  if (use_modal == false)  {
    int numLocSub = dR.numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double *dcv = dCtrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      double (*dr)[dim] = dR.subData(iSub);
      double (*drm)[dim] = (*dRm).subData(iSub);
      for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
        double invcv = 1.0 / cv[i];
        double dInvcv = ( (-1.0) / ( cv[i] * cv[i] ) ) * dcv[i];
        for (int j=0; j<dim; ++j)
          dr[i][j] = ( ( dr[i][j] * invcv ) + ( r[i][j] * dInvcv ) );
      }
    }
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
    if (use_modal == false)  {
      int numLocSub = R.numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
        for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
          double invcv = 1.0 / cv[i];
          for (int j=0; j<dim; ++j)
            r[i][j] *= invcv;
        }
      }
    }
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
    if (use_modal == false)  {
      int numLocSub = R.numLocSub();
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        double *cv = ctrlVol.subData(iSub);
        double (*r)[dim] = R.subData(iSub);
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

template<int dim>
// Included (MB)
void SpaceOperator<dim>::computeResidual(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistVec<double> &Phi,
                                         DistSVec<double,dim> &R, 
                                         DistExactRiemannSolver<dim> *riemann, int it)
{

  R = 0.0;
  DistSVec<double,1> PhiS(Phi.info(), reinterpret_cast<double (*)[1]>(Phi.data()));
  varFcn->conservativeToPrimitive(U, *V, &Phi);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0){
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, Phi, *V);
    timer->addNodalGradTime(t0);
  }

  if (dynamic_cast<RecFcnConstant<1> *>(recFcnLS) == 0)
    ngradLS->compute(geoState->getConfig(), X, ctrlVol, PhiS);

  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (xpol)
    xpol->compute(geoState->getConfig(),geoState->getInletNodeNorm(), X);

  if (fet) {
    domain->computeGalerkinTerm(fet, *bcData, *geoState, X, *V, R);
    if (!dvms) bcData->computeNodeValue(X);
  }

  if (volForce)
    domain->computeVolumicForceTerm(volForce, ctrlVol, *V, R);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);
  //if (dynamic_cast<RecFcnConstant<1> *>(recFcnLS) == 0)
  //  ngradLS->limit(recFcnLS, X, ctrlVol, PhiS);

  domain->computeFiniteVolumeTerm(ctrlVol, *riemann, fluxFcn, recFcn, *bcData,
                                  *geoState, X, *V, Phi, *ngrad, egrad,
                                  *ngradLS, R, it, failsafe,rshift);  

  if (use_modal == false)  {
    int numLocSub = R.numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r)[dim] = R.subData(iSub);
      for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
        double invcv = 1.0 / cv[i];
        for (int j=0; j<dim; ++j)
          r[i][j] *= invcv;
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::computeResidualLS(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                           DistVec<double> &Phi,  DistSVec<double,dim> &U,
                                           DistVec<double> &PhiF)
{
  PhiF = 0.0;
  DistSVec<double,1> PhiS(Phi.info(), reinterpret_cast<double (*)[1]>(Phi.data()));

  varFcn->conservativeToPrimitive(U, *V, &Phi);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0){
    double t0 = timer->getTime();
    ngrad->compute(geoState->getConfig(), X, ctrlVol, Phi, *V);
    timer->addLSNodalWeightsAndGradTime(t0);
  }

  if (dynamic_cast<RecFcnConstant<1> *>(recFcnLS) == 0)
    ngradLS->compute(geoState->getConfig(), X, ctrlVol, PhiS);
  if (egrad)
    egrad->compute(geoState->getConfig(), X);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *V);

  if (dynamic_cast<RecFcnConstant<1> *>(recFcnLS) == 0)
    ngradLS->limit(recFcnLS, X, ctrlVol, PhiS);

  domain->computeFiniteVolumeTermLS(fluxFcn, recFcn, recFcnLS, *bcData, *geoState, X, *V,
                                    *ngrad, *ngradLS, egrad, PhiS, PhiF);

  if (use_modal == false)  {
    int numLocSub = PhiF.numLocSub();
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *cv = ctrlVol.subData(iSub);
      double (*r) = PhiF.subData(iSub);
      for (int i=0; i<ctrlVol.subSize(iSub); ++i) {
        double invcv = 1.0 / cv[i];
        r[i] *= invcv;
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

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::storePreviousPrimitive(DistSVec<double,dim> &U, 
                                DistSVec<double,dim> &Vg, DistVec<double> &Phi,
                                DistSVec<double,dim> *Vgf, DistVec<double> *weight)
{

  varFcn->conservativeToPrimitive(U, Vg, &Phi);
  //if(riemann->RiemannUpdatePhase())
    //nothing to do, everything has been done when computing the fluxes
  if(Vgf && weight){
    //domain->storePrimitive(Vg,*Vgf,*weight,Phi);
    storeGhost(Vg,Phi,*Vgf,*weight);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::updatePhaseChange(DistSVec<double,dim> &Vg, 
                             DistSVec<double,dim> &U,
                             DistVec<double> &Phi, DistVec<double> &Phin,
                             DistSVec<double,dim> *Vgf, DistVec<double> *weight,
                             DistExactRiemannSolver<dim> *riemann)
{

  if (riemann->DoUpdatePhase()){
    // the solution of the riemann problem is used to replace values of a node
    // that changed nature (fluid1 to fluid2 or vice versa)
    // **** GFMPAR-like ****
    // checkWeights is to suppress 'cavitation'
    domain->checkWeights(Phi, Phin, *(riemann->getRiemannUpdate()), *(riemann->getRiemannWeight()));
    varFcn->updatePhaseChange(Vg, U, Phi, Phin, riemann->getRiemannUpdate(), riemann->getRiemannWeight());
  }

  else if (Vgf && weight)
    // an extrapolation is used to replace values of a node
    // that changed nature (fluid1 to fluid2 or vice versa)
    // **** GFMPAR-variation ****
    varFcn->updatePhaseChange(Vg, U, Phi, Phin, Vgf, weight);

  else
    // no solution of the riemann problem was computed and we just use
    // the values that we have to convert back to conservative variables
    // **** GFMP-like ****
    varFcn->primitiveToConservative(Vg, U, &Phi);

}

//------------------------------------------------------------------------------

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

  if (use_modal)  {
    DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
    unitCtrlVol = 1.0;
    if (fet)
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, unitCtrlVol, *V, A);

    domain->computeJacobianFiniteVolumeTerm(fluxFcn, *bcData, *geoState, *irey, X, unitCtrlVol, *V, A);

    if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, unitCtrlVol, *V, A);
  }
  else  {

    if (fet)
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);

    domain->computeJacobianFiniteVolumeTerm(fluxFcn, *bcData, *geoState, *irey, X, ctrlVol, *V, A);

    if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, ctrlVol, *V, A);
  }
  irey = 0;
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistSVec<double,dim> &U, DistMat<Scalar,neq> &A,
                                         DistVec<double> &Phi, DistExactRiemannSolver<dim> *riemann)
{

  //fprintf(stdout, "going through computeJacobian for two-phase flows\n");
#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V, &Phi);
#endif

  A = 0.0;

  if (use_modal)  {
    fprintf(stderr, "**Error: no modal for multiphase flows.. Exiting\n");
    exit(1);
  }
  else  {
    if (fet)
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);
    domain->computeJacobianFiniteVolumeTerm(*riemann, fluxFcn, *bcData, *geoState, *ngrad, *ngradLS, X, ctrlVol, *V, A, Phi);
    if (volForce)
      domain->computeJacobianVolumicForceTerm(volForce, ctrlVol, *V, A);
  }
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
                                      DistVec<double> &Phi, DistSVec<double,dim> &rhs)
{

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V, &Phi);
#endif

  if(xpol)
    domain->recomputeRHS(varFcn, *V, Phi, rhs, xpol, *bcData, *geoState, X);

}

                                      
//------------------------------------------------------------------------------                                                                                      
                                                                                                      
template<int dim>
template<class Scalar, int neq>
void SpaceOperator<dim>::computeViscousJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                         DistMat<Scalar,neq> &A)  {

  A = 0.0;
  if (fet)  {
    if (use_modal)  {
      DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
      unitCtrlVol = 1.0;
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, unitCtrlVol, *V, A);
      domain->finishJacobianGalerkinTerm(unitCtrlVol, A);
    }
    else  {
      domain->computeJacobianGalerkinTerm(fet, *bcData, *geoState, X, ctrlVol, *V, A);
      domain->finishJacobianGalerkinTerm(ctrlVol, A);
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
void SpaceOperator<dim>::applyBCsToSolutionVector(DistSVec<double,dim> &U)
{

  if (bcFcn)
    domain->applyBCsToSolutionVector(bcFcn, *bcData, U);

}

//------------------------------------------------------------------------------

template<int dim>
void SpaceOperator<dim>::applyBCsToResidual(DistSVec<double,dim> &U, DistSVec<double,dim> &R)
{

  if (bcFcn)
    domain->applyBCsToResidual(bcFcn, *bcData, U, R);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SpaceOperator<dim>::applyBCsToDerivativeOfResidual(DistSVec<double,dim> &U, DistSVec<double,dim> &dR)
{

//Remark: Error mesage for pointers
  if (dU == 0) {
    fprintf(stderr, "*** Error: Varible dU does not exist!\n");
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

  if (use_modal)  {
    DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
    unitCtrlVol = 1.0;
    domain->computeH1(fluxFcn, *bcData, *geoState, unitCtrlVol, *V, H1);
  }
  else
    domain->computeH1(fluxFcn, *bcData, *geoState, ctrlVol, *V, H1);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void SpaceOperator<dim>::computeH2(DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
				   DistSVec<double,dim> &U, DistMat<Scalar,dim> &H2,
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
template<class Scalar>
void SpaceOperator<dim>::computeH2LS(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                    DistVec<double> &Q, DistSVec<double,dim> &U,DistMat<Scalar,1> &H2)
{
                                                                                                                    
#ifdef DOUBLE_CHECK
//  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0) {
    ngrad->compute(geoState->getConfig(), X, ctrlVol, U);
    ngrad->limit(recFcn, X, ctrlVol, U);
//  }
#endif
                                                                                                                    
  H2 = 0.0;
                                                                                                                    
  varFcn->conservativeToPrimitive(U, *V);
  domain->computeH2LS(*geoState, X, *V, *ngrad, H2);
                                                                                                                    
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
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
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

  if (use_modal)  {
    DistVec<double> unitCtrlVol(domain->getNodeDistInfo());
    unitCtrlVol = 1.0;
    domain->computeMatVecProdH2(recFcn, X, unitCtrlVol, H2, aij, aji, bij, bji, V2, *distNodalGrad, prod);
  }
  else
    domain->computeMatVecProdH2(recFcn, X, ctrlVol, H2, aij, aji, bij, bji, V2, *distNodalGrad, prod);

// Included (MB*)
  if (bcFcn)
    if ((iod->eqs.type == EquationsData::NAVIER_STOKES) && (iod->eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
      if ((iod->bc.wall.integration == BcsWallData::WALL_FUNCTION) && (iod->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS))
        domain->applyBCsToProduct(bcFcn, *bcData, U, prod);

}

//------------------------------------------------------------------------------
                                                                                                                    
template<int dim>
template<class Scalar1, class Scalar2>
void SpaceOperator<dim>::applyH2LS(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                         DistVec<double> &U, DistMat<Scalar1,1> &H2,
                         DistSVec<double,1> &aij, DistSVec<double,1> &aji,
                         DistSVec<double,1> &bij, DistSVec<double,1> &bji,
                         DistVec<Scalar2> &p, DistVec<Scalar2> &prod)
{
                                                                                                                    
  int numLocSub = p.numLocSub();
                                                                                                                    
  domain->computeMatVecProdH2LS(recFcn, X, ctrlVol, H2, aij, aji, bij, bji, p, prod);
                                                                                                                    
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

  domain->computeMatVecProdH2T(recFcn, X, ctrlVol, H2, aij, aji, bij, bji, p,  prod, prod2, prod3, prod4);

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    distNodalGrad->computeT(geoState->getConfig(), X, ctrlVol, prod, prod3, prod4);

  domain->computeMatVecProdH2Tb(recFcn, X, ctrlVol, H2, *distNodalGrad, p, prod, prod2);

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    double (*locV)[dim] = V->subData(iSub);
    Scalar2 (*locp)[dim] = prod.subData(iSub);

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
void SpaceOperator<dim>::computeDerivativeOfGradP(DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &U, DistSVec<double,dim> &dU)
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

template<int dim>
void SpaceOperator<dim>::computeDerivativeOfGradP(DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &U)
{

  varFcn->conservativeToPrimitive(U, *V);
  *dV = 0.0;

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  {
    ngrad->compute(geoState->getConfig(), X, ctrlVol, *V);
    ngrad->computeDerivative(geoState->getConfigSA(), X, dX, ctrlVol, dCtrlVol, *V, *dV);
    ngrad->limitDerivative(recFcn, X, dX, ctrlVol, dCtrlVol, *V, *dV);
  }

  domain->getDerivativeOfGradP(*ngrad);

}

//------------------------------------------------------------------------------
