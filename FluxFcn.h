#ifndef _FLUX_FCN_H_
#define _FLUX_FCN_H_

#include "FluxFcnBase.h"
#include "FluxFcnDescSG.h"
#include "FluxFcnDescTait.h"
#include "FluxFcnDescJwl.h"

#include <IoData.h>
#include "BcDef.h"
#include <VarFcn.h>
#include <VarFcnJwl.h>
#include <VarFcnSGSA.h>
#include <VarFcnSGKE.h>
#include <VarFcnTait.h>
#include <VarFcnSGEuler.h>

#include <cassert>
#include <cmath>
#include <DebugTools.h>

//#define NDEBUG // if commented, assert statements are evaluated
//--------------------------------------------------------------------------
//
// This class is mostly a collection of FluxFcn used during the simulation
// and is created dynamically according to the needs of the simulation.
// The function of a specific FluxFcn can be accessed if a tag integer
// is provided for multiphase flows. A default tag value of zero
// always calls the member functions of the first FluxFcn.
//
//--------------------------------------------------------------------------
class FluxFcn {

private:
  int numPhases_;
  FluxFcnBase **ff_;

  VarFcn *vf_;
  
  void check(int tag=0) const{ 
    if(tag>=numPhases_){
      fprintf(stdout, "*** Error: An unknown fluid model with FluidID = %d is detected. Could be a software bug!\n", tag);
      DebugTools::PrintBacktrace();
      fflush(stdout);
      exit(1);
    }
    assert(tag<numPhases_);
  }

  FluxFcnBase *createFluxFcn(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, VarFcnBase *vfb, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE);

  //for Implicit Segregated Navier-Stokes solver ONLY!
  FluxFcnBase *createFluxFcnSeg1(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE);
  FluxFcnBase *createFluxFcnSeg2(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, VarFcnBase *vfb, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE);

public:
  FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE); 
  //for Implicit Segregated Navier-Stokes solver ONLY! (segPart = 1 or 2)
  FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, int segPart, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE); 

  FluxFcn() {}
  ~FluxFcn() {
    for(int i=0; i<numPhases_; i++)
      delete ff_[i];
    delete [] ff_;
  }

  VarFcn *getVarFcn() { return vf_; }


  //----- General Functions -----//
  void compute(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *flux, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->compute(length, irey, normal, normalVel, VL, VR, flux, useLimiter);
  }
  void computeJacobian(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *jacL, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->computeJacobian(length, irey, normal, normalVel, VL, VR, jacL, useLimiter);
  }
  void computeJacobians(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *jacL, double *jacR, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->computeJacobians(length, irey, normal, normalVel, VL, VR, jacL, jacR, useLimiter);
  }

  //----- Sensitivity Analysis Functions -----//
  void computeDerivative
  (
     double ire, double dIre, double *n, double *dn, double nv, double dnv, 
     double *vl, double *dvl, double *vr, double *dvr, double dmach, 
     double *f, double *df, int tag=0
  )
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->computeDerivative(ire,dIre,n,dn,nv,dnv,vl,dvl,vr,dvr,dmach,f,df);
  }

  void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *v, double *ub, double *dub, double *f, double *df,
    int tag=0
  )
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->computeDerivative(ire,dIre,n,dn,nv,dnv,v,ub,dub,f,df);
  }

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline
FluxFcn::FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, FluxFcnBase::Type typeJac) : vf_(vf){

  numPhases_ = iod.eqs.numPhase;
  ff_ = new FluxFcnBase *[numPhases_];

  if(vf_==0 || vf_->varFcn==0){
    fprintf(stderr, "*** Error: VarFcn is NULL in FluxFcn constructor\n");
    exit(1);
  }
  if(numPhases_>0){
    if(vf_->varFcn[0] == 0){
      fprintf(stderr, "*** Error: member varFcn[0] of VarFcn has not been initialized correctly\n");
      exit(1);
    }
    ff_[0] = createFluxFcn(rshift, ffType, iod.eqs.fluidModel, iod, vf_->varFcn[0], typeJac);
  }

  for(int iPhase=1; iPhase<numPhases_; iPhase++){
    map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(iPhase);
    if(it == iod.eqs.fluidModelMap.dataMap.end()){
      fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", iPhase);
      exit(1);
    }
    ff_[iPhase] = createFluxFcn(rshift, ffType, *it->second, iod, vf_->varFcn[iPhase], typeJac);
  }

}

//------------------------------------------------------------------------------

inline
FluxFcn::FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, int segPart, FluxFcnBase::Type typeJac) : vf_(vf){

  numPhases_ = iod.eqs.numPhase;
  ff_ = new FluxFcnBase *[numPhases_];

  if(vf_==0 || vf_->varFcn==0 || vf_->varFcn[0]==0 || numPhases_!=1 || ((segPart!=1)&&(segPart!=2)) ){
    fprintf(stderr, "*** Error: Unable to construct FluxFcn for Implicit Segregated NS solver!\n");
    exit(1);
  }
  if(segPart==1)  
    ff_[0] = createFluxFcnSeg1(rshift, ffType, iod.eqs.fluidModel, iod, typeJac); //an Euler (inviscid) VarFcn will be created
  if(segPart==2) 
    ff_[0] = createFluxFcnSeg2(rshift, ffType, iod.eqs.fluidModel, iod, vf_->varFcn[0], typeJac);
}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcn(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                    VarFcnBase *vfb, FluxFcnBase::Type typeJac){

  FluxFcnBase *localff;
  double gamma = iod.schemes.ns.gamma;

  if(fmodel.fluid == FluidModelData::GAS){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        VarFcnSGSA *vfsgsa = dynamic_cast<VarFcnSGSA *>(vfb);
        if(vfsgsa == 0){
          fprintf(stderr, "*** Error: a VarFcnSGSA is expected to create the associated FluxFcn\n");
          exit(-1);
        }
// Spalart-Allmaras for Stiffened Gas
        switch(ffType){

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnSGOutflowSA3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalOutflowSA3D(iod, vfsgsa, typeJac);
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnSGOutflowSA3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalInflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_INTERNAL:
            if (iod.schemes.ns.flux == SchemeData::ROE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacRoeSA3D(gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
                localff = new FluxFcnSGApprJacRoeSA3D(rshift, gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
                localff = new FluxFcnSGExactJacRoeSA3D(gamma, iod, vfsgsa, typeJac);
            }
            else if (iod.schemes.ns.flux == SchemeData::HLLE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLESA3D(gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLESA3D(rshift, gamma, iod, vfsgsa, typeJac);
              }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            break;

        }
      }
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
// k-epsilon turbulent model for Stiffened Gas
        VarFcnSGKE *vfsgke = dynamic_cast<VarFcnSGKE *>(vfb);
        if(vfsgke == 0){
          fprintf(stderr, "*** Error: a VarFcnSGKE is expected to create the associated FluxFcn\n");
          exit(-1);
        }
        switch(ffType){

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            localff = new FluxFcnSGOutflowKE3D(iod, vfsgke, typeJac);
            break;
  
          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallKE3D(iod, vfsgke, typeJac);
            break;
  
          case BC_INTERNAL:
            if (iod.schemes.ns.flux == SchemeData::ROE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacRoeKE3D(gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
                localff = new FluxFcnSGApprJacRoeKE3D(rshift, gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
                localff = new FluxFcnSGExactJacRoeKE3D(gamma, iod, vfsgke, typeJac);
            }
            else if (iod.schemes.ns.flux == SchemeData::HLLE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLEKE3D(gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLEKE3D(rshift, gamma, iod, vfsgke, typeJac);
              }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            break;

        }
      } // end - k-epsilon turbulent model for Stiffened Gas
    } // end - turbulence
    else{
// Euler or Navier-Stokes for Stiffened Gas 
      VarFcnSGEuler *vfsgeuler = dynamic_cast<VarFcnSGEuler *>(vfb);
      if(vfsgeuler == 0){
        fprintf(stderr, "*** Error: a VarFcnSGEuler is expected to create the associated FluxFcn\n");
        exit(-1);
      }
      switch(ffType){

        case BC_OUTLET_FIXED:
        case BC_OUTLET_MOVING:
          if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnSGInternalOutflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
            localff = new FluxFcnSGOutflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else{
            fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Stiffened Gas\n");
            exit(-1);
          }
          break;
     
        case BC_INLET_FIXED:
        case BC_INLET_MOVING:
          if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnSGInternalInflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
            localff = new FluxFcnSGInflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else{
            fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Stiffened Gas\n");
            exit(-1);
          }
          break;

        case BC_ADIABATIC_WALL_MOVING:
        case BC_ADIABATIC_WALL_FIXED:
        case BC_SLIP_WALL_MOVING:
        case BC_SLIP_WALL_FIXED:
        case BC_SYMMETRY:
        case BC_ISOTHERMAL_WALL_MOVING:
        case BC_ISOTHERMAL_WALL_FIXED:
          localff = new FluxFcnSGWallEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_INTERNAL:
          if (iod.schemes.ns.flux == SchemeData::VANLEER)
          {
            localff = new FluxFcnSGVanLeerEuler3D(iod, vfsgeuler, typeJac);
          }
          else if (iod.schemes.ns.flux == SchemeData::ROE) 
          {
            if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
            {
              localff = new FluxFcnSGFDJacRoeEuler3D(gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
            {
              localff = new FluxFcnSGApprJacRoeEuler3D(rshift, gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
            {
              localff = new FluxFcnSGExactJacRoeEuler3D(gamma, iod, vfsgeuler, typeJac);
            }
          }
          else if (iod.schemes.ns.flux == SchemeData::HLLE) 
          {
            if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
              localff = new FluxFcnSGFDJacHLLEEuler3D(gamma, iod, vfsgeuler, typeJac);
            else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
              localff = new FluxFcnSGApprJacHLLEEuler3D(rshift, gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
              fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
              exit(1);
            }
          }
          break;

      }
    } // end - Euler or Navier-Stokes for Stiffened Gas 
  } // end - All Stiffened Gas
  else if (fmodel.fluid == FluidModelData::LIQUID){
// Euler or Navier-Stokes for Tait EOS
    VarFcnTait *vftait = dynamic_cast<VarFcnTait *>(vfb);
    if(vftait == 0){
      fprintf(stderr, "*** Error: a VarFcnTait is expected to create the associated FluxFcn\n");
      exit(-1);
    }
    switch(ffType){

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
        if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
           iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
        else if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnTaitInternalOutflowEuler3D(iod, vftait, typeJac);
        else{
          fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Tait\n");
          exit(-1);
        }
        break;
    
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
           iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
        else if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnTaitInternalInflowEuler3D(iod, vftait, typeJac);
        else{
          fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Tait\n");
          exit(-1);
        }
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
        localff = new FluxFcnTaitWallEuler3D(iod, vftait, typeJac);
        break;

      case BC_INTERNAL:
        if (iod.schemes.ns.flux == SchemeData::ROE &&
            iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
          localff = new FluxFcnTaitApprJacRoeEuler3D(rshift, gamma, iod, vftait, typeJac);
        else{
          fprintf(stderr, "*** Error: only the Roe flux is available for Tait\n");
          exit(-1);
        }
        break;

    }
  } // end - Euler or Navier-Stokes for Tait EOS
  else if (fmodel.fluid == FluidModelData::JWL){
// Euler or Navier-Stokes for JWL EOS
    VarFcnJwl *vfjwl = dynamic_cast<VarFcnJwl *>(vfb);
    if(vfjwl == 0){
      fprintf(stderr, "*** Error: a VarFcnJwl is expected to create the associated FluxFcn\n");
      exit(-1);
    }
    switch(ffType){

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        localff = new FluxFcnJwlGhidagliaEuler3D(iod, vfjwl, typeJac);
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
        localff = new FluxFcnJwlWallEuler3D(iod, vfjwl, typeJac);
        break;

      case BC_INTERNAL:
        if (iod.schemes.ns.flux == SchemeData::ROE &&
            iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE){
          localff = new FluxFcnJwlApprJacRoeEuler3D(rshift, gamma, iod, vfjwl, typeJac);
        }else{
          fprintf(stderr, "*** Error: only the Roe flux is available for JWL\n");
          exit(-1);
        }
        break;

    }
    
  } // end - Euler or Navier-Stokes for JWL EOS

  return localff;

}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcnSeg1(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                        FluxFcnBase::Type typeJac){

  FluxFcnBase *localff;
  double gamma = iod.schemes.ns.gamma;

  if(fmodel.fluid == FluidModelData::GAS){
// Euler or Navier-Stokes for Stiffened Gas 
    VarFcnSGEuler *vfsgeuler = new VarFcnSGEuler(fmodel);

    switch(ffType){

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
        if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnSGInternalOutflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
          localff = new FluxFcnSGOutflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else{
          fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Stiffened Gas\n");
          exit(-1);
        }
        break;
     
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnSGInternalInflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
          localff = new FluxFcnSGInflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else{
          fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Stiffened Gas\n");
          exit(-1);
        }
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
        localff = new FluxFcnSGWallEuler3D(iod, vfsgeuler, typeJac);
        break;

      case BC_INTERNAL:
        localff = new FluxFcnSGApprJacRoeEuler3D(0, gamma, iod, vfsgeuler, typeJac);
        break;
    }
  } // end - All Stiffened Gas
  else {
    fprintf(stderr, "Exiting: No turbulence model for Tait, JWL, or multiphase simulations\n");
    exit(-1);
  }

  return localff;

}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcnSeg2(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                        VarFcnBase *vfb, FluxFcnBase::Type typeJac){

  FluxFcnBase *localff;
  double gamma = iod.schemes.ns.gamma;

  if(fmodel.fluid == FluidModelData::GAS){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        VarFcnSGSA *vfsgsa = dynamic_cast<VarFcnSGSA *>(vfb);
        if(vfsgsa == 0){
          fprintf(stderr, "*** Error: a VarFcnSGSA is expected to create the associated FluxFcn\n");
          exit(-1);
        }
// Spalart-Allmaras for Stiffened Gas
        switch(ffType){

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnSGOutflowSAturb3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalOutflowSAturb3D(iod, vfsgsa, typeJac);
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnSGOutflowSAturb3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalInflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_INTERNAL:
            localff = new FluxFcnSGRoeSAturb3D(gamma, iod, vfsgsa, typeJac);
            break;
        }
      }
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
// k-epsilon turbulent model for Stiffened Gas
        VarFcnSGKE *vfsgke = dynamic_cast<VarFcnSGKE *>(vfb);
        if(vfsgke == 0){
          fprintf(stderr, "*** Error: a VarFcnSGKE is expected to create the associated FluxFcn\n");
          exit(-1);
        }
        switch(ffType){

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            localff = new FluxFcnSGOutflowKEturb3D(iod, vfsgke, typeJac);
            break;
  
          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallKEturb3D(iod, vfsgke, typeJac);
            break;
  
          case BC_INTERNAL:
            localff = new FluxFcnSGRoeKEturb3D(gamma, iod, vfsgke, typeJac);
            break;

        }
      } // end - k-epsilon turbulent model for Stiffened Gas
      else {
        fprintf(stderr,"Error: Seg. solver is only implemented for SA and KE models.\n");
        exit(-1);
      }
    } // end - turbulence
    else {
      fprintf(stderr,"Error: Seg. solver is only implemented for Navier-Stokes Eqs.\n");
      exit(-1);
    }
  } else {
    fprintf(stderr, "Exiting: No turbulence model for Tait, JWL, or multiphase simulations\n");
    exit(-1);
  }
 

  return localff;

}
//------------------------------------------------------------------------------

#endif
