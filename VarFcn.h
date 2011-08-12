#ifndef _VAR_FCN_H_
#define _VAR_FCN_H_

#include <IoData.h>
#include <DistVector.h>
#include <Vector3D.h>
#include <VarFcnJwl.h>
#include <VarFcnSGEuler.h>
#include <VarFcnSGSA.h>
#include <VarFcnSGKE.h>
#include <VarFcnTait.h>
#include <VarFcnTaitSA.h>
#include <VarFcnTaitKE.h>

#include <cassert>
#include <cmath>
#include <complex>
typedef std::complex<double> bcomp;
#include <iostream>

using std::cout;
using std::endl;


//#define NDEBUG // if commented, assert statements are evaluated
//--------------------------------------------------------------------------
//
// This class is mostly a collection of VarFcn used during the simulation
// and is created dynamically according to the needs of the simulation.
// The function of a specific VarFcn can be accessed if a tag integer
// is provided for multiphase flows. A default tag value of zero
// always calls the member functions of the first VarFcn.
//
// In addition, gravity related and mesh velocity terms can be found
// in this class. (***It might be better to move them somewhere else?***)
//
//
//--------------------------------------------------------------------------
class VarFcn {

friend class FluxFcn;

private:
  int numPhases;
  int ghostPhase;
  int ghostId;
  VarFcnBase **varFcn;
  
  Vec3D meshVel;

  double gravity[3];
  double gravity_norm;
  double depth;

  int size()                {return numPhases;}
  int structPhase()          {return ghostPhase;}
  int structId()             {return ghostId;}
  void check(int tag) const {assert(0<=tag && tag<numPhases);} 

  VarFcnBase *createVarFcnBase(IoData &iod, FluidModelData& fluidModel);

public:
  VarFcn(IoData &iod); 
  VarFcn() { meshVel = 0.0; fprintf(stderr, "MeshVel initialized: %e %e %e\n", meshVel[0], meshVel[1], meshVel[2]);}
  ~VarFcn() {
    for(int i=0; i<numPhases; i++)
      delete varFcn[i];
    delete [] varFcn;
  }


  //----- Transformation Operators -----//
  /* Distributed Operators */
  template<int dim>
  void conservativeToPrimitive(DistSVec<double,dim> &U, DistSVec<double,dim> &V, DistVec<int> *tag = 0);
  template<int dim>
  void primitiveToConservative(DistSVec<double,dim> &V, DistSVec<double,dim> &U, DistVec<int> *tag = 0);
  template<int dim>
  void conservativeToPrimitiveDerivative(DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistVec<int> *tag = 0);
  template<int dim>
  void primitiveToConservativeDerivative(DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistVec<int> *tag = 0);
  /* Non-Distributed Operators */
  template<int dim>
  void conservativeToPrimitive(SVec<double,dim> &U, SVec<double,dim> &V, Vec<int> *tag = 0);
  template<int dim>
  void primitiveToConservative(SVec<double,dim> &V, SVec<double,dim> &U, Vec<int> *tag = 0);
  template<int dim>
  void conservativeToPrimitiveDerivative(SVec<double,dim> &U, SVec<double,dim> &dU, SVec<double,dim> &V, SVec<double,dim> &dV, Vec<int> *tag = 0);
  template<int dim>
  void primitiveToConservativeDerivative(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<int> *tag = 0);
  /* Node-Level Operators */
  void conservativeToPrimitive(double *U, double *V, int tag=0) { check(tag); varFcn[tag]->conservativeToPrimitive(U,V); }
  void primitiveToConservative(double *V, double *U, int tag=0) { check(tag); varFcn[tag]->primitiveToConservative(V,U); }
  void primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU, int tag=0) { check(tag); varFcn[tag]->primitiveToConservativeDerivative(V,dV,U,dU); }
  void conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV, int tag=0) { check(tag); varFcn[tag]->conservativeToPrimitiveDerivative(U,dU,V,dV); }
  int  conservativeToPrimitiveVerification(int glob, double *U, double *V, int tag=0) { check(tag); return varFcn[tag]->conservativeToPrimitiveVerification(glob,U,V); }

  void multiplyBydVdU(double *V, double *vec, double *res, int tag=0) { check(tag); varFcn[tag]->multiplyBydVdU(V,vec,res); }
  void multiplyBydVdU(double *V, bcomp *vec, bcomp *res, int tag=0) { check(tag); varFcn[tag]->multiplyBydVdU(V,vec,res); }
  void multiplyBydVdUT(double *V, double *vec, double *res, int tag=0) { check(tag); varFcn[tag]->multiplyBydVdUT(V,vec,res); }
  void multiplyBydVdUT(double *V, bcomp *vec, bcomp *res, int tag=0) { check(tag); varFcn[tag]->multiplyBydVdUT(V,vec,res); }
  void postMultiplyBydUdV(double *V, double *mat, double *res, int tag=0) { check(tag); varFcn[tag]->postMultiplyBydUdV(V,mat,res); }
  void postMultiplyBydVdU(double *V, double *mat, double *res, int tag=0) { check(tag); varFcn[tag]->postMultiplyBydVdU(V,mat,res); }
  void postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res, int tag=0) { check(tag); varFcn[tag]->postMultiplyBydUdV(V,mat,res); }
  void preMultiplyBydUdV(double *V, double *mat, double *res, int tag=0) { check(tag); varFcn[tag]->preMultiplyBydUdV(V,mat,res); }

  void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V, int tag=0) { check(tag); varFcn[tag]->extrapolatePrimitive(un,c,Vb,Vinter,V); }
  void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV, int tag=0) { check(tag); varFcn[tag]->extrapolateCharacteristic(n,un,c,Vb,dV); }
  void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW, int tag=0) { check(tag); varFcn[tag]->primitiveToCharacteristicVariations(n,V,dV,dW); }
  void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV, int tag=0) { check(tag); varFcn[tag]->characteristicToPrimitiveVariations(n,V,dW,dV); }

  //----- General Functions -----//
  int getType(int tag=0) const{ 
    check(tag); 
    return varFcn[tag]->getType(); }
  const char * pname(int dim, int tag=0) const{ return varFcn[tag]->pname[dim]; }

  bool doVerification() const{
    for(int iPhase=0; iPhase<numPhases; iPhase++)
      if(varFcn[iPhase]->doVerification())
        return true;
    return false;
  }

  double getDensity(double *V, int tag=0)   const{ check(tag); return varFcn[tag]->getDensity(V); }
  Vec3D  getVelocity(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->getVelocity(V); }
  double getVelocityX(double *V, int tag=0) const{ check(tag); return varFcn[tag]->getVelocityX(V); }
  double getVelocityY(double *V, int tag=0) const{ check(tag); return varFcn[tag]->getVelocityY(V); }
  double getVelocityZ(double *V, int tag=0) const{ check(tag); return varFcn[tag]->getVelocityZ(V); }
  double getPressure(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->getPressure(V); }

  void setDensity(const double density, double *V, int tag=0) { check(tag); return varFcn[tag]->setDensity(density,V); }
  void setPressure(const double p, double *V, int tag=0)  { check(tag); return varFcn[tag]->setPressure(p,V); }
  void setDensity(double *V, double *Vorig, int tag=0) { check(tag); return varFcn[tag]->setDensity(V,Vorig); }
  void setPressure(double *V, double *Vorig, int tag=0){ check(tag); return varFcn[tag]->setPressure(V,Vorig); }

  //checks that the Euler equations are still hyperbolic
  double checkPressure(double *V, int tag=0) const{check(tag); return varFcn[tag]->checkPressure(V); }

  double computeTemperature(double *V, int tag=0) const{check(tag); return varFcn[tag]->computeTemperature(V); }
  double computeRhoEnergy(double *V, int tag=0)   const{check(tag); return varFcn[tag]->computeRhoEnergy(V); }
  //this function computes the internal energy (=rho*e-0.5*rho*u^2)
  double computeRhoEpsilon(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->computeRhoEpsilon(V); }
  double getVelocitySquare(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->getVelocitySquare(V); }
  double getVelocityNorm(double *V, int tag=0)    const{ check(tag); return varFcn[tag]->getVelocityNorm(V); }

  double computeMachNumber(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->computeMachNumber(V); }
  double computeSoundSpeed(double *V, int tag=0)  const{ check(tag); return varFcn[tag]->computeSoundSpeed(V); }
  double computeSoundSpeed(const double density, const double entropy, int tag=0) const{check(tag); return varFcn[tag]->computeSoundSpeed(density,entropy); }
  double computeEntropy(const double density, const double pressure, int tag=0) const{check(tag); return varFcn[tag]->computeEntropy(density,pressure); }
  double computeIsentropicPressure(const double entropy, const double density, int tag=0) const{check(tag); return varFcn[tag]->computeIsentropicPressure(entropy,density); }
  double computeTotalPressure(double machr, double *V, int tag=0) const{check(tag); return varFcn[tag]->computeTotalPressure(machr,V); }

  double getTurbulentNuTilde(double *V, int tag=0)         const{check(tag); return varFcn[tag]->getTurbulentNuTilde(V); }
  double getTurbulentKineticEnergy(double *V, int tag=0)   const{check(tag); return varFcn[tag]->getTurbulentKineticEnergy(V); }
  double getTurbulentDissipationRate(double *V, int tag=0) const{check(tag); return varFcn[tag]->getTurbulentDissipationRate(V); }

  void rstVar(IoData &iod) { assert(numPhases==1); varFcn[0]->rstVar(iod); }
  Vec3D getDerivativeOfVelocity(double *dV, int tag=0) { check(tag); return varFcn[tag]->getDerivativeOfVelocity(dV); }
  double computeDerivativeOfTemperature(double *V, double *dV, int tag=0) { check(tag); return varFcn[tag]->computeDerivativeOfTemperature(V,dV); }
  double computeDerivativeOfMachNumber(double *V, double *dV, double dMach, int tag=0) { check(tag); return varFcn[tag]->computeDerivativeOfMachNumber(V,dV,dMach); }
  double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach, int tag=0) { check(tag); return varFcn[tag]->computeDerivativeOfSoundSpeed(V,dV,dMach); }
  double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dmach, int tag=0) { check(tag); return varFcn[tag]->computeDerivativeOfTotalPressure(machr,dmachr,V,dV,dmach); }
  double getDerivativeOfVelocityNorm(double *V, double *dV, int tag=0) { check(tag); return varFcn[tag]->getDerivativeOfVelocityNorm(V,dV); }
  double getDerivativeOfPressureConstant(int tag=0) { check(tag); return varFcn[tag]->getDerivativeOfPressureConstant(); }

  double specificHeatCstPressure(int tag=0) const{ check(tag); return varFcn[tag]->specificHeatCstPressure(); }
  double computePressureCoefficient(double *V, double pinfty, double mach, bool dimFlag, int tag=0) const{ check(tag); return varFcn[tag]->computePressureCoefficient(V,pinfty,mach,dimFlag); }

  //----- Mesh Velocity Related Functions -----//
  void setMeshVel(Vec3D &v)  { meshVel = v; }
  double *getMeshVel() { return meshVel.v; }
  double computeU2(double *V) const{ return V[1]*V[1]+V[2]*V[2]+V[3]*V[3]; }
  double computeWtU2(double *V) const{ 
    double v1 = V[1] - meshVel[0]; 
    double v2 = V[2] - meshVel[1]; 
    double v3 = V[3] - meshVel[2];
    return v1*v1+v2*v2+v3*v3;
  }
  double computeWtMachNumber(double *V, int tag=0) {
    check(tag);
    double m = sqrt(computeWtU2(V)) / computeSoundSpeed(V,tag);
    if (isnan(m))  {
      fprintf(stderr, "Nan: %e %e %e\n", meshVel[0], meshVel[1], meshVel[2]);
      exit(-1);
    }
    return m;
  }

  //----- Gravity Related Functions -----//
  double gravity_value() const { return gravity_norm; }
  double hydrostaticPressure(const double density, double *X) const{ 
    return density*(gravity_norm*depth + gravity[0]*X[0]+gravity[1]*X[1]+gravity[2]*X[2]); }
  double hydrodynamicPressure(double *V, double *X, int tag=0) const{ 
    check(tag);
    double localP = varFcn[tag]->getPressure(V);
    return localP - hydrostaticPressure(V[0],X); }
  double DerivativeHydrostaticPressure(const double density, const double ddensity, 
                                       double *X, double *dX) const{
    return ddensity*(gravity_norm*depth + gravity[0]*X[0]+gravity[1]*X[1]+gravity[2]*X[2])
         + density*(gravity[0]*dX[0]+gravity[1]*dX[1]+gravity[2]*dX[2]);
  }

  //----- Equation of State Parameters -----//
  /* PERFECT/STIFFENED GAS EOS */
  double getGamma(int tag=0)  const{ check(tag); return varFcn[tag]->getGamma(); }
  double getGamma1(int tag=0) const{ check(tag); return varFcn[tag]->getGamma1(); }
  double getPressureConstant(int tag=0) const{ check(tag); return varFcn[tag]->getPressureConstant(); }
  /* TAIT EOS */
  double getCv(int tag=0)         const{ check(tag); return varFcn[tag]->getCv(); }
  double getAlphaWater(int tag=0) const{ check(tag); return varFcn[tag]->getAlphaWater(); }
  double getBetaWater(int tag=0)  const{ check(tag); return varFcn[tag]->getBetaWater(); }
  double getPrefWater(int tag=0)  const{ check(tag); return varFcn[tag]->getPrefWater(); }
  /* JWL EOS */
  double getOmega(int tag=0) const{ check(tag); return varFcn[tag]->getOmega(); }
  double getOmegap1(int tag=0) const{ check(tag); return varFcn[tag]->getOmegap1(); }
  double getInvOmega(int tag=0) const{ check(tag); return varFcn[tag]->getInvOmega(); }
  double getInvOmegap1(int tag=0) const{ check(tag); return varFcn[tag]->getInvOmegap1(); }
  double getA1(int tag=0) const{ check(tag); return varFcn[tag]->getA1(); }
  double getA2(int tag=0) const{ check(tag); return varFcn[tag]->getA2(); }
  double getR1(int tag=0) const{ check(tag); return varFcn[tag]->getR1(); }
  double getR2(int tag=0) const{ check(tag); return varFcn[tag]->getR2(); }
  double getRhoref(int tag=0) const{ check(tag); return varFcn[tag]->getRhoref(); }
  double getR1r(int tag=0) const{ check(tag); return varFcn[tag]->getR1r(); }
  double getR2r(int tag=0) const{ check(tag); return varFcn[tag]->getR2r(); }
  double getPmin(int tag=0) const { check(tag); return varFcn[tag]->pmin; }

  //----- EOS related functions -----//
  double computeExponentials(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeExponentials(density); }
  double computeDerivativeOfExponentials(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeDerivativeOfExponentials(density); }
  double computeExponentials2(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeExponentials2(density); }
  double computeDerivativeOfExponentials2(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeDerivativeOfExponentials2(density); }
  double computeFrho(double *V, int tag=0) const{ check(tag); return varFcn[tag]->computeFrho(V); }
  double computeFrho(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeFrho(density); }
  double computeFrhop(double *V, int tag=0) const{ check(tag); return varFcn[tag]->computeFrhop(V); }
  double computeFrhop(const double density, int tag=0) const{ check(tag); return varFcn[tag]->computeFrhop(density); }

  VarFcnBase* getVarFcnBase(int tag = 0) const { check(tag); return varFcn[tag]; }

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline
VarFcn::VarFcn(IoData &iod)
{
  numPhases = iod.eqs.numPhase;
  int numRealFluid = numPhases;

  // In the case of Embedded Framework, add a VarFcn for Ghost
  if(iod.problem.framework == ProblemData::EMBEDDED)
    numPhases++;

  varFcn = new VarFcnBase *[numPhases];
  for(int iPhase = 0; iPhase < numPhases; iPhase++) 
    varFcn[iPhase] = 0;

  if(numRealFluid>0){
    varFcn[0] = createVarFcnBase(iod, iod.eqs.fluidModel);
  } else {
    fprintf(stderr,"*** Error: No FluidModel specified!\n");
    exit(1);
  }

  for(int iPhase=1; iPhase<numRealFluid; iPhase++){
    map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(iPhase);
    if(it == iod.eqs.fluidModelMap.dataMap.end()){
      fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", iPhase);
      exit(1);
    }
    varFcn[iPhase] = createVarFcnBase(iod, *it->second);
  }

  // In the case of a Embedded Simulation, Structure VarFcn points on VarFcn[0]
  // by default. It is just convienient. This can be changed eventually.
  if(iod.problem.framework == ProblemData::EMBEDDED) {
    ghostId = numPhases - 1;
    varFcn[ghostId] = createVarFcnBase(iod, iod.eqs.fluidModel);
    ghostPhase = 0; //same as varFcn[0]
  }

  meshVel = 0.0;

  gravity[0] = iod.eqs.gravity_x;
  gravity[1] = iod.eqs.gravity_y;
  gravity[2] = iod.eqs.gravity_z;
  gravity_norm = sqrt(gravity[0]*gravity[0]+gravity[1]*gravity[1]+gravity[2]*gravity[2]);
  depth = iod.bc.hydro.depth;
}  
  
//------------------------------------------------------------------------------

inline
VarFcnBase *VarFcn::createVarFcnBase(IoData &iod, FluidModelData &fluidModel) { 
  
  VarFcnBase *vf_;
 
  if(fluidModel.fluid == FluidModelData::PERFECT_GAS || 
     fluidModel.fluid == FluidModelData::STIFFENED_GAS) {
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES){
        vf_ = new VarFcnSGSA(fluidModel);
        //fprintf(stdout, "Debug: VarFcnSGSA created\n");
      }else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE){
        vf_ = new VarFcnSGKE(fluidModel);
        //fprintf(stdout, "Debug: VarFcnSGKE created\n");
      }
    }
    else{
      vf_ = new VarFcnSGEuler(fluidModel);
      //fprintf(stdout, "Debug: VarFcnSGEuler created\n");
    }
  }
  else if(fluidModel.fluid == FluidModelData::LIQUID){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES){
        vf_ = new VarFcnTaitSA(fluidModel);
        fprintf(stdout, "Debug: VarFcnTaitSA created\n");
      }else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE){
        vf_ = new VarFcnTaitKE(fluidModel);
        fprintf(stdout, "Debug: VarFcnTaitKE created\n");
      }
    }
    else{
      vf_ = new VarFcnTait(fluidModel);
      fprintf(stdout, "Debug: VarFcnTaitEuler created\n");
    }
  }else if(fluidModel.fluid == FluidModelData::JWL){
    vf_ = new VarFcnJwl(fluidModel);
    //fprintf(stdout, "Debug: VarFcnSGJwl created\n");
  }else{
    fprintf(stdout, "No VarFcn created\n");
    fflush(stdout);
    exit(1);
  }
  return vf_;
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitive(SVec<double,dim> &U, SVec<double,dim> &V, Vec<int> *tag)
{
  if(tag)
    for (int i=0; i<U.size(); ++i)
      varFcn[(*tag)[i]]->conservativeToPrimitive(U[i], V[i]);
  else
    for (int i=0; i<U.size(); ++i)
      varFcn[0]->conservativeToPrimitive(U[i], V[i]);
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitive(DistSVec<double,dim> &U, DistSVec<double,dim> &V, DistVec<int> *tag)
{
 
  int numLocSub = U.numLocSub();
 
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if(tag){
      int    *loctag   = tag->subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i) {
        varFcn[loctag[i]]->conservativeToPrimitive(u[i], v[i]);
	/*if (loctag[i] == 0) {
	  for (int k = 0; k < dim; ++k)
	    std::cout << u[i][k] << " " << v[i][k] << std::endl;
	}
	if (loctag[i] == 0) {
	  std::cout << std::endl;
	  }*/
      }

    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[0]->conservativeToPrimitive(u[i], v[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(SVec<double,dim> &V, SVec<double,dim> &U, Vec<int> *tag)
{
  if(tag)
    for (int i=0; i<U.size(); ++i)
      varFcn[(*tag)[i]]->primitiveToConservative(V[i], U[i]);
  else
    for (int i=0; i<U.size(); ++i)
      varFcn[0]->primitiveToConservative(V[i], U[i]);
}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservative(DistSVec<double,dim> &V, DistSVec<double,dim> &U, DistVec<int> *tag)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    if(tag){
      int *loctag = tag->subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[loctag[i]]->primitiveToConservative(v[i], u[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[0]->primitiveToConservative(v[i], u[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitiveDerivative(SVec<double,dim> &U, SVec<double,dim> &dU, SVec<double,dim> &V, SVec<double,dim> &dV, Vec<int> *tag)
{

 if (tag){
    for (int i=0; i<U.size(); ++i)
      varFcn[(*tag)[i]]->conservativeToPrimitiveDerivative(U[i], dU[i], V[i], dV[i]);
  }else{
    for (int i=0; i<U.size(); ++i){
      varFcn[0]->conservativeToPrimitiveDerivative(U[i], dU[i], V[i], dV[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::conservativeToPrimitiveDerivative(DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistVec<int> *tag)
{

  int numLocSub = U.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double (*dv)[dim] = dV.subData(iSub);
    if(tag){
      int *loctag = tag->subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[loctag[i]]->conservativeToPrimitiveDerivative(u[i], du[i], v[i], dv[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[0]->conservativeToPrimitiveDerivative(u[i], du[i], v[i], dv[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservativeDerivative(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<int> *tag)
{

 if (tag){
    for (int i=0; i<U.size(); ++i)
      varFcn[(*tag)[i]]->primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i]);
  }else{
    for (int i=0; i<U.size(); ++i){
      varFcn[0]->primitiveToConservativeDerivative(V[i], dV[i], U[i], dU[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void VarFcn::primitiveToConservativeDerivative(DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistVec<int> *tag)
{

  int numLocSub = U.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*u)[dim] = U.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    double (*v)[dim] = V.subData(iSub);
    double (*dv)[dim] = dV.subData(iSub);
    if(tag){
      int *loctag = tag->subData(iSub);
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[loctag[i]]->primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i]);
    }else{
      for (int i=0; i<U.subSize(iSub); ++i)
        varFcn[0]->primitiveToConservativeDerivative(v[i], dv[i], u[i], du[i]);
    }
  }

}

//------------------------------------------------------------------------------
#endif
