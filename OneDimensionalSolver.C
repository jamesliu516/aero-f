#include "OneDimensionalSolver.h"

#include <fstream>
#include <iostream>
using namespace std;
#include <cmath>

#include "FluxFcn.h"
#include "VarFcn.h"
#include "LocalRiemannDesc.h"
#include "FluidSelector.h"
#include "OneDimensionalSourceTerm.h"
#include "RefVal.h"

#include "IoData.h"
#include "Domain.h"

//------------------------------------------------------------------------------
OneDimensional::OneDimensional(IoData &ioData, Domain *domain) : 
  numPoints(ioData.oneDimensionalInfo.numPoints),
  ctrlVol(numPoints), ctrlSurf(numPoints+1), X(numPoints), Y(numPoints+1), 
  U(numPoints), V(numPoints), R(numPoints),
  gradV(numPoints), Phi(numPoints), Rphi(numPoints), gradPhi(numPoints),
  fluidId(numPoints),fluidIdn(numPoints),
  fluidSelector(ioData.eqs.numPhase,ioData), refVal(ioData.ref.rv)
{
  // equation modelling
  coordType  = ioData.oneDimensionalInfo.coordType;
  volumeType = ioData.oneDimensionalInfo.volumeType;

  // time and space domain definition
  maxDistance = ioData.oneDimensionalInfo.maxDistance;
  finalTime = ioData.ts.maxTime;
  cfl = ioData.ts.cfl0;

  // output
  frequency = ioData.output.transient.frequency;
  int sp = strlen(ioData.output.transient.prefix) + 1;
  if(ioData.output.transient.oneDimensionalRes[0] != 0) {
    outfile = new char[sp + strlen(ioData.output.transient.oneDimensionalRes)];
    sprintf(outfile, "%s%s", ioData.output.transient.prefix, ioData.output.transient.oneDimensionalRes);
  }
  else
    outfile = 0;


  // necessary for computation: varFcn to compute different state quantities
  //                            fluxFcn to compute fluxes at interfaces
  //                            riemann to solve riemann problem between two fluids
  //                            source to compute source term(s) for spherical and cylindrical problems
  varFcn = new VarFcn(ioData);

  fluxFcn = new FluxFcn *[3];
  fluxFcn[0] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn);
  fluxFcn[1] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn);
  fluxFcn[2] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn);

  //riemann = new LocalRiemannGfmparGasJWL(varFcn,0,1,0,MultiFluidData::RK2);
  riemann = new LocalRiemannGfmpGasJWL(varFcn,0,1);

  source = 0;
  if(volumeType == OneDimensionalInfo::CONSTANT_VOLUME){
    if(coordType == OneDimensionalInfo::CARTESIAN)
      source = new CartesianOneDSourceTerm(varFcn);
    else if(coordType == OneDimensionalInfo::CYLINDRICAL)
      source = new CylindricalOneDSourceTerm(varFcn);
    else if(coordType == OneDimensionalInfo::SPHERICAL)
      source = new SphericalOneDSourceTerm(varFcn);
  }else if(volumeType == OneDimensionalInfo::REAL_VOLUME){
    if(coordType == OneDimensionalInfo::SPHERICAL)
      source = new SphericalOneDSourceTerm2(varFcn);
    else if(coordType == OneDimensionalInfo::CYLINDRICAL)
      source = new CylindricalOneDSourceTerm2(varFcn);
  }

}
//------------------------------------------------------------------------------
OneDimensional::~OneDimensional(){

  delete varFcn;
  delete riemann;
  for(int i=0; i<3; i++) delete fluxFcn[i];
  delete [] fluxFcn;
  delete source;

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void OneDimensional::spatialSetup(){

  cout << "creating uniform mesh"<<endl;
  // at first uniform mesh
  X = 0.0;
  double dx = maxDistance/(numPoints-1);
  for(int i=0; i<numPoints; i++) X[i][0] = i*dx;

  cout << "computing cell boundaries" <<endl;
  Y = 0.0;
  Y[0][0] = X[0][0];
  for(int i=1; i<numPoints; i++) Y[i][0] = 0.5*(X[i-1][0]+X[i][0]);
  Y[numPoints][0] = X[numPoints-1][0];

  cout << "computing control volumes"<<endl;
  // computation of control volumes, for volumeType == CONSTANT_VOLUME
  ctrlVol[0][0] = 0.5*(X[1][0]-X[0][0]);
  for(int i=1; i<numPoints-1; i++)
    ctrlVol[i][0] = 0.5*(X[i+1][0]-X[i-1][0]);
  ctrlVol[numPoints-1][0] = 0.5*(X[numPoints-1][0]-X[numPoints-2][0]);

  // computation of control surfaces, for volumeType == CONSTANT_VOLUME
  cout << "computing control surfaces"<<endl;
  for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = 1.0;


  // in case a real Finite Volume approach is considered
  // with real cylindrical/spherical volumes, the control volumes and
  // surfaces must be computed differently as follows:
  if(volumeType == OneDimensionalInfo::REAL_VOLUME){

    if(coordType == OneDimensionalInfo::SPHERICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = (Y[i+1][0]-Y[i][0])*
                        (Y[i+1][0]*Y[i+1][0]+Y[i+1][0]*Y[i][0]+Y[i][0]*Y[i][0])/3.0;
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0]*Y[i][0];

    }else if(coordType == OneDimensionalInfo::CYLINDRICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = 0.5*(Y[i+1][0]-Y[i][0])*(Y[i+1][0]+Y[i][0]);
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0];
    }
  }

}
//------------------------------------------------------------------------------
void OneDimensional::temporalSetup(){
}
//------------------------------------------------------------------------------
void OneDimensional::stateInitialization(OneDimensionalInfo &data){
  V = 0.0;
  // initialize V
  for(int i=0; i<numPoints; i++){
    if(X[i][0]<data.interfacePosition){
      V[i][0] = data.density1;
      V[i][1] = data.velocity1;
      V[i][4] = data.pressure1;
    }else{
      V[i][0] = data.density2;
      V[i][1] = data.velocity2;
      V[i][4] = data.pressure2;
    }
  }

  // initialize Phi
  for(int i=0; i<numPoints; i++)
    Phi[i][0]  = -V[i][0]*(X[i][0]-data.interfacePosition);

  // initialize fluidId and U
  fluidSelector.getFluidId(fluidId,Phi);
  varFcn->primitiveToConservative(V,U,&fluidId);

  // compute boundary states
  double temp0[5] = {data.density1, data.velocity1, 0.0, 0.0, data.pressure1};
  double temp1[5] = {data.density2, data.velocity2, 0.0, 0.0, data.pressure2};
  varFcn->primitiveToConservative(temp0,BC[0],fluidId[0]);
  varFcn->primitiveToConservative(temp1,BC[1],fluidId[numPoints-1]);
  BCphi[0] = Phi[0][0]/V[0][0];
  BCphi[1] = Phi[numPoints-1][0]/V[numPoints-1][0];

  // output
  cout<<"**primitive boundary conditions are:"<<endl;
  cout<<"*    "<<temp0[0]<<" "<<temp0[1]<<" "<<temp0[2]<<" "<<temp0[3]<<" "<<temp0[4]<<endl;
  cout<<"*    "<<temp1[0]<<" "<<temp1[1]<<" "<<temp1[2]<<" "<<temp1[3]<<" "<<temp1[4]<<endl;
  cout<<"**conservative boundary conditions are:"<<endl;
  cout<<"*    "<<BC[0][0]<<" "<<BC[0][1]<<" "<<BC[0][2]<<" "<<BC[0][3]<<" "<<BC[0][4]<<endl;
  cout<<"*    "<<BC[1][0]<<" "<<BC[1][1]<<" "<<BC[1][2]<<" "<<BC[1][3]<<" "<<BC[1][4]<<endl;
  
}
//------------------------------------------------------------------------------
void OneDimensional::totalTimeIntegration(){

  // messages
  cout<<"***************************************************************"<<endl;
  cout<<"***  ctrlVol[0] = "<<ctrlVol[0][0]<<" --- ctrlVol[end] = "<<ctrlVol[numPoints-1][0]<<endl;
  cout<<"***  ctrlSur[0] = "<<ctrlSurf[0][0]<<" --- ctrlSur[end] = "<<ctrlSurf[numPoints][0]<<endl;
  cout<<"***************************************************************"<<endl;
  SVec<double,1> timeSteps(numPoints);

  double time = 0.0;
  double dt   = 0.0;
  int iteration = 0;

  resultsOutput(time,iteration);

  while(time<finalTime){
    // determine how much to advance in time
    computeTimeSteps(timeSteps);
    dt = cfl*timeSteps.min();

    if(time+dt>finalTime) dt = finalTime-time;
    if(iteration % frequency == 0)
      cout <<endl<<"*** Iteration " << iteration <<": Time = "<<time*refVal.time<<", and dt = "<<dt*refVal.time<<endl;
    time += dt;
    iteration++;

    // advance one iteration
    singleTimeIntegration(dt);

    if(iteration % frequency == 0)
      resultsOutput(time,iteration);
  }
  resultsOutput(time,iteration);

}
//------------------------------------------------------------------------------
void OneDimensional::computeTimeSteps(SVec<double,1> &timeSteps){

  // very crude CFL law
  double c = 0.0;
  varFcn->conservativeToPrimitive(U,V,&fluidId);

  for(int i=0; i<numPoints; i++){
    c = varFcn->computeSoundSpeed(V[i],fluidId[i]);
    timeSteps[i][0] = 0.5*(Y[i+1][0]-Y[i][0])/c;
  }

}
//------------------------------------------------------------------------------
void OneDimensional::singleTimeIntegration(double dt){
// for now, assume forward Euler

  fluidSelector.getFluidId(fluidId,Phi);

  R = 0.0;
  computeEulerFluxes();

  Rphi = 0.0;
  computeLevelSetFluxes();

  //preliminary update of U and Phi
  for(int i=0; i<numPoints; i++){
    for(int idim=0; idim<dim; idim++)
      U[i][idim] -= dt*R[i][idim]/ctrlVol[i][0];
    Phi[i][0]  -= dt*Rphi[i][0]/ctrlVol[i][0];
  }

  // store previous primitive with old fluidId
  varFcn->conservativeToPrimitive(U,V,&fluidId);

  // update fluidId
  fluidIdn = fluidId;
  fluidSelector.getFluidId(fluidId,Phi);

  //update PhaseChange with new fluidId
  varFcn->primitiveToConservative(V,U,&fluidId);

}
//------------------------------------------------------------------------------
void OneDimensional::computeEulerFluxes(){

  double normal[3] = {1.0, 0.0, 0.0};
  double length = 1.0;
  double normalVel = 0.0; // no ALE
  double flux[dim];
  int i,j;

  varFcn->conservativeToPrimitive(U,V,&fluidId);
  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    i = iEdge;
    j = iEdge+1;

    if(fluidId[i] == fluidId[j]){

      fluxFcn[0]->compute(length, 0.0, normal, normalVel, V[i], V[j], flux, fluidId[i]);
      for(int k=0; k<dim; ++k) {
        R[i][k] += ctrlSurf[iEdge+1][0]*flux[k];
        R[j][k] -= ctrlSurf[iEdge+1][0]*flux[k];
      }

    }else{
      double gradphi[3] = {1.0, 0.0, 0.0};
      double Vi[2*dim] = {V[i][0],V[i][1],V[i][2],V[i][3],V[i][4],V[i][0],V[i][1],V[i][2],V[i][3],V[i][4]};
      double Vj[2*dim] = {V[j][0],V[j][1],V[j][2],V[j][3],V[j][4],V[j][0],V[j][1],V[j][2],V[j][3],V[j][4]};
      double Wi[2*dim], Wj[2*dim];
      double Wir[dim], Wjr[dim];
      double wi, wj;
      double dx[3] = {X[j][0]-X[i][0], 0.0, 0.0};
      int iteration = 0;
      double fluxi[dim], fluxj[dim];

      riemann->computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,V[i],V[j],
                                      Wi,Wj,Wir,Wjr,wi,wj,dx,iteration);

      fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Wi, fluxi, fluidId[i]);
      fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wj, Vj, fluxj, fluidId[j]);
      for (int k=0; k<dim; k++){
        R[i][k] += ctrlSurf[iEdge+1][0]*fluxi[k];
        R[j][k] -= ctrlSurf[iEdge+1][0]*fluxj[k];
      }

    }
  }

  double dummy[dim];
  // flux at maxDistance
  fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[numPoints-1], BC[1], flux, fluidId[numPoints-1]);
  for (int k=0; k<dim; ++k)
    R[numPoints-1][k] += ctrlSurf[numPoints][0]*flux[k];

  // flux at left (for cartesian) - use of non-reflecting BC
  // flux at center (for cylindrical and spherical) - use of wall=symmetry
  normal[0] = ctrlSurf[0][0];
  if(coordType == OneDimensionalInfo::CARTESIAN)
    fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[0], BC[0], flux, fluidId[0]);
  else if(coordType == OneDimensionalInfo::CYLINDRICAL)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);
  else if(coordType == OneDimensionalInfo::SPHERICAL)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);

  for (int k=0; k<dim; ++k)
    R[0][k] -= ctrlSurf[0][0]*flux[k];

  // source term
  if(source){
    for(int i=0; i<numPoints; i++){
      source->computeSourceTerm(V[i],Y[i],Y[i+1],flux,fluidId[i]);
      for (int k=0; k<dim; ++k)
        R[i][k] += flux[k];
    }
  }

  // for debug
  //cout<<"flux[0] = "<<flux[0]<<" "<<flux[1]<<" "<<flux[2]<<" "<<flux[3]<<" "<<flux[4]<<endl;
  //cout<<"ctrlSurf[0] = "<<ctrlSurf[0][0]<<endl;
}
//------------------------------------------------------------------------------
void OneDimensional::computeLevelSetFluxes(){

  double uroe,flux;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    int i = iEdge;
    int j = iEdge+1;

    if(Phi[i][0] == Phi[j][0]){
      uroe = 0.5*(V[i][1]+V[j][1]);
    }
    else{
      uroe = (Phi[j][0]*V[j][1] - Phi[i][0]*V[i][1])/(Phi[j][0] - Phi[i][0]);
    }
    flux = V[i][1]*Phi[i][0] + V[j][1]*Phi[j][0] - fabs(uroe) *(Phi[j][0]-Phi[i][0]);
    Rphi[i][0] += 0.5*flux*ctrlSurf[iEdge+1][0];
    Rphi[j][0] -= 0.5*flux*ctrlSurf[iEdge+1][0];
  }

  // flux at maxDistance
  flux = (V[numPoints-1][1] > 0) ? Phi[numPoints-1][0]*V[numPoints-1][1] : V[numPoints-1][0]*BCphi[1]*V[numPoints-1][1];
  Rphi[numPoints-1][0] += flux*ctrlSurf[numPoints][0];

  // flux at center is zero since u=0 for radial and spherical
  // but for cartesian, flux needs to be computed
  if(coordType == OneDimensionalInfo::CARTESIAN){
    flux = (V[0][1] < 0) ? Phi[0][0]*V[0][1] : V[0][0]*BCphi[0]*V[0][1];
    Rphi[0][0] -= flux*ctrlSurf[0][0];
  }

  // source term
  if(source){
    for(int i=0; i<numPoints; i++){
      source->computeLevelSetSourceTerm(Phi[i], V[i], Y[i], Y[i+1], &flux, fluidId[i]);
      Rphi[i][0] += flux;
    }
  }

}
//------------------------------------------------------------------------------
void OneDimensional::resultsOutput(double time, int iteration){

  if(iteration==0) cout << "outputting results in file "<<outfile<<endl;
  fstream output;
  int sp = strlen(outfile)+1;
  char str[10];
  sprintf(str,"%d",iteration);
  char *currentfile = new char[sp + strlen(str) ];
  sprintf(currentfile, "%s%s", outfile, str);
  output.open(currentfile, fstream::out | fstream::trunc);

  if(iteration%frequency != 0) cout << "outputting last results in file "<<currentfile<<endl;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  output << "# time = " << time*refVal.time << endl;
  output << "# " << numPoints << endl;
  for(int i=0; i<numPoints; i++)
    output << X[i][0]*refVal.length <<" "<< V[i][0]*refVal.density <<" "<<V[i][1]*refVal.velocity <<" "<<V[i][4]*refVal.pressure <<" "<<Phi[i][0]/V[i][0]*refVal.length <<endl;
  output.close();

}
//------------------------------------------------------------------------------
