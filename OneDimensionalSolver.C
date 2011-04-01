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
OneDimensional::OneDimensional(int np,double* mesh,IoData &ioData, Domain *domain) : 
  numPoints(np),
  ctrlVol(numPoints), ctrlSurf(numPoints+1), X(numPoints), Y(numPoints+1), 
  U(numPoints), V(numPoints), R(numPoints),
  gradV(numPoints), Phi(numPoints), Rphi(numPoints), gradPhi(numPoints),
  fluidId(numPoints),fluidIdn(numPoints),
  fluidSelector(ioData.eqs.numPhase,ioData), refVal(ioData.ref.rv),
  Wr(numPoints),Vslope(numPoints),Phislope(numPoints),
  rupdate(numPoints), weight(numPoints), interfacialWi(numPoints),
  interfacialWj(numPoints), riemannStatus(numPoints), Phin(numPoints),
  programmedBurn(NULL)
{
  // equation modelling
  coordType  = ioData.oneDimensionalInfo.coordType;
  volumeType = ioData.oneDimensionalInfo.volumeType;

  // time and space domain definition
  maxDistance = mesh[np-1];
  finalTime = ioData.ts.maxTime;
  cfl = ioData.ts.cfl0;

  // Copy 1D mesh to X
  X = 0.0;
  for (int i = 0; i < np; ++i)
    X[i][0] = mesh[i];
  delete [] mesh;
  

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
  //riemann = new LocalRiemannGfmpGasJWL(varFcn,0,1);

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


  recFcn = createRecFcn(ioData);
  recFcnLS = createRecFcnLS(ioData);

  if (ioData.ts.type != TsData::EXPLICIT) {
    fprintf(stderr,"Only explcit integration available for the 1D solver!\n");
    exit(1);
  }

  switch (ioData.ts.expl.type) {

  case ExplicitData::RUNGE_KUTTA_4:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK4, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK4, numPoints);
    break;
  case ExplicitData::RUNGE_KUTTA_2:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK2, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK2, numPoints);
    break;
  case ExplicitData::FORWARD_EULER:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::FE, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::FE, numPoints);
    break;

  default:
    fprintf(stderr,"Unavailable explicit integration for the 1D solver!\n");
    exit(1);
    break;
  }

  loadSparseGrid(ioData);

  riemann = new ExactRiemannSolver<5>(ioData,rupdate,weight, interfacialWi,
				      interfacialWj, varFcn,
				      tabulationC);
  
  if (ioData.oneDimensionalInfo.programmedBurn.unburnedEOS >= 0) {
    programmedBurn = new ProgrammedBurn(ioData,&this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }
}
//------------------------------------------------------------------------------
OneDimensional::~OneDimensional(){

  delete varFcn;
  delete riemann;
  for(int i=0; i<3; i++) delete fluxFcn[i];
  delete [] fluxFcn;
  delete source;

  if (tabulationC)
    delete tabulationC;

}
//------------------------------------------------------------------------------
void OneDimensional::load1DMesh(IoData& ioData,int& numPts,double* &meshPoints) {

  if (ioData.input.oneDimensionalMesh[0] != 0) {
    char mesh1d[256];
    
    sprintf(mesh1d,"%s%s",ioData.input.prefix,ioData.input.oneDimensionalMesh);
    FILE* fin = fopen(mesh1d,"r");
    fscanf(fin, "%i",&numPts);
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i) {

      fscanf(fin,"%lf",&meshPoints[i]);
    }
  } else {

    numPts = ioData.oneDimensionalInfo.numPoints;
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i)
      meshPoints[i] = (double)i / (numPts-1)*ioData.oneDimensionalInfo.maxDistance;
  }
}
//------------------------------------------------------------------------------
void OneDimensional::spatialSetup(){

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
      if (varFcn->getType(1) != VarFcnBase::TAIT) {
	V[i][0] = data.density1;
	V[i][1] = data.velocity1;
	V[i][4] = data.pressure1;
      } else {
	V[i][0] = data.density1;
	V[i][1] = data.velocity1;
	V[i][4] = data.temperature1;
      }
    }else{
      if (varFcn->getType(0) != VarFcnBase::TAIT) {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.pressure2;
      } else {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.temperature2;
      }
    }
  }

  // initialize Phi
  for(int i=0; i<numPoints; i++)
    Phi[i][0]  = -V[i][0]*(X[i][0]-data.interfacePosition);

  // initialize fluidId and U
  fluidSelector.getFluidId(fluidId,Phi);
  varFcn->primitiveToConservative(V,U,&fluidId);

  fluidIdn = fluidId;

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

  time = 0.0;
  double dt   = 0.0;
  int iteration = 0;

  resultsOutput(time,iteration);

  while(time<finalTime){
    // determine how much to advance in time
    computeTimeSteps(timeSteps);
    dt = cfl*timeSteps.min();

    if (programmedBurn)
      programmedBurn->setCurrentTime(time,varFcn, U,fluidId,fluidIdn);

    if(time+dt>finalTime) dt = finalTime-time;
    if(iteration % frequency == 0)
      cout <<"*** Iteration " << iteration <<": Time = "<<time*refVal.time<<", and dt = "<<dt*refVal.time<<endl;
    iteration++;

    // advance one iteration
    singleTimeIntegration(dt);
    time += dt;

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
    //std::cout << "c = " << c <<  " " << fluidId[i] <<  std::endl;
    timeSteps[i][0] = 0.5*(Y[i+1][0]-Y[i][0])/c;
  }

}
//------------------------------------------------------------------------------
void OneDimensional::singleTimeIntegration(double dt){
// for now, assume forward Euler

  double Vtemp[5];

  fluidSelector.getFluidId(fluidId,Phi);

  riemannStatus = 0;

  Vintegrator->integrate(this,&OneDimensional::EulerF,
			 U,time,dt);
  //R = 0.0;
  //computeEulerFluxes();

  Phin = Phi;
  Phiintegrator->integrate(this,&OneDimensional::PhiF,
			   Phi,time,dt);
  //Rphi = 0.0;
  //computeLevelSetFluxes();

  //preliminary update of U and Phi
  for(int i=0; i<numPoints; i++){
    if (Phi[i][0]*Phin[i][0] < 0.0 &&
	!riemannStatus[i])
      Phi[i][0] = Phin[i][0];
  }

  // store previous primitive with old fluidId
  varFcn->conservativeToPrimitive(U,V,&fluidId);

  // update fluidId
  fluidIdn = fluidId;
  fluidSelector.getFluidId(fluidId,Phi);
  for(int i=0; i<numPoints; i++){

    if (fluidId[i] != fluidIdn[i]) { // Phase change
      
      if (!riemannStatus[i])
	std::cout << "Have a problem!" << std::endl;
      //fprintf(stderr,"Node %d changes phase!\n", i);
      //fprintf(stderr,"Wr[i] = %lf %lf %lf %lf %lf\n", Wr[i][0],Wr[i][1],Wr[i][2],Wr[i][3],Wr[i][4]);
      memcpy(V[i],Wr[i],sizeof(double)*5);
    }
  }

  double interfaceLocation;
  // initialize Phi
  /*for(int i=0; i<numPoints-1; i++) {

    if (Phi[i][0]*Phi[i+1][0] < 0.0) {
      interfaceLocation = (double)i+(Phi[i][0]/V[i][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
      break;
    }
  }
    
  for(int i=0; i<numPoints-1; i++) {
    Phi[i][0]  = V[i][0]*(interfaceLocation-i);
    }*/

  //update PhaseChange with new fluidId
  varFcn->primitiveToConservative(V,U,&fluidId);

  // Check solution (clip pressure, that is), if necessary
  if (varFcn->doVerification()) {
    for(int i=0; i<numPoints; i++){
      varFcn->conservativeToPrimitiveVerification(i+1, U[i], Vtemp, fluidId[i]);
    }
  }

  if (programmedBurn) {

    programmedBurn->setFluidIds(time, fluidId,U);
  }
}   

void OneDimensional::EulerF(double t, SVec<double,5>& y,SVec<double,5>& k) {

  R = 0.0;
  computeEulerFluxes(y);
  for(int i=0; i<numPoints; i++){
    for(int idim=0; idim<dim; idim++)
      k[i][idim] = -R[i][idim] / ctrlVol[i][0];
  }
}

//------------------------------------------------------------------------------
void OneDimensional::computeEulerFluxes(SVec<double,5>& y){

  double normal[3] = {1.0, 0.0, 0.0};
  double length = 1.0;
  double normalVel = 0.0; // no ALE
  double flux[dim];
  double Udummy[dim];
  int i,j,k;

  varFcn->conservativeToPrimitive(y,V,&fluidId);
  computeSlopes(V,Vslope,fluidId,true);
  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    i = iEdge;
    j = iEdge+1;

    double Vi[dim*2],Vj[dim*2],Vsi[dim],Vsj[dim];
    for (k = 0; k < dim; ++k) {
      Vsi[k] = Vslope[i][k]*(X[j][0]-X[i][0]);
      Vsj[k] = Vslope[j][k]*(X[j][0]-X[i][0]);
			  
    }
    
    recFcn->compute(V[i], Vsi,V[j], Vsj, Vi, Vj);

    //std::cout << "Hello" << std::endl;
    varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
    varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);

    if(fluidId[i] == fluidId[j]){

      fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Vj, flux, fluidId[i]);
      for(int k=0; k<dim; ++k) {
        R[i][k] += ctrlSurf[iEdge+1][0]*flux[k];
        R[j][k] -= ctrlSurf[iEdge+1][0]*flux[k];
      }

    }else{
      double gradphi[3] = {1.0, 0.0, 0.0};
      double Wi[2*dim], Wj[2*dim];
      double Wir[dim], Wjr[dim];
      double wi, wj;
      double dx[3] = {X[j][0]-X[i][0], 0.0, 0.0};
      int iteration = 0;
      double fluxi[dim], fluxj[dim];

      riemann->computeRiemannSolution(Vi,Vj,fluidId[i],fluidId[j],gradphi,varFcn,
				      Wi,Wj,i,j,i,dx);

      
      memcpy(Wr[i], Wj, sizeof(double)*5);
      memcpy(Wr[j], Wi, sizeof(double)*5);
      riemannStatus[i] = riemannStatus[j] = 1;

      //fprintf(stderr,"Edge %d-%d crosses interface\n",i,j);

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
void OneDimensional::PhiF(double t, SVec<double,1>& y,SVec<double,1>& k) {

  Rphi = 0.0;
  computeLevelSetFluxes(y);
  for(int i=0; i<numPoints; i++){
    k[i][0] = -Rphi[i][0] / ctrlVol[i][0];
  }
}

void OneDimensional::computeLevelSetFluxes(SVec<double,1>& y){

  double uroe,flux;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  computeSlopes(y,Phislope,fluidId,false);
  double Phii[1],Phij[1],Phisi[1],Phisj[1];
  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    int i = iEdge;
    int j = iEdge+1;

    Phisi[0] = Phislope[i][0]*(X[j][0]-X[i][0]);
    Phisj[0] = Phislope[j][0]*(X[j][0]-X[i][0]);

    recFcnLS->compute(y[i], Phisi,y[j], Phisj, Phii, Phij);

    double uav = 0.5*(V[i][1]+V[j][1]);

    if(uav > 0.0){
      flux = Phii[0]*uav;
    }
    else{
      flux = Phij[0]*uav;
    }
    Rphi[i][0] += flux*ctrlSurf[iEdge+1][0];
    Rphi[j][0] -= flux*ctrlSurf[iEdge+1][0];
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
      source->computeLevelSetSourceTerm(y[i], V[i], Y[i], Y[i+1], &flux, fluidId[i]);
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
    output << X[i][0]*refVal.length <<" "<< V[i][0]*refVal.density <<" "<<V[i][1]*refVal.velocity <<" "<<varFcn->getPressure(V[i],fluidId[i])*refVal.pressure <<" "<<Phi[i][0]/V[i][0]*refVal.length << " " << fluidId[i] << " " << varFcn->computeTemperature(V[i],fluidId[i])*refVal.temperature << endl;
  output.close();

}
//------------------------------------------------------------------------------

template <int neq>
void OneDimensional::computeSlopes(SVec<double,neq>& VV, SVec<double,neq>& slopes,
				   Vec<int>& fid,bool crossInterface) {
  
  int i,j;
  double r,sig;
  int stat = 0;
  slopes = 0.0;
  for(i=1; i<numPoints-1; ++i) {
    
    if (crossInterface) {
      if (fid[i] == fid[i+1] &&
	  fid[i] == fid[i-1])
	stat = 0;
      else if (fid[i] == fid[i+1])
	stat = 1;
      else if (fid[i] == fid[i-1])
	stat = 2;
      else
	stat = 3;
    }
    for (j = 0; j < neq; ++j) {

      if (stat == 0)
	slopes[i][j] = (VV[i+1][j]-VV[i-1][j])/(X[i+1][0]-X[i-1][0]);
      else if (stat == 1)
	slopes[i][j] = (VV[i+1][j]-VV[i][j])/(X[i+1][0]-X[i][0]);
      else if (stat == 2)
	slopes[i][j] = (VV[i][j]-VV[i-1][j])/(X[i][0]-X[i-1][0]);
      else
	slopes[i][j] = 0.0;

      //slopes[i][j] *= 0.5;
    }
  }

}

RecFcn* OneDimensional::createRecFcn(IoData &ioData)
{

  RecFcn *rf = 0;

  double beta = ioData.schemes.ns.beta;
  double eps = ioData.schemes.ns.eps;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    /*if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
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
      }*/
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
    fprintf(stderr, "*** Error: no valid choice for the reconstruction\n");
    exit(1);
  }

  return rf;

}

RecFcn* OneDimensional::createRecFcnLS(IoData &ioData)
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

void OneDimensional::loadSparseGrid(IoData& ioData) {

  if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
    // only the ioData.eqs.fluidModel is considered since only the Riemann invariant of one EOS is tabulated!
    // (no need to specify two different EOS)

    double *refIn  = new double[2];
    double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = pow(ioData.ref.rv.density,-ioData.eqs.fluidModel.jwlModel.omega)*ioData.ref.rv.velocity*ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.velocity;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);

  }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){

    double *refIn = new double[5]; double *refOut = new double[2];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = ioData.ref.rv.pressure;
    refIn[2] = ioData.ref.rv.density;
    refIn[3] = ioData.ref.rv.pressure;
    refIn[4] = ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.density;
    refOut[1] = ioData.ref.rv.density;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);
  }else{
    tabulationC = 0;
  }
}
