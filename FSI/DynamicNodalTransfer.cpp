/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel
 */
#include <iostream>
#include <FSI/DynamicNodalTransfer.h>
#include <IoData.h>
#include <Vector3D.h>
#include <list>
//------------------------------------------------------------------------------

DynamicNodalTransfer::DynamicNodalTransfer(IoData& iod, Communicator &c): com(c) , F(1), 
                           fScale(iod.ref.rv.tforce), XScale(iod.ref.rv.tlength),
                           tScale(iod.ref.rv.time), structure(iod,c)

{
  com.fprintf(stderr,"fscale = %e, XScale = %e, tScale = %e.\n", fScale, XScale, tScale);
{  Communication::Window<double> window(com, 1, &dts);
  window.fence(true);
  structure.sendTimeStep(&window);
  window.fence(false);
}
  dts /= tScale;
  com.barrier();
  com.fprintf(stderr,"dt = %e\n", dts);
 
}

//------------------------------------------------------------------------------

DynamicNodalTransfer::~DynamicNodalTransfer() {
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendForce() {
  std::pair<double *, int> embedded = structure.getTargetData();
  double *embeddedData = embedded.first;
  int length = embedded.second;
  Communication::Window<double> window(com, length, embeddedData);
  window.fence(true);
  window.accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  window.fence(false);
  structure.processReceivedForce();
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::getDisplacement(SVec<double,3>& structU) {
  Communication::Window<double> window(com, 3*structU.size()*sizeof(double), (double *)structU.data());
  window.fence(true);
  structure.sendDisplacement(&window);
  window.fence(false);

  structU = 1.0/XScale*structU;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateOutputToStructure(double dt, double dtLeft, SVec<double,3> &fs)
{
  if(F.size() != fs.size())
    F.resize(fs.size());
  F = fScale * fs;
}

//------------------------------------------------------------------------------

EmbeddedStructure::EmbeddedStructure(IoData& iod, Communicator &comm) : com(comm), 
                                             meshFile(iod.embeddedStructure.surfaceMeshFile)

{
  // read the input.
//  meshFile = iod.embeddedStructure.surfaceMeshFile;
  mode = iod.embeddedStructure.mode;
  tMax = iod.embeddedStructure.tMax;
  dt = iod.embeddedStructure.dt;
  omega = iod.embeddedStructure.omega;
  dx = iod.embeddedStructure.dx;
  dy = iod.embeddedStructure.dy;
  dz = iod.embeddedStructure.dz;

  com.fprintf(stderr,"Structure Information read from inputfile ...\n");
  com.fprintf(stderr,"mode = %d, tMax = %f, dt = %f, omega = %f, dx = %f, dy = %f, dz = %f.\n", mode, tMax, dt, omega, dx, dy, dz);
  
  // defaults
  nNodes = 0;
  X = 0;
  U = 0;
  F = 0;
  it = 0;
  
  // load structure nodes (from file).
  FILE *nodeFile = 0;
  nodeFile = fopen(meshFile,"r");
  if(!nodeFile) fprintf(stderr,"top file not found (in structure codes)\n");
  char c1[200], c2[200];
  int num0 = 0, num1 = 0, count, nInputs;
  double x1,x2,x3;
  fscanf(nodeFile, "%s %s\n", c1, c2);
  char debug[6]="Nodes";
  for (int i=0; i<5; i++) 
    if(debug[i]!=c1[i]) {fprintf(stderr,"Failed in reading file: %s\n", meshFile); exit(-1);}
  
  std::list<Vec3D> nodeList;
  std::list<Vec3D>::iterator it;

  while(1) {
    nInputs = fscanf(nodeFile,"%s", c1);
    if(nInputs!=1) break;
    num1 = atoi(c1);
    if(num0+1!=num1) break;

    fscanf(nodeFile,"%lf %lf %lf\n", &x1, &x2, &x3);
    nodeList.push_back(Vec3D(x1,x2,x3));
    num0 = num1;
  }
  nNodes = nodeList.size();

  fprintf(stderr,"In EmbeddedStructure: nNodes = %d\n",nNodes);  
  X = new (com) double[nNodes][3];
  
  count = 0;
  for (it=nodeList.begin(); it!=nodeList.end(); it++) {
    X[count][0] = (*it)[0]; 
    X[count][1] = (*it)[1]; 
    X[count][2] = (*it)[2]; 
    count++;
  }
  if(count!=nNodes) {fprintf(stderr,"WRONG!!\n"); exit(-1);}

  U = new (com) double[nNodes][3];
  F = new (com) double[nNodes][3];

  fclose(nodeFile);
    
/*
  if(coupled) {
    matchNodes = new MatchNodeSet *
    structExc = new StructExc(iod, mns, 6, strCom, com, 1);
  }
*/
}

//------------------------------------------------------------------------------

EmbeddedStructure::~EmbeddedStructure()
{
  if(X) delete[] X;
  if(U) delete[] U;
  if(F) delete[] F;
}

//------------------------------------------------------------------------------

pair<double*, int>
EmbeddedStructure::getTargetData() 
{
  //clear the force.
  for (int i=0; i<nNodes; i++)
    F[i][0] = F[i][1] = F[i][2] = 0.0;

  return pair<double*,int>((double*)F,nNodes);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendTimeStep(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc #1 sends the time.
//  Communication::Window<double> win0(com, 1, &dt);
{
  for(int i = 0; i < com.size(); ++i) {
    std::cout << "Sending the timestep (" << dt << ") to fluid " << i << std::endl;
    window->put(&dt, 0, 1, i, 0);
  }
}
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendDisplacement(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  it++;
  double time = dt*(double)it;

    for(int i = 0; i < nNodes; ++i) {
      U[i][0] = (1-cos(omega*time))*dx;
      U[i][1] = (1-cos(omega*time))*dy;
      U[i][2] = (X[i][1]*X[i][1])*(1-cos(omega*time))*dz;
    }
//    Communication::Window<double> win2(*intracom, 3*nNodes, (double *)data);
    std::cout << "Sending a displacement to fluid ";
    for(int i = 0; i < com.size(); ++i) {
      std::cout << i;
      window->put((double *)U, 0, 3*nNodes, i, 0);
      if(i!=com.size()-1) std::cout << ", ";
     }
    std::cout << std::endl;

}

//------------------------------------------------------------------------------

void
EmbeddedStructure::processReceivedForce()
{ 
  if(com.cpuNum()>0) return;

  double fx=0, fy=0, fz=0; //write the total froce.
  for(int i = 0; i < nNodes; ++i) {
    fx += F[i][0];  
    fy += F[i][1]; 
    fz += F[i][2];
  }
  std::cout << "Total force: " << fx << " " << fy << " " << fz << std::endl;
}

//------------------------------------------------------------------------------

