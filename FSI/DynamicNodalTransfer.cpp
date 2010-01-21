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
#include <MatchNode.h>
#include <StructExc.h>
#include <DistVector.h>
#include <list>
#include <map>
//------------------------------------------------------------------------------

DynamicNodalTransfer::DynamicNodalTransfer(IoData& iod, Communicator &c, Communicator &sc, Timer *tim): com(c) , F(1), 
                           fScale(iod.ref.rv.tforce), XScale(iod.ref.rv.tlength), UScale(iod.ref.rv.tvelocity),
                           tScale(iod.ref.rv.time), structure(iod,c,sc,tim)

{
  timer = tim;
  com.fprintf(stderr,"fscale = %e, XScale = %e, tScale = %e.\n", fScale, XScale, tScale);

{  Communication::Window<double> window(com, 1, &dts);
  window.fence(true);
  structure.sendTimeStep(&window);
  window.fence(false);
}

{  Communication::Window<double> window2(com, 1, &tMax);
  window2.fence(true);
  structure.sendMaxTime(&window2);
  window2.fence(false);
}

  algNum = structure.getAlgorithmNumber();
  dts /= tScale;
  tMax /= tScale;
  com.barrier();
  com.fprintf(stderr,"dts = %e, tMax = %e\n", dts, tMax);
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

  com.barrier();
  double t0 = timer->getTime();
  Communication::Window<double> window(com, 3*length*sizeof(double), embeddedData); 
  window.fence(true);
  window.accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  window.fence(false);
  timer->addEmbedComTime(t0);

  structure.processReceivedForce();
  
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::getDisplacement(SVec<double,3>& structU, SVec<double,3>& structUdot) {

  double UandUdot[2*structU.size()*3];
  Communication::Window<double> window(com, 2*3*structU.size()*sizeof(double), (double *)UandUdot);

  window.fence(true);
  structure.sendDisplacement(&window);
  window.fence(false);

  for(int i=0; i<structU.size(); i++) 
    for(int j=0; j<3; j++) {
      structU[i][j] = UandUdot[i*3+j];
      structUdot[i][j] = UandUdot[(structU.size()+i)*3+j];
    }
//  com.fprintf(stderr,"norm of received disp = %e.\n", structU.norm());

  structU = 1.0/XScale*structU;
  structUdot = 1.0/UScale*structUdot;

}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateOutputToStructure(double dt, double dtLeft, SVec<double,3> &fs)
{
  if(F.size() != fs.size()) {
//    fprintf(stderr,"force vector resized (from %d to %d)!\n", F.size(), fs.size());
    F.resize(fs.size());
  }
  F = fs;
}

//------------------------------------------------------------------------------

EmbeddedStructure::EmbeddedStructure(IoData& iod, Communicator &comm, Communicator &strCom, Timer *tim) : com(comm), 
                                             meshFile(iod.embeddedStructure.surfaceMeshFile),
                                             matcherFile(iod.embeddedStructure.matcherFile),
                                             tScale(iod.ref.rv.time), XScale(iod.ref.rv.tlength),
                                             UScale(iod.ref.rv.tvelocity)
{
  timer = tim;

  // read the input.
  coupled = ((iod.embeddedStructure.type == EmbeddedStructureInfo::TWOWAY) ||
             (iod.embeddedStructure.type == EmbeddedStructureInfo::ONEWAY)) ? true : false;

  // ---- for 2-way coupling only -----
  dim2Treatment = (iod.embeddedStructure.dim2Treatment == EmbeddedStructureInfo::YES) ? true : false;

  // ---- for 1-way coupling (not forced-motion) only -----
  oneWayCoupling = (iod.embeddedStructure.type == EmbeddedStructureInfo::ONEWAY) ? true : false;

  // ---- for forced-motion only ------
  if (iod.embeddedStructure.forcedMotionMode == EmbeddedStructureInfo::HEAVING)
    mode = 1;
  else if (iod.embeddedStructure.forcedMotionMode == EmbeddedStructureInfo::CONSTHEAVING)
    mode = 2;
  else 
    mode = 99; //for debugging use.

  tMax = iod.embeddedStructure.tMax;
  dt = iod.embeddedStructure.dt;
  omega = 2.0*acos(-1.0)*iod.embeddedStructure.omega;
  dx = iod.embeddedStructure.dx;
  dy = iod.embeddedStructure.dy;
  dz = iod.embeddedStructure.dz;
  t0= iod.embeddedStructure.t0;
  // ----------------------------------

  com.fprintf(stderr,"*** Embedded structure information read from inputfile ***\n");
  com.fprintf(stderr,"  coupled = %d,  dim2Treatment = %d, oneWayCoupling = %d\n", coupled, dim2Treatment, oneWayCoupling);
  if(!coupled)
    com.fprintf(stderr,"  (forced motion) mode = %d, tMax = %f, dt = %f, omega = %f, dx = %f, dy = %f, dz = %f.\n", mode, tMax, dt, omega, dx, dy, dz);

  // defaults
  nNodes = 0;
  X = 0;
  U = 0;
  Udot = 0;
  UandUdot = 0;
  F = 0;
  it = 0;
  structExc = 0;
  algNum = 6; //the default is A6. For 2-way coupling, algNum will be modified by the info from structure.

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
  Udot = new (com) double[nNodes][3];
  UandUdot = new (com) double[nNodes*2][3];
  F = new (com) double[nNodes][3];

  fclose(nodeFile);
   

  //for 2-D simulation only: pairing nodes in cross-section direction
  if(dim2Treatment) {
    double pairTol = 1.0e-6;
    for (int i=0; i<nNodes; i++) 
      for(int j=i; j<nNodes; j++)
        if(std::abs(X[i][0]+X[j][0])<pairTol &&
           std::abs(X[i][1]-X[j][1])<pairTol &&
           std::abs(X[i][2]-X[j][2])<pairTol)
          pairing[i] = j;
  }

  if(coupled) {
    MatchNodeSet **mns = new MatchNodeSet *[1];
    if(com.cpuNum() == 0)
      mns[0] = new MatchNodeSet(matcherFile);
    else 
      mns[0] = new MatchNodeSet();
    structExc = new StructExc(iod, mns, 6, &strCom, &com, 1);
    structExc->negotiate();
    structExc->getInfo();
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
//    tMax = 0.02402e-6; //TODO: FIX THIS

    int *locToGlob = new int[1];
    locToGlob[0] = 0;
    di = new DistInfo(1, 1, 1, locToGlob, &com);

    di->setLen(0,nNodes);
    di->finalize(false);
  
    algNum = structExc->getAlgorithmNumber();
  }

}

//------------------------------------------------------------------------------

EmbeddedStructure::~EmbeddedStructure()
{
  if(X) delete[] X;
  if(U) delete[] U;
  if(Udot) delete[] Udot;
  if(UandUdot) delete[] UandUdot;
  if(F) delete[] F;
  if(structExc) delete structExc;  
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
  std::cout << "Sending the timestep (" << dt << ") to fluid " << std::endl;
{
  for(int i = 0; i < com.size(); ++i) {
    window->put(&dt, 0, 1, i, 0);
  }
}
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendMaxTime(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc #1 sends the time.
  std::cout << "Sending the max time (" << tMax << ") to fluid " << std::endl;
{
  for(int i = 0; i < com.size(); ++i) {
    window->put(&tMax, 0, 1, i, 0);
  }
}
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendDisplacement(Communication::Window<double> *window)
{
  if(coupled) {
     DistSVec<double,3> Y0(*di, X);
     DistSVec<double,3> V(*di, U);
     DistSVec<double,3> Ydot(*di, Udot);
     DistSVec<double,3> Y(*di);
     Y = Y0;

     structExc->getDisplacement(Y0, Y, Ydot, V);
     V = XScale*V;
     Ydot = UScale*Ydot;
  }
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  it++;
  if(!coupled) {
    double time = t0 + dt*(double)it;
    
    if (mode==1) //heaving
      for(int i=0; i < nNodes; ++i) {
        U[i][0] = (1-cos(omega*time))*dx;
        U[i][1] = (1-cos(omega*time))*dy;
        U[i][2] = (1-cos(omega*time))*dz;
        Udot[i][0] = dx*omega*sin(omega*time);
        Udot[i][1] = dy*omega*sin(omega*time); 
        Udot[i][2] = dz*omega*sin(omega*time); 
      }
    else if (mode==2) //heaving with a constant velocity (in this case dx dy dz are velocity
      for(int i=0; i < nNodes; ++i) {
        U[i][0] = time*dx;
        U[i][1] = time*dy;
        U[i][2] = time*dz;
        Udot[i][0] = dx;
        Udot[i][1] = dy;
        Udot[i][2] = dz;
      }
    else if (mode==99) // for debugging use.
      for(int i=0; i < nNodes; ++i) { // expand / shrink the structure in y-z plane w.r.t. the origin
        double cosTheta = X[i][1]/sqrt(X[i][1]*X[i][1]+X[i][2]*X[i][2]);
        double sinTheta = X[i][2]/sqrt(X[i][1]*X[i][1]+X[i][2]*X[i][2]);
        U[i][0] = 0.0;
        U[i][1] = dy*time*cosTheta;
        U[i][2] = dz*time*sinTheta;
        Udot[i][0] = 0.0;
        Udot[i][1] = dy*cosTheta;
        Udot[i][2] = dz*sinTheta;
      }
  }

  for(int i=0; i<nNodes; i++)
    for(int j=0; j<3; j++) {
      UandUdot[i][j] = U[i][j];
      UandUdot[(i+nNodes)][j] = Udot[i][j];
    }

  double t0 = timer->getTime();
  for(int i = 0; i < com.size(); ++i) 
    window->put((double*)UandUdot, 0, 2*3*nNodes, i, 0);
  timer->addEmbedComTime(t0);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::processReceivedForce()
{ 
  if(coupled) {

    if(dim2Treatment) { // averaging the force on paired nodes
      std::map<int,int>::iterator it;
      for(it=pairing.begin(); it!=pairing.end(); it++) {
        double Fave;
        for(int iDim=0; iDim<3; iDim++) {
          Fave = 0.5*(F[it->first][iDim] + F[it->second][iDim]);
          F[it->first][iDim] = F[it->second][iDim] = Fave;
        }
      }
    } 

    DistSVec<double,3> f(*di, F);
  
    if(oneWayCoupling) // send a 0 force to structure
      f = 0.0;

    structExc->sendForce(f);
  }

  if(com.cpuNum()>0) return;
  double fx=0, fy=0, fz=0; //write the total froce.
  for(int i = 0; i < nNodes; ++i) {
    fx += F[i][0];  
    fy += F[i][1]; 
    fz += F[i][2];
  }
//  std::cout << "Total force (from AERO-F): " << fx << " " << fy << " " << fz << std::endl;

}

//------------------------------------------------------------------------------

