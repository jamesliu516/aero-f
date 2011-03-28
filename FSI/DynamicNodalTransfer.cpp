/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel, kevin
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
//  com.fprintf(stderr,"fscale = %e, XScale = %e, tScale = %e.\n", fScale, XScale, tScale);

{  Communication::Window<double> window(com, 1*sizeof(double), &dts);
  window.fence(true);
  structure.sendTimeStep(&window);
  window.fence(false);
}

{  Communication::Window<double> window2(com, 1*sizeof(double), &tMax);
  window2.fence(true);
  structure.sendMaxTime(&window2);
  window2.fence(false);
}

  algNum = structure.getAlgorithmNumber();
//  com.fprintf(stderr,"--- Received initial structure Time-step: %e, Final Time: %e\n", dts, tMax);
  dts /= tScale;
  tMax /= tScale;
  com.barrier();

  wintime  = 0;
  winForce = 0;
  winDisp  = 0;
  dt_tmax  = 0;
  UandUdot = 0;
}

//------------------------------------------------------------------------------

DynamicNodalTransfer::~DynamicNodalTransfer() {
  if(wintime)  delete   wintime;
  if(winForce) delete   winForce;
  if(winDisp)  delete   winDisp;
  if(UandUdot) delete[] UandUdot;
  if(dt_tmax)  delete[] dt_tmax;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendForce() {

  std::pair<double *, int> embedded = structure.getTargetData();
  double *embeddedData = embedded.first;
  int length = embedded.second;

  if(!winForce) winForce = new Communication::Window<double> (com, 3*length*sizeof(double), embeddedData);
    //(KW) WARNING: if "length" changes, winForce needs to be re-initialized !! 

  com.barrier(); //for timing purpose
  winForce->fence(true);
  winForce->accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  winForce->fence(false);

  structure.processReceivedForce();
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateInfo() {
  if(!dt_tmax) dt_tmax = new double[2];
  if(!wintime) wintime = new Communication::Window<double> (com, 2*sizeof(double), (double*)dt_tmax);
  com.barrier(); //for timing purpose
  wintime->fence(true);
  structure.sendInfo(wintime);
  wintime->fence(false);
  dts   = dt_tmax[0];
  tMax  = dt_tmax[1];
//  fprintf(stderr,"*** CPU %d received dts = %e, tMax = %e (dimensional)\n", com.cpuNum(), dts, tMax);
  dts  /= tScale;
  tMax /= tScale;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::getDisplacement(SVec<double,3>& structU, SVec<double,3>& structUdot) {

  if(!UandUdot) 
    UandUdot = new double[2*structU.size()*3];

  if(!winDisp) winDisp = new  Communication::Window<double> (com, 2*3*structU.size()*sizeof(double), (double *)UandUdot);
    //(KW) WARNING: if "structU.size" changes, UandUdot and winDisp needs to be re-initialized !! 

  com.barrier(); //for timing purpose
  winDisp->fence(true);
  structure.sendDisplacement(winDisp);
  winDisp->fence(false);

  for(int i=0; i<structU.size(); i++) 
    for(int j=0; j<3; j++) {
      structU[i][j] = UandUdot[i*3+j];
      structUdot[i][j] = UandUdot[(structU.size()+i)*3+j];
    }
//  com.fprintf(stderr,"norm of received disp/velo = %.12e/%.12e.\n", structU.norm(), structUdot.norm());

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
                                             tScale(iod.ref.rv.time), XScale(iod.ref.rv.tlength),
                                             UScale(iod.ref.rv.tvelocity),  nNodes(0), nElems(0), X(0), Tria(0), U(0),
                                             Udot(0), UandUdot(0), F(0), it(0), structExc(0), algNum(6) //A6
{
  timer = tim;

  // ----------------------------------
  //       User-Provided Info
  // ----------------------------------
  if(iod.problem.type[ProblemData::AERO])
    coupled = true;
  else if(iod.problem.type[ProblemData::FORCED])
    coupled = false;
  else {
    com.fprintf(stderr,"ERROR: Simulation type is not supported by the embedded framework.\n");
    exit(-1);
  }

  // ---- input files ----
  int sp = strlen(iod.input.prefix) + 1;
  meshFile = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(meshFile,"%s%s", iod.input.prefix, iod.input.embeddedSurface);
  matcherFile = new char[sp + strlen(iod.input.match)];
  sprintf(matcherFile,"%s%s", iod.input.prefix, iod.input.match); 

  getSurfFromFEM = coupled && (!strlen(iod.input.embeddedSurface));
  if(getSurfFromFEM)
    com.fprintf(stderr,"- Using the embedded surface provided by structure code.\n");

  // ---- for debug ----
  dim2Treatment = (iod.embed.dim2Treatment == EmbeddedFramework::YES) ? true : false; //by default it's false
  oneWayCoupling = (iod.embed.coupling == EmbeddedFramework::ONEWAY) ? true : false; //by default it's false

  // ---- for forced-motion only ------
  if(!coupled)
    if(iod.forced.type==ForcedData::HEAVING)
      mode = 1;
    else {
      com.fprintf(stderr,"ERROR: Forced motion type is not supported by the embedded framework.\n");
      exit(-1);
    }
  else
    mode = -1;

  // NOTE: All variables stored in EmbeddedStructure must be dimensional! (unless the simulation itself is non-dim)
  tMax  = tScale*iod.ts.maxTime; //iod.ts.maxTime is already non-dimensionalized in IoData
  dt    = tScale*iod.forced.timestep;
  t0    = tScale*iod.restart.etime;
  omega = 2.0*acos(-1.0)*(1.0/tScale)*iod.forced.frequency;
  dx    = iod.ref.rv.length*iod.forced.hv.ax; //tlength = length / aero.displacementScaling;
  dy    = iod.ref.rv.length*iod.forced.hv.ay;
  dz    = iod.ref.rv.length*iod.forced.hv.az;  
/*  com.fprintf(stderr,"*** User Specified Embedded Structure Info ***\n");
  com.fprintf(stderr,"  coupled = %d,  dim2Treatment = %d\n", coupled, dim2Treatment);
  if(!coupled){
    com.fprintf(stderr,"  (forced motion) mode = %d, t0 = %e, tMax = %e, dt = %e\n", mode, t0, tMax, dt);
    com.fprintf(stderr,"                  omega = %e, dx = %e, dy = %e, dz = %e.\n", omega, dx, dy, dz);
  }*/
  // ----------------------------------
  //               End
  // ----------------------------------

  // ---------------------------------
  //       F-S Communication
  // ---------------------------------
  if(coupled) {
    MatchNodeSet **mns = new MatchNodeSet *[1];
    if(com.cpuNum() == 0 && !getSurfFromFEM)
      mns[0] = new MatchNodeSet(matcherFile);
    else 
      mns[0] = new MatchNodeSet();
    structExc = new StructExc(iod, mns, 6, &strCom, &com, 1);

    if(getSurfFromFEM) {//receive embedded surface from FEM
      structExc->getEmbeddedWetSurfaceInfo(nNodes, nElems);
      X = new (com) double[nNodes][3];
      Tria = new (com) int[nElems][3];
      structExc->getEmbeddedWetSurface(nNodes, (double*)X, nElems, (int*)Tria);
      if(com.cpuNum()==0) {
        mns[0]->autoInit(nNodes);
        structExc->updateMNS(mns);
      }
    }


    structExc->negotiate();
    structExc->getInfo();
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
    algNum = structExc->getAlgorithmNumber();
  }
  // ----------------------------------
  //               End
  // ----------------------------------


  // ----------------------------------
  //    Load Structure Mesh (at t=0)
  // ----------------------------------
  if(!getSurfFromFEM) {//otherwise mesh is already loaded.
    // load structure nodes (from file).
    FILE *nodeFile = 0;
    nodeFile = fopen(meshFile,"r");
    if(!nodeFile) com.fprintf(stderr,"ERROR: Embedded surface mesh file could not be found!\n");
    char c1[200], c2[200];
    int num0 = 0, num1 = 0, count, nInputs;
    double x1,x2,x3;
    int toto = fscanf(nodeFile, "%s %s\n", c1, c2);
    char debug[6]="Nodes";
    for (int i=0; i<5; i++) 
      if(debug[i]!=c1[i]) {com.fprintf(stderr,"ERROR: The embedded surface file (%s) must begin with keyword `Nodes'!\n", meshFile); exit(-1);}
  
    std::list<Vec3D> nodeList;
    std::list<int> indexList;
    std::list<Vec3D>::iterator it1;
    std::list<int>::iterator it2;
    int maxIndex = 0;

    while(1) {
      nInputs = fscanf(nodeFile,"%s", c1);
      if(nInputs!=1) break;
      if(c1[0]=='E') //done with the node set
        break;
      num1 = atoi(c1);
      if(num1<1) {com.fprintf(stderr,"ERROR: detected a node with index %d in the embedded surface file!\n",num1); exit(-1);}
      indexList.push_back(num1);
      if(num1>maxIndex)
        maxIndex = num1;

      toto = fscanf(nodeFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList.push_back(Vec3D(x1,x2,x3));
    }
    nNodes = nodeList.size();
    if(nNodes != maxIndex) {
      com.fprintf(stderr,"ERROR: The node set of the embedded surface have gap(s). \n");
      com.fprintf(stderr,"       Detected max index = %d, number of nodes = %d\n", maxIndex, nNodes);
      com.fprintf(stderr,"NOTE: Currently the node set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
      exit(-1);
    }
    X = new (com) double[nNodes][3];

    it2=indexList.begin();
    for (it1=nodeList.begin(); it1!=nodeList.end(); it1++) {
      X[(*it2)-1][0] = (*it1)[0]; 
      X[(*it2)-1][1] = (*it1)[1]; 
      X[(*it2)-1][2] = (*it1)[2];
      it2++; 
    }
    fclose(nodeFile);
  }

  // ----------------------------------
  //               End
  // ----------------------------------


  // allocate memory for other stuff...
  U = new (com) double[nNodes][3];
  Udot = new (com) double[nNodes][3];
  UandUdot = new (com) double[nNodes*2][3];
  F = new (com) double[nNodes][3];

  // prepare distinfo for struct nodes
  if(coupled) {
    int *locToGlob = new int[1];
    locToGlob[0] = 0;
    di = new DistInfo(1, 1, 1, locToGlob, &com);
    di->setLen(0,nNodes);
    di->finalize(false);
  }

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

}

//------------------------------------------------------------------------------

EmbeddedStructure::~EmbeddedStructure()
{
  if(X)         delete[] X;
  if(Tria)      delete[] Tria;
  if(U)         delete[] U;
  if(Udot)      delete[] Udot;
  if(UandUdot)  delete[] UandUdot;
  if(F)         delete[] F;
  if(structExc) delete structExc;  
  delete[] meshFile;
  delete[] matcherFile;
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
EmbeddedStructure::sendInfo(Communication::Window<double> *window)
{
  if(coupled)
    structExc->getInfo();
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  if(coupled) {
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
  } /* else: nothing to be done */

  dt_tmax[0] = dt; dt_tmax[1] = tMax;
  for(int i = 0; i < com.size(); ++i)
    window->put((double*)dt_tmax, 0, 2, i, 0);
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
    double time;
    if(it==1 && algNum==6)
      time = t0 + 0.5*dt*(double)it;
    else
      time = t0 + dt*(double)it;
     
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

  for(int i = 0; i < com.size(); ++i) 
    window->put((double*)UandUdot, 0, 2*3*nNodes, i, 0);
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
/*
  for(int i=0; i<nNodes; ++i) 
    fprintf(stderr,"%d %e %e %e\n", i+1, F[i][0], F[i][1], F[i][2]);
  sleep(2);
*/
//  std::cout << "Total force (from AERO-F): " << fx << " " << fy << " " << fz << std::endl;

}

//------------------------------------------------------------------------------

