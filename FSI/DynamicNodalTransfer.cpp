/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel, kevin
 */
#include <iostream>
#include <FSI/DynamicNodalTransfer.h>
#include <FSI/CrackingSurface.h>
#include <IoData.h>
#include <Vector3D.h>
#include <MatchNode.h>
#include <StructExc.h>
#include <DistVector.h>
#include <assert.h>
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

  //initialize windows
  dt_tmax = new double[2];
  wintime = new Communication::Window<double> (com, 2*sizeof(double), (double*)dt_tmax);

  XandUdot = new double[2*3*structure.totalNodes];
  winDisp = new  Communication::Window<double> (com, 2*3*structure.nNodes*sizeof(double), (double *)XandUdot);

  std::pair<double *, int> embedded = structure.getTargetData();
  double *embeddedData = embedded.first;
  int length = embedded.second;
  winForce = new Communication::Window<double> (com, 3*length*sizeof(double), embeddedData);

  //get structure position
  com.barrier(); //for timing purpose
  winDisp->fence(true);
  structure.sendInitialPosition(winDisp);
  winDisp->fence(false);

  int N = structure.nNodes;
  for(int i=0; i<N; i++)
    for(int j=0; j<3; j++) {
      XandUdot[i*3+j] *= 1.0/XScale;
      XandUdot[(N+i)*3+j] *= 1.0/UScale;
    }

  structureSubcycling = (algNum==22) ? getStructSubcyclingInfo() : 0;
    
}

//------------------------------------------------------------------------------

DynamicNodalTransfer::~DynamicNodalTransfer() {
  if(wintime)  delete   wintime;
  if(winForce) delete   winForce;
  if(winDisp)  delete   winDisp;
  if(XandUdot) delete[] XandUdot;
  if(dt_tmax)  delete[] dt_tmax;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendForce()
{
  structure.clearForceVector();

  com.barrier(); //for timing purpose
  winForce->fence(true);
  winForce->accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  winForce->fence(false);

  structure.processReceivedForce();
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendFluidSuggestedTimestep(double dtf0)
{
  dtf0 *= tScale;
  structure.sendFluidSuggestedTimestep(dtf0);
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateInfo() {
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

int
DynamicNodalTransfer::getNewCracking()
{
  if(!cracking()) {
    com.fprintf(stderr,"WARNING: Cracking is not considered in the structure simulation!\n");
    return 0;
  }
  int numConnUpdate = structure.getNewCracking();

  if(numConnUpdate) { //update "windows".
    int N = structure.nNodes;
    if((N*2*3*sizeof(double) == winDisp->size())) {com.fprintf(stderr,"WEIRD!\n");exit(-1);}

    delete winDisp;
    winDisp = new Communication::Window<double> (com, 2*3*N*sizeof(double), (double *)XandUdot);
   
    delete winForce;
    std::pair<double *, int> embedded = structure.getTargetData();
    double *embeddedData = embedded.first;
    int length = embedded.second;
    if(length!=N) {com.fprintf(stderr,"WEIRD TOO!\n");exit(-1);}
    winForce = new Communication::Window<double> (com, 3*length*sizeof(double), embeddedData);
  }

  return numConnUpdate;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::getDisplacement()
{
  int N = numStNodes();
  if(N*2*3*sizeof(double) != winDisp->size()) {
    com.fprintf(stderr,"SOFTWARE BUG: length of winDisp (i.e. # nodes) is %d (should be %d)!\n", 
                        winDisp->size()/(2*3*sizeof(double)), N);
    exit(-1);
  }

  com.barrier(); //for timing purpose
  winDisp->fence(true);
  structure.sendDisplacement(winDisp);
  winDisp->fence(false);

  for(int i=0; i<N; i++) { 
    for(int j=0; j<3; j++) {
      XandUdot[i*3+j] *= 1.0/XScale;
      XandUdot[(N+i)*3+j] *= 1.0/UScale;
    }
  }
}

//------------------------------------------------------------------------------

int DynamicNodalTransfer::getStructSubcyclingInfo()
{
  int subcyc = structure.sendSubcyclingInfo();
  return subcyc;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateOutputToStructure(double dt, double dtLeft, SVec<double,3> &fs)
{
  if(F.size() != fs.size()) {
    com.fprintf(stderr,"force vector in DynamicNodalTransfer resized (from %d to %d)!\n", F.size(), fs.size());
    F.resize(fs.size());
  }
  F = fs;
}

//------------------------------------------------------------------------------

EmbeddedStructure::EmbeddedStructure(IoData& iod, Communicator &comm, Communicator &strCom, Timer *tim) : com(comm), 
                                             tScale(iod.ref.rv.time), XScale(iod.ref.rv.tlength),
                                             UScale(iod.ref.rv.tvelocity),  nNodes(0), nElems(0), elemType(3),
                                             cracking(0), X(0), X0(0), Tria(0), U(0), Udot(0), XandUdot(0), F(0), it(0), 
                                             structExc(0), mns(0), algNum(6) //A6
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
  if(!coupled) {
    if(iod.forced.type==ForcedData::HEAVING)
      mode = 1;
    else if(iod.forced.type==ForcedData::PITCHING)
      mode = 2;
    else {
      com.fprintf(stderr,"ERROR: Forced motion type is not supported by the embedded framework.\n");
      exit(-1);
    }
  } else mode = -1;

  // NOTE: All variables stored in EmbeddedStructure must be dimensional! (unless the simulation itself is non-dim)
  tMax  = tScale*iod.ts.maxTime; //iod.ts.maxTime is already non-dimensionalized in IoData
  dt    = tScale*iod.forced.timestep;
  t0    = tScale*iod.restart.etime;
  omega = 2.0*acos(-1.0)*(1.0/tScale)*iod.forced.frequency;

  // for heaving
  dx    = iod.ref.rv.length*iod.forced.hv.ax; //tlength = length / aero.displacementScaling;
  dy    = iod.ref.rv.length*iod.forced.hv.ay;
  dz    = iod.ref.rv.length*iod.forced.hv.az;  

  // for pitching
  alpha_in  = (acos(-1.0)*iod.forced.pt.alpha_in) / 180.0;  // initial angle of rotation
  alpha_max = (acos(-1.0)*iod.forced.pt.alpha_max) / 180.0;  // maximum angle of rotation
  x1[0] = iod.forced.pt.x1;  x1[1] = iod.forced.pt.y1;  x1[2] = iod.forced.pt.z1;
  x2[0] = iod.forced.pt.x2;  x2[1] = iod.forced.pt.y2;  x2[2] = iod.forced.pt.z2;
  for(int i=0; i<3; i++) { //get back to user-specified coordinates.
    x1[i] *= iod.ref.rv.length;
    x2[i] *= iod.ref.rv.length;
  }
  u = x2[0]-x1[0];  v = x2[1]-x1[1];  w = x2[2]-x1[2];
  // unit normals of axis of rotation //
  ix = u/sqrt(u*u+v*v+w*w);  iy = v/sqrt(u*u+v*v+w*w);  iz = w/sqrt(u*u+v*v+w*w);

  // ----------------------------------
  //               End
  // ----------------------------------

  // ---------------------------------
  //       F-S Communication
  // ---------------------------------
  if(coupled) {
    mns = new MatchNodeSet *[1];
    if(com.cpuNum() == 0 && !getSurfFromFEM)
      mns[0] = new MatchNodeSet(matcherFile);
    else 
      mns[0] = new MatchNodeSet();
    structExc = new StructExc(iod, mns, 6, &strCom, &com, 1);

    if(getSurfFromFEM) {//receive embedded surface from the structure code.
      bool crack;
      int  nStNodes, nStElems, totalStNodes, totalStElems;
      structExc->getEmbeddedWetSurfaceInfo(elemType, crack, nStNodes, nStElems); 

      // initialize cracking information
      if(crack) {
        structExc->getInitialCrackingSetup(totalStNodes, totalStElems);
        cracking = new CrackingSurface(elemType, nStElems, totalStElems, nStNodes, totalStNodes);
      } else {
        totalStNodes = nStNodes;
        totalStElems = nStElems;
      }

      // allocate memory for the node list
      totalNodes = totalStNodes;
      nNodes     = nStNodes;
      X = new (com) double[totalNodes][3];
      X0 = new (com) double[totalNodes][3];
      for(int i=0; i<totalNodes; i++)
        X[i][0] = X[i][1] = X[i][2] = X0[i][0] = X0[i][1] = X0[i][2] = 0.0;

      // allocate memory for the element topology list
      int tmpTopo[nStElems][4];
      switch (elemType) {
        case 3: //all triangles
          totalElems = totalStElems;
          nElems     = nStElems; 
          Tria = new (com) int[totalElems][3];
          structExc->getEmbeddedWetSurface(nNodes, (double*)X0, nElems, (int*)Tria, elemType);
          break;
        case 4: //quadrangles include triangles represented as degenerated quadrangles.
          structExc->getEmbeddedWetSurface(nNodes, (double*)X0, nStElems, (int*)tmpTopo, elemType);
          if(cracking) {
            totalElems = totalStElems*2;
            Tria = new (com) int[totalElems][3];
            nElems = cracking->splitQuads((int*)tmpTopo, nStElems, Tria);
          } else
            splitQuads((int*)tmpTopo, nStElems); //memory for Tria will be allocated
          break;
        default:
          com.fprintf(stderr,"ERROR: Element type (%d) of the wet surface not recognized! Must be 3 or 4.\n", elemType);
          exit(-1);
      }

      for(int i=0; i<nNodes; i++)
        for(int j=0; j<3; j++)
          X[i][j] = X0[i][j];

      if(com.cpuNum()==0) {
        mns[0]->autoInit(totalNodes); //in case of cracking, match all the nodes including inactive ones.
        structExc->updateMNS(mns);
      }
    }

    structExc->negotiate();
    structExc->getInfo();
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
    algNum = structExc->getAlgorithmNumber();

    if(cracking) {
      getInitialCrack();
      cracking->setNewCrackingFlag(false);
    }
  }

  // ----------------------------------
  //               End
  // ----------------------------------


  // ----------------------------------
  //    Load Structure Mesh (at t=0)
  // ----------------------------------
  if(!getSurfFromFEM) {//otherwise mesh is already loaded.
    // load structure nodes and elements (from file).
    FILE *topFile = 0;
    topFile = fopen(meshFile,"r");
    if(!topFile) com.fprintf(stderr,"ERROR: Embedded surface mesh file could not be found!\n");
    char c1[200], c2[200], c3[200];
    int num0 = 0, num1 = 0, count, nInputs;
    double x1,x2,x3;
    int toto = fscanf(topFile, "%s %s\n", c1, c2);
    char debug[6]="Nodes";
    for (int i=0; i<5; i++) 
      if(debug[i]!=c1[i]) {com.fprintf(stderr,"ERROR: The embedded surface file (%s) must begin with keyword `Nodes'!\n", meshFile); exit(-1);}
  
    std::list<Vec3D> nodeList;
    std::list<int> indexList;
    std::list<Vec3D>::iterator it1;
    std::list<int>::iterator it2;
    int maxIndex = 0;

    while(1) {
      nInputs = fscanf(topFile,"%s", c1);
      if(nInputs!=1) break;
      if(c1[0]=='E') //done with the node set
        break;
      num1 = atoi(c1);
      if(num1<1) {com.fprintf(stderr,"ERROR: detected a node with index %d in the embedded surface file!\n",num1); exit(-1);}
      indexList.push_back(num1);
      if(num1>maxIndex)
        maxIndex = num1;

      toto = fscanf(topFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList.push_back(Vec3D(x1,x2,x3));
    }
    nNodes = totalNodes = nodeList.size();
    if(nNodes != maxIndex) {
      com.fprintf(stderr,"ERROR: The node set of the embedded surface have gap(s). \n");
      com.fprintf(stderr,"       Detected max index = %d, number of nodes = %d\n", maxIndex, nNodes);
      com.fprintf(stderr,"NOTE: Currently the node set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
      exit(-1);
    }
    X0 = new (com) double[totalNodes][3];
    X  = new (com) double[totalNodes][3];

    it2=indexList.begin();
    for (it1=nodeList.begin(); it1!=nodeList.end(); it1++) {
      X[(*it2)-1][0] = X0[(*it2)-1][0] = (*it1)[0]; 
      X[(*it2)-1][1] = X0[(*it2)-1][1] = (*it1)[1]; 
      X[(*it2)-1][2] = X0[(*it2)-1][2] = (*it1)[2];
      it2++; 
    }

    // now load the elements
    if(nInputs!=1) {
      com.fprintf(stderr,"ERROR: Failed reading embedded surface from file: %s\n", meshFile); exit(-1);}
    int nm = fscanf(topFile,"%s %s %s\n", c1,c2,c3);
    char debug2[6] = "using";
    for (int i=0; i<5; i++)
      if(debug2[i]!=c2[i]) {com.fprintf(stderr,"ERROR: Failed reading embedded surface from file: %s\n", meshFile); exit(-1);}

    std::list<int> elemIdList;
    std::list<int> elemList1;
    std::list<int> elemList2;
    std::list<int> elemList3;
    std::list<int>::iterator it_0;
    std::list<int>::iterator it_1;
    std::list<int>::iterator it_2;
    std::list<int>::iterator it_3;
    int node1, node2, node3;
    maxIndex = -1;

    while(1) {
      nInputs = fscanf(topFile,"%d", &num0);
      if(nInputs!=1) break;
      toto = fscanf(topFile,"%d %d %d %d\n", &num1, &node1, &node2, &node3);
      if(num0<1) {com.fprintf(stderr,"ERROR: Detected an element with Id %d in the embedded surface (%s)!\n", num0, meshFile); exit(-1);}
      elemIdList.push_back(num0-1);  //start from 0.
      elemList1.push_back(node1-1);
      elemList2.push_back(node2-1);
      elemList3.push_back(node3-1);
      if(num0-1>maxIndex)
        maxIndex = num0-1;
    }
    nElems = totalElems = elemList1.size();
    if(nElems != maxIndex+1) {
      com.fprintf(stderr,"ERROR: The element set of the embedded surface have gap(s). \n");
      com.fprintf(stderr,"       Detected max index = %d, number of elements = %d\n", maxIndex+1, nElems);
      com.fprintf(stderr,"NOTE: Currently the element set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
      exit(-1);
    }

    Tria = new (com) int[totalElems][3];

    it_0 = elemIdList.begin();
    it_1 = elemList1.begin();
    it_2 = elemList2.begin();
    it_3 = elemList3.begin();
    for (int i=0; i<nElems; i++) {
      Tria[*it_0][0] = *it_1;
      Tria[*it_0][1] = *it_2;
      Tria[*it_0][2] = *it_3;
      it_0++;
      it_1++;
      it_2++;
      it_3++;
    }

    fclose(topFile);
  }

  // ----------------------------------
  //               End
  // ----------------------------------

  // allocate memory for other stuff...
  U = new (com) double[totalNodes][3];
  Udot = new (com) double[totalNodes][3];
  XandUdot = new (com) double[totalNodes*2][3];
  F = new (com) double[totalNodes][3];

  // prepare distinfo for struct nodes
  if(coupled) {
    int *locToGlob = new int[1];
    locToGlob[0] = 0;
    di = new DistInfo(1, 1, 1, locToGlob, &com);
    di->setLen(0,totalNodes);
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
/* PJSA
  if(X)         delete[] X;
  if(Tria)      delete[] Tria;
  if(U)         delete[] U;
  if(Udot)      delete[] Udot;
  if(XandUdot)  delete[] XandUdot;
  if(F)         delete[] F;
*/
  if(structExc) delete structExc;  
  if(mns)      {delete mns[0]; delete [] mns;}
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
EmbeddedStructure::clearForceVector()
{
  //clear the force.
  for (int i=0; i<nNodes; i++)
    F[i][0] = F[i][1] = F[i][2] = 0.0;
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
EmbeddedStructure::sendInitialPosition(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  for(int i=0; i<nNodes; i++)
    for(int j=0; j<3; j++) {
      XandUdot[i][j] = X[i][j];
      XandUdot[(i+nNodes)][j] = Udot[i][j];
    }

  for(int i = 0; i < com.size(); ++i)
    window->put((double*)XandUdot, 0, 2*3*nNodes, i, 0);
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
     Y = Y0; //KW: as long as Y = Y0, it doesn't matter if Y0 is X or X0, or anything else...
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
    else if (mode==2) { //pitching
      double theta = alpha_in + alpha_max*sin(omega*time);
      double costheta = cos(theta);
      double sintheta = sin(theta);
      for(int i=0; i<nNodes; ++i) {
        U[i][0] = U[i][1] = U[i][2] = 0.0;
        double p[3];
        for(int j=0; j<3; j++)
          p[j] = X0[i][j] - x1[j];

        U[i][0] += (costheta + (1 - costheta) * ix * ix) * p[0];
        U[i][0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];
        U[i][0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];

        U[i][1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];
        U[i][1] += (costheta + (1 - costheta) * iy * iy) * p[1];
        U[i][1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];

        U[i][2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];
        U[i][2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];
        U[i][2] += (costheta + (1 - costheta) * iz * iz) * p[2];

        U[i][0] += x1[0];
        U[i][1] += x1[1];
        U[i][2] += x1[2];

        if(it==1) {
          if(algNum==6)
            for(int j=0; j<3; j++)
              Udot[i][j] = (U[i][j]-X0[i][j])/(0.5*dt);
          else
            for(int j=0; j<3; j++)
              Udot[i][j] = (U[i][j]-X0[i][j])/dt;
        } else {
          for(int j=0; j<3; j++)
            Udot[i][j] = (U[i][j]-X[i][j])/dt;
        }

        for(int j=0; j<3; j++)
          U[i][j] -= X0[i][j];
      }
    }
    else if (mode==3) //heaving with a constant velocity (in this case dx dy dz are velocity
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

  for(int i=0; i<nNodes; i++) {
//    fprintf(stderr,"U %d %e %e %e\n", i+1, U[i][0], U[i][1], U[i][2]);
//    fprintf(stderr,"Udot %d %e %e %e\n", i+1, Udot[i][0], Udot[i][1], Udot[i][2]);
    for(int j=0; j<3; j++) {
      X[i][j] = X0[i][j] + U[i][j];
      XandUdot[i][j] = X[i][j];
      XandUdot[(i+nNodes)][j] = Udot[i][j];
    }
  }

  for(int i = 0; i < com.size(); ++i) 
    window->put((double*)XandUdot, 0, 2*3*nNodes, i, 0);
}

//------------------------------------------------------------------------------

int
EmbeddedStructure::sendSubcyclingInfo(/*Communication::Window<int> *window*/)
{
/*  int subcyc = structExc->getSubcyclingInfo();
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  for(int i = 0; i < com.size(); ++i)
    window->put((int*)&subcyc, 0, 1, i, 0);
*/
  return structExc->getSubcyclingInfo();
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


/*
  if(com.cpuNum()>0) return;
  double fx=0, fy=0, fz=0; //write the total froce.
  for(int i = 0; i < nNodes; ++i) {
    fx += F[i][0];  
    fy += F[i][1]; 
    fz += F[i][2];
  }
  for(int i=0; i<nNodes; ++i) 
    fprintf(stderr,"%d %e %e %e\n", i+1, F[i][0], F[i][1], F[i][2]);
  sleep(2);
*/
//  std::cout << "Total force (from AERO-F): " << fx << " " << fy << " " << fz << std::endl;

}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendFluidSuggestedTimestep(double dtf0)
{
  structExc->sendFluidSuggestedTimestep(dtf0);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::splitQuads(int* quads, int nStElems)
{ 
  int nTrias = 0;
  for(int i=0; i<nStElems; i++)
    if(quads[i*4+2]==quads[i*4+3])
      nTrias += 1;
    else 
      nTrias += 2;

  Tria = new (com) int[nTrias][3];

  int count = 0;
  for(int i=0; i<nStElems; i++) { 
    Tria[count][0] = quads[i*4];
    Tria[count][1] = quads[i*4+1];
    Tria[count][2] = quads[i*4+2];
    count++;

    if(quads[i*4+2]==quads[i*4+3])
      continue;

    Tria[count][0] = quads[i*4];
    Tria[count][1] = quads[i*4+2];
    Tria[count][2] = quads[i*4+3];
    count++;
  }

  if(count!=nTrias) {com.fprintf(stderr,"Software bug in FSI/EmbeddedStructure/splitQuad!\n");exit(-1);}
  nElems = totalElems = nTrias;
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::getInitialCrack()
{
  int newNodes, numConnUpdate, numLSUpdate;
  bool need2update = structExc->getNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); //inputs will be modified.
  if(!need2update) return; //Nothing new :)

  // get initial phantom nodes.
  structExc->getInitialPhantomNodes(newNodes,X,nNodes);
    //NOTE: nNodes will be updated in "getNewCracking"

  // get initial phantom elements (topo change).
  getNewCracking(numConnUpdate, numLSUpdate, newNodes);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::getNewCracking(int numConnUpdate, int numLSUpdate, int newNodes)
{
  if(numConnUpdate<1)  return;

  int phantElems[5*numConnUpdate]; // elem.id and node id.
  double phi[4*numLSUpdate];
  int phiIndex[numLSUpdate];
  int new2old[newNodes*2];
 
  structExc->getNewCracking(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, new2old, newNodes); 

  if(elemType!=4) {com.fprintf(stderr,"ERROR: only support quadrangles for cracking!\n");exit(1);} 
  nNodes += newNodes;
  nElems += cracking->updateCracking(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, Tria, nNodes, new2old, newNodes);
  if(nElems!=cracking->usedTrias()) {
    com.fprintf(stderr,"ERROR: inconsistency in the number of used triangles. (Software bug.)\n");exit(-1);}

  if(com.cpuNum()==0)
    structExc->updateNumStrNodes(nNodes); //set numStrNodes to nNodes in structExc and mns
}

//------------------------------------------------------------------------------

int
EmbeddedStructure::getNewCracking()
{
  int newNodes, numConnUpdate, numLSUpdate;
  bool need2update = structExc->getNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); //inputs will be modified.
  if(!need2update) {assert(numConnUpdate==0); return 0;} //Nothing new :)
  getNewCracking(numConnUpdate, numLSUpdate, newNodes);
  return numConnUpdate;
}

//------------------------------------------------------------------------------

