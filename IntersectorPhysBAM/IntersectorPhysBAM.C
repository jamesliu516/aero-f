#include <iostream>
#include <stdio.h>
#include "IntersectorPhysBAM.h"
#include "FloodFill.h"
#include "Mpi_Utilities.h"
#include "Vector3D.h"
#include "Communicator.h"
#include <Domain.h>
#include <IoData.h>
#include <Vector.h>
#include <DistVector.h>
#include <Timer.h>
#include "LevelSet/LevelSetStructure.h"
#include "parser/Assigner.h"
#include <Connectivity.h>
#include <queue>

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>

using std::pair;
using std::map;
using std::list;

using PhysBAM::TRIPLE;
using PhysBAM::GLOBAL_SUBD_ID;

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

const int IntersectorPhysBAM::UNDECIDED, IntersectorPhysBAM::INSIDE, IntersectorPhysBAM::OUTSIDE;

//----------------------------------------------------------------------------

DistIntersectorPhysBAM::DistIntersectorPhysBAM(IoData &iod, Communicator *comm)
{
  this->numFluid = iod.eqs.numPhase;
  floodFill=new FloodFill();
  com = comm;
  com->fprintf(stderr,"Using Intersector PhysBAM\n");

  //get embedded structure surface mesh and restart pos
  char *struct_mesh, *struct_restart_pos;
  int sp = strlen(iod.input.prefix) + 1;

  struct_mesh        = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(struct_mesh,"%s%s", iod.input.prefix, iod.input.embeddedSurface);
  struct_restart_pos = new char[sp + strlen(iod.input.positions)];
  if(iod.input.positions[0] != 0)
    sprintf(struct_restart_pos,"%s%s", iod.input.prefix, iod.input.positions);
  else //no restart position file provided
    struct_restart_pos = ""; 
  interpolatedNormal = (iod.embed.structNormal==EmbeddedFramework::NODE_BASED) ? 
                        true : false;
  
  //nodal or facet normal?

  //initialize the following to 0(NULL)
  physInterface = 0;
  triNorms = 0;
  triSize = 0;
  nodalNormal = 0;
  status = 0;
  status0 = 0;
  occluded_node = 0;
  swept_node = 0;
  pseudoPhi = 0;
  boxMin = 0;
  boxMax = 0;

  //Load files. Compute structure normals. Initialize PhysBAM Interface
  init(struct_mesh, struct_restart_pos);
  comm->barrier();

  delete[] struct_mesh;
  if(iod.input.positions[0] != 0) delete[] struct_restart_pos;
}

//----------------------------------------------------------------------------

DistIntersectorPhysBAM::~DistIntersectorPhysBAM() 
{
  if(Xs)          delete[] Xs;
  if(Xs0)         delete[] Xs0;
  if(Xs_n)        delete[] Xs_n;
  if(Xs_np1)      delete[] Xs_np1;
  if(Xsdot)       delete[] Xsdot;
  if(status)      delete   status;
  if(status0)     delete   status0;
  if(triSize)     delete[] triSize;
  if(triNorms)    delete[] triNorms;
  if(nodalNormal) delete[] nodalNormal;
  if(boxMax)      delete   boxMax;
  if(boxMin)      delete   boxMin;

  for(int i = 0; i < numLocSub; ++i) {
      delete intersector[i];}
  delete []intersector;
  delete physInterface;
  delete floodFill;
}

//----------------------------------------------------------------------------

LevelSetStructure &
DistIntersectorPhysBAM::operator()(int subNum) const {
  return *intersector[subNum];
}

//----------------------------------------------------------------------------

/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorPhysBAM::init(char *solidSurface, char *restartSolidSurface) {

  // Read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(solidSurface, "r");
  if (topFile == NULL) {com->fprintf(stderr, "Embedded structure surface mesh doesn't exist: %s, %s\n",solidSurface,restartSolidSurface); exit(1); }

  // load the nodes and initialize all node-based variables.

  // load solid nodes at t=0
  char c1[200], c2[200], c3[200];
  int num0 = 0, num1 = 0, nInputs;
  double x1,x2,x3;
  fscanf(topFile, "%s %s\n", c1, c2);
  char debug[6]="Nodes";
  for (int i=0; i<5; i++)
    if(debug[i]!=c1[i]) {fprintf(stderr,"Failed in reading file: %s\n", solidSurface); exit(-1);}

  std::list<std::pair<int,Vec3D> > nodeList;
  std::list<std::pair<int,Vec3D> >::iterator it;

  int ndMax = 0;

  while(1) {
    nInputs = fscanf(topFile,"%s", c1);
    if(nInputs!=1) break;
    char *endptr;
    num1 = strtol(c1, &endptr, 10);
    if(endptr == c1) break;

    fscanf(topFile,"%lf %lf %lf\n", &x1, &x2, &x3);
    nodeList.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
    ndMax = std::max(num1, ndMax);
    num0 = num1;
  }
  numStNodes = ndMax;

  // feed data to Xss. 
  Xs      = new Vec3D[numStNodes];
  Xs0     = new Vec3D[numStNodes];
  Xs_n    = new Vec3D[numStNodes];
  Xs_np1  = new Vec3D[numStNodes];
  Xsdot   = new Vec3D[numStNodes];
  solidX  = new Vec<Vec3D>(numStNodes, Xs);
  solidX0 = new Vec<Vec3D>(numStNodes, Xs0);
  solidXn = new Vec<Vec3D>(numStNodes, Xs_n);
  
  for (it=nodeList.begin(); it!=nodeList.end(); it++) 
    Xs[it->first-1] = it->second;

  for (int k=0; k<numStNodes; k++) {
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);
  }

  // load the elements.
  if(nInputs!=1) {
    fprintf(stderr,"Failed reading elements from file: %s\n", solidSurface); exit(-1);}
  fscanf(topFile,"%s %s %s\n", c1,c2,c3);
  char debug2[6] = "using";
  for (int i=0; i<5; i++) 
    if(debug2[i]!=c2[i]) {fprintf(stderr,"Failed in reading file: %s\n", solidSurface); exit(-1);}
    
  std::list<int> elemList1;
  std::list<int> elemList2;
  std::list<int> elemList3;
  std::list<int>::iterator it1;
  std::list<int>::iterator it2;
  std::list<int>::iterator it3;
  int node1, node2, node3;

  while(1) {
    nInputs = fscanf(topFile,"%d", &num0);
    if(nInputs!=1) break;
    fscanf(topFile,"%d %d %d %d\n", &num1, &node1, &node2, &node3);
    elemList1.push_back(node1-1);
    elemList2.push_back(node2-1);
    elemList3.push_back(node3-1);
  }
  numStElems = elemList1.size();

  stElem = new int[numStElems][3];
  
  it1 = elemList1.begin();
  it2 = elemList2.begin();
  it3 = elemList3.begin();
  for (int i=0; i<numStElems; i++) {
    stElem[i][0] = *it1;
    stElem[i][1] = *it2;
    stElem[i][2] = *it3;
    it1++;
    it2++;
    it3++;
  } 

  fclose(topFile);

  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");
    if(resTopFile==NULL) {com->fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
    int ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    while(1) {
      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;    
      char *endptr;
      num1 = strtol(c1, &endptr, 10);
      if(endptr == c1) break;

      fscanf(resTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
      ndMax = std::max(num1, ndMax);
    }
    if (ndMax!=numStNodes) {
      com->fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong.\n");
      exit(1);
    }

    for (int k=0; k<numStNodes; k++)
      Xs[k] = Vec3D(0,0,0);
    
    for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++)
      Xs[it2->first-1] = it2->second;

    for (int k=0; k<numStNodes; k++) {
      Xs_n[k]         = Xs[k];
      Xs_np1[k]    = Xs[k];
    }
    fclose(resTopFile);
  }

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) 
    com->fprintf(stderr,"Ok.\n");
  else 
    exit(-1); 

  initializePhysBAM();
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::initializePhysBAM() { //NOTE: In PhysBAM array index starts from 1 instead of 0
// Initialize the Particles list
  PhysBAM::TRIANGULATED_SURFACE<double>& physbam_triangulated_surface=*PhysBAM::TRIANGULATED_SURFACE<double>::Create();

  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >& physbam_solids_particle = physbam_triangulated_surface.particles;
  physbam_solids_particle.array_collection.Resize(numStNodes);
  for (int i=0; i<numStNodes; i++) 
    physbam_solids_particle.X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0],Xs[i][1], Xs[i][2]);
  
  // Initialize the Triangle list
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > & physbam_stElem=physbam_triangulated_surface.mesh.elements;
  physbam_stElem.Resize(numStElems);
  for (int i=0; i<numStElems; i++)
    physbam_stElem(i+1) = PhysBAM::VECTOR<int,3>(stElem[i][0]+1,stElem[i][1]+1,stElem[i][2]+1);

  // Construct TRIANGLE_MESH triangle_mesh.
  physbam_triangulated_surface.Update_Number_Nodes();
  physbam_triangulated_surface.mesh.Initialize_Adjacent_Elements();

  // Construct TRIANGULATED_SURFACE.
  physbam_triangulated_surface.Update_Triangle_List();
  if(physInterface) delete physInterface;
  physInterface = new PhysBAMInterface<double>(physbam_triangulated_surface);
}

//----------------------------------------------------------------------------

EdgePair DistIntersectorPhysBAM::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

//----------------------------------------------------------------------------

bool DistIntersectorPhysBAM::checkTriangulatedSurface()
{
  map<iipair, ibpair> edgeMap;

  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = stElem[iTriangle][0];  to1 = stElem[iTriangle][1];  found1 = false;
    from2 = stElem[iTriangle][1];  to2 = stElem[iTriangle][2];  found2 = false;
    from3 = stElem[iTriangle][2];  to3 = stElem[iTriangle][0];  found3 = false;

    EdgePair ep[3];
    ep[0] = makeEdgePair(from1, to1, iTriangle);
    ep[1] = makeEdgePair(from2, to2, iTriangle);
    ep[2] = makeEdgePair(from3, to3, iTriangle);

    for(int i=0; i < 3; ++i) {
      map<iipair, ibpair>::iterator it = edgeMap.find(ep[i].first);
      if(it != edgeMap.end()) { // we found this edge
         if(it->second.second == ep[i].second.second)
           {com->fprintf(stderr,"ERROR: surface is not closed or a triangle orientation problem. exit.\n"); return false;}
         else {
             int oTriangle = it->second.first;
             int n1 = it->second.second ? ep[i].first.first : ep[i].first.second;
             int edgeIndex;
             if(stElem[oTriangle][0] == n1)
               edgeIndex = 0;
             else if(stElem[oTriangle][1] == n1)
               edgeIndex = 1;
             else
               edgeIndex = 2;
         }
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
  return true;
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::buildSolidNormals() {
  if(!triNorms) triNorms = new Vec3D[numStElems];
  if(!triSize)  triSize = new double[numStElems];
  if(interpolatedNormal) {
    if(!nodalNormal)
      nodalNormal = new Vec3D[numStNodes];
    for(int i=0; i<numStNodes; i++)
      nodalNormal[i] = 0.0;
  }

  // Also look to determine a point inside the solid but away from the structure.
  double nrmMax = 0;
  int trMaxNorm = -1;
  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int n1 = stElem[iTriangle][0];
    int n2 = stElem[iTriangle][1];
    int n3 = stElem[iTriangle][2];
    double x1 = Xs[n1][0];
    double y1 = Xs[n1][1];
    double z1 = Xs[n1][2];
    double dx2 = Xs[n2][0]-x1;
    double dy2 = Xs[n2][1]-y1;
    double dz2 = Xs[n2][2]-z1;
    double dx3 = Xs[n3][0]-x1;
    double dy3 = Xs[n3][1]-y1;
    double dz3 = Xs[n3][2]-z1;

   // first calculate the "size" of the triangle.
    double dx23 = Xs[n3][0] - Xs[n2][0];
    double dy23 = Xs[n3][1] - Xs[n2][1];
    double dz23 = Xs[n3][2] - Xs[n2][2];
    double size12 = Vec3D(dx2,dy2,dz2).norm();
    double size13 = Vec3D(dx3,dy3,dz3).norm();
    double size23 = Vec3D(dx23,dy23,dz23).norm();
    triSize[iTriangle] = min(size12,min(size13,size23));

    // now calculate the normal.
    triNorms[iTriangle] = Vec3D(dx2, dy2, dz2)^Vec3D(dx3,dy3,dz3);
    
    if(interpolatedNormal){ // compute nodal normal (weighted by area)
      nodalNormal[n1] += triNorms[iTriangle];
      nodalNormal[n2] += triNorms[iTriangle];
      nodalNormal[n3] += triNorms[iTriangle];
    }

    double nrm = triNorms[iTriangle].norm();
    if(nrm > nrmMax) {
      nrmMax = nrm;
      trMaxNorm = iTriangle;
    }
    // normalize the normal.
    if(nrm > 0.0)
       triNorms[iTriangle] /= nrm;
  }

  if(interpolatedNormal) //normalize nodal normals.
    for(int i=0; i<numStNodes; i++) {
      nodalNormal[i] /= nodalNormal[i].norm();
    }
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistIntersectorPhysBAM::initialize(Domain *d, DistSVec<double,3> &X, IoData &iod) {
  if(this->numFluid<1) {
    fprintf(stderr,"ERROR: numFluid = %d!\n", this->numFluid);
    exit(-1);
  }
  this->X = &X;
  domain = d;
  numLocSub = d->getNumLocSub();
  intersector = new IntersectorPhysBAM*[numLocSub];

  status = new DistVec<int>(domain->getNodeDistInfo());  
  status0 = new DistVec<int>(domain->getNodeDistInfo());  
  occluded_node = new DistVec<bool>(domain->getNodeDistInfo());  
  swept_node = new DistVec<bool>(domain->getNodeDistInfo());  
  pseudoPhi = new DistVec<double>(domain->getNodeDistInfo());  
  boxMin = new DistSVec<double,3>(domain->getNodeDistInfo());
  boxMax = new DistSVec<double,3>(domain->getNodeDistInfo());

  // for hasCloseTriangle
  DistVec<bool> tId(domain->getNodeDistInfo());

  buildSolidNormals();
  d->findNodeBoundingBoxes(X,*boxMin,*boxMax);

  int numIntersectedEdges=0;
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new IntersectorPhysBAM(*(d->getSubDomain()[i]), X(i), (*status)(i), (*status0)(i), (*occluded_node)(i), (*swept_node)(i), *this);
    intersector[i]->hasCloseTriangle(X(i), (*boxMin)(i), (*boxMax)(i), tId(i));
    numIntersectedEdges += intersector[i]->findIntersections(X(i),tId(i),*com);}

  list< pair<Vec3D,int> > points; //pair points with fluid model ID.
#if 1
  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
      pointIt != iod.embed.embedIC.pointMap.dataMap.end();
      pointIt ++){
      int myID = pointIt->second->fluidModelID;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);
      points.push_back(pair<Vec3D,int>(xyz, myID));

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
        exit(-1);}}}
  else
    com->fprintf(stderr, "Point-based initial conditions could not be found.  Assuming single-phase flow\n");
#else
  // TODO(jontg): Hack that only specifically works for IMP45 simulation
  points.push_back(pair<Vec3D,int>(Vec3D(1.000000e-03, 2.500000e+00, 0.000000e+00),IntersectorPhysBAM::INSIDE));
#endif

  findActiveNodesUsingFloodFill(tId,points);
  *status0=*status;

/*
  //Kevin's debug
  char myName[20] = "Xedges_A.top";
  myName[7] += com->cpuNum();
  FILE* myFile = fopen(myName,"w");
  int count = 0;
  myName[8] = '\0'; //get Xedges_A
  for(int i = 0; i < numLocSub; ++i) {
    fprintf(myFile, "Elements %s using FluidNodes\n", myName);
    Vec<bool>& edgeX = intersector[i]->edgeIntersections;
    int (*ptr)[2] = intersector[i]->edges.getPtr();
    int *locToGlob = intersector[i]->locToGlobNodeMap;
    for(int i=0; i<edgeX.size(); i++) {
      if(edgeX[i]) {
        int ni = ptr[i][0], nj = ptr[i][1];
        fprintf(myFile, "%d 1 %d %d\n", ++count, locToGlob[ni]+1, locToGlob[nj]+1);
      }
    }
  }
  fclose(myFile);
*/

}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::findActiveNodesUsingFloodFill(const DistVec<bool>& tId,const list<pair<Vec3D,int> >& points) {
  DistVec<int> nodeColors(X->info());

// Perform local floodFill
  int localColorCount[numLocSub];
#pragma omp parallel for
  for(int iSub=0;iSub<numLocSub;++iSub) {
      localColorCount[iSub]=FloodFill::floodFillSubDomain(*domain->getSubDomain()[iSub],intersector[iSub]->edges,
                                                       intersector[iSub]->edgeIntersections,
                                                       intersector[iSub]->occluded_node,nodeColors(iSub));}
  floodFill->generateConnectionsSet(*domain,*com,nodeColors);

  map<pair<GLOBAL_SUBD_ID,int>,int> localToGlobalColorMap; // only contains valid data for local SubDomains.
  floodFill->unionColors(*domain,*com,numLocSub,localColorCount,localToGlobalColorMap);

// Determine the status of local colors
  map<int,int> globalColorToGlobalStatus;
  for(list<pair<Vec3D,int> >::const_iterator iter = points.begin(); iter!=points.end(); iter++)
    com->fprintf(stderr,"found point (%e %e %e) with FluidModel %d\n", (iter->first)[0], (iter->first)[1], (iter->first)[2], iter->second);

  for(int iSub=0;iSub<numLocSub;++iSub){int ffNode=domain->getSubDomain()[iSub]->findFarfieldNode();
    if(ffNode >= 0){
      fprintf(stderr,"Setting global color %d to Fluid Model %d\n",localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(*domain->getSubDomain()[iSub]),nodeColors(iSub)[ffNode])],IntersectorPhysBAM::INSIDE);
      globalColorToGlobalStatus[localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(*domain->getSubDomain()[iSub]),nodeColors(iSub)[ffNode])]]=IntersectorPhysBAM::INSIDE;}}

  for(int iSub=0;iSub<numLocSub;++iSub){
    SubDomain& sub=*(domain->getSubDomain()[iSub]);
    for(int iElem=0; iElem<sub.numElems(); iElem++)
      for(list<pair<Vec3D,int> >::const_iterator iP=points.begin(); iP!=points.end(); iP++){
        if(sub.isINodeinITet(iP->first, iElem, (*X)(iSub))){ // TODO(jontg): Use a robust implementation of this routine
          fprintf(stderr,"Setting global color %d to Fluid Model %d\n",localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[sub.getElemNodeNum(iElem)[0]])],iP->second);
          globalColorToGlobalStatus[localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[sub.getElemNodeNum(iElem)[0]])]]=iP->second;}}}

  PHYSBAM_MPI_UTILITIES::syncMap(*domain,*com,globalColorToGlobalStatus);

// Compute node status (occluded nodes are OUTSIDE the fluid regime)
#pragma omp parallel for
  for(int iSub=0;iSub<numLocSub;++iSub){
      SubDomain& sub=*(domain->getSubDomain()[iSub]);
      for(int i=0;i<(*status)(iSub).size();++i){int color=localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[i])];
          if((*occluded_node)(iSub)[i] || globalColorToGlobalStatus.find(color)==globalColorToGlobalStatus.end())
              (*status)(iSub)[i]=IntersectorPhysBAM::OUTSIDE;
          else 
              (*status)(iSub)[i]=globalColorToGlobalStatus[color];}}
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::findActiveNodes(const DistVec<bool>& tId) {
#pragma omp parallel for
    for(int iSub=0;iSub<numLocSub;++iSub){
        SubDomain& sub=*(domain->getSubDomain()[iSub]);
        Connectivity &nToN = *(sub.getNodeToNode());
        for(int i=0;i<(*status)(iSub).size();++i){
            if((*occluded_node)(iSub)[i]) (*status)(iSub)[i]=IntersectorPhysBAM::OUTSIDE;
            else if(!(*swept_node)(iSub)[i]) (*status)(iSub)[i]=(*status0)(iSub)[i];
            else{
              int stat=-5;
              for(int n=0;n<nToN.num(i);++n){int neighborNode=nToN[i][n];
                if(i!=neighborNode && !intersector[iSub]->edgeIntersectsStructure(0.0,i,neighborNode) && !(*swept_node)(iSub)[neighborNode]){
                  if(stat != -5 && stat != (*status0)(iSub)[neighborNode]){
                    fprintf(stderr,"ERROR: Swept node (%d) detected inconsistent neighbors.\n", intersector[iSub]->locToGlobNodeMap[i]+1);
                    for(int j=0;j<nToN.num(i);++j){int tmp=nToN[i][j];
                      if(i==tmp) continue; fprintf(stderr,"\tNeighbor %d is (visible [%d], status0 [%d], swept [%d])\n",intersector[iSub]->locToGlobNodeMap[tmp]+1,!intersector[iSub]->edgeIntersectsStructure(0.0,i,tmp),(*status0)(iSub)[tmp],(*swept_node)(iSub)[tmp]);}
                    exit(-1);}
                  stat = (*status0)(iSub)[neighborNode];}}
              (*status)(iSub)[i]=stat;}}}

  operMax<int> maxOp;
  domain->assemble(domain->getLevelPat(),*status,maxOp);

 //Debug
#pragma omp parallel for
    for(int iSub=0;iSub<numLocSub;++iSub) {
        SubDomain& sub=*(domain->getSubDomain()[iSub]);
        Connectivity &nToN = *(sub.getNodeToNode());
        for(int i=0;i<(*status)(iSub).size();++i) 
//          if(intersector[iSub]->locToGlobNodeMap[i]==246960-1)
//            fprintf(stderr,"****** status0 [%d], status [%d], occluded [%d], swept [%d]\n",
//			    (*status0)(iSub)[i],(*status)(iSub)[i],(*occluded_node)(iSub)[i],(*swept_node)(iSub)[i]);
          if((*status)(iSub)[i]==-5) {
            (*status)(iSub)[i]=-2;
            (*occluded_node)(iSub)[i]=true;
          }
    }
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::updateStructure(Vec3D *xs, Vec3D *Vs, int nNodes) {
  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);}

  for (int i=0; i<nNodes; i++) {
    Xs_n[i] = Xs_np1[i];
    Xs_np1[i] = xs[i];
    Xsdot[i] = Vs[i];}
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::updatePhysBAMInterface(Vec3D *particles, int size, const DistSVec<double,3>& fluid_nodes) {
  physInterface->SaveOldState();
  for (int i=0; i<size; i++)
    physInterface->triangulated_surface->particles.X(i+1) = PhysBAM::VECTOR<double,3>(particles[i][0],
                                                                     particles[i][1], particles[i][2]);

  physInterface->Update(true);
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistIntersectorPhysBAM::recompute(double dtf, double dtfLeft, double dts) {
  if (dtfLeft<-1.0e-6) {
    fprintf(stderr,"There is a bug in time-step!\n");
    exit(-1);
  }
  //get current struct coordinates.
  double alpha = 1.0;
  //double alpha = (dts - dtfLeft + dtf)/dts;
  for (int i=0; i<numStNodes; i++) 
    Xs[i] = (1.0-alpha)*Xs_n[i] + alpha*Xs_np1[i];

  // for hasCloseTriangle
  DistVec<bool> tId(X->info());
  
  updatePhysBAMInterface(Xs, numStNodes,*X);

/*
  map<int,int> glob2loc;
  glob2loc[902] = 1;
  glob2loc[906] = 2;
  glob2loc[908] = 3; 
  glob2loc[909] = 4; 
  glob2loc[762] = 5; 
  glob2loc[907] = 6; 
  glob2loc[946] = 7; 
  glob2loc[901] = 8; 
  glob2loc[926] = 9; 
  glob2loc[1035] = 10; 
  glob2loc[1042] = 11; 
  glob2loc[1045] = 12; 
  glob2loc[1044] = 13; 
  glob2loc[1094] = 14; 
  glob2loc[1043] = 15; 
  glob2loc[1238] = 16; 
  glob2loc[1036] = 17; 
  glob2loc[1205] = 18; 


  // Debug
  static int count = 0;
  if(com->cpuNum()==0) {
    FILE *myFile = fopen("particles.top","a");
    fprintf(myFile, "Its %d\n", ++count);
    for(int i=0; i<numStNodes; i++) {
      map<int,int>::iterator it = glob2loc.find(i+1);
      if(it!=glob2loc.end())
        fprintf(myFile,"%d %e %e %e\n", it->second, Xs[i][0], Xs[i][1], Xs[i][2]);
    }
    fclose(myFile);
  }


  if(com->cpuNum()==0 && count==1) {
    FILE *myFile = fopen("topology.top","w");
    fprintf(myFile,"#Id TypeOfElement Node1 Node2 Node3\n");
    for(int i=0; i<numStElems; i++) {
      int id = 0;
      int A = stElem[i][0]+1;
      int B = stElem[i][1]+1;
      int C = stElem[i][2]+1;
      map<int,int>::iterator it1 = glob2loc.find(A);
      map<int,int>::iterator it2 = glob2loc.find(B);
      map<int,int>::iterator it3 = glob2loc.find(C);
      if(it1==glob2loc.end()||it2==glob2loc.end()||it3==glob2loc.end())
        continue;
      fprintf(myFile, "%d %d %d %d %d\n", ++id, 4, it1->second, it2->second, it3->second);
    }
  }
*/


  buildSolidNormals();

  for(int iSub = 0; iSub < numLocSub; ++iSub){
      intersector[iSub]->reset();
      intersector[iSub]->hasCloseTriangle((*X)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), tId(iSub));}

  for(int iSub = 0; iSub < numLocSub; ++iSub) intersector[iSub]->findIntersections((*X)(iSub),tId(iSub),*com);

  for(int iSub = 0; iSub < numLocSub; ++iSub) intersector[iSub]->computeSweptNodes((*X)(iSub),tId(iSub),*com);

  findActiveNodes(tId);
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::printFirstLayer(SubDomain& sub, SVec<double,3>&X, int TYPE)
{
  int mySub = sub.getGlobSubNum();
  int myLocSub = sub.getLocSubNum();
  int (*ptr)[2] = edges.getPtr();
  Connectivity &nToN = *(sub.getNodeToNode());
  std::string fileName = PhysBAM::STRING_UTILITIES::string_sprintf("firstLayer/%d/firstLayer%d.top",mySub);
  std::string nodesName = PhysBAM::STRING_UTILITIES::string_sprintf("Nodes InsideNodes%d\n",mySub);

  FILE* firstLayer = fopen(fileName.c_str(),"w");
  fprintf(firstLayer, "Nodes InsideNodes%s\n", nodesName.c_str());
  for (int i=0; i<sub.numNodes(); i++) fprintf(firstLayer,"%d %e %e %e\n", i+1, X[i][0], X[i][1], X[i][2]);
  fprintf(firstLayer, "Elements FirstLayer%s using InsideNodes%s\n", nodesName.c_str(), nodesName.c_str());
  for (int l=0; l<edges.size(); l++){
    int x1 = ptr[l][0], x2 = ptr[l][1];
    if(edgeIntersectsStructure(0.0,x1,x2)) fprintf(firstLayer,"%d %d %d %d\n", l+1, (int)1, x1+1, x2+1);}
  fclose(firstLayer);
}

//----------------------------------------------------------------------------

IntersectorPhysBAM::IntersectorPhysBAM(SubDomain &sub, SVec<double,3> &X,
                    Vec<int> &stat, Vec<int> &stat0, Vec<bool>& occ_node, Vec<bool>& swe_node, DistIntersectorPhysBAM &distInt) :
                      edgeIntersections(sub.numEdges()),distIntersector(distInt),
                      status(stat),status0(stat0), occluded_node(occ_node), swept_node(swe_node), edges(sub.getEdges()), globIndex(sub.getGlobSubNum())
{
  int numEdges = edges.size();

  status = UNDECIDED;
  status0 = UNDECIDED;

  locToGlobNodeMap = sub.getNodeMap();

  nFirstLayer = 0;
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::reset()
{
  status0 = status;
  status = UNDECIDED;
  occluded_node = false;
  swept_node = false;
  edgeIntersections = false;
  nFirstLayer = 0;

  reverse_mapping.Remove_All();
  forward_mapping.Remove_All();
  xyz.Remove_All();
}

//----------------------------------------------------------------------------

/** Find the closest structural triangle for each node. If no triangle intersect the bounding box of the node,
* no closest triangle exists
*/
int IntersectorPhysBAM::hasCloseTriangle(SVec<double,3> &X,SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<bool> &tId) 
{
  const double TOL=1e-4;
  PhysBAMInterface<double>& physbam_interface=*distIntersector.physInterface;
  int numCloseNodes=0;
  for(int i=0;i<boxMin.size();++i){
    VECTOR<double,3> min_corner(boxMin[i][0],boxMin[i][1],boxMin[i][2]), max_corner(boxMax[i][0],boxMax[i][1],boxMax[i][2]);
    tId[i]=physbam_interface.HasCloseTriangle(VECTOR<double,3>(X[i][0],X[i][1],X[i][2]),min_corner,max_corner,i,TOL,&occluded_node[i]);
    if(tId[i]) ++numCloseNodes;}
  nFirstLayer += numCloseNodes;

  reverse_mapping.Resize(nFirstLayer);
  forward_mapping.Resize(X.size()+1);  PhysBAM::ARRAYS_COMPUTATIONS::Fill(forward_mapping,-1);
  xyz.Preallocate(nFirstLayer);
  for(int i=0;i<X.size();++i) if(tId[i]){
      int shrunk_index=xyz.Append(VECTOR<double,3>(X[i][0],X[i][1],X[i][2]));
      reverse_mapping(shrunk_index) = i;forward_mapping(i+1) = shrunk_index;}

  return numCloseNodes;
}

//----------------------------------------------------------------------------

int IntersectorPhysBAM::findIntersections(SVec<double,3>&X,Vec<bool>& tId,Communicator& com)
{
  int (*ptr)[2] = edges.getPtr();
  const double TOL=1e-4;
  ARRAY<bool> occludedNode(nFirstLayer);
  for(int i=1;i<=nFirstLayer;++i) occludedNode(i) = occluded_node[reverse_mapping(i)];

#if 0 // Debug output
  for(int i=0;i<tId.size();++i){
      if(tId[i] && (i != reverse_mapping(forward_mapping(i+1)))) fprintf(stderr,"BAH %d, %d, %d\n",i+1,forward_mapping(i+1),reverse_mapping(forward_mapping(i+1)));
      else if(!tId[i] && forward_mapping(i+1) != -1) fprintf(stderr,"forward_mapping %d should be -1, is %d (tId = %d)\n",i,forward_mapping(i+1),tId[i]);}
  for(int i=1;i<=nFirstLayer;++i)
      if(i != forward_mapping(reverse_mapping(i)+1)) fprintf(stderr,"BLAH %d, %d, %d\n",i,reverse_mapping(i)+1,forward_mapping(reverse_mapping(i)+1));

  for(int i=0;i<tId.size();++i)
      if(tId[i] && occludedNode(forward_mapping(i+1)) != occluded_node[i]) fprintf(stderr,"OCC ERR, %d -> %d, %d -> %d\n",i,forward_mapping(i+1),occluded_node[i],occludedNode(forward_mapping(i+1)));
#endif

  ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<double>,IntersectionResult<double> > > edgeRes;
  for (int l=0; l<edges.size(); l++) {
    int p = ptr[l][0], q = ptr[l][1];
    if(tId[p] && tId[q])
        edgeRes.Append(TRIPLE<VECTOR<int,3>,IntersectionResult<double>,IntersectionResult<double> >(VECTOR<int,3>(forward_mapping(p+1),forward_mapping(q+1),l),
                                                                                                    IntersectionResult<double>(),IntersectionResult<double>()));}
  distIntersector.getInterface().Intersect(xyz,occludedNode,edgeRes,TOL);

  int intersectedEdgeCount=0;
  edgeIntersections=false;
  for(int i=1; i<=edgeRes.Size(); ++i){
      if(edgeRes(i).y.triangleID != -1 && edgeRes(i).z.triangleID != -1){
          int l=edgeRes(i).x.z;
          edgeIntersections[l] = true;
          CrossingEdgeRes[l] = edgeRes(i).y;
          ReverseCrossingEdgeRes[l] = edgeRes(i).z;
          ++intersectedEdgeCount;}}

#if 1 // Debug output
  bool should_quit=false;
  for(int l=0;l<edges.size();++l){
      int p = ptr[l][0], q = ptr[l][1];
      if((!tId[p] || !tId[q]) && (edgeIntersections[l])){
          fprintf(stderr, "%03d encountered a bad edge, case 1, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n",globIndex,l,edgeIntersections[l],
          p,locToGlobNodeMap[p]+1,tId[p],occluded_node[p],q,locToGlobNodeMap[q]+1,tId[q],occluded_node[q]);should_quit=true;}
      if((occluded_node[p] || occluded_node[q]) && (!edgeIntersections[l])){
          fprintf(stderr, "%03d encountered a bad edge, case 2, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n",globIndex,l,edgeIntersections[l],
          p,locToGlobNodeMap[p]+1,tId[p],occluded_node[p],q,locToGlobNodeMap[q]+1,tId[q],occluded_node[q]);should_quit=true;}
      if(edgeIntersections[l] && (CrossingEdgeRes[l].triangleID == -1 || ReverseCrossingEdgeRes[l].triangleID == -1)){
          fprintf(stderr, "%03d encountered a bad edge, case 3, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n",globIndex,l,edgeIntersections[l],
          p,locToGlobNodeMap[p]+1,tId[p],occluded_node[p],q,locToGlobNodeMap[q]+1,tId[q],occluded_node[q]);should_quit=true;}}
  for(int i=1; i<=edgeRes.Size(); ++i){
      if(!edgeIntersections[edgeRes(i).x.z] && (occluded_node[reverse_mapping(edgeRes(i).x.x)] || occluded_node[reverse_mapping(edgeRes(i).x.y)]))
          fprintf(stderr,"Detected a bad edge, case 2, with the following results: %d (%d,%d, [%d])\n\t%d -> %d (%d), %e, (%e, %e, %e)\n\t%d -> %d (%d), %e, (%e, %e, %e)\n",
                  edgeRes(i).x.z, reverse_mapping(edgeRes(i).x.x),reverse_mapping(edgeRes(i).x.y),edges.find(reverse_mapping(edgeRes(i).x.x),reverse_mapping(edgeRes(i).x.y)),
                  ptr[edgeRes(i).x.z][0],ptr[edgeRes(i).x.z][1],edgeRes(i).y.triangleID, edgeRes(i).y.alpha, edgeRes(i).y.zeta[0], edgeRes(i).y.zeta[1], edgeRes(i).y.zeta[2],
                  ptr[edgeRes(i).x.z][1],ptr[edgeRes(i).x.z][0],edgeRes(i).z.triangleID, edgeRes(i).z.alpha, edgeRes(i).z.zeta[0], edgeRes(i).z.zeta[1], edgeRes(i).z.zeta[2]);
  }

  if(should_quit) exit(-1);
#endif
  return intersectedEdgeCount;
}

//----------------------------------------------------------------------------

int IntersectorPhysBAM::computeSweptNodes(SVec<double,3>& X, Vec<bool>& tId,Communicator& com)
{
  const double TOL=1e-4;
  ARRAY<bool> swept;

  swept.Resize(nFirstLayer); PhysBAM::ARRAYS_COMPUTATIONS::Fill(swept,false);
  distIntersector.getInterface().computeSweptNodes(xyz,swept,(double).1,TOL);

  int numSweptNodes=0;
  for(int i=1;i<=nFirstLayer;++i){swept_node[reverse_mapping(i)] = swept(i);
      if(swept(i)) ++numSweptNodes;}
  return numSweptNodes;
}

//----------------------------------------------------------------------------

LevelSetResult
IntersectorPhysBAM::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
  int edgeNum = edges.find(ni, nj);
  if (!edgeIntersectsStructure(0.0,edgeNum)) {
    fprintf(stderr,"%02d There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
                   globIndex,locToGlobNodeMap[ni]+1, status[ni],occluded_node[ni], locToGlobNodeMap[nj]+1, status[nj],occluded_node[nj],edgeNum);
    PHYSBAM_MPI_UTILITIES::dump_stack_trace();
    exit(-1);}

  const IntersectionResult<double>& result = (ni<nj) ? CrossingEdgeRes[edgeNum] : ReverseCrossingEdgeRes[edgeNum];
  double alpha0 = result.alpha;
  int trueTriangleID = result.triangleID-1;

  LevelSetResult lsRes;
  lsRes.alpha = alpha0;
  lsRes.xi[0] = result.zeta[0];
  lsRes.xi[1] = result.zeta[1];
  lsRes.xi[2] = 1-result.zeta[0]-result.zeta[1];
  lsRes.trNodes[0] = distIntersector.stElem[trueTriangleID][0];
  lsRes.trNodes[1] = distIntersector.stElem[trueTriangleID][1];
  lsRes.trNodes[2] = distIntersector.stElem[trueTriangleID][2];
  lsRes.normVel = lsRes.xi[0]*distIntersector.Xsdot[lsRes.trNodes[0]]
                + lsRes.xi[1]*distIntersector.Xsdot[lsRes.trNodes[1]]
                + lsRes.xi[2]*distIntersector.Xsdot[lsRes.trNodes[2]]; 

  if(!distIntersector.interpolatedNormal)
    lsRes.gradPhi = distIntersector.getSurfaceNorm(trueTriangleID);
  else { //use nodal normals.
    Vec3D ns0 = distIntersector.getNodalNorm(lsRes.trNodes[0]);
    Vec3D ns1 = distIntersector.getNodalNorm(lsRes.trNodes[1]);
    Vec3D ns2 = distIntersector.getNodalNorm(lsRes.trNodes[2]);
    lsRes.gradPhi = lsRes.xi[0]*ns0 + lsRes.xi[1]*ns1 + lsRes.xi[2]*ns2;
    lsRes.gradPhi /= lsRes.gradPhi.norm();
  }

  return lsRes;
}

//----------------------------------------------------------------------------

bool IntersectorPhysBAM::edgeIntersectsStructure(double t, int ni, int nj) const {
  return edgeIntersections[edges.find(ni,nj)];
}

//----------------------------------------------------------------------------

bool IntersectorPhysBAM::edgeIntersectsStructure(double t, int edge_num) const {
  return edgeIntersections[edge_num];
}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::isPointOnSurface(Vec3D pt, int N1, int N2, int N3) 
{
  Vec<Vec3D> &solidX = distIntersector.getStructPosition();
  Vec3D X1 = solidX[N1];
  Vec3D X2 = solidX[N2];
  Vec3D X3 = solidX[N3];

  Vec3D normal = (X2-X1)^(X3-X1);
  normal /=  normal.norm();

  return fabs((pt-X1)*normal);
}
