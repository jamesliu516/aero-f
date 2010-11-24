#include <stdio.h>
#include <iostream>
#include "IntersectorFRG.h"
#include "Vector3D.h"
#include "Communicator.h"
#include <Domain.h>
#include <IoData.h>
#include <Vector.h>
#include <DistVector.h>
#include <Timer.h>
#include "parser/Assigner.h"
#include "Geometry/KDTree.h"
#include <Connectivity.h>
#include <queue>

using std::pair;
using std::map;
using std::list;
using std::set;

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

const int IntersectorFRG::UNDECIDED, IntersectorFRG::INSIDE, IntersectorFRG::OUTSIDE;
int IntersectorFRG::OUTSIDECOLOR;

//----------------------------------------------------------------------------

/** Utility class to find and store bounding boxes for triangles */
class MyTriangle {
  int id;
  double x[3], w[3];
public:
  MyTriangle() {}
  MyTriangle(int i, Vec<Vec3D> &coord, int *nd) {
    id = i;

    for(int j = 0; j <3; ++j) {
      x[j] = std::min(std::min(coord[nd[0]][j], coord[nd[1]][j]), coord[nd[2]][j]);
      w[j] = std::max(std::max(coord[nd[0]][j], coord[nd[1]][j]), coord[nd[2]][j])-x[j];
    }

  }
  double val(int i) const { return x[i]; }
  double width(int i) const { return w[i]; }
  int trId() const { return id; }
};

//----------------------------------------------------------------------------

/** Utility class to find the signed distance to the surface of the structure */
class ClosestTriangle {
  int (*triNodes)[3];
  Vec3D *structX;
  Vec3D *structNorm;
public:
  bool isFirst;
  int bestTrId;
  int n1, n2; //!< if both ns are non negative, the best point is on an edge
  Vec3D n;
  double minDist; //!< Signed distance to the surface
  Vec3D x;
  map<int,int> nd2tri; //needed for the vertex case (only if isConsistent = false)
  map<int,int> nd2tester;

  static const int maxNPairs = 100;

  bool isConsistent, isPositive, hasBeenFixed;
  int nPairs;
  int periTri[maxNPairs];

  /** Check an edge
   * returns true if this edge is the new closest one.
   */
  bool checkEdge(int trId, int p1, int p2, int p3, double trDist);
  void checkVertex(int vn, int trId, double trDist);
  int registerNodes(int ip1, int trId);
  double project(Vec3D x0, int tria, double& xi1, double& xi2);
  double edgeProject(int n1, int n2, double &alpha);
  double getSignedVertexDistance() const;
public:
  ClosestTriangle(int (*triNodes)[3], Vec3D *structX, Vec3D *sN);
  void start(Vec3D x);
  void checkTriangle(int trId);
  double signedDistance() const {
    if(n1 < 0 || n2 >= 0) return minDist;
    else
      return getSignedVertexDistance();
  }

  int bestTriangle() const { return bestTrId; }

  int mode;
};

//----------------------------------------------------------------------------

ClosestTriangle::ClosestTriangle(int (*nd)[3], Vec3D *sX, Vec3D *sN) {
  triNodes = nd;
  structX = sX;
  structNorm = sN;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::start(Vec3D xp) {
  x = xp;
  isFirst = true;
  bestTrId = -1;
  n1 = n2 = -1;
  mode = -1;
  nPairs = 0;
  minDist = 1.0e10;
  nd2tri.clear();
  nd2tester.clear();
}
//----------------------------------------------------------------------------
double ClosestTriangle::edgeProject(int n1, int n2, double &alpha)
{
  Vec3D xA =   structX[n1];
  Vec3D xB =   structX[n2];
  Vec3D AB= xB-xA;
  Vec3D AX = x-xA;
  alpha = AB*AX/(AB*AB);
  Vec3D P = xA + alpha*AB;
  return (P-x).norm();
}
//----------------------------------------------------------------------------

double ClosestTriangle::project(Vec3D x0, int tria, double& xi1, double& xi2)
{
  int iA = triNodes[tria][0];
  int iB = triNodes[tria][1];
  int iC = triNodes[tria][2];
  Vec3D xA = structX[iA];
  Vec3D xB = structX[iB];
  Vec3D xC = structX[iC];

  Vec3D ABC = (xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  double dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (((xB-xp)^(xC-xp))*dir);
  double areaPCA = (((xC-xp)^(xA-xp))*dir);
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

  return dist;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::checkTriangle(int trId) {
  double dist;
  int *nd = triNodes[trId];
  double xi1, xi2;
  dist = project(x, trId, xi1, xi2);
  double xi3 = 1-xi1-xi2;
  const double eps = 0;
  if(xi1 >= -eps && xi2 >= -eps && xi3 >= -eps) {
    if(isFirst || std::abs(minDist) > std::abs(dist)) {
      isFirst = false;
      minDist = dist;
      bestTrId = trId;
      n1 = n2 = -1;
      mode = 0;
    }
  }else {
    if(xi1 < -eps)
      checkEdge(trId, nd[1], nd[2], nd[0], dist);
    if(xi2 < -eps)
      checkEdge(trId, nd[0], nd[2], nd[1], dist);
    if(xi3 < -eps)
      checkEdge(trId, nd[0], nd[1], nd[2], dist);
  }
}

//----------------------------------------------------------------------------

void ClosestTriangle::checkVertex(int ip1, int trId, double trDist) {
  // If this node is already our best solution
  if(n1 == ip1 && n2 < 0) {
    if(nPairs == maxNPairs) {
      std::cerr << "WARNING: On the embedded structure surface, detected too many (>" << maxNPairs << ") peripheral triangles to a node!" << std::endl;
      throw "Too many peripheral triangles";
    }
    periTri[nPairs++] = trId;
    int repeated_node = registerNodes(ip1,trId);

    bool thisSign = (trDist >= 0);

    if((thisSign != isPositive || !isConsistent) && !hasBeenFixed) { // need to figure out the correct sign...
      isConsistent = false;
      if(repeated_node != -1) {// this triangle and another traversed triangle share an edge that has this damn node.
        //Step 1. Determine if it is a "hill" or a "dent". 
        double xi1_temp, xi2_temp;
        double dist2 = project(structX[nd2tester[triNodes[trId][repeated_node]]], trId, xi1_temp, xi2_temp);
        int type = (dist2>0) ? 1 : 2; //1: dent, 2: hill
        //Step 2. Check the angles between (ip1->x) and the normal of the two triangles
        Vec3D ip1_x = x - structX[ip1];
        double alpha = ip1_x*structNorm[trId];
        double beta  = ip1_x*structNorm[nd2tri[triNodes[trId][repeated_node]]];
        if(type==1) 
          isPositive = (alpha>=0&&beta>=0) ? true : false;
        else
          isPositive = (alpha<0&&beta<0) ? false : true;
         
        hasBeenFixed = true; //the sign is determined.
      }
    }
    return;
  }
  // Compute the distance to the node. Determining the sign of the distance
  // is delayed until the very end
  double dist = (x-structX[ip1]).norm();
  if(isFirst || dist < std::abs(minDist)) {
    isFirst = false;
    n1 = ip1;
    n2 = -1;
    nPairs = 1;
    minDist = dist;
    periTri[0] = trId;
    bestTrId = trId;
    nd2tri.clear();
    nd2tester.clear();
    registerNodes(ip1,trId);
    mode = 2;
    hasBeenFixed = false;
    isConsistent = true;
    isPositive = trDist >= 0;
  }
}

//----------------------------------------------------------------------------

int ClosestTriangle::registerNodes(int ip1, int trId)
{
    map<int,int>::iterator it;
    int repeated_node = -1;
    if(ip1==triNodes[trId][0]) {
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][1]] = trId;
          nd2tester[triNodes[trId][1]] = triNodes[trId][2];
        } else
          repeated_node = 1;
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][2]] = trId;
           nd2tester[triNodes[trId][2]] = triNodes[trId][1];
        } else
          repeated_node = 2;
    }
    else if(ip1==triNodes[trId][1]) {
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][2]] = trId;
          nd2tester[triNodes[trId][2]] = triNodes[trId][0];
        } else
          repeated_node = 2;
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][0]] = trId;
           nd2tester[triNodes[trId][0]] = triNodes[trId][2];
        } else
          repeated_node = 0;
    }
    else if(ip1==triNodes[trId][2]) {
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][0]] = trId;
          nd2tester[triNodes[trId][0]] = triNodes[trId][1];
        } else
          repeated_node = 0;
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][1]] = trId;
           nd2tester[triNodes[trId][1]] = triNodes[trId][0];
        } else
          repeated_node = 1;
    } else
      fprintf(stderr,"ERROR (for debug): node %d doesn't belong to triangle %d!\n", ip1+1, trId+1);

    return repeated_node;
} 

//----------------------------------------------------------------------------

bool
ClosestTriangle::checkEdge(int trId, int ip1, int ip2, int p3, double trDist) {
  int p1, p2;
  if(ip1 < ip2) {
    p1 = ip1;
    p2 = ip2;
  } else {
    p1 = ip2;
    p2 = ip1;
  }
  if(n1 == p1 && n2 == p2) { //<! If we've already looked at this edge before
    // If we have already examined this edge and the sign opinions disagree, we need to make the right decision
    // Check if the p3 is in the positive or negative side of the other triangle
    // When signs disagree, the true distance has the opposite sign. because the triangle surface is the
    // end of point with the same sign as p3;
    if(trDist*minDist >= 0)
      return true;
    double xi1, xi2;
    double d2 = project(structX[p3], bestTrId, xi1, xi2);
    if(d2*minDist>0)
      minDist = -minDist;
  } else {
    double dist, alpha;
    double sign = trDist >= 0 ? 1 : -1;
    dist = sign*edgeProject(p1, p2, alpha);

    int cn1 = (alpha > 1) ? p2 : p1;
    int cn2 = p2;
    if(alpha < 0 || alpha > 1) {
      checkVertex(cn1, trId, trDist);
      return true;
    }

    if(isFirst || std::abs(dist) < std::abs(minDist)) {
      isFirst = false;
      bestTrId = trId;
      n = structNorm[bestTrId];
      minDist = dist;
      n1 = cn1;
      n2 = cn2;
      mode = 1;
    }
  }
}

//----------------------------------------------------------------------------

double ClosestTriangle::getSignedVertexDistance() const {
  if(isConsistent)
    return isPositive ? minDist : -minDist;
  else if(hasBeenFixed)
    return isPositive ? minDist : -minDist;
  return minDist;
}

//----------------------------------------------------------------------------

DistIntersectorFRG::DistIntersectorFRG(IoData &iod, Communicator *comm, int nNodes,
                                               double (*xyz)[3], int nElems, int (*abc)[3])
{
  this->numFluid = iod.eqs.numPhase;
  twoPhase = false;

  com = comm;
  com->fprintf(stderr,"- Using Intersector: FRG\n");

  //get embedded structure surface mesh and restart pos
  char *struct_mesh, *struct_restart_pos;
  int sp = strlen(iod.input.prefix) + 1;

  struct_mesh        = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(struct_mesh,"%s%s", iod.input.prefix, iod.input.embeddedSurface);
  struct_restart_pos = new char[sp + strlen(iod.input.positions)];
  if(iod.input.positions[0] != 0)
    sprintf(struct_restart_pos,"%s%s", iod.input.prefix, iod.input.positions);
  else //no restart position file provided
    struct_restart_pos[0] = '\0'; 
  
  //nodal or facet normal?
  interpolatedNormal = (iod.embed.structNormal==EmbeddedFramework::NODE_BASED) ? 
                        true : false;

  //initialize the following to 0(NULL)
  globPhysInterface = 0;
  triNorms = 0;
  nodalNormal = 0;
  status = 0;
  status0 = 0;
  boxMin = 0;
  boxMax = 0;
  distance = 0;
  tId = 0;
  poly = 0;

  //Load files. Compute structure normals. Initialize PhysBAM Interface
  if(nNodes && xyz && nElems && abc)
    init(nNodes, xyz, nElems, abc, struct_restart_pos);
  else
    init(struct_mesh, struct_restart_pos);

  delete[] struct_mesh;
  delete[] struct_restart_pos;
}

//----------------------------------------------------------------------------

DistIntersectorFRG::~DistIntersectorFRG() 
{
  delete [] stElem;
  if(Xs)          delete[] Xs;
  if(Xs0)         delete[] Xs0;
  if(Xs_n)        delete[] Xs_n;
  if(Xs_np1)      delete[] Xs_np1;
  if(Xsdot)       delete[] Xsdot;
  delete solidX;
  delete solidX0;
  delete solidXn;
  if(status)      delete   status;
  if(status0)     delete   status0;
  if(triNorms)    delete[] triNorms;
  if(nodalNormal) delete[] nodalNormal;
  if(boxMax)      delete   boxMax;
  if(boxMin)      delete   boxMin;
  if(distance)    delete   distance;
  if(tId)         delete   tId;
  if(poly)        delete   poly;
  delete globPhysInterface;
  for(int i=0;i<domain->getNumLocSub();++i)
    {
      delete intersector[i];
    }
  delete[] intersector;
}

//----------------------------------------------------------------------------

LevelSetStructure &
DistIntersectorFRG::operator()(int subNum) const {
  return *intersector[subNum];
}

//----------------------------------------------------------------------------

/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorFRG::init(char *solidSurface, char *restartSolidSurface) {

  // Read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(solidSurface, "r");
  if (topFile == NULL) {com->fprintf(stderr, "Embedded structure surface mesh doesn't exist :(\n"); exit(1); }

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
  com->fprintf(stderr,"- IntersectorFRG: Checking the embedded structure surface...   ");
  if (checkTriangulatedSurface()) 
    com->fprintf(stderr,"Ok.\n");
  else {
    com->fprintf(stderr,"\n");
    exit(-1); 
  }

//  getBoundingBox();
  initializePhysBAM();
}

//----------------------------------------------------------------------------
/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorFRG::init(int nNodes, double (*xyz)[3], int nElems, int (*abc)[3], char *restartSolidSurface) {

  // node set
  numStNodes = nNodes;
  // feed data to Xss. 
  Xs      = new Vec3D[numStNodes];
  Xs0     = new Vec3D[numStNodes];
  Xs_n    = new Vec3D[numStNodes];
  Xs_np1  = new Vec3D[numStNodes];
  Xsdot   = new Vec3D[numStNodes];
  solidX  = new Vec<Vec3D>(numStNodes, Xs);
  solidX0 = new Vec<Vec3D>(numStNodes, Xs0);
  solidXn = new Vec<Vec3D>(numStNodes, Xs_n);

  for (int k=0; k<numStNodes; k++) {
    Xs[k]     = Vec3D(xyz[k][0], xyz[k][1], xyz[k][2]);
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);
  }

  // elem set
  numStElems = nElems;
  stElem = new int[numStElems][3];
  for (int i=0; i<numStElems; i++)
    for (int k=0; k<3; k++)
      stElem[i][k] = abc[i][k];

  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");
    if(resTopFile==NULL) {com->fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
    int ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    int ndMax;
    while(1) {
      int nInputs, num1;
      double x1, x2, x3;
      char *endptr, c1[200];

      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;
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
      Xs_n[k]   = Xs[k];
      Xs_np1[k] = Xs[k];
    }
    fclose(resTopFile);
  }

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"- IntersectorFRG: Checking the embedded structure surface...   ");
  if (checkTriangulatedSurface())
    com->fprintf(stderr,"Ok.\n");
  else {
    com->fprintf(stderr,"\n");
    exit(-1);
  }

//  getBoundingBox();
  initializePhysBAM();
}

//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::getBoundingBox() {
  xMin = xMax = Xs[0][0];
  yMin = yMax = Xs[0][1];
  zMin = zMax = Xs[0][2];
  for(int i = 1; i < numStNodes; ++i) {
    xMin = std::min(xMin, Xs[i][0]);
    xMax = std::max(xMax, Xs[i][0]);
    yMin = std::min(yMin, Xs[i][1]);
    yMax = std::max(yMax, Xs[i][1]);
    zMin = std::min(zMin, Xs[i][2]);
    zMax = std::max(zMax, Xs[i][2]);
  }
}
*/
//----------------------------------------------------------------------------

void
DistIntersectorFRG::initializePhysBAM() { //NOTE: In PhysBAM array index starts from 1 instead of 0
  // Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle=new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(numStNodes);
  for (int i=0; i<numStNodes; i++) physbam_solids_particle->X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0],Xs[i][1], Xs[i][2]);
  
  // Initialize the Triangle list.
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > physbam_stElem(numStElems);
  for (int i=0; i<numStElems; i++) physbam_stElem(i+1) = PhysBAM::VECTOR<int,3>(stElem[i][0]+1,stElem[i][1]+1,stElem[i][2]+1);

  // Initialize the mesh.
  // fprintf(stderr,"Initializing the Mesh with %d particles and %d triangles\n",physbam_solids_particle->array_collection.Size(),physbam_stElem.Size());
  PhysBAM::TRIANGLE_MESH *mesh = new PhysBAM::TRIANGLE_MESH(numStNodes,physbam_stElem);
  mesh->Initialize_Adjacent_Elements();mesh->Set_Number_Nodes(numStNodes);

  // Construct TRIANGULATED_SURFACE.
  if(globPhysInterface) delete globPhysInterface;
  globPhysInterface = new PhysBAMInterface<double>(*mesh,*physbam_solids_particle);
}

//----------------------------------------------------------------------------

EdgePair DistIntersectorFRG::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

//----------------------------------------------------------------------------

bool DistIntersectorFRG::checkTriangulatedSurface()
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
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
  return true;
}

//----------------------------------------------------------------------------

void
DistIntersectorFRG::buildSolidNormals() {
  if(!triNorms) triNorms = new Vec3D[numStElems];
  if(interpolatedNormal) {
    if(!nodalNormal)
      nodalNormal = new Vec3D[numStNodes];
    for(int i=0; i<numStNodes; i++)
      nodalNormal[i] = 0.0;
  }

  // Also look to determine a point inside the solid but away from the structure.
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

    // calculate the normal.
    triNorms[iTriangle] = Vec3D(dx2, dy2, dz2)^Vec3D(dx3,dy3,dz3);
    
    if(interpolatedNormal){ // compute nodal normal (weighted by area)
      nodalNormal[n1] += triNorms[iTriangle];
      nodalNormal[n2] += triNorms[iTriangle];
      nodalNormal[n3] += triNorms[iTriangle];
    }

    double nrm = triNorms[iTriangle].norm();
    // normalize the normal.
    if(nrm > 0.0)
      triNorms[iTriangle] /= nrm;
    else
      fprintf(stderr,"ERROR: area of Triangle %d is %e.\n", iTriangle+1, nrm);
  }

  if(interpolatedNormal) //normalize nodal normals.
    for(int i=0; i<numStNodes; i++)
      nodalNormal[i] /= nodalNormal[i].norm();
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::expandScope()
{
  SubDomain **subs = domain->getSubDomain();

  // 1. setup communication pattern
  int numLocSub = domain->getNumLocSub();
  SubDTopo *subTopo = domain->getSubTopo();
  CommPattern<int> trader(subTopo, com, CommPattern<int>::CopyOnSend, CommPattern<int>::NonSym);
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++)
      trader.setLen(sndChannel[iNei], 1+intersector[iSub]->package[iNei].size());     
  }    
  trader.finalize();

  // 2. send packages to neighbour subdomains.
  set<int>::iterator it;
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.getSendBuffer(sndChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      buffer[0] = intersector[iSub]->package[iNei].size(); 
      int count = 0;
      for(it=intersector[iSub]->package[iNei].begin(); it!= intersector[iSub]->package[iNei].end(); it++)
        buffer[++count] = *it;
    }
  }

  // 3. exchange
  trader.exchange();

  // 4. expand scope
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {

/*    char ch[20] = "debug_cpuX";
    ch[9] = 'A' + subs[iSub]->getGlobSubNum();
    FILE* myDebug = fopen(ch,"w");*/

    int *rcvChannel = subs[iSub]->getRcvChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.recData(rcvChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      for(int j=0; j<buffer[0]; j++)
        intersector[iSub]->scope.insert(buffer[j+1]);
    }
/*
    fprintf(myDebug,"Scope (%d) of Sub %d is:\n", intersector[iSub]->scope.size(), subs[iSub]->getGlobSubNum());
    for(set<int>::iterator itt=intersector[iSub]->scope.begin(); itt!=intersector[iSub]->scope.end(); itt++)
      fprintf(myDebug,"%d\n", *itt+1);
    fprintf(myDebug,"\n");
    fclose(myDebug);*/
  }
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistIntersectorFRG::initialize(Domain *d, DistSVec<double,3> &X, IoData &iod, DistVec<int> *point_based_id) {
  if(this->numFluid<1) {
    fprintf(stderr,"ERROR: number of fluid = %d!\n", this->numFluid);
    exit(-1);
  }
  this->X = &X;
  domain = d;
  numLocSub = d->getNumLocSub();
  intersector = new IntersectorFRG*[numLocSub];

  status = new DistVec<int>(domain->getNodeDistInfo());  
  status0 = new DistVec<int>(domain->getNodeDistInfo());  
  boxMin = new DistSVec<double,3>(domain->getNodeDistInfo());
  boxMax = new DistSVec<double,3>(domain->getNodeDistInfo());
  distance = new DistVec<double>(domain->getNodeDistInfo());
  tId = new DistVec<int>(domain->getNodeDistInfo());
  
  poly = new DistVec<bool>(domain->getNodeDistInfo());
  findPoly();

  buildSolidNormals();
  d->findNodeBoundingBoxes(X,*boxMin,*boxMax);

#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new IntersectorFRG(*(d->getSubDomain()[i]), X(i), (*status)(i), (*status0)(i), *this);
    intersector[i]->getClosestTriangles(X(i), (*boxMin)(i), (*boxMax)(i), (*tId)(i), (*distance)(i), false);
    intersector[i]->computeFirstLayerNodeStatus((*tId)(i), (*distance)(i));
  }
  findInAndOut();
  finishStatusByPoints(iod, point_based_id);   
 
#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i)  
    intersector[i]->findIntersections(X(i), false);

//  for(int iSub=0; iSub<numLocSub; iSub++)
//    intersector[iSub]->printFirstLayer(*(domain->getSubDomain()[iSub]), X(iSub), 1); 
}

//----------------------------------------------------------------------------

void 
DistIntersectorFRG::findPoly() {
  if(!poly) {
//    com->fprintf(stderr,"ERROR: poly not initialized.\n"); 
    exit(-1);
  }

  (*poly) = false;
  DistVec<int> tester(domain->getNodeDistInfo());
  tester = 1;
  domain->assemble(domain->getLevelPat(), tester);

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) 
    for(int i=0; i<tester(iSub).size(); i++)
      if(tester(iSub)[i]>2)
        (*poly)(iSub)[i] = true;
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updateStructure(Vec3D *xs, Vec3D *Vs, int nNodes) 
{
//  com->fprintf(stderr,"DistIntersectorFRG::updateStructure called!\n");
  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);
  }

  for (int i=0; i<nNodes; i++) {
    Xs_n[i] = Xs_np1[i];
    Xs_np1[i] = xs[i];
    Xsdot[i] = Vs[i];
  }
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updatePhysBAMInterface() 
{
  for (int i=0; i<numStNodes; i++)
    globPhysInterface->particles.X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0], Xs[i][1], Xs[i][2]);
  globPhysInterface->Update(false); //also rebuild the topology (not really needed for now).
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void DistIntersectorFRG::recompute(double dtf, double dtfLeft, double dts) 
{

  updateStructCoords(0.0, 1.0);
  buildSolidNormals();
  expandScope();

  int subdXing = 0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    intersector[iSub]->reset(); 
    intersector[iSub]->rebuildPhysBAMInterface(Xs, numStNodes, stElem, numStElems);
    intersector[iSub]->getClosestTriangles((*X)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), (*tId)(iSub), (*distance)(iSub), true);
    intersector[iSub]->computeFirstLayerNodeStatus((*tId)(iSub), (*distance)(iSub));
    subdXing += !(intersector[iSub]->finishStatusByHistory(*(domain->getSubDomain()[iSub])));
    if(twoPhase) {
//      intersector[iSub]->printFirstLayer(*(domain->getSubDomain()[iSub]), (*X)(iSub), 1);
      int error = intersector[iSub]->findIntersections((*X)(iSub), true);
      int nCalls = 0;
      while(error && nCalls<100) {
        nCalls++;
        com->fprintf(stderr,"Recomputing intersections (%d) ...\n", error);
        intersector[iSub]->CrossingEdgeRes.clear();
        intersector[iSub]->ReverseCrossingEdgeRes.clear();
        error = intersector[iSub]->findIntersections((*X)(iSub), true);
      }
    }
  }

  if(!twoPhase) {
    com->globalMax(1,&subdXing);
    if(subdXing) { //this happens if the structure is entering/leaving a subdomain. Rarely happens but needs to be delt with...
      com->fprintf(stderr,"FS Interface entering/leaving a subdomain...\n");
      finalizeStatus(); //contains a local communication
    }
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      int error = intersector[iSub]->findIntersections((*X)(iSub), true);
      while(error) {
        com->fprintf(stderr,"Recomputing intersections (%d) ...\n", error);
        intersector[iSub]->CrossingEdgeRes.clear();
        intersector[iSub]->ReverseCrossingEdgeRes.clear();
        error = intersector[iSub]->findIntersections((*X)(iSub), true);
      }
    }
  }

}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updateStructCoords(double c_n, double c_np1)
{
  for (int i=0; i<numStNodes; i++)
    Xs[i] = c_n*Xs_n[i] + c_np1*Xs_np1[i];
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::findInAndOut()
{
  int nUndecided[numLocSub], total;
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());
//  DistVec<int> one(domain->getNodeDistInfo());
//  one = 1;

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) 
    intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);

  while(1) { //get out only when all nodes are decided

    //1. check if all the nodes (globally) are determined
    total = 0;
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalMax(1,&total);
//    com->fprintf(stderr,"total = %d\n",total);
    if(total==0) //done!
      break; 

    //2. try to get a seed from neighbor subdomains, then floodFill
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<status_and_weight(iSub).size(); i++) {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i] + 1; 
          //status(temp) = 0(UNDECIDED),1(INSIDE),or -1(OUTSIDE).
        status_and_weight(iSub)[i][1] = ((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) ? 0 : 1;
      } 
         
    domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_and_weight(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);      
    }
  }
}

//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::findInAndOut()
{
  int nUndecided[numLocSub], total;
  DistVec<int> status_temp(domain->getNodeDistInfo());
  DistVec<int> one(domain->getNodeDistInfo());
  one = 1;

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);

  while(1) { //get out only when all nodes are decided

    //1. check if all the nodes (globally) are determined
    total = 0;
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalMax(1,&total);
//    com->fprintf(stderr,"total = %d\n",total);
    if(total==0) //done!
      break;

    //2. try to get a seed from neighbor subdomains, then floodFill
    status_temp = *status + one; //status_temp = 0(UNDECIDED),1(INSIDE),or -1(OUTSIDE).
    domain->assemble(domain->getLevelPat(),status_temp);
    status_temp -= one;

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_temp(iSub), (*poly)(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }
}
*/
//----------------------------------------------------------------------------

void IntersectorFRG::printFirstLayer(SubDomain& sub, SVec<double,3>&X, int TYPE)
{
  int mySub = sub.getGlobSubNum();
  int myLocSub = sub.getLocSubNum();
  int (*ptr)[2] = edges.getPtr();
  Connectivity &nToN = *(sub.getNodeToNode());
  char fileName[50] = "firstLayera.top";
  char nodesName[2] = "a";

  fileName[10] += mySub;
  nodesName[0] += mySub;

  FILE* firstLayer = fopen(fileName,"w");
  fprintf(firstLayer, "Nodes InsideNodes%s\n", nodesName);
  for (int i=0; i<sub.numNodes(); i++)
    if (status[i]==TYPE) fprintf(firstLayer,"%d %e %e %e\n", i+1, X[i][0], X[i][1], X[i][2]);
  fprintf(firstLayer, "Elements FirstLayer%s using InsideNodes%s\n", nodesName, nodesName);
  for (int l=0; l<edges.size(); l++){
    int x1 = ptr[l][0], x2 = ptr[l][1];
    if (status[x1]!=TYPE || status[x2]!=TYPE) continue;
    int crit = 0;
    for (int i=0; i<nToN.num(x1); i++)
      if (status[nToN[x1][i]]!=TYPE) {crit++; break;}
    for (int i=0; i<nToN.num(x2); i++)
      if (status[nToN[x2][i]]!=TYPE) {crit++; break;}
    if (crit==2)
      fprintf(firstLayer,"%d %d %d %d\n", l+1, (int)1, x1+1, x2+1);
  }
  fclose(firstLayer);

}

//----------------------------------------------------------------------------

void DistIntersectorFRG::finishStatusByPoints(IoData &iod, DistVec<int> *point_based_id)
{
  if((numFluid==1 || numFluid==2) && iod.embed.embedIC.pointMap.dataMap.empty()) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++)
        if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE)
          (*status)(iSub)[i]=1; 

    twoPhase = true; 
    return;
  }

  list< pair<Vec3D,int> > Points; //pair points with fluid model ID.
  map<int,int> pid2id;

  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    int count = 0;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
        pointIt != iod.embed.embedIC.pointMap.dataMap.end();
        pointIt ++){
      int myID = pointIt->second->fluidModelID;
      int myPID = ++count;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);

      if(point_based_id) {
        Points.push_back(pair<Vec3D,int>(xyz, myPID));
        pid2id[myPID] = myID;
      } else
        Points.push_back(pair<Vec3D,int>(xyz, myID));

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
        exit(-1);
      } 
    }
  } else {
    com->fprintf(stderr, "ERROR: (INTERSECTOR) Point-based initial conditions are required for multi-phase FSI.\n");
    exit(-1);
  }

  list< pair<Vec3D,int> >::iterator iter;

  int nUndecided[numLocSub], total;
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());

  // first round
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    nUndecided[iSub] = 0;

    // 1. move "OUTSIDE" nodes to "UNDECIDED".
    for(int i=0; i<(*status)(iSub).size(); i++)
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        (*status)(iSub)[i] = IntersectorFRG::UNDECIDED;
        nUndecided[iSub]++;
      }
    // 2. find seeds by points
    int nSeeds = intersector[iSub]->findSeedsByPoints(*(domain->getSubDomain()[iSub]), (*X)(iSub), Points, nUndecided[iSub]);
    // 3. flood fill if seeds are found
    if(nSeeds>0)
      intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
  }   

  int total0 = 0;
  while(1) { //get out only when all nodes are decided
    //1. check if all the nodes (globally) are determined
    total = 0;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalSum(1,&total);
    if(total==0 || total==total0)
      break; //done
    total0 = total;

    //2. try to get a seed from neighbor subdomains
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<status_and_weight(iSub).size(); i++) {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i] + 1;
          //status(temp) = 0(UNDECIDED),1(INSIDE),or 2,3, ...
        status_and_weight(iSub)[i][1] = ((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) ? 0 : 1;
      }

    domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_and_weight(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }

  if(point_based_id) {
    *point_based_id = -999;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {
        int myPID = (*status)(iSub)[i];
        map<int,int>::iterator pit = pid2id.find(myPID);
        if(pit!=pid2id.end()) {
          (*point_based_id)(iSub)[i] = myPID;
          (*status)(iSub)[i] = pit->second;
        }
      }
  }

  if(total) {//still have undecided nodes. They must be ghost nodes (i.e. covered by solid).
    com->fprintf(stderr,"- IntersectorFRG: Ghost node(s) detected...\n");
    twoPhase = false;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {
        if((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) {
          (*status)(iSub)[i] = IntersectorFRG::OUTSIDECOLOR;
        }
      }
  } else 
    twoPhase = (numFluid<3) ? true : false;
}

//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::finishStatusByPoints(IoData &iod)
{
  if((numFluid==1 || numFluid==2) && iod.embed.embedIC.pointMap.dataMap.empty()) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++)
        if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE)
          (*status)(iSub)[i]=1; 

    twoPhase = true; 
    return;
  }

  list< pair<Vec3D,int> > Points; //pair points with fluid model ID.
  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
        pointIt != iod.embed.embedIC.pointMap.dataMap.end();
        pointIt ++){
      int myID = pointIt->second->fluidModelID;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);
      Points.push_back(pair<Vec3D,int>(xyz, myID));

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
        exit(-1);
      } 
    }
  } else {
    com->fprintf(stderr, "ERROR: (INTERSECTOR) Point-based initial conditions are required for multi-phase FSI.\n");
    exit(-1);
  }

  list< pair<Vec3D,int> >::iterator iter;
  for(iter = Points.begin(); iter!=Points.end(); iter++)
    com->fprintf(stderr,"  - Detected point (%e %e %e) with FluidModel %d\n", (iter->first)[0], (iter->first)[1], (iter->first)[2], iter->second);
  

  int nUndecided[numLocSub], total;
  DistVec<int> status_temp(domain->getNodeDistInfo());
  DistVec<int> one(domain->getNodeDistInfo());
  one = 1;

  // first round
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    nUndecided[iSub] = 0;

    // 1. move "OUTSIDE" nodes to "UNDECIDED".
    for(int i=0; i<(*status)(iSub).size(); i++)
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        (*status)(iSub)[i] = IntersectorFRG::UNDECIDED;
        nUndecided[iSub]++;
      }
    // 2. find seeds by points
    int nSeeds = intersector[iSub]->findSeedsByPoints(*(domain->getSubDomain()[iSub]), (*X)(iSub), Points, nUndecided[iSub]);
    // 3. flood fill if seeds are found
    if(nSeeds>0)
      intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
  }   

  int total0 = 0;
  while(1) { //get out only when all nodes are decided
    //1. check if all the nodes (globally) are determined
    total = 0;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalSum(1,&total);
    com->fprintf(stderr,"Total number of undecided nodes = %d\n", total);    
    if(total==0 || total==total0)
      break; //done

    total0 = total;
    //2. try to get a seed from neighbor subdomains
    status_temp = *status + one; //status_temp = 0,1,2,3,....
    domain->assemble(domain->getLevelPat(),status_temp);
    status_temp -= one;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_temp(iSub), (*poly)(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }

  if(total) {//still have undecided nodes. They must be ghost nodes (i.e. covered by solid).
    twoPhase = false;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {

        if(com->cpuNum()==154&&intersector[iSub]->locToGlobNodeMap[i]==700111-1) {
          fprintf(stderr,"Node 700111 belongs to %d subdomains.\n", intersector[iSub]->nodeToSubD.num(i));
          Connectivity &nToN = *(intersector[iSub]->subD->getNodeToNode());
          for(int j=0; j<nToN.num(i); j++)
            fprintf(stderr,"  -- Neighbour: %d\n", intersector[iSub]->locToGlobNodeMap[nToN[i][j]]+1);
        }

        if(intersector[iSub]->locToGlobNodeMap[i]==700111-1)
          fprintf(stderr,"CPU %d has ndoe 700111 with status %d.\n", com->cpuNum(), (*status)(iSub)[i]);


        if((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) {
          fprintf(stderr,"CPU %d: Node %d is undecided!!!\n", com->cpuNum(), intersector[iSub]->locToGlobNodeMap[i]+1);
          (*status)(iSub)[i] = IntersectorFRG::OUTSIDECOLOR;
        }
      }
  } else 
    twoPhase = (numFluid<3) ? true : false;

  exit(-1);
}
*/
//----------------------------------------------------------------------------

void DistIntersectorFRG::finalizeStatus()
{
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    for(int i=0; i<status_and_weight(iSub).size(); i++) {
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {//this means its real status hasn't been decided.
        status_and_weight(iSub)[i][0] = 0;
        status_and_weight(iSub)[i][1] = 0;
      } else {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i];
        status_and_weight(iSub)[i][1] = 1;
      }
    
//      if(intersector[iSub]->locToGlobNodeMap[i]+1==151870)
//        fprintf(stderr,"Before Assemble: CPU%d: status/weight of 151870 is %d/%d...\n", com->cpuNum(), status_and_weight(iSub)[i][0],status_and_weight(iSub)[i][1]);
    }  


  domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    for(int i=0; i<(*status)(iSub).size(); i++) {
//      if(intersector[iSub]->locToGlobNodeMap[i]+1==151870)
//        fprintf(stderr,"After Assemble: CPU%d: status/weight of 151870 is %d/%d...\n", com->cpuNum(), status_and_weight(iSub)[i][0],status_and_weight(iSub)[i][1]);
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        (*status)(iSub)[i] = status_and_weight(iSub)[i][0] / status_and_weight(iSub)[i][1];

        if((*status)(iSub)[i]<=0 || (*status)(iSub)[i]>numFluid)
          fprintf(stderr,"ERROR: failed at determining status for node %d. Structure may have crossed two or more layers of nodes in one time-step!\n", 
                  (intersector[iSub]->locToGlobNodeMap)[i]+1);
      }
//      if((intersector[iSub]->locToGlobNodeMap)[i]+1==81680)
//        fprintf(stderr,"Node 81680 has status %d.....\n", (*status)(iSub)[i]);
    }
}

//----------------------------------------------------------------------------

bool IntersectorFRG::finishStatusByHistory(SubDomain& sub)
{
  bool good = true;
  bool two_phase = distIntersector.twoPhase;

  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());
  for(int i=0; i<numNodes; i++) {
    if(status[i]==OUTSIDE){
      if(two_phase) {
        status[i] = 1;
        continue;
      }

      if(status0[i]!=INSIDE)
        status[i] = status0[i];
      else {//status0[i]==INSIDE
        bool done = false;    
        for(int iNei=0; iNei<nToN.num(i); iNei++)
          if(status0[nToN[i][iNei]]!=INSIDE) {
            status[i] = status0[nToN[i][iNei]];
            done = true;
            break;
          }
        if(!done) {
          good = false; 
//          fprintf(stderr,"ERROR(CPU %d): F-S Interface moved over 2 (or more) layers of nodes in 1 time-step!\n", distIntersector.com->cpuNum()); 
//          fprintf(stderr,"Node %d: status = %d, status0 = %d\n",locToGlobNodeMap[i]+1,status[i],status0[i]);
//          for(int iNei=0; iNei<nToN.num(i); iNei++) 
//            fprintf(stderr,"Neighbor %d: status = %d, status0 = %d\n", locToGlobNodeMap[nToN[i][iNei]]+1, status[nToN[i][iNei]], status0[nToN[i][iNei]]);
        }
      }
    } 
    else if(status[i]==UNDECIDED)
      if(status0[i]!=UNDECIDED)
        status[i] = status0[i];
      else  
        fprintf(stderr,"ERROR: Unable to determine node status for Node %d\n", locToGlobNodeMap[i]+1);

//    if(locToGlobNodeMap[i]+1==151870)
//      fprintf(stderr,"CPU%d: status of 151870 is %d...\n", distIntersector.com->cpuNum(), status[i]);

  }

  return good;
}

//----------------------------------------------------------------------------

void IntersectorFRG::floodFill(SubDomain& sub, int& nUndecided)
{
  nUndecided = 0;
  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());

  // Propogate status to the entire subdomain.
  // List is used as a FIFO queue
  Vec<int> seed(numNodes); //list of decided nodes  
  Vec<int> level(numNodes);
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i) {
    if(status[i] != UNDECIDED) {
      seed[lead++] = i;
      level[i] = 0;
    } else
      nUndecided++;
  }

  while(next < lead) { //still have seeds not used
    int cur = seed[next++];
    int curStatus = status[cur];
    int curLevel = level[cur];
    for(int i = 0; i < nToN.num(cur); ++i) {
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        level[nToN[cur][i]] = curLevel+1;
        seed[lead++] = nToN[cur][i]; 
        nUndecided--;
      } else 
        if(status[nToN[cur][i]] != curStatus && ( curLevel != 0 || level[nToN[cur][i]] != 0))
          std::cerr << "Incompatible nodes have met: " << locToGlobNodeMap[cur]+1 << "("<< status[cur]
                    << ") and " << locToGlobNodeMap[nToN[cur][i]]+1 << "(" << status[nToN[cur][i]] << ") "
                    << " " << curLevel << " " << level[nToN[cur][i]] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------

IntersectorFRG::IntersectorFRG(SubDomain &sub, SVec<double,3> &X,
                    Vec<int> &stat, Vec<int> &stat0, DistIntersectorFRG &distInt) :
                      subD(&sub), distIntersector(distInt), status(stat), status0(stat0),
                      edges(sub.getEdges()), globIndex(sub.getGlobSubNum()), nodeToSubD(*sub.getNodeToSubD())
{
  physInterface = 0;

  status = UNDECIDED;
  status0 = UNDECIDED;

  OUTSIDECOLOR = distInt.numOfFluids();

  locToGlobNodeMap = sub.getNodeMap();

  iscope = 0;
  package = new set<int>[sub.getNumNeighb()];
  int *neighb = sub.getNeighb();
  for(int i=0; i<sub.getNumNeighb(); i++)
    sub2pack[neighb[i]] = i;

  nFirstLayer = 0;
}

//----------------------------------------------------------------------------

IntersectorFRG::~IntersectorFRG()
{
  delete[] package;
  if(iscope) delete[] iscope;

  if(physInterface) delete physInterface;
}

//----------------------------------------------------------------------------

void IntersectorFRG::reset()
{
  for(int i=0; i<subD->getNumNeighb(); i++)
    package[i].clear();

  if(iscope) {delete[] iscope; iscope = 0;}
  n2p.clear();
  particle.clear();

  CrossingEdgeRes.clear();
  ReverseCrossingEdgeRes.clear();

  status0 = status;
  status = UNDECIDED;
  if(physInterface) {delete physInterface; physInterface = 0;}

  nFirstLayer = 0;
}

//----------------------------------------------------------------------------
//TODO: discuss with Jon
void IntersectorFRG::rebuildPhysBAMInterface(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElems)
{ //IMPORTANT: In PhysBAM array index starts fr*om 1 instead of 0
  int nPar, nTri = scope.size();
  if(nTri==0)
    return;
  nPar = buildScopeTopology(sElem, nsElems);

  int count = 0;
  // Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle=new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(nPar);
  for (list<int>::iterator lit=particle.begin(); lit!=particle.end(); lit++)
    physbam_solids_particle->X(++count) = PhysBAM::VECTOR<double,3>(Xs[*lit][0],Xs[*lit][1], Xs[*lit][2]);

  // Initialize the Triangle list
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > physbam_stElem;
  physbam_stElem.Preallocate(nTri);
  for (set<int>::iterator it=scope.begin(); it!=scope.end(); it++)
    physbam_stElem.Append(PhysBAM::VECTOR<int,3>(n2p[sElem[*it][0]]+1, n2p[sElem[*it][1]]+1, n2p[sElem[*it][2]]+1));

  // Initialize the mesh.
  // fprintf(stderr,"Initializing the Mesh with %d particles and %d triangles\n",physbam_solids_particle->array_collection.Size(),physbam_stElem.Size());
  PhysBAM::TRIANGLE_MESH *mesh = new PhysBAM::TRIANGLE_MESH(count,physbam_stElem);
  mesh->Initialize_Adjacent_Elements();mesh->Set_Number_Nodes(count);

  // Construct TRIANGULATED_SURFACE.
  if(physInterface) delete physInterface;
  physInterface = new PhysBAMInterface<double>(*mesh,*physbam_solids_particle);





/*







  // Initialize the Particles list 
  int count = 0;
  physbam_solids_particle = new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(nPar);
  for (list<int>::iterator lit=particle.begin(); lit!=particle.end(); lit++) 
    physbam_solids_particle->X(++count) = PhysBAM::VECTOR<double,3>(Xs[*lit][0],Xs[*lit][1], Xs[*lit][2]);

  // Initialize the Triangle list
  physbam_stElem = new PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> >();
  for (set<int>::iterator it=scope.begin(); it!=scope.end(); it++)
    physbam_stElem->Append(PhysBAM::VECTOR<int,3>(n2p[sElem[*it][0]]+1, n2p[sElem[*it][1]]+1, n2p[sElem[*it][2]]+1));

  // Construct TRIANGLE_MESH triangle_mesh (it stores its own copy of physbam_stElem).
  physbam_triangle_mesh = new PhysBAM::TRIANGLE_MESH(physbam_solids_particle->array_collection.Size(), *physbam_stElem);
  physbam_triangle_mesh->Initialize_Adjacent_Elements(); //need this?

  // Construct TRIANGULATED_SURFACE. (only stores the reference to physbam_triangle_mesh and physbam_solids_particle)
  physbam_triangulated_surface = new PhysBAM::TRIANGULATED_SURFACE<double>(*physbam_triangle_mesh, *physbam_solids_particle);
  physbam_triangulated_surface->Update_Triangle_List();

  // Construct PhysBAMInterface. (only stores a reference to physbam_triangulated_surface)
  physInterface = new PhysBAMInterface<double>(*physbam_triangulated_surface);
*/
/*  char ch[20] = "physbam_cpuX";
  ch[11] = '0' + subD->getGlobSubNum();
  FILE* myDebug = fopen(ch,"w");
  fprintf(myDebug,"Particles: \n");
  for(int i=0; i<nPar; i++)
    fprintf(myDebug, "%d %e %e %e\n", i+1, physbam_solids_particle.X(i+1)[1], physbam_solids_particle.X(i+1)[2], physbam_solids_particle.X(i+1)[3]);
  fprintf(myDebug,"Node to Particle map:\n");
  for(map<int,int>::iterator it=n2p.begin(); it!=n2p.end(); it++)
    fprintf(myDebug,"%d -> %d\n", it->first, it->second);
  fprintf(myDebug,"Elements:\n");
  for(int i=0; i<nTri; i++)
    fprintf(myDebug,"%d: %d %d %d\n", i+1, physbam_stElem(i+1)[1], physbam_stElem(i+1)[2], physbam_stElem(i+1)[3]);
  fclose(myDebug);
*/
}

//----------------------------------------------------------------------------
/*
void IntersectorFRG::buildSolidNormals(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElem)
{
  int nTr = scope.size();
  if(nTr==0)
    return;
  int nPt = particle.size();
  bool interp = distIntersector.interpolatedNormal;

  // allocate memory.
  if(trNormal) delete[] trNormal;
  trNormal = new Vec3D[nTr];
  if(interp) {
    if(ndNormal) delete[] ndNormal;
    ndNormal = new Vec3D[nPt];
    for(int i=0; i<nPt; i++)
      ndNormal[i] = 0.0;
  }

  for(int iTr=0; iTr<nTr; iTr++) {
    int n1 = sElem[iscope[iTr]][0];
    int n2 = sElem[iscope[iTr]][1];
    int n3 = sElem[iscope[iTr]][2];
    double x1 = Xs[n1][0];
    double y1 = Xs[n1][1];
    double z1 = Xs[n1][2];
    double dx2 = Xs[n2][0]-x1;
    double dy2 = Xs[n2][1]-y1;
    double dz2 = Xs[n2][2]-z1;
    double dx3 = Xs[n3][0]-x1;
    double dy3 = Xs[n3][1]-y1;
    double dz3 = Xs[n3][2]-z1;

    // calculate the normal.
    trNormal[iTr] = Vec3D(dx2,dy2,dz2)^Vec3D(dx3,dy3,dz3);

    if(interp){ // compute nodal normal (weighted by 2*area)
      ndNormal[n2p[n1]] += trNormal[iTr];
      ndNormal[n2p[n2]] += trNormal[iTr];
      ndNormal[n2p[n3]] += trNormal[iTr];
    } 

    // normalize the normal.
    double nrm = trNormal[iTr].norm();
    if(nrm > 0.0)
       trNormal[iTr] /= nrm;
    else
      fprintf(stderr,"ERROR: Area of Triangle %d is %e.\n", iscope[iTr]+1, 0.5*nrm); 
  }

  if(interp) //normalize nodal normals.
    for(int i=0; i<nPt; i++) 
      ndNormal[i] /= ndNormal[i].norm(); //TODO: need a local communication
}
*/
//----------------------------------------------------------------------------

int IntersectorFRG::buildScopeTopology(int (*sElem)[3], int nsElem)
{ // construct iscope, n2p, and particles.
  iscope = new int[scope.size()]; //it's memory has been released in "reset".
  int nd, newID = 0, newTrID = 0;
  for(set<int>::iterator it=scope.begin(); it!=scope.end(); it++){
    iscope[newTrID++] = *it;
    for(int k=0; k<3; k++) {
      nd = sElem[*it][k];
      if(n2p.find(nd)==n2p.end()) {
        n2p[nd] = newID++;
        particle.push_back(nd);
      }
    }
  }
  return newID;
}

//----------------------------------------------------------------------------

/** Find the closest structural triangle for each node. If no triangle intersect the bounding box of the node,
* no closest triangle exists
*/
void IntersectorFRG::getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist, bool useScope) 
{
  int (*triNodes)[3] = distIntersector.stElem;
  Vec3D *structX = (distIntersector.getStructPosition()).data();
  int ntri;
  MyTriangle *myTris;

  // build the KDTree
  if(useScope) {
    ntri = scope.size();
    if(ntri==0) {
      tId = -1;
      return;
    }
    myTris = new MyTriangle[ntri];
    int count = 0;
    set<int>::iterator it;
    for(it=scope.begin(); it!=scope.end(); it++)
      myTris[count++] = MyTriangle(*it, distIntersector.getStructPosition(), triNodes[*it]);
  } else {
    ntri = distIntersector.getNumStructElems();
    myTris = new MyTriangle[ntri];
    for(int i = 0; i < ntri; ++i)
      myTris[i] = MyTriangle(i, distIntersector.getStructPosition(), triNodes[i]);
  }
  KDTree<MyTriangle> structureTree(ntri, myTris);

  // clear scope for refill 
  scope.clear();

  // find candidates
  ClosestTriangle closestTriangle(triNodes, structX, distIntersector.triNorms);
  int nMaxCand = 500;
  MyTriangle *candidates = new MyTriangle[nMaxCand];
  for(int i = 0; i < X.size(); ++i) {
    int nFound = structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates, nMaxCand);
    if(nFound > nMaxCand) {
//      std::cerr << "For Fluid node " << locToGlobNodeMap[i]+1 << ", number of candidates: " << nFound << std::endl;
      nMaxCand = nFound;
      delete [] candidates;
      candidates = new MyTriangle[nMaxCand];
      structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates, nMaxCand);
    }
    closestTriangle.start(X[i]);
    for(int j = 0; j < nFound; ++j) {
      int myId = candidates[j].trId();
      scope.insert(myId);
      addToPackage(i, myId);    
      closestTriangle.checkTriangle(myId);
    }
    
    if(nFound <= 0) {
      tId[i] = -1;
      dist[i] = 0.0;
    }
    else {
      tId[i] = closestTriangle.bestTriangle();
      if(tId[i] < 0)
        std::cout << "Horror!!!" << std::endl;
      dist[i] = closestTriangle.signedDistance();
    }
  }

  delete [] candidates;
  delete [] myTris;
}

//----------------------------------------------------------------------------

void IntersectorFRG::computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist)
{
  for (int i=0; i<tId.size(); i++) {
    if (tId[i]<0)
      continue;
    status[i] = (dist[i]>=0.0) ? INSIDE : OUTSIDE;
    nFirstLayer++;
  }
}

//----------------------------------------------------------------------------

int IntersectorFRG::findNewSeedsAfterMerging(SVec<int,2>& status_and_weight, int& nUndecided)
{
  int numNewSeeds = 0;
  for(int i=0; i<status.size(); i++){
    if(status[i]!=UNDECIDED || status_and_weight[i][1]==0)
      continue; //already decided || didn't get anything from neighbours
    status[i] = status_and_weight[i][0] / status_and_weight[i][1] - 1;//back to normal conventional
    nUndecided--;
    numNewSeeds++;
  }
  return numNewSeeds;
}

//----------------------------------------------------------------------------

int IntersectorFRG::findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided)
{
  int numNodes = status.size();
  int myStatus, numNewSeeds = 0;

  for(int i=0; i<numNodes; i++){
    if(status_temp[i]==status[i] || poly[i])
      continue; // inside or on the boundary but the neighbor hasn't decided it, or lies on n>2 subdomains.

    myStatus = status_temp[i];
//    if(myStatus!=INSIDE && myStatus!=OUTSIDE){ 
//      fprintf(stderr,"horror!\n"); exit(-1);}
    if (status[i]==UNDECIDED) {
      status[i] = myStatus;
      nUndecided--;
      numNewSeeds++;
    }
//    else if (status[i]!=myStatus)
//      fprintf(stderr,"ERROR: Node %d got different statuses from different subdomains.\n", locToGlobNodeMap[i]+1);
  }
  return numNewSeeds;
}

//----------------------------------------------------------------------------

int IntersectorFRG::findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list<pair<Vec3D,int> > P, int& nUndecided)
{
  int *myNodes, nSeeds = 0;
  list<pair<Vec3D,int> >::iterator iP;

  for(int iElem=0; iElem<sub.numElems(); iElem++) {
    myNodes = sub.getElemNodeNum(iElem);
    if(status[myNodes[0]]==INSIDE && status[myNodes[1]]==INSIDE &&
       status[myNodes[2]]==INSIDE && status[myNodes[3]]==INSIDE)//this tet is inside farfield fluid
      continue;

    for(iP=P.begin(); iP!=P.end(); iP++)
      if(sub.isINodeinITet(iP->first, iElem, X)) 
        for(int i=0; i<4; i++)
          if(status[myNodes[i]]==UNDECIDED) {
            status[myNodes[i]] = iP->second;
            nUndecided--;
            nSeeds++;
          }
  }

  return nSeeds;
}

//----------------------------------------------------------------------------

void IntersectorFRG::noCheckFloodFill(SubDomain& sub, int& nUndecided)
{
  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());

  // Propogate status to the entire subdomain.
  // List is used as a FIFO queue
  Vec<int> seed(numNodes); //list of decided nodes
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i)
    if(status[i] != UNDECIDED && status[i] != INSIDE) {
      seed[lead++] = i;
    }

  while(next < lead) { //still have seeds not used
    int cur = seed[next++];
    int curStatus = status[cur];
    for(int i = 0; i < nToN.num(cur); ++i)
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        seed[lead++] = nToN[cur][i];
        nUndecided--;
      } 
  }
}

//----------------------------------------------------------------------------

int IntersectorFRG::findIntersections(SVec<double,3>&X, bool useScope)
{
  int error = 0;

  int (*ptr)[2] = edges.getPtr();
  const double TOL = 1.0e-4;
  int MAX_ITER = 20;
  int max_iter = 0;

  bool GlobPhysBAMUpdated = false;

  for (int l=0; l<edges.size(); l++) {
    int p = ptr[l][0], q = ptr[l][1];
    if(status[p]==status[q]) continue;

    //now need to find an intersection for this edge .
    Vec3D xp(X[p]), xq(X[q]);
    Vec3D dir = xq - xp;
    Vec3D xpPrime, xqPrime;

    IntersectionResult<double> res1;

    for (int iter=0; iter<MAX_ITER; iter++) {
      double coeff = iter*iter*TOL;
      Vec3D xpPrime = xp - coeff*dir;
      Vec3D xqPrime = xq + coeff*dir;
      VECTOR<double,3> xyz1(xpPrime[0],xpPrime[1],xpPrime[2]), 
                       xyz2(xqPrime[0],xqPrime[1],xqPrime[2]);

      if(useScope) {
        res1 = physInterface->Intersect(xyz1, xyz2, coeff*dir.norm());
        if(res1.triangleID>0)
          res1.triangleID = iscope[res1.triangleID-1] + 1;
      } else 
        res1 = distIntersector.getInterface().Intersect(xyz1,xyz2, coeff*dir.norm());
      // the triangle Id stored in edgeRes starts from 1, i.e. using PhysBAM convention.

      if (res1.triangleID>0 /*|| edgeRes(2).y.triangleID>0*/) {
        CrossingEdgeRes[l] = res1;
        ReverseCrossingEdgeRes[l] = res1; //TODO: redundent!
        if (iter>max_iter) max_iter = iter;
        break;
      }
    }

    if (res1.triangleID<0 && useScope) { // try global intersector as a 'fail-safe'.
      if(!GlobPhysBAMUpdated) {
        distIntersector.updatePhysBAMInterface();
        GlobPhysBAMUpdated = true;
//        fprintf(stderr,"Now using the global intersector...\n");
      }

      for (int iter=0; iter<MAX_ITER; iter++) {
        double coeff = iter*iter*TOL;
        Vec3D xpPrime = xp - coeff*dir;
        Vec3D xqPrime = xq + coeff*dir;
        VECTOR<double,3> xyz1(xpPrime[0],xpPrime[1],xpPrime[2]), 
                         xyz2(xqPrime[0],xqPrime[1],xqPrime[2]);
        res1 = distIntersector.getInterface().Intersect(xyz1,xyz2, coeff*dir.norm()); 

        if (res1.triangleID>0 /*|| edgeRes(2).y.triangleID>0*/) {
          CrossingEdgeRes[l] = res1;
          ReverseCrossingEdgeRes[l] = res1; //TODO: redundent!
          if (iter>max_iter) max_iter = iter;
          break;
        }
      }
 
      if(res1.triangleID<0) {
         error++;
//         fprintf(stderr,"WARNING: No intersection between node %d(status = %d, status0 = %d) and %d(status = %d, status0 = %d). \n",
//                        locToGlobNodeMap[p]+1,status[p], status0[p], locToGlobNodeMap[q]+1,status[q], status0[q]);
         if(status[p]!=status0[p] || status[q]!=status0[q]) {
           status[p] = status0[p];
           status[q] = status0[q];
         } else
           status[p] = status[q] = 0; //TODO: this is not a fix!
      }
    }
  }

  return error;
}

//----------------------------------------------------------------------------

void IntersectorFRG::addToPackage(int iNode, int trID)
{
  int nSub = nodeToSubD.num(iNode);
  for(int iSub=0; iSub<nSub; iSub++) {
    if(nodeToSubD[iNode][iSub]==globIndex) 
      continue;
    package[sub2pack[nodeToSubD[iNode][iSub]]].insert(trID);  
  }
}

//----------------------------------------------------------------------------

LevelSetResult
IntersectorFRG::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
  int edgeNum = edges.find(ni, nj);
  if (!edgeIntersectsStructure(0.0,ni,nj)) {
    fprintf(stderr,"There is no intersection between node %d(status:%d) and %d(status:%d)! Abort...\n",
                   locToGlobNodeMap[ni]+1, status[ni], locToGlobNodeMap[nj]+1, status[nj]);
    exit(-1);
  }

  IntersectionResult<double> result; //need to determine which result to choose.
  IntersectionResult<double> rij = CrossingEdgeRes[edgeNum];
  IntersectionResult<double> rji = ReverseCrossingEdgeRes[edgeNum];
  double alpha0 = 0.0;

  if (rij.triangleID>0 && rji.triangleID>0 && rij.triangleID==rji.triangleID) {
    result = rij;
    alpha0 = (ni<nj)? result.alpha : 1.0-result.alpha;
  }  
  else if (rij.triangleID>0 && rji.triangleID>0 && rij.triangleID!=rji.triangleID) {
    result = (ni<nj) ? rij : rji;
    alpha0 = result.alpha;
  }
  else if (rij.triangleID>0 && rji.triangleID<0) {
    result = rij;
    alpha0 = (ni<nj)? result.alpha : 1.0-result.alpha; 
  }
  else if (rij.triangleID<0 && rji.triangleID>0) {
    result = rji;
    alpha0 = (ni<nj)? 1.0-result.alpha : result.alpha;
  }
  else //we really have no intersection for this edge!
    fprintf(stderr,"ERROR: intersection between node %d and node %d can not be detected.\n", locToGlobNodeMap[ni]+1, locToGlobNodeMap[nj]+1);

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

bool IntersectorFRG::edgeIntersectsStructure(double t, int eij) const       
{
  return status[edges.getPtr()[eij][0]] != status[edges.getPtr()[eij][1]];
}

//----------------------------------------------------------------------------

void IntersectorFRG::projection(Vec3D x0, int tria, double& xi1, double& xi2, double& dist)
{
  int iA = distIntersector.stElem[tria][0];
  int iB = distIntersector.stElem[tria][1];
  int iC = distIntersector.stElem[tria][2];
  Vec3D xA = distIntersector.Xs[iA];
  Vec3D xB = distIntersector.Xs[iB];
  Vec3D xC = distIntersector.Xs[iC];

  Vec3D ABC = 0.5*(xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (0.5*(xB-xp)^(xC-xp))*dir;
  double areaPCA = (0.5*(xC-xp)^(xA-xp))*dir;
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

}

//----------------------------------------------------------------------------

double IntersectorFRG::isPointOnSurface(Vec3D pt, int N1, int N2, int N3) 
{
  Vec<Vec3D> &solidX = distIntersector.getStructPosition();
  Vec3D X1 = solidX[N1];
  Vec3D X2 = solidX[N2];
  Vec3D X3 = solidX[N3];

  Vec3D normal = (X2-X1)^(X3-X1);
  normal /=  normal.norm();

  return fabs((pt-X1)*normal);
}


