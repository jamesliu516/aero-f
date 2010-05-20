#include "../LevelSet/LevelSetStructure.C"
#include <stdio.h>
#include <iostream>
#include "PhysBAMIntersect.h"
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

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

const int PhysBAMIntersector::UNDECIDED, PhysBAMIntersector::INSIDE, PhysBAMIntersector::OUTSIDE;
//CAUTION: INSIDE/OUTSIDE means inside/outside fluid #0.

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

  static const int maxNPairs = 100;

  bool isConsistant, isPositive;
  int nPairs;
  int periTri[maxNPairs];

  /** Check an edge
   * returns true if this edge is the new closest one.
   */
  bool checkEdge(int trId, int p1, int p2, int p3, double trDist);
  void checkVertex(int vn, int trId, double trDist);
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
      std::cerr << "Too many peripheral triangles to a node!" << std::endl;
      throw "Too many peripheral triangles";
    }
    periTri[nPairs++] = trId;
    bool thisSign = (trDist >= 0);
    if(thisSign != isPositive)
      isConsistant = false;
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
    mode = 2;
    isConsistant = true;
    isPositive = trDist >= 0;
  }
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
  if(isConsistant)
    return isPositive ? minDist : -minDist;
  return minDist;
}

//----------------------------------------------------------------------------

DistPhysBAMIntersector::DistPhysBAMIntersector(IoData &iod, Communicator *comm) 
{
  this->numFluid = iod.eqs.numPhase;
  com = comm;

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
  
  //nodal or facet normal?
  interpolatedNormal = (iod.embed.structNormal==EmbeddedFramework::NODE_BASED) ? 
                        true : false;

  //initialize the following to 0(NULL)
  physInterface = 0;
  triNorms = 0;
  triSize = 0;
  nodalNormal = 0;
  status = 0;
  status0 = 0;
  boxMin = 0;
  boxMax = 0;
  poly = 0;

  //Load files. Compute structure normals. Initialize PhysBAM Interface
  init(struct_mesh, struct_restart_pos);

  delete[] struct_mesh;
  delete[] struct_restart_pos;
}

//----------------------------------------------------------------------------

DistPhysBAMIntersector::~DistPhysBAMIntersector() 
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
  if(poly)  delete   poly;
}

//----------------------------------------------------------------------------

LevelSetStructure &
DistPhysBAMIntersector::operator()(int subNum) const {
  return *intersector[subNum];
}

//----------------------------------------------------------------------------

/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistPhysBAMIntersector::init(char *solidSurface, char *restartSolidSurface) {

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
  com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) 
    com->fprintf(stderr,"Ok.\n");
  else 
    exit(-1); 

  getBoundingBox();
  initializePhysBAM();
}

//----------------------------------------------------------------------------

void DistPhysBAMIntersector::getBoundingBox() {
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

//----------------------------------------------------------------------------

void
DistPhysBAMIntersector::initializePhysBAM() { //NOTE: In PhysBAM array index starts from 1 instead of 0
// Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >& physbam_solids_particle = *new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle.array_collection.Resize(numStNodes);
  for (int i=0; i<numStNodes; i++) 
    physbam_solids_particle.X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0],Xs[i][1], Xs[i][2]);
  
  // Initialize the Triangle list
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > & physbam_stElem=*new PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> >();
  for (int i=0; i<numStElems; i++){
    int nx, ny, nz;
    nx = stElem[i][0] + 1;  ny = stElem[i][1] + 1;  nz = stElem[i][2] + 1;
    physbam_stElem.Append(PhysBAM::VECTOR<int,3>(nx, ny, nz));
  }

  // Construct TRIANGLE_MESH triangle_mesh.
  PhysBAM::TRIANGLE_MESH& physbam_triangle_mesh=*new PhysBAM::TRIANGLE_MESH(physbam_solids_particle.array_collection.Size(), physbam_stElem);
  physbam_triangle_mesh.Initialize_Adjacent_Elements();

  // Construct TRIANGULATED_SURFACE.
  PhysBAM::TRIANGULATED_SURFACE<double>& physbam_triangulated_surface=*new PhysBAM::TRIANGULATED_SURFACE<double>(physbam_triangle_mesh, physbam_solids_particle);
  physbam_triangulated_surface.Update_Triangle_List();
  if(physInterface) delete physInterface;
  physInterface = new PhysBAMInterface<double>(physbam_triangulated_surface);
}

//----------------------------------------------------------------------------

EdgePair DistPhysBAMIntersector::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

//----------------------------------------------------------------------------

bool DistPhysBAMIntersector::checkTriangulatedSurface()
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
DistPhysBAMIntersector::buildSolidNormals() {
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
DistPhysBAMIntersector::initialize(Domain *d, DistSVec<double,3> &X, IoData &iod) {
  if(this->numFluid<1) {
    fprintf(stderr,"ERROR: numFluid = %d!\n", this->numFluid);
    exit(-1);
  }
  this->X = &X;
  domain = d;
  numLocSub = d->getNumLocSub();
  intersector = new PhysBAMIntersector*[numLocSub];
  pseudoPhi = new DistVec<double>(X.info()); //TODO: not needed at all!

  status = new DistVec<int>(domain->getNodeDistInfo());  
  status0 = new DistVec<int>(domain->getNodeDistInfo());  
  boxMin = new DistSVec<double,3>(domain->getNodeDistInfo());
  boxMax = new DistSVec<double,3>(domain->getNodeDistInfo());

  poly = new DistVec<bool>(domain->getNodeDistInfo());
  findPoly();

  // for getClosestTriangles
  DistVec<double> distance(X.info());
  DistVec<int> tId(X.info());

  buildSolidNormals();
  d->findNodeBoundingBoxes(X,*boxMin,*boxMax);

  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new PhysBAMIntersector(*(d->getSubDomain()[i]), X(i), (*status)(i), (*status0)(i), *this);
    intersector[i]->getClosestTriangles(X(i), (*boxMin)(i), (*boxMax)(i), tId(i), distance(i));
    intersector[i]->computeFirstLayerNodeStatus(tId(i), distance(i));
  }
  findInAndOut();
  finishStatusByPoints(iod);   
 
  for(int i = 0; i < numLocSub; ++i) 
    intersector[i]->findIntersections(X(i));

//  for(int iSub=0; iSub<numLocSub; iSub++)
//    intersector[iSub]->printFirstLayer(*(domain->getSubDomain()[iSub]), X(iSub), 1); 
}

//----------------------------------------------------------------------------

void 
DistPhysBAMIntersector::findPoly() {
  if(!poly) {
    com->fprintf(stderr,"ERROR: poly not initialized.\n"); 
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

void
DistPhysBAMIntersector::updateStructure(Vec3D *xs, Vec3D *Vs, int nNodes) {

//  com->fprintf(stderr,"DistPhysBAMIntersector::updateStructure called!\n");
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

void
DistPhysBAMIntersector::updatePhysBAMInterface(Vec3D *particles, int size) {
  for (int i=0; i<size; i++)
    physInterface->triangulated_surface->particles.X(i+1) = PhysBAM::VECTOR<double,3>(particles[i][0],
                                                                     particles[i][1], particles[i][2]);
  physInterface->Update(true); //also rebuild the topology (not really needed for now).
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistPhysBAMIntersector::recompute(double dtf, double dtfLeft, double dts) {

  if (dtfLeft<-1.0e-8) {
    fprintf(stderr,"There is a bug in time-step!\n");
    exit(-1);
  }
  //get current struct coordinates.
  double alpha = 1.0;
  //double alpha = (dts - dtfLeft + dtf)/dts;
  for (int i=0; i<numStNodes; i++) 
    Xs[i] = (1.0-alpha)*Xs_n[i] + alpha*Xs_np1[i];

  // for getClosestTriangles
  DistVec<double> distance(X->info());
  DistVec<int> tId(X->info());
  
  updatePhysBAMInterface(Xs, numStNodes);
  getBoundingBox();
  buildSolidNormals();

  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    intersector[iSub]->reset();
    intersector[iSub]->getClosestTriangles((*X)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), tId(iSub), distance(iSub));
    intersector[iSub]->computeFirstLayerNodeStatus(tId(iSub), distance(iSub));
  }

  findInAndOut();
 
  for(int iSub = 0; iSub < numLocSub; ++iSub){ 
    intersector[iSub]->finishStatusByHistory(*(domain->getSubDomain()[iSub]));   
    intersector[iSub]->findIntersections((*X)(iSub));
  }
}

//----------------------------------------------------------------------------

void DistPhysBAMIntersector::findInAndOut()
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
    status_temp = *status + one; //status_temp = 0,1,or 2.
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

//----------------------------------------------------------------------------

void PhysBAMIntersector::printFirstLayer(SubDomain& sub, SVec<double,3>&X, int TYPE)
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

void DistPhysBAMIntersector::finishStatusByPoints(IoData &iod)
{
  if(numFluid<3) //no need to look at points
    return;

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
    com->fprintf(stderr, "ERROR: (INTERSECTOR) Point-based initial conditions could not be found.\n");
    exit(-1);
  }

  list< pair<Vec3D,int> >::iterator iter;
  for(iter = Points.begin(); iter!=Points.end(); iter++)
    com->fprintf(stderr,"found point (%e %e %e) with FluidModel %d\n", (iter->first)[0], (iter->first)[1], (iter->first)[2], iter->second);
  

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
      if((*status)(iSub)[i]==PhysBAMIntersector::OUTSIDE) {
        (*status)(iSub)[i] = PhysBAMIntersector::UNDECIDED;
        nUndecided[iSub]++;
      }
    // 2. find seeds by points
    int nSeeds = intersector[iSub]->findSeedsByPoints(*(domain->getSubDomain()[iSub]), (*X)(iSub), Points, nUndecided[iSub]);
    // 3. flood fill if seeds are found
    if(nSeeds>0)
      intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
  }   

  while(1) { //get out only when all nodes are decided

    //1. check if all the nodes (globally) are determined
    total = 0;
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalMax(1,&total);
//    com->fprintf(stderr,"max of total = %d\n", total);
    if(total==0) //done!
      break;

    //2. try to get a seed from neighbor subdomains
    status_temp = *status + one; //status_temp = 0,1,2,3,....
    domain->assemble(domain->getLevelPat(),status_temp);
    status_temp -= one;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_temp(iSub), (*poly)(iSub), nUndecided[iSub]);

//      fprintf(stderr,"CPU %d: nUndecided = %d, newSeeds = %d\n",com->cpuNum(), nUndecided[iSub], nNewSeed);

      if(nNewSeed)
        intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }

  }
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::finishStatusByHistory(SubDomain& sub)
{
  if(numOfFluids()<3) //no need to do the following 
    return;

  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());
  for(int i=0; i<numNodes; i++) {
    if(status[i]==OUTSIDE){
      if(status0[i]!=INSIDE)
        status[i] = status0[i];
      else //status0[i]==INSIDE    
        for(int iNei=0; iNei<nToN.num(i); iNei++)
          if(status0[nToN[i][iNei]]!=INSIDE)
            status[i] = status0[nToN[i][iNei]];
    } else if(status[i]==UNDECIDED)
        fprintf(stderr,"Still have undecided nodes...\n");
  }
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::floodFill(SubDomain& sub, int& nUndecided)
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
  for(int i = 0; i < numNodes; ++i)
    if(status[i] != UNDECIDED) {
      seed[lead++] = i;
      level[i] = 0;
    } else
      nUndecided++;

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

PhysBAMIntersector::PhysBAMIntersector(SubDomain &sub, SVec<double,3> &X,
                    Vec<int> &stat, Vec<int> &stat0, DistPhysBAMIntersector &distInt) :
                      distIntersector(distInt), status(stat), status0(stat0),
                      edges(sub.getEdges()), globIndex(sub.getGlobSubNum())
{
  int numEdges = edges.size();

  status = UNDECIDED;
  status0 = UNDECIDED;

  locToGlobNodeMap = sub.getNodeMap();

  nFirstLayer = 0;
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::reset()
{
  status0 = status;
  status = UNDECIDED;
  nFirstLayer = 0;
}

//----------------------------------------------------------------------------

/** Find the closest structural triangle for each node. If no triangle intersect the bounding box of the node,
* no closest triangle exists
*/
void PhysBAMIntersector::getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist) 
{
  int ntri = distIntersector.getNumStructElems();
  MyTriangle *myTris = new MyTriangle[ntri];
  int (*triNodes)[3] = distIntersector.stElem;
  Vec3D *structX = (distIntersector.getStructPosition()).data();
  for(int i = 0; i < ntri; ++i)
    myTris[i] = MyTriangle(i, distIntersector.getStructPosition(), triNodes[i]);

  int nClose = 0;
  double maxDist = 0, maxSize = 0;

  int trBased = 0, edgeBased = 0, vertexBased = 0;

  KDTree<MyTriangle> structureTree(ntri, myTris);
  ClosestTriangle closestTriangle(triNodes, structX, distIntersector.triNorms);

  double maxErr = 0;
  double exactPhi = -1000;
  double numPhi = -1000;
  int errNode = -1;
  int maxCand = 0;
  int ndMaxCand = -1;

  int nMaxCand = 500;
  MyTriangle *candidates = new MyTriangle[nMaxCand];
  for(int i = 0; i < X.size(); ++i) {
    int nFound = structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates, nMaxCand);
    if(nFound > nMaxCand) {
      std::cerr << "For Fluid node " << locToGlobNodeMap[i]+1 << ", there were more candidates than we can handle: " << nFound << std::endl;
      nMaxCand = nFound;
      delete [] candidates;
      candidates = new MyTriangle[nMaxCand];
      structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates, nMaxCand);
    }
    if(nFound > maxCand) {
      maxCand = nFound;
      ndMaxCand = i;
    }
    closestTriangle.start(X[i]);
    for(int j = 0; j < std::min(nMaxCand,nFound); ++j) {
      double xi1, xi2, dist;

      closestTriangle.checkTriangle(candidates[j].trId());
    }
    if(nFound <= 0)
      tId[i] = -1;
    else {
      tId[i] = closestTriangle.bestTriangle();
      if(tId[i] < 0)
        std::cout << "Horror!!!" << std::endl;
      dist[i] = closestTriangle.signedDistance();
      if(closestTriangle.mode == 0 && dist[i] != 0)
        trBased++;
      if(closestTriangle.mode == 1 && dist[i] < 0)
              edgeBased++;
      if(closestTriangle.mode == 2 && dist[i] != 0)
              vertexBased++;
      nClose++;
      maxDist = std::max(maxDist, std::abs(dist[i]));
      maxSize = std::max(maxSize, sqrt(
          (boxMin[i][0]-boxMax[i][0])*(boxMin[i][0]-boxMax[i][0]) +
          (boxMin[i][1]-boxMax[i][1])*(boxMin[i][1]-boxMax[i][1]) +
          (boxMin[i][2]-boxMax[i][2])*(boxMin[i][2]-boxMax[i][2]) ));
    }
  }

  delete [] candidates;
  delete [] myTris;
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist)
{
  const double TOL = 0;
  for (int i=0; i<tId.size(); i++) {
    if (tId[i]<0)
      continue;
    status[i] = (dist[i]>=TOL) ? INSIDE : OUTSIDE;
    nFirstLayer++;
  }
}

//----------------------------------------------------------------------------

int PhysBAMIntersector::findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided)
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

int PhysBAMIntersector::findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list<pair<Vec3D,int> > P, int& nUndecided)
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

void PhysBAMIntersector::noCheckFloodFill(SubDomain& sub, int& nUndecided)
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

void PhysBAMIntersector::findIntersections(SVec<double,3>&X)
{
  int (*ptr)[2] = edges.getPtr();
  const double TOL = 1.0e-4;
  int MAX_ITER = 20;
  int max_iter = 0;
  double maxEdgeSize = 0;
  for (int l=0; l<edges.size(); l++) {
    int p = ptr[l][0], q = ptr[l][1];
    if(status[p]==status[q]) continue;

    //now need to find an intersection for this edge .
    Vec3D xp(X[p]), xq(X[q]);
    Vec3D dir = xq - xp;
    Vec3D xpPrime, xqPrime;
    maxEdgeSize = std::max(maxEdgeSize, dir.norm());

    ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(2);
    edgeRes(1).x[1] = 1;  edgeRes(1).x[2] = 2;
    edgeRes(2).x[1] = 2;  edgeRes(2).x[2] = 1;

    for (int iter=0; iter<MAX_ITER; iter++) {
      double coeff = iter*iter*TOL;
      Vec3D xpPrime = xp - coeff*dir;
      Vec3D xqPrime = xq + coeff*dir;
      ARRAY<VECTOR<double,3> > xyz(2);
      xyz(1)[1] = xpPrime[0];  xyz(1)[2] = xpPrime[1];  xyz(1)[3] = xpPrime[2];
      xyz(2)[1] = xqPrime[0];  xyz(2)[2] = xqPrime[1];  xyz(2)[3] = xqPrime[2];

      distIntersector.getInterface().Intersect(xyz, edgeRes, coeff*dir.norm());

      if (edgeRes(1).y.triangleID>=0 || edgeRes(2).y.triangleID>=0) {
        CrossingEdgeRes[l] = edgeRes(1).y;
        ReverseCrossingEdgeRes[l] = edgeRes(2).y;
        if (iter>max_iter) max_iter = iter;
        break;
      }
    }

    if (edgeRes(1).y.triangleID<0 && edgeRes(2).y.triangleID<0)
    fprintf(stderr,"ERROR: failed to get an intersection between node %d(%d) and %d(%d). \n",
                    locToGlobNodeMap[p]+1,status[p],locToGlobNodeMap[q]+1,status[q]);
  }
}

//----------------------------------------------------------------------------

LevelSetResult
PhysBAMIntersector::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
  int edgeNum = edges.find(ni, nj);
  if (!edgeIntersectsStructure(0.0,ni,nj)) {
    fprintf(stderr,"There is no intersection between node %d(status:%d) and %d(status:%d)! Abort...\n",
                   ni, status[ni], nj, status[nj]);
    exit(-1);
  }

  IntersectionResult<double> result; //need to determine which result to choose.
  IntersectionResult<double> rij = CrossingEdgeRes[edgeNum];
  IntersectionResult<double> rji = ReverseCrossingEdgeRes[edgeNum];
  double alpha0 = 0.0;

  if (rij.triangleID>=0 && rij.triangleID>=0 && rij.triangleID==rji.triangleID) {
    result = rij;
    alpha0 = (ni<nj)? result.alpha : 1.0-result.alpha;
  }  
  else if (rij.triangleID>=0 && rji.triangleID>=0 && rij.triangleID!=rji.triangleID) {
    result = (ni<nj) ? rij : rji;
    alpha0 = result.alpha;
  }
  else if (rij.triangleID>=0 && rji.triangleID<0) {
    result = rij;
    alpha0 = (ni<nj)? result.alpha : 1.0-result.alpha; 
  }
  else if (rij.triangleID<0 && rji.triangleID>=0) {
    result = rji;
    alpha0 = (ni<nj)? 1.0-result.alpha : result.alpha;
  }
  else //we really have no intersection for this edge!
    fprintf(stderr,"ERROR: intersection between %d and %d can not be detected.\n", ni, nj);

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

bool PhysBAMIntersector::isActive(double t, int n, int phase) const{
  return status[n] == phase;
}

//----------------------------------------------------------------------------

bool PhysBAMIntersector::wasActive(double t, int n, int phase) const{
  return status0[n] == phase;
}

//----------------------------------------------------------------------------

bool PhysBAMIntersector::edgeIntersectsStructure(double t, int ni, int nj) const {
  return status[ni] != status[nj];
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::projection(Vec3D x0, int tria, double& xi1, double& xi2, double& dist)
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

double PhysBAMIntersector::isPointOnSurface(Vec3D pt, int N1, int N2, int N3) 
{
  Vec<Vec3D> &solidX = distIntersector.getStructPosition();
  Vec3D X1 = solidX[N1];
  Vec3D X2 = solidX[N2];
  Vec3D X3 = solidX[N3];

  Vec3D normal = (X2-X1)^(X3-X1);
  normal /=  normal.norm();

  return fabs((pt-X1)*normal);
}


