#include <stdio.h>
#include <iostream>
#include "PhysBAMIntersect.h"
#include "Vector3D.h"
#include "Communicator.h"
#include <Domain.h>
#include <Vector.h>
#include <DistVector.h>
#include <Timer.h>
#include "LevelSet/IntersectionFactory.h"
#include "parser/Assigner.h"
#include "Geometry/KDTree.h"
#include <Connectivity.h>
#include <vector>
#include <queue>

using std::vector;
using std::pair;

static Timer *timer;
const int PhysBAMIntersector::UNDECIDED, PhysBAMIntersector::INSIDE, PhysBAMIntersector::OUTSIDE;

class PhysBAMIntersectorConstructor : public IntersectorConstructor {
    const char *structureFile;
    double tolIntersect;

  public:
    PhysBAMIntersectorConstructor() {
      structureFile = 0;
      tolIntersect = 1e-3;
    }

    DistLevelSetStructure *getIntersector(IntersectProblemData&) {
      DistPhysBAMIntersector *inter = new DistPhysBAMIntersector(tolIntersect);
      std::string solidSurface = structureFile;
      inter->init(solidSurface);

      return inter;
     // return 0;
    }

    int print();

    void init(ParseTree &dataTree) {
        std::cout << "inside the init function" << std::endl;
        ClassAssigner *ca = new ClassAssigner("PhysBAMIntersectorConstructor", 2, 0);
        new ClassStr<PhysBAMIntersectorConstructor>(ca, "structureFile", this, &PhysBAMIntersectorConstructor::structureFile);
        new ClassDouble<PhysBAMIntersectorConstructor>(ca, "tolerance", this, &PhysBAMIntersectorConstructor::tolIntersect);
        dataTree.implement(ca);
    }


};

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

  static const int maxNPairs = 32;

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
static int neg = 0, pos = 0;
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
    if(dist < 0) neg++;
    else pos++;
   if(xi1 < -eps)
      checkEdge(trId, nd[1], nd[2], nd[0], dist);
    if(xi2 < -eps)
      checkEdge(trId, nd[0], nd[2], nd[1], dist);
    if(xi3 < -eps)
      checkEdge(trId, nd[0], nd[1], nd[2], dist);
  }
}

FILE *eNodes = 0;
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
int agree=0, disagree=0;
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

double ClosestTriangle::getSignedVertexDistance() const {
  /* if(eNodes == 0)
    eNodes = fopen("EdgeNodes", "w");
  fprintf(eNodes, "%d %f %f %f %f %f %f\n", n1+1, x[0], x[1], x[2], structX[n1][0], structX[n1][1], structX[n1][2]);
  fflush(eNodes);*/
  if(isConsistant)
    return isPositive ? minDist : -minDist;
  return minDist;
}


IntersectorConstructor *myIntersect =
  IntersectionFactory::registerClass("PhysBAM", new PhysBAMIntersectorConstructor());


int PhysBAMIntersectorConstructor::print() {
  return 0;
}



DistPhysBAMIntersector::DistPhysBAMIntersector(double tol) {
  com = IntersectionFactory::getCommunicator();
  tolerance = tol;
}

LevelSetStructure &
DistPhysBAMIntersector::operator()(int subNum) const {
  return *intersector[subNum];
}
/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistPhysBAMIntersector::init(std::string solidSurface) {
  // Read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(solidSurface.c_str(), "r");
  if (topFile == NULL) {fprintf(stderr, "topFile doesn't exist at all :(\n"); exit(1); }
  int len;
  len = fscanf(topFile,"%d %d", &length_solids_particle_list, &length_triangle_list);
  triangle_list = new int[length_triangle_list][3];
  solids_particle_list = new Vec3D[length_solids_particle_list];
  solidX = new Vec<Vec3D>(length_solids_particle_list, solids_particle_list);

  com->fprintf(stderr,"solid surface: %d nodes, %d elements.\n", length_solids_particle_list, length_triangle_list);

  int thisNode;
  for (int iNode=0; iNode<length_solids_particle_list; iNode++)
    len = fscanf(topFile, "%d %lf %lf %lf", &thisNode, &(solids_particle_list[iNode][0]),
        &(solids_particle_list[iNode][1]), &(solids_particle_list[iNode][2]));
  if (thisNode!=length_solids_particle_list) {fprintf(stderr,"error in loading surface from file *!\n"); exit(1);}

  int nothing;
  for (int iElem=0; iElem<length_triangle_list; iElem++) {
    len = fscanf(topFile, "%d %d %d %d %d", &thisNode, &nothing, &(triangle_list[iElem][0]), &(triangle_list[iElem][1]),
        &(triangle_list[iElem][2]));
    triangle_list[iElem][0]--; triangle_list[iElem][1]--; triangle_list[iElem][2]--;
  }
  if (thisNode!=length_triangle_list) {fprintf(stderr,"error in loading surface from file **!\n", thisNode); exit(1);}
  fclose(topFile);

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  /*com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) com->fprintf(stderr,"OK.\n");
  else exit(-1);*/

  getBoundingBox();

  initializePhysBAM();
}

void DistPhysBAMIntersector::getBoundingBox() {
   xMin = xMax = solids_particle_list[0][0];
   yMin = yMax = solids_particle_list[0][1];
   zMin = zMax = solids_particle_list[0][2];
   for(int i = 1; i < length_solids_particle_list; ++i) {
      xMin = std::min(xMin, solids_particle_list[i][0]);
      xMax = std::min(xMax, solids_particle_list[i][0]);
      yMin = std::min(yMin, solids_particle_list[i][0]);
      yMax = std::max(yMax, solids_particle_list[i][0]);
      zMin = std::max(zMin, solids_particle_list[i][0]);
      zMax = std::max(zMax, solids_particle_list[i][0]);
   }
}

bool DistPhysBAMIntersector::checkTriangulatedSurface() {
  if ((length_triangle_list==0) || (length_solids_particle_list==0))
    {fprintf(stderr,"Solid surface not loaded.\n"); return false;}
  int numOfEdges = length_triangle_list*3/2,iEdge=0;
  if ((numOfEdges-((double)length_triangle_list)*3.0/2.0>0.1) ||
      (numOfEdges-((double)length_triangle_list)*3.0/2.0<-0.1))
    {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
  int edgeOrientation[numOfEdges][2];//edgeOrientation[i][0,1] stores from, to nodes.
  for (int iTriangle=0; iTriangle<length_triangle_list; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = triangle_list[iTriangle][0];  to1 = triangle_list[iTriangle][1];  found1 = false;
    from2 = triangle_list[iTriangle][1];  to2 = triangle_list[iTriangle][2];  found2 = false;
    from3 = triangle_list[iTriangle][2];  to3 = triangle_list[iTriangle][0];  found3 = false;
    for (int j=0; j<iEdge; j++) {
      if ((edgeOrientation[j][0]==from1) && (edgeOrientation[j][1]==to1))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to1) && (edgeOrientation[j][1]==from1))
        found1 = true;
      if ((edgeOrientation[j][0]==from2) && (edgeOrientation[j][1]==to2))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to2) && (edgeOrientation[j][1]==from2))
        found2 = true;
      if ((edgeOrientation[j][0]==from3) && (edgeOrientation[j][1]==to3))
        {fprintf(stderr,"Not all triangles on SolidSurface point outward. exit.\n"); return false; }
      else if ((edgeOrientation[j][0]==to3) && (edgeOrientation[j][1]==from3))
        found3 = true;
    }
    if (!found1) {
      if (iEdge>=numOfEdges)
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from1;  edgeOrientation[iEdge][1] = to1; iEdge++;
    }
    if (!found2) {
      if (iEdge>=numOfEdges)
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from2;  edgeOrientation[iEdge][1] = to2; iEdge++;
    }
    if (!found3) {
      if (iEdge>=numOfEdges)
        {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
      edgeOrientation[iEdge][0] = from3;  edgeOrientation[iEdge][1] = to3; iEdge++;
    }
  }
  return true;
}

void
DistPhysBAMIntersector::buildSolidNormals() {
  triNorms = new Vec3D[length_triangle_list];
  triSize = new double[length_triangle_list];
  // Also look to determine a point inside the solid but away from the structure.
  double nrmMax = 0;
  int trMaxNorm = -1;
  for (int iTriangle=0; iTriangle<length_triangle_list; iTriangle++) {
    int n1 = triangle_list[iTriangle][0];
    int n2 = triangle_list[iTriangle][1];
    int n3 = triangle_list[iTriangle][2];
    double x1 = solids_particle_list[n1][0];
    double y1 = solids_particle_list[n1][1];
    double z1 = solids_particle_list[n1][2];
    double dx2 = solids_particle_list[n2][0]-x1;
    double dy2 = solids_particle_list[n2][1]-y1;
    double dz2 = solids_particle_list[n2][2]-z1;
    double dx3 = solids_particle_list[n3][0]-x1;
    double dy3 = solids_particle_list[n3][1]-y1;
    double dz3 = solids_particle_list[n3][2]-z1;

   // first calculate the "size" of the triangle.
    double dx23 = solids_particle_list[n3][0] - solids_particle_list[n2][0];
    double dy23 = solids_particle_list[n3][1] - solids_particle_list[n2][1];
    double dz23 = solids_particle_list[n3][2] - solids_particle_list[n2][2];
    double size12 = Vec3D(dx2,dy2,dz2).norm();
    double size13 = Vec3D(dx3,dy3,dz3).norm();
    double size23 = Vec3D(dx23,dy23,dz23).norm();
    triSize[iTriangle] = min(size12,min(size13,size23));

    // now calculate the normal.
    triNorms[iTriangle] = Vec3D(dx2, dy2, dz2)^Vec3D(dx3,dy3,dz3);
    double nrm = triNorms[iTriangle].norm();
    if(nrm > nrmMax) {
      nrmMax = nrm;
      trMaxNorm = iTriangle;
    }
    if(nrm != 0)
       triNorms[iTriangle] /= nrm;
  }
  if(trMaxNorm >= 0) {
    int n1 = triangle_list[trMaxNorm][0];
    int n2 = triangle_list[trMaxNorm][1];
    int n3 = triangle_list[trMaxNorm][2];
    Vec3D trCenter =
      Vec3D(solids_particle_list[n1][0]+solids_particle_list[n2][0]+solids_particle_list[n3][0],
            solids_particle_list[n1][1]+solids_particle_list[n2][1]+solids_particle_list[n3][1],
            solids_particle_list[n1][2]+solids_particle_list[n2][2]+solids_particle_list[n3][2])/3;
    // offset the center point by a small amount, but bigger than the epsilon used.
     Vec3D p1 = trCenter - 2*tolerance*triNorms[trMaxNorm];
     double maxDist = Vec3D(xMax-xMin, yMax-yMin, zMax-zMin).norm();
     Vec3D p2 = p1 - maxDist*triNorms[trMaxNorm];

     LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(2);
     edgeRes(1).x[1] = 1;
     edgeRes(1).x[2] = 2;
     edgeRes(2).x[1] = 2;
     edgeRes(2).x[2] = 1;

     LIST_ARRAY<VECTOR<double,3> > xyz(2);
     xyz(1)[1] = p1[0];
     xyz(1)[2] = p1[1];
     xyz(1)[3] = p1[2];
     xyz(2)[1] = p2[0];
     xyz(2)[2] = p2[1];
     xyz(2)[3] = p2[2];
     fprintf(stderr, "P1: %f %f %f, P2: %f %f %f\n", p1[0], p1[1], p1[2],
           p2[0], p2[1], p2[2]);
     getInterface().Intersect(xyz, edgeRes,getTolerance());
     fprintf(stderr, "Starting from %d, reached %d and %d\n",
              trMaxNorm+1, edgeRes(1).y.triangleID, edgeRes(2).y.triangleID);
     if(edgeRes(1).y.triangleID < 0)
       fprintf(stderr, "ERROR: OPEN SURFACE\n");
     insidePoint = 0.5*((1+edgeRes(1).y.alpha)*p1+(1-edgeRes(1).y.alpha)*p2);
     fprintf(stderr, "Distance is: %f vs max %f\n", (insidePoint-p1).norm(), maxDist);
     fprintf(stderr, "%f Inside point: %f %f %f\n", edgeRes(1).y.alpha, insidePoint[0], insidePoint[1], insidePoint[2]);
  } else
    fprintf(stderr, "All triangles are degenerate!!\n");
}


void
DistPhysBAMIntersector::initializePhysBAM() {
// Initialize the Particles list
  PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >& physbam_solids_particle=*new PhysBAM::SOLIDS_PARTICLE<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle.Add_Particles(length_solids_particle_list);
  for (int i=0; i<length_solids_particle_list; i++) {
    physbam_solids_particle.X(i+1) = PhysBAM::VECTOR<double,3>(solids_particle_list[i][0],
        solids_particle_list[i][1], solids_particle_list[i][2]);
  }
//  std::cout<<physbam_solids_particle.X.Size()<<std::endl;

  // Initialize the Triangle list
  PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> > & physbam_triangle_list=*new PhysBAM::LIST_ARRAY<PhysBAM::VECTOR<int,3> >();
  for (int i=0; i<length_triangle_list; i++){
    int nx, ny, nz;
    nx = triangle_list[i][0] + 1;  ny = triangle_list[i][1] + 1;  nz = triangle_list[i][2] + 1;
    physbam_triangle_list.Append(PhysBAM::VECTOR<int,3>(nx, ny, nz));
  }

  // Construct TRIANGLE_MESH triangle_mesh.
  PhysBAM::TRIANGLE_MESH& physbam_triangle_mesh=*new PhysBAM::TRIANGLE_MESH(physbam_solids_particle.number, physbam_triangle_list);
  physbam_triangle_mesh.Initialize_Adjacent_Elements();
  // Construct TRIANGULATED_SURFACE.
  PhysBAM::TRIANGULATED_SURFACE<double>& physbam_triangulated_surface=*new PhysBAM::TRIANGULATED_SURFACE<double>(physbam_triangle_mesh, physbam_solids_particle);
  physbam_triangulated_surface.Update_Triangle_List();
  std::cout <<"Going to make PhysBAMInterface" << std::endl;
  physInterface = new PhysBAMInterface<double>(physbam_triangulated_surface);
  std::cout <<"Done making PhysBAMInterface" << std::endl;
  // Compute the normals of each of the structural triangles
  buildSolidNormals();
}

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistPhysBAMIntersector::initialize(Domain *d, DistSVec<double,3> &X) {
  this->X = &X;
  domain = d;
  timer = domain->getTimer();
  numLocSub = d->getNumLocSub();
  intersector = new PhysBAMIntersector*[numLocSub];
  pseudoPhi = new DistVec<double>(X.info());

  // for getClosestTriangles
  DistSVec<double,3> boxMax(X.info());
  DistSVec<double,3> boxMin(X.info());
  DistVec<double> distance(X.info());
  DistVec<int> tId(X.info());

  d->findNodeBoundingBoxes(X,boxMin,boxMax);

  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new PhysBAMIntersector(*(d->getSubDomain()[i]), X(i), *this);
    intersector[i]->getClosestTriangles(X(i), boxMin(i), boxMax(i), tId(i), distance(i));

    intersector[i]->computeFirstLayerNodeStatus(tId(i), distance(i));
    intersector[i]->fixUntouchedSubDomain(X(i));
    intersector[i]->finishNodeStatus(*(d->getSubDomain()[i]), X(i));
    intersector[i]->findIntersections(X(i));
  }
}

//----------------------------------------------------------------------------

PhysBAMIntersector::PhysBAMIntersector(SubDomain &sub, SVec<double,3> &X,
                    DistPhysBAMIntersector &distInt) :
 distIntersector(distInt), status(sub.numNodes()),
 edges(sub.getEdges()), globIndex(sub.getGlobSubNum())
{
  int numEdges = edges.size();

  status = UNDECIDED;

  nodeMap = sub.getNodeMap();

  locToGlobNodeMap = sub.getNodeMap();

  nFirstLayer = 0;
}

bool verb =  false;
//----------------------------------------------------------------------------

/** Find the closest structural triangle for each node. If no triangle intersect the bounding box of the node,
* no closest triangle exists
*/
void PhysBAMIntersector::getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist) {
  int ntri = distIntersector.length_triangle_list;
  MyTriangle *myTris = new MyTriangle[ntri];
  int (*triNodes)[3] = distIntersector.triangle_list;
  Vec3D *structX = distIntersector.solids_particle_list;
  for(int i = 0; i < ntri; ++i)
    myTris[i] = MyTriangle(i, *distIntersector.solidX, triNodes[i]);

  int nClose = 0;
  double maxDist = 0, maxSize = 0;

  int trBased = 0, edgeBased = 0, vertexBased = 0;

  KDTree<MyTriangle> structureTree(ntri, myTris);
  ClosestTriangle closestTriangle(distIntersector.triangle_list, structX, distIntersector.triNorms);

  double maxErr = 0;
  double exactPhi = -1000;
  double numPhi = -1000;
  int errNode = -1;
  int maxCand = 0;
  int ndMaxCand = -1;
  for(int i = 0; i < X.size(); ++i) {
    MyTriangle candidates[500];
    verb = i == 2410;
    int nFound = structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates, 500);
    if(nFound > 500)
      std::cerr << "There were more candidates than we can handle: " << nFound << std::endl;
    if(nFound > maxCand) {
      maxCand = nFound;
      ndMaxCand = i;
    }
    closestTriangle.start(X[i]);
    /*if(i == 9130)
      std::cout << "Box: " << boxMin[i][0] << " " << boxMin[i][1] << " "<< boxMin[i][2] << " "
         << boxMax[i][0] << " " << boxMax[i][1] << " "<< boxMax[i][2] << std::endl;*/
    for(int j = 0; j < std::min(500,nFound); ++j) {
      double xi1, xi2, dist;

      closestTriangle.checkTriangle(candidates[j].trId());
    /*  if(i == 9130) {
        std::cout << " " << candidates[j].trId()<< " box " << candidates[j].val(0) << " " << candidates[j].val(1) << " " << candidates[j].val(2) <<
                " " << candidates[j].val(0)+candidates[j].width(0) << " " << candidates[j].val(1)+candidates[j].width(1) << " " << candidates[j].val(2)+candidates[j].width(2) << std::endl;
      }*/
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
      double a = 0.685;
      double b = a/6;
      double f = (X[i][0]-a)*(X[i][0]-a)/(a*a)+X[i][1]*X[i][1]/(b*b)+X[i][2]*X[i][2]/(b*b)-1;
      Vec3D gradF(2*(X[i][0]-a)/(a*a), 2*X[i][1]/(b*b), 2*X[i][2]/(b*b));
      double phiExact = f/gradF.norm();
      if(std::abs(phiExact-dist[i]) > maxErr && phiExact*dist[i] < 0) {
        maxErr= std::abs(phiExact-dist[i]);
        exactPhi = phiExact;
        numPhi = dist[i];
        errNode = i;
      }
    }
  }
  std::cout << "Close nodes: " << nClose << " Maximum distance: " << maxDist << " Max size: " << maxSize << std::endl;
  std::cout << "Neg probs: " << neg << " Pos probs: " << pos << std::endl;
  std::cout << "triangle Based: " << trBased << " edge based: " << edgeBased << " vertex based: " << vertexBased << std::endl;
  std::cout << "Maximum error: " << maxErr << " exact: " << exactPhi << " numerical: " << numPhi << " node: " << errNode <<
  " " << X[errNode][0] << " " << X[errNode][1] << " " << X[errNode][2] << std::endl;
  std::cout << "Fluid box: " << boxMin[errNode][0] << " " << boxMin[errNode][1] << " "<< boxMin[errNode][2] << " "
  << boxMax[errNode][0] << " " << boxMax[errNode][1] << " "<< boxMax[errNode][2] << std::endl;

  std::cout << "Max candidates: " << maxCand << " at " << ndMaxCand << std::endl;
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist)
{
  const double TOL = 0;
  int nInside=0, nOutside=0;
  for (int i=0; i<tId.size(); i++) {
    if (tId[i]<0)
      continue;
    status[i] = (dist[i]>=TOL) ? INSIDE : OUTSIDE;
    if(status[i] == INSIDE)
      nInside++;
    else
      nOutside++;
    nFirstLayer++;
  }
  std::cout << "Number inside: " << nInside << " outside: " << nOutside << " out of " << tId.size() << std::endl;
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::fixUntouchedSubDomain(SVec<double,3>&X)
{
  if(nFirstLayer == 0) {
     LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(1);
     edgeRes(1).x[1] = 1;
     edgeRes(1).x[2] = 2;

     Vec3D insidePoint = distIntersector.getInsidePoint();
     LIST_ARRAY<VECTOR<double,3> > xyz(2);
     xyz(1)[1] = X[0][0];
     xyz(1)[2] = X[0][1];
     xyz(1)[3] = X[0][2];
     xyz(2)[1] = insidePoint[0];
     xyz(2)[2] = insidePoint[1];
     xyz(2)[3] = insidePoint[2];
     fprintf(stderr, "Test edge is %f %f %f, %f %f %f\n",
           X[0][0], X[0][1], X[0][2],
           insidePoint[0], insidePoint[1], insidePoint[2]);
     distIntersector.getInterface().Intersect(xyz, edgeRes,distIntersector.getTolerance());
     if(edgeRes(1).y.triangleID >= 0) {
       Vec3D edgeVec(insidePoint[0]-X[0][0], insidePoint[1]-X[0][1], insidePoint[2]-X[0][2]);
       const Vec3D &trNorm = distIntersector.getSurfaceNorm(edgeRes(1).y.triangleID-1);
       if(edgeVec*trNorm < 0) {
         fprintf(stderr, "This subdomain is outside the structure\n");
         status[0] = INSIDE;
       } else {
         fprintf(stderr, "This subdomain is inside the structure\n");
         status[0] = OUTSIDE;
       }
     } else {
       fprintf(stderr, "This subdomain is trivially inside the structure\n");
       status[0] = OUTSIDE;
     }
  }
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::finishNodeStatus(SubDomain& sub, SVec<double,3>&X)
{
  int numNodes = status.size();
  Connectivity &nToN = *(sub.createEdgeBasedConnectivity());

  // Propogate status to the entire subdomain.
  // List is used as a FIFO queue
  Vec<int> list(numNodes);
  Vec<int> level(numNodes);
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i)
    if(status[i] != UNDECIDED) {
      list[lead++] = i;
      level[i] = 0;
    }
  std::cout << "Initial lead (# decided nodes): " << lead << std::endl;
  while(next < lead) {
    int cur = list[next++];
    int curStatus = status[cur];
    int curLevel = level[cur];
    for(int i = 0; i < nToN.num(cur); ++i) {
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        level[nToN[cur][i]] = curLevel+1;
        list[lead++] = nToN[cur][i];
      } else {
        if(status[nToN[cur][i]] != curStatus && ( curLevel != 0 || level[nToN[cur][i]] != 0))
          std::cerr << "Incompatible nodes have met: " << cur << " and " << nToN[cur][i] <<
             " " << curLevel << " " << level[nToN[cur][i]] << std::endl;
      }
    }
  }

  // calculate the number of inside/outside/undecided nodes.
  int nInside = 0, nOutside = 0, nUndecided = 0;
  for(int i = 0; i < numNodes; ++i) {
    if(status[i] == INSIDE)        nInside++;
    else if(status[i] == OUTSIDE)  nOutside++;
    else                           nUndecided++;
  }

  if (nUndecided)
    fprintf(stderr,"WARNING: In Subdomain %d, there are %d nodes that are UNDECIDED!\n", sub.getGlobSubNum(), nUndecided);

  //debug only: Plot the first layer of inside nodes.
/*  int (*ptr)[2] = edges.getPtr();
  FILE* firstLayer = fopen("firstLayer.top","w");
  fprintf(firstLayer, "Nodes InsideNodes\n");
  for (int i=0; i<sub.numNodes(); i++)
    if (status[i]==OUTSIDE) fprintf(firstLayer,"%d %e %e %e\n", i+1, X[i][0], X[i][1], X[i][2]);
  fprintf(firstLayer, "Elements FirstLayer using InsideNodes\n");
  for (int l=0; l<edges.size(); l++){
    int x1 = ptr[l][0], x2 = ptr[l][1];
    if (status[x1]!=OUTSIDE || status[x2]!=OUTSIDE) continue;
    int crit = 0;
    for (int i=0; i<nToN.num(x1); i++)
      if (status[nToN[x1][i]]==INSIDE) {crit++; break;}
    for (int i=0; i<nToN.num(x2); i++)
      if (status[nToN[x2][i]]==INSIDE) {crit++; break;}
    if (crit==2)
      fprintf(firstLayer,"%d %d %d %d\n", l+1, (int)1, x1+1, x2+1);
  }
  fclose(firstLayer);
*/
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::findIntersections(SVec<double,3>&X)
{
  int (*ptr)[2] = edges.getPtr();
  const double TOL = 1.0e-3;
  int MAX_ITER = 50;
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

    LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(2);
    edgeRes(1).x[1] = 1;  edgeRes(1).x[2] = 2;
    edgeRes(2).x[1] = 2;  edgeRes(2).x[2] = 1;

    for (int iter=0; iter<MAX_ITER; iter++) {
      double coeff = iter*iter*TOL;
      Vec3D xpPrime = xp - coeff*dir;
      Vec3D xqPrime = xq + coeff*dir;
      LIST_ARRAY<VECTOR<double,3> > xyz(2);
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
    fprintf(stderr,"ERROR: failed to get an intersection between node %d and %d. \n",
                    locToGlobNodeMap[p]+1,locToGlobNodeMap[q]+1);
  }
 
  std::cout << "Maximum edge distance: " << maxEdgeSize << std::endl;
  fprintf(stderr,"In subdomain %d: maximum iteration: %d. maximum tolerance: %e.\n", globIndex, max_iter, (max_iter+1)*(max_iter+1)*TOL);
}

//----------------------------------------------------------------------------

/*
void PhysBAMIntersector::computeLocalPseudoPhi(SVec<double,3> &X, Vec<int> &tId, Vec<double> &dist,
                                               SVec<double,3> &normApprox, Vec<double> &weightSum)
{
  bool* edgeMasterFlag = edges.getMasterFlag();

  int numEdges = edges.size();
  int numNodes = status.size();

  phi = 0;
  normApprox = 0;
  weightSum = 0;

  // Compute the contribution of each intersection to pseudoPhi, the normals and the weights.
  // Add intersecting edges to the list of reverse edges
  for(int i = 0; i < edgeRes.Size(); ++i) {
    //if (!edgeMasterFlag[i]) continue;
    // check if this edge intersects the structure
    if(edgeRes(i+1).y.triangleID >= 0) {
      int p = edgeRes(i+1).x[1]-1, q = edgeRes(i+1).x[2]-1;
      // Add the reverse edge to the list of edges to compute
      reverseEdges.push_back(pair<int,int>(q,p));

      if(edgeMasterFlag[i])
        updatePhi(p, q, edgeRes(i+1).y, X, phi, normApprox, weightSum);
      //trIntersectCount[edgeRes(i+1).y.triangleID-1]++;
      int trId = edgeRes(i+1).y.triangleID-1;
      int *nd = distIntersector.triangle_list[trId];
      nIntersect++;
    }
  }
  std::cout << "Number of intersections: " << nIntersect << " vs " << numEdges << " in " << (t-t0) << std::endl;
  //std::cout << "Nd min: " << trNdIsectCount.min() << std::endl;

  LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > reverseEdgeRes(nIntersect);
  for(int i = 0; i < nIntersect; ++i) {
      reverseEdgeRes(i+1).x[1] = reverseEdges[i].first+1;
      reverseEdgeRes(i+1).x[2] = reverseEdges[i].second+1;
  }

  t0 = timer->getTime();
  // First call to intersect to find all intersections
  distIntersector.getInterface().Intersect(xyz, reverseEdgeRes,distIntersector.getTolerance());
  t = timer->getTime();

  nIntersect = 0;
  for(int i = 0; i < reverseEdges.size(); ++i) {
    int p = reverseEdgeRes(i+1).x[1]-1, q = reverseEdgeRes(i+1).x[2]-1;
    secondIntersection[edges.find(p,q)] = reverseEdgeRes(i+1).y;

    if(reverseEdgeRes(i+1).y.triangleID >= 0) {
      if(edgeMasterFlag[edges.find(q,p)])
        updatePhi(p, q, reverseEdgeRes(i+1).y, X, phi, normApprox, weightSum);
      //trIntersectCount[edgeRes(i+1).y.triangleID-1]++;
      nIntersect++;
    } else {
      std::cout << "Reverse between " << p << " and " << q << " has no intersection" << std::endl;
      LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(1);
      edgeRes(1).x[1] = 1;
      edgeRes(1).x[2] = 2;

      LIST_ARRAY<VECTOR<double,3> > xyz(2);
      xyz(1)[1] = X[p][0];
      xyz(1)[2] = X[p][1];
      xyz(1)[3] = X[p][2];
      xyz(2)[1] = X[q][0];
      xyz(2)[2] = X[q][1];
      xyz(2)[3] = X[q][2];
      distIntersector.getInterface().Intersect(xyz, edgeRes,distIntersector.getTolerance()+1e-9);
      std::cout << "redoing it gives: " << edgeRes(1).y.triangleID << std::endl;
      if(edgeRes(1).y.triangleID < 0)
        std::cerr << "Reverse cut could not be fixed! This will probably lead to a crash or incorrect results!" << std::endl;
      else {
        if(edgeMasterFlag[edges.find(q,p)])
          updatePhi(p, q, edgeRes(1).y, X, phi, normApprox, weightSum);
        secondIntersection[edges.find(p,q)] = edgeRes(1).y;
      }
    }
  }
  int nZeros = 0;
}

void PhysBAMIntersector::finishPseudoPhi(SubDomain &sub, SVec<double,3> &X, SVec<double,3> &normApprox,
                                         Vec<double> &weightSum)
 {
  int numNodes = status.size();

  for(int i = 0; i < numNodes; ++i) {
    if(weightSum[i] > 0) {
      phi[i] /= weightSum[i];
      status[i] = (phi[i] >= 0) ? INSIDE : OUTSIDE;
      double nl = Vec3D(normApprox[i]).norm();
      if(nl != 0)
        locNorm[i] = Vec3D(normApprox[i]) / nl;
      else
        std::cout << "Really null norm" << std::endl;
    }
  }
  if(nIntersect == 0) {
     LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(1);
     edgeRes(1).x[1] = 1;
     edgeRes(1).x[2] = 2;

     Vec3D insidePoint = distIntersector.getInsidePoint();
     LIST_ARRAY<VECTOR<double,3> > xyz(2);
     xyz(1)[1] = X[0][0];
     xyz(1)[2] = X[0][1];
     xyz(1)[3] = X[0][2];
     xyz(2)[1] = insidePoint[0];
     xyz(2)[2] = insidePoint[1];
     xyz(2)[3] = insidePoint[2];
     fprintf(stderr, "Test edge is %f %f %f, %f %f %f\n",
           X[0][0], X[0][1], X[0][2],
           insidePoint[0], insidePoint[1], insidePoint[2]);
     distIntersector.getInterface().Intersect(xyz, edgeRes,distIntersector.getTolerance());
     if(edgeRes(1).y.triangleID >= 0) {
       Vec3D edgeVec(insidePoint[0]-X[0][0], insidePoint[1]-X[0][1], insidePoint[2]-X[0][2]);
       const Vec3D &trNorm = distIntersector.getSurfaceNorm(edgeRes(1).y.triangleID-1);
       if(edgeVec*trNorm < 0) {
         fprintf(stderr, "This subdomain is outside the structure\n");
         status[0] = INSIDE;
       } else {
         fprintf(stderr, "This subdomain is inside the structure\n");
         status[0] = OUTSIDE;
       }
     } else {
       fprintf(stderr, "This subdomain is trivially inside the structure\n");
       status[0] = OUTSIDE;
     }
  }

  Connectivity &nToN = *(sub.createEdgeBasedConnectivity());

  // Find problem nodes: undecided nodes that have connections to both inside and outside nodes
  std::vector<int> problemNodes;
  for(int cur = 0; cur < sub.numNodes(); ++cur) {
    if(status[cur] != UNDECIDED)
      continue;
    int v = 0;
    for(int i = 0; i < nToN.num(cur); ++i)
      if (status[nToN[cur][i]] != UNDECIDED)
        v |= (status[nToN[cur][i]] == INSIDE) ? 1 : 2;
    if(v == 3)
      problemNodes.push_back(cur);
  }


  for(int j = 0; j < problemNodes.size(); ++j) {
    int nIn = 0, nOut = 0;
    int nd = problemNodes[j];
    Vec3D thisNode(X[nd]);
    double phiThis, phiMin, phiMax;
    std::cerr << "Node " << nd << " is causing problems." << std::endl;
    double maxDist = 0, minDot, maxDot;
    double phiSum = 0;
    double wSum = 0;
    bool isInit = false;
    for(int i = 0; i < nToN.num(nd); ++i) {
      int oNd = nToN[nd][i];
      if(status[oNd] == UNDECIDED)
        continue;
      Vec3D otherNode(X[oNd]);
      double phiONd = phi[oNd];
      Vec3D oNdNormal = locNorm[oNd];
      double dot = (thisNode-otherNode)*oNdNormal;
      double phiEst = phiONd + dot;
      double dist = (thisNode-otherNode).norm();
      double weight = std::abs(dot/dist);
      phiSum += weight*phiEst;
      wSum += weight;
      if(!isInit) {
        phiThis = phiMin = phiMax = phiEst;
        maxDist = dist;
        maxDot = minDot = (thisNode-otherNode)*oNdNormal/dist;
        isInit = true;
      }
      if(std::abs(phiEst) < std::abs(phiThis))
        phiThis = phiEst;
      if(phiEst < phiMin) { phiMin = phiEst; minDot = (thisNode-otherNode)*oNdNormal/dist; }
      if(phiEst > phiMax) { phiMax = phiEst; maxDot = (thisNode-otherNode)*oNdNormal/dist; }
      maxDist = std::max(maxDist, dist);
    }
    if(phiMin * phiMax <= 0 || nd >= 0) {
      std::cerr << "WOW! " << nd << " " << phiSum/wSum << std::endl;
      std::cout << "Phi: "<<phiThis << " Min: " << phiMin << " Max : " << phiMax  << " minDot: " <<minDot << " maxDot: " << maxDot << std::endl;
      if(nd >= 0 ) {
        std::cout << "My status: " << status[nd] << " and weight " << weightSum[nd] << std::endl;
        for(int i = 0; i < nToN.num(nd); ++i) {
          int oNd = nToN[nd][i];
          if(oNd == nd)
            continue;
          int edgeNum = edges.find(nd, oNd);
          std::cout << oNd << " Status: " << status[oNd] << " triangle: " << edgeRes(edgeNum+1).y.triangleID << std::endl;
          if(status[oNd] == UNDECIDED)
            continue;
          Vec3D otherNode(X[oNd]);
          double phiONd = phi[oNd];
          Vec3D oNdNormal = locNorm[oNd];
          double dot = (thisNode-otherNode)*oNdNormal;
          double phiEst = phiONd + dot;
          double dist = (thisNode-otherNode).norm();
          double weight = std::abs(dot/dist);
        }
      }
    }

  }

  // List is used as a FIFO queue
  Vec<int> list(numNodes);
  Vec<int> level(numNodes);
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i)
    if(status[i] != UNDECIDED) {
      list[lead++] = i;
      level[i] = 0;
    }
  std::cout << "Initial lead: " << lead << std::endl;
  while(next < lead) {
    int cur = list[next++];
    int curStatus = status[cur];
    int curLevel = level[cur];
    for(int i = 0; i < nToN.num(cur); ++i) {
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        level[nToN[cur][i]] = curLevel+1;
        list[lead++] = nToN[cur][i];
      } else {
        if(status[nToN[cur][i]] != curStatus && ( curLevel != 0 || level[nToN[cur][i]] != 0))
          std::cerr << "Incompatible nodes have met: " << cur << " and " << nToN[cur][i] <<
             " " << curLevel << " " << level[nToN[cur][i]] << std::endl;
      }
    }
  }

  int nInside = 0, nOutside = 0, nUndecided = 0;
  for(int i = 0; i < sub.numNodes(); ++i) {
    if(status[i] == INSIDE)
      nInside++;
    else if(status[i] == OUTSIDE)
      nOutside++;
    else
      nUndecided++;
  }

  int iter = 0;
  while (nUndecided) {
    for (int i=0; i<sub.numNodes(); ++i)
      if (status[i]==UNDECIDED) {
        LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > edgeRes(1);
        edgeRes(1).x[1] = 1;
        edgeRes(1).x[2] = 2;

        Vec3D insidePoint = distIntersector.getInsidePoint();
        LIST_ARRAY<VECTOR<double,3> > xyz(2);
        xyz(1)[1] = X[i][0];
        xyz(1)[2] = X[i][1];
        xyz(1)[3] = X[i][2];
        xyz(2)[1] = insidePoint[0];
        xyz(2)[2] = insidePoint[1];
        xyz(2)[3] = insidePoint[2];
        fprintf(stderr, "Test edge is %f %f %f, %f %f %f\n",
              X[0][0], X[0][1], X[0][2],
              insidePoint[0], insidePoint[1], insidePoint[2]);
        distIntersector.getInterface().Intersect(xyz, edgeRes,distIntersector.getTolerance());

        if(edgeRes(1).y.triangleID >= 0) {
          Vec3D edgeVec(insidePoint[0]-X[0][0], insidePoint[1]-X[0][1], insidePoint[2]-X[0][2]);
          const Vec3D &trNorm = distIntersector.getSurfaceNorm(edgeRes(1).y.triangleID-1);
          if(edgeVec*trNorm < 0) {
            status[i] = INSIDE; nInside++;}
          else {status[i] = OUTSIDE; nOutside++;}
        } else {status[i] = OUTSIDE; nOutside++;}
        nUndecided--;

        next = 0, lead = 0;
        list[lead++] = i;
        level[i] = 0;
        while(next < lead) {
          int cur = list[next++];
          int curStatus = status[cur];
          int curLevel = level[cur];
          for(int k = 0; k < nToN.num(cur); ++k) {
            if(status[nToN[cur][k]] == UNDECIDED) {
              nUndecided--;
              if (curStatus==INSIDE) nInside++; else nOutside++;
              status[nToN[cur][k]] = curStatus;
              level[nToN[cur][k]] = curLevel+1;
              list[lead++] = nToN[cur][k];
            }
          }
        }

      }
  }

  // Supplemental check
  int numWeird = 0;
  for(int cur = 0; cur < sub.numNodes(); ++cur)
    for(int i = 0; i < nToN.num(cur); ++i)
      if(status[nToN[cur][i]] != status[cur])  {
        int edgeNum = edges.find(cur, nToN[cur][i]);
        if(edgeRes(edgeNum+1).y.triangleID < 0)
          numWeird++;
        if(weightSum[cur] == 0 && weightSum[nToN[cur][i]] == 0)
          std::cout << "Got an edge with no life! while tol = " <<distIntersector.getTolerance() << std::endl;
          //std::cout << "Found conflicting edge " << edgeNum << " between " << cur <<
          //   " and " << nToN[cur][i] << std::endl;
      }
  std::cout << "We have " << numWeird << " intersections to revisit" << std::endl;

}

void PhysBAMIntersector::updatePhi(int p, int q, IntersectionResult<double> &res, SVec<double, 3> &X, Vec<double> &phi,
                                   SVec<double, 3> &normApprox, Vec<double> &weightSum) {
  double alpha = res.alpha;
  int trID = res.triangleID-1;
  const Vec3D &trNorm = distIntersector.getSurfaceNorm(trID);
  Vec3D edgeVec(X[q][0]-X[p][0], X[q][1]-X[p][1], X[q][2]-X[p][2]);
  if(edgeVec.norm() == 0)
    return;
  Vec3D edgeDir = edgeVec / edgeVec.norm();
  double dot = trNorm*edgeVec;
  double weight = std::abs(trNorm*edgeDir)+1e-3;
  double phip = (alpha-1)*dot;
  phi[p] += weight*phip;
  for(int j = 0; j < 3; ++j)
    normApprox[p][j] += weight*trNorm[j];
  weightSum[p] += weight;

}
*/

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
  Vec3D nrm = distIntersector.getSurfaceNorm(trueTriangleID);

  LevelSetResult lsRes(nrm[0], nrm[1], nrm[2], 0, 0, 0);
  lsRes.alpha = alpha0;
  lsRes.xi[0] = result.zeta[0];
  lsRes.xi[1] = result.zeta[1];
  lsRes.xi[2] = 1-result.zeta[0]-result.zeta[1];
  lsRes.trNodes[0] = distIntersector.triangle_list[trueTriangleID][0];
  lsRes.trNodes[1] = distIntersector.triangle_list[trueTriangleID][1];
  lsRes.trNodes[2] = distIntersector.triangle_list[trueTriangleID][2];

  return lsRes;
}

//----------------------------------------------------------------------------

bool PhysBAMIntersector::isActive(double t, int n) {
  return status[n] == INSIDE;
}

//----------------------------------------------------------------------------

bool PhysBAMIntersector::edgeIntersectsStructure(double t, int ni, int nj) const {
  return status[ni] != status[nj];
}

//----------------------------------------------------------------------------

void PhysBAMIntersector::projection(Vec3D x0, int tria, double& xi1, double& xi2, double& dist)
{
  int iA = distIntersector.triangle_list[tria][0];
  int iB = distIntersector.triangle_list[tria][1];
  int iC = distIntersector.triangle_list[tria][2];
  Vec3D xA = distIntersector.solids_particle_list[iA];
  Vec3D xB = distIntersector.solids_particle_list[iB];
  Vec3D xC = distIntersector.solids_particle_list[iC];

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

