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
#include <Connectivity.h>
#include <vector>
#include <queue>

using std::vector;
using std::pair;
using std::map;

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

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
  triToTri = new int[length_triangle_list][3];

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
  com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) com->fprintf(stderr,"OK.\n");
  else exit(-1);

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

EdgePair DistPhysBAMIntersector::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

bool DistPhysBAMIntersector::checkTriangulatedSurface()
{
  map<iipair, ibpair> edgeMap;

  for (int iTriangle=0; iTriangle<length_triangle_list; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = triangle_list[iTriangle][0];  to1 = triangle_list[iTriangle][1];  found1 = false;
    from2 = triangle_list[iTriangle][1];  to2 = triangle_list[iTriangle][2];  found2 = false;
    from3 = triangle_list[iTriangle][2];  to3 = triangle_list[iTriangle][0];  found3 = false;

    EdgePair ep[3];
    ep[0] = makeEdgePair(from1, to1, iTriangle);
    ep[1] = makeEdgePair(from2, to2, iTriangle);
    ep[2] = makeEdgePair(from3, to3, iTriangle);
    for(int i=0; i < 3; ++i) {
      map<iipair, ibpair>::iterator it = edgeMap.find(ep[i].first);
      if(it != edgeMap.end()) { // we found this edge
         if(it->second.second == ep[i].second.second)
           {fprintf(stderr,"triangulated surface is not closed. exit.\n"); return false;}
         else {
             int oTriangle = it->second.first;
             triToTri[iTriangle][i] = oTriangle;
             int n1 = it->second.second ? ep[i].first.first : ep[i].first.second;
             int edgeIndex;
             if(triangle_list[oTriangle][0] == n1)
               edgeIndex = 0;
             else if(triangle_list[oTriangle][1] == n1)
               edgeIndex = 1;
             else
               edgeIndex = 2;
             triToTri[oTriangle][edgeIndex] = iTriangle;
         }
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
}
/*bool DistPhysBAMIntersector::checkTriangulatedSurface() {
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
*/

void
DistPhysBAMIntersector::buildSolidNormals() {
  triNorms = new Vec3D[length_triangle_list];
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
//  std::cout <<"Going to make PhysBAMInterface" << std::endl;
  physInterface = new PhysBAMInterface<double>(physbam_triangulated_surface);
//  std::cout <<"Done making PhysBAMInterface" << std::endl;
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
  DistVec<double> weight(X.info());
  DistSVec<double,3> normApprox(X.info());
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new PhysBAMIntersector(*(d->getSubDomain()[i]), X(i), (*pseudoPhi)(i), *this);
    intersector[i]->computeLocalPseudoPhi(X(i), normApprox(i), weight(i));
  }
  // Addup the neighboring phis and weights
  d->assemble(*pseudoPhi);
  // Make sure the value for pseudoPhi is independent of summation order.
  pseudoPhi->zeroNonMaster();
  d->assemble(*pseudoPhi);

  d->assemble(weight);

  d->assemble(normApprox);
  // Use the values of phi and weights and finish the status
  for(int i = 0; i < numLocSub; ++i)
    intersector[i]->finishPseudoPhi(*(d->getSubDomain()[i]), X(i), normApprox(i), weight(i));
}


PhysBAMIntersector::PhysBAMIntersector(SubDomain &sub, SVec<double,3> &X,
              Vec<double> &pPhi, DistPhysBAMIntersector &distInt) :
 distIntersector(distInt), status(sub.numNodes()), phi(pPhi), locNorm(sub.numNodes()),
 edges(sub.getEdges()), edgeRes(sub.getEdges().size()), globIndex(sub.getGlobSubNum())
{
  int numEdges = edges.size();
  int (*ptr)[2] = edges.getPtr();

  for(int i = 0; i < numEdges; ++i) {
    edgeRes(i+1).x[1] = ptr[i][0]+1;
    edgeRes(i+1).x[2] = ptr[i][1]+1;
  }
  status = UNDECIDED;

  nodeMap = sub.getNodeMap();

  locToGlobNodeMap = sub.getNodeMap();
}

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

// check. (TODO:to be deleted.)
  double areaPAB = (0.5*(xA-xp)^(xB-xp))*dir;
  double xi3 = areaPAB/areaABC;
  if (xi1+xi2+xi3-1.0>1e-10) fprintf(stderr,"Oh no!\n");
}

int PhysBAMIntersector::closestTriangle(Vec3D x0, int tria, double& dist)
{
  const int MAX_ITER = 100;
  const double TOL = 1.0e-8;
  int iter = 0;
  int curTri = tria;
  double dst, xi1, xi2;
  pair<int,double> prevTri[MAX_ITER];

  while(iter<MAX_ITER) {
    projection(x0, curTri, xi1, xi2, dst);
    if (xi1>=-TOL && xi1<1.0+TOL && xi2>=-TOL && xi2<1.0+TOL) {//projection inside the triangle
      double temp_dst, temp_xi1, temp_xi2;
      int closestTri = curTri;
      for (int iNei=0; iNei<3; iNei++) {
        projection(x0, distIntersector.triToTri[curTri][iNei], temp_xi1, temp_xi2, temp_dst);
        if (temp_xi1>-TOL && temp_xi1<1.0+TOL && temp_xi2>-TOL && temp_xi2<1.0+TOL)
          if (fabs(temp_dst)<fabs(dst)) {
            dst = temp_dst;
            closestTri = distIntersector.triToTri[curTri][iNei];
          }
      }
      dist = dst;
      return closestTri;
    }

    //check if "curTri" has already been visited. If yes, we run into a closed loop.
    for (int iPrev=0; iPrev<iter; iPrev++)
      if (prevTri[iPrev].first==curTri) {
        //check. (TODO: to be deleted.)
        int closestTri = curTri;
        bool isPos = (dst+1.0e-16>0.0);
        for (int j=iPrev; j<iter; j++) {
          if ((prevTri[j].second+1.0e-16>0.0)!=isPos)
            {fprintf(stderr,"Cannot proceed. Stop.\n"); exit(-1);}
          if (fabs(prevTri[j].second+1.0e-16)<fabs(dst+1.0e-16)) {
            dst = prevTri[j].second;  closestTri = prevTri[j].first;
          }
        }

        dist = dst;
        return closestTri;
      }

    //put the last traversed triangle onto list.
    prevTri[iter].first = curTri;
    prevTri[iter].second = dst;

    //go to the correct neighbor.
    if(xi1<-TOL && (1.0-xi1-xi2)>=-2.0*TOL)  curTri = distIntersector.triToTri[curTri][1];
    else if(xi2<-TOL && (xi1>=-TOL))     curTri = distIntersector.triToTri[curTri][2];
    else                                 curTri = distIntersector.triToTri[curTri][0];

    iter++;
  }

  //shouldn't reach here if the closest triangle is found.
  fprintf(stderr,"failed in finding the closest triangle to (%e, %e, %e). Traversed triangles include:\n");
  for (int i=0; i<MAX_ITER-1; i++)
    fprintf(stderr,"TRIANGLE # %d (dist = %e).\n", prevTri[i].first, prevTri[i].second);
  exit(-1);
}

void PhysBAMIntersector::computeLocalPseudoPhi(SVec<double,3> &X, SVec<double,3> &normApprox,
                                               Vec<double> &weightSum)
{
  bool* edgeMasterFlag = edges.getMasterFlag();

  LIST_ARRAY<VECTOR<double,3> > xyz(X.size());
  for(int i = 0; i < X.size(); ++i)
    for(int j = 0; j < 3; ++j)
      xyz(i+1)[j+1] = X[i][j];

  double t0 = timer->getTime();
  distIntersector.getInterface().Intersect(xyz, edgeRes,distIntersector.getTolerance());
  double t = timer->getTime();
  nIntersect = 0;
  int numEdges = edges.size();
  int numNodes = status.size();

  phi = 0;
  normApprox = 0;
  weightSum = 0;

  LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > reverseEdgeRes(edges.size());
  for(int i = 0; i < edges.size(); ++i) {
      reverseEdgeRes(i+1).x[1] = edgeRes(i+1).x[2];
      reverseEdgeRes(i+1).x[2] = edgeRes(i+1).x[1];
  }
  distIntersector.getInterface().Intersect(xyz, reverseEdgeRes,distIntersector.getTolerance());

  // fix conflicts between edgeRes and reverseEdgeRes.

  for (int i=0; i< edges.size(); i++) {
    int trIDij = edgeRes(i+1).y.triangleID;
    int trIDji = reverseEdgeRes(i+1).y.triangleID;
    if (trIDij<0 && trIDji>=0) {
      edgeRes(i+1).y = reverseEdgeRes(i+1).y;
      edgeRes(i+1).y.alpha = 1.0 - edgeRes(i+1).y.alpha;
    }
    if (trIDji<0 && trIDij>=0) {
      reverseEdgeRes(i+1).y = edgeRes(i+1).y;
      reverseEdgeRes(i+1).y.alpha = 1.0 - reverseEdgeRes(i+1).y.alpha;
    }
/*
    if (trIDij>=0 && trIDji>=0) {
      if (trIDij==trIDji) {
        reverseEdgeRes(i+1).y = edgeRes(i+1).y;
        reverseEdgeRes(i+1).y.alpha = 1.0 - reverseEdgeRes(i+1).y.alpha;
      } else {
        if (edgeRes(i+1).y.alpha + reverseEdgeRes(i+1).y.alpha > 1.0) {
          IntersectionResult<double> temp = edgeRes(i+1).y;
          edgeRes(i+1).y = reverseEdgeRes(i+1).y;
          edgeRes(i+1).y.alpha = 1.0 - edgeRes(i+1).y.alpha;
          reverseEdgeRes(i+1).y = temp;
          reverseEdgeRes(i+1).y.alpha = 1.0 - reverseEdgeRes(i+1).y.alpha;
        }
      }
    }
*/
  }

  //Guanyuan's debug
  FILE* ftemp = fopen("temp","w");
  int (*ptr)[2] = edges.getPtr();
  for (int i=0; i<edges.size(); i++) {
    int nA = ptr[i][0], nB = ptr[i][1];
    if (nA==1092158 || nB==1092158)
      fprintf(ftemp,"%d %d %d %d\n", i+1, (int)1, nA+1, nB+1);
  }
  for (int i=0; i<edges.size(); i++) {
    int nA = ptr[i][0], nB = ptr[i][1];
    if (nA==1160679 || nB==1160679)
      fprintf(ftemp,"%d %d %d %d\n", i+1, (int)1, nA+1, nB+1);
  }
  for (int i=0; i<edges.size(); i++) {
    int nA = ptr[i][0], nB = ptr[i][1];
    if (nA==1083022 || nB==1083022)
      fprintf(ftemp,"%d %d %d %d\n", i+1, (int)1, nA+1, nB+1);
  }
  fclose(ftemp);


/*  int (*ptr)[2] = edges.getPtr();
  FILE* con1 = fopen("conflict1.top","w");
  FILE* con2 = fopen("conflict2.top","w");
  FILE* con3 = fopen("conflict3.top","w");
  for (int i=0; i<edges.size(); i++) {
    int trIDij = edgeRes(i+1).y.triangleID;
    int trIDji = reverseEdgeRes(i+1).y.triangleID;
    if (trIDij < 0 && trIDji <0) continue;
    if (trIDij==trIDji)
      if (edgeRes(i+1).y.alpha + reverseEdgeRes(i+1).y.alpha > 1+1e-10) {
        fprintf(stderr,"edge (%d, %d), trID = %d, alpha = %e, reverse alpha = %e.\n", ptr[i][0]+1, ptr[i][1]+1, trIDij, edgeRes(i+1).y.alpha, reverseEdgeRes(i+1).y.alpha);
        fprintf(con1, "%d %d %d %d\n", i+1, (int)1, ptr[i][0]+1, ptr[i][1]+1);
      }
    if (trIDij!=trIDji && trIDij>=0 && trIDji>=0)
      if (edgeRes(i+1).y.alpha + reverseEdgeRes(i+1).y.alpha > 1+1e-10) {
        fprintf(stderr,"edge (%d, %d), trIDij = %d, trIDji = %d, alpha = %e, reverse alpha = %e.\n", ptr[i][0]+1, ptr[i][1]+1, trIDij, trIDji, edgeRes(i+1).y.alpha, reverseEdgeRes(i+1).y.alpha);
        fprintf(con2, "%d %d %d %d\n", i+1, (int)1, ptr[i][0]+1, ptr[i][1]+1);
      }
    if (trIDij*trIDji<0) {
      fprintf(stderr,"edge (%d, %d), trIDij = %d, trIDji = %d, alpha = %e, reverse alpha = %e.\n", ptr[i][0]+1, ptr[i][1]+1, trIDij, trIDji, edgeRes(i+1).y.alpha, reverseEdgeRes(i+1).y.alpha);
      fprintf(con3, "%d %d %d %d\n", i+1, (int)1, ptr[i][0]+1, ptr[i][1]+1);
    }
  }
*/
  // Compute the contribution of each intersection to pseudoPhi, the normals and the weights.
  for(int i = 0; i < edgeRes.Size(); ++i)
    // check if this edge intersects the structure
    if(edgeRes(i+1).y.triangleID >= 0) {
      int p = edgeRes(i+1).x[1]-1, q = edgeRes(i+1).x[2]-1;
      if(edgeMasterFlag[i])
        updatePhi(p, q, edgeRes(i+1).y, X, phi, normApprox, weightSum);
      nIntersect++;
    }


  // Compute the contribution of each reverse intersection.
  for(int i = 0; i < edgeRes.Size(); ++i) {
    int p = reverseEdgeRes(i+1).x[1]-1, q = reverseEdgeRes(i+1).x[2]-1;
    secondIntersection[edges.find(p,q)] = reverseEdgeRes(i+1).y;
    if(reverseEdgeRes(i+1).y.triangleID >= 0) {
      if(edgeMasterFlag[edges.find(q,p)])
        updatePhi(p, q, reverseEdgeRes(i+1).y, X, phi, normApprox, weightSum);
      nIntersect++;
    }
  }

  t = timer->getTime();
  std::cout << "Number of intersections: " << nIntersect << " vs " << numEdges << " in " << (t-t0) << std::endl;
}

void PhysBAMIntersector::finishPseudoPhi(SubDomain &sub, SVec<double,3> &X, SVec<double,3> &normApprox,
                                         Vec<double> &weightSum)
 {
  int numNodes = status.size();
  int (*ptr)[2] = edges.getPtr();

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

  // find missed intersections.
  for (int i=0; i<numNodes; ++i) {
    if (weightSum[i] >0) { //look at its neighbors.
      Vec3D X0(X[i]);
      for (int iNei=0; iNei<nToN.num(i); ++iNei) {
        int me = nToN[i][iNei];
        if (weightSum[me]>0) continue;
        Vec3D X1(X[me]);
        double phiPrime = phi[i] + locNorm[i]*(X1-X0);
        if (phi[i]*phiPrime<0) { // we missed an intersection !
          phi[me] = -phi[i];
          locNorm[me] = locNorm[i];

          int edgeNum = edges.find(me, i);
          int n1 = ptr[edgeNum][0], n2 = ptr[edgeNum][1];
          if (edgeRes(edgeNum+1).y.triangleID>=0)
            fprintf(stderr,"Nothing is impossible... (especially for edge %d.)\n", edgeNum);

          Vec3D Xn1(X[n1]), Xn2(X[n2]);
          Vec3D dir(Xn2-Xn1);
          Xn1 = Xn1 - 0.1*dir;
          Xn2 = Xn2 + 0.1*dir;

          LIST_ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<double> > > myEdgeRes(1);
          myEdgeRes(1).x[1] = 1;  myEdgeRes(1).x[2] = 2;
          LIST_ARRAY<VECTOR<double,3> > myxyz(2);
          myxyz(1)[1] = Xn1[0];
          myxyz(1)[2] = Xn1[1];
          myxyz(1)[3] = Xn1[2];
          myxyz(2)[1] = Xn2[0];
          myxyz(2)[2] = Xn2[1];
          myxyz(2)[3] = Xn2[2];

          distIntersector.getInterface().Intersect(myxyz, myEdgeRes,distIntersector.getTolerance());
          if (myEdgeRes(1).y.triangleID<0) fprintf(stderr,"You have no hope... (node %d: phi = %e, locNorm = %e %e %e.  node %d.\n", i+1, phi[i], locNorm[i][0], locNorm[i][1], locNorm[i][2], me+1);
          else fprintf(stderr,"one missed intersection found.\n");
          edgeRes(edgeNum+1).y = myEdgeRes(1).y;

        }
      }
    }
  }

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

  // Kevin's debug
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
  exit(-1);

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
//  double dot = trNorm*edgeVec;
  double weight = std::abs(trNorm*edgeDir);

//  double phip = (alpha-1)*dot;
  int trNode = distIntersector.triangle_list[trID][0];
  Vec3D Xp(X[p][0],X[p][1],X[p][2]);
  double phip = trNorm*(Xp-distIntersector.solids_particle_list[trNode]);

  weight *= exp(-std::abs(phip)/edgeVec.norm());

  phi[p] += weight*phip;
  for(int j = 0; j < 3; ++j)
    normApprox[p][j] += weight*trNorm[j];
  weightSum[p] += weight;

}

LevelSetResult
PhysBAMIntersector::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
  int edgeNum = edges.find(ni, nj);
  int triangleID = edgeRes(edgeNum+1).y.triangleID-1;
  if(triangleID < 0) {
    Vec3D nrm = (std::abs(phi[ni]) < std::abs(phi[nj])) ? locNorm[ni] : locNorm[nj];
    if(nrm.norm() == 0)
      std::cerr << "Norm (1) is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
    return LevelSetResult(nrm[0], nrm[1], nrm[2], 0, 0, 0);
  }

  int triangle2ID = secondIntersection[edgeNum].triangleID-1;
  Vec3D nrm;
  IntersectionResult<double> result = edgeRes(edgeNum+1).y;
  if (triangle2ID>=0 && ni>nj)
    result = secondIntersection[edgeNum];
  int trueTriangleID = result.triangleID-1;
  nrm = distIntersector.getSurfaceNorm(trueTriangleID);

  if(nrm.norm() == 0)
    std::cerr << "Norm (2) is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;

  LevelSetResult lsRes(nrm[0], nrm[1], nrm[2], 0, 0, 0);
  lsRes.alpha = result.alpha;
  lsRes.xi[0] = result.zeta[0];
  lsRes.xi[1] = result.zeta[1];
  lsRes.xi[2] = 1-result.zeta[0]-result.zeta[1];
  lsRes.trNodes[0] = distIntersector.triangle_list[trueTriangleID][0];
  lsRes.trNodes[1] = distIntersector.triangle_list[trueTriangleID][1];
  lsRes.trNodes[2] = distIntersector.triangle_list[trueTriangleID][2];

  return lsRes;
}

bool PhysBAMIntersector::isActive(double t, int n) {
  return status[n] == INSIDE;
}

bool PhysBAMIntersector::edgeIntersectsStructure(double t, int ni, int nj) const {
  return status[ni] != status[nj];
}
