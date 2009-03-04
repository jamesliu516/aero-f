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

static Timer *timer;
const int PhysBAMIntersector::UNDECIDED, PhysBAMIntersector::INSIDE, PhysBAMIntersector::OUTSIDE;

class PhysBAMIntersectorConstructor : public IntersectorConstructor {
    const char *structureFile;

  public:
    PhysBAMIntersectorConstructor() {
      structureFile = 0;
    }

    DistLevelSetStructure *getIntersector(IntersectProblemData&) {
      DistPhysBAMIntersector *inter = new DistPhysBAMIntersector();
      std::string solidSurface = structureFile;
      inter->init(solidSurface);
      return inter;
     // return 0;
    }

    int print();

    void init(ParseTree &dataTree) {
        std::cout << "inside the init function" << std::endl;
        ClassAssigner *ca = new ClassAssigner("PhysBAMIntersectorConstructor", 1, 0);
        new ClassStr<PhysBAMIntersectorConstructor>(ca, "structureFile", this, &PhysBAMIntersectorConstructor::structureFile);
        dataTree.implement(ca);
    }
};


IntersectorConstructor *myIntersect =
  IntersectionFactory::registerClass("PhysBAM", new PhysBAMIntersectorConstructor());


int PhysBAMIntersectorConstructor::print() {
  return 0;
}



DistPhysBAMIntersector::DistPhysBAMIntersector() {
  com = IntersectionFactory::getCommunicator();
}

LevelSetStructure &
DistPhysBAMIntersector::operator()(int subNum) const {
  return *intersector[subNum];
}

void DistPhysBAMIntersector::init(std::string solidSurface) {
  //2.read data from "prolateSurface.top".
  FILE *topFile;
  topFile = fopen(solidSurface.c_str(), "r");
  if (topFile == NULL) {fprintf(stderr, "topFile doesn't exist at all :(\n"); exit(1); }

  fscanf(topFile,"%d %d", &length_solids_particle_list, &length_triangle_list);
  triangle_list = new int[length_triangle_list][3];
  solids_particle_list = new Vec3D[length_solids_particle_list];

  com->fprintf(stderr,"solid surface: %d nodes, %d elements.\n", length_solids_particle_list, length_triangle_list);

  int thisNode;
  for (int iNode=0; iNode<length_solids_particle_list; iNode++)
    fscanf(topFile, "%d %lf %lf %lf", &thisNode, &(solids_particle_list[iNode][0]),
        &(solids_particle_list[iNode][1]), &(solids_particle_list[iNode][2]));
  if (thisNode!=length_solids_particle_list) {fprintf(stderr,"error in loading surface from file *!\n"); exit(1);}

  int nothing;
  for (int iElem=0; iElem<length_triangle_list; iElem++) {
    fscanf(topFile, "%d %d %d %d %d", &thisNode, &nothing, &(triangle_list[iElem][0]), &(triangle_list[iElem][1]),
        &(triangle_list[iElem][2]));
    triangle_list[iElem][0]--; triangle_list[iElem][1]--; triangle_list[iElem][2]--;
  }
  if (thisNode!=length_triangle_list) {fprintf(stderr,"error in loading surface from file **!\n", thisNode); exit(1);}
  fclose(topFile);

  //3. verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"Checking the solid surface...\n");
  if (checkTriangulatedSurface()) com->fprintf(stderr,"OK.\n");
  else exit(-1);

  initializePhysBAM();
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
    if(nrm != 0)
       triNorms[iTriangle] /= nrm;
  }
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
  buildSolidNormals();
}

void
DistPhysBAMIntersector::initialize(Domain *d, DistSVec<double,3> &X) {
  this->X = &X;
  domain = d;
  timer = domain->getTimer();
  numLocSub = d->getNumLocSub();
  intersector = new PhysBAMIntersector*[numLocSub];
  pseudoPhi = new DistVec<double>(X.info());
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new PhysBAMIntersector(*(d->getSubDomain()[i]), X(i), *this);
  }
}


PhysBAMIntersector::PhysBAMIntersector(SubDomain &sub, SVec<double,3> &X, DistPhysBAMIntersector &distInt) :
 distIntersector(distInt), status(sub.numNodes()), phi(sub.numNodes()), locNorm(sub.numNodes()),
 edges(sub.getEdges()), edgeRes(sub.getEdges().size())
{
  int numEdges = edges.size();
  int (*ptr)[2] = edges.getPtr();

  for(int i = 0; i < numEdges; ++i) {
    edgeRes(i+1).x[1] = ptr[i][0]+1;
    edgeRes(i+1).x[2] = ptr[i][1]+1;
  }
  LIST_ARRAY<VECTOR<double,3> > xyz(X.size());
  for(int i = 0; i < X.size(); ++i)
    for(int j = 0; j < 3; ++j)
      xyz(i+1)[j+1] = X[i][j];
  double t0 = timer->getTime();
  distIntersector.getInterface().Intersect(xyz, edgeRes,1e-3);
  double t = timer->getTime();
  int nIntersect = 0;

  status = UNDECIDED;

  Vec<bool> isVisited(sub.numNodes());
  Vec<double> weightSum(sub.numNodes());
  Vec<Vec3D> normApprox(sub.numNodes());
  phi = 0;
  weightSum = 0;
  isVisited = false;

  for(int i = 0; i < numEdges; ++i)
    if(edgeRes(i+1).y.triangleID >= 0) {
      int p = edgeRes(i+1).x[1]-1, q = edgeRes(i+1).x[2]-1;
      const Vec3D &trNorm = distIntersector.getSurfaceNorm(edgeRes(i+1).y.triangleID-1);
      Vec3D edgeVec(X[q][0]-X[p][0], X[q][1]-X[p][1], X[q][2]-X[p][2]);
      if(edgeVec.norm() == 0)
        continue;
      Vec3D edgeDir = edgeVec / edgeVec.norm();
      double weight = std::abs(edgeDir*trNorm)+1e-3;
      double dot = trNorm*edgeVec;
      double alpha = edgeRes(i+1).y.alpha;
      double phiq = alpha*dot;
      double phip = phiq-dot;

      phi[p] += weight*phip;
      phi[q] += weight*phiq;
      normApprox[p] += weight*trNorm;
      normApprox[q] += weight*trNorm;
      weightSum[p] += weight;
      weightSum[q] += weight;
      /*bool orientation = trNorm*edgeVec < 0;
      int pp = status[p];
      int pq = status[q];
      status[p] = orientation ? INSIDE : OUTSIDE;
      status[q] = orientation ? OUTSIDE : INSIDE;
      if(pp != UNDECIDED && pp != status[p])
        std::cout << "p Change of heart on " << p << " " << pp << status[p] << std::endl;
      if(pq != UNDECIDED && pq != status[q])
        std::cout << "q Change of heart on " << q << " " << pq << status[q] << std::endl;*/
      isVisited[p] = isVisited[q] = true;
      nIntersect++;
    }

  for(int i = 0; i < sub.numNodes(); ++i)
    if(weightSum[i] > 0) {
      phi[i] /= weightSum[i];
      status[i] = (phi[i] >= -1e-3) ? INSIDE : OUTSIDE;
      double nl = normApprox[i].norm();
      if(nl != 0)
        locNorm[i] = (normApprox[i] /= nl);
      else
        std::cout << "Really null norm" << std::endl;
    }

  Vec<int> list(sub.numNodes());
  Connectivity &nToN = *(sub.createEdgeBasedConnectivity());
  // Look for a start point
  int next = 0, lead = 0;
  for(int i = 0; i < sub.numNodes(); ++i)
    if(status[i] != UNDECIDED)
      list[lead++] = i;
  std::cout << "Initial lead: " << lead << std::endl;
  while(next < lead) {
    int cur = list[next++];
    int curStatus = status[cur];
    if(curStatus == UNDECIDED)
      std::cout << "Weird!!!" << std::endl;
    for(int i = 0; i < nToN.num(cur); ++i) {
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        list[lead++] = nToN[cur][i];
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
  std::cout << "Number of intersections: " << nIntersect << " vs " << numEdges << " in " << (t-t0) << std::endl;
  std::cout << "Check: " << lead << " vs "<< sub.numNodes() << std::endl;
  std::cout << "Inside: " << nInside << " outside: " << nOutside << " undecided: " << nUndecided << std::endl;
  // Supplemental check
  int numWeird = 0;
  for(int cur = 0; cur < sub.numNodes(); ++cur)
    for(int i = 0; i < nToN.num(cur); ++i)
      if(status[nToN[cur][i]] != status[cur])  {
        int edgeNum = edges.find(cur, nToN[cur][i]);
        if(edgeRes(edgeNum+1).y.triangleID < 0)
          numWeird++;
        if(!isVisited[cur] && !isVisited[nToN[cur][i]])
          std::cout << "Got an edge with no life!" << std::endl;
          //std::cout << "Found conflicting edge " << edgeNum << " between " << cur <<
          //   " and " << nToN[cur][i] << std::endl;
      }
  std::cout << "We have " << numWeird << " intersections to revisit" << std::endl;
}

LevelSetResult
PhysBAMIntersector::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
  int edgeNum = edges.find(ni, nj);
  int triangleID = edgeRes(edgeNum+1).y.triangleID-1;
 //  std:cerr << "Going for the norm " << triangleID << endl;
  if(triangleID < 0) {
    Vec3D nrm = (std::abs(phi[ni]) < std::abs(phi[nj])) ? locNorm[ni] : locNorm[nj];
    // std::cerr << "Norm is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
    if(nrm.norm() == 0)
      std::cerr << "Norm (1) is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
    return LevelSetResult(nrm[0], nrm[1], nrm[2], 0, 0, 0);
  }
  Vec3D nrm = distIntersector.getSurfaceNorm(triangleID);
  // std::cerr << "Norm is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
  if(nrm.norm() == 0)
    std::cerr << "Norm (2) is " << nrm[0] << " " << nrm[1] << " " << nrm[2] << std::endl;
  return LevelSetResult(nrm[0], nrm[1], nrm[2], 0, 0, 0);
}

bool PhysBAMIntersector::isActive(double t, int n) {
  return status[n] == INSIDE;
}

bool PhysBAMIntersector::edgeIntersectsStructure(double t, int ni, int nj) const {
  return status[ni] != status[nj];
}
