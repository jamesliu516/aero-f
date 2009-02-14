#include <stdio.h>
#include <iostream>
#include "PhysBAMIntersect.h"
#include "Vector3D.h"
#include "Communicator.h"

#include "LevelSet/IntersectionFactory.h"
#include "parser/Assigner.h"



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


PhysBAMIntersector::PhysBAMIntersector() {

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

#include "PHYSBAM_INTERFACE.h"
#include "Particles/SOLIDS_PARTICLE.h"
#include "Grids/TRIANGLE_MESH.h"
#include "Geometry/TRIANGULATED_SURFACE.h"

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
}

LevelSetResult
PhysBAMIntersector::getLevelSetDataAtEdgeCenter(double t, int ni, int nj) {
 return LevelSetResult(0,0,0,1,0,0);
}

bool PhysBAMIntersector::isActive(double t, int n) {
  return true;
}

bool PhysBAMIntersector::edgeIntersectsStructure(double t, int ni, int nj) {
  return false;
}
