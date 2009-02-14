#ifndef _PHYSBAMINTERSECT_H_
#define _PHYSBAMINTERSECT_H_

#include <string>

#include "LevelSet/LevelSetStructure.h"
#include "LevelSet/IntersectionFactory.h"
#include "PHYSBAM_INTERFACE.h"

using PhysBAM::PhysBAMInterface;

class Vec3D;
class Communicator;
class PhysBAMIntersector;

class DistPhysBAMIntersector : public DistLevelSetStructure {
  protected:
    int length_solids_particle_list, length_triangle_list;
    int (*triangle_list)[3];
    Vec3D *solids_particle_list;
    Communicator *com;
    PhysBAMIntersector **intersector;

    PhysBAMInterface<double> *physInterface;
  public:
    DistPhysBAMIntersector();
    void init(std::string structureFileName);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    LevelSetStructure & operator()(int subNum) const;
};

class PhysBAMIntersector : public LevelSetStructure {


   public:
     PhysBAMIntersector();
     LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
     bool isActive(double t, int n);
     bool edgeIntersectsStructure(double t, int ni, int nj);



};

#endif
