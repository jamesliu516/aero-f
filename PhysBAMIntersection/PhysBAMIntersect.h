#ifndef _PHYSBAMINTERSECT_H_
#define _PHYSBAMINTERSECT_H_

#include <string>

#include "LevelSet/LevelSetStructure.h"
#include "LevelSet/IntersectionFactory.h"

class Vec3D;
class Communicator;

class PhysBAMIntersector : public LevelSetStructure {
     int length_solids_particle_list, length_triangle_list;
     int (*triangle_list)[3];
     Vec3D *solids_particle_list;
     Communicator *com;

   public:
     PhysBAMIntersector();
     LevelSetResult
       getLevelSetDataAtEdgeCenter(double t, int ni, int nj);
     bool isActive(double t, int n);
     bool edgeIntersectsStructure(double t, int ni, int nj);


     void init(std::string structureFileName);
     bool checkTriangulatedSurface();
     void initializePhysBAM();
};

#endif
