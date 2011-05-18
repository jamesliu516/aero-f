#ifndef _SURF_MESH_GEN_H_
#define _SURF_MESH_GEN_H_

#include <GappyOffline.h>
template <int dim>
class SurfMeshGen : public GappyOffline<dim> {

	void setUpPodResJac() { ;}
	void setUpGreedy() { ;}
	void setUpPseudoInverse() { ;}

	void determineSampleNodes() { ;}
	void addSampleNodesAndNeighbors() { ;}
	void outputSampleNodes() { ;}
	void buildGappyMatrices() { ;}

	public:
	SurfMeshGen(Communicator *, IoData &, Domain &, DistGeoState *);
};
#include "SurfMeshGen.C"
#endif
