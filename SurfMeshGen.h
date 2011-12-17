#ifndef _SURF_MESH_GEN_H_
#define _SURF_MESH_GEN_H_

#include <GnatPreprocessing.h>
template <int dim>
class SurfMeshGen : public GnatPreprocessing<dim> {

	void setUp() { ;}

	void determineSampleNodes() { ;}

	void computePseudoInverse() { ;}

	void assembleOnlineMatrices() { ;}

	void outputOnlineMatrices() { ;}

	void outputSampleNodes() { ;}

	void addSampleNodesAndNeighbors() { ;}

	public:
	SurfMeshGen(Communicator *, IoData &, Domain &, DistGeoState *);
};
#include "SurfMeshGen.C"
#endif