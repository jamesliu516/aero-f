#include <SurfMeshGen.h>

template<int dim>
SurfMeshGen<dim>::SurfMeshGen(Communicator *_com, IoData &_ioData, Domain &dom, DistGeoState *_geoState) : 
GappyOffline<dim>(_com, _ioData, dom, _geoState) {
	this->nSampleNodes = 1;
}

