#ifndef _GAPPY_ONLY_PSEUDO_H_
#define _GAPPY_ONLY_PSEUDO_H_

#include <GappyOffline.h>

template <int dim>
class GappyOfflineOnlyPseudo : public GappyOffline<dim> {

	typedef VecSet< DistSVec<double,dim> > SetOfVec;
	SetOfVec podHatTmp;

	int notRestricted;	// true if reading in full pod basis`

	virtual void setUpGreedy();

	virtual void readInPodResJac(int *);

	virtual void computePodTPod();

	virtual void determineSampleNodes();

	virtual void buildRemainingMesh();

	virtual void outputTopFile();

	virtual void outputReducedToFullNodes();

	virtual void outputSampleNodes();

	virtual void outputStateReduced();

	virtual void outputWallDistanceReduced();


public:

	GappyOfflineOnlyPseudo(Communicator *, IoData &, Domain &, DistGeoState *);

};
#include "GappyOfflineOnlyPseudo.C"
#endif
