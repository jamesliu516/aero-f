#ifndef _GAPPY_STEP_2_H_
#define _GAPPY_STEP_2_H_

#include <GnatPreprocessing.h>

template <int dim>
class GnatPreprocessingStep2 : public GnatPreprocessing<dim> {

	typedef VecSet< DistSVec<double,dim> > SetOfVec;
	SetOfVec podHatTmp;

	bool backupPlan;	// true if reading in full pod basis`

	virtual void setUpGreedy();

	virtual void readInPodResJac();

	virtual void computePodTPod();

	virtual void determineSampleNodes();

	virtual void buildRemainingMesh();

	virtual void outputTopFile();

	virtual void outputReducedToFullNodes();

	virtual void outputSampleNodes();

	virtual void outputStateReduced();

	virtual void outputWallDistanceReduced();


public:

	GnatPreprocessingStep2(Communicator *, IoData &, Domain &, DistGeoState *);

};
#include "GnatPreprocessingStep2.C"
#endif
