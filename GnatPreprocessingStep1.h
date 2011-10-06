#ifndef _GAPPY_STEP_1_H_
#define _GAPPY_STEP_1_H_

#include <GnatPreprocessing.h>

template <int dim>
class GnatPreprocessingStep1 : public GnatPreprocessing<dim> {

	int computeGappyRes;
	virtual void computePseudoInverse();	// do nothing
	virtual void assembleOnlineMatrices();	// do nothing
	virtual void computePodTPod();
	virtual void outputOnlineMatricesGeneral(const char *onlineMatrix, 
			int numNodes, const std::map<int,int> &sampleNodeMap, const
			std::vector<int> &sampleNodeVec);

public:
	GnatPreprocessingStep1(Communicator *, IoData &, Domain &, DistGeoState *);
	//~GnatPreprocessingStep1();

};
#include "GnatPreprocessingStep1.C"
#endif
