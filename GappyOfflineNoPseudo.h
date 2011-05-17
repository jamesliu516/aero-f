#ifndef _GAPPY_NO_PSEUDO_H_
#define _GAPPY_NO_PSEUDO_H_

#include <GappyOffline.h>

template <int dim>
class GappyOfflineNoPseudo : public GappyOffline<dim> {

	virtual void computePseudoInverse();
	virtual void assembleOnlineMatrices();
	void computePodTPod();
	virtual void outputOnlineMatricesGeneral(const char *onlineMatrix, 
			int numNodes, const std::map<int,int> &sampleNodeMap, const
			std::vector<int> &sampleNodeVec);

public:
	GappyOfflineNoPseudo(Communicator *, IoData &, Domain &, DistGeoState *);
	//~GappyOfflineNoPseudo();

};
#include "GappyOfflineNoPseudo.C"
#endif
