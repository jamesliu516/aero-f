#include <GappyOfflineOnlyPseudo.h>
#include <Domain.h>

template<int dim>
GappyOfflineOnlyPseudo<dim>::GappyOfflineOnlyPseudo(Communicator *_com, IoData
		&_ioData, Domain &dom, DistGeoState *_geoState) : GappyOffline<dim>(_com,
			_ioData, dom, _geoState), podHatTmp(0, dom.getNodeDistInfo() ) { 

			if (this->ioData->gnat.sampleMeshUsed == GNATData::SAMPLE_MESH_NOT_USED)
				usingSampleMesh = false;
			else
				usingSampleMesh = true;

			if (usingSampleMesh) {
				this->outputOnlineMatricesFull = true;
				this->outputOnlineMatricesSample = false;
			}

}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::setUpGreedy() {

}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::readInPodResJac() {

	for (int i = 0 ; i < this->nPodBasis ; ++i){	// only do for number of required bases
		this->podHat[i].resize(this->nPod[i]);
	}

	//	read in both bases
	this->com->fprintf(stderr, " ... Reading POD bases for the residual and/or Jacobian ...\n");

	//this->domain.readMultiPodBasis(this->input->podFileResJacHat, this->podHat.a,
	//		this->nPod, this->nPodBasis, podFiles);

	//this->domain.readPodBasis(this->input->podFileResHat, this->nPod[0], this->podHat[0]);

	//if (this->nPodBasis == 2) {
	//	this->domain.readPodBasis(this->input->podFileJacHat, this->nPod[1], this->podHat[1]);
	//}
	this->domain.readPodBasis(this->input->podFileRes, this->nPod[0], this->podHat[0]);

	if (this->nPodBasis == 2) {
		this->domain.readPodBasis(this->input->podFileJac, this->nPod[1], this->podHat[1]);
	}

	if (!usingSampleMesh) {
		podHatTmp.resize(this->nPod[0]);
		for (int iPod = 0; iPod < this->nPod[0]; ++iPod) {
			podHatTmp[iPod] = 0.0;
		}
	}
}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::computePodTPod() {
	// read in PodTPod

	const char *onlineMatExtension = {".PodJacTPodRes"};
	FILE *onlineMatrix;
	int sp = strlen(this->ioData->output.transient.prefix);
	char *onlineMatrixFile = new char[sp +
		strlen(this->ioData->output.rom.onlineMatrix)+strlen(onlineMatExtension)+1];
	onlineMatrix = fopen(onlineMatrixFile, "r");

	if (this->thisCPU == 0) { 
		this->podTpod = new double * [this->nPod[1]];
		for (int i = 0; i < this->nPod[1]; ++i)
			this->podTpod[i] = new double [this->nPod[0]];
		for (int i = 0; i < this->nPod[1]; ++i)
			for (int j = 0; j < this->nPod[0]; ++j) { 
				fscanf(onlineMatrix, "%8.15e",&(this->podTpod[i][j]));
			}
	}
	delete [] onlineMatrixFile;
}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::determineSampleNodes() {

	// set globalSampleNodes, nSampleNodes

	this->nSampleNodes = 0;
	this->domain.readSampleNodes(this->globalSampleNodes, this->nSampleNodes,
			this->input->sampleNodes);
	for (int iSampleNode = 0; iSampleNode < this->nSampleNodes; ++iSampleNode) {
		int cpuTmp = 0;
		int subDTmp = 0;
		int locNodeTmp = 0;
		int globalSampleNode = this->globalSampleNodes[iSampleNode];

		this->globalSampleNodeRankMap.insert(pair<int, int > (globalSampleNode, iSampleNode));

		SubDomainData<dim> locPodHat, locPodHatTmp;
		for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
			int nLocNodes = this->nodeDistInfo.subSize(iSub);	// number of nodes in this subdomain
			int *locToGlobNodeMap = this->subD[iSub]->getNodeMap();
			bool *locMasterFlag = this->nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
			for (int iLocNode = 0; iLocNode < this->subD[iSub]->numNodes(); ++iLocNode) {	// all local globalNodes in subdomain
				if (locToGlobNodeMap[iLocNode] == globalSampleNode && locMasterFlag[iLocNode]) {
					cpuTmp = this->thisCPU;
					subDTmp = iSub;
					locNodeTmp = iLocNode;
					if (!usingSampleMesh) {
						assert(this->nPodBasis == 1);
						for (int iPod = 0; iPod < this->nPod[0]; ++iPod) {
							locPodHat = this->podHat[0][iPod].subData(iSub);
							locPodHatTmp = podHatTmp[iPod].subData(iSub);
							for (int iDim = 0; iDim < dim ; ++iDim) {
								locPodHatTmp[iLocNode][iDim] = locPodHat[iLocNode][iDim];	// zeros everywhere except at the chosen sample nodes
							}
						}
					}
					break;
				}
			}
		}
		this->com->barrier();
		this->com->globalSum(1,&cpuTmp);
		this->com->globalSum(1,&subDTmp);
		this->com->globalSum(1,&locNodeTmp);
		this->globalNodeToCpuMap.insert(pair<int, int > (globalSampleNode, cpuTmp));
		this->globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalSampleNode, subDTmp));
		this->globalNodeToLocalNodesMap.insert(pair<int, int > (globalSampleNode, locNodeTmp));
	}

	if (!usingSampleMesh) 
		this->podHat.a[0] = &podHatTmp;	// set podHatTmp to be the right one
}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::buildRemainingMesh() {
	// do nothing
}


template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputTopFile() {
	// do nothing
}

/*
template<int dim>
void GappyOfflineOnlyPseudo<dim>::assembleOnlineMatrices() {
	// do SAME
}
*/

/*
template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputOnlineMatrices() {
	// do SOMETHING
	// NEED numFullNodes, globalSampleNodeRankMap
}
*/

template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputReducedToFullNodes() {
	// do nothing
}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputSampleNodes() {
	// do nothing
}

template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputStateReduced() {
	// do nothing
}
template<int dim>
void GappyOfflineOnlyPseudo<dim>::outputWallDistanceReduced() {
	// do nothing
}
