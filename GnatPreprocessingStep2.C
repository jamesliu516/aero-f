#include <GnatPreprocessingStep2.h>
#include <Domain.h>

template<int dim>
GnatPreprocessingStep2<dim>::GnatPreprocessingStep2(Communicator *_com, IoData
		&_ioData, Domain &dom, DistGeoState *_geoState) : GnatPreprocessing<dim>(_com,
			_ioData, dom, _geoState), podHatTmp(0, dom.getNodeDistInfo() ) { 

			if (this->ioData->gnat.sampleMeshUsed == GNATData::SAMPLE_MESH_NOT_USED)
				backupPlan = true;
			else
				backupPlan = false;

			this->outputOnlineMatricesFull = true;
			this->outputOnlineMatricesSample = false;	// already using sample mesh, so output in `full' coordinates

}

template<int dim>
void GnatPreprocessingStep2<dim>::setUpGreedy() {
// do nothing
}

template<int dim>
void GnatPreprocessingStep2<dim>::readInPodResJac() {

	for (int iPodBasis = 0; iPodBasis < this->nPodBasis; ++iPodBasis){
		this->podHat[iPodBasis].resize(this->nPod[iPodBasis]);
		for (int i = 0; i < this->nPod[iPodBasis]; ++i) this->podHat[iPodBasis][i] = 0.0;
	}

	//	read in both bases
	this->com->fprintf(stderr, " ... Reading POD bases for the residual and/or Jacobian ...\n");

	this->domain.readPodBasis(this->input->podFileRes, this->nPod[0], this->podHat[0]);

	if (this->nPodBasis == 2) {
		this->domain.readPodBasis(this->input->podFileJac, this->nPod[1], this->podHat[1]);
	}

	if (backupPlan) {
		podHatTmp.resize(this->nPod[0]);
		for (int iPod = 0; iPod < this->nPod[0]; ++iPod) {
			podHatTmp[iPod] = 0.0;
		}
	}
}

template<int dim>
void GnatPreprocessingStep2<dim>::computePodTPod() {
	// read in PodTPod

	const char *onlineMatExtension = {".PodJacTPodRes"};
	FILE *onlineMatrix;
	int sp = strlen(this->ioData->output.rom.prefix);
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
void GnatPreprocessingStep2<dim>::determineSampleNodes() {

	// set globalSampleNodes, nSampleNodes

	//int dbgWait = 0;
	//if (this->thisCPU == 1)
	//	dbgWait = 1;
	//while (dbgWait==1);
	this->com->barrier();
	this->nSampleNodes = 0;
	this->domain.readSampleNodes(this->globalSampleNodes, this->nSampleNodes,
			this->input->sampleNodes);
	this->com->barrier();
	for (int iSampleNode = 0; iSampleNode < this->nSampleNodes; ++iSampleNode) {
		int cpuTmp = 0;
		int subDTmp = 0;
		int locNodeTmp = 0;
		int globalSampleNode = this->globalSampleNodes[iSampleNode];

		this->globalSampleNodeRankMap.insert(pair<int, int > (globalSampleNode, iSampleNode));

		SubDomainData<dim> locPodHat, locPodHatTmp;
		for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
			bool foundNode = false;
			int nLocNodes = this->nodeDistInfo.subSize(iSub);	// number of nodes in this subdomain
			int *locToGlobNodeMap = this->subD[iSub]->getNodeMap();
			bool *locMasterFlag = this->nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
			for (int iLocNode = 0; iLocNode < this->subD[iSub]->numNodes(); ++iLocNode) {	// all local globalNodes in subdomain
				if (locToGlobNodeMap[iLocNode] == globalSampleNode && locMasterFlag[iLocNode]) {
					cpuTmp = this->thisCPU;
					subDTmp = iSub;
					locNodeTmp = iLocNode;
					if (backupPlan) {
						assert(this->nPodBasis == 1);
						for (int iPod = 0; iPod < this->nPod[0]; ++iPod) {
							locPodHat = this->podHat[0][iPod].subData(iSub);
							locPodHatTmp = podHatTmp[iPod].subData(iSub);
							for (int iDim = 0; iDim < dim ; ++iDim) {
								locPodHatTmp[iLocNode][iDim] = locPodHat[iLocNode][iDim];	// zeros everywhere except at the chosen sample nodes
							}
						}
					}
					foundNode = true;
					break;
				}
			}
			if (foundNode == true) {
				break;
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

	if (backupPlan) 
		this->podHat.a[0] = &podHatTmp;	// set podHatTmp to be the right one
}

template<int dim>
void GnatPreprocessingStep2<dim>::buildRemainingMesh() {
	// do nothing
}


template<int dim>
void GnatPreprocessingStep2<dim>::outputTopFile() {
	// do nothing
}

/*
template<int dim>
void GnatPreprocessingStep2<dim>::assembleOnlineMatrices() {
	// do same
}
*/

/*
template<int dim>
void GnatPreprocessingStep2<dim>::outputOnlineMatrices() {
	// do same
	
}
*/

template<int dim>
void GnatPreprocessingStep2<dim>::outputReducedToFullNodes() {
	// do nothing
}

template<int dim>
void GnatPreprocessingStep2<dim>::outputSampleNodes() {
	// do nothing
}

template<int dim>
void GnatPreprocessingStep2<dim>::outputStateReduced() {
	// do nothing
}
template<int dim>
void GnatPreprocessingStep2<dim>::outputWallDistanceReduced() {
	// do nothing
}
