#include <GappyOfflineNoPseudo.h>

template<int dim>
GappyOfflineNoPseudo<dim>::GappyOfflineNoPseudo(Communicator *_com, IoData
		&_ioData, Domain &dom, DistGeoState *_geoState) : GappyOffline<dim>(_com,
			_ioData, dom, _geoState) {

	computeGappyRes = this->ioData->gnat.computeGappyRes;
}

template<int dim>
void GappyOfflineNoPseudo<dim>::computePseudoInverse() {

	// do not compute

}

template<int dim>
void GappyOfflineNoPseudo<dim>::assembleOnlineMatrices() {

	// do not compute

}

template<int dim>
void GappyOfflineNoPseudo<dim>::computePodTPod() {

	if (computeGappyRes==0)
		return;

	GappyOffline<dim>::computePodTPod();

	// output podTpod separately

	const char *onlineMatExtension = {".PodJacTPodRes"};
	FILE *onlineMatrix;

	int sp = strlen(this->ioData->output.transient.prefix);
	char *onlineMatrixFile = new char[sp +
		strlen(this->ioData->output.rom.onlineMatrix)+strlen(onlineMatExtension)+1];
	if (this->thisCPU ==0){
		sprintf(onlineMatrixFile, "%s%s%s",
				this->ioData->output.transient.prefix, this->ioData->output.rom.onlineMatrix,
				onlineMatExtension); 
		onlineMatrix = fopen(onlineMatrixFile, "wt");
	}

	for (int iPodJac = 0; iPodJac < this->nPod[1]; ++iPodJac) {
		for (int iPodRes = 0; iPodRes < this->nPod[0]; ++iPodRes) {
			this->com->fprintf(onlineMatrix,"%8.16e ", this->podTpod[iPodJac][iPodRes]);
		}
		this->com->fprintf(onlineMatrix,"\n");
	}

	for (int i = 0; i < this->nPod[1]; ++i)
		if (this->podTpod[i]) delete [] this->podTpod[i];
	if (this->podTpod) delete [] this->podTpod;

}

template<int dim>
void GappyOfflineNoPseudo<dim>::outputOnlineMatricesGeneral(const
		char *onlineMatricesName, int numNodes, const std::map<int,int>
		&sampleNodeMap, const std::vector<int> &sampleNodeVec) {

	// prepare files

	char *onlineMatrixFile;
	FILE *onlineMatrix;
	int sp = strlen(this->ioData->output.transient.prefix);
	const char *(onlineMatExtension [2]) = {".gappyResSample",".gappyJacSample"};

	int nPodBasisMax = this->nPodBasis;
	if (computeGappyRes==0)
		nPodBasisMax = 1;	// only output one of the matrices

	const char *onlineMatricesNameUse;
	const char *onlineMatricesNameUseExtension;

	for (int iPodBasis = 0; iPodBasis < nPodBasisMax; ++iPodBasis) {	
		GappyOffline<dim>::determineFileName(onlineMatricesName,
				onlineMatExtension[iPodBasis],onlineMatricesNameUse,onlineMatricesNameUseExtension);
		onlineMatrixFile = new char[sp +
			strlen(onlineMatricesNameUse)+strlen(onlineMatricesNameUseExtension)+1];
		if (this->thisCPU ==0){
			sprintf(onlineMatrixFile, "%s%s%s",
					this->ioData->output.rom.prefix, onlineMatricesNameUse,
					onlineMatricesNameUseExtension); 
			onlineMatrix = fopen(onlineMatrixFile, "wt");
		}

		this->com->fprintf(onlineMatrix,"Vector abMatrix under load for FluidNodes\n");
		this->com->fprintf(onlineMatrix,"%d\n", numNodes);
		outputReducedSVec(this->podHat[iPodBasis][0],onlineMatrix,this->nPod[iPodBasis]); // dummy output
		for (int iPod = 0; iPod < this->nPod[iPodBasis]; ++iPod) {	// # rows in A and B
			this->com->fprintf(stderr," ... writing vector %d of %d ... \n", iPod, this->nPod[iPodBasis]);
			outputReducedSVec(this->podHat[iPodBasis][iPod],onlineMatrix,iPod);
		}
		delete [] onlineMatrixFile;
	}
}
