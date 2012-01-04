#include <GnatPreprocessingStep1.h>

template<int dim>
GnatPreprocessingStep1<dim>::GnatPreprocessingStep1(Communicator *_com, IoData
		&_ioData, Domain &dom, DistGeoState *_geoState) : GnatPreprocessing<dim>(_com,
			_ioData, dom, _geoState) {

	computeGappyRes = this->ioData->gnat.computeGappyRes;
}

template<int dim>
void GnatPreprocessingStep1<dim>::computePseudoInverse() {

	// do not compute

}

template<int dim>
void GnatPreprocessingStep1<dim>::assembleOnlineMatrices() {

	// do not compute

}

template<int dim>
void GnatPreprocessingStep1<dim>::computePodTPod() {

	if (computeGappyRes==0)
		return;

	GnatPreprocessing<dim>::computePodTPod();

	// output podTpod separately

	const char *podTpodExtension = {".PodJacTPodRes"};
	FILE *PODTPOD;

	int sp = strlen(this->ioData->output.rom.prefix);
	char *podTpodFile = new char[sp +
		strlen(this->ioData->output.rom.podNonlinRed)+strlen(podTpodExtension)+1];
	if (this->thisCPU ==0){
		sprintf(podTpodFile, "%s%s%s",	// TODO fix this
				this->ioData->output.rom.prefix, this->ioData->output.rom.podNonlinRed,
				podTpodExtension); 
		PODTPOD = fopen(podTpodFile, "wt");
	}

	for (int iPodJac = 0; iPodJac < this->nPod[1]; ++iPodJac) {
		for (int iPodRes = 0; iPodRes < this->nPod[0]; ++iPodRes) {
			this->com->fprintf(PODTPOD,"%8.16e ", this->podTpod[iPodJac][iPodRes]);
		}
		this->com->fprintf(PODTPOD,"\n");
	}

	for (int i = 0; i < this->nPod[1]; ++i)
		if (this->podTpod[i]) delete [] this->podTpod[i];
	if (this->podTpod) delete [] this->podTpod;
	if (podTpodFile) delete [] podTpodFile; 
	if (this->thisCPU ==0) fclose(PODTPOD);

}

template<int dim>
void GnatPreprocessingStep1<dim>::outputOnlineMatricesGeneral(const
		char *onlineMatricesName, int numNodes, const std::map<int,int>
		&sampleNodeMap, const std::vector<int> &sampleNodeVec) {

	// prepare files

	char *onlineMatrixFile;
	FILE *onlineMatrix;
	int sp = strlen(this->ioData->output.rom.prefix);
	const char *(onlineMatExtension [2]) = {".robResSample",".robJacSample"};

	int nPodBasisMax = this->nPodBasis;
	if (computeGappyRes==0)
		nPodBasisMax = 1;	// only output one of the matrices

	const char *onlineMatricesNameUse;
	const char *onlineMatricesNameUseExtension;

	for (int iPodBasis = 0; iPodBasis < nPodBasisMax; ++iPodBasis) {	
		GnatPreprocessing<dim>::determineFileName(this->ioData->output.rom.podNonlinRed,
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
		if (onlineMatrixFile) delete [] onlineMatrixFile;
		if (this->thisCPU ==0) fclose(onlineMatrix);
	}
	this->com->barrier();
}
