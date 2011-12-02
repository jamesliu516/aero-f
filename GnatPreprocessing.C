#include <GnatPreprocessing.h>
#include <TsDesc.h>

template<int dim>
GnatPreprocessing<dim>::GnatPreprocessing(Communicator *_com, IoData &_ioData, Domain
		&dom, DistGeoState *_geoState) : 
	domain(dom), 	com(_com), ioData(&_ioData), 
	podRes(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	podJac(0, dom.getNodeDistInfo() ),
	podHatRes(0, dom.getNodeDistInfo() ),
	podHatJac(0, dom.getNodeDistInfo() ),
	errorRes(0, dom.getNodeDistInfo() ),
	errorJac(0, dom.getNodeDistInfo() ),
	pseudoInvRhs(0, dom.getNodeDistInfo() ),
	handledNodes(0), nPodBasis(0),
	debugging(true),
	// distribution info
	numLocSub(dom.getNumLocSub()), nTotCpus(_com->size()), thisCPU(_com->cpuNum()),
	nodeDistInfo(dom.getNodeDistInfo()), subD(dom.getSubDomain()),parallelRom(2),
	geoState(_geoState), X(_geoState->getXn()),
	residual(0), jacobian(1), outputOnlineMatricesFull(false), outputOnlineMatricesSample(true),
	initializeLeastSquaresDone(false)
{
	// create temporary objects to build postOp, which is needed for surfaces

	twoLayers = ioData->gnat.layers == 2;
	geoSourceTmp = new GeoSource(*ioData);
	tsDescTmp = new TsDesc<dim>(*ioData, *geoSourceTmp, &domain);
  bcDataTmp = tsDescTmp->createBcData(*ioData);
  varFcnTmp = new VarFcn(*ioData);
  postOp = new PostOperator<dim>(*ioData, varFcnTmp, bcDataTmp, geoState, &domain);

	input = new TsInput(_ioData);
	includeLiftFaces = ioData->gnat.includeLiftFaces;

	// initialize vectors to point to the approprite bases

	for(int i=0; i<2; ++i) parallelRom[i] = new ParallelRom<dim>(dom,_com);
	pod.a[0] = &podRes;	// make pod point to res and jac
	pod.a[1] = &podJac;
	podHat.a[0] = &podHatRes;	// make pod point to res and jac
	podHat.a[1] = &podHatJac;
	error.a[0] = &errorRes;	// make pod point to res and jac
	error.a[1] = &errorJac;

	handledVectors[0] = 0;	// have not yet handled any vectors
	handledVectors[1] = 0;
	nRhs[0] = 0;
	nRhs[1] = 0;

	int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE)+1;	// there is a zero condition
	for (int i = 0; i < 2; ++i){
		for (int j = 0; j < 3; ++j)
			bcFaces[i][j] = new std::vector< int > [BC_CODE_EXTREME];
		bcFaceSurfID[i] = new std::vector< int > [BC_CODE_EXTREME];
	}

  numFullNodes = domain.getNumGlobNode();	// # globalNodes in full mesh

	com->fprintf(stderr," ... Number of full nodes in domain is %d ...\n",numFullNodes);
	
	nodesToHandle = NULL;
	globalNodes = NULL;
	cpus = NULL;
	locSubDomains = NULL;
	localNodes = NULL;
	totalNodesCommunicated = NULL;
	totalEleCommunicated = NULL;


	for (int i = 0; i < 3; ++i)
		nodesXYZ[i] = NULL;
	elements = NULL;
	for (int i = 0; i < 4; ++i)
		elemToNode[i] = NULL;

	for (int i = 0; i < 2; ++i) {
		podHatPseudoInv[i] = NULL;
		nRhsGreedy[i] = NULL;
		onlineMatrices[i] = NULL;
	}

	onlyInletOutletBC = true;	// first node should be on inlet/outlet
}

//----------------------------------------------

template<int dim>
GnatPreprocessing<dim>::~GnatPreprocessing() 
{
  if (input) delete input;
	if (postOp) delete postOp;
	if (varFcnTmp) delete varFcnTmp;
	if (bcDataTmp) delete bcDataTmp;
	if (tsDescTmp) delete tsDescTmp;
	if (geoSourceTmp) delete geoSourceTmp;


	if (globalNodes) delete [] globalNodes;


}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::buildReducedModel() {

	setUp();

	// compute the sample mesh used by Gappy POD

	//======================================
	// PURPOSE
	// 	build sample mesh by selecting sample globalNodes
	// INPUTS
	// 	nPod[0], nPod[1], nSampleNodes
	// 	full domain decomposition: pod[0], pod[1]
	// OUTPUTS
	// 	reduced domain decomposition: mesh, podHat[0], podHat[1],
	//======================================

	//======================================
	// NOTES
	// nRhsMax = number of POD vectors for each interpolation node that is
	// selected. It can be 1 to dim, where dim corresponds to interpolation.
	//
	// Fix the size of PhiHat using trick 1
	//
	// First compute a node, then compute the associated mask
	//
	//======================================

	//======================================
	// compute nRhsMax, nGreedyIt, nodesToHandle, possibly fix nSampleNodes
	//======================================

	determineSampleNodes();	// use greedy algorithm to determine sample nodes

	buildRemainingMesh();	// need two node layers b/c 2nd order flux

	computePseudoInverse();	// only requires sample nodes

	if (thisCPU == 0) {

		outputTopFile();

		assembleOnlineMatrices(); // handle online matrices so you can free up pseudo inverse matrix

	}

	outputOnlineMatrices();

	if (thisCPU == 0)
		outputSampleNodes();

	// from here on, only require:
	// geoState, ioData, nReducedNodes, domain
	// globalNodes, globalNodeToCpuMap, globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

	outputStateReduced();	// distributed info (parallel)

	outputWallDistanceReduced();	// distributed info (parallel)

	com->fprintf(stderr," ... finished with buildGappy ...\n");

	// STRATEGY
	// compute the mesh for the masked domain
	// put zeros in the pod[0], pod[1] on unsampled nodes
	// don't evaluate rhat and jphihat for unsampled nodes

	// TRICK: adding zero rows to matrices has no effect on the qr decomposition

} 

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::setUp() {

	// determine whether one or two pod bases are used

	setUpPodResJac();
	setUpGreedy();
	setUpPseudoInverse();

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::setUpPodResJac() {

	// use one basis if 1) same file name, or 2) podFileJac unspecified

	if (strcmp(input->podFileRes,input->podFileJac)==0 ||
			strcmp(input->podFileJac,"")==0) 
		nPodBasis = 1;
	else 
		nPodBasis = 2;
	
	nPod[0] = ioData->gnat.nRobRes; 
	if (nPod[0] == -1)
		nPod[0] = ioData->gnat.nRobNonlin;

	nPod[1] = ioData->gnat.nRobJac;
	if (nPod[1] == -1)
		nPod[1] = ioData->gnat.nRobNonlin;

	int podFiles [2] = {0, 1};	// which pod files should be read

	nPodMax = max(nPod[0],nPod[1]);	// compute maximum nPod

	if (nPodBasis == 1) {
		pod.a[1] = &podRes;
		podHat.a[1] = &podHatRes;	// make pod point to res and jac
		error.a[1] = &errorRes;	
		podFiles[1]= -1;	// do not read the file for the Jacobian
	}

	if (ioData->gnat.robGreedy == GNATData::UNSPECIFIED_GREEDY ||
			ioData->gnat.robGreedy == GNATData::BOTH_GREEDY) {
		errorBasis[0] = 0;
		errorBasis[1] = nPodBasis-1;
	}
	else if (ioData->gnat.robGreedy == GNATData::RESIDUAL_GREEDY) {
		errorBasis[0] = 0;
		errorBasis[1] = 0;
	}
	else if (ioData->gnat.robGreedy == GNATData::JACOBIAN_GREEDY) {
		errorBasis[0] = 1;
		errorBasis[1] = 1;
	}

	readInPodResJac();

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::setUpPseudoInverse() {

	// compute pod[0]^Tpod[1] (so you can delete these from memory sooner)
	if (nPodBasis == 2)
		computePodTPod();

	initializeLeastSquares();	// no least squares the first greedy it
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::setUpGreedy() {

	//==================================================================
	// PURPOSE: compute number of POD basis vectors (nRhsMax) and globalNodes and handled by each
	// iteration of the greedy algorithm
	// OUTPUTS
	// 	nRhsMax, nGreedyIt, nodesToHandle
	// STRATEGY:
	// nothing special happens when nSampleNodes < nPodGreedy < nSampleNodes * dim 
	// 	1) require nPodGreedy < nSampleNodes * dim to avoid underdetermined system
	// 	2) if nSampleNodes > nPodGreedy, need to treat more globalNodes per iteration
	// 	(nRhsMax is 1)
	//==================================================================
	
	nPodGreedy = ioData->gnat.nRobGreedy;
	if (nPodGreedy == 0 || nPodGreedy > nPodMax)
		nPodGreedy = nPodMax;
	nSampleNodes = ioData->gnat.nSampleNodes;

	if (nSampleNodes == 0) {
		int sampleNodeFactor = (int) ioData->gnat.sampleNodeFactor;
		if (sampleNodeFactor == -1)
			sampleNodeFactor = 2;
		nSampleNodes = static_cast<int>(ceil(double(nPodMax *
						sampleNodeFactor)/double(dim)));
			// this will give interpolation or the smallest possible least squares
	}

	if (nSampleNodes * dim < nPodGreedy) {	
		int nSampleNodesOld = nSampleNodes; 
		nSampleNodes = static_cast<int>(ceil(double(nPodGreedy)/double(dim))); 
		com->fprintf(stderr,"Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d",nSampleNodesOld,nSampleNodes);
	}

	nRhsMax = static_cast<int>(ceil(double(nPodGreedy)/double(nSampleNodes))); // nSampleNodes * nRhsMax >= max(nPod[0],nPod[1])

	// the following should always hold because of the above fix (safeguard)
	
	if (nRhsMax > dim) {
		com->fprintf(stderr,"Warning: nRhsMax > dim. More nodes should have been added.");
	}

	//==================================================================
	// 2) if nSampleNodes > nPodGreedy, need to treat more nodes per iteration
	// strategy: fill more nodes at the earlier iterations because POD basis vectors are optimally ordered
	//==================================================================

	nGreedyIt = min(nPodGreedy, nSampleNodes);	// number of greedy iterations (at most nPodGreedy; if nSampleNodes > nPodGreedy need to take care of more nodes per iteration)
	nodesToHandle = new int[nGreedyIt];	// number of nodes for each greedy iteration

	for (int iGreedyIt = 0; iGreedyIt < nGreedyIt; ++iGreedyIt)	{
		nodesToHandle[iGreedyIt] = (nSampleNodes * nRhsMax) / nPodGreedy;
		if (iGreedyIt < nSampleNodes % nPodGreedy && nRhsMax ==1)	// only in the dangerous case with nRhsMax = 1
			++nodesToHandle[iGreedyIt];
	}

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
		nRhsGreedy[iPodBasis] = new int [nGreedyIt];
		int nRhsMin = min(nPod[iPodBasis],nPodGreedy) / nGreedyIt;
		int nRhsExtra = min(nPod[iPodBasis],nPodGreedy) % nGreedyIt;
		for (int iGreedyIt = 0; iGreedyIt < nGreedyIt; ++iGreedyIt) {
			nRhsGreedy[iPodBasis][iGreedyIt] = nRhsMin;
			if (iGreedyIt < nRhsExtra) ++nRhsGreedy[iPodBasis][iGreedyIt];
		}
	}

	// initialize sampled pod basis and error vectors

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
		podHat[iPodBasis].resize(nPod[iPodBasis]);
		for (int i = 0; i < nPod[iPodBasis]; ++i) podHat[iPodBasis][i] = 0.0;
		error[iPodBasis].resize(nRhsMax);
		for (int i = 0; i < nRhsMax; ++i) error[iPodBasis][i] = pod[iPodBasis][i];	// for the first iteration, just pick out largest element
	}

	//===============================================
	// initialize the least squares problems
	//===============================================
	
	initializeLeastSquares();	// no least squares the first greedy it
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::findMaxAndFillPodHat(const double myMaxNorm, const int
		locSub, const int locNode, const int globalNode) {

	//===============================================
	// PURPOSE: fill locPodHat
	// INPUTS
	// 	Local: myMaxNorm, locSub, locNode, globalNode
	// OUTPUTS
	// 	Global: locPodHat for maximum entry, globally summed cpuSample, locSubSample,
	// 	locNodeSample, globalSampleNodes, xyz
	//===============================================
		
	double globalMaxNorm = myMaxNorm;
	com->barrier();
	com->globalMax(1, &globalMaxNorm);  // find the maximum value over all cpus

	// define variables to be summed globally

	int cpuTemp = 0;
	int locSubTemp = 0;
	int locNodeTemp = 0;
	int globalNodeTemp = 0;
	double xyz [3];
	for (int i=0; i<3; ++i)
		xyz[i]=0.0;

	// ensure only one cpu enters this loop

	int cpuHasMaxVal = 0;	// indicates if CPU has max value
	int cpuNumWithMaxVal = nTotCpus;
	if (myMaxNorm == globalMaxNorm) {
		cpuHasMaxVal = 1;
		cpuNumWithMaxVal = thisCPU;
	}
	com->globalSum(1, &cpuHasMaxVal);	// total CPUs with maximum value
	if (cpuHasMaxVal > 1) { 
		com->globalMin(1, &cpuNumWithMaxVal);	// take CPU with smallest number that has max val
	}

	if (thisCPU == cpuNumWithMaxVal) {  // if this CPU has the maximum value

		// save the global subdomain and local node indices (sum at the very end of
		// algorithm)

		cpuTemp = thisCPU;
		locSubTemp = locSub;
		locNodeTemp = locNode;
		globalNodeTemp = globalNode;
		assert(locSubTemp!=-1 && locNodeTemp !=-1 && globalNodeTemp != -1);
		computeXYZ(locSub, locNode, xyz);

		// fill out sampled matrices (all columns for the current rows)

		SubDomainData<dim> locPod, locPodHat;

		for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
			for (int iPod = 0 ; iPod < nPod[iPodBasis]; ++iPod) {
				locPod = pod[iPodBasis][iPod].subData(locSub);	// cannot access iDim entry
				locPodHat = podHat[iPodBasis][iPod].subData(locSub);
				for (int iDim = 0; iDim < dim ; ++iDim) {
					locPodHat[locNode][iDim] = locPod[locNode][iDim];
						// zeros everywhere except at sample nodes
				}
			}
		}
	}

	// make sure all cpus have the same copy
	
	com->barrier();
	com->globalSum(1, &cpuTemp);
	com->globalSum(1, &locSubTemp);
	com->globalSum(1, &locNodeTemp);
	com->globalSum(1, &globalNodeTemp);
	com->globalSum(3, xyz);

	if (debugging){
	 com->fprintf(stderr, "CPU %d has sample node: globalNode = %d, locNode = %d, locSub = %d  \n",cpuTemp,globalNodeTemp,locNodeTemp,locSubTemp);
	}

	// add information to all cpus copies

	cpuSample.push_back(cpuTemp);
	locSubSample.push_back(locSubTemp);
	locNodeSample.push_back(locNodeTemp);
	globalSampleNodes.push_back(globalNodeTemp);

	// define maps for SAMPLE nodes
	
	globalSampleNodeRankMap.insert(pair<int, int > (globalNodeTemp,
				handledNodes)); globalNodeToCpuMap.insert(pair<int, int >
				(globalNodeTemp, cpuTemp));
	globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalNodeTemp,
				locSubTemp)); globalNodeToLocalNodesMap.insert(pair<int, int >
				(globalNodeTemp, locNodeTemp)); StaticArray<double, 3> XYZ(xyz);
	nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalNodeTemp,
				XYZ));

	++handledNodes;

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::determineSampleNodes() {

	for (int greedyIt = 0; greedyIt < nGreedyIt; ++greedyIt)  {
		com->fprintf(stderr,"... greedy iteration %d ...\n", greedyIt);
		greedyIteration(greedyIt);
	}

	if (nodesToHandle) delete [] nodesToHandle;
	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		if (nRhsGreedy[iPodBasis]) delete [] nRhsGreedy[iPodBasis];
	}

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
		// pod no longer needed
		pod[iPodBasis].resize(0);
		error[iPodBasis].resize(0);
	}

	if (debugging){
		com->fprintf(stderr,"globalSampleNodes are:");
		for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes)
			com->fprintf(stderr,"%d ",globalSampleNodes[iSampleNodes]);
		com->fprintf(stderr,"\n");
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::greedyIteration(int greedyIt) {

	// Differences for 1st iteration compared with other greedy iterations:
	// 1) no least squares problem is solved (just take the maximum entry)
	// 2) look at the inlet face to ensure boundary condition is handled

	double myMaxNorm;
	int locSub, locNode, globalNode;  // temporary global subdomain and local node

	bool doLeastSquares = true;

	if (greedyIt == 0) {	
		doLeastSquares = false;	// don't do least squares if first iteration
	}

	// determine number of rhs for each
	// typically have nRhs = nRhsMax. exception: there are less than nRhsMax remaining vectors

	for (int iPodBasis = 0; iPodBasis  < nPodBasis; ++iPodBasis)	
		nRhs[iPodBasis] = nRhsGreedy[iPodBasis][greedyIt];
	assert(max(nRhs[0],nRhs[1]) > 0);		// must have at least one RHS

	if (doLeastSquares) {
		leastSquaresReconstruction();		// solve the least-squares reconstruction
	}

	for (int iFillNode = 0; iFillNode < nodesToHandle[greedyIt]; ++iFillNode) {	// fill up the appropriate number of nodes

		// initialize parameters
		myMaxNorm = 0.0;	// initial maximum is zero
		locSub = -1; locNode = -1; globalNode = -1;	// where the maximum is located

		// loop over nodes, and add to set if it is the maximum
		// subdomains -> nodes
		for (int iSub = 0; iSub < numLocSub; ++iSub) {

			// get subdomain info for all RHS

			getSubDomainError(iSub);

			// find maximum error on the subdomain
			
			subDFindMaxError(iSub, onlyInletOutletBC, myMaxNorm, locSub, locNode, globalNode);
			
		}

		// find global subdomain number and local node number for node with maximum norm

		findMaxAndFillPodHat(myMaxNorm, locSub, locNode, globalNode);


		if (onlyInletOutletBC == true)	// only add one 
			onlyInletOutletBC = false;
	}

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)	
		handledVectors[iPodBasis] += nRhs[iPodBasis];

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::initializeLeastSquares() {

	// initialize least squares problems
	// TODO: only allocate memory for required columns!
	if (initializeLeastSquaresDone == true) return;

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		parallelRom[iPodBasis]->parallelLSMultiRHSInit(podHat[iPodBasis], error[iPodBasis], nPodGreedy);
	}

	initializeLeastSquaresDone = true;
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::makeNodeMaxIfUnique(double nodeError, double
		&myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int
		&globalNode) {

	// PURPOSE: make the node the current maximum on this cpu if it hasn't been added already
	// INPUT
	// 	Local: nodeError, myMaxNorm (can change), iSub, locNodeNum (in subdomain node numbering
	// system)
	// 	Global: handledNodes, globalSampleNodes,
	// OUTPUT: 
	// 	Local: myMaxNorm (can change), locSub, locNode, globalNode
	
	// only do if the node could be the maximum

	if (nodeError >= myMaxNorm) {

		bool newNode = true;
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();
		int thisGlobalNode = locToGlobNodeMap[locNodeNum];

		// check that the node hasn't already been added (look at globalSampleNodes)

		for (int iIslandCheck = 0; iIslandCheck < handledNodes; ++iIslandCheck) {
			if (thisGlobalNode == globalSampleNodes[iIslandCheck]) { 
				newNode = false;
				break;
			}
		}

		// make local maximum if it is maximum and not already in the set

		if (newNode) {
			myMaxNorm = nodeError;
			locSub = iSub;
			locNode = locNodeNum;
			globalNode = thisGlobalNode;
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError) {

	// PURPOSE: compute the sum of squares of node error for all RHS, both bases
	// 	at the locNodeNum node
	// INPUT
	// 	Local: locMasterFlag, locNodeNum 
	// 	Global: nRhs
	// OUTPUT
	// 	Global: nodeError

	 nodeError = 0.0; // initialize normed error to zero

	 if (locMasterFlag[locNodeNum]) { 	// use subdomain node number
		 for (int iPodBasis = errorBasis[0]; iPodBasis  <= errorBasis[1] ; ++iPodBasis)
			 for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs){  // add components of error from all vectors on this node
				 for (int k = 0; k < dim ; ++k){	// add contributions from residual and jacobian reconstruction errors where possible
					 double componentError = locError[iPodBasis][iRhs][locNodeNum][k];
					 nodeError += componentError * componentError;
			 }
		 }
	 }
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::getSubDomainError(int iSub) {

	// INPUT 
	// 	Passed: iSub
	// 	Global: nPodBasis, error
	// OUTPUT
	// 	Global: locError

	locError.resize(nRhs);	// create locError entity of maximal size nRhsMax

	for (int iPodBasis = 0; iPodBasis < nPodBasis ; ++iPodBasis) {
		for (int iRhs = 0; iRhs < nRhs[iPodBasis] ; ++iRhs){  // add components of error from all vectors on this node
			locError[iPodBasis][iRhs] = error[iPodBasis][iRhs].subData(iSub);	// first iteration, it is just the pod vectors themselves
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::leastSquaresReconstruction() {

	// PURPOSE: compute least squares reconstruction error of the pod basis
	// vectors
	// INPUT
	// 	Global: nPodBasis, nRhs, error, podHat, error

	double ** (lsCoeff[2]);
	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		lsCoeff[iPodBasis] = new double * [ nRhs[iPodBasis] ];
		for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs)  {
			// temporarily fill the error vector with the RHS (solving nRhs[iPodBasis] problems)
			error[iPodBasis][iRhs] = podHat[iPodBasis][handledVectors[iPodBasis] + iRhs];	// NOTE: PODHAT
			lsCoeff[iPodBasis][iRhs] = new double [handledVectors[iPodBasis]];
		}

		parallelLSMultiRHSGap(iPodBasis,lsCoeff[iPodBasis]);

		// compute reconstruction error

		for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs) { 
			error[iPodBasis][iRhs] = pod[iPodBasis][handledVectors[iPodBasis] + iRhs];	// NOTE: POD
			for (int jPod = 0; jPod < handledVectors[iPodBasis]; ++jPod) {
				error[iPodBasis][iRhs] -= pod[iPodBasis][jPod] * lsCoeff[iPodBasis][iRhs][jPod];
			}
		}
		for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs)  
			if (lsCoeff[iPodBasis][iRhs]) delete [] lsCoeff[iPodBasis][iRhs];
		if (lsCoeff[iPodBasis]) delete [] lsCoeff[iPodBasis];

	}
}

//----------------------------------------------

template<int dim> void GnatPreprocessing<dim>::subDFindMaxError(int iSub, bool
		onlyInletOutletBC, double &myMaxNorm, int &locSub, int &locNode, int
		&globalNode) {

	// PURPOSE: Search the iSub subdomain for possible maximum error
	// Inputs:
	// 	Passed: iSub, onlyInletOutletBC, myMaxNorm, locSub, locNode, globalNode
	// Outputs:
	// 	Passed: (all can change) myMaxNorm, locSub, locNode, globalNode

	bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
	int nLocNodes = nodeDistInfo.subSize(iSub);	// number of nodes in this subdomain
	double nodeError; // the error at the node

	if (onlyInletOutletBC) {	// consider only nodes on inlet BC
		FaceSet& currentFaces = subD[iSub]->getFaces();
		for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {	
			// only consider inlet boundary conditions
			if (currentFaces[iFace].getCode() == BC_INLET_MOVING ||
					currentFaces[iFace].getCode() == BC_INLET_FIXED ||
					currentFaces[iFace].getCode() == BC_OUTLET_MOVING ||
					currentFaces[iFace].getCode() == BC_OUTLET_FIXED) {
			 locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
			 nLocNodes = currentFaces[iFace].numNodes(); // number of nodes on this face
			 for (int iLocNode = 0; iLocNode < nLocNodes ; ++iLocNode){
				 int locNodeNum = currentFaces[iFace][iLocNode];	// local node number (in subdomain numbering) from the face
				 computeNodeError(locMasterFlag, locNodeNum, nodeError);	// compute the error at the node
				 // make the node max if it is the maximum so far and has not yet been added
				 makeNodeMaxIfUnique(nodeError, myMaxNorm, iSub, locNodeNum, locSub, locNode, globalNode);
			 }
			}
		}
	}

	else {	// consider all nodes
		for (int locNodeNum = 0; locNodeNum < nLocNodes ; ++locNodeNum){
			computeNodeError(locMasterFlag, locNodeNum, nodeError);	
			// make the node max if it is the maximum so far and has not yet been added
			makeNodeMaxIfUnique(nodeError, myMaxNorm, iSub, locNodeNum, locSub, locNode, globalNode);
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff) {

	bool lsCoeffAllCPU = true; // all cpus need solution
	parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis],error[iPodBasis],
			handledVectors[iPodBasis], nRhs[iPodBasis], lsCoeff, lsCoeffAllCPU);
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::buildRemainingMesh() {

	// globalNodes[iIsland][0] is the sample node itself 
	// globalNodes[iIsland][iNode] is the iNode neighbor of iIsland.

	// quantities for entire sample mesh (not just sample nodes)

	globalNodes = new std::vector <int> [nSampleNodes];
	cpus = new std::vector <int> [nSampleNodes];
	locSubDomains = new std::vector <int> [nSampleNodes];
	localNodes = new std::vector <int> [nSampleNodes];
	for (int i = 0; i < 3; ++i)
		nodesXYZ[i] = new std::vector <double> [nSampleNodes];
	elements = new std::vector <int> [nSampleNodes];
	for (int i = 0; i < 4; ++i)
		elemToNode[i] = new std::vector <int> [nSampleNodes];
	totalNodesCommunicated = new int [nSampleNodes];
	totalEleCommunicated = new int [nSampleNodes];
	for (int i = 0; i < nSampleNodes; ++i) {
		totalNodesCommunicated[i] = 0;
		totalEleCommunicated[i] = 0;
	}

	// compute two layers of globalNodes
	
	nodeOffset = new int [nSampleNodes];
	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) 
		nodeOffset[iSampleNodes] = 0;

	com->fprintf(stderr," ... adding sample nodes and neighbors ...\n");
	addSampleNodesAndNeighbors();	

	com->fprintf(stderr," ... computing BC faces ...\n");
	computeBCFaces(true);	// compute BC faces and add nodes/elements for lift surfaces

	// add all neighbor globalNodes and elements of the sample node's neighbors
	// and neighbor nodes of the lift surface nodes (need 1 layer of nodes for lift)

	communicateAll();	// all cpus need all nodes/elements for adding neighbors

	if (twoLayers) {
		com->fprintf(stderr," ... adding second layer of neighbors ...\n");
		for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
			com->fprintf(stderr," ... adding neighbors for sample node %d of %d...\n",iSampleNodes,nSampleNodes);
			addNeighbors(iSampleNodes, nodeOffset[iSampleNodes]);
		}
		communicateAll();	// all cpus need all nodes/elements to define maps
	}

	com->fprintf(stderr," ... defining maps ...\n");
	defineMaps();	// define element and node maps
		// deletes cpus, locSubDomains, locNodes, totalNodesCommunicated,
		// totalEleCommunicated, nodesXYZ, elemToNode
		// (want to keep globalNodes and elements: see below)

	// from now on, use maps and not these vectors

	cpuSample.erase(cpuSample.begin(),cpuSample.end());
	locSubSample.erase(locSubSample.begin(),locSubSample.end());
	locNodeSample.erase(locNodeSample.begin(),locNodeSample.end());

	// remove redundant entries from the globalNodes and elements

	domain.makeUnique(globalNodes, nSampleNodes);
	domain.makeUnique(elements, nSampleNodes);

	computeBCFaces(false);	// compute BC faces already in mesh
	communicateBCFaces();

	nReducedNodes = globalNodes[0].size();	// number of nodes in the sample mesh

	if (nodeOffset) delete [] nodeOffset; 

} 

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::addSampleNodesAndNeighbors() {

	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
		com->fprintf(stderr," ... adding neighbors for sample node %d of %d...\n",iSampleNodes,nSampleNodes);

		// add the sample node itself to all processors

		globalNodes[iSampleNodes].push_back(globalSampleNodes[iSampleNodes]);
		cpus[iSampleNodes].push_back(cpuSample[iSampleNodes]);
		locSubDomains[iSampleNodes].push_back(locSubSample[iSampleNodes]);
		localNodes[iSampleNodes].push_back(locNodeSample[iSampleNodes]);

		StaticArray<double, 3> xyzVals = nodesXYZmap.find(globalSampleNodes[iSampleNodes])->second;
		for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
			nodesXYZ[iXYZ][iSampleNodes].push_back(xyzVals[iXYZ]);
		}

		// add all neighbor globalNodes and elements of the sample node itself
		nodeOffset[iSampleNodes] = globalNodes[iSampleNodes].size();
		addNeighbors(iSampleNodes);
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::addNeighbors(int iIslands, int startingNodeWithNeigh = 0) {

	// add all global neighbor globalNodes/elements in the iIslands row of and elements to the iIsland node set
	
	Connectivity *nodeToNode, *nodeToEle, *eleToNode;
	double xyz [3] = {0.0, 0.0, 0.0};
	int globalNodeNum;	// NOTE: must work with global node numbers because the greedy selection only operated on master globalNodes. Here, we care about the node even if it wasn't a master node.
	int locNodeNum, locEleNum;	// local node number of the node/element to be added
	int nNodesToAdd, nEleToAdd;	// number of globalNodes that should be added
	bool elemBasedConnectivityCreated;	// whether createElemBasedConnectivity has been created for the current subdomain
	int *nToNpointer, *nToEpointer;	// pointers to connectivity information

	int nNodeWithNeigh = globalNodes[iIslands].size();	// number of globalNodes whose neighbors will be added
	bool *locMasterFlag;

	for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {	// all subdomains
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();
		int *locToGlobElemMap = subD[iSub]->getElemMap();
		locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain

		elemBasedConnectivityCreated = false;
		for (int iLocNode = 0; iLocNode < subD[iSub]->numNodes(); ++iLocNode) {	// all local globalNodes in subdomain
			//only add master nodes
			globalNodeNum = locToGlobNodeMap[iLocNode];	// global node number
			for (int iNodeWithNeigh = startingNodeWithNeigh; iNodeWithNeigh < nNodeWithNeigh; ++iNodeWithNeigh) {	// check if this local node is in the current row
				if (globalNodeNum == globalNodes[iIslands][iNodeWithNeigh]) {// add all neighbors of this node
					if (!elemBasedConnectivityCreated) {
						// taken from createElemBasedConnectivity
						nodeToNode = subD[iSub]->createElemBasedConnectivity();
						nodeToEle = subD[iSub]->createNodeToElementConnectivity();
						eleToNode = subD[iSub]->createElementToNodeConnectivity();
						elemBasedConnectivityCreated = true;
					}

					// add neighbors of iLocNode

					nToNpointer = nodeToNode->ptr();
					nToEpointer = nodeToEle->ptr();
					nNodesToAdd = nToNpointer[iLocNode+1] - nToNpointer[iLocNode];	// number of neighbors this node has
					nEleToAdd = nToEpointer[iLocNode+1] - nToEpointer[iLocNode];	// number of neighbors this node has

					for (int iNodeToAdd = 0; iNodeToAdd < nNodesToAdd; ++iNodeToAdd) { 
						locNodeNum = *((*nodeToNode)[iLocNode]+iNodeToAdd);
						globalNodes[iIslands].push_back(locToGlobNodeMap[locNodeNum]);
						cpus[iIslands].push_back(thisCPU);
						locSubDomains[iIslands].push_back(iSub);
						localNodes[iIslands].push_back(locNodeNum);
						for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) xyz[iXYZ] = 0.0; // KTCREMOVE
						computeXYZ(iSub, locNodeNum, xyz);
						for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
							nodesXYZ[iXYZ][iIslands].push_back(xyz[iXYZ]);
						}
					}
					for (int iEleToAdd = 0; iEleToAdd < nEleToAdd; ++iEleToAdd) { 
						// add global element
						locEleNum = *((*nodeToEle)[iLocNode]+iEleToAdd);
						elements[iIslands].push_back(locToGlobElemMap[locEleNum]);
						
						// determine globalNodes connected to the element
						for (int iNodesConn = 0; iNodesConn < 4; ++iNodesConn){
							int locNodeNumTmp = *((*eleToNode)[locEleNum] + iNodesConn);
							int globalNodeNumTmp = locToGlobNodeMap[locNodeNumTmp];
							elemToNode[iNodesConn][iIslands].push_back(globalNodeNumTmp);
						}
					}
				}
			}
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computeBCFaces(bool liftContribution) {

	// PURPOSE: determine which faces of the sample mesh are on the boundary
	// INPUT: liftContribution (true: adding faces contributing to lift; false:
	//   adding other BC faces)
	// ASSUME: call with liftContribution = true first!
	// METHOD: loop over ALL FACES
	// NOTE: this is done in parallel; thus, a communication is needed to make
	// sure all cpus have the same copy
		// determine if any faces are boundary conditions

	int faceBCCode = 0;
	bool includeFace;
	int globalFaceNodes [3];	// global node numbers of current face
	for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {	// all subdomains
		FaceSet& currentFaces = subD[iSub]->getFaces();	// faces on subdomain
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();	// global node numbering
		for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {	// check all faces	
			faceBCCode = currentFaces[iFace].getCode();
			int codeIsPos = faceBCCode >0;
			if (faceBCCode != 0) {
				
				// check if the face is in the sample mesh
				 if (liftContribution && includeLiftFaces > 0) {	// including lift faces
					 includeFace = checkFaceContributesToLift(currentFaces,
							 iFace, iSub, locToGlobNodeMap);
					 if (includeFace) {
						 addFaceNodesElements(currentFaces, iFace, iSub, locToGlobNodeMap);// add nodes
					 }
				 }
				 else {// include faces already in the mesh
					 includeFace = checkFaceInMesh(currentFaces, iFace, iSub, locToGlobNodeMap);
					 if (includeFace)
						 includeFace *= checkFaceAlreadyAdded(thisCPU, iSub, iFace);
				 }

				if (includeFace) { // add face to extra set
					for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) { 
						int localFaceNode = currentFaces[iFace][iFaceNode];
						globalFaceNodes[iFaceNode] = locToGlobNodeMap[localFaceNode];
					}
					for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) {
						bcFaces[codeIsPos][iFaceNode][abs(faceBCCode)].push_back(globalFaceNodes[iFaceNode]);
					}
					int surfaceID = currentFaces[iFace].getSurfaceID();
					bcFaceSurfID[codeIsPos][abs(faceBCCode)].push_back(surfaceID);
					if (liftContribution)	{// only do when considering lift surfaces
						// NOTE: each face is on one subdomain, so no communication of
						// bcFacesInfo is required
						int faceLoc [3]	= {thisCPU, iSub, iFace};
						StaticArray <int,3> faceLocation(faceLoc);
						bcFacesInfo.insert(faceLocation);
					}
				}
			}
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::communicateAll() {

	domain.communicateMesh(globalNodes, nSampleNodes, totalNodesCommunicated);
	domain.communicateMesh(cpus, nSampleNodes, totalNodesCommunicated);
	domain.communicateMesh(locSubDomains, nSampleNodes, totalNodesCommunicated);
	domain.communicateMesh(localNodes, nSampleNodes, totalNodesCommunicated);
	for (int i = 0; i < 3; ++i)
		domain.communicateMesh(nodesXYZ[i], nSampleNodes, totalNodesCommunicated);

	domain.communicateMesh(elements, nSampleNodes, totalEleCommunicated);
	for (int i = 0; i < 4; ++i)
		domain.communicateMesh(elemToNode[i], nSampleNodes, totalEleCommunicated);

	for (int i = 0; i < nSampleNodes; ++i) { 
		totalNodesCommunicated[i] = globalNodes[i].size();
		totalEleCommunicated[i] = elements[i].size();
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::defineMaps() {

	// defines nodesXYZmap and elemToNodeMap

	int globalNodeNumTmp;
	StaticArray<double, 3> nodesXYZTmp;
	int globalEleNumTmp;
	StaticArray<int, 4> elemToNodeTmp;
	std::map<int,int>::const_iterator sampleNodeRank; 

	// first time, establish that no maps have been defined

	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {

		// define nodesXYZmap
		for (int iNeighbor = 0; iNeighbor < globalNodes[iSampleNodes].size(); ++iNeighbor) {

			// do not re-define map for sample nodes (already defined as master) 
			globalNodeNumTmp = globalNodes[iSampleNodes][iNeighbor];
			sampleNodeRank = globalSampleNodeRankMap.find(globalNodeNumTmp);
			if (sampleNodeRank != globalSampleNodeRankMap.end()) {continue;}	// skip if a sample node

			int cpuTemp = cpus[iSampleNodes][iNeighbor];
			int locSubDomTemp = locSubDomains[iSampleNodes][iNeighbor];
			int localNodesTemp = localNodes[iSampleNodes][iNeighbor];
			for (int iXYZ = 0 ; iXYZ < 3; ++iXYZ)
				nodesXYZTmp[iXYZ] = nodesXYZ[iXYZ][iSampleNodes][iNeighbor];

			globalNodeToCpuMap.insert(pair<int, int > (globalNodeNumTmp, cpuTemp));
			globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalNodeNumTmp, locSubDomTemp));
			globalNodeToLocalNodesMap.insert(pair<int, int > (globalNodeNumTmp, localNodesTemp));
			nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalNodeNumTmp, nodesXYZTmp));
		}

		// define elemToNodeMap

		for (int iEle = 0; iEle < elements[iSampleNodes].size(); ++iEle) {
			globalEleNumTmp = elements[iSampleNodes][iEle];
			for (int iNodesConn = 0 ; iNodesConn  < 4; ++iNodesConn)
				elemToNodeTmp[iNodesConn] = elemToNode[iNodesConn][iSampleNodes][iEle];
			elemToNodeMap.insert(pair<int, StaticArray <int, 4> > (globalEleNumTmp, elemToNodeTmp));
		}
	}

	// no longer need this information (use in the form of maps)

	if (cpus) delete [] cpus;
	if (locSubDomains) delete [] locSubDomains;
	if (localNodes) delete [] localNodes;
	if (totalNodesCommunicated) delete [] totalNodesCommunicated;
	if (totalEleCommunicated) delete [] totalEleCommunicated;
for (int i = 0; i < 3; ++i)
		if (nodesXYZ) delete [] nodesXYZ[i];
	for (int i = 0; i < 4; ++i)
		if (elemToNode[i]) delete [] elemToNode[i];

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::communicateBCFaces(){
	
	int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE)+1;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j) {
			domain.communicateMesh(bcFaces[i][j], BC_CODE_EXTREME, NULL);
		}
		domain.communicateMesh(bcFaceSurfID[i], BC_CODE_EXTREME, NULL);
	}
}

//----------------------------------------------

template<int dim>
bool GnatPreprocessing<dim>::checkFaceInMesh(FaceSet& currentFaces, const int iFace, const int iSub, const int *locToGlobNodeMap){

	// PURPOSE: determine wheteher or not currentFace is in the sample mesh
	// OUTPUT: faceInMesh = true if the face is in the sample mesh

	bool faceInMesh = false;
	bool faceSomewhereInMesh;
	bool *nodeSomewhereInMesh = new bool [currentFaces[iFace].numNodes()];	// one bool for each face node
	int * globalNodeNum = new int [currentFaces[iFace].numNodes()];
	for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
		nodeSomewhereInMesh[iNodeFace] = false;
	}
	for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
		globalNodeNum[iNodeFace] = locToGlobNodeMap[currentFaces[iFace][iNodeFace]];
	}
	for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
		// check to see if the iNodeFace is in iIsland
		for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
			for (int iReducedMeshNode = 0; iReducedMeshNode < globalNodes[iSampleNodes].size(); ++iReducedMeshNode) { // globalNodes on island
				if (globalNodeNum[iNodeFace] == globalNodes[iSampleNodes][iReducedMeshNode]){ 
					nodeSomewhereInMesh[iNodeFace] = true;
					break;
				}
			}
			if (nodeSomewhereInMesh[iNodeFace]) break;
		}
	}

	faceSomewhereInMesh = true;
	for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace)
		faceSomewhereInMesh *= nodeSomewhereInMesh[iNodeFace];

	//check if a given element has all the globalNodes

	bool faceInElement = true;
	StaticArray <int, 4> globalNodesInElem;
	if (faceSomewhereInMesh) {	// if the face is in the island
		for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
			for (int iEle = 0; iEle < elements[iSampleNodes].size();++iEle) {	// check all elements
				faceInElement = true;
				int globalEleNum = elements[iSampleNodes][iEle];
				globalNodesInElem = elemToNodeMap.find(globalEleNum)->second;
				for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace){
					bool nodeInElement = false;
					for (int iNode = 0; iNode < 4; ++iNode) {
						if (globalNodeNum[iNodeFace] == globalNodesInElem[iNode]){
							nodeInElement = true;
							break;
						}
					}
					faceInElement *=nodeInElement;
				}
				if (faceInElement)
					break;
			}
			if (faceInElement)
				break;

		}
	}

	if (faceSomewhereInMesh && faceInElement) {
		faceInMesh = true;
	}
	
	if (globalNodeNum) delete [] globalNodeNum;
	if (nodeSomewhereInMesh) delete [] nodeSomewhereInMesh;

	return faceInMesh;
}

//----------------------------------------------

template<int dim>
bool GnatPreprocessing<dim>::checkFaceAlreadyAdded(const int cpuNum, const int
		iSub, const int iFace){
	bool includeFace;
	std::set<StaticArray <int, 3> >::iterator bcFacesInfoIt;	// {iCPU,iSub,iFace}
	int faceLoc [3] = {cpuNum, iSub, iFace};
	StaticArray <int,3> faceLocation(faceLoc);
	bcFacesInfoIt = bcFacesInfo.find(faceLocation);
	if (bcFacesInfoIt == bcFacesInfo.end())	// new face
		includeFace = true;
	else 	// already added face
		includeFace = false;
	return includeFace;
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::addFaceNodesElements(FaceSet&
		currentFaces, const int iFace, const int iSub, const int
		*locToGlobNodeMap){
	int *locNodeNums = new int [currentFaces[iFace].numNodes()];
	addNodesOnFace(currentFaces, iFace, iSub, locToGlobNodeMap, locNodeNums);// add nodes
	addElementOfFace(currentFaces, iFace, iSub, locToGlobNodeMap, locNodeNums);
	if (locNodeNums) delete [] locNodeNums ;
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::addNodesOnFace(FaceSet&
		currentFaces, const int iFace, const int iSub, const int
		*locToGlobNodeMap, int *locNodeNums = NULL){
	// output: locNodeNums

	int globalNodeNum;
	int localNodeNum;
  double xyz [3] = {0.0, 0.0, 0.0};

	for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
		localNodeNum = currentFaces[iFace][iNodeFace];
		globalNodeNum = locToGlobNodeMap[localNodeNum];
		for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) xyz[iXYZ] = 0.0;
		computeXYZ(iSub, localNodeNum, xyz);

		// add to the global set
		globalNodes[0].push_back(globalNodeNum);
		cpus[0].push_back(thisCPU);
		locSubDomains[0].push_back(iSub);
		localNodes[0].push_back(localNodeNum);
		for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) 
			nodesXYZ[iXYZ][0].push_back(xyz[iXYZ]);
		locNodeNums[iNodeFace] = localNodeNum;
	} 

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::addElementOfFace(FaceSet&
		currentFaces, const int iFace, const int iSub, const int
		*locToGlobNodeMap, const int *locNodeNums){


	bool faceInElement;
	int *locToGlobElemMap = subD[iSub]->getElemMap();

	Connectivity *nodeToEle = subD[iSub]->createNodeToElementConnectivity();
	int *nToEpointer = nodeToEle->ptr();
	int nEleConnected;
	int locEleNum, globEleNum;
	std::set<int> locEleNumsCurrent, locEleNums;
	std::vector<int> intersection;

	// loop on globalNodeNums and add element they all have in common
	for (int iLocNode = 0; iLocNode < currentFaces[iFace].numNodes(); ++iLocNode){
		locEleNumsCurrent.clear();
		int locNodeNum = locNodeNums[iLocNode];
		nEleConnected = nToEpointer[ locNodeNum + 1 ] - nToEpointer[ locNodeNum ];
			//number of elements this node belongs to

		// local element numbers connected to node
		for (int iEle = 0; iEle < nEleConnected; ++iEle) { 
			// add global element
			locEleNumsCurrent.insert(*((*nodeToEle)[locNodeNum]+iEle));
		}

		// which elements the nodes have in common
		if (iLocNode == 0) {
			locEleNums = locEleNumsCurrent;
		}
		else {
			intersection.clear();
			std::set_intersection(locEleNumsCurrent.begin(),locEleNumsCurrent.end(),
					locEleNums.begin(), locEleNums.end(),
					std::back_inserter(intersection));
			locEleNums.clear();
			for (int iIntersect = 0; iIntersect < intersection.size(); ++iIntersect)
				locEleNums.insert(intersection[iIntersect]);
		}
	}
	assert(locEleNums.size() == 1);	// must only have one element in common
	// find most common element (shared by all nodes)
	locEleNum = *(locEleNums.begin());
	elements[0].push_back(locToGlobElemMap[locEleNum]);

	// add extra node

	// loop on nodes of element, and add if it is not in the set

	Connectivity *eleToNode;
	eleToNode = subD[iSub]->createElementToNodeConnectivity();
	int locNodeNum;
	for (int iNodesConn = 0; iNodesConn < 4; ++iNodesConn){
		locNodeNum = *((*eleToNode)[locEleNum] + iNodesConn);
		int globalNodeNumTmp = locToGlobNodeMap[locNodeNum];
		elemToNode[iNodesConn][0].push_back(globalNodeNumTmp);
		bool locNodeNumNew = true;
		for (int iLocNode = 0; iLocNode < currentFaces[iFace].numNodes(); ++iLocNode) {
			if (locNodeNum == locNodeNums[iLocNode] ) {
				locNodeNumNew = false;
				break;
			}
		}
		if (locNodeNumNew)	{// must add node if it is new
			double xyz [3] = {0.0, 0.0, 0.0};
			int *locToGlobNodeMap = subD[iSub]->getNodeMap();
			int globalNodeNum = locToGlobNodeMap[locNodeNum];
			computeXYZ(iSub, locNodeNum, xyz);

			globalNodes[0].push_back(globalNodeNum);
			cpus[0].push_back(thisCPU);
			locSubDomains[0].push_back(iSub);
			localNodes[0].push_back(locNodeNum);
			for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) 
				nodesXYZ[iXYZ][0].push_back(xyz[iXYZ]);
		}
	}

	if (nodeToEle) delete nodeToEle;
}

//----------------------------------------------

template<int dim>
bool GnatPreprocessing<dim>::checkFaceContributesToLift(FaceSet& faces, const int iFace, const int iSub, const int *locToGlobNodeMap ){

	bool faceContributesToLift;
	map<int, int> surfOutMap = postOp->getSurfMap();
	int idx;
	map<int,int>::iterator it = surfOutMap.find(faces[iFace].getSurfaceID());
	if(it != surfOutMap.end() && it->second != -2)
		idx = it->second;
	else if (includeLiftFaces == 2){	// include face if it is on any moving wall (perhaps remove)
		if(faces[iFace].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
				faces[iFace].getCode() == BC_ADIABATIC_WALL_MOVING  ||
				faces[iFace].getCode() == BC_SLIP_WALL_MOVING)
			idx = 0;
		else
			idx = -1;
	}
	else
		idx = -1;

	if(idx >= 0) 
		faceContributesToLift = true;
	else
		faceContributesToLift = false;
	return faceContributesToLift; 

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputTopFile() {

	com->fprintf(stderr," ... Writing TOP file ...\n");

	// initialize file

	int sp = strlen(ioData->output.rom.prefix);
	const char *fileName;
	const char *fileNameExtension;
	determineFileName(ioData->output.rom.mesh, ".top",fileName,fileNameExtension);
	char *outMeshFile = new char[sp + strlen(fileName)+ strlen(fileNameExtension) +1];
	if (thisCPU == 0) sprintf(outMeshFile, "%s%s%s", ioData->output.rom.prefix, fileName, fileNameExtension);
	FILE *reducedMesh;
	if (thisCPU == 0) reducedMesh = fopen(outMeshFile, "wt");

	// write out globalNodes

	com->fprintf(reducedMesh,"Nodes FluidNodesRed\n");

	int globalNodeNum, globalEleNum;
	StaticArray <double, 3> xyzVals;
	StaticArray <int, 4> globalNodesTmp, reducedNodes;
	std::map<int, int > globalToReducedNodeNumbering;	// one mapping for each island

	reducedSampleNodes.resize(nSampleNodes);

	// save the reduced node number for the sample node

	for (int j = 0; j < nReducedNodes; ++j) {

		// compute xyz position of the node
		globalNodeNum = globalNodes[0][j];	// global node numbers on sample mesh have been sorted in increasing order
		xyzVals = nodesXYZmap.find(globalNodeNum)->second;
		com->fprintf(reducedMesh, "%d %8.15e %8.15e %8.15e \n", j+1, xyzVals[0], xyzVals[1], xyzVals[2]);	

		// associate global node to the current reduced node
		globalToReducedNodeNumbering.insert(pair<int, int> (globalNodeNum, j));

		// determine if a sample node
		std::map<int, int>::const_iterator sampleNodeMapLoc = globalSampleNodeRankMap.find(globalNodeNum);
		if ( sampleNodeMapLoc != globalSampleNodeRankMap.end()) {
			int globalNodeRank = sampleNodeMapLoc->second;
			reducedSampleNodes[globalNodeRank] = j;	// the globalNodeRank sample node is node j in sample mesh
			reducedSampleNodeRankMap.insert(pair<int,int>( j , globalNodeRank));
		}
	}

	// write elements

	com->fprintf(reducedMesh,"Elements FluidMeshRed using FluidNodesRed\n");

	for (int iEle = 0; iEle < elements[0].size(); ++iEle) {
		globalEleNum = elements[0][iEle];
		globalNodesTmp = elemToNodeMap.find(globalEleNum)->second;
		for (int k = 0; k < 4; ++k) {
			reducedNodes[k] = globalToReducedNodeNumbering.find(globalNodesTmp[k])->second;
		}
		com->fprintf(reducedMesh, "%d 5 %d %d %d %d \n", iEle + 1,
				reducedNodes[0]+1, reducedNodes[1]+1, reducedNodes[2]+1,
				reducedNodes[3]+1);	
	}

	// write out boundary faces

	// from sower user manual
	boundaryConditionsMap.insert(pair<int, std::string > (-5, "OutletMoving"));
	boundaryConditionsMap.insert(pair<int, std::string > (-4, "InletMoving"));
	boundaryConditionsMap.insert(pair<int, std::string > (-3, "StickMoving"));
	boundaryConditionsMap.insert(pair<int, std::string > (-2, "SlipMoving"));
	boundaryConditionsMap.insert(pair<int, std::string > (-1, "IsothermalMoving"));	// not sure (not in manual)
	boundaryConditionsMap.insert(pair<int, std::string > (0, "Internal"));
	boundaryConditionsMap.insert(pair<int, std::string > (1, "IsothermalFixed"));	// not sure (not in manual)
	boundaryConditionsMap.insert(pair<int, std::string > (2, "SlipFixed"));
	boundaryConditionsMap.insert(pair<int, std::string > (3, "StickFixed"));
	boundaryConditionsMap.insert(pair<int, std::string > (4, "InletFixed"));
	boundaryConditionsMap.insert(pair<int, std::string > (5, "OutletFixed"));
	boundaryConditionsMap.insert(pair<int, std::string > (6, "Symmetry"));

	for (int iSign = 0; iSign < 2; ++iSign) {
		for (int iBCtype = 0; iBCtype <= max(BC_MAX_CODE, -BC_MIN_CODE); ++iBCtype) {
			int maxID = 0;
			if (bcFaceSurfID[iSign][iBCtype].size() > 0) {
				std::vector<int>::iterator maxIDit = max_element(bcFaceSurfID[iSign][iBCtype].begin(),bcFaceSurfID[iSign][iBCtype].end());
				maxID = *maxIDit;
			}
			
			for (int iID = 0; iID < maxID; ++iID) {
				bool firstTime = true;
				int faceCounter = 0;
				for (int iFace = 0; iFace < bcFaces[iSign][0][iBCtype].size() ; ++iFace) {
					int currentID = bcFaceSurfID[iSign][iBCtype][iFace];
					if (currentID == iID + 1) {	// face is in the current set
						if (firstTime) {	// only output the first time (and only if size > 0)
							int boundaryCondNumber = iBCtype * ((iSign > 0)*2 - 1) ;	// returns boundary condition number
							std::string boundaryCond = boundaryConditionsMap.find(boundaryCondNumber)->second;
							com->fprintf(reducedMesh,"Elements %sSurface_%d using FluidNodesRed\n",boundaryCond.c_str(), bcFaceSurfID[iSign][iBCtype][iFace]);
							firstTime = false;
						}
						for (int k = 0; k < 3; ++k) {
							int globalNode = bcFaces[iSign][k][iBCtype][iFace];
							reducedNodes[k] = globalToReducedNodeNumbering[globalNode];
						}
						com->fprintf(reducedMesh, "%d 4 %d %d %d \n", ++faceCounter, reducedNodes[0]+1, reducedNodes[1]+1, reducedNodes[2]+1);	
					}
				}
			}
		}
	}

	if (thisCPU == 0) fclose(reducedMesh);
	if (outMeshFile) delete [] outMeshFile;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (bcFaces[i][j]) delete [] bcFaces[i][j];
		}
		if (bcFaceSurfID[i]) delete [] bcFaceSurfID[i];
	}
	if (elements) delete [] elements;
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputSampleNodes() {

	// write out sample node numbers in sample mesh node numbering system

	com->fprintf(stderr," ... writing sample node file with respect to sample mesh...\n");
	outputSampleNodesGeneral(reducedSampleNodes,ioData->output.rom.sampleNodes,".sampledNodes");

	com->fprintf(stderr," ... writing sample node file with respect to full mesh...\n");
	outputSampleNodesGeneral(globalSampleNodes,ioData->output.rom.sampleNodesFull,".sampledNodesFull");

	reducedSampleNodes.resize(0);
	globalSampleNodes.resize(0);
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputSampleNodesGeneral(const std::vector<int> &sampleNodes, const char *sampleNodeFile,const char *sampleNodeFileExtension) {

	const char *fileName;
	const char *fileNameExtension;
	determineFileName(sampleNodeFile, sampleNodeFileExtension,fileName,fileNameExtension);
	char *outSampleNodeFile = new char[strlen(ioData->output.rom.prefix) +
		strlen(fileName)+strlen(fileNameExtension)+1];
	if (thisCPU ==0) sprintf(outSampleNodeFile, "%s%s%s",
			ioData->output.rom.prefix, fileName,fileNameExtension);
	FILE *writingFile;
	if (thisCPU ==0) writingFile = fopen(outSampleNodeFile, "wt");

	com->fprintf(writingFile, "%d", nSampleNodes);	// first print number of sample nodes
	for (int i = 0; i < nSampleNodes; ++i) {
		com->fprintf(writingFile, "\n");
		com->fprintf(writingFile, "%d %d", i+1, sampleNodes[i]+1);	
	}

	delete [] outSampleNodeFile;
	if (thisCPU == 0) fclose(writingFile);

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computeXYZ(int iSub, int iLocNode, double *xyz) {
	
	// input: iSub, iLocNode
	// output: x,y,z coordinates of the iLocNode located on iSub
	DistSVec<double,3>& X = geoState->getXn();
	SVec<double,3>& Xsub = X(iSub);	// X is of type DistSVec<double,3>
	for (int i = 0; i < 3; ++i) xyz[i] = Xsub[iLocNode][i];	
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computePseudoInverse(int iPodBasis) {

//======================================
// Purpose
// 	compute pseudo inverses of podHatRes or podHatJac
// Inputs
// 	iPodBasis (0 or 1)
// Outputs
// 	the pseudo inverse
// Approach
// 	the ith column of the pseudo inverse solves min || podHat xi = ei ||_2
// 	where ei is zero except 1 at the node and dim of interest
//======================================

	// generate pseudoInvRhs in chunks

	int nNodesAtATime = ioData->gnat.pseudoInverseNodes;	
	bool lastTime = false;
	nNodesAtATime = min(nSampleNodes,nNodesAtATime);	// fix if needed
	int numRhs = nNodesAtATime * dim;

	int nHandledVectors = 0;

	pseudoInvRhs.resize(numRhs);	// make correct size (possible memory problem!)
	for (int iRhs = 0; iRhs < numRhs; ++iRhs)
		pseudoInvRhs[iRhs] = 0.0;

	if (thisCPU == 0) {
		podHatPseudoInv[iPodBasis] = new double * [nSampleNodes * dim] ;
		for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs)  
			podHatPseudoInv[iPodBasis][iRhs] = new double [nPod[iPodBasis] ] ;
	}

	// all processors store the temporary solution
	double **podHatPseudoInvTmp;

	int iVector = 0;
	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
		// compute unit vectors
		int currentGlobalNode = globalSampleNodes[iSampleNodes];
		int currentCPU = globalNodeToCpuMap.find(currentGlobalNode)->second;
		if (thisCPU == currentCPU) {
			int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
			int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
			for (int iDim = 0; iDim < dim; ++iDim) {
				SubDomainData<dim> locValue = pseudoInvRhs[iVector+iDim].subData(iSubDomain);
				locValue[iLocalNode][iDim] = 1.0;
			}
		}
		iVector += dim;
		com->barrier();	// temporary debugging
		// compute part of pseudo-inverse

		if (thisCPU == 0)	// directly store computed value on cpu 0
			podHatPseudoInvTmp = podHatPseudoInv[iPodBasis]+nHandledVectors;

		if (((iSampleNodes + 1)  % nNodesAtATime == 0 && !lastTime) || iSampleNodes == nSampleNodes - 1) {
			com->fprintf(stderr," computing pseudo inverse at sample node %d of %d \n", iSampleNodes+1, nSampleNodes);
			parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis], pseudoInvRhs,
					nPod[iPodBasis], numRhs, podHatPseudoInvTmp, false);

			nHandledVectors+=numRhs;
			if (nNodesAtATime > nSampleNodes - iSampleNodes - 1) {
				nNodesAtATime = nSampleNodes - iSampleNodes - 1;
				lastTime = true;
			}
			
			int numRhsOld = numRhs;
			numRhs = nNodesAtATime * dim;
			if (numRhs != numRhsOld)
				pseudoInvRhs.resize(numRhs);
			for (int iRhs = 0; iRhs < numRhs; ++iRhs)
				pseudoInvRhs[iRhs] = 0.0;
			iVector = 0;	// re-filling rhs
		}
	}

	if (thisCPU == 0 && iPodBasis == 0) {
		podHatPseudoInv[1] = podHatPseudoInv[0];
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computePseudoInverse() {

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)
		computePseudoInverse(iPodBasis);
	//checkConsistency();	// debugging check

	// podHat no longer needed

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)
		podHat[iPodBasis].resize(0);

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::computePodTPod() {

	// podTpod is nPod[1] x nPod[0] array
	//
	// since podTpod = I for nPodBasis = 1, only call this function for
	// nPodBasis = 2

	if (thisCPU == 0) { 
		podTpod = new double * [nPod[1]];
		for (int i = 0; i < nPod[1]; ++i)
			podTpod[i] = new double [nPod[0]];
	}

	double podTpodTmp;
	for (int i = 0; i < nPod[1]; ++i) {
		com->fprintf(stderr," ... computing podJac^TpodRes row %d of %d ...\n",i,nPod[1]);
		for (int j = 0; j < nPod[0]; ++j) { 
			podTpodTmp = pod[1][i]*pod[0][j];
			if (thisCPU == 0)
				podTpod[i][j] = podTpodTmp;
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::assembleOnlineMatrices() {

	// Purpose: assemble matrices that are used online
	// Inputs: podHatPseudoInv, podTpod
	// Outputs: onlineMatrices

	int numCols = nSampleNodes * dim;

	if (nPodBasis == 1) {	// only need to use podTpod if there are two bases (otherwise, it is identity)
		onlineMatrices[0] = podHatPseudoInv[0];	// related to the jacobian
	}
	else {	// nPod[1] != 0 because nPodBasis == 2
		int nPodJac = (nPod[1] == 0) ? nPod[0] : nPod[1];
		onlineMatrices[1] = podHatPseudoInv[1];	// related to the jacobian

		onlineMatrices[0] = new double * [numCols] ;
		for (int i = 0; i < numCols; ++i) {
			onlineMatrices[0][i] = new double [nPod[1] ] ;
			for (int iPod = 0; iPod < nPod[1]; ++iPod) 
				onlineMatrices[0][i][iPod] = 0.0;
		}
		for (int iPod = 0; iPod < nPod[1]; ++iPod) {
			for (int jPod = 0; jPod < nPod[0]; ++jPod) {
				for (int k = 0; k < numCols; ++k) { 
					onlineMatrices[0][k][iPod] += podTpod[iPod][jPod] * podHatPseudoInv[0][k][jPod];
				}
			}
		}
		// no longer need podTpod

		for (int i = 0; i < nPod[1]; ++i) {
			if (podTpod[i]) delete [] podTpod[i];
		}
		if (podTpod) delete [] podTpod;
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputOnlineMatrices() {

	// output matrices A and B in ASCII form as VecSet< DistSVec> with the
	// DistSVec defined on the sample mesh. Each column in this VecSet
	// corresponds to a row of A or B.

	outputReducedToFullNodes();

	com->fprintf(stderr," ... Writing online matrices ...\n");

	// sample mesh
	if (ioData->output.rom.onlineMatrix[0] != 0 || outputOnlineMatricesSample) 
	outputOnlineMatricesGeneral(ioData->output.rom.onlineMatrix,
			nReducedNodes, reducedSampleNodeRankMap, reducedSampleNodes);
	
	if (ioData->output.rom.onlineMatrixFull[0] != 0 || outputOnlineMatricesFull) 
		outputOnlineMatricesGeneral(ioData->output.rom.onlineMatrixFull,
				numFullNodes, globalSampleNodeRankMap, globalSampleNodes);

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		if (podHatPseudoInv[iPodBasis]) {
			for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs) {
				if (podHatPseudoInv[iPodBasis][iRhs]) delete [] podHatPseudoInv[iPodBasis][iRhs];
			}
			delete [] podHatPseudoInv[iPodBasis];
		}
	}

	if (nPodBasis == 2) {	// only need to delete memory if 2 POD bases 
												// (see assembleOnlineMatrices)
		if (onlineMatrices[0]) {
			for (int i = 0; i < nSampleNodes * dim; ++i) {
				if (onlineMatrices[0][i]) delete [] onlineMatrices[0][i];
			}
			if (onlineMatrices[0]) delete [] onlineMatrices[0];
		}
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputStateReduced() {
	//INPUTS
	// ioData, nReducedNodes, domain
	// needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
	// globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

	com->fprintf(stderr," ... Writing POD state and initial condition in sample mesh coordinates ...\n");
	int sp = strlen(ioData->output.rom.prefix);
	const char *fileName;
	const char *fileNameExtension;
	determineFileName(ioData->output.rom.podStateRed, ".sampledROBState",fileName,fileNameExtension);
	char *outPodStateFile= new char[sp + strlen(fileName) + strlen(fileNameExtension)+1];
	const char *fileName2;
	const char *fileNameExtension2;
	determineFileName(ioData->output.rom.solution, ".sol",fileName2,fileNameExtension2);
	char *outInitialConditionFile= new char[sp + strlen(fileName2) + strlen(fileNameExtension2)+1];
	if (thisCPU ==0) sprintf(outPodStateFile, "%s%s%s", ioData->output.rom.prefix, fileName,fileNameExtension);
	if (thisCPU ==0) sprintf(outInitialConditionFile, "%s%s%s", ioData->output.rom.prefix, fileName2, fileNameExtension2);

	FILE *outPodState;
	if (thisCPU ==0) outPodState = fopen(outPodStateFile, "wt");
	FILE *outInitialCondition;
	if (thisCPU ==0) outInitialCondition = fopen(outInitialConditionFile, "wt");

	com->fprintf(outInitialCondition,"Vector InitialCondition under load for FluidNodesRed\n");
	com->fprintf(outInitialCondition,"%d\n", nReducedNodes);
	
	// read in initial condition

	DistSVec<double,dim> *initialCondition = new DistSVec<double,dim>( domain.getNodeDistInfo() );
	double tmp;
	domain.readVectorFromFile(ioData->input.solutions, 0, &tmp, *initialCondition);
	outputReducedSVec(*initialCondition,outInitialCondition,0);
	if (initialCondition) delete initialCondition;

	// note: need to output the first POD basis vector twice

	com->fprintf(stderr, " ... Reading POD basis for the state ...\n");
	nPodState = ioData->gnat.nRobState;
	SetOfVec *podState = new SetOfVec(0, domain.getNodeDistInfo() );	// need to put podState in reduced coordinates
	domain.readPodBasis(input->podFileState, nPodState, *podState);	// want to read in all (nPodState should be all)
	com->fprintf(outPodState,"Vector PodState under load for FluidNodesRed\n");
	com->fprintf(outPodState,"%d\n", nReducedNodes);
	outputReducedSVec((*podState)[0],outPodState,nPodState);	// dummy output

	for (int iPod = 0; iPod < nPodState; ++iPod) {	// # rows in A and B
		outputReducedSVec((*podState)[iPod],outPodState,iPod);
	}
	if (podState) delete podState;

	if (outInitialConditionFile) delete [] outInitialConditionFile;
	if (outPodStateFile) delete [] outPodStateFile;
	if (thisCPU == 0) fclose(outPodState);
	if (thisCPU == 0) fclose(outInitialCondition );
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputReducedSVec(const DistSVec<double,dim>
		&distSVec, FILE* outFile , int iVector) {

	// INPUTS: vector, output file name, vector index
	// globalNodes, globalNodeToCpuMap, globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap

	com->barrier();
	com->fprintf(outFile,"%d\n", iVector);

	// save the reduced node number for the sample node
	for (int j = 0; j < nReducedNodes; ++j) {
		int currentGlobalNode = globalNodes[0][j];
		int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
		int iSubDomain, iLocalNode;
		if (thisCPU == iCpu){
			iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
			iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
		}

		// find value for each dim
		for (int iDim = 0; iDim < dim; ++iDim) {
			double value = 0.0; // initialize value to zero
			if (thisCPU == iCpu) {
				SubDomainData<dim> locValue = distSVec.subData(iSubDomain);
				value = locValue[iLocalNode][iDim];
			}
			com->globalSum(1, &value);
			com->fprintf(outFile,"%8.15e ", value);
		}
		com->fprintf(outFile,"\n");
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputWallDistanceReduced() {
	// INPUTS: geoState, ioData, nReducedNodes
	// needed by outputReducedVec: globalNodes, globalNodeToCpuMap, globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap

	com->fprintf(stderr," ... Writing wall distance for sample mesh ...\n");

	// load in wall distance

	DistVec<double> *d2wall = geoState->getd2wall();
	DistVec<double> d2wallOutput(d2wall->info());
	d2wallOutput = -ioData->bc.wall.delta;
		// must subtract off for input files (see DistGeoState.C constructor)
	d2wallOutput += *d2wall;

	int sp = strlen(ioData->output.rom.prefix);
	const char *fileName;
	const char *fileNameExtension;
	determineFileName(ioData->output.rom.wallDistanceRed, ".dwall",fileName,fileNameExtension);
	char *outWallDistFile = new char[sp + strlen(fileName)+ strlen(fileNameExtension) +1];
	if (thisCPU ==0) sprintf(outWallDistFile, "%s%s%s", ioData->output.rom.prefix,
			fileName, fileNameExtension);

	FILE *outWallDist;
	if (thisCPU ==0) outWallDist = fopen(outWallDistFile, "wt");

	// note: need to output the first POD basis vector twice

	com->fprintf(outWallDist,"Scalar walldist under load for FluidNodesRed\n");
	com->fprintf(outWallDist,"%d\n", nReducedNodes);
	outputReducedVec(d2wallOutput,outWallDist,0);

	if (outWallDistFile) delete [] outWallDistFile;
	if (thisCPU == 0) fclose(outWallDist);
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector) {

	com->fprintf(outFile,"%d\n", iVector);

	for (int j = 0; j < nReducedNodes; ++j) {
		double value = 0.0; // initialize value to zero
		int currentGlobalNode = globalNodes[0][j];
		int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;

		if (thisCPU == iCpu){
			int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
			int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
			double *locValue = distVec.subData(iSubDomain);
			value = locValue[iLocalNode];
		}

		com->globalSum(1, &value);
		com->fprintf(outFile,"%8.15e ", value);
		com->fprintf(outFile,"\n");
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputReducedToFullNodes() {

	int sp = strlen(ioData->output.rom.prefix);

	const char *fileName;
	const char *fileNameExtension;
	determineFileName(ioData->output.rom.reducedfullnodemap, ".reducedFullNodeMap",fileName,fileNameExtension);
	char *outMeshFile = new char[sp + strlen(fileName)+strlen(fileNameExtension)+1];
	if (thisCPU ==0) sprintf(outMeshFile, "%s%s%s", ioData->output.rom.prefix, fileName, fileNameExtension);
	FILE *outMesh;
	if (thisCPU ==0) outMesh = fopen(outMeshFile, "wt");

	// save the reduced node number for the sample node
	for (int j = 0; j < globalNodes[0].size(); ++j) {
		com->fprintf(outMesh,"%d \n", globalNodes[0][j]+1);
	}
	if (thisCPU == 0) fclose(outMesh);
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::checkConsistency() {

	// PURPOSE: debugging

	int numRhs = nSampleNodes * dim;	// number of RHS treated
	double **consistency = new double * [numRhs];
	for (int i = 0; i < numRhs; ++i)
		consistency[i] = new double [numRhs];

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
		for (int i = 0; i < numRhs; ++i){
			for (int j = 0; j < numRhs; ++j) {
				consistency[i][j] = 0.0;
			}
		}

		for (int i = 0; i < nPod[iPodBasis]; ++i){
			for (int j = 0; j < nPod[iPodBasis]; ++j) {
				int counterk = 0;
				for (int k = 0; k < numRhs; ++k) {
					int sampleNodeIndex = k/dim;
					int iDim = k%dim;
					int currentSampleNode = globalSampleNodes[sampleNodeIndex];
					int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second; 
					if (thisCPU == currentCpu){
						int currentSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second; 
						int currentLocNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second; 
						assert(domain.getSubDomain()[currentSub]->getNodeMap()[currentLocNode] == currentSampleNode);
						SubDomainData<dim> locValue;
					 	locValue	= podHat[iPodBasis][j].subData(currentSub);
						consistency[i][j] += podHatPseudoInv[iPodBasis][k][i]*locValue[currentLocNode][iDim];
						++counterk;
					}

				}
				com->globalSum(1, &counterk);
				assert(counterk==numRhs);
			}
			com->globalSum(nPod[iPodBasis], consistency[i]);
		}
		int asdf = 0;
	}

	for (int i = 0; i < numRhs; ++i) {
		if (consistency[i]) delete [] consistency[i];
	}
	if (consistency) delete [] consistency;

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::outputOnlineMatricesGeneral(const char
		*onlineMatricesName, int numNodes,
		const std::map<int,int> &sampleNodeMap, const std::vector<int>
		&sampleNodeVec) {

	if (thisCPU != 0)	return;	// only execute for cpu 0

	// prepare files

	char *onlineMatrixFile;
	FILE *onlineMatrix;
	int sp = strlen(ioData->output.rom.prefix);
	const char *(onlineMatExtension [2]) = {".gappyRes",".gappyJac"};
	int nPodJac = (nPod[1] == 0) ? nPod[0] : nPod[1];

	const char *onlineMatricesNameUse;
	const char *onlineMatricesNameUseExtension;

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {	
		determineFileName(onlineMatricesName,
				onlineMatExtension[iPodBasis],onlineMatricesNameUse,onlineMatricesNameUseExtension);
		onlineMatrixFile = new char[sp +
			strlen(onlineMatricesNameUse)+strlen(onlineMatricesNameUseExtension)+1];
		if (thisCPU ==0){
			sprintf(onlineMatrixFile, "%s%s%s",
					ioData->output.rom.prefix, onlineMatricesNameUse,
					onlineMatricesNameUseExtension); 
			onlineMatrix = fopen(onlineMatrixFile, "wt");
		}

		com->fprintf(onlineMatrix,"Vector abMatrix under load for FluidNodes\n");
		com->fprintf(onlineMatrix,"%d\n", numNodes);
		com->fprintf(onlineMatrix,"%d\n", nPodJac);	// # rows in A and B

		// dummy loop needed by POD
		for (int iNode = 0; iNode < numNodes; ++iNode) {
			for (int iDim = 0; iDim < dim; ++iDim) { 
				com->fprintf(onlineMatrix,"%8.15e ", 0.0);
			}
			com->fprintf(onlineMatrix,"\n");
		}

		for (int iPod = 0; iPod < nPodJac; ++iPod) {	// # rows in A and B
			com->fprintf(onlineMatrix,"%d\n", iPod);
			for (int iNode = 0; iNode < numNodes; ++iNode) {
				bool isSampleNode = false;
				int sampleNodeRank = -1;

				// determine if it is a sample node
				std::map<int,int>::const_iterator sampleNodeSetLoc = sampleNodeMap.find(iNode);
				if ( sampleNodeSetLoc != sampleNodeMap.end()){
					sampleNodeRank = sampleNodeMap.find(iNode)->second;
					isSampleNode = true;
				}

				for (int iDim = 0; iDim < dim; ++iDim) { 
					if (isSampleNode) {
						com->fprintf(onlineMatrix,"%8.15e ",
							onlineMatrices[iPodBasis][sampleNodeRank * dim + iDim][iPod]);
					}
					else // must output zeros if it is not a sample node
						com->fprintf(onlineMatrix,"%8.15e ", 0.0);
				}
				com->fprintf(onlineMatrix,"\n");
			}
		}
		if (onlineMatrixFile) delete [] onlineMatrixFile;
		if (thisCPU == 0) fclose(onlineMatrix);
	}
}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::readInPodResJac() {

	for (int i = 0 ; i < nPodBasis ; ++i){	// only do for number of required bases
		pod[i].resize(nPod[i]);
	}

	//	read in both bases
	com->fprintf(stderr, " ... Reading POD bases for the residual and/or Jacobian ...\n");

	domain.readPodBasis(input->podFileRes, nPod[0],pod[0]);	// always read in basis for residual

	if (nPodBasis == 2) {
		domain.readPodBasis(input->podFileJac, nPod[1],pod[1]);
	}

}

//----------------------------------------------

template<int dim>
void GnatPreprocessing<dim>::determineFileName(const char *fileNameInput, const char
		*currentExtension, const char *(&fileNameBase), const char
		*(&fileNameExtension)) {

	if (strcmp(fileNameInput,"") == 0) {
		fileNameBase = ioData->output.rom.gnatPrefix;
		fileNameExtension = currentExtension;
	}
	else {
		fileNameBase = fileNameInput;
		fileNameExtension = "";
	}
}
