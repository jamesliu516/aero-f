#include <GappyOffline.h>

template<int dim>
GappyOffline<dim>::GappyOffline(Communicator *_com, IoData &_ioData, Domain &dom, DistGeoState *_geoState) : 
	domain(dom), 	com(_com), ioData(&_ioData), 
	podRes(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	podJac(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	podState(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	initialCondition( dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	podHatRes(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	podHatJac(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	errorRes(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	errorJac(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	pseudoInvRhs(0, dom.getNodeDistInfo() ),	// two pod bases (residual + jacobian)
	handledNodes(0), nPodBasis(2),
	debugging(1),
	// distribution info
	numLocSub(dom.getNumLocSub()), nTotCpus(_com->size()), thisCPU(_com->cpuNum()),
	nodeDistInfo(dom.getNodeDistInfo()), subD(dom.getSubDomain()),parallelRom(2),
	geoState(_geoState), X(_geoState->getXn()),
	residual(0), jacobian(1)
{
	// initialize vectors to point to the approprite bases

        for(int i=0; i<2; ++i) parallelRom[i] = new ParallelRom<dim>(dom,_com);
	pod.a[0] = &podRes;	// make pod point to res and jac
	pod.a[1] = &podJac;
	podHat.a[0] = &podHatRes;	// make pod point to res and jac
	podHat.a[1] = &podHatJac;
	error.a[0] = &errorRes;	// make pod point to res and jac
	error.a[1] = &errorJac;

	handledVectors[0] = 0;	// have not yet handled any vectors
	handledVectors[1] = 0;	// have not yet handled any vectors
	nRhs[0] = 0;
	nRhs[1] = 0;

	int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE)+1;	// there is a zero condition
	for (int i = 0; i < 2; ++i){
		bcIsland[i] = new std::vector< int > [BC_CODE_EXTREME];
		for (int j = 0; j < 3; ++j)
			bcFaces[i][j] = new std::vector< int > [BC_CODE_EXTREME];
	}
	
}
template<int dim>
GappyOffline<dim>::~GappyOffline() 
{
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j)
			delete [] bcFaces[i][j];
		delete [] bcIsland[i];
	}

	for (int i = 0; i < nSampleNodes * dim; ++i)  
		delete [] onlineMatrices[0][i];
	delete [] onlineMatrices[0];
	
	delete [] cpuSet;
	delete [] locSubSet;
	delete [] locNodeSet;
	delete [] globalNodeSet;
	delete [] nodesToHandle;
	delete [] nodes;

	for (int i = 0; i < 3; ++i)
		delete [] nodesXYZ[i];
	delete [] elements;
	for (int i = 0; i < 4; ++i)
		delete [] elemToNode[i];

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs)
			delete [] podHatPseudoInv[iPodBasis][iRhs];
		delete [] podHatPseudoInv[iPodBasis];
	}
	// ktc: failing in destruction of vecset<distsvec>
}
	//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::buildGappy() {

	// KTC DEBUG
	int dbgwait = 0;
	int thisCPU = com->cpuNum();
	if (thisCPU == 0)
		while (dbgwait);
	
	setUpPodBases();

	// compute the reduced mesh used by Gappy POD

   buildGappyMesh();

	// compute matries A and B required online

	 com->barrier();
   buildGappyMatrices();

   com->fprintf(stderr," ... finished with buildGappy ...\n");

// STRATEGY
// compute the mesh for the masked domain
// 	Advantages: can use existing data structures. Can use typical decomposition techniques (e.g. DistSVec), which will make integration with ScaLapack much easier. For connectivity, assume a bunch of "islands" which do not connect.
// put zeros in the pod[0], pod[1] on the masked quantities that aren't restricted
// don't evaluate rhat and jphihat for masked quantities that aren't restricted
// the association between the restricted and full domain is through the ordering of the POD basis vectors! This is all hidden in matrices A,B

// TRICKS
// 1) you can add zero rows to matrices and it will have no effect on the qr decomposition!

// ?
// -how to build a mesh here? can you create a new object with connectivity, decomposition, etc.?
// -how to specify PODData?

// INPUT NOTES
// nPod[1] = 0 indicates that the same basis is used for the jacobian and residual

} 

//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::setUpPodBases() {

	// determine whether one or two pod bases are used

	nPodState = ioData->Rob.numROB;	// TODO want to read in all pod state vectors
	nPod[0] = ioData->Rob.numROBRes;	// number of POD basis vectors 
	nPod[1] = ioData->Rob.numROBJac;
	int podFiles [2] = {0, 1};	// which pod files should be read

	nPodMax = max(nPod[0],nPod[1]);	// compute maximum nPod

	if (nPod[1] == 0) {	// only need one basis (shared between jacobian and residual)
		nPodBasis = 1;
		pod.a[1] = &podRes;
		podHat.a[1] = &podHatRes;	// make pod point to res and jac
		error.a[1] = &errorRes;	
		podFiles[1]= -1;	// do not read the file for the Jacobian
	}
	else if (nPod[0] ==0) {
		nPodBasis = 1;
		pod.a[0] = &podJac;
		podHat.a[0] = &podHatJac;
		error.a[0] = &errorJac;	
		podFiles[0]= 1;	// only read file for the Jacobian
		podFiles[1]= -1;	// do not read the file for the Jacobian
		nPod[0] = nPod[1];	// from now on, only deal with the first basis
	}

	com->fprintf(stderr, " ... Reading POD bases for Gappy POD contruction\n");//DA
	// XXX: nSampleNodes will be an input

	int dimDbg = dim;	// KTC debug
	nSampleNodes = static_cast<int>(ceil(double(nPodMax)/double(dim)));	// this will give interpolation or the smallest possible least squares

	// require nSampleNodes * dim >= max(nPod[0],nPod[1]); nSampleNodes >= ceil(double(max(nPod[0],nPod[1]))/double(dim))
	assert(nSampleNodes * dim >= max(nPod[0],nPod[1])); // KTC debug

	for (int i = 0 ; i < nPodBasis ; ++i){	// only do for number of required bases
		pod[i].resize(nPod[i]);
		podHat[i].resize(nPod[i]);
		error[i].resize(nPod[i]);
	}


	//	read in both Pod bases
  domain.readMultiPodBasis(ioData->input.podFileResJac, pod.a, nPod, nPodBasis, podFiles);

	// read in POD basis for state and initial condition
  //ioData->output.restart.solutions
	domain.readPodBasis(ioData->input.podFile, nPodState, podState);	// want to read in all (nPodState should be all)
	double tmp;
	domain.readVectorFromFile(ioData->input.solutions, 0, &tmp, initialCondition);
}

template<int dim>
void GappyOffline<dim>::buildGappyMesh() {

	//======================================
	// PURPOSE
	// 	build reduced mesh by selecting sample nodes
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
	// Mesh will be entirely disconnected, which makes decomposition trivial
	// (tradeoff: possibly more (redundant) unknowns, but no communication
	//
	// Fix the size of PhiHat using trick 1
	//
	// First compute a node, then compute the associated mask
	//
	// First node captures the inlet boundary conditions
	//======================================

	//======================================
	// compute nRhsMax, nGreedyIt, nodesToHandle, possibly fix nSampleNodes
	//======================================

	computeGreedyIterationInfo();

	// initialize restricted mesh quantities

	cpuSet = new int[nSampleNodes];  // set of cpus containing local nodes
	locSubSet = new int[nSampleNodes];  // set of local subdomains
	locNodeSet = new int[nSampleNodes];  // set of local nodes
	globalNodeSet = new int[nSampleNodes];  // set of global nodes
	for (int i = 0; i < nSampleNodes; ++i){ cpuSet[i] = 0; locSubSet[i] = 0; locNodeSet[i] = 0; globalNodeSet[i] = 0; } // initialize
		// must be zero because do a global sum later

	//==================
	// Option 1: build on a reduced mesh, which is currently unknown
	// *Option 2*: build on the full mesh, which is currently known (easier, but slower computations)
	//==================

	// initialize restricted pod basis and error vectors

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
		podHat[iPodBasis].resize(nPod[iPodBasis]);
		for (int i = 0; i < nPod[iPodBasis]; ++i) podHat[iPodBasis][i] = 0.0;
		error[iPodBasis].resize(nRhsMax);
		for (int i = 0; i < nRhsMax; ++i) error[iPodBasis][i] = pod[iPodBasis][i];	// for the first iteration, just pick out largest element
	}

	//===============================================
	// initialize the least squares problems
	//===============================================
	
	if (nGreedyIt > 1)
		initializeGappyLeastSquares();	// no least squares the first greedy it

	//===============================================
	// run the rest of the greedy method
	//===============================================

	for (int greedyIt = 0; greedyIt < nGreedyIt; ++greedyIt) 
		greedy(greedyIt);

	// print global nodes
	if (debugging){
		com->fprintf(stderr,"globalNodeSet is:");
		for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes)
			com->fprintf(stderr,"%d ",globalNodeSet[iSampleNodes]);
	}

	//==============================================
	// add required node neighbors to the node set
	// NOTE: adding two layers of neighbors!
	//==============================================

	//KTC OK to here
	
	addTwoNodeLayers();	// add two layers of nodes 

	//==============================================
	// output TOP file
	//==============================================
	
	if (thisCPU == 0) outputTopFile();
}


//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::computeGreedyIterationInfo() {
	// KTC: verified in c++

	//==================================================================
	// PURPOSE: compute number of POD basis vectors (nRhsMax) and nodes and handled by each
	// iteration of the greedy algorithm
	// OUTPUTS
	// 	nRhsMax, nGreedyIt, nodesToHandle
	// STRATEGY:
	// nothing special happens when nSampleNodes < nPodMax < nSampleNodes * dim 
	// 	1) require nPodMax < nSampleNodes * dim to avoid underdetermined system
	// 	2) if nSampleNodes > nPodMax, need to treat more nodes per iteration
	// 	(nRhsMax is 1)
	//==================================================================
	
	//==================================================================
	// 	1) require nPodMax < nSampleNodes * dim to avoid underdetermined system
	//==================================================================

	if (nSampleNodes * dim < nPodMax) {	
		int nSampleNodesOld = nSampleNodes; 
		nSampleNodes = static_cast<int>(ceil(double(nPodMax)/double(dim))); 
		com->fprintf(stderr,"Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d",nSampleNodesOld,nSampleNodes);
	}

	nRhsMax = static_cast<int>(ceil(double(nPodMax)/double(nSampleNodes))); // nSampleNodes * nRhsMax >= max(nPod[0],nPod[1])

	// the following should always hold !this should always be true because of the above fix (safeguard)!
	
	if (nRhsMax > dim) {
		com->fprintf(stderr,"Warning: nRhsMax > dim. More nodes should have been added.");
	}

	//==================================================================
	// 2) if nSampleNodes > nPodMax, need to treat more nodes per iteration
	// strategy: fill more nodes at the earlier iterations because POD basis vectors are optimally ordered
	//==================================================================

	nGreedyIt = min(nPodMax, nSampleNodes);	// number of greedy iterations (at most nPodMax; if nSampleNodes > nPodMax need to take care of more nodes per iteration)
	nodesToHandle = new int[nGreedyIt];	// number of nodes for each greedy iteration

	for (int iGreedyIt = 0; iGreedyIt < nGreedyIt; ++iGreedyIt)	{
		nodesToHandle[iGreedyIt] = (nSampleNodes * nRhsMax) / nPodMax;
		if (iGreedyIt < nSampleNodes % nPodMax && nRhsMax ==1)	// only in the dangerous case with nRhsMax = 1
			++nodesToHandle[iGreedyIt];
	}

}

//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::findMaxAndFillPodHat(double myMaxNorm, int locSub, int locNode, int globalNode) {

	//===============================================
	// PURPOSE: fill locPodHat
	// INPUTS
	// 	Local: myMaxNorm, locSub, locNode, globalNode
	// OUTPUTS
	// 	Global: locPodHat for maximum entry, globally summed cpuSet, locSubSet,
	// 	locNodeSet, globalNodeSet, xyz
	//===============================================
		
	int thisCPU = com->cpuNum();
	double globalMaxNorm = myMaxNorm;
	com->barrier();
	com->globalMax(1, &globalMaxNorm);  // find the maximum value over all cpus
	double *xyz = new double [3];
	for (int i=0; i<3; ++i)
		xyz[i]=0.0;

	if (myMaxNorm == globalMaxNorm) {  // if this CPU has the maximum value

		if (debugging){
		 com->fprintf(stderr, "CPU %d has myMaxNorm with locNode %d, locSub %d \n",thisCPU,locNode,locSub);
		}

		// save the global subdomain and local node indices (sum at the very end of
		// algorithm)

		cpuSet[handledNodes] = thisCPU;
		locSubSet[handledNodes] = locSub;
		locNodeSet[handledNodes] = locNode;
		globalNodeSet[handledNodes] = globalNode;
		computeXYZ(locSub, locNode, xyz);

		// fill out restricted matrices (all columns for the current rows)

		SubDomainData<dim> locPod, locPodHat;

		for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
			for (int iPod = 0 ; iPod < nPod[iPodBasis]; ++iPod) {
				locPod = pod[iPodBasis][iPod].subData(locSub);	// cannot access iDim entry
				locPodHat = podHat[iPodBasis][iPod].subData(locSub);
				for (int iDim = 0; iDim < dim ; ++iDim) {
					for (int iPod = 0; iPod < nPod[0] ; ++iPod) {
						locPodHat[locNode][iDim] = locPod[locNode][iDim];	// zeros everywhere except at the chosen sample nodes
					}
				}
			}
		}
	}

	// make sure all cpus have the same copy
	
	com->barrier();
	com->globalSum(1, cpuSet + handledNodes);
	com->globalSum(1, locSubSet + handledNodes);
	com->globalSum(1, locNodeSet + handledNodes);
	com->globalSum(1, globalNodeSet + handledNodes);
	com->globalSum(3, xyz);
	StaticArray<double, 3> XYZ(xyz); 
	delete [] xyz;
	nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalNodeSet[handledNodes], XYZ));

	++handledNodes;
}

template<int dim>
void GappyOffline<dim>::greedy(int greedyIt) {

	// Differences for 1st iteration compared with other greedy iterations:
	// 1) no least squares problem is solved (just take the maximum entry)
	// 2) look at the inlet face to ensure boundary condition is handled
	//

	double myMaxNorm;
	int locSub, locNode, globalNode;  // temporary global subdomain and local node

	bool doLeastSquares = true;
	bool onlyInletOutletBC = false;

	if (greedyIt == 0) {	
		doLeastSquares = false;// don't do least squares if first iteration
		onlyInletOutletBC = true;// first iteration, only consider inlet BC
	}

	// determine number of rhs for each
	// typically have nRhs = nRhsMax. exception: there are less than nRhsMax remaining vectors

	for (int iPodBasis = 0; iPodBasis  < nPodBasis; ++iPodBasis)	
		nRhs[iPodBasis] = min(nRhsMax, max(nPod[iPodBasis] - handledVectors[iPodBasis], 0));	

	if (doLeastSquares) {
		leastSquaresReconstruction();		// solve the least squares reconstruction
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

	}

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)	
		handledVectors[iPodBasis] += nRhs[iPodBasis];
}


template<int dim>
void GappyOffline<dim>::initializeGappyLeastSquares() {

	// initialize least squares problems

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		parallelRom[iPodBasis]->parallelLSMultiRHSInit(podHat[iPodBasis], error[iPodBasis]);
	}

}

template<int dim>
void GappyOffline<dim>::makeNodeMaxIfUnique(double nodeError, double &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int &globalNode) {

	// PURPOSE: make the node the current maximum on this cpu if it hasn't been added already
	// INPUT
	// 	Local: nodeError, myMaxNorm (can change), iSub, locNodeNum (in subdomain node numbering
	// system)
	// 	Global: handledNodes, globalNodeSet,
	// OUTPUT: 
	// 	Local: myMaxNorm (can change), locSub, locNode, globalNode
	
	// only do if the node could be the maximum

	if (nodeError > myMaxNorm) {

		bool newNode = true;
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();
		int thisGlobalNode = locToGlobNodeMap[locNodeNum];

		// check that the node hasn't already been added (look at globalNodeSet)

		 for (int iIslandCheck = 0; iIslandCheck < handledNodes; ++iIslandCheck) {
			 if (thisGlobalNode == globalNodeSet[iIslandCheck]) { 
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

template<int dim>
void GappyOffline<dim>::computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError) {

	// PURPOSE: compute the sum of squares of node error for all RHS, both bases
	// 	at the locNodeNum node
	// INPUT
	// 	Local: locMasterFlag, locNodeNum 
	// 	Global: nRhs
	// OUTPUT
	// 	Global: nodeError

	 nodeError = 0.0; // initialize normed error to zero

	 if (locMasterFlag[locNodeNum]) { 	// use subdomain node number
		 for (int iPodBasis = 0; iPodBasis  < nPodBasis ; ++iPodBasis)
			 for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs){  // add components of error from all vectors on this node
				 for (int k = 0; k < dim ; ++k){	// add contributions from residual and jacobian reconstruction errors where possible
					 nodeError += pow(locError[iPodBasis][iRhs][locNodeNum][k],2);
			 }
		 }
	 }
}

template<int dim>
void GappyOffline<dim>::getSubDomainError(int iSub) {
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

template<int dim>
void GappyOffline<dim>::leastSquaresReconstruction() {

	// PURPOSE: compute least squares reconstruction error of the pod basis
	// vectors
	// INPUT
	// 	Global: nPodBasis, nRhs, error, podHat, error
	//
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
			delete [] lsCoeff[iPodBasis][iRhs];
		delete [] lsCoeff[iPodBasis];

	}
}

template<int dim>
void GappyOffline<dim>::subDFindMaxError(int iSub, bool onlyInletOutletBC, double &myMaxNorm, int &locSub, int &locNode, int &globalNode) {

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
		for (int locNodeNum = 0; locNodeNum < nLocNodes ; ++locNodeNum){// local node number (in subdomain numbering) from the face
			computeNodeError(locMasterFlag, locNodeNum, nodeError);	// computes nodeError
			// make the node max if it is the maximum so far and has not yet been added
			 makeNodeMaxIfUnique(nodeError, myMaxNorm, iSub, locNodeNum, locSub, locNode, globalNode);
		}
	}
}

template<int dim>
void GappyOffline<dim>::parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff) {
		parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis],error[iPodBasis],
				handledVectors[iPodBasis], nRhs[iPodBasis], lsCoeff);
}
template<int dim>
void GappyOffline<dim>::addTwoNodeLayers() {

	// another option: compute H2, find which nodes the nonzeros correspond to, then find which elements have all their nodes in that set

	// nodes[iIsland][iNode] is the iNode neighbor of iIsland.
	// nodes[iIsland][0] is the sample node itself 
	// NOTE: globalNodeSet has already been communicated

	nodes = new std::vector <int> [nSampleNodes];
	for (int i = 0; i < 3; ++i)
		nodesXYZ[i] = new std::vector <double> [nSampleNodes];
	elements = new std::vector <int> [nSampleNodes];
	for (int i = 0; i < 4; ++i)
		elemToNode[i] = new std::vector <int> [nSampleNodes];
	
	// compute two layers of nodes
	
	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
		
		// add the sample node itself
		nodes[iSampleNodes].push_back(globalNodeSet[iSampleNodes]);
		StaticArray<double, 3> xyzVals = nodesXYZmap.find(globalNodeSet[iSampleNodes])->second;
		for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
			nodesXYZ[iXYZ][iSampleNodes].push_back(xyzVals[iXYZ]);
		}

		// add all neighbor nodes and elements of the sample node itself
		addNeighbors(iSampleNodes);
	}

	communicateAll();

	// add all neighbor nodes and elements of the sample node's neighbors 
	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) 
		addNeighbors(iSampleNodes, 1);

	communicateAll();

	// globally sum them up (make consistent across all cpus)

	defineMaps();	// must be done before makeUnique (this will sort in a different order) and after communication (all procs have same info) 

	// remove redundant entries from the nodes and elements

	makeUnique(nodes);
	makeUnique(elements);

	// compute faces on boundary of reduced mesh

	computeBCFaces();
	communicateBCFaces();
} 

template<int dim>
void GappyOffline<dim>::addNeighbors(int iIslands, int startingNodeWithNeigh = 0) {

	// add all global neighbor nodes/elements in the iIslands row of and elements to the iIsland node set
	
	Connectivity *nodeToNode, *nodeToEle, *eleToNode;
	double xyz [3] = {0.0, 0.0, 0.0};
	int globalNodeNum;	// NOTE: must work with global node numbers because the greedy selection only operated on master nodes. Here, we care about the node even if it wasn't a master node.
	int locNodeNum, locEleNum;	// local node number of the node/element to be added
	int nNodesToAdd, nEleToAdd;	// number of nodes that should be added
	bool elemBasedConnectivityCreated;	// whether createElemBasedConnectivity has been created for the current subdomain
	int *nToNpointer, *nToEpointer, *eToNpointer;	// pointers to connectivity information

	int nNodeWithNeigh = nodes[iIslands].size();	// number of nodes whose neighbors will be added

	for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {	// all subdomains
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();
		int *locToGlobElemMap = subD[iSub]->getElemMap();
		elemBasedConnectivityCreated = false;
		for (int iLocNode = 0; iLocNode < subD[iSub]->numNodes(); ++iLocNode) {	// all local nodes in subdomain
			globalNodeNum = locToGlobNodeMap[iLocNode];	// global node number
			for (int iNodeWithNeigh = startingNodeWithNeigh; iNodeWithNeigh < nNodeWithNeigh; ++iNodeWithNeigh) {	// check if this local node is in the current row
				if (globalNodeNum == nodes[iIslands][iNodeWithNeigh]) {// add all neighbors of this node
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
						nodes[iIslands].push_back(locToGlobNodeMap[locNodeNum]);
						for (int iCrap = 0; iCrap < 3 ; ++iCrap) xyz[iCrap] = 0.0; // KTCREMOVE
						computeXYZ(iSub, locNodeNum, xyz);
						for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
							nodesXYZ[iXYZ][iIslands].push_back(xyz[iXYZ]);
						}
					}
					for (int iEleToAdd = 0; iEleToAdd < nEleToAdd; ++iEleToAdd) { 
						// add global element
						locEleNum = *((*nodeToEle)[iLocNode]+iEleToAdd);
						elements[iIslands].push_back(locToGlobElemMap[locEleNum]);
						
						// determine nodes connected to the element
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
	// delete nodeToNode; 
	// delete nodeToEle; 
	// delete eleToNode; 
}

template<int dim>
void GappyOffline<dim>::computeBCFaces() {

	// PURPOSE: determine which faces of the reduced mesh are on the boundary
	// NOTE: this is done in parallel; thus, a communication is needed to make
	// sure all cpus have the same copy
		// determine if any faces are boundary conditions
		// KTC: how to check which faces are attached to an element???

	int faceBCCode = 0, whichIsland;
	bool faceInMesh;
	int globalFaceNodes [3];	// global node numbers of current face
	for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {	// all subdomains
		FaceSet& currentFaces = subD[iSub]->getFaces();	// faces on subdomain
		int *locToGlobNodeMap = subD[iSub]->getNodeMap();	// global node numbering
		for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {	// check all faces	
			faceBCCode = currentFaces[iFace].getCode();
			int codeIsPos = faceBCCode >0;
			if (faceBCCode != 0) {	// only do for boundary faces
				
				// check if the face is in the reduced mesh
				 checkFaceInMesh(currentFaces, iFace, iSub, locToGlobNodeMap, faceInMesh, whichIsland);

				if (faceInMesh) {
					// add face to extra set
					for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) { 
						int localFaceNode = currentFaces[iFace][iFaceNode];
						globalFaceNodes[iFaceNode] = locToGlobNodeMap[localFaceNode];
					}
					for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) {
						bcFaces[codeIsPos][iFaceNode][abs(faceBCCode)].push_back(globalFaceNodes[iFaceNode]);
					}
					bcIsland[codeIsPos][abs(faceBCCode)].push_back(whichIsland);
				}
			}
		}
	}
}

template<int dim>
template<typename Scalar>
void GappyOffline<dim>::communicateMesh(std::vector <Scalar> * nodeOrEle, int arraySize){
	
	// loop over iIslands
		// figure out how many total entries each cpu has for the iIsland
		// 	numNeigh = {0 9 15 58} means that the first cpu has 9, second has 6, etc.
		// initiate memory for the total number of nodes (using last entry in above vector)
		// fill out entries [iCpu] to [iCpu+1] in above array using vector
		// do a global sum
		// overwrite node and element vectors with this global data

	int* numNeigh = new int [nTotCpus + 1];
	for (int iArraySize = 0; iArraySize < arraySize; ++iArraySize) {
		for (int i = 0; i <=nTotCpus; ++i)	// initialize
			numNeigh[i] = 0;
		numNeigh[thisCPU+1] = nodeOrEle[iArraySize].size();	// number of entries on this cpu
		com->globalSum(nTotCpus+1,numNeigh);
		for (int i = 1; i <=nTotCpus; ++i)	// accumulate
			numNeigh[i] += numNeigh[i-1];
		int totalNodeOrEle = numNeigh[nTotCpus];	// total across all processors

		Scalar *nodeOrEleArray = new Scalar [totalNodeOrEle];
		for (int iNeighbor = 0; iNeighbor < totalNodeOrEle; ++iNeighbor) {
			if (iNeighbor >= numNeigh[thisCPU] && iNeighbor < numNeigh[thisCPU+1]) 
				nodeOrEleArray[iNeighbor] = nodeOrEle[iArraySize][iNeighbor - numNeigh[thisCPU]];	// fill in this cpu's contribution
			else
				nodeOrEleArray[iNeighbor] = 0;
		}

		com->globalSum(totalNodeOrEle,nodeOrEleArray);

		// fill in the array with all global entries
		nodeOrEle[iArraySize].clear();
		for (int iNeighbor = 0; iNeighbor < totalNodeOrEle; ++iNeighbor) 
			nodeOrEle[iArraySize].push_back(nodeOrEleArray[iNeighbor]);

		delete [] nodeOrEleArray;
		
	}
	delete [] numNeigh;
}

template<int dim>
void GappyOffline<dim>::communicateAll() {

	communicateMesh(nodes, nSampleNodes);
	for (int i = 0; i < 3; ++i)
		communicateMesh(nodesXYZ[i], nSampleNodes);
	communicateMesh(elements, nSampleNodes);
	for (int i = 0; i < 4; ++i)
		communicateMesh(elemToNode[i], nSampleNodes);

}
template<int dim>
void GappyOffline<dim>::defineMaps() {

	// defines nodesXYZmap and elemToNodeMap

	int globalNodeNumTmp;
	StaticArray<double, 3> nodesXYZTmp;
	int globalEleNumTmp;
	StaticArray<int, 4> elemToNodeTmp;

	for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {

		// define nodesXYZmap
		for (int iNeighbor = 0; iNeighbor < nodes[iSampleNodes].size(); ++iNeighbor) {
			globalNodeNumTmp = nodes[iSampleNodes][iNeighbor];
			for (int iXYZ = 0 ; iXYZ < 3; ++iXYZ)
				nodesXYZTmp[iXYZ] = nodesXYZ[iXYZ][iSampleNodes][iNeighbor];
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
}

template<int dim>
void GappyOffline<dim>::makeUnique( std::vector <int> * nodeOrEle, int length = -1) {
	// remove redundant entries from a vector <int> nodeOrEle *
	// apply to nodeOrEle and elements

	if (length == -1) length = nSampleNodes;
	vector<int>::iterator it;

	for (int iNodes = 0; iNodes < length; ++iNodes) {
		sort(nodeOrEle[iNodes].begin(), nodeOrEle[iNodes].end());	// sort: puts in order
		it = unique(nodeOrEle[iNodes].begin(), nodeOrEle[iNodes].end()); // remove duplicate consecutive elements (reason for sort)
		nodeOrEle[iNodes].resize(it - nodeOrEle[iNodes].begin());	// remove extra entries
	}

}

template<int dim>
void GappyOffline<dim>::communicateBCFaces(){
	
	int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE);

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j) {
			communicateMesh(bcFaces[i][j],BC_CODE_EXTREME );
		}
		communicateMesh(bcIsland[i],BC_CODE_EXTREME );
	}
}


template<int dim>
void GappyOffline<dim>::checkFaceInMesh(FaceSet& currentFaces, int iFace, int iSub, int *locToGlobNodeMap , bool &faceInMesh, int &whichIsland){

	// PURPOSE: determine wheteher or not currentFace is in the reduced mesh
	// OUTPUT: faceInMesh = true if the face is in the reduced mesh; whichIsland

	faceInMesh = false;
	bool faceInIsland;
	bool *nodeInIsland = new bool [currentFaces[iFace].numNodes()];	// one bool for each face node
	int * globalNodeNum = new int [currentFaces[iFace].numNodes()];
	for (int iIsland = 0; iIsland < nSampleNodes ; ++iIsland) { 	// islands in reduced mesh
		for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // nodes on face
			nodeInIsland[iNodeFace] = false;
		}
		for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // nodes on face
			globalNodeNum[iNodeFace] = locToGlobNodeMap[currentFaces[iFace][iNodeFace]];
		}
		for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // nodes on face
			// check to see if the iNodeFace is in iIsland
			for (int iIslandNode = 0; iIslandNode < nodes[iIsland].size(); ++iIslandNode) { // nodes on island
				if (globalNodeNum[iNodeFace] == nodes[iIsland][iIslandNode]){ 
					nodeInIsland[iNodeFace] = true;
					break;
				}
			}
		}

		faceInIsland = true;
		for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace)
			faceInIsland *= nodeInIsland[iNodeFace];

		//check if a given element has all the nodes

		bool faceInElement = true;

		// KTC: sometimes getting faces that are not part of an element
		// 	this checks to make sure an element owns part of the face
		
		//if (faceInIsland) {	// if the face is in the island
		//	for (int iEle = 0; iEle < elements[iIsland].size();++iEle) {
		//		faceInElement = true;
		//		for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace){
		//			bool nodeInElement = false;
		//			for (int iNode = 0; iNode < 4; ++iNode) {
		//				if (globalNodeNum[iNodeFace] == elemToNode[iNode][iIsland][iEle]){
		//					nodeInElement = true;
		//					break;
		//				}
		//			}
		//			faceInElement *=nodeInElement;
		//		}
		//		if (faceInElement)
		//			break;
		//	}
		//}
		if (faceInIsland && faceInElement) {
			faceInMesh = true;
			whichIsland = iIsland;
			break;
		}
	}

	delete [] nodeInIsland;
	delete [] globalNodeNum;
}

//
//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::outputTopFile() {

	// KTC CURRENT

   com->fprintf(stderr," ... Printing TOP file ...\n");
   int sp = strlen(ioData->output.transient.prefix);
   char *outMeshFile = new char[sp + strlen(ioData->output.transient.mesh)+1];
   sprintf(outMeshFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.mesh);
   FILE *reducedMesh = fopen(outMeshFile, "w");

	 // write out nodes

   com->fprintf(reducedMesh,"Nodes FluidNodesRed\n");
	 reducedNodeCount = 0;
	 int  reducedEleCount = 0, reducedFaceCount = 0, globalNodeNum, globalEleNum;
	 StaticArray <double, 3> xyzVals;
	 StaticArray <int, 4> globalNodes, reducedNodes;
	 std::map<int, int > globalToReducedNodeNumbering [nSampleNodes];	// one mapping for each island
	 sampleToReducedNodeNumbering = new int [nSampleNodes];

   for (int i = 0; i < nSampleNodes; ++i) {
		 // save the reduced node number for the sample node
     for (int j = 0; j < nodes[i].size(); ++j) {
			 ++reducedNodeCount;
			 // compute xyz position of the node
			 globalNodeNum = nodes[i][j];
			 if (globalNodeNum == globalNodeSet[i])
				 sampleToReducedNodeNumbering[i] = reducedNodeCount;
			 for (int iXYZ = 0; iXYZ < 3; ++iXYZ) xyzVals[iXYZ] = 0.0;	//KTCREMOVE
			 xyzVals = nodesXYZmap.find(globalNodeNum)->second;
       com->fprintf(reducedMesh, "%d %e %e %e \n", reducedNodeCount, xyzVals[0], xyzVals[1], xyzVals[2]);	
       //com->fprintf(reducedMesh, "%d %d %e %e %e \n", reducedNodeCount, globalNodeNum+1, xyzVals[0], xyzVals[1], xyzVals[2]);	
			 // associate global node to the current reduced node
			 globalToReducedNodeNumbering[i].insert(pair<int, int> (globalNodeNum, reducedNodeCount));
     }
   }

	 // write out elements

   com->fprintf(reducedMesh,"Elements FluidMeshRed using FluidNodesRed\n");

   for (int i = 0; i < nSampleNodes; ++i) {
     for (int j = 0; j < elements[i].size(); ++j) {
			 ++reducedEleCount;
			 globalEleNum = elements[i][j];
			 globalNodes = elemToNodeMap.find(globalEleNum)->second;
			 for (int k = 0; k < 4; ++k) {
				reducedNodes[k] = globalToReducedNodeNumbering[i].find(globalNodes[k])->second;
			 }
       com->fprintf(reducedMesh, "%d 5 %d %d %d %d \n", reducedEleCount, reducedNodes[0], reducedNodes[1], reducedNodes[2], reducedNodes[3]);	
       //com->fprintf(reducedMesh, "%d 5 %d %d %d %d \n", globalEleNum + 1, globalNodes[0]+1, globalNodes[1]+1, globalNodes[2]+1, globalNodes[3]+1);	
		 }
	 }

	 // write out boundary faces

		// from sower user manual
		
		boundaryConditionsMap.insert(pair<int, std::string > (-5, "OutletMoving"));
		boundaryConditionsMap.insert(pair<int, std::string > (-4, "InletMoving"));
		boundaryConditionsMap.insert(pair<int, std::string > (-3, "StickMoving"));
		boundaryConditionsMap.insert(pair<int, std::string > (0, "Internal"));
		boundaryConditionsMap.insert(pair<int, std::string > (2, "SlipFixed"));
		boundaryConditionsMap.insert(pair<int, std::string > (3, "StickFixed"));
		boundaryConditionsMap.insert(pair<int, std::string > (4, "InletFixed"));
		boundaryConditionsMap.insert(pair<int, std::string > (5, "OutletFixed"));
		boundaryConditionsMap.insert(pair<int, std::string > (6, "Symmetry"));

		for (int iSign = 0; iSign < 2; ++iSign) {
			for (int iBCtype = 0; iBCtype <= max(BC_MAX_CODE, -BC_MIN_CODE); ++iBCtype) {
				reducedFaceCount = 0;	// reduced faces for this bc type
				for (int iFace = 0; iFace < bcFaces[iSign][0][iBCtype].size() ; ++iFace) {
					++reducedFaceCount;
					if (iFace == 0) {	// only output the first time (and only if size > 0)
						int boundaryCondNumber = iBCtype * ((iSign > 0)*2 - 1) ;	// returns boundary condition number
						std::string boundaryCond = boundaryConditionsMap.find(boundaryCondNumber)->second;
						com->fprintf(reducedMesh,"Elements %sSurfaceRed using FluidNodesRed\n",boundaryCond.c_str());
					}
					int whichIsland = bcIsland[iSign][iBCtype][iFace];	// island defines global to reduced node mapping
					for (int k = 0; k < 3; ++k) {
						int globalNode = bcFaces[iSign][k][iBCtype][iFace];
						reducedNodes[k] = globalToReducedNodeNumbering[whichIsland][globalNode];
					}
					com->fprintf(reducedMesh, "%d 4 %d %d %d \n", reducedFaceCount, reducedNodes[0], reducedNodes[1], reducedNodes[2]);	
				}
			}
		}

	 // write out sample node numbers in reduced mesh node numbering system
	 
   com->fprintf(stderr," ... Printing sample node mapping file ...\n");
   sp = strlen(ioData->output.transient.prefix);
   char *outSampleNodeFile = new char[sp + strlen(ioData->output.transient.sampleNodes)+1];
   sprintf(outSampleNodeFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.sampleNodes);
   FILE *sampleNodeFile = fopen(outSampleNodeFile, "wt");

	 com->fprintf(sampleNodeFile, "%d", nSampleNodes);	// first print number of sample nodes
   for (int i = 0; i < nSampleNodes; ++i) {
		 com->fprintf(sampleNodeFile, "\n");	// first print number of sample nodes
		 com->fprintf(sampleNodeFile, "%d %d", i+1, sampleToReducedNodeNumbering[i]);	
   }
	 
	 // write out sample node numbers in full mesh node numbering system

   com->fprintf(stderr," ... Printing sample node file with respect to full mesh...\n");
   sp = strlen(ioData->output.transient.prefix);
   char *outSampleNodeGlobFile = new char[sp + strlen(ioData->output.transient.sampleNodesGlob)+1];
   sprintf(outSampleNodeGlobFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.sampleNodesGlob);
   FILE *sampleNodeGlobFile = fopen(outSampleNodeGlobFile, "wt");

	 com->fprintf(sampleNodeGlobFile, "%d", nSampleNodes);	// first print number of sample nodes
   for (int i = 0; i < nSampleNodes; ++i) {
		 com->fprintf(sampleNodeGlobFile, "%\n");	// first print number of sample nodes
		 com->fprintf(sampleNodeGlobFile, "%d %d", i+1, globalNodeSet[i]+1);	
   }
	 delete [] outMeshFile;
	 delete [] outSampleNodeFile;
	 delete [] outSampleNodeGlobFile;

	 //delete [] sampleToReducedNodeNumbering;
}

//---------------------------------------------------------------------------------------
template<int dim>
void GappyOffline<dim>::computeXYZ(int iSub, int iLocNode, double *xyz) {
	
	// input: iSub, iLocNode
	// output: x,y,z coordinates of the iLocNode located on iSub
	DistSVec<double,3>& X = geoState->getXn();
	SVec<double,3>& Xsub = X(iSub);	// X is of type DistSVec<double,3>
	for (int i = 0; i < 3; ++i) xyz[i] = Xsub[iLocNode][i];	
}

//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::buildGappyMatrices() {

//======================================
// Inputs
// 	full domain decomposition: pod[0], pod[1]
// 	reduced domain decomposition: podHat[0], podHat[1]
// Outputs
// 	reduced domain decomposition: A, B
//======================================

	// compute podHatPseudoInv (pseudo-inverses of podHatRes and podHatJac)
	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)
		computePseudoInverse(iPodBasis);
	
	// compute pod[0]^Tpod[1]
	computePodTPod();

	// assemble onlineMatrices (A and B)
	assembleOnlineMatrices();

	// output matrices A and B in ASCII form as VecSet< DistSVec> with the
	// DistSVec defined on the reduced mesh. Each column in this VecSet
	// corresponds to a row of A or B.
  numTotNodes = domain.getNumGlobNode();	// # nodes in full mesh
	if (thisCPU == 0) outputOnlineMatrices();
	if (thisCPU == 0)  outputStateReduced();
	if (thisCPU == 0)  outputWallDistanceReduced();
}

//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::computePseudoInverse(int iPodBasis) {

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

	int numRhs = nSampleNodes * dim;	// number of RHS treated

	// compute rhs matrix pseudoInvRhs
	// 	for each column, all zeros but one if it is the [iSampleNode][iDim]
	// 	location

	computePseudoInverseRHS(); // pseudoInvRhs.resize(numRhs);

	// allocate memory for pseudo-inverse with dimension: (numRhs) x nPod 

	podHatPseudoInv[iPodBasis] = new double * [numRhs] ;
	for (int iRhs = 0; iRhs < numRhs; ++iRhs)  
		podHatPseudoInv[iPodBasis][iRhs] = new double [nPod[iPodBasis] ] ;

	// compute pseudo-inverse
	// numRhs = nSampleNodes * dim
	// A: (nSampleNodes * dim) x nPod
	parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis], pseudoInvRhs,
			nPod[iPodBasis], numRhs, podHatPseudoInv[iPodBasis]);
	//parallelRom[iPodBasis].parallelLSMultiRHS(podHat[iPodBasis],error[iPodBasis],
//				handledVectors[iPodBasis], nRhs[iPodBasis], lsCoeff);

	if (iPodBasis == 0)	// set them equal by default
		podHatPseudoInv[1] = podHatPseudoInv[0];
}

//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::computePseudoInverseRHS() {

	int numRhs = nSampleNodes * dim;

	SubDomainData<dim> locData;

	// initialize pseudoInvRhs to have numRhs columns of zeros
	pseudoInvRhs.resize(numRhs);	// make correct size
	for (int i = 0; i < numRhs; ++i) pseudoInvRhs[i] = 0.0;
	
	// determine which entries should be one
	for (int iSampleNode = 0; iSampleNode < nSampleNodes; ++iSampleNode) { 
		for (int iDim = 0; iDim < dim ; ++iDim) {
			if (thisCPU == cpuSet[iSampleNode]) {
				int iRhs = iSampleNode * dim + iDim;
				// pseudoInvRhs[iVec][iSub][iLocNode][iDim] = 1.0
				locData = pseudoInvRhs[iRhs].subData(locSubSet[iSampleNode]);
				locData[locNodeSet[iSampleNode]][iDim] = 1.0;
			}
		}
	}
}

//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::computePodTPod() {

	// podTpod is nPod[1] x nPod[0] array

	podTpod = new double * [nPod[1]];
	for (int i = 0; i < nPod[1]; ++i)
		podTpod[i] = new double [nPod[0]];

	for (int i = 0; i < nPod[1]; ++i)
		for (int j = 0; j < nPod[0]; ++j)
			podTpod[i][j] = pod[1][i]*pod[0][j];
}

template<int dim>
void GappyOffline<dim>::assembleOnlineMatrices() {

	// Purpose: assemble matrices that are used online
	// Inputs: podHatPseudoInv, podTpod
	// Outputs: onlineMatrices

	int numCols = nSampleNodes * dim;

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

	for (int i = 0; i < nPod[1]; ++i)
		delete [] podTpod[i];
	delete [] podTpod;

}

template<int dim>
void GappyOffline<dim>::outputOnlineMatrices() {

   com->fprintf(stderr," ... Printing online matrices file ...\n");
   int sp = strlen(ioData->output.transient.prefix);
   char *(onlineMatrixFile [2]);
   FILE *(onlineMatrix [2]);
	 onlineMatrixFile[0] = new char[sp + strlen(ioData->output.transient.bMatrix)+1];
	 onlineMatrixFile[1] = new char[sp + strlen(ioData->output.transient.aMatrix)+1];

	 sprintf(onlineMatrixFile[0], "%s%s", ioData->output.transient.prefix, ioData->output.transient.bMatrix);
	 sprintf(onlineMatrixFile[1], "%s%s", ioData->output.transient.prefix, ioData->output.transient.aMatrix);

	 // write each file

	 // KTC: could make the upper limit nPodBasis because the files should be
	 // the same if phir and phij are the same (leave it like this for now)
	 
	 // loop over reduced nodes and output the correct value if it is a sample
	 // node. Otherwise, output a zero. This is to make it in DistSVec form for
	 // the online part.

	 // ----------------------
	 // output in reduced mesh
	 // ----------------------
	 
	 bool isSampleNode;	// is the current node is a sample node?
	 int iSampleNode = 0;	// counter for passing sample nodes
	 for (int iPodBasis = 0; iPodBasis < 2; ++iPodBasis) {	// always write out 2 files for now

		 onlineMatrix[iPodBasis] = fopen(onlineMatrixFile[iPodBasis], "wt");

		 com->fprintf(onlineMatrix[iPodBasis],"Vector abMatrix under load for FluidNodesRed\n");
		 com->fprintf(onlineMatrix[iPodBasis],"%d\n", reducedNodeCount);
		 com->fprintf(onlineMatrix[iPodBasis],"%d\n", nPod[1]);	// # rows in A and B

		// dummy loop needed by POD

		 for (int iReducedNode = 0; iReducedNode < reducedNodeCount; ++iReducedNode) {
			 for (int iDim = 0; iDim < dim; ++iDim) { 
				 com->fprintf(onlineMatrix[iPodBasis],"%e ", 0.0);
			 }
			 com->fprintf(onlineMatrix[iPodBasis],"\n");
		 }

		 for (int iPod = 0; iPod < nPod[1]; ++iPod) {	// # rows in A and B
			 com->fprintf(onlineMatrix[iPodBasis],"%d\n", iPod);
			 iSampleNode = 0;	// reset counter for sample nodes
			 for (int iReducedNode = 0; iReducedNode < reducedNodeCount; ++iReducedNode) {
				 isSampleNode = false;
				 // determine if it is a sample node
																 // TODO check the first condition on the following line
				 if (iSampleNode < nSampleNodes && iReducedNode + 1 == sampleToReducedNodeNumbering[iSampleNode]) {
					 isSampleNode = true;
					 ++iSampleNode;	// we have passed a sample node
				 }
				 for (int iDim = 0; iDim < dim; ++iDim) { 
					 if (isSampleNode)
						 com->fprintf(onlineMatrix[iPodBasis],"%e ",
								 onlineMatrices[iPodBasis][(iSampleNode - 1) * dim + iDim][iPod]);
					 else // must output zeros if it is not a sample node
						 com->fprintf(onlineMatrix[iPodBasis],"%e ", 0.0);
				 }
				 com->fprintf(onlineMatrix[iPodBasis],"\n ");
			 }
		 }
	 }

	 delete [] onlineMatrixFile[0];
	 delete [] onlineMatrixFile[1];
	 
	 // ----------------------
	 // output in full mesh
	 // ----------------------

	 if (ioData->output.transient.bMatrixFull[0] != 0 && ioData->output.transient.aMatrixFull[0] != 0) {
		 onlineMatrixFile[0] = new char[sp + strlen(ioData->output.transient.bMatrixFull)+1];
		 onlineMatrixFile[1] = new char[sp + strlen(ioData->output.transient.aMatrixFull)+1];
		 sprintf(onlineMatrixFile[0], "%s%s", ioData->output.transient.prefix, ioData->output.transient.bMatrixFull);
		 sprintf(onlineMatrixFile[1], "%s%s", ioData->output.transient.prefix, ioData->output.transient.aMatrixFull);
		 int sampleNodeNum;

		 for (int iPodBasis = 0; iPodBasis < 2; ++iPodBasis) {	// always write out 2 files (need A,B)

			 onlineMatrix[iPodBasis] = fopen(onlineMatrixFile[iPodBasis], "wt");
			 com->fprintf(onlineMatrix[iPodBasis],"Vector abMatrix under load for FluidNodes\n");
			 com->fprintf(onlineMatrix[iPodBasis],"%d\n", numTotNodes);
			 com->fprintf(onlineMatrix[iPodBasis],"%d\n", nPod[1]);	// # rows in A and B

				// dummy loop needed by POD

				 for (int iFullNode = 0; iFullNode < numTotNodes; ++iFullNode) {
					 for (int iDim = 0; iDim < dim; ++iDim) { 
						 com->fprintf(onlineMatrix[iPodBasis],"%e ", 0.0);
					 }
					 com->fprintf(onlineMatrix[iPodBasis],"\n");
				 }

			 for (int iPod = 0; iPod < nPod[1]; ++iPod) {	// # rows in A and B
				 com->fprintf(onlineMatrix[iPodBasis],"%d\n", iPod);
				 for (int iFullNode = 0; iFullNode < numTotNodes; ++iFullNode) {
					 isSampleNode = false;
					 // determine if it is a sample node
					 for (int iSampleNode = 0; iSampleNode < nSampleNodes; ++iSampleNode) {
						 if (globalNodeSet[iSampleNode] == iFullNode){
							 isSampleNode = true;
							 sampleNodeNum = iSampleNode;
							 break;
						 }
					 }
					 for (int iDim = 0; iDim < dim; ++iDim) { 
						 if (isSampleNode)
							 com->fprintf(onlineMatrix[iPodBasis],"%e ",
									 onlineMatrices[iPodBasis][sampleNodeNum * dim + iDim][iPod]);
						 else // must output zeros if it is not a sample node
							 com->fprintf(onlineMatrix[iPodBasis],"%e ", 0.0);
					 }
					 com->fprintf(onlineMatrix[iPodBasis],"\n ");
				 }
			 }
		 }
	 }

	 delete [] onlineMatrixFile[0];
	 delete [] onlineMatrixFile[1];
}

template<int dim>
void GappyOffline<dim>::outputStateReduced() {

   com->fprintf(stderr," ... Printing POD state and initial condition in reduced mesh coordinates ...\n");
   int sp = strlen(ioData->output.transient.prefix);
   char *outPodStateFile= new char[sp + strlen(ioData->output.transient.podStateRed)+1];
   char *outInitialConditionFile= new char[sp + strlen(ioData->output.restart.solutions)+1];
	 sprintf(outPodStateFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podStateRed);
	 sprintf(outInitialConditionFile, "%s%s", ioData->output.restart.prefix, ioData->output.restart.solutions);

	 FILE *outPodState = fopen(outPodStateFile, "wt");
	 FILE *outInitialCondition = fopen(outInitialConditionFile, "wt");

	 com->fprintf(outInitialCondition,"Vector InitialCondition under load for FluidNodesRed\n");
	 com->fprintf(outInitialCondition,"%d\n", reducedNodeCount);
	 outputReducedSVec(initialCondition,outInitialCondition,0);

	 // note: need to output the first pod basis vector twice

	 com->fprintf(outPodState,"Vector PodState under load for FluidNodesRed\n");
	 com->fprintf(outPodState,"%d\n", reducedNodeCount);
	 outputReducedSVec(podState[0],outPodState,nPodState);	// dummy output

	 for (int iPod = 0; iPod < nPodState; ++iPod) {	// # rows in A and B
		 outputReducedSVec(podState[iPod],outPodState,iPod);
	 }

	 delete [] outInitialConditionFile;
	 delete [] outPodStateFile;
}

template<int dim>
void GappyOffline<dim>::outputReducedSVec(const DistSVec<double,dim> &distSVec, FILE* outFile , int iVector) {
	 com->fprintf(outFile,"%d\n", iVector);
	 for (int i = 0; i < nSampleNodes; ++i) {
		 // save the reduced node number for the sample node
		 for (int j = 0; j < nodes[i].size(); ++j) {
			 for (int iDim = 0; iDim < dim; ++iDim) { 
				 int currentGlobalNode = nodes[i][j];
				 com->fprintf(outFile,"%e ", distSVec[currentGlobalNode][iDim]);
			 }
			 com->fprintf(outFile,"\n ");
		 }
	 }
}



//template<int dim>
//void GappyOffline<dim>::outputMinimalTopFile() {
//
//
//	std::vector <int> minimalNodes;	// nodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
//	
//
//
//}
//
//template<int dim>
//void GappyOffline<dim>::makeUnique( std::vector<int> minimalNodes,std::vector <int> * nodeOrEle, int length = -1) {
//	for (int i = 0; i < nSampleNodes; ++i) {
//	 // save the reduced node number for the sample node
//	 for (int j = 0; j < nodes[i].size(); ++j) {
//		 minimalNodes.push_back(nodes[i][j]);
//	 }
//	}
//
//	sort(minimalNodes.begin(), minimalNodes.end());	// sort: puts in order
//	it = unique(minimalNodes.begin(), minimalNodes.end()); // remove duplicate consecutive elements (reason for sort)
//	minimalNodes.resize(it - minimalNodes.begin());	// remove extra entries
//}
template<int dim>
void GappyOffline<dim>::outputWallDistanceReduced() {

   com->fprintf(stderr," ... Printing wall distance for reduced mesh ...\n");

	// load in wall distance

	 DistVec<double> *d2wall = geoState->getd2wall();

   int sp = strlen(ioData->output.transient.prefix);
   char *outWallDistFile= new char[sp + strlen(ioData->output.transient.wallDistanceRed)+1];
	 sprintf(outWallDistFile, "%s%s", ioData->output.transient.prefix,
			 ioData->output.transient.wallDistanceRed);

	 FILE *outWallDist = fopen(outWallDistFile, "wt");

	 // note: need to output the first pod basis vector twice

	 com->fprintf(outWallDist,"Scalar walldist under load for FluidNodesRed\n");
	 com->fprintf(outWallDist,"%d\n", reducedNodeCount);
	 outputReducedVec(*d2wall,outWallDist,0);

	 delete [] outWallDistFile;
}

template<int dim>
void GappyOffline<dim>::outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector) {
	 com->fprintf(outFile,"%d\n", iVector);
	 for (int i = 0; i < nSampleNodes; ++i) {
		 // save the reduced node number for the sample node
		 for (int j = 0; j < nodes[i].size(); ++j) {
			 int currentGlobalNode = nodes[i][j];
			 com->fprintf(outFile,"%e ", distVec[currentGlobalNode]);
			 com->fprintf(outFile,"\n ");
		 }
	 }
}
