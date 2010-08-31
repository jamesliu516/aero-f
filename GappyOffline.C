#include <GappyOffline.h>

template<int dim>
GappyOffline<dim>::GappyOffline(Communicator *_com, IoData &_ioData, Domain &dom, TsInput *_tInput) : 
	domain(dom), 	com(_com), ioData(&_ioData), tInput(_tInput),
	podRes(0, dom.getNodeDistInfo() ), podJac(0, dom.getNodeDistInfo() ),
	podHatRes(0, dom.getNodeDistInfo() ),	podHatJac(0, dom.getNodeDistInfo() ),	
	errorRes(0, dom.getNodeDistInfo() ),	errorJac(0, dom.getNodeDistInfo() ),	
	handledNodes(0), nPodBasis(2),
	debugging(1),
	// distribution info
	numLocSub(dom.getNumLocSub()), nTotCpus(_com->size()), thisCPU(_com->cpuNum()),
	nodeDistInfo(dom.getNodeDistInfo()), subD(dom.getSubDomain()),parallelRom(2,ParallelRom<dim>(dom,_com)),
	residual(0), jacobian(1)
{
	// initialize vectors to point to the approprite bases

	pod.a[0] = &podRes;	// make pod point to res and jac 
	pod.a[1] = &podJac; 
	podHat.a[0] = &podHatRes;	// make pod point to res and jac 
	podHat.a[1] = &podHatJac; 
	error.a[0] = &errorRes;	// make pod point to res and jac 
	error.a[1] = &errorJac; 

	handledVectors[0] = 0;	// have not yet handled any vectors
	handledVectors[1] = 0;	// have not yet handled any vectors
	// scalapack least-squares info

		
  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		parallelRom[iPodBasis].parallelLSMultiRHSInit(podHat[iPodBasis], error[iPodBasis]);
	} 
	
}
template<int dim>
GappyOffline<dim>::~GappyOffline() 
{
}
	//---------------------------------------------------------------------------------------

template<int dim>
void GappyOffline<dim>::buildGappy() {

	// determine whether one or two pod bases are used

	nPod[0] = ioData->Rob.numROB;	
	nPod[1] = ioData->Rob2.numROB;

	nPodMax = max(nPod[0],nPod[1]);	// compute maximum nPod

	if (nPod[1] == 0) {	// only need one basis (shared between jacobian and residual)
		nPodBasis = 1;
		pod.a[1] = &podRes;
		podHat.a[1] = &podHatRes;	// make pod point to res and jac
		error.a[1] = &errorRes;	
	}
	else if (nPod[0] ==0) {
		nPodBasis = 1;
		pod.a[0] = &podJac;
		podHat.a[0] = &podHatJac;
		error.a[0] = &errorJac;	
	}

	com->fprintf(stderr, " ... Reading POD bases for Gappy POD contruction\n");
	int nSampleNodes = ceil(double(nPodMax)/double(dim)); // this will give interpolation or the smallest possible least squares
	// require nSampleNodes * dim >= max(nPod[0],nPod[1]); nSampleNodes >= ceil(double(max(nPod[0],nPod[1]))/double(dim))

	for (int i = 0 ; i < nPodBasis ; ++i)// only do for number of required bases
		pod[i].resize(nPod[i]);

	//	read in both Pod bases
	//	tInput->podFile: file containing file names of multiple bases
  domain.readMultiPodBasis(tInput->podFile, pod.a, nPod, nPodBasis);

	// compute the reduced mesh used by Gappy POD

   //buildGappyMesh(); 

	// compute matries A and B required online

   //buildGappyMatrices();

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

	// STOPPING HERE! David doesn't care about anything after this point
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
		nSampleNodes = ceil(double(nPodMax)/double(dim)); 
		com->fprintf(stderr,"Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d",nSampleNodesOld,nSampleNodes);
	}

	nRhsMax = ceil(double(nPodMax)/double(nSampleNodes)); // nSampleNodes * nRhsMax >= max(nPod[0],nPod[1])

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
	// 	Global: locPodHat for maximum entry, globally summed cpuSet etc.
	//===============================================
		
	int thisCPU = com->cpuNum();
	double globalMaxNorm = myMaxNorm;
	com->barrier();
	com->globalMax(1, &globalMaxNorm);  // find the maximum value over all cpus

	if (myMaxNorm == globalMaxNorm) {  // if this CPU has the maximum value

		if (debugging){
		 fprintf(stderr, "CPU %d has myMaxNorm with locNode %d, locSub %d \n",thisCPU,locNode,locSub);
		}

		// save the global subdomain and local node indices (sum at the very end of
		// algorithm)

		cpuSet[handledNodes] = thisCPU;
		locSubSet[handledNodes] = locSub;
		locNodeSet[handledNodes] = locNode;
		globalNodeSet[handledNodes] = globalNode;

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
	
	com->globalSum(1, cpuSet + handledNodes);
	com->globalSum(1, locSubSet + handledNodes);
	com->globalSum(1, locNodeSet + handledNodes);
	com->globalSum(1, globalNodeSet + handledNodes);

	++handledNodes;
}

template<int dim>
void GappyOffline<dim>::greedy(int greedyIt) {

	// Differences for 1st iteration compared with other greedy iterations:
	// 1) no least squares problem is solved (just take the maximum entry)
	// 2) look at the inlet face to ensure boundary condition is handled

	double myMaxNorm;
	int locSub, locNode, globalNode;  // temporary global subdomain and local node

	bool doLeastSquares = true;
	bool onlyInletBC = false;

	if (greedyIt == 0) {	
		doLeastSquares = false;// don't do least squares if first iteration
		onlyInletBC = true;// first iteration, only consider inlet BC
	}
		

	// locError[ nPodBasis ][ nRhs[iPodBasis] ][ locNodeNum ][ iDim ]
	// 																				 [  SUBDOMAIN INFO    ]

	locError.resize(nRhs);	// create locError entity of maximal size nRhsMax

	if (doLeastSquares) {
		leastSquaresReconstruction();		// solve the least squares reconstruction
	}

	for (int iFillNode ; iFillNode < nodesToHandle[0]; ++iFillNode) {	// fill up the appropriate number of nodes

		// initialize parameters
		myMaxNorm = 0.0;	// initial maximum is zero
		locSub = -1; locNode = -1; globalNode = -1;	// where the maximum is located

		// determine number of rhs for each
		// typically have nRhs = nRhsMax. exception: there are less than nRhsMax remaining vectors

		for (int iPodBasis = 0; iPodBasis  < nPodBasis; ++iPodBasis)	
			nRhs[iPodBasis] = min(nRhsMax, max(nPod[iPodBasis] - handledVectors[iPodBasis], 0));	

		// loop over nodes, and add to set if it is the maximum
		// subdomains -> nodes
		for (int iSub = 0; iSub < numLocSub; ++iSub) {

			// get subdomain info for all RHS

			getSubDomainError(iSub);

			// find maximum error on the subdomain
			
			subDFindMaxError(iSub, onlyInletBC, myMaxNorm, locSub, locNode, globalNode);
			
		}

		// find global subdomain number and local node number for node with maximum norm

		findMaxAndFillPodHat(myMaxNorm, locSub, locNode, globalNode);

	}

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)	
		handledVectors[iPodBasis] +=nRhs [iPodBasis];
}


template<int dim>
void GappyOffline<dim>::initializeGappyLeastSquares() {

	// initialize least squares problems

	for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
		parallelRom[iPodBasis].parallelLSMultiRHSInit(podHat[iPodBasis], error[iPodBasis]); 
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
	// 	Local: nodeError

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
		lsCoeff[iPodBasis] = new double * [nRhs[iPodBasis]];
		for (int iPod = 0; iPod < nRhs[iPodBasis]; ++iPod)  {
			// temporarily fill the error vector with the RHS (solving nRhs[iPodBasis] problems)
			error[iPodBasis][iPod] = podHat[iPodBasis][handledVectors[iPodBasis] + iPod];	// NOTE: PODHAT
      lsCoeff[iPodBasis][nRhs[iPodBasis]] = new double [handledVectors[iPodBasis]];
		}

		parallelLSMultiRHSGap(iPodBasis,lsCoeff[iPodBasis]);	

		// compute reconstruction error

		for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs) { 
			error[iPodBasis][iRhs] = pod[iPodBasis][handledVectors[iPodBasis] + iRhs];	// NOTE: POD
			for (int jPod = 0; jPod < handledVectors[iPodBasis]; ++jPod) {
				error[iPodBasis][iRhs] -= error[iPodBasis][jPod] * lsCoeff[iPodBasis][iRhs][jPod];
			}
		}
	}
}

template<int dim>
void GappyOffline<dim>::subDFindMaxError(int iSub, bool onlyInletBC, double &myMaxNorm, int &locSub, int &locNode, int &globalNode) {

	// PURPOSE: Search the iSub subdomain for possible maximum error
	// Inputs:
	// 	Passed: iSub, onlyInletBC, myMaxNorm, locSub, locNode, globalNode
	// Outputs:
	// 	Passed: (all can change) myMaxNorm, locSub, locNode, globalNode

	bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
	int nLocNodes = nodeDistInfo.subSize(iSub);	// number of nodes in this subdomain
	double nodeError; // the error at the node

	if (onlyInletBC) {	// consider only nodes on inlet BC
		FaceSet& currentFaces = subD[iSub]->getFaces();
		for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {	
			// only consider inlet boundary conditions
			if (currentFaces[iFace].getCode() == BC_INLET_MOVING || currentFaces[iFace].getCode() == BC_INLET_FIXED) {
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
		parallelRom[iPodBasis].parallelLSMultiRHS(podHat[iPodBasis],error[iPodBasis], handledVectors[iPodBasis], nRhs[iPodBasis], lsCoeff);
}
