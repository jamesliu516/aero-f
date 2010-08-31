#ifndef _GAPPY_H_
#define _GAPPY_H_

class IoData;
class Domain;
#include <Communicator.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <ParallelRom.h>

template <int dim>
class ArrayVecDist {
	public:
		// use a to access pointer to vecset, use [] to access the vecset itself
		VecSet< DistSVec<double,dim> > *(a[2]);	// B[0] = *(a[0])
		VecSet< DistSVec<double,dim> >& operator [](int i) { return *(a[i]); };
};

template <int dim>
struct SubDomainData { 
	// need this struct because we are creating vectors of them (unknown number of
	// entries
	double (*a) [dim];
	double * operator [](int i) {return a[i];};
	SubDomainData& operator=(double (*b)[dim] ) {
		a = b;
		return *this;
	};
};

template <int dim>
struct VecSubDomainData { 
	// need this struct because the error vector is a 4-D array and complicated
	// 
	// locError[ nPodBasis ][ nRhs[iPodBasis] ][ locNodeNum ][ iDim ]

	std::vector< SubDomainData<dim> > a [2];	// array of two vectors, each vector is nRhs long

	std::vector< SubDomainData<dim> > operator [](int i) {return a[i];};
	void resize (int size) {a[0].resize(size); a[1].resize(size); };	// within constructor, specify maximum size
	void resize (int *size) {a[0].resize(size[0]); a[1].resize(size[1]); };	// within constructor, specify maximum size
};

template <int size>
class StaticArray {	//used for the value of a map
	private:
		double a [size];
	public:
		// use a to access pointer to vecset, use [] to access the vecset itself
		StaticArray() { for (int i=0; i<size; ++i) a[i] = 0.0;}
		StaticArray(double b [size]){ for (int i=0; i<size; ++i) a[i] = b[i];}
		StaticArray(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i];}
		StaticArray& operator=(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i];}
		double& operator[] (int i){ return a[i];}
		const double &operator[] (int i) const{ return a[i];}
};

template <int dim>
class GappyOffline {

private:

	typedef VecSet< DistSVec<double,dim> > SetOfVec;
	bool debugging; 	// debugging flag

	std::vector< ParallelRom<dim> > parallelRom;	// object for all parallel operations
	void buildGappyMesh();	// build reduced mesh offline
	void buildGappyMatrices();	// build matrices A and B

	int nSampleNodes;	// number of parent sample nodes
	void newNeighbors();	// add unique neighbors of a node

	Domain &domain;
	Communicator *com;
	IoData *ioData;	
	TsInput *tInput;
	DistGeoState *geoState;
	DistSVec<double,3> X;

	const int residual;	// refer to residual as 0
	const int jacobian;
	int nPod [2];	// nPod[0] = nPodRes, nPod[1] = nPodJac
	int nPodMax;
	int * cpuSet, * locSubSet, * locNodeSet, * globalNodeSet;	// info for master sample nodes

	int nPodBasis;	// either 1 or 2
	ArrayVecDist<dim> pod;	// pod bases for residual and jacobian
	SetOfVec podRes, podJac;
	ArrayVecDist<dim> podHat;	// restricted pod bases for residual and jacobian
	SetOfVec podHatRes, podHatJac;
	ArrayVecDist<dim> error;	// restricted pod bases for residual and jacobian
	SetOfVec errorRes, errorJac;

	// greedy data

	int nRhsMax, nGreedyIt, handledNodes;
	int handledVectors [2];
	int nRhs [2];	// nRhs at a given greedy iteration
	int *nodesToHandle;	// how many nodes are handled at each greedy iteration
	VecSubDomainData<dim> locError;

	void parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff);

	// greedy functions

	void greedy(int greedyIt);
	void computeGreedyIterationInfo();
	void computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError);
	void findMaxAndFillPodHat(double myMaxNorm, int locSub, int locNode, int globalNode);
	void makeNodeMaxIfUnique(double nodeError, double &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int &globalNode);
	void getSubDomainError(int iSub);	// computes locError for iSub subdomain
	void leastSquaresReconstruction();
	void subDFindMaxError(int iSub, bool onlyInletBC, double &myMaxNorm, int &locSub, int &locNode, int &globalNode);	// computes locError for iSub subdomain

	// least squares parameters
	
	void initializeGappyLeastSquares();

	// distribution info

	int numLocSub, nTotCpus, thisCPU;
	DistInfo &nodeDistInfo;
	SubDomain** subD; 
	
	// mesh construction
	// each of these arrays has nSampleNodes elements
	std::vector <int> *nodes;	// nodes[iObsNode][iNode] is the global node number of the iNode in the iObsNode island 
	std::vector <int> *elements;		// elements[iObsNode][iEle] is the global element number of the iNode in the iObsNode island 
	std::vector< int > *(bcFaces [2][3]);	// boundary faces. faces[iSign][whichNode][BCtype][iFace] returns the global node number of whichNode on the iFace face corresponding to iSign/BCtype. iSign = 0 if the BC definition is negative, and iSign = 1 if positive. BCtype can be found in BcDefs.h
	std::vector< int > *(bcIsland [2]);	// bcIsland[iSign][BCtype][iFace] returns the island to which the face belongs
	double *sampleToReducedNodeNumbering;	// nodes[iObsNode][iNode] is the global node number of the iNode in the iObsNode island 

	std::map<int, StaticArray <3> > nodesXYZ;	// key: global node #, values: x, y, z
	std::map <int, StaticArray <4> > elemToNode;	// key: global elem #, values: global node #s
	std::map <int, string > boundaryConditions;	// mapping between BC numbers in BcDef.h and Sower's identification
	int reducedNodeCount;	// total number of nodes in the reduced mesh

		// KTC!!! then, when outputting the TOP file, need another key that maps global
		// node # to reduced mesh node #... this mapping will be different for
		// each island!

		// also create a pointer of vectors to handle the faces that might be
		// boundary conditions

	void computeXYZ(int iSub, int iLocNode, double *xyz);


	void addTwoNodeLayers();
	void computeBCFaces();
	void checkFaceInMesh(FaceSet& currentFaces, int iFace, int iSub, int *locToGlobNodeMap , bool &faceInMesh, int &whichIsland);
	void addNeighbors(int iSampleNodes, int startingNodeWithNeigh);
	void communicateMesh( std::vector <int> *nodeOrEle , int arraySize);
	void communicateBCFaces();
	void makeUnique( std::vector <int>  *nodeOrEle, int length);
	void outputTopFile();

	// A and B matrices functions

	// pseudo-inverse functions
	double **(podHatPseudoInv [2]);	// dimension: (nSampleNode*dim) x nPod[i]
	void computePseudoInverse(int iPodBasis);
	void computePseudoInverseRHS();	// computes the RHS matrix pseudoInvRhs
	SetOfVec pseudoInvRhs;

	// podTpod
	double **podTpod;	// stores phiJ^TphiR
	void computePodTPod();
	double **(onlineMatrices [2]);	// dimension: (nSampleNode*dim) x nPod[1]
		// onlineMatrices[0] is related to the residual: 
		// 		pod[1]^Tpod[0] * podHatPseudoInv[0]^T
		// onlineMatrices[1] is related to the jacobian: 
		// 		podHatPseudoInv[1]^T
	void assembleOnlineMatrices();
	void outputOnlineMatrices();

public:
	GappyOffline(Communicator *, IoData &, Domain &, TsInput *, DistGeoState *);
	~GappyOffline();
	void buildGappy();	// build all offline info (do everything)

};
#include "GappyOffline.C"
#endif
