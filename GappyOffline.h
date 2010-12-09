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

	std::vector< SubDomainData<dim> > &operator [](int i) {return a[i];};
	void resize (int size) {a[0].resize(size); a[1].resize(size); };	// within constructor, specify maximum size
	void resize (int *size) {a[0].resize(size[0]); a[1].resize(size[1]); };	// within constructor, specify maximum size
};

template <typename Scalar, int size>
class StaticArray {	//used for the value of a map
	private:
		Scalar a [size];
	public:
		// use a to access pointer to vecset, use [] to access the vecset itself
		StaticArray() { for (int i=0; i<size; ++i) a[i] = static_cast<Scalar>(0);}
		StaticArray(Scalar b [size]){ for (int i=0; i<size; ++i) a[i] = b[i];}
		StaticArray(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i];}
		StaticArray& operator=(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i];}
		Scalar& operator[] (int i){ return a[i];}
		const Scalar &operator[] (int i) const{ return a[i];}
};

template <int dim>
class GappyOffline {

private:

	typedef VecSet< DistSVec<double,dim> > SetOfVec;
	bool debugging; 	// debugging flag

	std::vector< ParallelRom<dim> *> parallelRom;	// object for all parallel operations
	void setUpPodBases();	// read in POD bases
	void buildGappyMesh();	// build reduced mesh offline
	void buildGappyMatrices();	// build matrices A and B

	int nSampleNodes;	// number of parent sample nodes
	void newNeighbors();	// add unique neighbors of a node

	Domain &domain;
	Communicator *com;
	IoData *ioData;	
	DistGeoState *geoState;
	DistSVec<double,3> &X;

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

	int nPodState;	// TODO want to treat all of these vectors
	SetOfVec podState;	// need to put podState in reduced coordinates
	DistSVec<double,dim> initialCondition;
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
	std::vector <int> *nodes;	// nodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
	std::vector <double> *(nodesXYZ [3]);	// nodesXYZ[iXYZ][iSampleNode][iNode] is the iXYZ coordinate of the iNode in the iSampleNode island
	std::vector <int> *elements;		// elements[iSampleNode][iEle] is the global element number of the iEle element in the iSampleNode island 
	std::vector <int> *(elemToNode [4]);	// elemToNode[iNode][iSampleNode][iEle] is the global node number of the iNode attached to the iEle element of the iSampleNode island 
	std::vector< int > *(bcFaces [2][3]);	// boundary faces. bcfaces[iSign][whichNode][BCtype][iFace] returns the global node number of whichNode on the iFace face corresponding to iSign/BCtype. iSign = 0 if the BC definition is negative, and iSign = 1 if positive. BCtype can be found in BcDefs.h
	std::vector< int > *(bcIsland [2]);	// bcIsland[iSign][BCtype][iFace] returns the island to which the face belongs
	int * sampleToReducedNodeNumbering;

	std::map<int, StaticArray <double, 3> > nodesXYZmap;	// key: global node #, values: x, y, z
	std::map <int, StaticArray <int, 4> > elemToNodeMap;	// key: global elem #, values: global node #s
	std::map <int, std::string > boundaryConditionsMap;	// mapping between BC numbers in BcDef.h and Sower's identification
	int reducedNodeCount;	// total number of nodes in the reduced mesh
	int numTotNodes;	// total number of nodes in the full mesh

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
	template<typename Scalar> void communicateMesh( std::vector <Scalar> *nodeOrEle , int arraySize);
	void communicateAll();
	void defineMaps();
	void communicateBCFaces();
	void makeUnique( std::vector <int>  *nodeOrEle, int length);
	void outputTopFile();

	// A and B matrices functions

	// pseudo-inverse functions
	double **(podHatPseudoInv [2]);	// each dimension: (nSampleNode*dim) x nPod[i]
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
	void outputStateReduced();
	void outputWallDistanceReduced();
	void outputReducedSVec(const DistSVec<double,dim> &distSVec, FILE* outFile , int iVector);
	void outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector);

public:
	GappyOffline(Communicator *, IoData &, Domain &, DistGeoState *);
	~GappyOffline();
	void buildGappy();	// build all offline info (do everything)

};
#include "GappyOffline.C"
#endif
