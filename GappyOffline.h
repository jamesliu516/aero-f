#ifndef _GAPPY_H_
#define _GAPPY_H_

class IoData;
class Domain;
#include <Communicator.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <ParallelRom.h> // KTC

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

	void parallelLSMultiRHSGap(int iPodBasis, double *lsCoeff);

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
	std::map<int, int [3] > nodesXYZ;	// key: global node #, values: x, y, z
	std::map <int, int [4] > elemToNode;	// key: global elem #, values: global node #s
		// then, when outputting the TOP file, need another key that maps global
		// node # to reduced mesh node #... this mapping will be different for
		// each island!
	void computeXYZ(int iSub, int iLocNode, double *xyz);


	void addTwoNodeLayers();
	void addNeighbors(int iSampleNodes, int startingNodeWithNeigh);
	void communicateMesh( std::vector <int> *nodeOrEle );
	void makeUnique( std::vector <int>  *nodeOrEle, int length);
	void outputTopFile();

public:
	GappyOffline(Communicator *, IoData &, Domain &, TsInput *);
	~GappyOffline();
	void buildGappy();	// build all offline info (do everything)

};
#include "GappyOffline.C"
#endif
