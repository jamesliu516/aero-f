#ifndef _GAPPY_H_
#define _GAPPY_H_

class IoData;
class Domain;
#include <Communicator.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <ParallelRom.h> // KTC

		// needed? (from Modal.h)

		// #ifdef TYPE_PREC
		// #define PreScalar TYPE_PREC
		// #else
		// #define PreScalar double
		// #endif

		// #include <SpaceOperator.h>
		// #include <KspSolver.h>
		// #include <KspPrec.h>
		// #include <TsInput.h>
		// #include <TsOutput.h>
		// #include <TsRestart.h>
		// 
		// #ifdef DO_MODAL
		//   #include <arpack++/include/ardnsmat.h>
		//   #include <arpack++/include/ardssym.h>
		// #endif
		// 
		// #include <complex.h>
		// typedef complex<double> bcomp;


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

	int nSampleNodes;	// number of parent sample nodes
	void newNeighbors();	// add unique neighbors of a node

	Domain &domain;
	Communicator *com;
	IoData *ioData;	
	TsInput *tInput;

	const int Residual;	// refer to Residual as 0
	const int Jacobian;
	int nPod [2];	// nPod[0] = nPodRes, nPod[1] = nPodJac
	int nPodMax;
	int * cpuSet, * locSubSet, * locNodeSet, * globalNodeSet;	// info for master sample nodes

	int nPodBasis;	// either 1 or 2
	ArrayVecDist<dim> pod;	// pod bases for Residual and Jacobian
	SetOfVec podRes, podJac;
	ArrayVecDist<dim> podHat;	// restricted pod bases for Residual and Jacobian
	SetOfVec podHatRes, podHatJac;
	ArrayVecDist<dim> error;	// restricted pod bases for Residual and Jacobian
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
	
public:
	GappyOffline(Communicator *, IoData &, Domain &, TsInput *);
	~GappyOffline();
	void buildGappy();	// build all offline info (do everything)

};
#include "GappyOffline.C"
#endif
