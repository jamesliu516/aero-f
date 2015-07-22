#ifndef _PARALLEL_ROM_H_
#define _PARALLEL_ROM_H_

//#include <Elem.h>	// use ElemSet
class Domain;
class Communicator;
class SubDomain;
template<typename Scalar> class GenFullM;
typedef GenFullM<double> FullM;

template <int dim>
class ParallelRom {

	private:
	Domain &domain;
	Communicator *com; 
	SubDomain **subDomain;
	int nTotCpus;
	int thisCPU;

	int *cpuMasterNodes;	//number of master nodes on cpu for given mesh decomposition
	int *cpuNodes;	//number of master nodes needed by cpu for scalapack block cyclic decomposition
	int maxCpuBlocks;	//	number of maximum blocks per cpu for block cyclic decomposition
	int *locSendReceive;

	int rowsPerBlock;
	int nodesPerBlock;

	// least squares
	int desc_a [9], desc_b [9];

	// private functions
	void scalapackCpuDecomp(const int nCol);	// transfer data from VecSet<DistSVec> to scalapack format
	template<class VecContainer> void transferData (VecContainer &snaps, double* subMat, int nSnaps);
	template<class VecContainer> void transferDataBack(double *U, VecContainer &Utrue , int nSnaps);
	void transferDataBackLS(double *subMatB, int n, double **lsSol, int nRhs, int subMatLLD, bool); 
	void setTransfer();

	public:
	ParallelRom(Domain &, Communicator *);
	~ParallelRom();

	// Parallel operations
	template<class VecContainer1, class VecContainer2> void parallelSVD(VecContainer1 &snaps,
			VecContainer2 &Utrue, double *S, FullM &Vtrue, int nSnaps, bool computeV=true);
	template<class VecContainer1, class VecContainer2> void parallelLSMultiRHSInit(const VecContainer1
			&A, const VecContainer2 &B, const int nA=0); 	// initialization
	template<class VecContainer1, class VecContainer2> void parallelLSMultiRHS(const VecContainer1 &A,
			const VecContainer2 &B, int n, int nRhs, double **lsSol, bool=true); // least-squares
			// via QR with multiple RHS
};
#endif
