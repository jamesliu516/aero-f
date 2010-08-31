#ifndef _PARALLEL_ROM_H_
#define _PARALLEL_ROM_H_

//#include <Elem.h>	// use ElemSet

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
	void scalapackCpuDecomp(int nCol);	// transfer data from VecSet<DistSVec> to scalapack format
	void transferData(VecSet< DistSVec<double, dim> > &snaps, double* subMat, int nSnaps);
	void transferDataBack(double *U, VecSet< DistSVec<double, dim> > &Utrue , int nSnaps);
	void transferDataBackLS(double *subMatB, int n, double **lsSol, int nRhs, int subMatLLD); 
	void setTransfer();

	public:
	ParallelRom(Domain &, Communicator *);
	~ParallelRom();

	// Parallel operations
	void parallelSVD(VecSet< DistSVec<double, dim> > &snaps, VecSet<DistSVec<double, dim> > &Utrue, double *S, FullM &Vtrue, int nSnaps);
	void parallelLSMultiRHSInit(VecSet< DistSVec<double, dim> > &A, VecSet<DistSVec<double, dim> > &B); 	// initialization
	void parallelLSMultiRHS(VecSet< DistSVec<double, dim> > &A, VecSet<DistSVec<double, dim> > &B, int n, int nRhs, double **lsSol); // least-squares
			// via QR with multiple RHS
};
#include "ParallelRom.C"
#endif
