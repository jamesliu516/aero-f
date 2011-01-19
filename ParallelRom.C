#include <LinkF77.h>
#include <ParallelRom.h>

extern "C"      {

  void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                       double *, double *, int &, double *, int &,
                       double *, const int &, int &);

  void F77NAME(thinsvd)(int &, int &, int &, int &, int&, int&, int&, int&, int&, double *, int &,
                       int &, int &, int &, int &, int &, int &, double *U, double *S, double *V,
                       int &, double *, int &);
  void F77NAME(lworksizesvd)(int &, int &, int &, int &, int&, int&, int&, int&, int&, int &);
	void F77NAME(globalmatrices)(int &, int &, int &, int &, int &, int&, int &, int &,int &, int *, int *);	

	void F77NAME(pdgels)(char &, int &, int &, int &, double *, int &, int &, int *,
											 double *, int &, int &, int *, double *, int &, int & );
};

template<int dim> 
ParallelRom<dim>::ParallelRom(Domain & _domain, Communicator *_com) : 
domain(_domain), com(_com), subDomain(_domain.getSubDomain()), maxCpuBlocks(0), nTotCpus(_com->size()), thisCPU(_com->cpuNum())
{
	cpuNodes = new int[nTotCpus];
	cpuMasterNodes = new int[nTotCpus];
	locSendReceive = new int[nTotCpus];
	for (int i = 0; i < 9; ++i) {
		desc_a[i] = -1;
		desc_b[i] = -1;
	}
}

template<int dim> 
ParallelRom<dim>::~ParallelRom() 
{

	delete [] cpuNodes;
	delete [] cpuMasterNodes;
	delete [] locSendReceive;
	
}
template<int dim> 
void ParallelRom<dim>::scalapackCpuDecomp(int nCol) {

	//===============================
	// PURPOSE: determine cpu decomposition for ScaLAPACK of a matrix which
	// is described as VecSet< DistSVec< >>
	// INPUTS: nCol
	// OUTPUTS: rowsPerBlock, cpuMasterNodes, cpuNodes, maxCpuBlocks
	//===============================

#ifdef DO_SCALAPACK

 // specify the block size (in terms of nodes), and require nodesPerBlock > nCol 

 nodesPerBlock = int(ceil(double(nCol)/double(dim)));	// ensures that square blocks can be used with one processor containing all the columns for its nodes 
 rowsPerBlock = nodesPerBlock * dim;	// number of entries per block (there are dim entries per node)


 //===============================
 // set transfer parameters
 //===============================

 setTransfer();	// determines cpuMasterNodes,cpuNodes,maxCpuBlocks from nodesPerBlock

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH SCALAPACK and DO_SCALAPACK Flag\n");
 exit(-1);
#endif

}

//----------------------------------------------------------------------------------

template<int dim>
void ParallelRom<dim>::parallelSVD(VecSet< DistSVec<double, dim> > &snaps, VecSet<DistSVec<double, dim> > &Utrue, double *S, FullM &Vtrue, int nSnaps) {

#ifdef DO_SCALAPACK
 	int numLocSub = domain.getNumLocSub();

 // specify the block size (in terms of nodes)
 // blockfactor > nSnaps => the right singular vector matrix V is computed by cpu 0

 nodesPerBlock = int(floor(2.0*nSnaps/dim));


 //set the transfer parameters

 setTransfer();

 //create the matrix of snapshots with the right distribution
 double *subMat = new double[nodesPerBlock*dim*maxCpuBlocks*nSnaps];
 for (int i=0; i < nodesPerBlock*dim*maxCpuBlocks*nSnaps; ++i)
   subMat[i] =0.0;
 
 //transfer the extra nodes where needed and fill subMat
 transferData(snaps,subMat,nSnaps);

 com->barrier();

 // Allocate for svd call
 // allocate for eigenvectors and eigenvalues
 int locLLD = dim * cpuNodes[thisCPU];
 int maxLLD = locLLD;
 com->globalMax(1, &maxLLD);
 int maxLocC = nSnaps;
 int locLLD_V = 1;
 if (thisCPU == 0)
   locLLD_V = nSnaps;
 
 double *U = new double[locLLD*maxLocC];
 double *V = new double[locLLD_V*maxLocC];
 for (int iSnaps=0; iSnaps<nSnaps; ++iSnaps)
   S[iSnaps] = 0.0;
 int nprow = nTotCpus;
 int npcol = 1;
 int globNumRows = dim * snaps[0].nonOverlapSize();
 // call svd
 int lwork = 0;
 int myrow, mycol, ictxt, info;
 int rowIndex = 1;
 int colIndex = 1;
 nodesPerBlock *= dim;
 
 F77NAME(lworksizesvd)(ictxt, nprow, npcol, myrow, mycol, globNumRows, nSnaps,
                  nodesPerBlock, nodesPerBlock, lwork);

 com->barrier();
 lwork *= 3;//debug
 double *work = new double[lwork + 2];
 work[0] = -99.0;

 com->barrier();
	// scalapack routine!
 F77NAME(thinsvd)(ictxt, nprow, npcol, myrow, mycol, globNumRows, nSnaps,
                  nodesPerBlock, nodesPerBlock, subMat, maxLLD, locLLD, nSnaps,
                  locLLD_V, thisCPU, rowIndex, colIndex, U, S, V, lwork, work,
                  info);

 com->barrier();
 
 for (int i = 0; i < nSnaps; i++){
   for (int j = 0; j < nSnaps; j++) {
     Vtrue[i][j] = 0.0;
     if (thisCPU==0)
       Vtrue[i][j] = V[nSnaps*i+j];
   }
 }
 delete[] V;

 com->globalSum(nSnaps*nSnaps, Vtrue.data());
 
 //Permute back matrix U
 com->barrier();
 for (int i = 0; i < nSnaps; ++i)
   Utrue[i] = 0.0;
 transferDataBack(U, Utrue , nSnaps);

 com->barrier();
 delete[] work; 
 delete[] U; 

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH SCALAPACK and DO_SCALAPACK Flag\n");
 exit(-1);
#endif
}

//----------------------------------------------------------------------------------

template<int dim>
void ParallelRom<dim>::parallelLSMultiRHSInit(VecSet< DistSVec<double, dim> > &A, VecSet<DistSVec<double, dim> > &B) {

	//====================================
	// PURPOSE: initialize the least squares problems for scalapack when
	// matrices are in the form VecSet< DistSVec <> > 
	// INPUTS: A, B OUTPUTS:
	// cpuMasterNodes, cpuNodes, desc_a, desc_b
	//====================================

	//====================================
	// CONSTRAINTS 
	// CTXT_A = CTXT_B (same contexts for both matrices)
	// A, B, and work must have no common elements
	// need A(ia:ia+min(m,n)-1, ja:ja+min(m,n)-1) to be in a single block
	// for m >= n (least squares)
	// 	MB_A = MB_B (number of rows in a block must be equal)
	// 	mod(ia-1, MB_A) = mod(ib-1,MB_B) (row block offset of A must be equal to the row block offset of B)
	// 	process row with first row of A must contain the first row of submatrix B
	//
	//	these are fixed by always setting MB_A = MB_B > = n_A, ia = ib = 0, RSRC_A = RSRC_B 
	//====================================

	//===============================
	// NOTES
	// desc_a contains ictxt (context), m_a, n_a, mb_a, nb_a, rsrc_a, csrc_a, lld_a
	// 	then, from ictxt, can get nprow, npcol, myrow, mycol via BLACS_GRIDINFO
	// Assume you can store an nCol x nCol matrix on a single processor
	// ScaLAPACK's driver routines operate on a (continuous) submatrix.  So make an unnecessarily large matrix (more columns than needed) to fit some of ScaLAPACK's constraints, and end up using only the first nCol columns.
	//===============================

	//====================================
	//	STEPS
	// 	1) figure out decomposition
	// 	2) define process grid (ictxt, nprow, npcol, myrow, mycol by running SL_INIT and BLACS_GRIDINFO)
	// 	3) define description vector for each global matrix
	//====================================

	//===============================
	// Step 1 - depends on # cols in A
	//===============================

	// decomposition information

	int numLocSub = domain.getNumLocSub();
	nTotCpus = com->size(); 
	thisCPU = com->cpuNum(); 

	maxCpuBlocks=0;	// maximum number of blocks per cpu 

	// inputs: nA
	// outputs: rowsPerBlock, cpuMasterNodes, cpuNodes (number of nodes in scalapack block cyclic decomp), maxCpuBlocks
	int nA = A.numVectors();
	scalapackCpuDecomp(nA); 

 //===============================
 // Step 2: define process grid - depends on # cols in A
 //===============================
 
 int locLLD = max(1,dim * cpuNodes[thisCPU]);	// number of rows in domain!
 int maxLLD = locLLD;	// find the maximum over all LLD
 com->globalMax(1, &maxLLD); 
 
 // define process (cpu) grid: nTotCpus x 1 (fit all data from a single node on a processor)
 // NOTE: can have NB_A > N_A! Your block size can be bigger than the number of rows
 
 int nprow = nTotCpus;	// number of processor rows
 int npcol = 1;	// number of processor columns
 int globNumRows = dim * A[0].nonOverlapSize();
	
	//===============================
	// Step 3 - depends on # cols in A and # cols in B
	//===============================

	// input: nprow, npcol, m_a, n_a, n_b, rowblock, colblock, locLLD
	// output: desc_a, desc_b
	
	int nB = B.numVectors();
  F77NAME(globalmatrices)(nprow, npcol, globNumRows, nA, nB, rowsPerBlock, rowsPerBlock,locLLD, thisCPU, desc_a, desc_b);	

}

//----------------------------------

template<int dim>
void ParallelRom<dim>::parallelLSMultiRHS(VecSet< DistSVec<double, dim> > &A,
		VecSet<DistSVec<double, dim> > &B, int n, int nRhs, double **lsSol) {

	// each time
  // 	ScaLAPACK only operates on submatrices; these are the submatrices of interest
	// 	define a new submat with new number of columns and correct data
	//  determine lwork (depends on size of submatrix)
	// 	solve least squares problem with this information
	// 	assume: use all rows!

	//====================================
	// PURPOSE: solve the least squares system(s) with submatrix A and sub-matrices for the RHS B 
	// INPUTS: A (matrix), B (RHS), n (# columns in submatrix), nRhs (# columns in B)
	//		cpuMasterNodes,cpuNodes, desc_a, desc_b 
	// OUTPUTS: lsSol (solutions)
	//====================================

#ifdef DO_SCALAPACK

	//===============================
	// TRANSFER DATA (WHEN NEEDED) TO MAKE IT CONSISTENT WITH SCALAPACK QUANTITIES
	//===============================

	//===============================
	// check for any changes in the dimensions of the least squares problem
	//===============================
	
	if (n != desc_a[3] || nRhs != desc_b[3])
		parallelLSMultiRHSInit(A, B);	// (re)-initialization needed

	// KTC: can make this more elegant if needed
	
	int globNumRows = desc_a[2];	// the third entry has number of global rows. This is m in the submatrix
	int nSubMatA = n;	// only choose number of columns actually used
	int nSubMatB = nRhs;

	//===============================
	// initialize the matrices with the right distribution
	//===============================

	int locLLD = dim * cpuNodes[thisCPU];
	int subMatLLD = locLLD;
	thisCPU = com->cpuNum();
	//int subMatLLD = maxLLD = rowsPerBlock * maxCpuBlocks;	// this is the safe way to do it

 double *subMatA = new double[n * subMatLLD];
 double *subMatB = new double[nRhs * subMatLLD];
 for (int i=0; i < n * subMatLLD; ++i)
   subMatA[i] =0.0;
 for (int i=0; i < nRhs * subMatLLD; ++i)	
   subMatB[i] =0.0;

 //===============================
 // transfer the extra nodes where needed and fill subMat
 // this goes from physical node-based decomposition to the computationally balanced decomposition
 // 	Note that A and B use the same decomposition
 //===============================

 locSendReceive = new int[com->size()];	// how many nodes the current cpu sends to other cpu (not needed after this 

 // input: B, cpuMasterNodes, cpuNodes, nB
 // output: locSendReceive, subMatB
 transferData(B,subMatB,nSubMatB);	// determines locSendReceive(), subMat
 transferData(A,subMatA,nSubMatA);	// determines locSendReceive(), subMat

 com->barrier();

 //===============================
 // compute local workspace array needed for least squares problem
 //===============================


 // always start with the top left corner of the matrix
 int rowIndex = 1;	// IA in PDGELS.f: row index in global array A indicating first row of submatrix
 int colIndex = 1;  // JA in PDGELS.f: col index in global array A indicating first col of submatrix

 // inputs: desc_a, desc_b, rowIndex, colIndex, info for submatrix problem (m, n, nRhs)
 // output: lwork
 // determine size of array

 int lwork = -1;	// size of local workspace array
 double worktmp;
 int info;
 char normalChar = 'N';

 // when lwork = -1, pdgels returns the optimal lworksize in worktmp

 F77NAME(pdgels)(normalChar, globNumRows, n, nRhs, subMatA, rowIndex, colIndex, desc_a, 
		 								subMatB, rowIndex, colIndex, desc_b, &worktmp, lwork, info);

 lwork = static_cast<int>(worktmp);

 com->barrier();
 // allocate memory for local work array
 double *work = new double[lwork];

 com->barrier();

 //===============================
 // solve least squares problem in scalapack
 //===============================

 // output: subMatB has been changed to contain the solution
 F77NAME(pdgels)(normalChar, globNumRows, n, nRhs, subMatA, rowIndex, colIndex,
	desc_a, subMatB, rowIndex, colIndex, desc_b, work, lwork, info);
	//what about rowIndex and colIndex for A and B???

 if (info != 0) {
	 com->fprintf(stderr, "  ... ERROR IN SCALAPACK ROUTINE PDGELS!\n");
	 exit(-1);
 }

 com->barrier();

 //===============================
 // put least squares solution (contained in subMatB) back in correct form
 //===============================

 // output: lsSol

 transferDataBackLS(subMatB, n, lsSol, nRhs, subMatLLD);

 com->barrier();

 //===============================
 // clean up
 //===============================

 delete[] work;
 delete[] locSendReceive;
 delete[] subMatA;
 delete[] subMatB;

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH SCALAPACK and DO_SCALAPACK Flag\n");
 exit(-1);
#endif

}

template<int dim>
void ParallelRom<dim>::transferData(VecSet< DistSVec<double, dim> > &snaps, double* subMat, int nSnaps) {

	//====================================
	// Purpose: fill subMat by sending/receiving nodes where needed
	// Inputs: snaps, cpuMasterNodes, cpuNodes, nSnaps
	// Outputs: locSendReceive (locSendReceive[jCpu] = n: this cpu sends n nodes to jCpu), subMat
	//====================================
	//====================================
	// Strategy:
	//  for sender cpus, the first nodes encountered (first master nodes of first subdomains) are the ones that are sent
	//  subMat is arranged as a vector with the first column first, second column second, etc.  The first nodes are those received from senders.
	//====================================

 // distribution info

 int numLocSub = domain.getNumLocSub();
 DistInfo &nodeDistInfo = domain.getNodeDistInfo();
 nTotCpus = com->size(); 
 thisCPU = com->cpuNum(); 

 int *cpuMasterNodesCopy = new int[nTotCpus];
 for (int iCpu=0; iCpu<nTotCpus; ++iCpu)
   cpuMasterNodesCopy[iCpu] = cpuMasterNodes[iCpu];

 int subMatRowsNum = cpuNodes[thisCPU]*dim;	// number of rows for this cpu
 double *recData;	// array for received data
 double *buffer;	// buffer is filled in for all nodes first, then columns [node1col1, node2col1, ...]^T

 // send/receive for each cpu (if positive, it has received, if negative, it has sent)
 for (int iCpu = 0; iCpu < nTotCpus; iCpu++)
   locSendReceive[iCpu] = 0;

 // loop over all domain nodes and compute the cpu to send to
 // -1 indicates that it's a slave node with no sending
 DistVec<int> cpuDestination(domain.getNodeDistInfo());
 cpuDestination = -1;

 int *totalSentNodes = new int[nTotCpus];	// how many nodes have been sent to previous cpus
 int nSendNodes,buffLen,index;
 int nRecNodes = 0;
 int masterNodeCount;	// counter for master nodes

 //============================
 // send/receive nodes where needed (rich give to the poor)
 // the received nodes become the first nodes of subMat
 //============================

 for (int iCpu = 0; iCpu < nTotCpus; ++iCpu) {
   totalSentNodes[iCpu] = 0;
   if (cpuMasterNodesCopy[iCpu] - cpuNodes[iCpu] > 0){	// if more master nodes than needed, iCpu is rich and must send to the poor!
     for (int jCpu = 0; jCpu < nTotCpus; ++jCpu){
       if ( cpuMasterNodesCopy[jCpu] - cpuNodes[jCpu] < 0) {	// if fewer master nodes than needed, jCpu is poor and must receive some nodes from the rich
         nSendNodes = min(cpuMasterNodesCopy[iCpu] - cpuNodes[iCpu],-(cpuMasterNodesCopy[jCpu] - cpuNodes[jCpu]));	// either get all nodes from iCpu or fill all missing nodes for jCpu (whichever happens first)
         //send from iCpu to jCpu
         if (thisCPU==iCpu)
           locSendReceive[jCpu] = -nSendNodes;
         else if(thisCPU==jCpu)
           locSendReceive[iCpu] = nSendNodes;

         buffLen = nSendNodes*dim*nSnaps;	// send all rows for chosen nodes and all columns (nSnaps)
         // create array for received data
         if (thisCPU==jCpu)	// thisCpu is receiver
           recData = new double[buffLen];

         if (thisCPU==iCpu) {	// thisCpu is sender
           // create buffer
           buffer = new double[buffLen];
           //fill buffer
           for (int iSnap = 0; iSnap < nSnaps; ++iSnap)  {	// loop over columns
             masterNodeCount = 0;	// counter for master nodes
             for (int iSub = 0; iSub < numLocSub; ++iSub) {	// loop over subdomains
               // create buffer
               double (*locSnap)[dim] = snaps[iSnap].subData(iSub);
               bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
               for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {	// loop over master nodes
                 if (locMasterFlag[iNode]) {
                   if (masterNodeCount == nSendNodes+totalSentNodes[iCpu])  break;	// only send nSendNodes starting at totalSentNodes (other cpus have received the first 0 to totalSentNodes-1 nodes!)
                   masterNodeCount++;
                   if (masterNodeCount-1 >= totalSentNodes[iCpu])  {
                     for (int j = 0; j < dim; ++j)
                       buffer[(masterNodeCount-totalSentNodes[iCpu]-1)*dim+j+iSnap*nSendNodes*dim] = locSnap[iNode][j];
                   }
                 }
               }
             }
           }
           totalSentNodes[iCpu] = masterNodeCount;

           // send data in buffer
           //fprintf(stderr, "*** CPU #%d sending to CPU #%d: %d nodes = %d entries\n", thisCPU, jCpu, nSendNodes, buffLen);
           com->sendTo(jCpu, thisCPU*nTotCpus+jCpu, buffer, buffLen);
           com->waitForAllReq();
         }//endstuff

         // receive data and populate submatrix

         if (thisCPU==jCpu) {
           com->recFrom(iCpu, iCpu*nTotCpus+thisCPU, recData, buffLen);
           for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
             for (int k = nRecNodes; k < nSendNodes+nRecNodes; ++k)  {
               for (int j = 0; j < dim; ++j)
                 subMat[iSnap*subMatRowsNum + k*dim+j] = recData[(k-nRecNodes)*dim+j+iSnap*nSendNodes*dim];	// arrange by column with the first nodes associated with the received nodes
             }
           }
           nRecNodes += nSendNodes;
         }
         cpuMasterNodesCopy[jCpu] += nSendNodes;	// update number of master nodes on receiver
         cpuMasterNodesCopy[iCpu] -= nSendNodes;	// update number of master nodes on sender
         if (thisCPU==iCpu)
           delete[] buffer;
         else if (thisCPU==jCpu)
           delete[] recData;
       }
     }
     com->barrier();
   }
 }//endfor

 //==================================
 // fill the last nodes of subMat with its own data
 //==================================

 int k;
 for (int iCpu = 0; iCpu < nTotCpus; ++iCpu) {
   if (thisCPU == iCpu) {	// execute for current cpu
     for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
       k = nRecNodes;	// initial nodes are the ones already received
       masterNodeCount = 0;
       for (int iSub = 0; iSub < numLocSub; ++iSub) {
         double (*locSnap)[dim] = snaps[iSnap].subData(iSub);
         bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
         for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
           if (locMasterFlag[iNode]) {
             //if not already sent
             if (masterNodeCount >= totalSentNodes[iCpu])  {
               for (int j = 0; j < dim; ++j)
                 subMat[iSnap*subMatRowsNum + k*dim+j] = locSnap[iNode][j];
               if (iSnap==nSnaps-1)
                 cpuMasterNodesCopy[iCpu]--;
               k++;
             }
             masterNodeCount++;
           }
         }
       }
     }
   }
 }

 delete[] totalSentNodes;
 delete[] cpuMasterNodesCopy;
	// NOTE: max size of subMat is (nSnap-1)* subMatRowsNum * + K*dim
}

//----------------------------------------------------------------------------------
template<int dim>
void ParallelRom<dim>::transferDataBack(double *U, VecSet< DistSVec<double, dim> > &Utrue , int nSnaps) {

 int numLocSub = domain.getNumLocSub();
 DistInfo &nodeDistInfo = domain.getNodeDistInfo();
 nTotCpus = com->size(); 
 thisCPU = com->cpuNum();

 double *recData;
 double *buffer;

 int locLLD = dim * cpuNodes[thisCPU];
 // loop over all domain nodes and compute the cpu to send to
 // -1 indicates that it's a slave node and no sending
 DistVec<int> cpuDestination(domain.getNodeDistInfo());
 cpuDestination = -1;

 int *totalSentNodes = new int[nTotCpus];
 int nSendNodes,buffLen,index;
 int transf = 0;
 int nRecNodes = 0;
 com->barrier();
 for (int iCpu=0; iCpu < nTotCpus; ++iCpu) {
   totalSentNodes[iCpu] = 0;
   for (int jCpu=0; jCpu < nTotCpus; ++jCpu) {
     nSendNodes = 0;
     buffLen = 0;
     transf = 0;
     if (thisCPU==iCpu) {
       if (locSendReceive[jCpu] > 0) {
         nSendNodes = locSendReceive[jCpu];
         buffLen = nSendNodes*dim*nSnaps;
         transf = 1;
       }
     }
     else if (thisCPU==jCpu) {
       if (locSendReceive[iCpu] < 0) {
         nSendNodes = -locSendReceive[iCpu];
         buffLen = nSendNodes*dim*nSnaps;
         transf = 1;
       }
     }
     com->barrier();
     com->globalSum(1, &transf);
     if (transf > 1) {
       // create array for received data
       if (thisCPU==jCpu)
         recData = new double[buffLen];
       com->barrier();
       //create buffer
       if (thisCPU==iCpu) {
         buffer = new double[buffLen];
          //fill buffer
         for (int iSnap = 0; iSnap < nSnaps; ++iSnap)  {
           index = 0;
           for (int iSub = 0; iSub < numLocSub; ++iSub) {
             // create buffer
             bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
             for (int iNode = 0; iNode < cpuNodes[iCpu]; ++iNode)  {
               if (index == nSendNodes+totalSentNodes[iCpu])  break;
               index++;
               if (index-1 >= totalSentNodes[iCpu])  {
                 for (int j = 0; j < dim; ++j)
                   buffer[(index-totalSentNodes[iCpu]-1)*dim+j+iSnap*nSendNodes*dim] = U[locLLD*iSnap+dim*iNode+j];
               }
             }
           }
         }
         totalSentNodes[iCpu] = index;
         // send data in buffer
         //fprintf(stderr, "*** CPU #%d sending back to CPU #%d: %d nodes = %d entries\n", thisCPU, jCpu, nSendNodes, buffLen);
         com->sendTo(jCpu, 10*thisCPU*nTotCpus+jCpu, buffer, buffLen);
       }
       com->barrier();
       // receive data and populate submatrix
       if (thisCPU==jCpu) {
         com->recFrom(iCpu, 10*iCpu*nTotCpus+thisCPU, recData, buffLen);
         for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
           index = 0;
           int k=0;
           for (int iSub = 0; iSub < numLocSub; ++iSub) {
             double (*locU)[dim] = Utrue[iSnap].subData(iSub);
             bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
             for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
               if (locMasterFlag[iNode]) {
                 if (k>=nRecNodes && k <nSendNodes+nRecNodes) {
                   for (int j = 0; j < dim; ++j)
                     locU[iNode][j] = recData[(k-nRecNodes)*dim+j+iSnap*nSendNodes*dim];
                 }
                 k++;
               }
             }
           }
         }
         nRecNodes += nSendNodes;
         delete[] recData;

       }
       com->barrier();
       if (thisCPU==iCpu)
         delete[] buffer;
     }
     com->barrier();
   }
 }
 //Completing the rest of the matrices
 for (int iCpu = 0; iCpu < nTotCpus; ++iCpu) {
   if (thisCPU == iCpu) {
     int k;
     for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
       k = totalSentNodes[iCpu];
       index = 0;
       for (int iSub = 0; iSub < numLocSub; ++iSub) {
         double (*locU)[dim] = Utrue[iSnap].subData(iSub);
         bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);
         for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); ++iNode)  {
           if (locMasterFlag[iNode]) {
             //if not already sent
             if (index >= nRecNodes)  {
               for (int j = 0; j < dim; ++j)
                 locU[iNode][j] = U[iSnap*locLLD+k*dim+j];
               k++;
             }
             index++;
           }
         }
       }
     }
   }
 }

 //complete the slave nodes
 CommPattern<double> *vpat = domain.getVecPat();
 for (int iSnap=0; iSnap < nSnaps; ++iSnap)
   domain.assemble(vpat, Utrue[iSnap]);
}

//----------------------------------------------------------------------------------

template<int dim>
void ParallelRom<dim>::transferDataBackLS (double *subMatB, int n, double **lsSol, int nRhs, int subMatLLD) {

	//====================================
	// Purpose: fill lsSol by using appropriate data from subMatB
	// Inputs: subMatB (scalapack fills in B with solution), n (number of columns in A), subMatLLD (so you know where the data for various RHS resides on subMatB)
	// Outputs: lsSol
	//====================================

	// distribution info

	int numLocSub = domain.getNumLocSub();
	DistInfo &nodeDistInfo = domain.getNodeDistInfo();
	nTotCpus = com->size();
	thisCPU = com->cpuNum();

	// info to keep track of the subMatrixEntry for each cpu and the current cpu that is filling lsSol

	int currentCpu = 0;	// start with cpu 0 then move up
	int *subMatEntry = new int[nTotCpus];	// current entry of submatrices (seen by all cpus)
	int *subMatLLDset = new int[nTotCpus];	// current entry of submatrices (seen by all cpus)

	// initialize lsSol

	for (int i = 0; i < nRhs; ++i){
		for (int j = 0; j < n; ++j)
		  lsSol[i][j] = 0.0;
  }

	for (int i = 0; i < nTotCpus; ++i)
		subMatLLDset[i] = 0;
	subMatLLDset[thisCPU] = subMatLLD;

	com->barrier();

	com->globalSum(nTotCpus, subMatLLDset);

	// compute lsSol from distributed data

  for (int iRhs = 0; iRhs < nRhs; ++iRhs) {
		
		// initialize subMatEntry for this RHS
	  for (int iCpu = 0; iCpu < nTotCpus; ++iCpu)
		  subMatEntry[iCpu] = iRhs * subMatLLDset[iCpu];
	  for (int lsSolCount = 0; lsSolCount < n; ++lsSolCount) {	// loop on entries of the solution
			// currentCpu: cpu which is filling up lsSol
		  if (thisCPU == currentCpu)	// fill lsSol with appropriate subMatB entry
			  lsSol[iRhs][lsSolCount] = subMatB[subMatEntry[currentCpu]];
		  ++subMatEntry[currentCpu];	// current cpu has moved on
			//KTC CHECK 
			if (subMatEntry[currentCpu] % rowsPerBlock == 0) {	// reached the end of a block, so go to the next cpu 
				++currentCpu;
			  if (currentCpu == nTotCpus) 	// should cycle back to zero 
				  currentCpu = 0;
		  }
		}
		currentCpu = 0;

	}


	com->barrier();

	for (int iRhs = 0; iRhs < nRhs; ++iRhs)
		com->globalSum(n, lsSol[iRhs]);

	delete[] subMatEntry;
	delete[] subMatLLDset;

}

//----------------------------------------------------------------------------------

template<int dim>
void ParallelRom<dim>::setTransfer() {

 //======================================
 // Purpose: determine number of nodes needed for block cyclic decomposition of scalapack
 // Inputs: nodesPerBlock
 // Outputs: cpuMasterNodes (number of master nodes on cpu for given mesh decomposition), cpuNodes (number of master nodes needed by cpu for scalapack block cyclic decomposition), maxCpuBlocks (number of maximum blocks per cpu for block cyclic decomposition)
 //======================================

 int numLocSub = domain.getNumLocSub();

 // loop over all domain nodes and compute the cpu to send to
 // -1 indicates that it's a slave node and no sending
 DistVec<int> cpuDestination(domain.getNodeDistInfo());
 cpuDestination = -1;

 subDomain = domain.getSubDomain();

 nTotCpus = com->size(); 
 thisCPU = com->cpuNum(); 

 //======================================
 // compute master nodes (non-overlapping subdomain sizes)
 //======================================

 for (int iCpu = 0; iCpu < nTotCpus; iCpu++)
   cpuMasterNodes[iCpu] = 0;

 #pragma omp parallel for
 for (int iSub = 0; iSub < numLocSub; ++iSub) {
   bool *locMasterFlag = cpuDestination.getMasterFlag(iSub);
   for (int iNode = 0; iNode < subDomain[iSub]->numNodes(); iNode++)  {
     if (locMasterFlag[iNode])
       cpuMasterNodes[thisCPU]++;	// add if it is a master node
   }
 }
 
 // maxDomainSize = maximum number of master nodes on a cpu
 int maxDomainSize = cpuMasterNodes[thisCPU];
 com->globalMax(1, &maxDomainSize);
 //make consistent globally
 com->globalSum(nTotCpus, cpuMasterNodes);

 //======================================
 // compute cpuNodes = the number of nodes needed by cpu for block cyclic decomposition
 //======================================

 int numTotalNodes = 0;	// number of total nodes in domain
 for (int iCpu=0; iCpu < nTotCpus; ++iCpu)
   numTotalNodes += cpuMasterNodes[iCpu];
 //com->fprintf(stderr, "There are %d nodes in total \n",numTotalNodes);

 int numBlocks = int(ceil(double(numTotalNodes)/double(nodesPerBlock))); //total number of blocks required (last may not be full)
 int nNodeLastBlock = nodesPerBlock - (numBlocks*nodesPerBlock - numTotalNodes);	// number of nodes in last block
 maxCpuBlocks = int(ceil(double(numBlocks)/double(nTotCpus)));	// maximum number of blocks for a cpu
 int *cpuBlocks = new int[nTotCpus];	// actual number of blocks per cpu
 int notMaxBlockCpus= nTotCpus*maxCpuBlocks-numBlocks;	// number of cpus without the full number of blocks
 for (int iCpu=0; iCpu < nTotCpus; ++iCpu)
   cpuBlocks[iCpu] = maxCpuBlocks;
 for (int i=1; i<=notMaxBlockCpus; ++i)	// correct number of blocks: the last cpus overestimated this
   cpuBlocks[nTotCpus-i] -= 1;
 int lastCpu;
 if (nNodeLastBlock>0) {	// this should ALWAYS hold
   lastCpu = nTotCpus - notMaxBlockCpus -1;	// the -1 goes from number of cpus (starts at 1) to index (starts at 0) // ???? THIS SEEMS DANGEROUS! COULD BE NEGATIVE ????
 }
 else {	// should never hold
	com->fprintf(stderr, "ERROR: nNodeLastBlock not positive!\n");
   lastCpu = nTotCpus - notMaxBlockCpus;
 }

 // compute the number of nodes needed by cpu

 for (int iCpu = 0; iCpu < nTotCpus; ++iCpu)
   cpuNodes[iCpu] = cpuBlocks[iCpu]*nodesPerBlock;	// blocks * number of nodes per block assuming full blocks
 if (nNodeLastBlock < nodesPerBlock)	// last block is not full
   cpuNodes[lastCpu] += (nNodeLastBlock-nodesPerBlock);	// subtract off the missing nodes in not full block
 
 delete[] cpuBlocks;

}

