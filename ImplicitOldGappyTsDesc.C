//------------------------------------------------------------------------------

template<int dim>
ImplicitOldGappyTsDesc<dim>::ImplicitOldGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom) {
  
		globalSubSet = 0;
		locNodeSet = 0;

		dom->readInterpNode(ioData.input.sampleNodes, nIntNodes, globalSubSet, locNodeSet);
		if (ioData.input.aMatrix) dom->readInterpMatrix(ioData.input.aMatrix, dimInterpMat, interpMat1);
		if (ioData.input.bMatrix) dom->readInterpMatrix(ioData.input.bMatrix, dimInterpMat, interpMat2);

		computeRestrictInfo();
}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitOldGappyTsDesc<dim>::solveNewtonSystem()  {

		computeAJGappy(it, U, F, rhs, AJ);

		// saving residual vectors (for GappyPOD)
		//writeBinaryVectorsToDisk1(false, _it, 0.0, F, Dummy);


    res = rhs*rhs;

    if (res < 0.0){
      fprintf(stderr, "*** negative residual: %e\n", res);
      exit(1);
    }
    res = sqrt(res);

    if (it == 0) {
      target = epsNewton*res;
      res0 = res;
    }

    if (res == 0.0 || res <= target) break;

    solveLinearSystem(it, rhs, dUrom);
}


//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::computeRestrictInfo() {

  // computes 1) myNNodeInt: the number of interpolation nodes contained on this CPU
  //          2) myLocalSubSet: the local subdomain numbers
  //          3) myInterpNodes: the global interpolation nodes it contains
  // NOTE: nodes are lumped according to subdomain for efficiency

  // local CPU information

  int numLocSub = this->domain->getNumLocSub();
  DistInfo &nodeDistInfo = this->domain->getNodeDistInfo();
  SubDomain** subD = this->domain->getSubDomain();
  myNNodeInt = 0;
  int globalSubNum;
  int localNode;

  int *myInterpNodesTemp = new int[nIntNodes];
  int *myLocalSubSetTemp = new int[nIntNodes];

  for (int iSub = 0; iSub < numLocSub; ++iSub) { // loop on subdomains on this CPU
    globalSubNum = subD[iSub]->getGlobSubNum(); 
    bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub);  // array of locMasterFlag
    for (int iNodeInt = 0; iNodeInt < nIntNodes; ++iNodeInt) {  // loop over global interpolation nodes (zero to nIntNodes-1)
      if (globalSubNum == globalSubSet[iNodeInt]){ // if the current interpolation node is on this subdomain
        localNode = locNodeSet[iNodeInt]; // the local node of interest
        if (locMasterFlag[localNode]){ // if this node is the master on this domain
          myInterpNodesTemp[myNNodeInt] = iNodeInt;
          myLocalSubSetTemp[myNNodeInt]= iSub;
          ++myNNodeInt;
        }
      }
    }
    //delete [] locMasterFlag;
  }

  myInterpNodes = new int[myNNodeInt];
  myLocalSubSet = new int[myNNodeInt];
  for (int i = 0; i < myNNodeInt; ++i) {
    myInterpNodes[i] = myInterpNodesTemp[i];
    myLocalSubSet[i] = myLocalSubSetTemp[i];
  }
/*
  delete [] myInterpNodesTemp;
  delete [] myLocalSubSetTemp;
  delete [] subD;
*/
}
template<int dim>
void ImplicitRomTsDesc<dim>::computeAJGappy(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F, 
                                                  Vec<double> &rhs, VecSet<DistSVec<double, dim> > &AJ)  {
  
  // The goal of this function is to compute rhs and jac for the newton step
  // nIntNodes = # nodes at which interpolation occurs
  // dimension of interpMat1 = nIntNodes * dim, interpMat2 = nIntNodes * dim

  int debugging = 0;

  mvpfd->evaluate(it, *this->X, *this->A, Q, F);  // prepares mvpfd->apply

  jac.setNewSize(nPod,nPod); // actual Jacobian to be used in Newton iterations

  for (int iPod = 0; iPod < nPod; ++iPod)
    mvpfd->apply(pod[iPod], AJ[iPod]);  // AJ is AJ: the VecSet of DistSVec with AJ[i] = dr/dw*Phi[i]

  FullM AJRestrict(nIntNodes*dim,nPod); // the matrix of dr/dw*Phi[i] restricted to interp nodes
  AJRestrict = 0.0;
  FullM resRestrict(nIntNodes*dim,1); // residual restricted to interp nodes
  resRestrict = 0.0;

  // info for parallel operations
  int numLocSub = this->domain->getNumLocSub();
  int nTotCpus = this->com->size();
  int thisCPU = this->com->cpuNum();
  DistInfo &nodeDistInfo = this->domain->getNodeDistInfo();
  SubDomain** subD = this->domain->getSubDomain();

  // compute AJ and F restricted to the specified nodes
  int oldLocalSub = -1; 
  int localNode,localSub,currentNodeIndex;
  for (int iPod = 0; iPod < nPod; ++iPod) { // loop over all POD vectors
    for (int iMyInterpNode = 0; iMyInterpNode < myNNodeInt; ++iMyInterpNode) { // loop over local interpolation nodes
      currentNodeIndex = myInterpNodes[iMyInterpNode];  // global index of current interpolation node
      localNode = locNodeSet[currentNodeIndex];  // the local node number
      localSub = myLocalSubSet[iMyInterpNode];  // the local subdomain number
      //if (iMyInterpNode == 0 || localSub != oldLocalSub) { // only reload subdomain part of locF, locAJ if needed
      //if (iPod == 0){
      double (*locF)[dim] = F.subData(localSub); // compute local F (should only be loaded for iPod == 0)
      //}
      double (*locAJ)[dim] = AJ[iPod].subData(localSub); // compute local AJ
      //}

      // KTC test: are the CPU #, global Sub #, local Node #, and myNNodeInt correct for a certain global currentNodeIndex?
      if (debugging) {
        fprintf(stderr,"currentNodeIndex %d, CPU %d, globalSub %d, localNode %d, localSub %d, myNNodeInt %d\n",currentNodeIndex,thisCPU,globalSubSet[currentNodeIndex],localNode,localSub, myNNodeInt);
      }
      for (int iDim = 0 ; iDim < dim; ++iDim) {
        // compute restricted residual and jacobian
        if (iPod == 0) resRestrict[currentNodeIndex*dim+iDim][0] = locF[localNode][iDim];  // compute for RHS (only one)
        AJRestrict[currentNodeIndex*dim+iDim][iPod] = locAJ[localNode][iDim]; // fill in the matrix of dr/dw*Phi[i] restricted to interp
      }
      oldLocalSub = localSub;
    }
  }
/*
  delete [] locAJ;
  delete [] locF;
  delete [] subD;
*/

  this->com->globalSum(nIntNodes*dim*nPod,AJRestrict.data()); // ensure all CPUs have the same copy
  this->com->globalSum(nIntNodes*dim,resRestrict.data()); // ensure all CPUs have the same copy

/*
  FullM Arhat(nIntNodes*dim,1);
  Arhat = interpMat*resRestrict;
  FullM rhatArhatmat = resRestrict^Arhat;
  double rhatArhat = rhatArhatmat[0][0];
  this->com->fprintf(stderr,"R^T A R = %e\n",rhatArhat);
*/
  // parallel implementation of AJRestrictTMat1 = AJRestrict^T * interpMat1
  // parallel implementation of AJRestrictTMat2 = AJRestrict^T * interpMat2

  int maxIndex, loadBal, loadBalMod, myMinIndex, myMaxIndex;  // indices used for parallel operations

  FullM AJRestrictTMat1(AJRestrict.numCol(),interpMat1.numCol());
  AJRestrictTMat1 = 0.0;  // initialize to zero

  rowPartition(myMinIndex,myMaxIndex,AJRestrictTMat1.numRow());
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < AJRestrictTMat1.numCol(); ++j){
      for (int k = 0; k < interpMat1.numRow(); ++k) {
        AJRestrictTMat1[i][j] += AJRestrict[k][i]*interpMat1[k][j];
      }
    }
  }
  this->com->globalSum(AJRestrictTMat1.numRow()*AJRestrictTMat1.numCol(),AJRestrictTMat1.data()); // ensure all CPUs have the same copy

  FullM AJRestrictTMat2(AJRestrict.numCol(),interpMat2.numCol());
  AJRestrictTMat2 = 0.0;  // initialize to zero

  rowPartition(myMinIndex,myMaxIndex,AJRestrictTMat2.numRow());
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < AJRestrictTMat2.numCol(); ++j){
      for (int k = 0; k < interpMat2.numRow(); ++k) {
        AJRestrictTMat2[i][j] += AJRestrict[k][i]*interpMat2[k][j];
      }
    }
  }
  this->com->globalSum(AJRestrictTMat2.numRow()*AJRestrictTMat2.numCol(),AJRestrictTMat2.data()); // ensure all CPUs have the same copy

  // compute jac = AJRestrictTMat1*AJRestrict

  jac = 0.0;  // initialize to zero
  rowPartition(myMinIndex,myMaxIndex,jac.numRow(),1);  // symmetric
       
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j <= i; ++j){
      for (int k = 0; k < AJRestrict.numRow(); ++k) {
        jac[i][j] += AJRestrictTMat1[i][k]*AJRestrict[k][j];
      }
      if (i > j)
        jac[j][i] = jac[i][j];
    }
  }

  this->com->globalSum(jac.numRow()*jac.numCol(),jac.data()); // ensure all CPUs have the same copy
  // compute rhs = -1.0*AJRestrictTMat2*resRestrict

  rowPartition(myMinIndex,myMaxIndex,rhs.size());
  rhs = 0.0; 
  for (int i = myMinIndex; i < myMaxIndex; ++i){
    for (int j = 0; j < resRestrict.numRow(); ++j)
      rhs[i] -= AJRestrictTMat2[i][j]*resRestrict[j][0]; // NOTE: the RHS is NEGATIVE!
  }
  this->com->globalSum(rhs.size(),rhs.data()); // ensure all CPUs have the same copy
 
  // KTC test: given AJRestrict, resRestrict, interpMat1, and interpMat2, are the computation of jacobian and rhs correct?
  FullM AJRestrictTMat1test, AJRestrictTMat2test, jacTest,rhsTest;
  if (debugging){
    AJRestrictTMat1test= AJRestrict^interpMat1;
    AJRestrictTMat2test= AJRestrict^interpMat2;
    jacTest = AJRestrictTMat1test*AJRestrict;
    rhsTest = AJRestrictTMat2test*resRestrict;
    rhsTest *= -1.0; 
    this->com->fprintf(stderr,"jac is \n");
    if (thisCPU == 0 ) {jac.print();}
    this->com->fprintf(stderr,"jacTest is \n");
    if (thisCPU == 0 ) {jacTest.print();}
    if (thisCPU == 0 ) {fprintf(stderr, "rhs is \n");}
    for (int j = 0; j < rhs.size(); ++j) { 
      if (thisCPU == 0 ) {fprintf(stderr, "%e \n", rhs[j]);}
    }
    this->com->fprintf(stderr,"rhsTest is \n");
    if (thisCPU == 0 ) {rhsTest.print();}
  }

}


//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::rowPartition(int &myMinIndex, int &myMaxIndex, int nRow, int sym ) {
 
 // this function partitions the rows of a matrix across processors (load balancing)

 long int loadBal;
 int loadBalMod;
 int nTotCpus = this->com->size();
 int thisCPU = this->com->cpuNum();
 switch (sym) {  // is the operation symmetric?
   case 0:  // evenly split rows across processors
     loadBal = nRow/nTotCpus;
     loadBalMod = nRow % nTotCpus;
     myMinIndex = loadBal*thisCPU + (thisCPU < loadBalMod)*thisCPU + (thisCPU >= loadBalMod)*loadBalMod;
     myMaxIndex = loadBal*(thisCPU+1) + ((thisCPU+1) < loadBalMod)*(thisCPU+1) + ((thisCPU +1) >= loadBalMod)*loadBalMod;
   break;
   case 1:  // split rows in a staggered manner
     loadBal = (nRow*nRow)/nTotCpus;
     int *minVal = new int[nTotCpus];
     int *maxVal = new int[nTotCpus];
     minVal[0] = 0;
     maxVal[0] = (int) sqrt(double(loadBal));
     for (int i = 1; i < nTotCpus; ++i){
       minVal[i] = maxVal[i-1];
       maxVal[i] = (int) sqrt(double(loadBal+ minVal[i]*minVal[i]));
     }    
     maxVal[nTotCpus - 1] = nRow;
     myMinIndex = minVal[thisCPU];
     myMaxIndex = maxVal[thisCPU];
     delete [] minVal;
     delete [] maxVal;
   break;
 }
 if (myMaxIndex > nRow || myMinIndex > nRow || myMinIndex > myMaxIndex) {
   fprintf(stderr, "*** Problem with rowPartition!!!");
 }
 // KTC: output minVal, maxVal?
}
