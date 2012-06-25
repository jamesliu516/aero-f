#include <MultiGridLevel.h>

#include <Connectivity.h>
#include <Domain.h>
#include <Edge.h>

#include <set>

template<class Scalar>
MultiGridLevel<Scalar>::MultiGridLevel(Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo)
  : nodeIdPattern(0), nodeVolPattern(0), nodePosnPattern(0),nodeVecPattern(0), domain(domain),
    nodeDistInfo(new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com)),
    edgeDistInfo(new DistInfo(refinedEdgeDistInfo.numLocThreads, refinedEdgeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedEdgeDistInfo.locSubToGlobSub, refinedEdgeDistInfo.com)),
    numLocSub(refinedNodeDistInfo.numLocSub), ownsData(true), connectivity(new Connectivity*[numLocSub]),
    edges(new EdgeSet*[numLocSub]), faces(new FaceSet*[numLocSub]), nodeMapping(refinedNodeDistInfo), edgeMapping(refinedEdgeDistInfo),elems(new ElemSet*[numLocSub]),
    distGeoState(0), lineMap(refinedNodeDistInfo),lineIDMap(refinedNodeDistInfo),lineLocIDMap(refinedNodeDistInfo)
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = 0;
    edges[iSub] = 0;
    faces[iSub] = 0;
    elems[iSub] = 0;
  }
  lineIDMap = -1;

  numLines = new int[numLocSub];
  memset(numLines,0,sizeof(int)*numLocSub);
  lineids = new std::vector<int>[numLocSub];

  lineLengths = new int*[numLocSub];
 
}

template<class Scalar>
void MultiGridLevel<Scalar>::WriteTopFile(const std::string& fileName) {

  std::ofstream outfile(fileName.c_str());
  int nodenum = 1;
  int* node_starts = new int[numLocSub];
  outfile << "Nodes FluidNodes\n";
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    node_starts[iSub] = nodenum;
    SVec<double,3>& Xn = (distGeoState->getXn())(iSub);
    for (int i = 0; i < Xn.size(); ++i,++nodenum)
      outfile << nodenum << " " << Xn[i][0] << " " << Xn[i][1] << " " << Xn[i][2] << "\n";
  }

  int elem_num = 1;
  outfile << "Elements FluidMesh using FluidNodes\n";
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    ElemSet& elem_set = *elems[iSub];
    for (int i = 0; i < elem_set.size(); ++i,++elem_num) {

      outfile << elem_num << " 5 ";
      for (int j = 0; j < 4; ++j)
        outfile << node_starts[iSub] + elem_set[i][j] << " ";
      outfile << "\n";
    }

  }

  outfile << "Elements StickMovingSurface_0 using FluidNodes\n";
  elem_num = 1;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    FaceSet& elem_set = *faces[iSub];
    for (int i = 0; i < elem_set.size(); ++i,++elem_num) {

      outfile << elem_num << " 4 ";
      for (int j = 0; j < 3; ++j)
        outfile << node_starts[iSub] + elem_set[i][j] << " ";
      outfile << "\n";
    }
  }
  

  outfile.close();
}

template<class Scalar>
template<int dim>
void MultiGridLevel<Scalar>::writeXpostFile(const std::string& fileName,
                                            DistSVec<Scalar,dim>& val,int id) {

  std::ofstream outfile(fileName.c_str());
  outfile << "Scalar Pressure under load for FluidNodes\n";
  int nodes_inc_overlap = 0;
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    nodes_inc_overlap += val.subSize(iSub);
  }
  outfile << nodes_inc_overlap << "\n0\n";
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    for (int i = 0; i < val.subSize(iSub); ++i) {

      outfile << val(iSub)[i][id] << "\n";
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
MultiGridLevel<Scalar>::~MultiGridLevel()
{
  if(ownsData){
    delete nodeIdPattern;
    delete nodeVolPattern;
    delete nodeVecPattern;
    delete nodePosnPattern;
  }

  delete nodeDistInfo;
  delete edgeDistInfo;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub){
    if(ownsData) {
      delete connectivity[iSub];
      delete edges[iSub];
      delete faces[iSub];
    }
  }
  delete []sharedNodes;
  delete []connectivity;
  delete []edges;
  delete []faces;
  
  delete [] lineids;
  delete [] numLines;

  if(ownsData) delete distGeoState;
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo,
                                              DistGeoState& refinedGeoState, Domain& domain)
{
  ownsData = false;
  nodeIdPattern = domain.getLevelPat();
  nodeVolPattern = domain.getVolPat();
  nodeVecPattern = domain.getVecPat();
  nodePosnPattern = domain.getVec3DPat();
  sharedNodes = new Connectivity*[numLocSub];

  sharedEdges = new EdgeDef**[numLocSub];
  numSharedEdges = new int*[numLocSub]; 

  distGeoState = &refinedGeoState;
 
  finestNodeMapping = NULL;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = domain.getSubDomain()[iSub]->getNodeToNode();
    sharedNodes[iSub] = domain.getSubDomain()[iSub]->getSharedNodes();
    edges[iSub] = &domain.getSubDomain()[iSub]->getEdges();
    faces[iSub] = &domain.getSubDomain()[iSub]->getFaces();
    elems[iSub] = &domain.getSubDomain()[iSub]->getElems();
    nodeDistInfo->setLen(iSub, refinedNodeDistInfo.subSize(iSub));
    edgeDistInfo->setLen(iSub, refinedEdgeDistInfo.subSize(iSub));

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    numSharedEdges[iSub] = new int[subD.getNumNeighb()];
    sharedEdges[iSub] = new EdgeDef*[subD.getNumNeighb()];
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
 
      numSharedEdges[iSub][jSub] = subD.getNumSharedEdges()[jSub];
      sharedEdges[iSub][jSub] = new EdgeDef[numSharedEdges[iSub][jSub]];
      memcpy(sharedEdges[iSub][jSub], subD.getSharedEdges()[jSub], 
             sizeof(EdgeDef)*numSharedEdges[iSub][jSub]);
    }
  }
  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < refinedNodeDistInfo.subSize(iSub); ++i) {
      nodeDistInfo->getMasterFlag(iSub)[i] = refinedNodeDistInfo.getMasterFlag(iSub)[i];
      nodeDistInfo->getInvWeight(iSub)[i] = refinedNodeDistInfo.getInvWeight(iSub)[i];
    }
}

//------------------------------------------------------------------------------

namespace { // MPI helper functions
template<class T,class OperType>
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, DistVec<T> & data, const OperType & op) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      T * buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        buffer[iNode] = data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ];
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      T * buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ] = OperType::apply(data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ], buffer[iNode]);
    }
  }
}

template<class T,class OperType,int dim>
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, DistSVec<T,dim> & data, const OperType & op,std::map<int,std::map<int,int> >* locToGlobMap = NULL) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        for (int j = 0; j < dim; ++j)
          buffer[iNode][j] = data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][j];
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        /*if (locToGlobMap && (*locToGlobMap)[iSub][(*sharedNodes[iSub])[jSub][iNode]] == 3911)
          std::cout << iSub << " " << data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][0] <<
                    " " << buffer[iNode][0] << std::endl;
        */
        for (int j = 0; j < dim; ++j)
          data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][j] = OperType::apply(data(iSub)[ (*sharedNodes[iSub])[jSub][iNode]][j], buffer[iNode][j]);
      }
    }
  }
}

template<class T,int dim>
void assemble(Domain & domain, CommPattern<T> & commPatDiag, Connectivity ** sharedNodes, DistMat<T,dim> & A) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = commPatDiag.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A(iSub).getElem_ii((*sharedNodes[iSub])[jSub][iNode]);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
  }

  commPatDiag.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPatDiag.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A(iSub).getElem_ii((*sharedNodes[iSub])[jSub][iNode]);
        for (int j = 0; j < dim*dim; ++j) a[j] += buffer[iNode][j];
      }
    }
  }

}

template<class T,int dim>
void assemble(Domain & domain, CommPattern<T> & commPatOffDiag, int **numSharedEdges, EdgeDef ***sharedEdges, DistMat<T,dim> & A) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = commPatOffDiag.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[2][dim*dim] = reinterpret_cast<T (*)[2][dim*dim]>(sInfo.data);

      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]/*sharedNodes[iSub]->num(jSub)*/; ++iEdge) {
       
        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;

        T *aij, *aji;

        if (sharedEdges[iSub][jSub][iEdge].sign > 0) {
  	  aij = A(iSub).getElem_ij(edgeNum);
	  aji = A(iSub).getElem_ji(edgeNum);
        }
        else {
 	  aij = A(iSub).getElem_ji(edgeNum);
	  aji = A(iSub).getElem_ij(edgeNum);
        }

        if (aij && aji) {
  	  for (int k=0; k<dim*dim; ++k) {
	    buffer[iEdge][0][k] = aij[k];
	    buffer[iEdge][1][k] = aji[k];
	  }
        }
        else {
 	  for (int k=0; k<dim*dim; ++k) {
	    buffer[iEdge][0][k] = 0.0;
	    buffer[iEdge][1][k] = 0.0;
	  }
        }      
      }
    }
  }

  commPatOffDiag.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPatOffDiag.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[2][dim*dim] = reinterpret_cast<T (*)[2][dim*dim]>(sInfo.data);
      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {

        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;

        T *aij, *aji;

        if (sharedEdges[iSub][jSub][iEdge].sign > 0) {
  	  aij = A(iSub).getElem_ij(edgeNum);
	  aji = A(iSub).getElem_ji(edgeNum);
        }
        else {
	  aij = A(iSub).getElem_ji(edgeNum);
	  aji = A(iSub).getElem_ij(edgeNum);
        }

        if (aij && aji) {
	  for (int k=0; k<dim*dim; ++k) {
	    aij[k] += buffer[iEdge][0][k];
	    aji[k] += buffer[iEdge][1][k];
 	  }
        }

      }

    }
  }
}
};

template <class Scalar>
void MultiGridLevel<Scalar>::identifyEdges(CommPattern<int> &edgeNumPat,
                                           int mySub) {
 
  int iSub, jSub;
  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  Connectivity *nodeToNeighb = sharedNodes[mySub]->reverse();

  int numNeighb = subD.getNumNeighb();

  numSharedEdges[mySub] = new int[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[mySub][iSub] = 0;

  int (*edgePtr)[2] = edges[mySub]->getPtr();

  int l;
  for (l=0; l<edges[mySub]->size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    /*if (mySub == 1 && leftNode == 384 && rightNode == 429)
      std::cout << "Hello" << std::endl;
*/
    if (nodeToNeighb->num(leftNode) == 0 ||
	nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub)
	if (subI == (*nodeToNeighb)[rightNode][jSub]) numSharedEdges[mySub][subI] += 1;

    }

  }

  sharedEdges[mySub] = new EdgeDef *[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub)
    sharedEdges[mySub][iSub] = new EdgeDef[ numSharedEdges[mySub][iSub] ];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[mySub][iSub] = 0;

  for (l=0; l<edges[mySub]->size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    if (nodeToNeighb->num(leftNode) == 0 ||
       nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub) {

	if (subI == (*nodeToNeighb)[rightNode][jSub]) {

	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].edgeNum = l;
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].glLeft = locToGlobMap[mySub][leftNode];
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].glRight = locToGlobMap[mySub][rightNode];
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].order();

	  numSharedEdges[mySub][subI] += 1;

	}

      }

    }
  }

  for (iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef OLD_STL
    sort(sharedEdges[mySub][iSub], sharedEdges[mySub][iSub]+numSharedEdges[mySub][iSub]);
#else
    stable_sort(sharedEdges[mySub][iSub], sharedEdges[mySub][iSub]+numSharedEdges[mySub][iSub]);
#endif
    edgeNumPat.setLen(subD.getSndChannel()[iSub], 2*numSharedEdges[mySub][iSub]);
  }

  delete nodeToNeighb;

}

template <class Scalar>
void MultiGridLevel<Scalar>::sndEdgeInfo(CommPattern<int> &edgeNumPat,int mySub)
{

  int iSub, iEdge;
  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  int numNeighb = subD.getNumNeighb();

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.getSendBuffer(subD.getSndChannel()[iSub]);

    for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge) {
      sInfo.data[2*iEdge]   = sharedEdges[mySub][iSub][iEdge].glLeft;
      sInfo.data[2*iEdge+1] = sharedEdges[mySub][iSub][iEdge].glRight;
    }

  }
}

template <class Scalar>
void MultiGridLevel<Scalar>::rcvEdgeInfo(CommPattern<int> &edgeNumPat,int mySub,int dim)
{

  int iSub, iEdge;

  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  // We receive the neighbor's list of edges and then also build the edgeMasterFlag

  bool *edgeMasterFlag = new bool[edges[mySub]->size()];

  for (iEdge = 0; iEdge < edges[mySub]->size(); ++iEdge) edgeMasterFlag[iEdge] = true;

  for (iSub = 0; iSub < subD.getNumNeighb(); ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.recData(subD.getRcvChannel()[iSub]);

    int nIndex = 0;
    int myIndex = 0;

    for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge) {

      int glLeft  = sharedEdges[mySub][iSub][iEdge].glLeft;
      int glRight = sharedEdges[mySub][iSub][iEdge].glRight;

      while (2*nIndex < sInfo.len &&
	     (sInfo.data[2*nIndex] < glLeft ||
	      (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] < glRight)))
	nIndex++;

      if (2*nIndex < sInfo.len &&
	  (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] == glRight)) {
	sharedEdges[mySub][iSub][myIndex] = sharedEdges[mySub][iSub][iEdge];
	myIndex++;
      }

    }

    numSharedEdges[mySub][iSub] = myIndex;
    offDiagMatPattern->setLen(subD.getSndChannel()[iSub],2*dim*dim*numSharedEdges[mySub][iSub]);
    offDiagMatPattern->setLen(subD.getRcvChannel()[iSub],2*dim*dim*numSharedEdges[mySub][iSub]);

    if (subD.getNeighb()[iSub] < subD.getGlobSubNum()) // I cannot be the master of these edges
      for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge)
	edgeMasterFlag[ sharedEdges[mySub][iSub][iEdge].edgeNum ] = false;

  }

  int num_edges = 0;
  for (iEdge = 0; iEdge < edges[mySub]->size(); ++iEdge) 
    if (edgeMasterFlag[iEdge]) num_edges++;
 // std::cout << num_edges << std::endl;
  
  edges[mySub]->setMasterFlag(edgeMasterFlag);
}


template<class Scalar>
void MultiGridLevel<Scalar>::agglomerate(const DistInfo& refinedNodeDistInfo,
                                         const DistInfo& refinedEdgeDistInfo,
                                         CommPattern<int>& refinedNodeIdPattern,
                                         DistGeoState& refinedDistGeoState,
                                         Connectivity** refinedSharedNodes,
                                         Connectivity ** nToN, EdgeSet ** refinedEdges,
                                         FaceSet ** refinedFaces,
                                         ElemSet ** refinedElems,
                                         EdgeDef*** refinedSharedEdges,
                                         int** refinedNumSharedEdges,
                                         Domain& domain,int dim,
                                         DistVec<int>* finestNodeMapping_p1)
{
  static int cnt = 0;
  ++cnt;
  int * numNodes = new int[numLocSub];
  nodeIdPattern = new CommPattern<int>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<int>::CopyOnSend);
  nodeVolPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodeVecPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodePosnPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  matPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  offDiagMatPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);

  Connectivity* subToSub = domain.getSubToSub();

  nodeMapping = -1;
  edgeMapping = -1;

  lineMap = -1;

  sharedNodes = new Connectivity*[numLocSub];
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    // We have to ensure that any agglomerated nodes are shared by the
    // correct subdomains (that they can talk to each other)
    Connectivity& rSharedNodes(*refinedSharedNodes[iSub]);
    std::map<int,std::set<int> > nodeSharedData;
    for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      for(int i = 0; i < rSharedNodes.num(jSub); ++i) {
        nodeSharedData[rSharedNodes[jSub][i]].insert(subD.getNeighb()[jSub]);
      }
    }


    for(int i = 0; i < nodeMapping(iSub).size(); ++i) // Disable shared slave nodes
      if(!refinedNodeDistInfo.getMasterFlag(iSub)[i]) nodeMapping(iSub)[i]=0;

    // Perform local agglomeration
    int seed_node = 0, agglom_id = 0;
    lineLengths[iSub] = new int[nodeMapping(iSub).size()];
    while(seed_node < nodeMapping(iSub).size() && nodeMapping(iSub)[seed_node] >= 0) ++seed_node;
    while(seed_node < nodeMapping(iSub).size()) {
      std::vector<int> agglom;
      std::map<int,double> weights;
      agglom.reserve(8);

      agglom.push_back(seed_node);
      std::set<int> agglomSubDs;
      agglomSubDs.insert(subD.getGlobSubNum());
      for (std::set<int>::iterator it = nodeSharedData[seed_node].begin();
           it != nodeSharedData[seed_node].end(); ++it)
        agglomSubDs.insert(*it);

      while(agglom.size() < 8) {
        const int current_node=agglom.back();
        for(int n = 0; n < nToN[iSub]->num(current_node); ++n) {
          const int neighbor = (*nToN[iSub])[current_node][n];
          if(current_node == neighbor || nodeMapping(iSub)[neighbor] >= 0) continue;
          bool ok = true;
          bool has_new = false;
          bool is_superset = true;
          for (std::set<int>::iterator it = nodeSharedData[neighbor].begin();
               it != nodeSharedData[neighbor].end(); ++it) {
          
            if (agglomSubDs.find(*it) == agglomSubDs.end()) {

              has_new = true;
              break;
            }
          }

          for (std::set<int>::iterator it = agglomSubDs.begin();
               it != agglomSubDs.end(); ++it) {

            if (nodeSharedData[neighbor].find(*it) == nodeSharedData[neighbor].end()) {

              is_superset = false;
              break;
            }
          }
          if (has_new && !is_superset) continue;

          if(std::find(agglom.begin(), agglom.end(), neighbor) == agglom.end()) {
            assert(refinedEdges[iSub]->length(refinedEdges[iSub]->findOnly(current_node, neighbor)) > 0.0);
            if (cnt < 4)
              weights[neighbor] += 1.0/refinedEdges[iSub]->length(refinedEdges[iSub]->findOnly(current_node, neighbor));
          }
        }

        int node_to_add = -1;
        double node_weight = -FLT_MAX;
        for(std::map<int,double>::const_iterator iter=weights.begin(); iter != weights.end(); iter++) {
          if(iter->second > node_weight) {
            node_to_add = iter->first;
            node_weight = iter->second;
          }
        }
        if(node_to_add != -1) {
          weights.erase(node_to_add);
          agglom.push_back(node_to_add);
          for (std::set<int>::iterator it = nodeSharedData[node_to_add].begin();
               it != nodeSharedData[node_to_add].end(); ++it) {
            agglomSubDs.insert(*it);
          }
        } else break;
      }

      int j = 0;
      bool isLine = true;
      for(std::vector<int>::const_iterator iter = agglom.begin(); iter != agglom.end(); iter++,++j) {
        const int current_node = *iter;
        //nodeMapping(iSub)[current_node] = agglom_id;
        for(int n = 0; n < nToN[iSub]->num(current_node) && isLine; ++n) {
          const int neighbor = (*nToN[iSub])[current_node][n];
          if (current_node == neighbor) continue;
  
          for (int k = 0; k < agglom.size(); ++k) {
            if ((k > j+1 || k < j-1) && neighbor == agglom[k]) { 
              isLine = false;
              break;
            }
          }
        }

        nodeMapping(iSub)[*iter] = agglom_id;
      }
      
      if (isLine && agglom.size() > 1) {
        int k;
        for (k = 0; k < agglom.size(); ++k) {
  
          assert(refinedNodeDistInfo.getMasterFlag(iSub)[agglom[k]]);
          lineLocIDMap(iSub)[agglom[k]] = k;         
 
          lineIDMap(iSub)[agglom[k]] = numLines[iSub];
          if (k > 0)
            lineMap(iSub)[agglom[k]][0] = agglom[k-1];
          else
            lineMap(iSub)[agglom[k]][0] = -1;
          if (k < agglom.size()-1)
            lineMap(iSub)[agglom[k]][1] = agglom[k+1];  
          else
            lineMap(iSub)[agglom[k]][1] = -1;
          lineids[iSub].push_back(agglom[k]);
        }
        lineLengths[iSub][numLines[iSub]] = k;

        for (; k < 8; ++k)
          lineids[iSub].push_back(-1);

        ++numLines[iSub];
      } else {
        for (int k = 0; k < agglom.size(); ++k) {
          lineMap(iSub)[agglom[k]][0] = lineMap(iSub)[agglom[k]][1] = -1;
        }
      }

      ++agglom_id;
      while(seed_node < nodeMapping(iSub).size() && nodeMapping(iSub)[seed_node] >= 0) ++seed_node;
    }

    numNodes[iSub] = agglom_id;
  }

  domain.getCommunicator()->fprintf(stdout,"hello1\n");
  fflush(stdout);

// ------------------------- MPI PARALLEL CRAP -----------------------------------------------------------------
  DistVec<int> ownerSubDomain(refinedNodeDistInfo);
  operMax<int> maxOp;

  maxNodesPerSubD = 0;
  for(int iSub = 0; iSub < numLocSub; ++iSub) maxNodesPerSubD = max(maxNodesPerSubD, nodeMapping(iSub).size());
  domain.getCommunicator()->fprintf(stdout,"hello6");
  fflush(stdout);
  domain.getCommunicator()->globalMax(1, &maxNodesPerSubD);

  domain.getCommunicator()->fprintf(stdout,"hello5n");
  fflush(stdout);
 
  int tot_nodes = 0;
  for(int iSub = 0; iSub < numLocSub; ++iSub) tot_nodes += numNodes[iSub];
  domain.getCommunicator()->globalSum(1,&tot_nodes);
  //std::cout << "Num Nodes = " << tot_nodes << std::endl;
  domain.getCommunicator()->fprintf(stdout,"hello4\n");
  fflush(stdout);
 
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < ownerSubDomain(iSub).size(); ++i)
      if(!refinedNodeDistInfo.getMasterFlag(iSub)[i])
        nodeMapping(iSub)[i] = ownerSubDomain(iSub)[i] = -1;
      else
        ownerSubDomain(iSub)[i] = domain.getSubDomain()[iSub]->getGlobSubNum();

  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, ownerSubDomain, maxOp);
  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, nodeMapping, maxOp);
  
  domain.getCommunicator()->fprintf(stdout,"hello3\n");
  fflush(stdout);

  
  std::set<int> globNodeIds;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    Connectivity& rSharedNodes(*refinedSharedNodes[iSub]);

    std::map<int,std::set<int> > toTransfer;
    std::map<int,std::map<int,int> > newlyInsertedNodes;
    std::map<int,int> newNodeIds;

    int * numNeighbors = new int[subD.getNumNeighb()];
    for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      for(int i = 0; i < rSharedNodes.num(jSub); ++i) {
        const int index = rSharedNodes[jSub][i];
        assert(nodeMapping(iSub)[index] >= 0);
        const int uniqueNodeID = ownerSubDomain(iSub)[index] * maxNodesPerSubD + nodeMapping(iSub)[index];
        if (uniqueNodeID == 3911) {
          std::cout << subD.getGlobSubNum() << " " << subD.getNeighb()[jSub] << std::endl;
        }
        if(refinedNodeDistInfo.getMasterFlag(iSub)[index]) {
          toTransfer[jSub].insert(nodeMapping(iSub)[index]);
          newlyInsertedNodes[jSub][uniqueNodeID] = nodeMapping(iSub)[index];
          //locToGlobMap[iSub][nodeMapping(iSub)[index]] = uniqueNodeID;
        } else /*if(ownerSubDomain(iSub)[index] == subD.getNeighb()[jSub])*/ {
          if(newNodeIds.find(uniqueNodeID) == newNodeIds.end())
            newNodeIds[uniqueNodeID] = numNodes[iSub]++;
          //locToGlobMap[iSub][newNodeIds[uniqueNodeID]] = uniqueNodeID;
          toTransfer[jSub].insert(newNodeIds[uniqueNodeID]);
          newlyInsertedNodes[jSub][uniqueNodeID] = newNodeIds[uniqueNodeID];
        }
      }
    }

    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {
      const int uniqueNodeID = ownerSubDomain(iSub)[i] * maxNodesPerSubD + nodeMapping(iSub)[i];
      globNodeIds.insert(uniqueNodeID);
      if(ownerSubDomain(iSub)[i] != subD.getGlobSubNum()) {
        nodeMapping(iSub)[i] = newNodeIds[ownerSubDomain(iSub)[i] * maxNodesPerSubD + nodeMapping(iSub)[i]];
        locToGlobMap[iSub][nodeMapping(iSub)[i]] = uniqueNodeID;
      }
      else
        locToGlobMap[iSub][nodeMapping(iSub)[i]] = uniqueNodeID;
    }

    for(std::map<int,std::set<int> >::const_iterator iter=toTransfer.begin(); iter != toTransfer.end(); ++iter) {
      numNeighbors[iter->first] = iter->second.size();
      nodeIdPattern->setLen(subD.getSndChannel()[iter->first],iter->second.size());
      nodeIdPattern->setLen(subD.getRcvChannel()[iter->first],iter->second.size());
      nodeVolPattern->setLen(subD.getSndChannel()[iter->first],iter->second.size());
      nodeVolPattern->setLen(subD.getRcvChannel()[iter->first],iter->second.size());
      nodePosnPattern->setLen(subD.getSndChannel()[iter->first],3*iter->second.size());
      nodePosnPattern->setLen(subD.getRcvChannel()[iter->first],3*iter->second.size());
      nodeVecPattern->setLen(subD.getSndChannel()[iter->first],dim*iter->second.size());
      nodeVecPattern->setLen(subD.getRcvChannel()[iter->first],dim*iter->second.size());
      matPattern->setLen(subD.getSndChannel()[iter->first],dim*dim*iter->second.size());
      matPattern->setLen(subD.getRcvChannel()[iter->first],dim*dim*iter->second.size());
    }

    sharedNodes[iSub] = new Connectivity(subD.getNumNeighb(),numNeighbors);
    int* j = new int[subD.getNumNeighb()];
    memset(j,0,sizeof(int)*subD.getNumNeighb());
    for(std::map<int,std::map<int,int> >::const_iterator neighbor = newlyInsertedNodes.begin(); neighbor != newlyInsertedNodes.end(); ++neighbor) {
      for(std::map<int,int>::const_iterator iter = neighbor->second.begin(); iter != neighbor->second.end(); ++iter) {
        assert(subD.getNeighb()[neighbor->first] != subD.getGlobSubNum());
        assert(j[neighbor->first] < numNeighbors[neighbor->first]);
        //if (locToGlobMap[iSub][iter->second] == 9)
        //  std::cout << "hello 3" << std::endl;
        //if (locToGlobMap[iSub][iter->second] == 58)
        //  std::cout << "hello 4" << std::endl;
        /*if (iSub == 1 && (iter->second == 429 || iter->second == 384))
          std::cout << "Hello2" <<std::endl;
        */
        (*sharedNodes[iSub])[neighbor->first][j[neighbor->first]++] = iter->second;
      }
    }
    for (int ii = 0; ii < subD.getNumNeighb(); ++ii)
      assert(j[ii] == numNeighbors[ii]);
  }
  nodeIdPattern->finalize();
  nodeVolPattern->finalize();
  nodeVecPattern->finalize();
  nodePosnPattern->finalize();
  matPattern->finalize();

  

#if 0
// #pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    fprintf(stderr,"SubDomain %d has... [%d nodes]\n", domain.getSubDomain()[iSub]->getGlobSubNum(), numNodes[iSub]);
    for(int jSub = 0; jSub < domain.getSubDomain()[iSub]->getNumNeighb(); ++jSub) {
      fprintf(stderr,"\tSubDomain %d as a neighbor with %d shared fine nodes; %d shared coarse nodes\n",
              domain.getSubDomain()[iSub]->getNeighb()[jSub],
              refinedSharedNodes[iSub]->num(jSub),
              sharedNodes[iSub]->num(jSub));
    }
  }
#endif
// ------------------------- END MPI PARALLEL CRAP -------------------------------------------------------------

  //std::cout << "gl node size " << globNodeIds.size() << std::endl;
  
  domain.getCommunicator()->fprintf(stdout,"hello2\n");
  fflush(stdout);


  std::set< std::pair<int,int> > globalEdges;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    // After agglomeration, allocate and build all the auxilliary data structures we'll be using
    int *numNodeNeighbors = new int[numNodes[iSub]];
    std::fill(numNodeNeighbors, numNodeNeighbors+numNodes[iSub], 1);

    int num_edges_skipped = 0;
    int (*edgePtr)[2] = refinedEdges[iSub]->getPtr();

    // Build the coarse EdgeSet, NodeToNode, and fine-to-coarse edge mapping
    edges[iSub] = new EdgeSet();
    int num_edges = 0;
    std::multimap<int,int> reverseEdgeMapping;
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      int i = edgePtr[l][0], j = edgePtr[l][1];
      int new_i = nodeMapping(iSub)[i], new_j = nodeMapping(iSub)[j];
#if 1 // Debug output
      if(new_i < 0 || new_j < 0 || new_i > numNodes[iSub] || new_j > numNodes[iSub]) {
          fprintf(stderr,"ERROR; on SubD %d, i = %d, new_i = %d, j = %d, new_j = %d, but total num nodes = %d!\n",
            domain.getSubDomain()[iSub]->getGlobSubNum(), i, new_i, j, new_j, numNodes[iSub]);
          assert(false);
      }
#endif
      if(new_i != new_j) {

        assert(locToGlobMap[iSub].find(new_i) != locToGlobMap[iSub].end());
        assert(locToGlobMap[iSub].find(new_j) != locToGlobMap[iSub].end()); 
        int i1 = locToGlobMap[iSub][new_i], j1 = locToGlobMap[iSub][new_j];
        if (i1 > j1) { 
          int tmp = i1;
          i1 = j1;
          j1 = tmp;
        }
        //if (i1 == 9 && j1 == 58)
        //  std::cout << "!!!" << std::endl; 
        globalEdges.insert( std::pair<int,int>( i1,j1));
        int new_edge = edges[iSub]->find(new_i, new_j);
        edgeMapping(iSub)[l] = new_edge;
        reverseEdgeMapping.insert(std::pair<int,int>(new_edge,l));
        
        if(num_edges < edges[iSub]->size()) {
          numNodeNeighbors[new_i]++;
          numNodeNeighbors[new_j]++;
          num_edges = edges[iSub]->size();
        }
      } else {
        ++num_edges_skipped;
        edgeMapping(iSub)[l] = -1;
      }
    }

    Vec<int> newNum(edges[iSub]->size());
    edges[iSub]->createPointers(newNum);
    edges[iSub]->setMasterFlag(new bool[edges[iSub]->size()]);

//    Vec<int> tmp(edgeMapping(iSub));
    // Edges are renumbered by create pointers!  Update our edge mapping
    // accordingly

    for (std::multimap<int,int>::iterator itr = reverseEdgeMapping.begin();
         itr != reverseEdgeMapping.end(); ++itr) {

      edgeMapping(iSub)[(*itr).second] = newNum[(*itr).first];
    }

    std::vector<Elem*> newElems;
    int newMap[4];
    bool ok = true;
    for (int i = 0; i < refinedElems[iSub]->size(); ++i) {

      Elem& elem((*refinedElems[iSub])[i]);
      switch (elem.type()) {
        
      case Elem::TET:
        ok = true;
        for (int k = 0; k < 4; ++k) {
          int new_node = nodeMapping(iSub)[elem[k]];
          for (int j = 0; j < k; ++j) {
            if (new_node == newMap[j])
              ok = false;
          }
          newMap[k] = new_node;
        }
        if (ok) {
          Elem* new_elem = new(refinedElems[iSub]->getBlockAllocator()) ElemTet;
          memcpy(new_elem->nodeNum(), newMap,sizeof(int)*4);
          for (int l = 0; l < 6; ++l) {
            new_elem->edgeNum(l) = edgeMapping(iSub)[elem.edgeNum(l)];
          }
          newElems.push_back(new_elem);
        }
      }
    }
    
    elems[iSub] = new ElemSet(newElems.size());
    for (int i = 0; i < newElems.size(); ++i)
      elems[iSub]->addElem(i, newElems[i]); 

    std::vector<Face*> newFaces;
    for (int i = 0; i < refinedFaces[iSub]->size(); ++i) {

      Face& face((*refinedFaces[iSub])[i]);
      switch (face.type()) {
        
      case Face::TRIA:
        ok = true;
        for (int k = 0; k < 3; ++k) {
          int new_node = nodeMapping(iSub)[face[k]];
          for (int j = 0; j < k; ++j) {
            if (new_node == newMap[j])
              ok = false;
          }
          newMap[k] = new_node;
        }
        if (ok) {
          Face* new_face = new(refinedFaces[iSub]->getBlockAllocator()) FaceTria;
          new_face->setup(face.getCode(),newMap, i, face.getSurfaceID());
          for (int l = 0; l < 3; ++l) {
            new_face->setEdgeNum(l, edgeMapping(iSub)[face.getEdgeNum(l)]);
          }
          newFaces.push_back(new_face);
        }
      }
    }
    
    faces[iSub] = new FaceSet(newFaces.size());
    for (int i = 0; i < newFaces.size(); ++i)
      faces[iSub]->addFace(i, newFaces[i]); 

    // Build the connectivity graph
    connectivity[iSub] = new Connectivity(numNodes[iSub], numNodeNeighbors);
    for(int i = 0; i < numNodes[iSub]; ++i) {
      (*connectivity[iSub])[i][0] = i;
      for(int n = 1; n < connectivity[iSub]->num(i);++n) (*connectivity[iSub])[i][n] = -1;
    }

    edgePtr = edges[iSub]->getPtr();
    for(int l = 0; l < edges[iSub]->size(); ++l) {
      int i=edgePtr[l][0], j=edgePtr[l][1];
      int index_i=0, index_j=0;
      while((*connectivity[iSub])[i][index_i] != -1) {
        ++index_i;
#if 1 // Debug output
        if(index_i >= connectivity[iSub]->num(i)) {
          fprintf(stderr,"ERROR: Trying to establish edge %d -> %d, but index_i = %d which is greater than %d\n",
                  i, j, index_i, connectivity[iSub]->num(i));
          assert(false);
        }
#endif
      }
      while((*connectivity[iSub])[j][index_j] != -1) {
          ++index_j;
#if 1 // Debug output
          if(index_j >= connectivity[iSub]->num(j)) {
          fprintf(stderr,"ERROR: Trying to establish edge %d -> %d, but index_j = %d which is greater than %d\n",
                  i, j, index_j, connectivity[iSub]->num(j));
          assert(false);
        }
#endif
      }
      (*connectivity[iSub])[i][index_i] = j;
      (*connectivity[iSub])[j][index_j] = i;
    }

    for (int l = 0; l < edges[iSub]->size(); ++l)
      edges[iSub]->getMasterFlag()[l] = true;

    // Create the "length" to be such that it's shorter between two agglomerated cells that share multiple edges
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      const int coarse_index = edgeMapping(iSub)[l];
      if(coarse_index < 0) continue;
      else {
        edges[iSub]->viewEdgeLength()[coarse_index] = std::max(refinedEdges[iSub]->length(l), edges[iSub]->viewEdgeLength()[coarse_index]);
//        edges[iSub]->getMasterFlag()[coarse_index] = refinedEdges[iSub]->getMasterFlag()[l] &&
//                                                     edges[iSub]->getMasterFlag()[coarse_index];
      }
    }

    /*fprintf(stderr, "\tAgglomeration finished with %d coarse cells (vs. %d fine cells) and %d edges vanished, %d coarse edges remaining (of %d)\n",
            numNodes[iSub], nodeMapping(iSub).size(), num_edges_skipped, edges[iSub]->size(), refinedEdges[iSub]->size());
*/
    nodeDistInfo->setLen(iSub, numNodes[iSub]);
    edgeDistInfo->setLen(iSub, edges[iSub]->size());
  }
//  fprintf(stderr, "\n");

  nodeDistInfo->finalize(true);
  //edgeDistInfo->finalize(true);

  //std::cout << "num Edges = " << globalEdges.size() << std::endl;
  
  sharedEdges = new EdgeDef**[numLocSub];
  numSharedEdges = new int*[numLocSub];
  
  CommPattern<int> edgeNumPat(domain.getSubTopo(), domain.getCommunicator(), CommPattern<int>::CopyOnSend,
			      CommPattern<int>::NonSym);
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    identifyEdges(edgeNumPat,iSub);

  edgeNumPat.finalize();
  
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    sndEdgeInfo(edgeNumPat,iSub);

  edgeNumPat.exchange();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    rcvEdgeInfo(edgeNumPat,iSub,dim);
  }

  edgeDistInfo->finalize(false);
/*
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    // We have to ensure that any agglomerated nodes are shared by the
    // correct subdomains (that they can talk to each other)
    std::map<int,std::set<int> > nodeSharedData;
    for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      for(int i = 0; i < sharedNodes[iSub]->num(jSub); ++i) {
        nodeSharedData[jSub].insert((*sharedNodes[iSub])[jSub][i]);
      }
    }

    int (*edgePtr)[2] = edges[iSub]->getPtr();
    for (int l = 0; l < edges[iSub]->size(); ++l) {

      int i = edgePtr[l][0],j = edgePtr[l][1];
      int gi = locToGlobMap[iSub][i], gj = locToGlobMap[iSub][j];
      if (gi == 9 && gj == 58 || gi == 58 && gj == 9)
        std::cout <<"ehllo" << std::endl;
      for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
        if (nodeSharedData[jSub].find(i) != nodeSharedData[jSub].end() &&
            nodeSharedData[jSub].find(j) != nodeSharedData[jSub].end() &&
            subD.getNeighb()[jSub] < subD.getGlobSubNum())
          edges[iSub]->getMasterFlag()[l] = false;
        
      }
    }
  }
*/
  std::set<std::pair<int,int> > new_edge_map;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    int (*edgePtr)[2] = edges[iSub]->getPtr();
    for (int l = 0; l < edges[iSub]->size(); ++l) {

      int i = edgePtr[l][0],j = edgePtr[l][1];
      int gi = locToGlobMap[iSub][i], gj = locToGlobMap[iSub][j];
      if (gi > gj) {

        int tmp = gi;
        gi = gj;
        gj = tmp;
      }
      assert(globalEdges.find(std::pair<int,int>(gi,gj)) != globalEdges.end());
      if (edges[iSub]->getMasterFlag()[l]) {

        assert(new_edge_map.find(std::pair<int,int>(gi,gj)) == new_edge_map.end());
        new_edge_map.insert(std::pair<int,int>(gi,gj));
      }
    }
  }
/*
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    int (*edgePtr)[2] = edges[iSub]->getPtr();
    for (int l = 0; l < edges[iSub]->size(); ++l) {
      int i = edgePtr[l][0],j = edgePtr[l][1];
      int gi = locToGlobMap[iSub][i], gj = locToGlobMap[iSub][j];
      if (gi > gj) {

        int tmp = gi;
        gi = gj;
        gj = tmp;
      }
      assert(new_edge_map.find(std::pair<int,int>(gi,gj)) != new_edge_map.end());
    }
  }
*/
  //std::cout << "New edge map size = " << new_edge_map.size() << std::endl;
 
/*
  std::vector<EdgeDef> newSharedEdges;  
  // now compute the shared edges for the new subdomain

  sharedEdges = new EdgeDef**[numLocSub];
  numSharedEdges = new int*[numLocSub];

  std::set<int> keptEdges;
  
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    
    Connectivity *nodeToNeighb = sharedNodes[iSub]->reverse();
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    sharedEdges[iSub] = new EdgeDef*[subD.getNumNeighb()];
    numSharedEdges[iSub] = new int[subD.getNumNeighb()];
  
    int (*edgePtr)[2] = edges[iSub]->getPtr(); 
    for (int l=0; l<edges[iSub]->size(); ++l) {

      int leftNode = edgePtr[l][0];
      int rightNode = edgePtr[l][1];
      if (nodeToNeighb->num(leftNode) == 0 ||
          nodeToNeighb->num(rightNode) == 0) continue;

      for (int iSub_ = 0; iSub_ < nodeToNeighb->num(leftNode); ++iSub_) {
        int subI = (*nodeToNeighb)[leftNode][iSub_];
        for (int jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub) {
          if (subI == (*nodeToNeighb)[rightNode][jSub]) {
            if (subD.getNeighb()[subI] < subD.getGlobSubNum())
              edges[iSub]->getMasterFlag()[l] = false;
          }
        }
      }
    }
    delete nodeToNeighb;

    int (*ptr)[2] = edges[iSub]->getPtr();
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
  
      newSharedEdges.clear();
      keptEdges.clear();
//      std::cout << jSub << std::endl;
      EdgeDef* refinedSharedEdgesLoc = refinedSharedEdges[iSub][jSub];
      int refinedNumSharedEdgesLoc = refinedNumSharedEdges[iSub][jSub];
      //std::cout << refinedNumSharedEdgesLoc << std::endl;
      for (int iEdge = 0; iEdge < refinedNumSharedEdgesLoc; ++iEdge) {

        EdgeDef& ed = refinedSharedEdgesLoc[iEdge];
        int emap = edgeMapping(iSub)[ed.edgeNum];
        if (emap >= 0) {

  //        if (subD.getGlobSubNum() > subD.getNeighb()[jSub])
  //          edges[iSub]->getMasterFlag()[emap] = false;        
        }
        if (emap >= 0 && keptEdges.find(emap) == keptEdges.end()) {
  //        EdgeDef ed;
          newSharedEdges.push_back(EdgeDef(ed.glLeft, ed.glRight, 
                                           emap,
                                           ed.sign));
          keptEdges.insert(emap);
        }
      }
      numSharedEdges[iSub][jSub] = newSharedEdges.size();
      //std::cout << iSub << " " << jSub <<  " " << numSharedEdges[iSub][jSub] << " " << std::endl;
      sharedEdges[iSub][jSub] = new EdgeDef[numSharedEdges[iSub][jSub]];
      memcpy(sharedEdges[iSub][jSub], &newSharedEdges[0], sizeof(EdgeDef)*newSharedEdges.size());
      offDiagMatPattern->setLen(subD.getSndChannel()[jSub],2*dim*dim*numSharedEdges[iSub][jSub]);
      offDiagMatPattern->setLen(subD.getRcvChannel()[jSub],2*dim*dim*numSharedEdges[iSub][jSub]);

    } 
    //std::cout << std::endl;
  }
  */
  offDiagMatPattern->finalize();

  finestNodeMapping = new DistVec<int>(*nodeDistInfo);

  // TODO(jontg): double-check this last part
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
   
    *nodeDistInfo->getMasterFlag(iSub) = true;
    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {

      if (nodeMapping(iSub)[i] >= 0) {
        if (finestNodeMapping_p1)
          (*finestNodeMapping)(iSub)[nodeMapping(iSub)[i]] = 
            (*finestNodeMapping_p1)(iSub)[i];
        else
          (*finestNodeMapping)(iSub)[nodeMapping(iSub)[i]] = i;

        // Nodes
        nodeDistInfo->getMasterFlag(iSub)[nodeMapping(iSub)[i]] = refinedNodeDistInfo.getMasterFlag(iSub)[i]/* && 
                                                                  nodeDistInfo->getMasterFlag(iSub)[nodeMapping(iSub)[i]]*/;
        nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]] = min(nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]],
                                                                     refinedNodeDistInfo.getInvWeight(iSub)[i]);
      
     // Edges
/*     edgeDistInfo->getMasterFlag(iSub)[edgeMapping(iSub)[i]] = refinedEdgeDistInfo.getMasterFlag(iSub)[i];
     edgeDistInfo->getInvWeight(iSub)[edgeMapping(iSub)[i]] = min(edgeDistInfo->getInvWeight(iSub)[edgeMapping(iSub)[i]],
                                                                   refinedEdgeDistInfo.getInvWeight(iSub)[i]);*/
      }
    }
  }

  distGeoState = new DistGeoState(refinedDistGeoState.getGeoData(),&domain, *nodeDistInfo, *edgeDistInfo);
  //distGeoState->getFaceNormal() = refinedDistGeoState.getFaceNormal();
  //distGeoState->getFaceNorVel() = refinedDistGeoState.getFaceNorVel();

  distGeoState->getFaceNorVel() = 0.0;

#if 1 // Debug output
  DistVec<bool> verifier(*nodeDistInfo);
  verifier = false;
// #pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < nodeMapping(iSub).size(); ++i)
      verifier(iSub)[nodeMapping(iSub)[i]] = true;

  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < verifier(iSub).size(); ++i)
      if(!verifier(iSub)[i]) fprintf(stderr, "Agglomerated cell %d on SubD %d has no refined nodes\n", i, iSub);
#endif
  delete []numNodes;

}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::computeRestrictedQuantities(const DistGeoState& refinedGeoState)
{
  distGeoState->getXn() *= 0.0;
  distGeoState->getCtrlVol() *= 0.0;
  const DistVec<double>& refinedVolume(refinedGeoState.getCtrlVol());
  const DistSVec<double, 3>& refinedX(refinedGeoState.getXn());

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < refinedVolume(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      distGeoState->getCtrlVol()(iSub)[coarseIndex] += refinedVolume(iSub)[i];
      for(int j = 0; j < 3; ++j)
        distGeoState->getXn()(iSub)[coarseIndex][j] += refinedVolume(iSub)[i] * refinedX(iSub)[i][j];
    }

    for(int i = 0; i < distGeoState->getXn()(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0/distGeoState->getCtrlVol()(iSub)[i];
      for(int j = 0; j < 3; ++j) distGeoState->getXn()(iSub)[i][j] *= one_over_volume;
    }
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < distGeoState->getXn()(iSub).size(); ++i)
      if(!nodeDistInfo->getMasterFlag(iSub)[i]) {
        for(int j = 0; j < 3; ++j) distGeoState->getXn()(iSub)[i][j] = 0.0;
        distGeoState->getCtrlVol()(iSub)[i] = 0.0;
      }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    (*distGeoState)(iSub).getEdgeNormal() = 0.0;
    for(int l = 0; l < edgeMapping(iSub).size(); ++l) {
      const int coarse_l = edgeMapping(iSub)[l];
      if(coarse_l < 0) continue;
      else {
        (*distGeoState)(iSub).getEdgeNormal()[coarse_l] += refinedGeoState(iSub).getEdgeNormal()[l];
      }
    }
  }
  operAdd<double> addOp;
  ::assemble(domain, *nodePosnPattern, sharedNodes, distGeoState->getXn(), addOp);
  ::assemble(domain, *nodeVolPattern, sharedNodes, distGeoState->getCtrlVol(), addOp);

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < faces[iSub]->size(); ++i) {

      (*faces[iSub])[i].computeNormal(distGeoState->getXn()(iSub), distGeoState->getFaceNormal()(iSub));
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2, dim>& fineData, DistSVec<Scalar2, dim>& coarseData) const
{
  coarseData *= 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i)  {
      for(int j = 0; j < dim; ++j)
        coarseData(iSub)[nodeMapping(iSub)[i]][j] += /*fineGrid.distGeoState->getCtrlVol()(iSub)[i] */ fineData(iSub)[i][j];
    }

    for(int i = 0; i < coarseData(iSub).size(); ++i) {
//      const Scalar one_over_volume = 1.0 / distGeoState->getCtrlVol()(iSub)[i];
      if (nodeDistInfo->getMasterFlag(iSub)[i]) {
//        for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] *= one_over_volume;
      } else {
        for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] = 0.0;
      }      
    }

  }

  operAdd<double> addOp;
  ::assemble(domain, *nodeVecPattern, const_cast<Connectivity **>(sharedNodes), coarseData, addOp);
  
  //std::cout << "Sum = " << sum << std::endl;
}

template<class Scalar> template<class Scalar2>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistVec<Scalar2>& fineData, DistVec<Scalar2>& coarseData) const
{
  coarseData *= 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i)
      coarseData(iSub)[nodeMapping(iSub)[i]] += /*fineGrid.distGeoState->getCtrlVol()(iSub)[i] */ fineData(iSub)[i];
/*
    for(int i = 0; i < coarseData(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0 / distGeoState->getCtrlVol()(iSub)[i];
      coarseData(iSub)[i] *= one_over_volume;
    } */
  }

}
//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::RestrictOperator(const MultiGridLevel<Scalar>& fineGrid,
                                              DistMat<Scalar2,dim>& fineOperator,
                                              DistMat<Scalar2,dim>& coarseOperator) {

  coarseOperator = 0.0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    bool* masterFlag = fineGrid.nodeDistInfo->getMasterFlag(iSub);
    for(int i = 0; i < fineGrid.nodeDistInfo->subSize(iSub); ++i) {

      if (!masterFlag[i])
        continue;

      int m = nodeMapping(iSub)[i];
      Scalar2* Aii = fineOperator(iSub).getElem_ii(i);
      Scalar2* Aii_coarse = coarseOperator(iSub).getElem_ii(m);
      Scalar2 rat = 1.0;/*fineGrid.distGeoState->getCtrlVol()(iSub)[i] / 
                   distGeoState->getCtrlVol()(iSub)[m];*/
      for (int k = 0; k < dim*dim; ++k) {
     
        Aii_coarse[k] += Aii[k] * rat * rat;
      }
    }

    int (*edgePtr)[2] = fineGrid.edges[iSub]->getPtr();
    int (*edgePtrCoarse)[2] = edges[iSub]->getPtr();

    masterFlag = fineGrid.edges[iSub]->getMasterFlag();    
    for (int l = 0; l < fineGrid.edges[iSub]->size(); ++l) {

      if (!masterFlag[l])
        continue;

      int i = edgePtr[l][0],j = edgePtr[l][1];
      int ci = nodeMapping(iSub)[i], cj = nodeMapping(iSub)[j];
      Scalar2* Aij = fineOperator(iSub).getElem_ij(l),
            * Aji = fineOperator(iSub).getElem_ji(l);
      Scalar2 rati = 1.0,/*fineGrid.distGeoState->getCtrlVol()(iSub)[i] /
                    distGeoState->getCtrlVol()(iSub)[ci],*/
             ratj = 1.0;/*fineGrid.distGeoState->getCtrlVol()(iSub)[j] /
                    distGeoState->getCtrlVol()(iSub)[cj];*/
      if (ci == cj) {
 
        Scalar2* Aii_coarse = coarseOperator(iSub).getElem_ii(ci);
        for (int k = 0; k < dim*dim; ++k) {

          Aii_coarse[k] += Aij[k]*rati*ratj+Aji[k]*ratj*rati;
        }
      } else {

        int cl = edgeMapping(iSub)[l];
        assert(cl >= 0);
        Scalar2* Aij_coarse = coarseOperator(iSub).getElem_ij(cl);
        Scalar2* Aji_coarse = coarseOperator(iSub).getElem_ji(cl);
        Scalar2* tmp;
        if (edgePtrCoarse[cl][0] != ci) {

          assert(edgePtrCoarse[cl][1] == ci);
          tmp = Aij_coarse; Aij_coarse = Aji_coarse; Aji_coarse = tmp;
        }
        for (int k = 0; k < dim*dim; ++k) { 
          Aij_coarse[k] += Aij[k]*rati*ratj;
          Aji_coarse[k] += Aji[k]*ratj*rati;
        }
      }
    }
    /*for (int i = 0; i < nodeDistInfo->subSize(iSub); ++i) {

      if (!nodeDistInfo->getMasterFlag(iSub)[i]) continue;

      Scalar2* Aii_coarse = coarseOperator(iSub).getElem_ii(i);
      for (int k = 0; k < dim*dim; ++k) {

        Aii_coarse[k] = (0.0);
      }
      for (int k = 0; k < dim; ++k) {

        Aii_coarse[k*dim+k] = 100.0;
      }
    }

    for (int l = 0; l < edges[iSub]->size(); ++l) {
      if (!edges[iSub]->getMasterFlag()[l]) continue;
      Scalar2* Aij_coarse = coarseOperator(iSub).getElem_ij(l);
      Scalar2* Aji_coarse = coarseOperator(iSub).getElem_ji(l);
      for (int k = 0; k < dim*dim; ++k) {

        Aij_coarse[k] = -1.0;
        Aji_coarse[k] = -1.0; 
      }
    }
  */
  }

  ::assemble(domain,*offDiagMatPattern, numSharedEdges,sharedEdges,coarseOperator);
  ::assemble(domain, *matPattern, sharedNodes, coarseOperator);
 /* 
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < nodeDistInfo->subSize(iSub); ++i) {

      if (locToGlobMap[iSub][i] == 3911) {

        std::cout << "isub = " << iSub << " i = " << i << std::endl;
        for (int k = 0; k <dim*dim; ++k) {

          if (k%dim == 0) std::cout << std::endl; 
          std::cout << coarseOperator(iSub).getElem_ii(i)[k] << " ";
        }
        std::cout << std::endl;
      }
    } 
  }
*/
}                                              

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Prolong(MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                     const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData) const
{

  //DistSVec<Scalar2,dim> correction(fineData);
  
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
/*      int gl = domain.getSubDomain()[iSub]->getNodeMap()[i];//locToGlobMap[iSub][coarseIndex];
      if (gl == 948) {

        std::cout << "Prolong: " << iSub << " " <<  i  << " " << fineData(iSub)[i][3] << " " << coarseData(iSub)[coarseIndex][3] << std::endl;
      }*/
    //  if (fineGrid.nodeDistInfo->getMasterFlag(iSub)[i])
        //for(int j = 0; j < dim; ++j) fineData(iSub)[i][j] += coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j];
        for(int j = 0; j < dim; ++j) {
          fineData(iSub)[i][j] += (coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j]);
    //      correction(iSub)[i][j] = (coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j]);
        }
/*fineGrid.distGeoState->getCtrlVol()(iSub)[i] / distGeoState->getCtrlVol()(iSub)[coarseIndex];*/
  //    else
  //      for(int j = 0; j < dim; ++j) fineData(iSub)[i][j] = 0.0;
    }
  }

//  fineGrid.assemble(fineData);
/*
  if (fineData(0).size() > 2000) {

    const_cast<MultiGridLevel<Scalar>*>(this)->writeXpostFile("correction.xpost",correction,0);
  }
*/
}

//------------------------------------------------------------------------------

template<class Scalar>
bool MultiGridLevel<Scalar>::isLine(int iSub,int edgei,int edgej,int* lineid,int* loci,int* locj)
{

  int idm1,idp1;
  idm1 = lineMap(iSub)[edgei][0];
  idp1 = lineMap(iSub)[edgei][1];

  if (edgei == edgej) {

    if (idm1 >= 0 || idp1 >= 0) {
      *lineid = lineIDMap(iSub)[edgei];
      *loci = lineLocIDMap(iSub)[edgei];
      *locj = *loci;
      return true;
    }
  }

  if (idm1 == edgej || idp1 == edgej) {
    *lineid = lineIDMap(iSub)[edgei];
    *loci = lineLocIDMap(iSub)[edgei];
    *locj = lineLocIDMap(iSub)[edgej];
    return true;
  } else {
    return false;
  }
  
}
//------------------------------------------------------------------------------
template<class Scalar>
template <class Scalar2,int dim,int neq> 
void MultiGridLevel<Scalar>::computeJacobian(DistSVec<Scalar2,dim>& U, DistSVec<Scalar2,dim>& V,
                                             DistVec<Scalar2>& irey,
                                             FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             FemEquationTerm* fet,
                                             DistMvpMatrix<Scalar2,neq>& matrices,
                                             DistTimeState<dim>* timeState) {

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    matrices(iSub) = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

//    matrices(iSub) = 0.0;
      
    if (fet) {

      elems[iSub]->computeJacobianGalerkinTerm(fet, (*distGeoState)(iSub), distGeoState->getXn()(iSub), distGeoState->getCtrlVol()(iSub), V(iSub), matrices(iSub));

      faces[iSub]->computeJacobianGalerkinTerm(*elems[iSub], fet, bcData(iSub), (*distGeoState)(iSub), distGeoState->getXn()(iSub), distGeoState->getCtrlVol()(iSub), V(iSub), matrices(iSub));
    }
  }

  ::assemble(domain,*offDiagMatPattern, numSharedEdges,sharedEdges,matrices);


#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    edges[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, (*distGeoState)(iSub),irey(iSub),
                                                 distGeoState->getXn()(iSub),
                                                 distGeoState->getCtrlVol()(iSub), 
                                                 V(iSub), matrices(iSub));
    faces[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn,bcData(iSub), (*distGeoState)(iSub),
                                                 V(iSub), matrices(iSub));
 
    Vec<double>& ctrlVol = distGeoState->getCtrlVol()(iSub); 
    for (int i=0; i<ctrlVol.size(); ++i) {
      Scalar voli = 1.0 / ctrlVol[i];
      Scalar2 *Aii = matrices(iSub).getElem_ii(i);
      for (int k=0; k<neq*neq; ++k)
        Aii[k] *= voli;
    }

/*
    for (int l=0; l < edges[iSub]->size(); ++l) {

      Scalar2 *Aij = matrices(iSub).getElem_ij(l);
      Scalar2 *Aji = matrices(iSub).getElem_ji(l);
      for (int k=0; k<neq*neq; ++k)
        Aij[k] = Aji[k] = 1.0;
    }
*/
  }

//  operAdd<double> addOp;

  ::assemble(domain, *matPattern, sharedNodes, matrices);  
  
  if (timeState) {

#pragma omp parallel for
    for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {
      Vec<double>& ctrlVol = distGeoState->getCtrlVol()(iSub); 
      for (int i=0; i<ctrlVol.size(); ++i) {
        (*timeState)(iSub).addToJacobianNoPrecLocal(i, ctrlVol[i],
                                            U(iSub),  matrices(iSub), 
                                            (*finestNodeMapping)(iSub)[i]);
      }
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar>   
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assemble(DistSVec<Scalar2,dim>& V)
{

  operAdd<double> addOp;
  ::assemble(domain, *nodeVecPattern, sharedNodes, V, addOp,&locToGlobMap);
}

template<class Scalar>   
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assembleMax(DistSVec<Scalar2,dim>& V)
{

//  operMax<double> maxOp;
//  ::assemble(domain, *nodeVecPattern, sharedNodes, V, maxOp);

  std::map<int,double> my_vals;
  int has_loc_to_glob_map = (locToGlobMap.size() > 0);
#pragma omp parallel for
  for (int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int i = 0; i < V(iSub).size(); ++i) {
    
      int gl_node = (has_loc_to_glob_map?locToGlobMap[iSub][i]:subD.getNodeMap()[i]);
      for (int j = 0; j < dim; ++j) {

        if (my_vals.find(gl_node*dim+j) != my_vals.end()) {
 
          double v1 = my_vals[gl_node*dim+j];
          double v2 = V(iSub)[i][j];
//          if (fabs(v1) >= 1.0e-10 || fabs(v2) >= 1.0e-10)
//            assert(v1 == v2 || fabs(v1-v2)/fabs(v1) < 1.0e-2);
        } else
          my_vals[gl_node*dim+j] = V(iSub)[i][j];
      }
    }
  }
}

template<class Scalar>
template <class Scalar2,int dim>
void MultiGridLevel<Scalar>::computeMatVecProd(DistMvpMatrix<Scalar2,dim>& mat,
                                               DistSVec<Scalar2,dim>& p,
                                               DistSVec<Scalar2,dim>& prod) {

  prod = 0.0;

/*  p = 1.0;

  std::cout << p*p << std::endl;
 */
  int i,j,l;

  int num_edges = 0;

#pragma omp parallel for 
  for (int iSub = 0; iSub < p.numLocSub(); ++iSub) {
    
    bool* nodeFlag = nodeDistInfo->getMasterFlag(iSub);
    bool* edgeFlag = edges[iSub]->getMasterFlag();

    GenMat<Scalar2,dim>& A = mat(iSub);
    Scalar2 (*a)[dim*dim] = A.data(); 
    int numNodes = nodeDistInfo->subSize(iSub);
    for (i=0; i<numNodes; ++i) {

      if (!nodeFlag[i])
        continue;

      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, i, p(iSub).v, i, prod(iSub).v, i);
    }

    int (*edgePtr)[2] = edges[iSub]->getPtr();
 
    for (l = 0; l < edgeDistInfo->subSize(iSub); ++l) {

      if (!edgeFlag[l])
        continue;

      ++num_edges;

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l, p(iSub).v, j, prod(iSub).v, i);
      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p(iSub).v, i, prod(iSub).v, j);
    }

  } 

  operAdd<double> addOp;
  ::assemble(domain, *nodeVecPattern, sharedNodes, prod, addOp,&locToGlobMap);

  //std::cout << "computed # edges = " << num_edges << std::endl;

//  std::cout << prod*prod << std::endl;
}

//------------------------------------------------------------------------------

#define INSTANTIATION_HELPER(T,dim) \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::Prolong(  MultiGridLevel<T> &, const DistSVec<double,dim> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::assemble(DistSVec<double,dim> &); \
  template void MultiGridLevel<T>::assembleMax(DistSVec<double,dim> &); \
  template void MultiGridLevel<T>::computeJacobian(DistSVec<double,dim>& U, DistSVec<double,dim>& V, \
                                             DistVec<double>& irey, \
                                             FluxFcn **fluxFcn, DistBcData<dim> &bcData, \
                                             FemEquationTerm*,DistMvpMatrix<double,dim> &A, \
                                             DistTimeState<dim>*); \
  template void MultiGridLevel<T>::computeMatVecProd(DistMvpMatrix<double,dim>& mat, \
                                                     DistSVec<double,dim>& p, \
                                                     DistSVec<double,dim>& prod); \
  template void MultiGridLevel<T>::RestrictOperator(const MultiGridLevel<T>& fineGrid, \
                                                    DistMat<double,dim>& fineOperator, \
                                                    DistMat<double,dim>& coarseOperator); \
  template void MultiGridLevel<T>::writeXpostFile(const std::string&, \
                                                    DistSVec<T,dim>&, \
                                                    int);

#define INSTANTIATION_HELPER2(T) \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistVec<double> &, DistVec<double> &) const; 

template class MultiGridLevel<double>;
INSTANTIATION_HELPER2(double);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
/*template class MultiGridLevel<float>;
INSTANTIATION_HELPER2(float);
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
INSTANTIATION_HELPER(float,7);
*/
