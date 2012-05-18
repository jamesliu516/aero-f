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
    edges(new EdgeSet*[numLocSub]), nodeMapping(refinedNodeDistInfo), edgeMapping(refinedEdgeDistInfo), distGeoState(0), lineMap(refinedNodeDistInfo),lineIDMap(refinedNodeDistInfo),lineLocIDMap(refinedNodeDistInfo)
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = 0;
    edges[iSub] = 0;
  }
  lineIDMap = -1;

  numLines = new int[numLocSub];
  memset(numLines,0,sizeof(int)*numLocSub);
  lineids = new std::vector<int>[numLocSub];

  lineLengths = new int*[numLocSub];
 
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
    }
  }
  delete []sharedNodes;
  delete []connectivity;
  delete []edges;
  
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

  distGeoState = &refinedGeoState;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = domain.getSubDomain()[iSub]->getNodeToNode();
    sharedNodes[iSub] = domain.getSubDomain()[iSub]->getSharedNodes();
    edges[iSub] = &domain.getSubDomain()[iSub]->getEdges();
    nodeDistInfo->setLen(iSub, refinedNodeDistInfo.subSize(iSub));
    edgeDistInfo->setLen(iSub, refinedEdgeDistInfo.subSize(iSub));
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
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, DistSVec<T,dim> & data, const OperType & op) {
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

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        for (int j = 0; j < dim; ++j)
          data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][j] = OperType::apply(data(iSub)[ (*sharedNodes[iSub])[jSub][iNode]][j], buffer[iNode][j]);
    }
  }
}

template<class T,int dim>
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, GenMat<T,dim> & A) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j = 0; j < dim*dim; ++j) a[j] += buffer[iNode][j];
      }
    }
  }
}
};

template<class Scalar>
void MultiGridLevel<Scalar>::agglomerate(const DistInfo& refinedNodeDistInfo,
                                         CommPattern<int>& refinedNodeIdPattern,
                                         DistGeoState& refinedDistGeoState,
                                         Connectivity** refinedSharedNodes,
                                         Connectivity ** nToN, EdgeSet ** refinedEdges,
                                         Domain& domain,int dim)
{
  int * numNodes = new int[numLocSub];
  nodeIdPattern = new CommPattern<int>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<int>::CopyOnSend);
  nodeVolPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodeVecPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodePosnPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);

  nodeMapping = -1;
  edgeMapping = -1;

  lineMap = -1;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
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
      while(agglom.size() < 8) {
        const int current_node=agglom.back();
        for(int n = 0; n < nToN[iSub]->num(current_node); ++n) {
          const int neighbor = (*nToN[iSub])[current_node][n];
          if(current_node == neighbor || nodeMapping(iSub)[neighbor] >= 0) continue;
          if(std::find(agglom.begin(), agglom.end(), neighbor) == agglom.end())
            weights[neighbor] += 1.0/refinedEdges[iSub]->length(refinedEdges[iSub]->findOnly(current_node, neighbor));
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

// ------------------------- MPI PARALLEL CRAP -----------------------------------------------------------------
  DistVec<int> ownerSubDomain(refinedNodeDistInfo);
  operMax<int> maxOp;

  int maxNodesPerSubD = 0;
  for(int iSub = 0; iSub < numLocSub; ++iSub) maxNodesPerSubD = max(maxNodesPerSubD, nodeMapping(iSub).size());
  domain.getCommunicator()->globalMax(1, &maxNodesPerSubD);

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < ownerSubDomain(iSub).size(); ++i)
      if(!refinedNodeDistInfo.getMasterFlag(iSub)[i])
        nodeMapping(iSub)[i] = ownerSubDomain(iSub)[i] = -1;
      else
        ownerSubDomain(iSub)[i] = domain.getSubDomain()[iSub]->getGlobSubNum();

  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, ownerSubDomain, maxOp);
  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, nodeMapping, maxOp);

  sharedNodes = new Connectivity*[numLocSub];
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
        const int uniqueNodeID = ownerSubDomain(iSub)[index] * maxNodesPerSubD + nodeMapping(iSub)[index];
        if(refinedNodeDistInfo.getMasterFlag(iSub)[index]) {
          toTransfer[jSub].insert(nodeMapping(iSub)[index]);
          newlyInsertedNodes[jSub][uniqueNodeID] = nodeMapping(iSub)[index];
        } else if(ownerSubDomain(iSub)[index] == subD.getNeighb()[jSub]) {
          if(newNodeIds.find(uniqueNodeID) == newNodeIds.end())
            newNodeIds[uniqueNodeID] = numNodes[iSub]++;
          toTransfer[jSub].insert(newNodeIds[uniqueNodeID]);
          newlyInsertedNodes[jSub][uniqueNodeID] = newNodeIds[uniqueNodeID];
        }
      }
    }

    for(int i = 0; i < nodeMapping(iSub).size(); ++i)
      if(ownerSubDomain(iSub)[i] != subD.getGlobSubNum())
        nodeMapping(iSub)[i] = newNodeIds[ownerSubDomain(iSub)[i] * maxNodesPerSubD + nodeMapping(iSub)[i]];

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
    }

    sharedNodes[iSub] = new Connectivity(subD.getNumNeighb(),numNeighbors);
    for(std::map<int,std::map<int,int> >::const_iterator neighbor = newlyInsertedNodes.begin(); neighbor != newlyInsertedNodes.end(); ++neighbor) {
      int j = 0;
      for(std::map<int,int>::const_iterator iter = neighbor->second.begin(); iter != neighbor->second.end(); ++iter)
        (*sharedNodes[iSub])[neighbor->first][j++] = iter->second;
    }
  }
  nodeIdPattern->finalize();
  nodeVolPattern->finalize();
  nodeVecPattern->finalize();
  nodePosnPattern->finalize();

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
        int new_edge = edges[iSub]->find(new_i, new_j);
        edgeMapping(iSub)[l] = new_edge;
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

    // Create the "length" to be such that it's shorter between two agglomerated cells that share multiple edges
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      const int coarse_index = edgeMapping(iSub)[l];
      if(coarse_index < 0) continue;
      else {
        edges[iSub]->viewEdgeLength()[coarse_index] = std::max(refinedEdges[iSub]->length(l), edges[iSub]->viewEdgeLength()[coarse_index]);
        edges[iSub]->getMasterFlag()[coarse_index] = refinedEdges[iSub]->getMasterFlag()[l];
      }
    }

    fprintf(stderr, "\tAgglomeration finished with %d coarse cells (vs. %d fine cells) and %d edges vanished, %d coarse edges remaining (of %d)\n",
            numNodes[iSub], nodeMapping(iSub).size(), num_edges_skipped, edges[iSub]->size(), refinedEdges[iSub]->size());

    nodeDistInfo->setLen(iSub, numNodes[iSub]);
    edgeDistInfo->setLen(iSub, edges[iSub]->size());
  }
  fprintf(stderr, "\n");

  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);

  // TODO(jontg): double-check this last part
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {
      nodeDistInfo->getMasterFlag(iSub)[nodeMapping(iSub)[i]] = refinedNodeDistInfo.getMasterFlag(iSub)[i];
      nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]] = min(nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]],
                                                                   refinedNodeDistInfo.getInvWeight(iSub)[i]);
    }
  }

  distGeoState = new DistGeoState(refinedDistGeoState.getGeoData(),&domain, *nodeDistInfo, *edgeDistInfo);

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
      if(!nodeDistInfo->getMasterFlag(iSub)[i])
        for(int j = 0; j < 3; ++j) distGeoState->getXn()(iSub)[i][j] = 0.0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int l = 0; l < edgeMapping(iSub).size(); ++l) {
      const int coarse_l = edgeMapping(iSub)[l];
      if(coarse_l < 0) continue;
      else {
        (*distGeoState)(iSub).getEdgeNormal()[coarse_l] += refinedGeoState(iSub).getEdgeNormal()[l];
      }
    }
  operAdd<double> addOp;
  ::assemble(domain, *nodePosnPattern, sharedNodes, distGeoState->getXn(), addOp);
  ::assemble(domain, *nodeVolPattern, sharedNodes, distGeoState->getCtrlVol(), addOp);
}

//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2, dim>& fineData, DistSVec<Scalar2, dim>& coarseData) const
{
  coarseData *= 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i) for(int j = 0; j < dim; ++j)
      coarseData(iSub)[nodeMapping(iSub)[i]][j] += fineGrid.distGeoState->getCtrlVol()(iSub)[i] * fineData(iSub)[i][j];

    for(int i = 0; i < coarseData(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0 / distGeoState->getCtrlVol()(iSub)[i];
      for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] *= one_over_volume;
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Prolong(const MultiGridLevel<Scalar>& coarseGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                     const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData) const
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      for(int j = 0; j < dim; ++j) fineData(iSub)[i][j] += coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j];
    }
  }
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
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assemble(DistSVec<Scalar2,dim>& V)
{

  operAdd<double> addOp;
  ::assemble(domain, *nodeVecPattern, sharedNodes, V, addOp);
}

//------------------------------------------------------------------------------

#define INSTANTIATION_HELPER(T,dim) \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<float,dim>  &, DistSVec<float,dim>  &) const; \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::Prolong( const MultiGridLevel<T> &, const DistSVec<float,dim>  &, const DistSVec<float,dim>  &, DistSVec<float,dim>  &) const; \
  template void MultiGridLevel<T>::Prolong( const MultiGridLevel<T> &, const DistSVec<double,dim> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::assemble(DistSVec<double,dim> &);

template class MultiGridLevel<double>;
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
template class MultiGridLevel<float>;
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
INSTANTIATION_HELPER(float,7);
