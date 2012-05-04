#include <MultiGridLevel.h>

#include <Connectivity.h>
#include <Domain.h>
#include <Edge.h>

template<class Scalar>
MultiGridLevel<Scalar>::MultiGridLevel(DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo)
  : nodeDistInfo(new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com)),
    edgeDistInfo(new DistInfo(refinedEdgeDistInfo.numLocThreads, refinedEdgeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedEdgeDistInfo.locSubToGlobSub, refinedEdgeDistInfo.com)),
    numLocSub(refinedNodeDistInfo.numLocSub), ownsData(true), connectivity(new Connectivity*[numLocSub]),
    edges(new EdgeSet*[numLocSub]), nodeMapping(refinedNodeDistInfo), edgeMapping(refinedEdgeDistInfo), X(0), volume(0)
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = 0;
    edges[iSub] = 0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
MultiGridLevel<Scalar>::~MultiGridLevel()
{
  delete nodeDistInfo;
  delete edgeDistInfo;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub){
    if(ownsData) {
      delete connectivity[iSub];
      delete edges[iSub];
    }
  }
  delete []connectivity;
  delete []edges;

  delete X;
  delete volume;
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo,
                                              const DistSVec<Scalar, 3>& Xn, Domain& domain)
{
  ownsData = false;

  X = new DistSVec<Scalar, 3>(refinedNodeDistInfo);
  volume = new DistVec<Scalar>(refinedNodeDistInfo);
  *X = Xn;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = domain.getSubDomain()[iSub]->getNodeToNode();
    edges[iSub] = &domain.getSubDomain()[iSub]->getEdges();
    nodeDistInfo->setLen(iSub, refinedNodeDistInfo.subSize(iSub));
    edgeDistInfo->setLen(iSub, refinedEdgeDistInfo.subSize(iSub));
  }
  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::setCtrlVolumes(const int iSub, Vec<Scalar>& ctrlVol)
{
  (*volume)(iSub) = ctrlVol;
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::agglomerate(Connectivity ** nToN, EdgeSet ** refinedEdges)
{
  nodeMapping = -1;
  edgeMapping = -1;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    int seed_node = 0, agglom_id = 0;

    // Perform agglomoration
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

      for(std::vector<int>::const_iterator iter = agglom.begin(); iter != agglom.end(); iter++)
        nodeMapping(iSub)[*iter] = agglom_id;

      ++agglom_id;
      while(seed_node < nodeMapping(iSub).size() && nodeMapping(iSub)[seed_node] >= 0) ++seed_node;
    }

    // After agglomoration, allocate and build all the auxilliary data structures we'll be using
    int numNodes = agglom_id;
    int *numNodeNeighbors = new int[numNodes];
    std::fill(numNodeNeighbors, numNodeNeighbors+numNodes, 0);

    int num_edges_skipped = 0;
    int (*edgePtr)[2] = refinedEdges[iSub]->getPtr();

    // Build the coarse EdgeSet, NodeToNode, and fine-to-coarse edge mapping
    std::map<std::pair<int,int>,int> newEdges;
    edges[iSub] = new EdgeSet();
    int num_edges = 0;
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      int i = edgePtr[l][0], j = edgePtr[l][1];
      if(nodeMapping(iSub)[i] != nodeMapping(iSub)[j]) {
        int new_i = nodeMapping(iSub)[i], new_j = nodeMapping(iSub)[j],
            new_edge = edges[iSub]->find(new_i, new_j);
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

    // Build the connectivity graph
    connectivity[iSub] = new Connectivity(numNodes, numNodeNeighbors);
    for(int i = 0; i < numNodes; ++i) for(int n = 0; n < connectivity[iSub]->num(i);++n) (*connectivity[iSub])[i][n] = -1;

    edgePtr = edges[iSub]->getPtr();
    for(int l = 0; l < edges[iSub]->size(); ++l) {
      int i=edgePtr[l][0], j=edgePtr[l][1];
      int index_i=0, index_j=0;
      while((*connectivity[iSub])[i][index_i] != -1) ++index_i;
      while((*connectivity[iSub])[j][index_j] != -1) ++index_j;
      (*connectivity[iSub])[i][index_i] = j;
      (*connectivity[iSub])[j][index_j] = i;
    }

    // Create the "length" to be such that it's shorter between two agglomorated cells that share multiple edges
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      const int coarse_index = edgeMapping(iSub)[l];
      if(coarse_index < 0) continue;
      else edges[iSub]->viewEdgeLength()[coarse_index] = std::max(refinedEdges[iSub]->length(l), edges[iSub]->viewEdgeLength()[coarse_index]);
    }

    fprintf(stdout, "\tAgglomeration finished with %d coarse cells (vs. %d fine cells) and %d edges vanished, %d coarse edges remaining (of %d)\n",
            numNodes, nodeMapping(iSub).size(), num_edges_skipped, edges[iSub]->size(), refinedEdges[iSub]->size());

    nodeDistInfo->setLen(iSub, numNodes);
    edgeDistInfo->setLen(iSub, edges[iSub]->size());
  }

  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);

  X = new DistSVec<Scalar, 3>(*nodeDistInfo);
  volume = new DistVec<Scalar>(*nodeDistInfo);

#if 0 // Debug output
  DistVec<bool> verifier(*nodeDistInfo);
  verifier = false;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < nodeMapping(iSub).size(); ++i)
      verifier(iSub)[nodeMapping(iSub)[i]] = true;

  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < verifier(iSub).size(); ++i)
      if(!verifier(iSub)[i]) fprintf(stderr, "Agglomerated cell %d on SubD %d has no refined nodes\n", i, iSub);
#endif
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::computeRestrictedQuantities(const DistVec<Scalar>& refinedVolume, const DistSVec<Scalar, 3>& refinedX)
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < refinedVolume(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      (*volume)(iSub)[coarseIndex] += refinedVolume(iSub)[i];
      for(int j = 0; j < 3; ++j)
        (*X)(iSub)[coarseIndex][j] += refinedVolume(iSub)[i] * refinedX(iSub)[i][j];
    }

    for(int i = 0; i < (*X)(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0/(*volume)(iSub)[i];
      for(int j = 0; j < 3; ++j) (*X)(iSub)[i][j] *= one_over_volume;
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2, dim>& fineData, DistSVec<Scalar2, dim>& coarseData) const
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i) for(int j = 0; j < dim; ++j)
      coarseData(iSub)[fineGrid.nodeMapping(iSub)[i]][j] += (*fineGrid.volume)(iSub)[i] * fineData(iSub)[i][j];

    for(int i = 0; i < coarseData(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0 / (*volume)(iSub)[i];
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
#define INSTANTIATION_HELPER(T,dim) \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<float,dim>  &, DistSVec<float,dim>  &) const; \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::Prolong( const MultiGridLevel<T> &, const DistSVec<float,dim>  &, const DistSVec<float,dim>  &, DistSVec<float,dim>  &) const; \
  template void MultiGridLevel<T>::Prolong( const MultiGridLevel<T> &, const DistSVec<double,dim> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const;

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
