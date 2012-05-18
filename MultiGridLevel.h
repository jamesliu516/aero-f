#ifndef _MULTIGRID_LEVEL_H_
#define _MULTIGRID_LEVEL_H_

#include <DistInfo.h>
#include <DistVector.h>
#include <Domain.h>

class Connectivity;
class EdgeSet;
template<class Scalar, int dim> class DistSVec;

template<class Scalar>
class MultiGridLevel {
  private:
    Domain& domain;
    CommPattern<int> * nodeIdPattern;
    CommPattern<double> * nodeVolPattern;
    CommPattern<double> * nodeVecPattern;
    CommPattern<double> * nodePosnPattern;
    Connectivity ** sharedNodes;

  protected:
    DistInfo * nodeDistInfo;
    DistInfo * edgeDistInfo;

    int numLocSub;
    bool ownsData;
    Connectivity ** connectivity;
    EdgeSet ** edges;

    DistVec<int> nodeMapping;
    DistVec<int> edgeMapping;

    DistGeoState * distGeoState;

    DistSVec<int, 2> lineMap;

    DistVec<int> lineIDMap;
    
    DistVec<int> lineLocIDMap;

    int* numLines;

    int** lineLengths;

    std::vector<int>* lineids;

  public:
    MultiGridLevel(Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo);
    ~MultiGridLevel();

    DistInfo& getNodeDistInfo()       { return *nodeDistInfo; }
    DistInfo& getEdgeDistInfo()       { return *edgeDistInfo; }
    Connectivity ** getConnectivity() { return connectivity; }
    EdgeSet ** getEdges()             { return edges; }
    DistGeoState& getDistGeoState()   { return *distGeoState; }
    CommPattern<int>& getIdPat()      { return *nodeIdPattern; }
    Connectivity ** getSharedNodes()  { return sharedNodes; }

    void copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo, DistGeoState& refinedGeoState, Domain& domain);

    void agglomerate(const DistInfo& refinedNodeDistInfo,
                     CommPattern<int>& refinedNodeIdPattern,
                     DistGeoState& refinedDistGeoState,
                     Connectivity** refinedSharedNodes,
                     Connectivity ** nToN, EdgeSet ** edges,
                     Domain& domain,int dim);

    void computeRestrictedQuantities(const DistGeoState& refinedGeoState);

    template<class Scalar2, int dim> void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                                   const DistSVec<Scalar2, dim>& fineData,
                                                   DistSVec<Scalar2, dim>& coarseData) const;
    template<class Scalar2, int dim> void Prolong(const MultiGridLevel<Scalar>& coarseGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                                  const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData) const;

    bool isLine(int iSub,int edgei,int edgej, int* lineid, int* loci, int* locj);

    int NumLines(int iSub) const { return numLines[iSub]; }

    int* getLineData(int iSub, int lineid) { return &lineids[iSub][lineid*8]; } 
    int lineLength(int iSub, int lineid) { return lineLengths[iSub][lineid]; } 

    template <class Scalar2,int dim>
    void assemble(DistSVec<Scalar2,dim>& V);

    template <class Scalar2,int dim,int neq> 
    void computeJacobian(DistSVec<Scalar2,dim>& V,
                         DistVec<Scalar2>& irey,
                         FluxFcn **fluxFcn, BcData<dim> &bcData,
                         DistMat<Scalar2,neq> &A);

};

#endif
