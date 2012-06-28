#ifndef _MULTIGRID_LEVEL_H_
#define _MULTIGRID_LEVEL_H_

#include <DistInfo.h>
#include <DistVector.h>
#include <Domain.h>
#include <MultiGridSmoothingMatrix.h>
#include <DistMvpMatrix.h>
#include <DistTimeState.h>

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
    CommPattern<double> * matPattern;
    CommPattern<double> * offDiagMatPattern;
    Connectivity ** sharedNodes;

    int maxNodesPerSubD;
    mutable std::map<int,std::map<int,int> > locToGlobMap;
  protected:
    DistInfo * nodeDistInfo;
    DistInfo * edgeDistInfo;

    int numLocSub;
    bool ownsData;
    Connectivity ** connectivity;
    EdgeSet ** edges;
    EdgeDef*** sharedEdges;
    int** numSharedEdges;

    DistVec<int> nodeMapping;
    DistVec<int> edgeMapping;

    DistSVec<int, 2> lineMap;

    DistVec<int> lineIDMap;
    
    DistVec<int> lineLocIDMap;

    int* numLines;

    int** lineLengths;

    std::vector<int>* lineids;

    void identifyEdges(CommPattern<int> &edgeNumPat, int mySub);

    void sndEdgeInfo(CommPattern<int> &edgeNumPat, int mySub);

    void rcvEdgeInfo(CommPattern<int> &edgeNumPat, int mySub,int dim);

  public:
    MultiGridLevel(Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo);
    ~MultiGridLevel();

    DistInfo& getNodeDistInfo()       { return *nodeDistInfo; }
    DistInfo& getEdgeDistInfo()       { return *edgeDistInfo; }
    Connectivity ** getConnectivity() { return connectivity; }
    EdgeSet ** getEdges()             { return edges; }
    CommPattern<int>& getIdPat()      { return *nodeIdPattern; }
    Connectivity ** getSharedNodes()  { return sharedNodes; }

    EdgeDef*** getSharedEdges() { return sharedEdges; }
    int** getNumSharedEdges() { return numSharedEdges; }

    void copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo, DistGeoState& refinedGeoState, Domain& domain);

    void agglomerate(const DistInfo& refinedNodeDistInfo,
                     const DistInfo& refinedEdgeDistInfo,
                     CommPattern<int>& refinedNodeIdPattern,
                     Connectivity** refinedSharedNodes,
                     Connectivity ** nToN, EdgeSet ** edges,
                     EdgeDef***, int**,
                     Domain& domain,int dim,
                     DistVec<int>*);

    void computeRestrictedQuantities(const DistGeoState& refinedGeoState);

    template<class Scalar2, int dim> void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                                   const DistSVec<Scalar2, dim>& fineData,
                                                   DistSVec<Scalar2, dim>& coarseData) const;
    template<class Scalar2>void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                         const DistVec<Scalar2>& fineData,
                                         DistVec<Scalar2>& coarseData) const;
    template<class Scalar2, int dim> void Prolong(MultiGridLevel<Scalar>& coarseGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                                  const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData) const;

    template<class Scalar2, int dim>
    void RestrictOperator(const MultiGridLevel<Scalar>& fineGrid,
                          DistMat<Scalar2,dim>& fineOperator,
                          DistMat<Scalar2,dim>& coarseOperator);

    bool isLine(int iSub,int edgei,int edgej, int* lineid, int* loci, int* locj);

    int NumLines(int iSub) const { return numLines[iSub]; }

    int* getLineData(int iSub, int lineid) { return &lineids[iSub][lineid*8]; } 
    int lineLength(int iSub, int lineid) { return lineLengths[iSub][lineid]; } 

    template <class Scalar2,int dim>
    void assemble(DistSVec<Scalar2,dim>& V);
    template <class Scalar2,int dim>
    void assembleMax(DistSVec<Scalar2,dim>& V);

    template <class Scalar2,int dim>
    void computeMatVecProd(DistMat<Scalar2,dim>& mat,
                           DistSVec<Scalar2,dim>& p,
                           DistSVec<Scalar2,dim>& prod);

    void WriteTopFile(const std::string& fileName);

    template <int dim>
    void writeXpostFile(const std::string& fileName,
                        DistSVec<Scalar,dim>& val,int id);
};

#endif
