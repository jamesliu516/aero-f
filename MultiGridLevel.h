#ifndef _MULTIGRID_LEVEL_H_
#define _MULTIGRID_LEVEL_H_

#include <DistInfo.h>
#include <DistVector.h>

class Connectivity;
class Domain;
class EdgeSet;
template<class Scalar, int dim> class DistSVec;

template<class Scalar>
class MultiGridLevel {
  protected:
    DistInfo * nodeDistInfo;
    DistInfo * edgeDistInfo;

    int numLocSub;
    bool ownsData;
    Connectivity ** connectivity;
    EdgeSet ** edges;

    DistVec<int> nodeMapping;
    DistVec<int> edgeMapping;

    DistSVec<Scalar, 3> * X;
    DistVec<Scalar> * volume;

    DistSVec<int, 2> lineMap;

    DistVec<int> lineIDMap;
    
    DistVec<int> lineLocIDMap;

    int* numLines;

    int** lineLengths;

    std::vector<int>* lineids;

  public:
    MultiGridLevel(DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo);
    ~MultiGridLevel();

    DistInfo& getNodeDistInfo()       { return *nodeDistInfo; }
    DistInfo& getEdgeDistInfo()       { return *edgeDistInfo; }
    Connectivity ** getConnectivity() { return connectivity; }
    EdgeSet ** getEdges()             { return edges; }
    DistSVec<Scalar,3>& getX()        { return *X; }
    DistVec<Scalar>& getVol()         { return *volume; }

    void copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo, const DistSVec<Scalar, 3>& Xn, Domain& dom);
    void setCtrlVolumes(const int iSub, Vec<Scalar>& ctrlVol);
    void agglomerate(Connectivity ** nToN, EdgeSet ** edges);

    void computeRestrictedQuantities(const DistVec<Scalar>& refinedVolume, const DistSVec<Scalar, 3>& refinedX);

    template<class Scalar2, int dim> void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                                   const DistSVec<Scalar2, dim>& fineData,
                                                   DistSVec<Scalar2, dim>& coarseData) const;
    template<class Scalar2, int dim> void Prolong(const MultiGridLevel<Scalar>& coarseGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                                  const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData) const;

    bool isLine(int iSub,int edgei,int edgej, int* lineid, int* loci, int* locj);

    int NumLines(int iSub) const { return numLines[iSub]; }

    int* getLineData(int iSub, int lineid) { return &lineids[iSub][lineid*8]; } 
    int lineLength(int iSub, int lineid) { return lineLengths[iSub][lineid]; } 
};

#endif
