#ifndef _MULTIGRID_LEVEL_H_
#define _MULTIGRID_LEVEL_H_

#include <DistInfo.h>
#include <DistVector.h>
#include <Domain.h>
#include <MultiGridSmoothingMatrix.h>
#include <DistMvpMatrix.h>
#include <DistTimeState.h>
#include <SparseMatrix.h>
#include <MultigridCommon.h>

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
    CommPattern<double> * edgeAreaPattern;
    Connectivity ** sharedNodes;
    Connectivity** nodeToNodeMaskILU;
    int maxNodesPerSubD;
    mutable std::map<int,std::map<int,int> > locToGlobMap;

    MultiGridMethod mgMethod;

    DistGeoState* myGeoState;

    IoData* myIoData;
    
  protected:
    DistInfo * nodeDistInfo;
    DistInfo * edgeDistInfo;
    DistInfo* faceDistInfo;
    DistInfo* faceNormDistInfo;
    DistInfo* inletNodeDistInfo;

    int numLocSub;
    bool ownsData;
    Connectivity ** connectivity;
    EdgeSet ** edges;
    FaceSet** faces;
    EdgeDef*** sharedEdges;
    int** numSharedEdges;

    DistVec<int> nodeMapping;
    DistVec<int> edgeMapping;
    DistVec<Vec3D>* edgeNormals;

    DistSVec<int, 2> lineMap;
    DistSVec<double, 3>* Xn;

    DistVec<int>* finestNodeMapping;

    DistVec<int> lineIDMap;
    
    DistVec<int> lineLocIDMap;

    DistVec<double>* ctrlVol; 

    int* numLines;

    int** lineLengths;

    int* numNodes;

    std::vector<int>* lineids;

    void identifyEdges(CommPattern<int> &edgeNumPat, int mySub);

    void sndEdgeInfo(CommPattern<int> &edgeNumPat, int mySub);

    void rcvEdgeInfo(CommPattern<int> &edgeNumPat, int mySub,int dim);

    MultiGridLevel* parent;
  
    DistSVec<int,2>* fv_comp_tag;


  public:
    MultiGridLevel(MultiGridMethod,MultiGridLevel*,Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo);
    ~MultiGridLevel();

    MultiGridLevel* getParent() { return parent; }

    int mapFineToCoarse(int iSub,int i) {
      if (parent)
        return nodeMapping(iSub)[parent->mapFineToCoarse(iSub,i)];
      else
        return i;
    }

    MultiGridLevel* getFinestLevel() {

      if (parent)
        return parent->getFinestLevel();
      else
        return this;
    }   

    DistGeoState& getGeoState() const { return *myGeoState; }

    DistSVec<double,3>& getXn() const { return myGeoState->getXn(); }
 
    void mapNodeList(int iSub,std::tr1::unordered_set<int>&);

    enum SeedNodeChoice { Random, Mavripilis };

    DistInfo& getNodeDistInfo()       { return *nodeDistInfo; }
    DistInfo& getEdgeDistInfo()       { return *edgeDistInfo; }
    DistInfo& getFaceDistInfo()       { return *faceDistInfo; }
    DistInfo& getInletNodeDistInfo()       { return *inletNodeDistInfo; }
    Connectivity ** getConnectivity() { return connectivity; }
    EdgeSet ** getEdges()             { return edges; }
    FaceSet ** getFaces()             { return faces; }
    CommPattern<int>& getIdPat()      { return *nodeIdPattern; }
    Connectivity ** getSharedNodes()  { return sharedNodes; }

    DistVec<int>& getNodeMapping() { return nodeMapping; }

    DistVec<Vec3D>& getEdgeNormals() { return *edgeNormals; }
    DistVec<double>& getCtrlVol() const { return *ctrlVol; }

    EdgeDef*** getSharedEdges() { return sharedEdges; }
    int** getNumSharedEdges() { return numSharedEdges; }
    
    DistSVec<int,2>& getFVCompTag() const { return *fv_comp_tag; }

    void copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo, const DistInfo& refinedFaceDistInfo,const DistInfo& refinedInletNodeDistInfo,DistGeoState& refinedGeoState, Domain& domain);

    void agglomerate(const DistInfo& refinedNodeDistInfo,
                     const DistInfo& refinedEdgeDistInfo,
                     CommPattern<int>& refinedNodeIdPattern,
                     Connectivity** refinedSharedNodes,
                     Connectivity ** nToN, EdgeSet ** edges,
                     EdgeDef***, int**,
                     Domain& domain,int dim,
                     DistVec<Vec3D>& refinedEdgeNormals,
                     DistVec<double>& refinedVol,
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
    void assemble(DistMat<Scalar2,dim>& V);
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
    
    void writeXpostFile(const std::string& fileName,
                        DistVec<Scalar>& val);

    Connectivity* createEdgeBasedConnectivity(int iSub);
    compStruct* createRenumbering(int iSub,Connectivity *nodeToNode,
 		 		  int typeRenum, int print);

    template <int dim>
    SparseMat<Scalar,dim>* createMaskILU(int iSub,int fill, 
                                         int renum, int *ndType);

};

#endif
