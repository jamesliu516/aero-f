#include <GnatPreprocessing.C>

#define INSTANTIATION_HELPER(dim)\
template \
GnatPreprocessing<dim>::GnatPreprocessing(Communicator *_com, IoData &_ioData, Domain\
    &dom, DistGeoState *_geoState);\
\
template \
GnatPreprocessing<dim>::~GnatPreprocessing(); \
\
template \
void GnatPreprocessing<dim>::initialize(); \
\
template \
void GnatPreprocessing<dim>::buildReducedModel();\
\
template \
void GnatPreprocessing<dim>::constructApproximatedMetric();\
\
template \
void GnatPreprocessing<dim>::computeCorrelationMatrixEVD();\
\
template \
void GnatPreprocessing<dim>::computePseudoInverseMaskedSnapshots();\
\
template \
void GnatPreprocessing<dim>::computeMaskedSnapshots();\
\
template \
void GnatPreprocessing<dim>::computePseudoInverseTranspose();\
\
template \
void GnatPreprocessing<dim>::computeApproximatedMetricLowRankFactor();\
\
template \
void GnatPreprocessing<dim>::outputApproxMetricLowRankFactorFullCoords();\
\
template \
void GnatPreprocessing<dim>::outputApproxMetricLowRankFactorReducedCoords();\
\
template \
void GnatPreprocessing<dim>::testInnerProduct(char *snapshotType);\
\
template \
void GnatPreprocessing<dim>::setSampleNodes(int iCluster);\
\
template \
void GnatPreprocessing<dim>::setUpPodResJac(int iCluster);\
\
template \
void GnatPreprocessing<dim>::setUpBasisBasisProducts();\
\
template \
void GnatPreprocessing<dim>::computeQROfWeightedPhiJ();\
\
template \
void GnatPreprocessing<dim>::setUpGreedy(int iCluster);\
\
template \
void GnatPreprocessing<dim>::findMaxAndFillPodHat(const double myMaxNorm, const int\
    locSub, const int locNode, const int globalNode);\
\
template \
void GnatPreprocessing<dim>::determineSampleNodes();\
\
template \
void GnatPreprocessing<dim>::greedyIteration(int greedyIt);\
\
template \
void GnatPreprocessing<dim>::initializeLeastSquares();\
\
template \
void GnatPreprocessing<dim>::initializeLeastSquaresPseudoInv(int numRhs);\
\
template \
void GnatPreprocessing<dim>::makeNodeMaxIfUnique(double nodeError, double\
    &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int\
    &globalNode);\
\
template \
void GnatPreprocessing<dim>::computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError);\
\
template \
void GnatPreprocessing<dim>::getSubDomainError(int iSub);\
\
template \
void GnatPreprocessing<dim>::leastSquaresReconstruction();\
\
template void GnatPreprocessing<dim>::subDFindMaxError(int iSub, bool\
    onlyInletOutletBC, double &myMaxNorm, int &locSub, int &locNode, int\
    &globalNode);\
\
template \
void GnatPreprocessing<dim>::parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff);\
\
template \
void GnatPreprocessing<dim>::buildRemainingMesh();\
\
template \
void GnatPreprocessing<dim>::addSampleNodesAndNeighbors();\
\
template \
void GnatPreprocessing<dim>::addNeighbors(int iIslands, int startingNodeWithNeigh);\
\
template \
void GnatPreprocessing<dim>::computeBCFaces(bool liftContribution);\
\
template \
void GnatPreprocessing<dim>::communicateAll();\
\
template \
void GnatPreprocessing<dim>::defineMaps();\
\
template \
void GnatPreprocessing<dim>::communicateBCFaces();\
  \
template \
bool GnatPreprocessing<dim>::checkFaceInMesh(FaceSet& currentFaces, const int iFace, const int iSub, const int *locToGlobNodeMap);\
\
template \
bool GnatPreprocessing<dim>::checkFaceAlreadyAdded(const int cpuNum, const int\
    iSub, const int iFace);\
  \
template \
void GnatPreprocessing<dim>::addFaceNodesElements(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
    *locToGlobNodeMap);\
\
template \
void GnatPreprocessing<dim>::addNodesOnFace(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
		*locToGlobNodeMap, int *locNodeNums);\
\
template \
void GnatPreprocessing<dim>::addElementOfFace(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
    *locToGlobNodeMap, const int *locNodeNums);\
\
template \
bool GnatPreprocessing<dim>::checkFaceContributesToLift(FaceSet& faces, const int iFace, const int iSub, const int *locToGlobNodeMap );\
\
template \
void GnatPreprocessing<dim>::outputTopFile(int iCluster);\
\
template \
void GnatPreprocessing<dim>::outputSampleNodes(int iCluster);\
\
template \
void GnatPreprocessing<dim>::outputSampleNodesGeneral(const std::vector<int> &sampleNodes, const char *outSampleNodeFile);\
\
template \
void GnatPreprocessing<dim>::computeXYZ(int iSub, int iLocNode, double *xyz);\
\
template \
void GnatPreprocessing<dim>::computePseudoInverse(int iPodBasis);\
\
template \
void GnatPreprocessing<dim>::computePseudoInverse();\
\
template \
void GnatPreprocessing<dim>::computePodTPod();\
\
template \
void GnatPreprocessing<dim>::assembleOnlineMatrices();\
\
template \
void GnatPreprocessing<dim>::outputOnlineMatrices(int iCluster);\
\
template \
void GnatPreprocessing<dim>::outputOnlineMatricesGeneral(int iCluster, int numNodes,\
    const std::map<int,int> &sampleNodeMap, const std::vector<int>\
    &sampleNodeVec);\
\
template \
void GnatPreprocessing<dim>::outputLocalReferenceStateReduced(int iCluster);\
\
template \
void GnatPreprocessing<dim>::outputInitialConditionReduced();\
\
template \
void GnatPreprocessing<dim>::outputClusterCentersReduced();\
\
template \
void GnatPreprocessing<dim>::outputLocalStateBasisReduced(int iCluster);\
\
template \
void GnatPreprocessing<dim>::outputReducedSVec(const DistSVec<double,dim>\
    &distSVec, FILE* outFile , double tag);\
\
template \
void GnatPreprocessing<dim>::outputWallDistanceReduced();\
\
template \
void GnatPreprocessing<dim>::outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector);\
\
template \
void GnatPreprocessing<dim>::outputReducedToFullNodes();\
\
template \
void GnatPreprocessing<dim>::checkConsistency();\
\
template \
void GnatPreprocessing<dim>::formMaskedNonlinearROBs();\
\
template \
void GnatPreprocessing<dim>::formReducedSampleNodeMap();

INSTANTIATION_HELPER(5);
INSTANTIATION_HELPER(6);

#undef INSTANTIATION_HELPER
