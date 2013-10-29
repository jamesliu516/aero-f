#include <NonlinearRomDatabaseConstruction.h>
#include <Modal.h>
#include <TsInput.h>
#include <cmath>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

using std::stable_sort;

template<int dim> 
NonlinearRomDatabaseConstruction<dim>::NonlinearRomDatabaseConstruction(Communicator* _com, IoData& _ioData, Domain& _domain, GeoSource& _geoSource)  : 
  NonlinearRomOnlineII<dim>(_com, _ioData, _domain), geoSource(_geoSource)
{ 
  // ioData->example, com->example, this->domain.example
  nSnapshotFiles = 0;
  stateSnapshotTags = NULL;
  stateSnapsFromFile = NULL;   
  stateSnapshotClustersAfterOverlap = NULL;

  projErrorLog = NULL;
  initialCondition = NULL; 

  arbitraryUniformIC = false;

  ioDataProjError = &(this->ioData->romOffline.rob.relativeProjectionError);
}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRomDatabaseConstruction<dim>::~NonlinearRomDatabaseConstruction() 
{
  if (stateSnapshotClustersAfterOverlap) {
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      for (int iSnap=0;iSnap<stateSnapsFromFile[iFile];++iSnap) {
        delete [] stateSnapshotClustersAfterOverlap[iFile][iSnap];
      }
    }
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      delete [] stateSnapshotClustersAfterOverlap[iFile];
    }
    delete [] stateSnapshotClustersAfterOverlap;
  }
  stateSnapshotClustersAfterOverlap = NULL;

  if (stateSnapshotTags) {
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      delete [] stateSnapshotTags[iFile];
    }
    delete [] stateSnapshotTags;
  }
  stateSnapshotTags = NULL;

  if (stateSnapsFromFile)  delete [] stateSnapsFromFile;
  stateSnapsFromFile = NULL;

  if (initialCondition) delete initialCondition;
  initialCondition = NULL;
}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::constructDatabase() {

// The functions kmeans(), localPod(), and localRelProjError() are called here, depending on inputs.

  // clustering
  if (this->ioData->romOffline.rob.clustering.useExistingClusters == ClusteringData::USE_EXISTING_CLUSTERS_FALSE) {
    readSnapshotFile("state", false);
    kmeans();
    this->outputClusteredSnapshots("state");

    // for snapshot collection method 0 (all snapshots from FOM)
    if (strcmp(this->ioData->input.residualSnapFile,"")!=0) this->placeNonStateSnapshotsInClusters("residual");
    if (strcmp(this->ioData->input.krylovSnapFile,"")!=0) this->placeNonStateSnapshotsInClusters("krylov");

    if (strcmp(this->ioData->input.sensitivitySnapFile,"")!=0 && 
         strcmp(this->sensitivityClusterName,"")!=0) {
       readSnapshotFile("sensitivity", false);
       this->outputClusteredSnapshots("sensitivity");
    } 
  }

  // local POD
  if (this->ioData->romOffline.rob.state.dataCompression.computePOD) localPod("state");
  if (this->ioData->romOffline.rob.krylov.dataCompression.computePOD) localPod("krylov");
  if (this->ioData->romOffline.rob.sensitivity.dataCompression.computePOD) localPod("sensitivity");
  if (this->ioData->romOffline.rob.residual.dataCompression.computePOD) localPod("residual");
  if (this->ioData->romOffline.rob.jacAction.dataCompression.computePOD) localPod("jacAction");

  // preprocessing for fast distance calculations
  if (this->ioData->romOffline.rob.distanceComparisons.preprocessForNoUpdates || 
      this->ioData->romOffline.rob.distanceComparisons.preprocessForExactUpdates || 
      this->ioData->romOffline.rob.distanceComparisons.preprocessForApproxUpdates) preprocessForDistanceComparisons();

  // preprocessing for basis updates
  // exact
  // approx

  // projection error
  if (ioDataProjError->relProjError!=RelativeProjectionErrorData::REL_PROJ_ERROR_OFF) localRelProjError();

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::readSnapshotFile(char* snapType, bool preprocess) {

  // Check for snapshot command file

  char *vecFile;
  bool typeIsState = false;

  if (strcmp(snapType, "state")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.stateSnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.stateSnapFile);
    typeIsState = true;
  } else if (strcmp(snapType,"sensitivity")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.sensitivitySnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.sensitivitySnapFile);
//  } else if (strcmp(snapType,"residual")==0) {
//    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.residualSnapFile) + 1];
//    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.residualSnapFile);
//  } else if (strcmp(snapType,"krylov")==0) {
//    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.krylovSnapFile) + 1];
//    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.krylovSnapFile);
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", snapType);
    exit (-1);
  }

  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    this->com->fprintf(stderr, "*** Error: No snapshots FILES in %s\n", vecFile);
    exit (-1);
  }

  delete [] vecFile;
  vecFile = NULL;

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  this->com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

  if (typeIsState) {
    nSnapshotFiles=nData;
    stateSnapsFromFile = new int[nSnapshotFiles];
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) 
      stateSnapsFromFile[iFile] = 0;
  }

  char** snapFile = new char*[nData];
  for (int iData = 0; iData < nData; ++iData)
    snapFile[iData] = new char[500];
  char snapFile1[500];
  int* numSnaps = new int[nData];
  int* startSnaps = new int[nData];
  int* endSnaps = new int[nData];
  int* sampleFreq = new int[nData];
  double* snapWeight = new double[nData];
  int nSnap, iStart, iEnd, iFreq;
  double weight;

  if (typeIsState && this->ioData->romOffline.rob.state.snapshots.incrementalSnaps)
    this->com->fprintf(stderr, "*** Warning: Incremental snapshots is not supported for multiple bases (yet) \n");

  if (typeIsState && (this->ioData->romOffline.rob.state.dataCompression.energyOnly == DataCompressionData::ENERGY_ONLY_TRUE))
    this->com->fprintf(stderr, "*** Warning: EnergyOnly is not supported for multiple bases\n");

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %lf", snapFile1,&iStart,&iEnd,&iFreq,&weight);
    if (iStart < 1) iStart = 1;
    if (iEnd < 0) iEnd = 0;
    if (iFreq < 1) iFreq = 1;
    //numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart - 1;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    this->com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);
  }

  // compute the total number of snapshots
  int nTotSnaps = 0;
  int dummyStep = 0;
  double dummyTag = 0.0;
  for (int iData = 0; iData < nData; ++iData) {
    bool status = this->domain.template readTagFromFile<double, dim>(snapFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));

    if (!status) {
      this->com->fprintf(stderr, "*** Error: could not read snapshotsfrom %s \n", snapFile[iData]);
      exit(-1);
    }

    if ((endSnaps[iData]==0) || (endSnaps[iData]>numSnaps[iData]))
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
        if (typeIsState) ++stateSnapsFromFile[iData]; 
      }
    }
  }

  bool incrementalSnaps = false;
  bool subtractRefSol = false;
  if (preprocess) {
    if (ioDataProjError->projectIncrementalSnaps) {
      incrementalSnaps = true;
      nTotSnaps -= nData;
    } else if (ioDataProjError->subtractRefSol) {
      subtractRefSol = true;
      if (!(this->ioData->input.stateSnapRefSolution)) {
        this->com->fprintf(stderr, "*** Error: Reference solution not found \n");
        exit (-1);
      }
      this->readReferenceState();
    }
  }

  this->snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufOld = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufNew = new DistSVec<double, dim>(this->domain.getNodeDistInfo());

  if (typeIsState) {
    stateSnapshotTags = new double*[nSnapshotFiles];
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      stateSnapshotTags[iFile] = new double[stateSnapsFromFile[iFile]];
      for (int iSnap=0;iSnap<stateSnapsFromFile[iFile];++iSnap)
        stateSnapshotTags[iFile][iSnap] = -1.0;
    }
  }

  *snapBufOld = 0.0;
  *snapBufNew = 0.0;

  double* tags = new double [nTotSnaps];
  double tagOld;
  double tagNew;

  int numCurrentSnapshots = 0;
  bool status;

  if (subtractRefSol) {
    *snapBufOld = *(this->snapRefState);
    delete (this->snapRefState);
    (this->snapRefState) = NULL;
  }

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        if ((iSnap == startSnaps[iData]) && incrementalSnaps) {
          status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tagOld, *snapBufOld);
        } else {
          status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tagNew, *snapBufNew);
          tags[numCurrentSnapshots] = tagNew;
          (*(this->snap))[numCurrentSnapshots] = *snapBufNew - *snapBufOld;  //snapBufOld = 0 if not using incremental snaps
          if (incrementalSnaps) *snapBufOld = *snapBufNew;
          if (snapWeight[iData]) (*(this->snap))[numCurrentSnapshots] *= snapWeight[iData]; //CBM--check
          ++numCurrentSnapshots;
        }
      }
    }
  }

  if (typeIsState) {
    int snapCount=0;
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      for (int iSnap=0;iSnap<stateSnapsFromFile[iFile];++iSnap) {
        stateSnapshotTags[iFile][iSnap] = tags[snapCount];
        ++snapCount;
      }
    }
  }

  delete snapBufOld;
  snapBufOld = NULL;
  delete snapBufNew;
  snapBufNew = NULL;

  for (int iData=0; iData < nData; ++iData) {
    delete [] snapFile[iData];
  }
  
  delete [] snapFile;
  snapFile = NULL;
  delete [] numSnaps;
  numSnaps = NULL;
  delete [] startSnaps;
  startSnaps = NULL;
  delete [] endSnaps;
  endSnaps = NULL;
  delete [] sampleFreq;
  sampleFreq = NULL;
  delete [] snapWeight;
  snapWeight = NULL;
  delete [] tags;
  tags = NULL;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::placeNonStateSnapshotsInClusters(char* snapType) {

  char *vecFile;

  if (strcmp(snapType,"residual")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.residualSnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.residualSnapFile);
  } else if (strcmp(snapType,"krylov")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.krylovSnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.krylovSnapFile);
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", snapType);
    exit (-1);
  }

  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    this->com->fprintf(stderr, "*** Error: No snapshots FILES in %s\n", vecFile);
    exit (-1);
  }

  delete [] vecFile;
  vecFile = NULL;

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  this->com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

  if (nData!=nSnapshotFiles) {
    this->com->fprintf(stderr, "*** Error: incorrect number of files listed in %s (these files must correspond to the state snapshot files)\n", vecFile);
    exit(-1);
  }

  char** snapFile = new char*[nData];
  for (int iData = 0; iData < nData; ++iData)
    snapFile[iData] = new char[500];
  char snapFile1[500];
  int* numSnaps = new int[nData];
  int* startSnaps = new int[nData];
  int* endSnaps = new int[nData];
  int* sampleFreq = new int[nData];
  double* snapWeight = new double[nData];
  int nSnap, iStart, iEnd, iFreq;
  double weight;

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %lf", snapFile1,&iStart,&iEnd,&iFreq,&weight);
    if (iStart < 1) iStart = 1;
    if (iEnd < 0) iEnd = 0;
    if (iFreq < 1) iFreq = 1;
    //numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart - 1;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    this->com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);
  }

  // compute the total number of snapshots
  int nTotSnaps = 0;
  int dummyStep = 0;
  double dummyTag = 0.0;
  for (int iData = 0; iData < nData; ++iData) {
    bool status = this->domain.template readTagFromFile<double, dim>(snapFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));
    if ((endSnaps[iData]==0) || (endSnaps[iData]>numSnaps[iData]))
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
      }
    }
  }

  DistSVec<double, dim>* snapBuf = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
  *snapBuf = 0.0;
  double tag;
  bool status;

  this->snapsInCluster = new int[this->nClusters]; 
  for (int iCluster=0;iCluster<(this->nClusters);++iCluster)
    this->snapsInCluster[iCluster] = 0;

  this->initializeClusteredOutputs();

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tag, *snapBuf);        
        if (snapWeight[iData]) *snapBuf *= snapWeight[iData]; //CBM--check
        // find associated state
        int stateSnap = -1;
        for (int iStateSnap=0; iStateSnap<stateSnapsFromFile[iData]; ++iStateSnap) {
          if (iStateSnap==(stateSnapsFromFile[iData]-1)) {
            stateSnap = iStateSnap;
          } else if (tag<stateSnapshotTags[iData][iStateSnap]) {
            // do nothing
          } else if ((tag>=stateSnapshotTags[iData][iStateSnap]) && (tag<stateSnapshotTags[iData][iStateSnap+1])) {
            stateSnap = iStateSnap;
            break;
          }
        }
        // store snapshot in all applicaple clusters
        if (strcmp(snapType,"residual")==0) {
          for (int iCluster=0;iCluster<(this->nClusters);++iCluster) {
            if (stateSnapshotClustersAfterOverlap[iData][stateSnap][iCluster]) { 
              this->com->fprintf(stdout, "... matched with state snapshot #%d from file #%d; writing to cluster %d\n",stateSnap,iData,iCluster); 
              this->writeClusteredBinaryVectors(iCluster, snapBuf, NULL, NULL);
              ++(this->snapsInCluster[iCluster]);
            }
          }
        } else if (strcmp(snapType,"krylov")==0) {
          for (int iCluster=0;iCluster<(this->nClusters);++iCluster) {
            if (stateSnapshotClustersAfterOverlap[iData][stateSnap][iCluster]) {
              this->com->fprintf(stdout, "... matched with state snapshot #%d from file #%d; writing to cluster %d\n",stateSnap,iData,iCluster);
              this->writeClusteredBinaryVectors(iCluster, NULL, NULL, snapBuf);
              ++(this->snapsInCluster[iCluster]);
            }
          }
        }
      }
    }
  }

  this->com->fprintf(stdout, "\n"); 
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    this->com->fprintf(stdout, "... %d %s vectors placed in cluster %d\n",this->snapsInCluster[iCluster], snapType, iCluster);
  }
  this->com->fprintf(stdout, "\n");


  delete [] (this->snapsInCluster);
  this->snapsInCluster = NULL;

  delete snapBuf;
  snapBuf = NULL;

  for (int iData=0; iData < nData; ++iData) {
    delete [] snapFile[iData];
  }
  
  delete [] snapFile;
  snapFile = NULL;
  delete [] numSnaps;
  numSnaps = NULL;
  delete [] startSnaps;
  startSnaps = NULL;
  delete [] endSnaps;
  endSnaps = NULL;
  delete [] sampleFreq;
  sampleFreq = NULL;
  delete [] snapWeight;
  snapWeight = NULL;

}


//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRomDatabaseConstruction<dim>::calcResidual(VecSet< DistSVec<double, dim> > &centers, VecSet< DistSVec<double, dim> > &centersOld) {

// Calculates the residual for the kmeans clustering.  The clustering converges when the cluster centers stop moving.

  double norm;
  double maxNorm = 0.0;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    norm = this->distanceFull( centers[iCluster], centersOld[iCluster]);
    if (norm > maxNorm) maxNorm = norm;
  }

  return maxNorm;

}


//----------------------------------------------------------------------------------

// this struct is used in the kmeans algorithm
struct sortStruct {
  int snapIndex; // snapshot #
  double dist; // distance to second closest cluster

  bool operator<(const sortStruct& a) const {
    return dist < a.dist;  
  }
};

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::kmeans() {

  this->com->fprintf(stdout, "\nUsing K-Means algorithm to cluster snapshots\n");

  // parameters that control the kmeans clustering
  int kMeansRandSeed = this->ioData->romOffline.rob.clustering.kMeansRandSeed;
  int minClusterSize = this->ioData->romOffline.rob.clustering.minClusterSize;
  double percentOverlap = this->ioData->romOffline.rob.clustering.percentOverlap;
  double kMeansTol = this->ioData->romOffline.rob.clustering.kMeansTol;

  int nTotSnaps = this->snap->numVectors();

  (this->clusterIndex) = new int[nTotSnaps];
  (this->clusterCenters) = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
  VecSet< DistSVec<double, dim> > clusterCentersOld((this->nClusters), this->domain.getNodeDistInfo());

  // pick random initial centers (use shuffle algorithm to ensure no duplicates)

  int randSeed;

  if (kMeansRandSeed == -1) {
    randSeed = time(NULL);
  } else {
    randSeed = kMeansRandSeed;
  }
 
  srand(randSeed);

  int shuffle[nTotSnaps];     

  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    shuffle[iSnap] = iSnap; 
  }

  for (int iSnap=0; iSnap<(this->nClusters); iSnap++) { // only need to shuffle first nClusters snapshots
    int randPosition = iSnap + (rand() % (nTotSnaps-iSnap));
    int temp = shuffle[iSnap];
    shuffle[iSnap] = shuffle[randPosition];
    shuffle[randPosition] = temp;
  }

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    (*(this->clusterCenters))[iCluster]=(*(this->snap))[shuffle[iCluster]];
  }

  int iterMax = this->ioData->romOffline.rob.clustering.maxIter;  // max number of kmeans iterations
  int iterSingle = this->ioData->romOffline.rob.clustering.maxIterSingle;  // number of aggressive kmeans iterations to use before switching to a more robust single update scheme

  int index1 = 0;
  int index2 = 0;

  (this->snapsInCluster) = new int[(this->nClusters)];

  // k-means algorithm
  int iter=0;
  double residual=1.0;
  while (((residual > kMeansTol) || (iter == 0)) && (iter<iterMax)) {
    
    this->com->fprintf(stdout, "Clustering iteration #%d \n", iter+1);

    double prevResidual = residual;
  
    if (iter<iterSingle) {
      this->com->fprintf(stdout, " ... updating all snapshots simultaneously\n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        this->closestCenterFull( (*(this->snap))[iSnap], &index1, &index2);
        (this->clusterIndex)[iSnap]=index1;
      }
  
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        (this->snapsInCluster)[iCluster]=0;
        clusterCentersOld[iCluster] = (*(this->clusterCenters))[iCluster];
        (*(this->clusterCenters))[iCluster] = 0.0;
      }

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        (*(this->clusterCenters))[(this->clusterIndex)[iSnap]] += (*(this->snap))[iSnap];
        ++((this->snapsInCluster)[(this->clusterIndex)[iSnap]]);  
      }

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
        double normalize = 0;
        if ((this->snapsInCluster)[iCluster] != 0) normalize = 1.0/double((this->snapsInCluster)[iCluster]);
        (*(this->clusterCenters))[iCluster] *= normalize;
      }

    residual = calcResidual(*(this->clusterCenters), clusterCentersOld);
  
    } else { // algorithm is likely stuck in a cycle -- begin updating updating one at a time

      this->com->fprintf(stdout, " ... updating one snapshot at a time\n");

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        clusterCentersOld[iCluster] = (*(this->clusterCenters))[iCluster];
      }
    
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        int oldIndex = (this->clusterIndex)[iSnap];
        int newIndex;

        this->closestCenterFull( (*(this->snap))[iSnap], &newIndex);

        if (oldIndex != newIndex) {
          (this->clusterIndex)[iSnap] = newIndex;
          (*(this->clusterCenters))[oldIndex] *= (double((this->snapsInCluster)[oldIndex])/(double((this->snapsInCluster)[oldIndex])-1.0));
          (*(this->clusterCenters))[oldIndex] -= (*(this->snap))[iSnap]*(1.0/(double((this->snapsInCluster)[oldIndex])-1.0));
          (*(this->clusterCenters))[newIndex] *= (double((this->snapsInCluster)[newIndex])/(double((this->snapsInCluster)[newIndex])+1.0));
          (*(this->clusterCenters))[newIndex] += (*(this->snap))[iSnap]*(1.0/(double((this->snapsInCluster)[newIndex])+1.0));
          --((this->snapsInCluster)[oldIndex]);
          ++((this->snapsInCluster)[newIndex]);
        } 
      }

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
      }

      residual = calcResidual(*(this->clusterCenters), clusterCentersOld);

    //  if (pow((prevResidual - residual),2)<pow(kMeansTol,2)) {
    //    this->com->fprintf(stderr, "*** Warning: Clustering algorithm is stuck. Exiting now.\n");
    //    break;
    //  }
    }
    this->com->fprintf(stdout, " ... absolute residual = %e (tolerance is set to %e)\n", residual, kMeansTol);
    ++iter;
  }


// after clustering, assimilate small clusters
  int emptyClusterCount = 0;
  if (minClusterSize>nTotSnaps) minClusterSize = nTotSnaps;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    if ((this->snapsInCluster)[iCluster]==0) {
      // if no snapshots in this cluster, skip
      ++emptyClusterCount;
    } else {
    // if smaller than tolerance, add to nearest cluster
      if ((this->snapsInCluster)[iCluster]<minClusterSize) {
        this->com->fprintf(stderr, "*** Warning: combining small cluster with nearest neighbor\n");
        int index1 = 0; // closest center (should be itself)
        int index2 = 0; // second closest center (the one we want)

        this->closestCenterFull((*(this->clusterCenters))[iCluster], &index1, &index2);
        (*(this->clusterCenters))[index2] =   (*(this->clusterCenters))[index2]*(double((this->snapsInCluster)[index2])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]))
                                            + (*(this->clusterCenters))[iCluster]*(double((this->snapsInCluster)[iCluster])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]));
 
        (*(this->clusterCenters))[iCluster] = 0.0;
 
        (this->snapsInCluster)[index2] = (this->snapsInCluster)[index2] + (this->snapsInCluster)[iCluster];
        (this->snapsInCluster)[iCluster] = 0;  

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if ((this->clusterIndex)[iSnap]==iCluster)  (this->clusterIndex)[iSnap]=index2;
        }
        
        ++emptyClusterCount;
      }
    }
  }

// remove any empty clusters, renumber existing clusters

  if (emptyClusterCount>0) {

    this->com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(((this->nClusters)-emptyClusterCount), this->domain.getNodeDistInfo()); 
    int* snapsInClusterNew = new int[((this->nClusters)-emptyClusterCount)];

    int renumberedIndices[(this->nClusters)];
    int clusterCount = 0;

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      if ((this->snapsInCluster)[iCluster]==0) {
        renumberedIndices[iCluster] = -1;
      } else {
        renumberedIndices[iCluster] = clusterCount; 
        (*clusterCentersNew)[clusterCount]=(*(this->clusterCenters))[iCluster];
        snapsInClusterNew[clusterCount]=(this->snapsInCluster)[iCluster];
        ++clusterCount;
      }
    }

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      (this->clusterIndex)[iSnap] = renumberedIndices[(this->clusterIndex)[iSnap]];
    }

    delete (this->clusterCenters);
    (this->clusterCenters) = clusterCentersNew;
    clusterCentersNew = NULL;

    delete (this->snapsInCluster);
    (this->snapsInCluster) = snapsInClusterNew; 
    snapsInClusterNew = NULL;  

    (this->nClusters) = (this->nClusters) - emptyClusterCount;
  }
  

// find snapshot closest to each center
  (this->nearestSnapsToCenters) = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      
    int sortSize = (this->snapsInCluster)[iCluster]; 
    sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

    int snapCount = 0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      if ((this->clusterIndex)[iSnap]==iCluster) { // only need to sort snapshots that are inside the current cluster
        snapDist[snapCount].snapIndex = iSnap;
        snapDist[snapCount].dist = this->distanceFull((*(this->snap))[iSnap],(*(this->clusterCenters))[iCluster]);
        ++snapCount;
      }
    }

    sort(snapDist, snapDist+sortSize);

    (*(this->nearestSnapsToCenters))[iCluster] = (*(this->snap))[snapDist[0].snapIndex];

    delete [] snapDist;
    snapDist = NULL;
  }

 
// add overlap if required
  if (percentOverlap>0) {

    this->com->fprintf(stdout, "Adding additional vectors to clusters (PercentOverlap = %2.1f%%)\n", percentOverlap);

    int index1 = 0;
    int index2 = 0;
    double dist1 = 0;
    double dist2 = 0;

    //determine which clusters are neighbors
    //approach: if a cluster is second closest to a snapshot in another cluster, then those two clusters are neighbors
    (this->clusterNeighbors) = new int*[(this->nClusters)];
    (this->clusterNeighborsCount) = new int[(this->nClusters)];

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      (this->clusterNeighborsCount)[iCluster] = 0;
      (this->clusterNeighbors)[iCluster] = new int[1];
    }

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) { // only consider snapshots that are inside of the current cluster
          this->closestCenterFull( (*(this->snap))[iSnap], &index1, &index2, &dist1, &dist2);
          //add the second closest cluster as a neighbor of current cluster (if this is the first time we've found it)
          bool unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
            if (index2 == (this->clusterNeighbors)[iCluster][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[iCluster])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[iCluster][iNeighbor];
            }
 
           clusterNeighborsNew[((this->clusterNeighborsCount)[iCluster])] = index2;
            delete (this->clusterNeighbors)[iCluster];
            (this->clusterNeighbors)[iCluster] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[iCluster];
         }
          //also add the current cluster as a neighbor of second closest cluster
          unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
            if (iCluster == (this->clusterNeighbors)[index2][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[index2])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[index2][iNeighbor];
            }
 
            clusterNeighborsNew[((this->clusterNeighborsCount)[index2])] = iCluster;
            delete (this->clusterNeighbors)[index2];
            (this->clusterNeighbors)[index2] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[index2];
         }
        }
      }
    }

    index1 = 0;
    index2 = 0;
    dist1 = 0;
    dist2 = 0;

    int origSnapsInCluster[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) origSnapsInCluster[iCluster] = (this->snapsInCluster)[iCluster];

    //share shapshots between neighboring clusters
    (this->clusterSnapshotMap) = new int*[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {  

      int numToAdd = int(ceil(double((this->snapsInCluster)[iCluster])*percentOverlap*.01/double((this->clusterNeighborsCount)[iCluster])));
      if (((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))>nTotSnaps) 
        numToAdd=int(floor(double((nTotSnaps-(this->snapsInCluster)[iCluster]))/double((this->clusterNeighborsCount)[iCluster])));
      (this->clusterSnapshotMap)[iCluster] = new int[((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))];

      // first add all snapshots that were originally in the cluster
      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) {
          (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }

      for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {

        int sortSize = origSnapsInCluster[(this->clusterNeighbors)[iCluster][iNeighbor]]; 
        sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

        int snapCount = 0;
        DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        *tmpDistVec = (*(this->clusterCenters))[iCluster] - (*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]];
        double offset = (pow(((*(this->clusterCenters))[iCluster]).norm(),2) - pow(((*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]]).norm(),2))*(-0.5);

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if ((this->clusterIndex)[iSnap]==(this->clusterNeighbors)[iCluster][iNeighbor]) { // only need to sort snapshots that are outside of the current cluster
            snapDist[snapCount].snapIndex = iSnap;
            snapDist[snapCount].dist = abs( (*tmpDistVec) * (*(this->snap))[iSnap] + offset );
           // (doesn't work well) snapDist[snapCount].dist = this->distance((*snap)[iSnap],(*(this->clusterCenters))[iCluster]);
            ++snapCount;
          }
        }
        delete tmpDistVec;
        tmpDistVec = NULL;
        sort(snapDist, snapDist+sortSize);

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          for (int jSnap=0; jSnap<numToAdd; ++jSnap) {
            if (snapDist[jSnap].snapIndex==iSnap) {
              (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
              ++mappedCount;
              ++((this->snapsInCluster)[iCluster]);
            }
          }
        }

        delete [] snapDist;
        snapDist = NULL;
      }

      this->com->fprintf(stdout, "... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);

    }
  }

  // create data structure needed for clustering non-state FOM information (explained in header file)

  stateSnapshotClustersAfterOverlap = new bool**[nSnapshotFiles];
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    stateSnapshotClustersAfterOverlap[iFile] = new bool*[stateSnapsFromFile[iFile]];
    for (int iSnap=0; iSnap<stateSnapsFromFile[iFile]; ++iSnap) {
      stateSnapshotClustersAfterOverlap[iFile][iSnap] = new bool[this->nClusters];
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        stateSnapshotClustersAfterOverlap[iFile][iSnap][iCluster] = false;
      }
    }
  }

  int** snapInfo = new int*[nTotSnaps];  // temporary - for constructing stateSnapshotClustersAfterOverlap
  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    snapInfo[iSnap] = new int[2];
  }

  int snapCount = 0;
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    for (int iSnap=0; iSnap<stateSnapsFromFile[iFile]; ++iSnap) {
      snapInfo[snapCount][0] = iFile;// file
      snapInfo[snapCount][1] = iSnap;// snapshot number from file
      ++snapCount;
    }
  }

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    for (int iSnap=0; iSnap<(this->snapsInCluster[iCluster]); ++iSnap) {
      int currentSnap = this->clusterSnapshotMap[iCluster][iSnap];
      stateSnapshotClustersAfterOverlap[snapInfo[currentSnap][0]][snapInfo[currentSnap][1]][iCluster] = true;
    }
  }

  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap)
    delete [] snapInfo[iSnap];

  delete [] snapInfo;

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::localPod(char* basisType) {

  // for limited memory SVD
  int maxVecStorage;

  //these four parameters control the size of the local bases
  double singValTolerance;
  double maxEnergyRetained;
  int maxBasisSize;
  int minBasisSize;
  
  if (strcmp(basisType, "state")==0) {
    if (strcmp(this->stateBasisName,"")==0) return;
    if (this->ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters) this->readNearestSnapsToCenters();
    maxVecStorage = this->ioData->romOffline.rob.state.dataCompression.maxVecStorage; 
    singValTolerance = this->ioData->romOffline.rob.state.dataCompression.singValTolerance;
    maxEnergyRetained = this->ioData->romOffline.rob.state.dataCompression.maxEnergyRetained;
    maxBasisSize = this->ioData->romOffline.rob.state.dataCompression.maxBasisSize;
    minBasisSize = this->ioData->romOffline.rob.state.dataCompression.minBasisSize;
  } else if (strcmp(basisType,"residual")==0) {
    if (strcmp(this->residualBasisName,"")==0) return;
    maxVecStorage = this->ioData->romOffline.rob.residual.dataCompression.maxVecStorage;
    singValTolerance = this->ioData->romOffline.rob.residual.dataCompression.singValTolerance;
    maxEnergyRetained = this->ioData->romOffline.rob.residual.dataCompression.maxEnergyRetained;
    maxBasisSize = this->ioData->romOffline.rob.residual.dataCompression.maxBasisSize;
    minBasisSize = this->ioData->romOffline.rob.residual.dataCompression.minBasisSize;
  } else if (strcmp(basisType,"jacAction")==0) {
    if (strcmp(this->jacActionBasisName,"")==0) return;
    maxVecStorage = this->ioData->romOffline.rob.jacAction.dataCompression.maxVecStorage;
    singValTolerance = this->ioData->romOffline.rob.jacAction.dataCompression.singValTolerance;
    maxEnergyRetained = this->ioData->romOffline.rob.jacAction.dataCompression.maxEnergyRetained;
    maxBasisSize = this->ioData->romOffline.rob.jacAction.dataCompression.maxBasisSize;
    minBasisSize = this->ioData->romOffline.rob.jacAction.dataCompression.minBasisSize;
  } else if (strcmp(basisType,"krylov")==0) {
    if (strcmp(this->krylovBasisName,"")==0) return;
    maxVecStorage = this->ioData->romOffline.rob.krylov.dataCompression.maxVecStorage;
    singValTolerance = this->ioData->romOffline.rob.krylov.dataCompression.singValTolerance;
    maxEnergyRetained = this->ioData->romOffline.rob.krylov.dataCompression.maxEnergyRetained;
    maxBasisSize = this->ioData->romOffline.rob.krylov.dataCompression.maxBasisSize;
    minBasisSize = this->ioData->romOffline.rob.krylov.dataCompression.minBasisSize;
  } else if (strcmp(basisType,"sensitivity")==0) {
    if (strcmp(this->sensitivityBasisName,"")==0) return;
    maxVecStorage = this->ioData->romOffline.rob.sensitivity.dataCompression.maxVecStorage;
    singValTolerance = this->ioData->romOffline.rob.sensitivity.dataCompression.singValTolerance;
    maxEnergyRetained = this->ioData->romOffline.rob.sensitivity.dataCompression.maxEnergyRetained;
    maxBasisSize = this->ioData->romOffline.rob.sensitivity.dataCompression.maxBasisSize;
    minBasisSize = this->ioData->romOffline.rob.sensitivity.dataCompression.minBasisSize;
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", basisType);
    exit (-1);
  }

  // read cluster centers
  this->readClusteredCenters();

  bool limitedMemorySVD = (maxVecStorage <= 0 ) ? false : true;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {

    if (strcmp(basisType,"sensitivity")==0)
      iCluster = -2;
 
    // compute POD
    #ifndef DO_SCALAPACK
      this->com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH SCALAPACK \n");
      exit(-1);
    #endif

    int nTotSnaps = 0;
    VecSet< DistSVec<double, dim> >* Utrue = NULL;
    Vec<double>* singVals = NULL;
    FullM* Vtrue = NULL;

    // read snapshots and preprocess
    this->readClusteredSnapshots(iCluster, true, basisType, 0, (maxVecStorage-1));  

    if (limitedMemorySVD && (this->snap->numVectors() == maxVecStorage)) {  // limited memory SVD algorithm

      this->com->fprintf(stdout, " ... beginning limited memory SVD algorithm \n");

      if (maxBasisSize <= 0) {
          this->com->fprintf(stderr, "*** Error: the limited memory SVD algorithm requires maxBasisSize to be specified \n");
          exit(-1);
      }

      VecSet< DistSVec<double, dim> >* fullSnaps = NULL; // matrix for final SVD

      int nStoredSnaps = 0;

      while (this->snap->numVectors() > maxBasisSize) {
        int nSnaps = this->snap->numVectors();

        VecSet< DistSVec<double, dim> > UTrueTmp(nSnaps, this->domain.getNodeDistInfo());
        Vec<double> singValsTmp(nSnaps);
        FullM VtrueTmp(nSnaps);

        int nKeep = maxBasisSize;
        this->com->fprintf(stdout, " ... performing SVD on matrix of size %d, storing first %d vectors to final SVD matrix \n", nSnaps, nKeep);

        ParallelRom<dim> parallelRom( this->domain, this->com);
        parallelRom.parallelSVD(*(this->snap), UTrueTmp, singValsTmp.data(), VtrueTmp, nSnaps, true);
       
        if (this->snap) delete (this->snap);
        this->snap = NULL;         
 
        VecSet< DistSVec<double, dim> >* fullSnapsNew = new VecSet< DistSVec<double, dim> >(nStoredSnaps + nKeep, this->domain.getNodeDistInfo());
        for (int iVec=0;iVec<nStoredSnaps;++iVec) {
          (*fullSnapsNew)[iVec] = (*fullSnaps)[iVec];
        }
        for (int iVec=nStoredSnaps;iVec<(nStoredSnaps+nKeep);++iVec) {
          (*fullSnapsNew)[iVec] = UTrueTmp[iVec-nStoredSnaps] * singValsTmp[iVec-nStoredSnaps];
        }
        nStoredSnaps += nKeep;
        if (fullSnaps) delete fullSnaps;
        fullSnaps = fullSnapsNew;
        fullSnapsNew = NULL;

        nTotSnaps += nSnaps;

        // read next chunk of snapshots
        this->readClusteredSnapshots(iCluster, true, basisType, nTotSnaps, (nTotSnaps + maxVecStorage - 1));  
      }

      // append extra (<maxBasisSize) snapshots to the fullSnaps matrix
      int nSnaps = this->snap->numVectors();
      if (nSnaps > 0) {
        this->com->fprintf(stdout, " ... appending trailing %d vectors to final SVD matrix \n", nSnaps);
        VecSet< DistSVec<double, dim> >* fullSnapsNew = new VecSet< DistSVec<double, dim> >(nStoredSnaps + nSnaps, this->domain.getNodeDistInfo());
        for (int iVec=0;iVec<nStoredSnaps;++iVec) {
          (*fullSnapsNew)[iVec] = (*fullSnaps)[iVec];
        }
        for (int iVec=nStoredSnaps;iVec<(nStoredSnaps+nSnaps);++iVec) {
          (*fullSnapsNew)[iVec] = (*this->snap)[iVec-nStoredSnaps];
        }
        nStoredSnaps += nSnaps;
        if (fullSnaps) delete fullSnaps;
        fullSnaps = fullSnapsNew;
        fullSnapsNew = NULL;
      }

      // SVD to compute final U
      nTotSnaps = fullSnaps->numVectors();
      this->com->fprintf(stdout, " ... performing final SVD on matrix of size %d\n", nTotSnaps);
      Utrue = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
      singVals = new Vec<double>(nTotSnaps);
      Vtrue = new FullM(nTotSnaps);
      ParallelRom<dim> parallelRom( this->domain, this->com);
      parallelRom.parallelSVD(*(fullSnaps), *Utrue, singVals->data(), *Vtrue, nTotSnaps, true);

      delete fullSnaps;
      fullSnaps = NULL;

    } else { // regular SVD
    
      nTotSnaps = this->snap->numVectors();
      Utrue = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
      singVals = new Vec<double>(nTotSnaps);
      Vtrue = new FullM(nTotSnaps);
      ParallelRom<dim> parallelRom( this->domain, this->com);
      parallelRom.parallelSVD(*(this->snap), *Utrue, singVals->data(), *Vtrue, nTotSnaps, true);
  
      delete (this->snap);
      (this->snap) = NULL;

    }
  
    //check how many vectors to keep
    double singValsTotal = 0;
    for(int i = 0; i < nTotSnaps; ++i){
      singValsTotal += (*singVals)[i];
    }
  
    int podSize = 0; 
    double singValsPartialSum = 0; 
    double target = maxEnergyRetained * singValsTotal; 
     
    for(int iSnap=0; iSnap<nTotSnaps; ++iSnap){ 
      singValsPartialSum += (*singVals)[iSnap];  
      podSize = iSnap+1; 
      if (podSize == maxBasisSize) { // setting maxBasisSize <= 0 guarantees that this is always false 
        this->com->barrier();
        this->com->fprintf(stdout, "Retaining the specified limit of %d vectors in basis %d \n", podSize, iCluster); 
        this->com->fprintf(stdout, "This basis retains %2.1f%% of the total energy\n", (singValsPartialSum/singValsTotal)*100); 
        break; 
      } else if (((singValsPartialSum/singValsTotal) >= target) && (podSize>=minBasisSize)){ 
        this->com->barrier();
        this->com->fprintf(stdout, "Reached specified energy (%2.1f%%) \n", (singValsPartialSum/singValsTotal)*100); 
        this->com->fprintf(stdout, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
        break; 
      } else if (((*singVals)[iSnap] <= singValTolerance) && (podSize>=minBasisSize)) { // singValTolerance is 1e-6 by default 
        this->com->barrier();
        this->com->fprintf(stderr, "*** Warning: Reached the singular value tolerance (%e)\n", singValTolerance); 
        this->com->fprintf(stderr, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
        break; 
      } 
    } 
   
    (this->basis) = new VecSet< DistSVec<double, dim> >(podSize, this->domain.getNodeDistInfo());
    for (int iSnap=0; iSnap<podSize; ++iSnap) { 
      (*(this->basis))[iSnap] = (*Utrue)[iSnap]; 
    }      
    delete Utrue;
    Utrue = NULL;
  
    (this->columnSumsV) = new std::vector<double>;
    this->columnSumsV->resize(nTotSnaps, 0.0);
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      for (int jSnap=0; jSnap<nTotSnaps; ++jSnap) {
        (*(this->columnSumsV))[iSnap] += (*Vtrue)[iSnap][jSnap];
      }
    }  
    delete Vtrue;
    Vtrue = NULL;
  
    (this->sVals) = new std::vector<double>;
    this->sVals->resize(nTotSnaps);
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      (*this->sVals)[iSnap] = (*singVals)[iSnap];
    } 
    delete singVals;    
    singVals = NULL;

    // output the basis and the update quantities
    this->outputClusteredBasis(iCluster, nTotSnaps, basisType); 
  
    this->com->barrier();

    if (strcmp(basisType,"sensitivity")==0)  break;
  }  

  delete [] (this->snapsInCluster);
  (this->snapsInCluster) = NULL;
  delete (this->clusterCenters);
  (this->clusterCenters) = NULL;

  if (this->nearestSnapsToCenters) {
    delete (this->nearestSnapsToCenters);
    (this->nearestSnapsToCenters) = NULL;
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::preprocessForDistanceComparisons() {

  this->com->fprintf(stdout, "\nPreprocessing for fast online cluster selection\n");

  readInitialCondition();
  
  this->readClusteredCenters();
  std::vector<std::vector<double> > clusterCenterNorms; // nClusters-by-1, or nClusters-by-dim
  clusterCenterNorms.resize(this->nClusters);
  
  if (arbitraryUniformIC) {  // preprocess for an arbitrary uniform initial condition

    // loop through all cluster centers, store squared norm of center and component-wise sums of center
    for (int iCenter=0; iCenter<(this->nClusters); ++iCenter) {
      clusterCenterNorms[iCenter].reserve(dim+1);
      DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
      double norm = ((*(this->clusterCenters))[iCenter]).norm();
      norm *= norm;
      (clusterCenterNorms[iCenter]).push_back(norm);
      double componentSums[dim];
      ((*(this->clusterCenters))[iCenter]).sum(componentSums);
      for (int iDim=0; iDim<dim; ++iDim) {
        (clusterCenterNorms[iCenter]).push_back(componentSums[iDim]);
      }
    } 

  } else { // preprocess for specified initial condition

    // loop through all cluster centers, store the squared norm of (center - IC)
    for (int iCenter=0; iCenter<(this->nClusters); ++iCenter) {
        DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
        difference = (*(this->clusterCenters))[iCenter] - *initialCondition;
        double norm = difference.norm();
        norm *= norm;  // norm squared
        (clusterCenterNorms[iCenter]).push_back(norm);
    }

  }
 
  this->outputCenterNorms(clusterCenterNorms);

  //--------- End (noUpdates || exactUpdates || approxUpdates) ---------------//


  if (this->ioData->romOffline.rob.distanceComparisons.preprocessForNoUpdates || 
      this->ioData->romOffline.rob.distanceComparisons.preprocessForExactUpdates) {

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      productOfBasisAndCenterDifferences(iCluster, "state");
      if (strcmp(this->krylovBasisName,"")!=0) productOfBasisAndCenterDifferences(iCluster, "krylov");
    }

    if (strcmp(this->sensitivityBasisName,"")!=0) productOfBasisAndCenterDifferences(-2, "sensitivity");

  }

  //--------- End (noUpdates || exactUpdates) ---------------//


  if (this->ioData->romOffline.rob.distanceComparisons.preprocessForExactUpdates) {
      this->com->fprintf(stdout, "\nPreprocessing for fast online cluster selection (For exact online ROB updates)\n");

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      productOfVectorAndCenterDifferences(iCluster, "refState");
    }
 
    productOfVectorAndCenterDifferences(-1, "initialCondition");

  }

  //--------- End (exactUpdates) ---------------//


  if (this->ioData->romOffline.rob.distanceComparisons.preprocessForApproxUpdates) {
    this->com->fprintf(stdout, "\nPreprocessing for fast online cluster selection (For approximate online ROB updates)\n");

    //

  }

  //--------- End (approxUpdates) ---------------//

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::productOfBasisAndCenterDifferences(int iCluster, char* basisType) {
  // helper function for fast distance calculation preprocessing:
  // computes w_(m,p) = 2 * basis^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that w_(m,p) = -w_(p,m)

  // note that in this case 0 <= p < m < nClusters, which differes from Amsallem et al., INJME 2012
  // (this allows for easy indexing [m][p] without any unnecessary storage)

  this->readClusteredBasis(iCluster, basisType);
  int nPodVecs = this->basis->numVectors();
  std::vector<std::vector<std::vector<double> > > w;  // collection of vectors: [m][p][iPodVec]
  w.resize(this->nClusters);

  // w[0][0] is a vector of zeros, no need to compute this
    
  for (int mCenter=1; mCenter<this->nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
      difference = (*(this->clusterCenters))[pCenter] - (*(this->clusterCenters))[mCenter]; 
      std::vector<double> product;
      product.reserve(nPodVecs);
      for (int iVec=0; iVec<nPodVecs; ++iVec) { 
        product.push_back( 2.0 * ((*(this->basis))[iVec] * difference));
      }
      w[mCenter].push_back(product); 
    }
  }
 
  this->outputClusteredInfoASCII(iCluster, basisType, NULL, NULL, &w);

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::productOfVectorAndCenterDifferences(int iCluster, char* vecType) {
  // helper function for fast distance calculation preprocessing:
  // computes w_(m,p) = 2 * vector^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that w_(m,p) = -w_(p,m)

  // note that in this case 0 <= p < m < nClusters, which differes from Amsallem et al., INJME 2012
  // (this allows for easy indexing [m][p] without any unnecessary storage)

  DistSVec<double,dim>* vec;
  vec = NULL;

  if (strcmp(vecType,"referenceState")==0) {
    this->readClusteredReferenceState(iCluster, "state");
    vec = this->Uref;
  } else if (strcmp(vecType, "initialCondition")) {
    if (initialCondition) {
      vec = initialCondition;
    } else { //need to precompute component-wise sums of (center_m - center_p)
      assert(arbitraryUniformIC);
      return;
    }
  } else {
    this->com->fprintf(stderr, "*** Error: unanticipated vector type '%s' encountered", vecType);
    exit(-1);
  }

  std::vector<std::vector<double> > result;  // collection of vectors: [m][p]
  result.resize(this->nClusters);

  // result[0][0] is a vector of zeros, no need to compute this

  for (int mCenter=1; mCenter<this->nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
      difference = (*(this->clusterCenters))[pCenter] - (*(this->clusterCenters))[mCenter]; 
      double tmp = 2.0 * ((*vec) * difference);
      result[mCenter].push_back(tmp); 
    }
  }

  this->outputClusteredInfoASCII(iCluster, vecType, NULL, &result);

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::readInitialCondition() {

  // check if this function has already been called
  if (initialCondition || arbitraryUniformIC) return;

  if (strcmp(this->ioData->input.solutions,"")!=0) {

    // read user-specified initial condition that will be used for the online ROM simulation
    initialCondition = new DistSVec<double,dim>( this->domain.getNodeDistInfo() );

    char* icFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.solutions) + 1];
    sprintf(icFile, "%s%s", this->ioData->input.prefix, this->ioData->input.solutions);

    double tmp;
    bool status = this->domain.readVectorFromFile(icFile, 0, &tmp, *initialCondition);

    if (!status) {
      this->com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", icFile);
      exit(-1);
    }

    delete [] icFile;

  } else {

    this->com->fprintf(stdout, "\n ... no initial condition specified; preprocessing for a uniform initial condition\n");
    arbitraryUniformIC = true;

  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::localRelProjError() {

  this->readClusteredCenters();
  delete (this->clusterCenters);
  (this->clusterCenters) = NULL;
 
  this->readSnapshotFile("state", true);
  int nTotSnaps = this->snap->numVectors();
 
  projErrorLog = new VecSet<Vec<double> >((this->nClusters),nTotSnaps);

  // if computing residuals
  //  if (ioDataProjError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF)
  //  ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {

    switch (ioDataProjError->relProjError) {
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_STATE):
        this->readClusteredBasis(iCluster, "state", true);
        if (ioDataProjError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF &&
            this->snapRefState!=NULL) 
          this->updateBasis(iCluster, *(this->snapRefState));
        if (ioDataProjError->krylov.include) this->appendNonStateDataToBasis(iCluster,"krylov",true);
        if (ioDataProjError->sensitivity.include) this->appendNonStateDataToBasis(iCluster,"sensitivity",true);
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_RESIDUAL):
        this->readClusteredBasis(iCluster, "residual", true);
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_JACACTION):
        this->readClusteredBasis(iCluster, "jacAction", true);
        break;
      default:
        exit (-1);
    }

    int nPodVecs = this->basis->numVectors();

    this->com->fprintf(stdout, "\nCalculating relative projection error for basis %d\n", iCluster);

    // basis^T * snapshots    
    this->com->fprintf(stdout, " ... tmp = basis^T * snapshots\n");
    Vec<double> tmpVec(nPodVecs);
    VecSet<Vec<double> >* tmpVecSet = new VecSet<Vec<double> >(nTotSnaps, nPodVecs);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      for (int iVec = 0; iVec < nPodVecs; iVec++){
        tmpVec[iVec] = (*(this->basis))[iVec] * (*(this->snap))[iSnap];
      }
      (*tmpVecSet)[iSnap] = tmpVec;
    }

    // basis * result 
    this->com->fprintf(stdout, " ... projected snaps = basis * tmp\n");
    VecSet<DistSVec<double, dim> >* projectedSnaps = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*projectedSnaps)[iSnap] = 0.0;
      for (int iVec = 0; iVec < nPodVecs; iVec++)
        (*projectedSnaps)[iSnap] += ((*tmpVecSet)[iSnap])[iVec] * (*(this->basis))[iVec];
    }
    delete tmpVecSet;
    tmpVecSet = NULL;

    // option to output postprocesed projected states (projectedSnaps + snapRefState if snapRefState is defined)
    
    if (ioDataProjError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON ) { //||
       // ioDataProjError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
  
      TsDesc<dim>* tsDesc = new TsDesc<dim>(*(this->ioData), geoSource, &(this->domain));    

      if (ioDataProjError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON) {

        DistSVec<double,dim> outVec(this->domain.getNodeDistInfo());

        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
          if (this->snapRefState) {
             outVec = *(this->snapRefState);
             outVec += (*projectedSnaps)[iSnap];
          } else {
             outVec = (*projectedSnaps)[iSnap];
          }
          tsDesc->performPostProForState(outVec);
        }
      }

      // option to output spatial residual of projected states (projectedSnaps + snapRefState if snapRefState is defined)
      //if (ioDataProjError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
       // this->spaceOp->computeResidual(*this->X, *this->A, Q, *R, this->timeState);
      //}
      if (tsDesc) delete tsDesc;
    }

    // snaps - projSnaps
    this->com->fprintf(stdout, " ... difference = originalSnaps - projectedSnaps\n");
    VecSet<DistSVec<double, dim> >* snapDifference = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*snapDifference)[iSnap] = (*(this->snap))[iSnap] - (*projectedSnaps)[iSnap];
    }

    delete projectedSnaps;
    projectedSnaps = NULL;
   
    // compute l2 norm for each 
    Vec<double> numeratorNorms(nTotSnaps);
    Vec<double> denominatorNorms(nTotSnaps);

    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      numeratorNorms[iSnap] = (*snapDifference)[iSnap].norm();
      denominatorNorms[iSnap] = (*(this->snap))[iSnap].norm();
    }

    delete snapDifference;
    snapDifference = NULL;

    this->com->fprintf(stdout, " ... relProjError = difference.norm / originalSnap.norm\n");
    Vec<double> relProjError(nTotSnaps);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
    if (denominatorNorms[iSnap] != 0) {
        relProjError[iSnap] = numeratorNorms[iSnap]/denominatorNorms[iSnap];
      } else {
        relProjError[iSnap] = 1;
      }
    }

    (*projErrorLog)[iCluster] = relProjError;

    delete (this->basis);
    (this->basis) = NULL;

  }  

  delete [] (this->snapsInCluster);
  (this->snapsInCluster) = NULL;
  delete (this->snap);
  (this->snap) = NULL;

  writeProjErrorToDisk();

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::writeProjErrorToDisk()  {

  // output projError as ASCII file in top level of ROM database
  FILE *projErrorFile;
  char *fullProjErrorName = new char[strlen(this->databasePrefix) + strlen(this->databaseName) + 1 + strlen(this->projErrorName) + 1];
  sprintf(fullProjErrorName, "%s%s/%s", this->databasePrefix, this->databaseName, this->projErrorName);
  this->com->fprintf(stdout, "\nWriting projection error to disk: %s \n", fullProjErrorName);

  int nTotSnaps = (*projErrorLog)[0].size();

  if (this->com->cpuNum() == 0) {  
    projErrorFile = fopen(fullProjErrorName, "wt");
    this->com->fprintf(projErrorFile,"Snapshot# ErrorUsingEachBasis\n");

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      this->com->fprintf(projErrorFile,"%d ", iSnap);
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      this->com->fprintf(projErrorFile,"%e ", ((*projErrorLog)[iCluster])[iSnap]);
      }
      this->com->fprintf(projErrorFile,"\n");
    }

    fclose (projErrorFile);
  }

  delete [] fullProjErrorName;
  fullProjErrorName = NULL;
  delete projErrorLog;
  projErrorLog = NULL;

  this->com->barrier();

}


