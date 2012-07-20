#include <LocalRom.h>
#include <Modal.h>
#include <TsInput.h>
#include <math.h>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

using std::stable_sort;

extern "C"      {
   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}


template<int dim> 
LocalRom<dim>::LocalRom(Communicator *_com, IoData &_ioData, Domain &_domain)  : 
com(_com), ioData(&_ioData), domain(_domain)
{ 

  // ioData->example, com->example, domain.example

  prefix = ioData->localRomDatabase.prefix;
  indexName = ioData->localRomDatabase.indexName;
  mapName = ioData->localRomDatabase.mapName;
  connName = ioData->localRomDatabase.connName;
  centersName = ioData->localRomDatabase.centersName;
  nearestName = ioData->localRomDatabase.nearestName;
  snapsName = ioData->localRomDatabase.snapsName;
  basisName = ioData->localRomDatabase.basisName;
  clusterName = ioData->localRomDatabase.clusterName;
  databaseName = ioData->localRomDatabase.databaseName;
  singValName = ioData->localRomDatabase.singValName;
  projErrorName = ioData->localRomDatabase.projErrorName;
  basisUpdateName = ioData->localRomDatabase.basisUpdateName;
  nClusters = ioData->snapshots.clustering.numClusters;  // overwritten later if there are fewer clusters
  maxBasisSize = ioData->snapshots.dataCompression.maxBasisSize;
}

//----------------------------------------------------------------------------------

template<int dim> 
LocalRom<dim>::~LocalRom() 
{
  //delete snap;
  //delete clusterCenters;
  //delete [] clusterIndex;
  
}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::createAndAnalyzeDatabase() {

// This function is called when the problem type is ROBConstruction and the number of clusters is greater than 1.
// The functions kmeans(), localPod(), and localRelProjError() are called here, depending on inputs.

  if (!(ioData->snapshots.clustering.useExistingClusters)) {
    percentOverlap = ioData->snapshots.clustering.percentOverlap;
    kMeansRandSeed = ioData->snapshots.clustering.kMeansRandSeed;
    kMeansTol = ioData->snapshots.clustering.kMeansTol;
    minClusterSize = ioData->snapshots.clustering.minClusterSize;
    readSnapshotFile(false);
    kmeans();
    outputClusterSnapshots();
  }

  if (ioData->snapshots.clustering.computePOD) {
    singValTolerance = ioData->snapshots.dataCompression.singValTolerance;
    maxPercentEnergyRetained = ioData->snapshots.dataCompression.maxPercentEnergyRetained;
    localPod();
  }

  if (ioData->snapshots.relativeProjectionError.relProjError) {
    localRelProjError();
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readSnapshotFile(bool preprocess) {

// This has been mostly borrowed from buildGlobalPod (in Modal.C).
// KMW: in the future, these portions of the codes should be merged  

  TsInput* tInput = new TsInput(*ioData);

  // Check for snapshot command file
  char *vecFile = tInput->snapFile;
  if (!vecFile) strcpy(vecFile,"snapshotFiles.in");
  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Warning: No snapshots FILES in %s\n", vecFile);
    exit (-1);
  }

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

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

  if (ioData->snapshots.incrementalSnaps)
    com->fprintf(stderr, "*** Warning: Incremental snapshots is not supported for multiple bases (yet) \n");

  if (ioData->snapshots.dataCompression.maxVecStorage > 0)
    com->fprintf(stderr, "*** Warning: Limited memory POD is not supported for multiple bases (yet) \n");

  if (ioData->snapshots.dataCompression.energyOnly == DataCompressionData::ENERGY_ONLY_TRUE)
    com->fprintf(stderr, "*** Warning: EnergyOnly is not supported for multiple bases\n");

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %d %lf", snapFile1,&nSnap,&iStart,&iEnd,&iFreq,&weight);
    if (iStart < 1) iStart = 1;
    if (iEnd < 0) iEnd = 0;
    if (iFreq < 1) iFreq = 1;
    numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);
  }

  // compute the total number of snapshots
  nTotSnaps = 0;
  for (int iData = 0; iData < nData; ++iData) {
    if (endSnaps[iData] == 0)
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = (startSnaps[iData]-1); iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
      }
    }
  }

  bool incrementalSnaps = false;
  bool subtractRefSol = false;
  if (preprocess) {
    if (ioData->snapshots.relativeProjectionError.projectIncrementalSnaps) {
      incrementalSnaps = true;
      nTotSnaps -= nData;
    } else if (ioData->snapshots.relativeProjectionError.subtractRefSol) {
      subtractRefSol = true;
      if (!(ioData->input.snapRefSolution)) {
        com->fprintf(stderr, "*** Error: Reference solution not found \n");
        exit (-1);
      }
      readReferenceSolution();
    }
  }

  int nStoredSnaps = nTotSnaps;
  
  snap = new VecSet< DistSVec<double, dim> >(nStoredSnaps, domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufOld = new DistSVec<double, dim>(domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufNew = new DistSVec<double, dim>(domain.getNodeDistInfo());

  *snapBufOld = 0.0;
  *snapBufNew = 0.0;

  double* eig = new double [nStoredSnaps];
  double eigBufOld;
  double eigBufNew;

  int numCurrentSnapshots = 0;
  bool endOfFile;

  if (subtractRefSol) {
    *snapBufOld = *snapRefSol;
    delete snapRefSol;
  }

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = (startSnaps[iData]-1); iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        if ((iSnap == (startSnaps[iData]-1)) && incrementalSnaps) {
          endOfFile = domain.readVectorFromFile(snapFile[iData], iSnap, &eigBufOld, *snapBufOld);
        } else {
          endOfFile = domain.readVectorFromFile(snapFile[iData], iSnap, &eigBufNew, *snapBufNew);
          eig[numCurrentSnapshots] = eigBufNew;
          (*snap)[numCurrentSnapshots] = *snapBufNew - *snapBufOld;  //snapBufOld = 0 if not using incremental snaps
          if (incrementalSnaps) *snapBufOld = *snapBufNew;
          if (snapWeight[iData]) (*snap)[numCurrentSnapshots] *= snapWeight[iData]; //CBM--check
          ++numCurrentSnapshots;
        }
      }
    }
  }

  delete snapBufOld;
  delete snapBufNew;

  for (int iData=0; iData < nData; ++iData) {
    delete [] snapFile[iData];
  }
  
  delete [] snapFile;
  delete [] numSnaps;
  delete [] startSnaps;
  delete [] endSnaps;
  delete [] sampleFreq;
  delete [] snapWeight;
  delete [] eig;
  //delete [] tInput;

}

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::closestCenter(DistSVec<double, dim> &vec, int &index1, int &index2, double &dist1, double &dist2) {

// Computes the closest cluster center to vec.  Returns index of closest center (index1), index of second closest center (index2)
// and associated distances (dist1 and dist2).

  double tmp;
  dist1 = 0;
  dist2 = 0;
  index1 = 0;
  index2 = 0;

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    tmp = distance( vec, (*clusterCenters)[iCluster]);
    if ((tmp<dist1) || (iCluster==0)) {
      dist2 = dist1;
      index2 = index1;
      dist1 = tmp;
      index1 = iCluster;
    } else if ((tmp<dist2) || (iCluster==1)) {
      dist2 = tmp;
      index2 = iCluster;
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
double LocalRom<dim>::distance(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {

// TODO: Implement several different distance functions

  DistSVec<double, dim> diff(domain.getNodeDistInfo());
  diff = (U1 - U2);
  double dist = diff.norm();

  return dist;

}

//----------------------------------------------------------------------------------

template<int dim>
double LocalRom<dim>::calcResidual(VecSet< DistSVec<double, dim> > &centers, VecSet< DistSVec<double, dim> > &centersOld) {

// Calculates the residual for the kmeans clustering.  The clustering converges when the cluster centers stop moving.

  double norm;
  double maxNorm = 0.0;

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    norm = distance( centers[iCluster], centersOld[iCluster]);
    if (norm > maxNorm) maxNorm = norm;
  }

  return maxNorm;

}


//----------------------------------------------------------------------------------

// this struct is used in the kmeans algorithm
struct sortStruct {
  int snapIndex; // snapshot #
  double distance; // distance to second closest cluster

  bool operator<(const sortStruct& a) const {
    return distance < a.distance;  
  }
};

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::kmeans() {

  com->fprintf(stdout, "\nUsing K-Means algorithm to cluster snapshots\n");

  clusterIndex = new int[nTotSnaps];
  clusterCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
  VecSet< DistSVec<double, dim> > clusterCentersOld(nClusters, domain.getNodeDistInfo());

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

  for (int iSnap=0; iSnap<nClusters; iSnap++) { // only need to shuffle first nClusters snapshots
    int randPosition = iSnap + (rand() % (nTotSnaps-iSnap));
    int temp = shuffle[iSnap];
    shuffle[iSnap] = shuffle[randPosition];
    shuffle[randPosition] = temp;
  }

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    (*clusterCenters)[iCluster]=(*snap)[shuffle[iCluster]];
  }

  int iterMax = ioData->snapshots.clustering.maxIter;  // max number of kmeans iterations
  int iterSingle = ioData->snapshots.clustering.maxIterSingle;  // number of aggressive kmeans iterations to use before switching to a more robust single update scheme

  int index1 = 0;
  int index2 = 0;
  double tmp1 = 0;
  double tmp2 = 0;

  snapsInCluster = new int[nClusters];

  // k-means algorithm
  int iter=0;
  double residual=1.0;
  while (((residual > kMeansTol) || (iter == 0)) && (iter<iterMax)) {
    
    com->fprintf(stdout, "Clustering iteration #%d \n", iter+1);

    double prevResidual = residual;
  
    if (iter<iterSingle) {
      com->fprintf(stdout, " ... updating all snapshots simultaneously\n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        closestCenter( (*snap)[iSnap], index1, index2, tmp1, tmp2);
        clusterIndex[iSnap]=index1;
      }
  
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        snapsInCluster[iCluster]=0;
        clusterCentersOld[iCluster] = (*clusterCenters)[iCluster];
        (*clusterCenters)[iCluster] = 0.0;
      }

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        (*clusterCenters)[clusterIndex[iSnap]] += (*snap)[iSnap];
        ++(snapsInCluster[clusterIndex[iSnap]]);  
      }

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, snapsInCluster[iCluster]);
        double normalize = 0;
        if (snapsInCluster[iCluster] != 0) normalize = 1.0/double(snapsInCluster[iCluster]);
        (*clusterCenters)[iCluster] *= normalize;
      }

    residual = calcResidual(*clusterCenters, clusterCentersOld);
  
    } else { // algorithm is likely stuck in a cycle -- begin updating updating one at a time

      com->fprintf(stdout, " ... updating one snapshot at a time\n");

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        clusterCentersOld[iCluster] = (*clusterCenters)[iCluster];
      }
    
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        int oldIndex = clusterIndex[iSnap];
        int newIndex = 1;
        int tmpIndex;
        double tmp1;
        double tmp2;

        closestCenter( (*snap)[iSnap], newIndex, tmpIndex, tmp1, tmp2);

        if (oldIndex != newIndex) {
          clusterIndex[iSnap] = newIndex;
          (*clusterCenters)[oldIndex] *= (double(snapsInCluster[oldIndex])/(double(snapsInCluster[oldIndex])-1.0));
          (*clusterCenters)[oldIndex] -= (*snap)[iSnap]*(1.0/(double(snapsInCluster[oldIndex])-1.0));
          (*clusterCenters)[newIndex] *= (double(snapsInCluster[newIndex])/(double(snapsInCluster[newIndex])+1.0));
          (*clusterCenters)[newIndex] += (*snap)[iSnap]*(1.0/(double(snapsInCluster[newIndex])+1.0));
          --(snapsInCluster[oldIndex]);
          ++(snapsInCluster[newIndex]);
        } 
      }

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, snapsInCluster[iCluster]);
      }

      residual = calcResidual(*clusterCenters, clusterCentersOld);

    //  if (pow((prevResidual - residual),2)<pow(kMeansTol,2)) {
    //    com->fprintf(stderr, "*** Warning: Clustering algorithm is stuck. Exiting now.\n");
    //    break;
    //  }
    }
    com->fprintf(stdout, " ... absolute residual = %e (tolerance is set to %e)\n", residual, kMeansTol);
    ++iter;
  }


// after clustering, assimilate small clusters
  int emptyClusterCount = 0;
  if (minClusterSize>nTotSnaps) minClusterSize = nTotSnaps;

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    if (snapsInCluster[iCluster]==0) {
      // if no snapshots in this cluster, skip
      ++emptyClusterCount;
    } else {
    // if smaller than tolerance, add to nearest cluster
      if (snapsInCluster[iCluster]<minClusterSize) {
        com->fprintf(stderr, "*** Warning: combining small cluster with nearest neighbor\n");
        int index1 = 0; // closest center (should be itself)
        int index2 = 0; // second closest center (the one we want)
        double tmp1 = 0;
        double tmp2 = 0;

        closestCenter((*clusterCenters)[iCluster], index1, index2, tmp1, tmp2);
        (*clusterCenters)[index2] = (*clusterCenters)[index2]*(double(snapsInCluster[index2])/double(snapsInCluster[index2]+snapsInCluster[iCluster]))
                                    + (*clusterCenters)[iCluster]*(double(snapsInCluster[iCluster])/double(snapsInCluster[index2]+snapsInCluster[iCluster]));
 
        (*clusterCenters)[iCluster] = 0.0;
 
        snapsInCluster[index2] = snapsInCluster[index2] + snapsInCluster[iCluster];
        snapsInCluster[iCluster] = 0;  

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if (clusterIndex[iSnap]==iCluster)  clusterIndex[iSnap]=index2;
        }
        
        ++emptyClusterCount;
      }
    }
  }

// remove any empty clusters, renumber existing clusters

  if (emptyClusterCount>0) {

    com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >((nClusters-emptyClusterCount), domain.getNodeDistInfo()); 
    int* snapsInClusterNew = new int[(nClusters-emptyClusterCount)];

    int renumberedIndices[nClusters];
    int clusterCount = 0;

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      if (snapsInCluster[iCluster]==0) {
        renumberedIndices[iCluster] = -1;
      } else {
        renumberedIndices[iCluster] = clusterCount; 
        (*clusterCentersNew)[clusterCount]=(*clusterCenters)[iCluster];
        snapsInClusterNew[clusterCount]=snapsInCluster[iCluster];
        ++clusterCount;
      }
    }

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      clusterIndex[iSnap] = renumberedIndices[clusterIndex[iSnap]];
    }

    delete clusterCenters;
    clusterCenters = clusterCentersNew;

    delete snapsInCluster;
    snapsInCluster = snapsInClusterNew;   

    nClusters = nClusters - emptyClusterCount;
  }
  

// find snapshot closest to each center
  nearestSnapsToCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      
    int sortSize = snapsInCluster[iCluster]; 
    sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

    int snapCount = 0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      if (clusterIndex[iSnap]==iCluster) { // only need to sort snapshots that are inside the current cluster
        snapDist[snapCount].snapIndex = iSnap;
        snapDist[snapCount].distance = distance((*snap)[iSnap],(*clusterCenters)[iCluster]);
        ++snapCount;
      }
    }

    sort(snapDist, snapDist+sortSize);

    (*nearestSnapsToCenters)[iCluster] = (*snap)[snapDist[0].snapIndex];

    delete [] snapDist;

  }

 
// add overlap if required
  if (percentOverlap>0) {

    com->fprintf(stdout, "Adding additional vectors to clusters (PercentOverlap = %2.1f%%)\n", percentOverlap);

    int index1 = 0;
    int index2 = 0;
    double dist1 = 0;
    double dist2 = 0;


    //determine which clusters are neighbors
    //approach: if a cluster is second closest to a snapshot in another cluster, then those two clusters are neighbors
    clusterNeighbors = new int*[nClusters];
    clusterNeighborsCount = new int[nClusters];

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      clusterNeighborsCount[iCluster] = 0;
      clusterNeighbors[iCluster] = new int[1];
    }

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if (clusterIndex[iSnap]==iCluster) { // only consider snapshots that are inside of the current cluster
          closestCenter( (*snap)[iSnap], index1, index2, dist1, dist2);
          //add the second closest cluster as a neighbor of current cluster (if this is the first time we've found it)
          bool unique = true;
          for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {
            if (index2 == clusterNeighbors[iCluster][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[(clusterNeighborsCount[iCluster])+1];
            for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = clusterNeighbors[iCluster][iNeighbor];
            }
 
           clusterNeighborsNew[(clusterNeighborsCount[iCluster])] = index2;
            delete clusterNeighbors[iCluster];
            clusterNeighbors[iCluster] = clusterNeighborsNew;
            ++clusterNeighborsCount[iCluster];
         }
          //also add the current cluster as a neighbor of second closest cluster
          unique = true;
          for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[index2]; ++iNeighbor) {
            if (iCluster == clusterNeighbors[index2][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[(clusterNeighborsCount[index2])+1];
            for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[index2]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = clusterNeighbors[index2][iNeighbor];
            }
 
           clusterNeighborsNew[(clusterNeighborsCount[index2])] = iCluster;
            delete clusterNeighbors[index2];
            clusterNeighbors[index2] = clusterNeighborsNew;
            ++clusterNeighborsCount[index2];
         }
        }
      }
    }

    index1 = 0;
    index2 = 0;
    dist1 = 0;
    dist2 = 0;

    int origSnapsInCluster[nClusters];
    for (int iCluster=0; iCluster<nClusters; ++iCluster) origSnapsInCluster[iCluster] = snapsInCluster[iCluster];

    //share shapshots between neighboring clusters
    clusterSnapshotMap = new int*[nClusters];
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {  

      int numToAdd = int(ceil(double(snapsInCluster[iCluster])*percentOverlap*.01/double(clusterNeighborsCount[iCluster])));
      if ((snapsInCluster[iCluster]+(numToAdd*clusterNeighborsCount[iCluster]))>nTotSnaps) 
        numToAdd=int(floor(double((nTotSnaps-snapsInCluster[iCluster]))/double(clusterNeighborsCount[iCluster])));
      clusterSnapshotMap[iCluster] = new int[(snapsInCluster[iCluster]+(numToAdd*clusterNeighborsCount[iCluster]))];

      // first add all snapshots that were originally in the cluster
      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if (clusterIndex[iSnap]==iCluster) {
          clusterSnapshotMap[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }

      for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {

        int sortSize = origSnapsInCluster[clusterNeighbors[iCluster][iNeighbor]]; 
        sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

        int snapCount = 0;
        DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(domain.getNodeDistInfo());
        *tmpDistVec = (*clusterCenters)[iCluster] - (*clusterCenters)[clusterNeighbors[iCluster][iNeighbor]];
        double offset = (pow(((*clusterCenters)[iCluster]).norm(),2) - pow(((*clusterCenters)[clusterNeighbors[iCluster][iNeighbor]]).norm(),2))*(-0.5);

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if (clusterIndex[iSnap]==clusterNeighbors[iCluster][iNeighbor]) { // only need to sort snapshots that are outside of the current cluster
            snapDist[snapCount].snapIndex = iSnap;
            snapDist[snapCount].distance = abs( (*tmpDistVec) * (*snap)[iSnap] + offset );
           // (doesn't work well) snapDist[snapCount].distance = distance((*snap)[iSnap],(*clusterCenters)[iCluster]);
            ++snapCount;
          }
        }
        delete tmpDistVec;
        sort(snapDist, snapDist+sortSize);

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          for (int jSnap=0; jSnap<numToAdd; ++jSnap) {
            if (snapDist[jSnap].snapIndex==iSnap) {
              clusterSnapshotMap[iCluster][mappedCount] = iSnap;
              ++mappedCount;
              ++(snapsInCluster[iCluster]);
            }
          }
        }


        delete [] snapDist;

      }

      com->fprintf(stdout, "... cluster %d has %d snaps \n", iCluster, snapsInCluster[iCluster]);

    }
  }

  /*if (percentOverlap>0) {

    com->fprintf(stdout, "Adding additional vectors to clusters (PercentOverlap = %2.1f%%)\n", percentOverlap);

    int index1 = 0;
    int index2 = 0;
    double dist1 = 0;
    double dist2 = 0;

    clusterSnapshotMap = new int*[nClusters];

      //determine which snapshots are shared
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      
      int sortSize = nTotSnaps-snapsInCluster[iCluster]; 
      sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 


      int snapCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if (clusterIndex[iSnap]!=iCluster) { // only need to sort snapshots that are outside of the current cluster
          snapDist[snapCount].snapIndex = iSnap;
          snapDist[snapCount].distance = distance((*snap)[iSnap],(*clusterCenters)[iCluster]);
          ++snapCount;
        }
      }

      sort(snapDist, snapDist+sortSize);

      int numToAdd = int(ceil(double(snapsInCluster[iCluster])*percentOverlap*.01));

      if ((snapsInCluster[iCluster]+numToAdd)>nTotSnaps) numToAdd=(nTotSnaps-snapsInCluster[iCluster]);

      clusterSnapshotMap[iCluster] = new int[(snapsInCluster[iCluster]+numToAdd)];
      
      snapCount = 0;
      bool addSnap = false;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        for (int jSnap=0; jSnap<numToAdd; ++jSnap) {
          if (snapDist[jSnap].snapIndex==iSnap) {
            addSnap=true;
            // adding a new vector to the cluster, so update cluster center
           (*clusterCenters)[iCluster] *= (double(snapsInCluster[iCluster])/(double(snapsInCluster[iCluster])+1.0));
           (*clusterCenters)[iCluster] += (*snap)[iSnap]*(1.0/(double(snapsInCluster[iCluster])+1.0));
           ++(snapsInCluster[iCluster]);
          }
        }
        if (addSnap || (clusterIndex[iSnap]==iCluster)) {
          clusterSnapshotMap[iCluster][snapCount] = iSnap;
          ++snapCount;
          addSnap=false;
        }
      }

      com->fprintf(stdout, "... cluster %d has %d snaps \n", iCluster, snapsInCluster[iCluster]);

      delete [] snapDist;

    }
  }*/
}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::outputClusterSnapshots()  { 

  // create top level of database
  int sp = strlen(prefix) + 1;
  char *fullDatabaseName = new char[sp + strlen(databaseName)];
  
  sprintf(fullDatabaseName, "%s%s", prefix, databaseName);

  if (com->cpuNum() == 0) {
    int status;
    status = mkdir(fullDatabaseName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  // output cluster index as ASCII file (original clusters before sharing)
  FILE *clustIndexFile;
  char *fullClustIndexName = new char[strlen(fullDatabaseName) + 1 + strlen(indexName) + 1];
  sprintf(fullClustIndexName, "%s/%s", fullDatabaseName, indexName);
  com->fprintf(stdout, "\nWriting cluster index to disk\n");

  if (com->cpuNum() == 0) {  
    clustIndexFile = fopen(fullClustIndexName, "wt");
    com->fprintf(clustIndexFile,"Snapshot# OriginalCluster\n");

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      com->fprintf(clustIndexFile,"%d %d\n", iSnap, clusterIndex[iSnap]);
    }

    fclose (clustIndexFile);
  }

  delete [] clusterIndex;

  // output cluster snapshot map as ASCII file (clusters after sharing)
  FILE *clustMapFile;
  char *fullClustMapName = new char[strlen(fullDatabaseName) + 1 + strlen(mapName) + 1];
  sprintf(fullClustMapName, "%s/%s", fullDatabaseName, mapName);
  com->fprintf(stdout, "\nWriting cluster-snapshot map to disk\n");

  if (com->cpuNum() == 0) {  
    clustMapFile = fopen(fullClustMapName, "wt");
    com->fprintf(clustMapFile,"Snapshot# Cluster\n");

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      for (int iSnap=0; iSnap<snapsInCluster[iCluster]; ++iSnap) {
        com->fprintf(clustMapFile,"%d %d\n", clusterSnapshotMap[iCluster][iSnap], iCluster);
      }
    }
    fclose (clustMapFile);
  }

  //TODO patch everything for the no-overlap scenario
  // output cluster connectivity as ASCII file
  FILE *clustConnFile;
  char *fullClustConnName = new char[strlen(fullDatabaseName) + 1 + strlen(connName) + 1];
  sprintf(fullClustConnName, "%s/%s", fullDatabaseName, connName);
  com->fprintf(stdout, "\nWriting cluster connectivity to disk\n");

  if (com->cpuNum() == 0) {  
    clustConnFile = fopen(fullClustConnName, "wt");
    com->fprintf(clustConnFile,"Cluster#  ...Neighbors...\n");
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->fprintf(clustConnFile,"%d", iCluster);
      for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {
        com->fprintf(clustConnFile," %d", clusterNeighbors[iCluster][iNeighbor]);
      }
      com->fprintf(clustConnFile,"\n");
    }
    fclose (clustConnFile);
  }

  com->barrier();

  // output cluster centers
  char *fullClustCentersName = new char[strlen(fullDatabaseName) + 1 + strlen(centersName) + 1];
  sprintf(fullClustCentersName, "%s/%s", fullDatabaseName, centersName);
  com->fprintf(stdout, "\nWriting cluster centers to disk\n");

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    com->barrier();
    domain.writeVectorToFile(fullClustCentersName, iCluster, double(snapsInCluster[iCluster]), (*clusterCenters)[iCluster]);
  }


  // output nearest snap to each cluster
  char *fullNearestSnapsName = new char[strlen(fullDatabaseName) + 1 + strlen(nearestName) + 1];
  sprintf(fullNearestSnapsName, "%s/%s", fullDatabaseName, nearestName);
  com->fprintf(stdout, "\nWriting nearest snapshot to each center to disk\n");

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    com->barrier();
    domain.writeVectorToFile(fullNearestSnapsName, iCluster, double(snapsInCluster[iCluster]), (*nearestSnapsToCenters)[iCluster]);
  }


  // iterate through clusters
  for (int iCluster=0; iCluster<nClusters; iCluster++) {

    int addedDigits = 1;
    if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
 
    char *fullClusterName = new char[strlen(fullDatabaseName) + 1 + strlen(clusterName) + addedDigits + 1];
    sprintf(fullClusterName, "%s/%s%d", fullDatabaseName, clusterName, iCluster);

    if (com->cpuNum() == 0) {
      int status;
      status = mkdir(fullClusterName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    com->barrier();
    sleep(5);

    char *fullSnapshotsName = new char[strlen(fullClusterName) + 1 + strlen(snapsName) + 1];
    sprintf(fullSnapshotsName, "%s/%s", fullClusterName, snapsName);

    com->fprintf(stdout, "\nWriting %d snapshots to cluster %d\n", snapsInCluster[iCluster], iCluster);

    /*  // Again, this would be a pain in the ass to implement
    bool incrementalSnaps = false;
    DistSVec<double, dim> > *snapBuf = new DistSVec<double, dim> >(domain.getNodeDistInfo());
    if ((ioData->snapshots.IncrementalSnaps) && !(ioData->snapshots.subtractCenters || ioData->snapshots.subtractNearestSnapsToCenters ||ioData->snapshots.subtractReference)) {
      com->fprintf(stderr, "*** Warning: Writing INCREMENTAL snapshots to disk \n");
      incrementalSnaps = true;  
      snapBuf = 0.0; 
    }*/

    int numWritten = 0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      for (int jSnap=0; jSnap<snapsInCluster[iCluster]; ++jSnap) {
        if (iSnap == clusterSnapshotMap[iCluster][jSnap]) {
          domain.writeVectorToFile(fullSnapshotsName, numWritten, double(numWritten), (*snap)[iSnap] );
          ++numWritten;
        }
      }
    } 

    //delete [] snapBuf;
    delete [] fullSnapshotsName;
    delete [] fullClusterName;

  }

  com->fprintf(stdout, "\nFreeing memory for parallel SVD; read in snapshots as needed\n");

  delete [] fullDatabaseName;
  delete [] fullClustIndexName;
  delete [] fullClustMapName;
  delete [] fullClustCentersName; 
  delete [] fullNearestSnapsName;

  delete snap;
  delete clusterCenters;
  delete nearestSnapsToCenters;

  for (int iCluster=0; iCluster<nClusters; ++iCluster) delete [] clusterSnapshotMap[iCluster];
  delete [] clusterSnapshotMap;
  delete [] snapsInCluster;

  for (int iCluster=0; iCluster<nClusters; ++iCluster) delete [] clusterNeighbors[iCluster];
  delete [] clusterNeighbors;
  delete [] clusterNeighborsCount;
  
}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::localPod() {

  // read cluster centers
  readClusterCenters();

  if (ioData->snapshots.subtractNearestSnapsToCenters) 
    readNearestSnapsToCenters();

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {

  // read snapshots and preprocess
  readClusterSnapshots(iCluster, true);  

  // compute POD
  #ifndef DO_SCALAPACK
    this->com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH SCALAPACK \n");
    exit(-1);
  #endif

  VecSet< DistSVec<double, dim> >* Utrue = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
  Vec<double> singVals(nTotSnaps);
  FullM Vtrue(nTotSnaps);
  ParallelRom<dim> parallelRom( domain, com);
  parallelRom.parallelSVD(*snap, *Utrue, singVals.data(), Vtrue, nTotSnaps, true);

  delete snap;

  //check how many vectors to keep
  double singValsTotal = 0;
  for(int i = 0; i < nTotSnaps; ++i){
    singValsTotal += singVals[i];
  }

  int podSize = 0; 
  double singValsPartialSum = 0; 
  double target = maxPercentEnergyRetained * singValsTotal * .01; 
     
  for(int iSnap=0; iSnap<nTotSnaps; ++iSnap){ 
    singValsPartialSum += singVals[iSnap];  
    if ((iSnap+1) == maxBasisSize) { // setting maxBasisSize = 0 guarantees that this is always false 
      podSize = iSnap+1; 
      com->barrier();
      com->fprintf(stdout, "Retaining the specified limit of %d vectors in basis %d \n", podSize, iCluster); 
      com->fprintf(stdout, "This basis retains %2.1f%% of the total energy\n", (singValsPartialSum/singValsTotal)*100); 
      break; 
    } else if (singValsPartialSum >= target){ 
      podSize = iSnap+1; 
      com->barrier();
      com->fprintf(stdout, "Reached specified energy \n"); 
      com->fprintf(stdout, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
      break; 
    } else if (singVals[iSnap] <= singValTolerance) { // 1e-6 by default 
      podSize = iSnap+1; 
      com->barrier();
      com->fprintf(stderr, "*** Warning: Reached the singular value tolerance (%e)\n", singValTolerance); 
      com->fprintf(stderr, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
      break; 
    } 
    ++podSize;
  } 
 
  basis = new VecSet< DistSVec<double, dim> >(podSize, domain.getNodeDistInfo());
  for (int iSnap=0; iSnap<podSize; ++iSnap) { 
    (*basis)[iSnap] = (*Utrue)[iSnap]; 
  }      
  delete Utrue;

  columnSumsV = new double[nTotSnaps];
  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    columnSumsV[iSnap] = 0.0;
    for (int jSnap=0; jSnap<nTotSnaps; ++jSnap) {
      columnSumsV[iSnap] += Vtrue[iSnap][jSnap]; //TODO check to make sure Vtrue is not actually V^T!
    }
  }  
   
  sVals = new double[nTotSnaps]; 
  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    sVals[iSnap] = singVals[iSnap];
  } 
  
  // output the basis and the update quantities
  outputClusterBasis(iCluster); 

  com->barrier();
  // free memory
  delete basis;
 
  }  

  delete [] snapsInCluster;
  delete clusterCenters;

  if (ioData->snapshots.subtractNearestSnapsToCenters)
    delete nearestSnapsToCenters;

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readClusterSnapshots(int iCluster, bool preprocess) {

  // reconstruct cluster name
  int addedDigits = 1;
  if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
  char *fullClusterName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(snapsName) + 1];
  sprintf(fullClusterName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, snapsName);

  // read snapshots from cluster
  nTotSnaps = snapsInCluster[iCluster];
  snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());

  if ((ioData->snapshots.subtractRefSol) && !(ioData->snapshots.subtractCenters || ioData->snapshots.subtractNearestSnapsToCenters)){
    if (!(ioData->input.snapRefSolution)) {
      com->fprintf(stderr, "*** Error: Reference solution not found \n");
      exit (-1);
    }
    readReferenceSolution();
  }

  com->fprintf(stdout, "\nReading %d snapshots from cluster %d \n", nTotSnaps, iCluster);

  double tmp;
  bool endOfFile; 

  // read in Snapshot Vectors
  for (int iSnap = 0; iSnap<nTotSnaps; ++iSnap) {
    endOfFile = domain.readVectorFromFile(fullClusterName, iSnap, &tmp, (*snap)[iSnap]);
  }

  if (preprocess) {
 
    if (ioData->snapshots.subtractNearestSnapsToCenters) {
      com->fprintf(stdout, " ... subtracting nearest snapshot to center from snapshots \n");
      if (ioData->snapshots.subtractCenters)
        com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring subtractCenters command \n");
      if (ioData->snapshots.subtractRefSol) 
        com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
      if (ioData->snapshots.incrementalSnaps) 
        com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
      --nTotSnaps; //
      VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
      DistSVec<double, dim>* tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());
      outputClusterReferenceSnapshot(iCluster, (*nearestSnapsToCenters)[iCluster]);
      com->barrier();
      int snapCount = 0;
      for (int iSnap=0; iSnap<=nTotSnaps; ++iSnap) {
        *tmpSnap = (*snap)[iSnap] - (*nearestSnapsToCenters)[iCluster];
        if ((*tmpSnap).norm() > 1e-6) { // one will be exactly zero; don't include this one
          (*snapNew)[snapCount] = *tmpSnap;
          ++snapCount;
        }
      }
      delete tmpSnap;
      delete snap;
      snap = snapNew;
    } else if (ioData->snapshots.subtractCenters) {
      com->fprintf(stdout, " ... subtracting cluster center from snapshots \n");
      if (ioData->snapshots.subtractRefSol) 
        com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
      if (ioData->snapshots.incrementalSnaps) 
        com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
      outputClusterReferenceSnapshot(iCluster,(*clusterCenters)[iCluster]);
      com->barrier();
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        (*snap)[iSnap] -= (*clusterCenters)[iCluster];
      }
    } else if (ioData->snapshots.subtractRefSol) {
      com->fprintf(stdout, " ... subtracting reference solution from snapshots \n");
      if (ioData->snapshots.incrementalSnaps) 
        com->fprintf(stderr, "*** Warning: Incompatible commands -- not performing incremental snapshots \n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        (*snap)[iSnap] -= *snapRefSol;
      }
      outputClusterReferenceSnapshot(iCluster,(*snapRefSol));
      com->barrier();
      delete snapRefSol;
    }

    if (ioData->snapshots.normalizeSnaps) {
      com->fprintf(stdout, " ... normalizing snapshots \n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        double tmpNorm = ((*snap)[iSnap]).norm();
        if (tmpNorm != 0) {
          double normalize = 1.0 / tmpNorm;
          (*snap)[iSnap] *= normalize;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::outputClusterReferenceSnapshot(int iCluster, DistSVec<double, dim>& ref) {
  // Automatically stores the reference snapshot for cluster iCluster.
  // This functionality is in place to reduce the user's workload.

  // form file name
  int addedDigits = 1;
  if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
  char *fullSnapRefName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisName) + strlen(".snapReference") + 1];
  sprintf(fullSnapRefName, "%s%s/%s%d/%s%s", prefix, databaseName, clusterName, iCluster, basisName, ".snapReference");

  // output vector
  com->fprintf(stdout, " ... storing reference snapshot to disk for later use\n", iCluster);
  domain.writeVectorToFile(fullSnapRefName, 0, double(iCluster), ref);

  delete [] fullSnapRefName;

}

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readClusterReferenceSnapshot(int iCluster) {
  // This function reads in the automatically stored reference snapshot for a cluster.
  // By storing these reference snapshots and reading them automatically it reduces the user's workload.

  // form file name
  int addedDigits = 1;
  if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
  char *fullSnapRefName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisName) + strlen(".snapReference") + 1];
  sprintf(fullSnapRefName, "%s%s/%s%d/%s%s", prefix, databaseName, clusterName, iCluster, basisName, ".snapReference");

  double tmp;
  bool endOfFile;
 
  com->fprintf(stdout, "Reading the reference snapshot for this cluster %s\n", fullSnapRefName);
  Uinit = new DistSVec<double, dim>(domain.getNodeDistInfo());
  endOfFile = domain.readVectorFromFile(fullSnapRefName, 0, &tmp, *Uinit);

  delete [] fullSnapRefName;

}

//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readReferenceSolution() {

  double tmp;
  bool endOfFile;
  char* fullRefName;
  const char* refSolName = ioData->input.snapRefSolution;
  fullRefName = new char[strlen(ioData->input.prefix) + strlen(refSolName) + 1];
  sprintf(fullRefName, "%s%s", ioData->input.prefix, refSolName);
  com->fprintf(stdout, "\nReading reference solution for snapshots from %s\n", fullRefName);
  snapRefSol = new DistSVec<double, dim>(domain.getNodeDistInfo());
  endOfFile = domain.readVectorFromFile(fullRefName, 0, &tmp, *snapRefSol);
  delete [] fullRefName;

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readClusterCenters() {

  // reconstruct cluster centers file name
  char *fullClusterCentersName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(centersName) + 1];
  sprintf(fullClusterCentersName, "%s%s/%s", prefix, databaseName, centersName);

  // read centers
  clusterCenters = new VecSet< DistSVec<double, dim> >(10, domain.getNodeDistInfo());
  snapsInCluster = new int[10];

  com->fprintf(stdout, "\nReading cluster centers\n");
  bool endOfFile;
  nClusters = 0;
  double tmp;
  DistSVec<double, dim>* tmpVec = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {

    endOfFile = domain.readVectorFromFile(fullClusterCentersName, nClusters, &tmp, *tmpVec);
    if (!endOfFile) break;
    
    ++nClusters;

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
    int* snapsInClusterNew = new int[nClusters];

    for (int iCluster=0; iCluster<(nClusters-1); ++iCluster) {
      (*clusterCentersNew)[iCluster] = (*clusterCenters)[iCluster];
      snapsInClusterNew[iCluster] = snapsInCluster[iCluster];
    }

    (*clusterCentersNew)[nClusters-1] = *tmpVec;
    snapsInClusterNew[nClusters-1] = int(tmp);

    delete clusterCenters;
    delete [] snapsInCluster;

    clusterCenters = clusterCentersNew;
    snapsInCluster = snapsInClusterNew;
 
  }

  com->fprintf(stdout, "\n ... found %d clusters in this database\n", nClusters);
  com->barrier();
  delete tmpVec;
  delete [] fullClusterCentersName;

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readNearestSnapsToCenters() {

  // NOTE: this function must only be run after readClusterCenters (this is because the
  // cluster centers file defines snapsInCluster and nClusters)

  // reconstruct file name
  char *fullNearestSnapsName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(nearestName) + 1];
  sprintf(fullNearestSnapsName, "%s%s/%s", prefix, databaseName, nearestName);

  // read centers
  nearestSnapsToCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  bool endOfFile;

  com->fprintf(stdout, "\nReading closest snapshot to each center\n");

  // read in Snapshot Vectors
  double tmp;
  for (int iCluster = 0; iCluster<nClusters; ++iCluster) {
    endOfFile = domain.readVectorFromFile(fullNearestSnapsName, iCluster, &tmp, (*nearestSnapsToCenters)[iCluster]);
    if (!endOfFile) {
       com->fprintf(stderr,"*** Error: readNearestSnapsToCenters attempted to read %d vecs from a file with only %d.\n", nClusters, iCluster-1);
       exit(-1);
    }
  }

  delete [] fullNearestSnapsName;

}
//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::readClusterBasis(int iCluster) {

  int basisSize;
  ( maxBasisSize >= ioData->rom.dimension ) ? basisSize = maxBasisSize : basisSize = ioData->rom.dimension;

  int addedDigits = 1;
  if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));

  // read singular values
  char *fullSingValName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(singValName) + 1];
  sprintf(fullSingValName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, singValName);
  com->fprintf(stdout, "\nReading singular values for basis %d\n", iCluster );

  FILE *singValFile = fopen(fullSingValName, "r");

  if (!singValFile)  {
    com->fprintf(stderr, "*** Error: unable to open file %s\n", fullSingValName);
    exit (-1);
  }

  //char dummyString;
  //int _n = fscanf(singValFile,"%s\n", &dummyString);
  int _n;
  int vecNumber;
  double percentEnergy;
  sVals = new double[basisSize];//[snapsInCluster[iCluster]];

  double energyTol = ioData->rom.energy;

 // for (int iVec=0; iVec<snapsInCluster[iCluster]; ++iVec) {
  for (int iVec=0; iVec<basisSize; ++iVec) {
   fscanf(singValFile,"%d %le %le", &vecNumber, &(sVals[iVec]), &percentEnergy);
   if ((percentEnergy>=energyTol)&&((iVec+1)>=(ioData->rom.minDimension))) {
     basisSize = (iVec+1);
     com->fprintf(stdout, " ... using %d vectors from basis %d (energy tolerance = %e)\n", basisSize, iCluster, energyTol);
     double* sValsNew = new double[basisSize];
     for (int iSVal = 0; iSVal<basisSize; ++iSVal)
       sValsNew[iSVal] = sVals[iSVal];
     delete sVals;
     sVals = sValsNew;
     break;
   }
  }

  fclose(singValFile);

  delete [] fullSingValName;


  // read in basis vectors
  char* fullBasisName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisName) + 1];
  sprintf(fullBasisName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, basisName);

  basis = new VecSet< DistSVec<double, dim> >(basisSize, domain.getNodeDistInfo());
  bool endOfFile;

  com->fprintf(stdout, "\nReading basis %d\n", iCluster);
  int numBasisVecs = 0;
  double tmpSVal;
  for (int iVec = 0; iVec<basisSize; ++iVec) {
    endOfFile = domain.readVectorFromFile(fullBasisName, iVec, &tmpSVal, (*basis)[iVec]);
    ++numBasisVecs; 
    if (!endOfFile) break;
  }

  delete [] fullBasisName;

  if (numBasisVecs != basisSize) {
    VecSet< DistSVec<double, dim> >* basisNew =  new VecSet< DistSVec<double, dim> >(numBasisVecs, domain.getNodeDistInfo());

    for (int iVec=0; iVec<numBasisVecs; ++iVec) {
      (*basisNew)[iVec]=(*basis)[iVec];
    }

    delete basis;
    basis = basisNew;
  }    

  // also read in the basis update quantities
  char *fullBasisUpdateName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisUpdateName) + 1];
  sprintf(fullBasisUpdateName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, basisUpdateName);
  com->fprintf(stdout, "\nReading update info for basis %d \n", iCluster);

  FILE  *basisUpdateFile = fopen(fullBasisUpdateName, "r");

  if (!basisUpdateFile)  {
    com->fprintf(stderr, "*** Error: unable to open file %s\n", fullBasisUpdateName);
    exit (-1);
  }

  columnSumsV = new double[snapsInCluster[iCluster]];
  for (int iVec = 0; iVec < snapsInCluster[iCluster]; ++iVec){
    _n = fscanf(basisUpdateFile, "%le", &(columnSumsV[iVec]));
    //com->fprintf(stderr, "*** colSumsV %e\n", columnSumsV[iVec]);
  }

  fclose(basisUpdateFile);

  delete [] fullBasisUpdateName;
  com->barrier();

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::outputClusterBasis(int iCluster) {

  int podSize = basis->numVectors();

  // form basis file name
  int addedDigits = 1;
  if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
  char *fullBasisName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisName) + 1];
  sprintf(fullBasisName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, basisName);

  // output vectors
  com->fprintf(stdout, "\nWriting basis %d to disk\n", iCluster);
  for (int iVec=0; iVec<podSize; ++iVec) {
      domain.writeVectorToFile(fullBasisName, iVec, double(iVec), (*basis)[iVec] );
  }

  delete [] fullBasisName;

  // write singular values to ASCII file
  FILE *singValFile;
  char *fullSingValName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(singValName) + 1];
  sprintf(fullSingValName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, singValName);
  com->fprintf(stdout, "\nWriting singular values to disk\n");

  double sValSum = 0;
  for (int iSnap = 0; iSnap<nTotSnaps; ++iSnap) sValSum += sVals[iSnap];
  double invSValSum = 1/sValSum;

  double sValPartialSum = 0;

  if (com->cpuNum() == 0) {
    singValFile = fopen(fullSingValName, "wt");
    //com->fprintf(singValFile,"Vector# SVal SValSum\n");

    for (int iVec=0; iVec<nTotSnaps; ++iVec) {
      sValPartialSum += (sVals[iVec]*invSValSum);
      com->fprintf(singValFile,"%d %e %e\n", iVec, sVals[iVec], sValPartialSum);
    }

    fclose (singValFile);
  }
 
  delete [] fullSingValName;
  delete [] sVals;


 // write thinSVD update info to ASCII file
  FILE *basisUpdateFile;
  char *fullBasisUpdateName = new char[strlen(prefix) + strlen(databaseName) + 1
                                   + strlen(clusterName) + addedDigits + 1 + strlen(basisUpdateName) + 1];
  sprintf(fullBasisUpdateName, "%s%s/%s%d/%s", prefix, databaseName, clusterName, iCluster, basisUpdateName);
  com->fprintf(stdout, "\nWriting thin SVD update info to disk\n");

  if (com->cpuNum() == 0) {
    basisUpdateFile = fopen(fullBasisUpdateName, "wt");
 //   com->fprintf(basisUpdateFile,"%e\n", normQuantity);

    for (int iVec=0; iVec<snapsInCluster[iCluster]; ++iVec) {
      com->fprintf(basisUpdateFile,"%le\n", columnSumsV[iVec]);
    }

    fclose (basisUpdateFile);
  }

  delete [] fullBasisUpdateName;
  delete [] columnSumsV;
  com->barrier();
}

//----------------------------------------------------------------------------------


template<int dim>
void LocalRom<dim>::localRelProjError() {

  readClusterCenters();
  delete clusterCenters;
  delete [] snapsInCluster;

  readSnapshotFile(true);  
  projErrorLog = new VecSet<Vec<double> >(nClusters,nTotSnaps);

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {

    readClusterBasis(iCluster);
    int nPodVecs = basis->numVectors();

    com->fprintf(stdout, "\nCalculating relative projection error for basis %d\n", iCluster);

    // basis^T * snapshots    
    com->fprintf(stdout, " ... tmp = basis^T * snapshots\n");
    Vec<double> tmpVec(nPodVecs);
    VecSet<Vec<double> >* tmpVecSet = new VecSet<Vec<double> >(nTotSnaps, nPodVecs);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      com->fprintf(stdout, " ***debugging: Currently working on snapshot %d\n", iSnap);
      for (int iVec = 0; iVec < nPodVecs; iVec++){
        tmpVec[iVec] = (*basis)[iVec] * (*snap)[iSnap];
      }
      (*tmpVecSet)[iSnap] = tmpVec;
    }

    // basis * result 
    com->fprintf(stdout, " ... projected snaps = basis * tmp\n");
    VecSet<DistSVec<double, dim> >* projectedSnaps = new VecSet<DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*projectedSnaps)[iSnap] = 0.0;
      for (int iVec = 0; iVec < nPodVecs; iVec++)
        (*projectedSnaps)[iSnap] += ((*tmpVecSet)[iSnap])[iVec] * (*basis)[iVec];
    }
    delete tmpVecSet;

    // snaps - projSnaps
    com->fprintf(stdout, " ... difference = originalSnaps - projectedSnaps\n");
    VecSet<DistSVec<double, dim> >* snapDifference = new VecSet<DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*snapDifference)[iSnap] = (*snap)[iSnap] - (*projectedSnaps)[iSnap];
    }

    delete projectedSnaps;

    // compute l2 norm for each 
    Vec<double> numeratorNorms(nTotSnaps);
    Vec<double> denominatorNorms(nTotSnaps);

    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      numeratorNorms[iSnap] = (*snapDifference)[iSnap].norm();
      denominatorNorms[iSnap] = (*snap)[iSnap].norm();
    }

    delete snapDifference;

    com->fprintf(stdout, " ... relProjError = difference.norm / originalSnap.norm\n");
    Vec<double> relProjError(nTotSnaps);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      if (denominatorNorms[iSnap] != 0) {
        relProjError[iSnap] = numeratorNorms[iSnap]/denominatorNorms[iSnap];
      } else {
        relProjError[iSnap] = 1;
      }
    }

    (*projErrorLog)[iCluster] = relProjError;

  }  

  delete snap;
  writeProjErrorToDisk();

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::writeProjErrorToDisk()  {

  // output projError as ASCII file in top level of ROM database
  FILE *projErrorFile;
  char *fullProjErrorName = new char[strlen(prefix) + strlen(databaseName) + 1 + strlen(projErrorName) + 1];
  sprintf(fullProjErrorName, "%s%s/%s", prefix, databaseName, projErrorName);
  com->fprintf(stdout, "\nWriting projection error to disk: %s \n", fullProjErrorName);

  if (com->cpuNum() == 0) {  
    projErrorFile = fopen(fullProjErrorName, "wt");
    com->fprintf(projErrorFile,"Snapshot# ErrorUsingEachBasis\n");

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      com->fprintf(projErrorFile,"%d ", iSnap);
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->fprintf(projErrorFile,"%e ", ((*projErrorLog)[iCluster])[iSnap]);
      }
      com->fprintf(projErrorFile,"\n");
    }

    fclose (projErrorFile);
  }

  delete [] fullProjErrorName;
  delete projErrorLog;

  com->barrier();

}


//----------------------------------------------------------------------------------

template<int dim>
void LocalRom<dim>::updateBasis(int iCluster, DistSVec<double, dim> &U) {

/* 
  When updateBasis is called the following quantities are available:
    - nClusters: number of clusters
    - snapsInCluster:  a vector containing the number of snapshots in each cluster, size nClusters
    - columnSumsV: a vector containing the sums of the columns of V, size ny + buffer
    - basis: the local ROB, ny + buffer
    - sVals: the singular values, size ny + buffer
  (As well as all of the ioData values)

*/
  readClusterReferenceSnapshot(iCluster); // reads Uinit

  int robSize = basis->numVectors();
  int kSize = robSize+1;

  DistSVec<double, dim> a(domain.getNodeDistInfo());
  a = *Uinit - U;
  
  if (a.norm() >= 1e-6) {  // only update if Uinit is different than U (this handles the case of time=0) 

    double m[robSize];
    for (int iVec=0; iVec<robSize; ++iVec) {
      m[iVec] = (*basis)[iVec] * a;
    }

    DistSVec<double, dim> p(domain.getNodeDistInfo());
    p = a;

    for (int iVec=0; iVec<robSize; ++iVec) {
      p -= (*basis)[iVec] * m[iVec];
    }

    double Ra = p.norm();
    double RaInv = 1/Ra; 
    p *= RaInv;

    double *K = new double[(kSize)*(kSize)];

    for (int iCol = 0; iCol < (kSize); ++iCol){
      for (int iRow = 0; iRow < (kSize); ++iRow) {
        if ((iCol == iRow) && (iCol < kSize-1)) {
          K[iCol*(kSize) + iRow] = sVals[iCol];
        } else {
          K[iCol*(kSize) + iRow] = 0.0;
        }
      }
    }

    double q = 0;
    for (int iVec=robSize; iVec<snapsInCluster[iCluster]; ++iVec) {
      q += pow(columnSumsV[iVec], 2);
    }
    q = pow(q, 0.5);

    for (int iRow = 0; iRow < (kSize-1); ++iRow) {
      for (int iCol = 0; iCol < (kSize-1); ++iCol){
        K[iCol*(kSize) + iRow] += m[iRow] * columnSumsV[iCol];
      }
      K[(kSize-1)*(kSize) + iRow] = m[iRow] * q;
    }

    for (int iCol = 0; iCol < kSize-1; ++iCol){
      K[iCol*(kSize) + (kSize-1)] += Ra * columnSumsV[iCol];
    }
    K[(kSize-1)*(kSize) + kSize-1] += Ra * q;

    double *sigma = new double[kSize];
    double *error = new double[kSize];
    double *work = new double[kSize];
    int info;
    double *zVec = new double[kSize*kSize]; // right singular vectors
    double *yVec = new double[kSize*kSize]; // left singular vectors

    com->fprintf(stdout, " ... computing rank one update to basis using current state\n");
    F77NAME(dsvdc)(K, kSize, kSize, kSize, sigma, error, yVec, kSize, zVec, kSize, work, 11, info);

    VecSet< DistSVec<double, dim> > basisOld(robSize, domain.getNodeDistInfo());

    for (int iVec=0; iVec<robSize; ++iVec)
        basisOld[iVec] = (*basis)[iVec];
 
    for (int iVec=0; iVec<robSize; ++iVec) {
      (*basis)[iVec] = p * yVec[(iVec*kSize) + kSize-1];
      for (int jVec=0; jVec<robSize; ++jVec) {
        (*basis)[iVec] += yVec[(iVec*kSize) + jVec] * basisOld[jVec];
      }
    }
  }
  delete Uinit;
}


