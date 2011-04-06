/* Multigrid.C

 */

#include <Multigrid.h>

#include <map>

#include <stdio.h>

#include <KDTree.h>

#include <fstream>

using std::pair;
using std::ifstream;

void Multigrid::constructMapping(SVec<double,3>& Xfine,
				 SVec<double,3>& Xcoarse,
				 int** decFine, int nSubDFine,int* subDsizeFine,
				 int** decCoarse,int nSubDCoarse,int* subDsizeCoarse,
				 const char* packageFile,
				 const char* collectionFile,
				 double radius0,int threshold,
				 double radiusf) {

 
  KDTree* myKDTrees = new KDTree[nSubDCoarse];

  FILE** packageFiles = new FILE*[nSubDCoarse];
    
  for (int i = 0; i < nSubDCoarse; ++i) {
    SVec<double,3> Xl(subDsizeCoarse[i]);
    for (int j = 0; j < subDsizeCoarse[i]; ++j) {
      
      memcpy(Xl[j], Xcoarse[ decCoarse[i][j] ], sizeof(double)*3);
    }
    myKDTrees[i].construct(Xl);

    char pfile[256];
    sprintf(pfile,"%s%i",packageFile,i);
    packageFiles[i] = fopen(pfile,"w");
  }

  typedef std::map<int, pair<int,int> > MyMap;
  for (int i = 0; i < nSubDFine; ++i) {
    
    char cfile[256];
    sprintf(cfile,"%s%i",collectionFile, i);
    FILE* collfile = fopen(cfile, "w");
    MyMap agglomeratedSubDMap;
    int cnt = 0;
    MyMap** pNeighborhoodNodes = new MyMap*[ subDsizeFine[i] ];
    for (int j = 0;  j < subDsizeFine[i]; ++j) {

      MyMap& neighborhoodNodes = *(pNeighborhoodNodes[j]);
      for (int k = 0; k < nSubDCoarse; ++k) {

	std::list<int> locids;
	int globNode = decFine[i][j];
	double radius = radius0;
	while (locids.size() < threshold) {
	  myKDTrees[k].collectPointsInRadius(Xfine[globNode], radius,
					     locids);
	  radius *= radiusf;
	}

	for (std::list<int>::iterator itr = locids.begin(); 
	     itr != locids.end(); ++itr) {

	  neighborhoodNodes[ decCoarse[k][*itr] ] = pair<int,int>(k, *itr);
	  agglomeratedSubDMap[ decCoarse[k][*itr] ] = pair<int,int>(k, *itr);
	}
      }     
    }
    
    // Now we know the nodes in the coarse neighborhood of
    // this fine node; now construct the tables for data 
    // transfer from coarse to fine

    // Step 1.  Reorder the data coming from the various subds. 
    // so that the subds. are contiguous.
    std::map<int,int>* reorder = new std::map<int,int>[nSubDCoarse];
    int* offsets = new int[nSubDCoarse];
    int* counts = new int[nSubDCoarse];
    int** packages = new int*[nSubDCoarse];

    memset(counts,0,sizeof(nSubDCoarse));
    for (MyMap::iterator itr = agglomeratedSubDMap.begin();
	 itr != agglomeratedSubDMap.end(); ++itr) {

      const pair<int,int>& p = itr->second;
      reorder[p.first][p.second] = counts[p.first];

      ++counts[p.first];
    }

    for (int k = 0; k < nSubDCoarse; ++k)
      packages[k] = new int[ counts[k] ];

    for (MyMap::iterator itr = agglomeratedSubDMap.begin();
	 itr != agglomeratedSubDMap.end(); ++itr) {

      const pair<int,int>& p = itr->second;
      packages[ p.first ][ reorder[p.first][p.second] ] = p.second;
    }

    for (int k = 0; k < nSubDCoarse; ++k) {

      fwrite(&counts[k], sizeof(int),1, packageFiles[k]);
      fwrite(packages[k], sizeof(int), counts[k], packageFiles[k]);
      delete [] packages[k];
    }
    delete [] packages;
    delete [] counts;
    
    offsets[0] = 0;
    for (int k = 1; i < nSubDCoarse; ++k)
      offsets[k] = offsets[k-1] + reorder[k-1].size();

    int* collection;

    for (int j = 0; j < subDsizeFine[i]; ++j) {
      MyMap& neighborhoodNodes = *(pNeighborhoodNodes[j]);
      int k = 0;
      int nss = neighborhoodNodes.size();
      fwrite(&nss, sizeof(nss),1, collfile);
      collection = new int[neighborhoodNodes.size()];
      // Rewrite itr (which is the local node 
      for (MyMap::iterator itr = neighborhoodNodes.begin(); 
	   itr != neighborhoodNodes.end(); ++itr) {
	const pair<int,int>& p = itr->second;
	collection[k] = offsets[p.first] + reorder[p.first][p.second];	
	++k;
      }

      fwrite(collection, sizeof(int), nss, collfile);
      delete [] collection;
    }
    
    for (int k = 0; k < subDsizeFine[i]; ++k) {
      delete [] pNeighborhoodNodes[k];
    }
    delete [] pNeighborhoodNodes;
    

    fclose(collfile);
  }
  
  for (int i = 0; i < nSubDCoarse; ++i) {
    fclose(packageFiles[i]);
  }
}

static void loadTopMesh(const char* fn,SVec<double,3>& nodes) {

  char word1[20],word2[20],word3[20];

  ifstream inMesh(fn);
  
  vector<double> X,Y,Z;
  int int1;
  double x,y,z;

  inMesh >> word1 >> word2;

  //load node coordinates
  while(!inMesh.eof()) {
    inMesh >> int1;
    if(inMesh.fail()) {
      inMesh.clear();
      break;
    }

    inMesh >> x >> y >> z;
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);
  }

  int nNodes = (int)X.size();
  //cout<<"Loaded "<< nNodes <<" nodes."<<endl;

  nodes.resize( nNodes );
  for (int i = 0; i < nNodes; ++i) {
    nodes[i][0] = X[i];
    nodes[i][1] = Y[i];
    nodes[i][2] = Z[i];
  }
}

static void loadDecomposition(const char* fn, int** &dec,
			      int &nSubD,int* &subDsize) {

  char line[512];
  ifstream decf(fn);
  decf.getline(line, 512);
  decf.getline(line, 512);
  decf.getline(line, 512);
  
  decf >> nSubD;

  subDsize = new int[nSubD];
  dec = new int*[nSubD];
  for (int i = 0; i < nSubD; ++i) {

    decf >> subDsize[i];
    dec[i] = new int[ subDsize[i] ];
    for (int j = 0; j < subDsize[i]; ++j)
      decf >> dec[i][j];
  }
  
}

void Multigrid::createMappingFromMeshes(IoData& ioData) {

  SVec<double,3> Xfine(5),Xcoarse(5);
  int** decFine,**decCoarse;
  int nSubDFine,nSubDCoarse;
  int* subDsizeFine,*subDsizeCoarse;

  loadTopMesh(ioData.multigrid.fineMesh,Xfine);
  loadTopMesh(ioData.multigrid.coarseMesh,Xcoarse);
  
  loadDecomposition(ioData.multigrid.fineDec,
		    decFine, nSubDFine, subDsizeFine);
  loadDecomposition(ioData.multigrid.coarseDec,
		    decCoarse, nSubDCoarse, subDsizeCoarse);
  
  constructMapping(Xfine,Xcoarse,
		   decFine,nSubDFine,subDsizeFine,
		   decCoarse,nSubDCoarse,subDsizeCoarse,
		   ioData.multigrid.packageFile,
		   ioData.multigrid.collectionFile,
		   ioData.multigrid.radius0,ioData.multigrid.threshold,
		   ioData.multigrid.radiusf); 
  
  delete [] subDsizeFine;
  delete [] subDsizeCoarse;

  for (int i = 0; i < nSubDFine; ++i)
    delete [] decFine[i];
 
  for (int i = 0; i < nSubDCoarse; ++i)
    delete [] decCoarse[i];
}
