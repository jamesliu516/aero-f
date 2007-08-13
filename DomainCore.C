#include <Domain.h>

#include <BcDef.h>
#include <TimeData.h>
#include <SubDomain.h>
#include <GeoSource.h>
#include <Vector3D.h>
#include <DistVector.h>
#include <Connectivity.h>
#include <KspPrec.h>
#include <BCApplier.h>
#include <MatchNode.h>
#include<LevelSet.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_CODES 4
#define FLUID_ID 0
#define STRUC_ID 1
#define HEAT_ID  2

//------------------------------------------------------------------------------

Domain::Domain()
{
  subDomain = 0;
  subTopo = 0;
  nodeType = 0;
  nodeFaceType = 0;

  nodeDistInfo = 0;
  edgeDistInfo = 0;
  faceDistInfo = 0;
  faceNormDistInfo = 0;
  inletNodeDistInfo = 0;

  vecPat = 0;
  compVecPat = 0;
  vec3DPat = 0;
  volPat = 0;
  levelPat = 0;
  weightPat = 0;
  edgePat = 0;
  momPat = 0;
  csPat = 0;
  fsPat = 0;
  inletVec3DPat = 0;
  inletCountPat = 0;
  inletRhsPat = 0;

  Delta = 0;
  CsDelSq = 0;
  PrT = 0;
  WCsDelSq = 0;
  WPrT = 0;
  tag = 0;
  tagBar = 0;
                                                                                                                      
  globCom = new Communicator;
  Communicator* allCom[MAX_CODES];
  globCom->split(FLUID_ID, MAX_CODES, allCom);

  com = allCom[FLUID_ID];
  timer = new Timer(com);
  com->setTimer(timer);

  strCom = allCom[STRUC_ID];
  if (strCom) {
    strTimer = new Timer(0);
    strCom->setTimer(strTimer);
  }
  else
    strTimer = 0;

  heatCom = allCom[HEAT_ID];
  if (heatCom) {
    heatTimer = new Timer(0);
    heatCom->setTimer(heatTimer);
  }
  else
    heatTimer = 0;

  meshMotionBCs = 0; //HB

}

//------------------------------------------------------------------------------

Domain::~Domain()
{
  return; //BUG omp

  if (subDomain) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) 
      if (subDomain[iSub]) delete subDomain[iSub];
    delete [] subDomain;
  }

  if (subTopo) delete subTopo;
  if (nodeType) delete [] nodeType;
  if (nodeFaceType) delete [] nodeFaceType;

  if (nodeDistInfo) delete nodeDistInfo;
  if (edgeDistInfo) delete edgeDistInfo;
  if (faceDistInfo) delete faceDistInfo;
  if (faceNormDistInfo) delete faceNormDistInfo;
  if (inletNodeDistInfo) delete inletNodeDistInfo;

  if (vecPat) delete vecPat;
  if (compVecPat) delete compVecPat;
  if (vec3DPat) delete vec3DPat;
  if (volPat) delete volPat;
  if (levelPat) delete levelPat;
  if (weightPat) delete weightPat;
  if (edgePat) delete edgePat;
  if (momPat) delete momPat;
  if (csPat) delete csPat;
  if (fsPat) delete fsPat;
  if (inletVec3DPat) delete inletVec3DPat;
  if (inletCountPat) delete inletCountPat;
  if (inletRhsPat) delete inletRhsPat;

  if (Delta) delete (Delta);
  if (CsDelSq) delete(CsDelSq);
  if (PrT) delete(PrT);
  if (WCsDelSq) delete(WCsDelSq);
  if (WPrT) delete(WPrT);
  if (tag) delete(tag);
  if (tagBar) delete(tagBar);

  //if (com) delete com;
  if(meshMotionBCs) delete meshMotionBCs; 
}

//------------------------------------------------------------------------------

void Domain::getGeometry(GeoSource &geoSource, IoData &ioData)
{
  int iSub;
  int numLocThreads = 0;
  numLocSub = geoSource.getNumLocSub();
  subDomain = new SubDomain *[numLocSub];

#ifdef _OPENMP
  numLocThreads = geoSource.getNumLocThreads();
  omp_set_num_threads(numLocThreads);
  //numLocThreads = omp_get_max_threads();
  /*
#ifdef sgi
  extern void setNumArenas(int);
  setNumArenas(numLocThreads);
#endif
  //omp_set_dynamic(0);
  */
#endif

  //set IoData in the Timer
  timer->setIoData(ioData); 
  if (strCom)
    strTimer->setIoData(ioData);

  if (com->getMaxVerbose() >= 8)
    fprintf(stdout, "CPU %d uses %d thread%s for %d subdomain%s\n", com->cpuNum(), 
	    numLocThreads, numLocThreads>1? "s":"", numLocSub, numLocSub>1? "s":"");
  else if (com->getMaxVerbose() >= 2 && com->cpuNum() == 0)
    fprintf(stdout, "%d MPI CPU for %d subdomain%s\n", com->size(),
	   geoSource.getNumGlobSub(), geoSource.getNumGlobSub()>1? "s":"");

#ifdef CXFS
  FILE *fp = fopen("subinfo.data", "w");
#endif

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub) {
    subDomain[iSub] = geoSource.getSubDomain(iSub);
    
#ifdef CXFS
    subDomain[iSub]->printInfo(fp);
#endif
  }

#ifdef CXFS
  fclose(fp);
  exit(1);
#endif

  subTopo = new SubDTopo(com->cpuNum(), geoSource.getSubToSub(), geoSource.getCpuToSub());
  volPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
  levelPat = new CommPattern<int>(subTopo, com, CommPattern<int>::CopyOnSend); // New Comm Pattern
  vec3DPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
  weightPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
  momPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
  csPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
  fsPat = new CommPattern<int>(subTopo, com, CommPattern<int>::CopyOnSend);

  int numGlobSub = geoSource.getNumGlobSub();
  int *locSubToGlobSub = (*geoSource.getCpuToSub())[com->cpuNum()];

  nodeDistInfo = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
  edgeDistInfo = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
  faceDistInfo      = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
  faceNormDistInfo  = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
  inletNodeDistInfo = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
  if (!(ioData.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
        ioData.schemes.bc.type != BoundarySchemeData::STEGER_WARMING &&
        ioData.schemes.bc.type != BoundarySchemeData::GHIDAGLIA)){
#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->markLenNull(*inletNodeDistInfo);
    }
                                                                                                  
    inletNodeDistInfo->finalize(true);
  }

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->markLenNodes(*nodeDistInfo);
    subDomain[iSub]->markLenFaces(*faceDistInfo);
    subDomain[iSub]->markLenFaceNorms(*faceNormDistInfo);
    subDomain[iSub]->setChannelNums(*subTopo);
    subDomain[iSub]->setComLenNodes(1, *volPat);
    subDomain[iSub]->setComLenNodes(1, *levelPat); // New Comm Pattern
    subDomain[iSub]->setComLenNodes(2, *csPat);
    subDomain[iSub]->setComLenNodes(2, *fsPat);
    subDomain[iSub]->setComLenNodes(3, *vec3DPat);
    subDomain[iSub]->setComLenNodes(6, *weightPat);
    subDomain[iSub]->setComLenNodes(16, *momPat);
  }

  nodeDistInfo->finalize(true);
  faceDistInfo->finalize(false);
  faceNormDistInfo->finalize(false);

  volPat->finalize();
  levelPat->finalize();
  csPat->finalize();
  fsPat->finalize();
  vec3DPat->finalize();
  weightPat->finalize();
  momPat->finalize();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->makeMasterFlag(*nodeDistInfo);

}

//------------------------------------------------------------------------------

void Domain::createVecPat(int dim, IoData *ioData)
{
  vecPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenNodes(dim, *vecPat);

  vecPat->finalize();

  if (ioData)
    if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
      compVecPat = new CommPattern<bcomp>(subTopo, com, CommPattern<bcomp>::CopyOnSend);
#pragma omp parallel for
      for (int iSub = 0; iSub<numLocSub; ++iSub)
        subDomain[iSub]->setComLenNodes(dim, *compVecPat);

      compVecPat->finalize();
    }
}

//------------------------------------------------------------------------------

void Domain::createRhsPat(int dim, IoData &ioData)
{
                                                                                                  
  if (ioData.schemes.bc.type == BoundarySchemeData::CONSTANT_EXTRAPOLATION ||
      ioData.schemes.bc.type == BoundarySchemeData::LINEAR_EXTRAPOLATION  ){

    inletRhsPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);

#pragma omp parallel for
    for ( int iSub = 0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->setComLenNodes(dim, *inletRhsPat);
    inletRhsPat->finalize();
  }
}

//------------------------------------------------------------------------------

void Domain::numberEdges()
{
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->numberEdges();

  CommPattern<int> edgeNumPat(subTopo, com, CommPattern<int>::CopyOnSend, 
			      CommPattern<int>::NonSym);

  edgePat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->identifyEdges(edgeNumPat);

  edgeNumPat.finalize();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->sndEdgeInfo(edgeNumPat);

  edgeNumPat.exchange();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->rcvEdgeInfo(edgeNumPat);
    subDomain[iSub]->setComLenEdges(4, *edgePat);
  }

  edgePat->finalize();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->markLenEdges(*edgeDistInfo);

  edgeDistInfo->finalize(false);
}

//------------------------------------------------------------------------------
// This routine makes each subdomain compute on their own the motion logic for each
// node and then exchange with the neighbors to unify this information

void Domain::setNodeType(IoData &ioData)
{
  int* facemap = reinterpret_cast<int *>(alloca(sizeof(int) * (BC_MAX_CODE - BC_MIN_CODE + 1)));
  facemap -= BC_MIN_CODE;

  int i;
  for (i=BC_MIN_CODE; i<=BC_MAX_CODE; ++i)
    facemap[i] = i;

  if (ioData.bc.wall.type == BcsWallData::ISOTHERMAL) {
    facemap[BC_ADIABATIC_WALL_MOVING] = BC_ISOTHERMAL_WALL_MOVING;
    facemap[BC_ADIABATIC_WALL_FIXED] = BC_ISOTHERMAL_WALL_FIXED;
  }

  int* bcpriority = reinterpret_cast<int *>(alloca(sizeof(int) * (BC_MAX_CODE - BC_MIN_CODE + 1)));
  bcpriority -= BC_MIN_CODE;

  bcpriority[BC_ADIABATIC_WALL_MOVING ] = 10;
  bcpriority[BC_ISOTHERMAL_WALL_MOVING] =  9;
  bcpriority[BC_SLIP_WALL_MOVING      ] =  8;
  bcpriority[BC_INLET_MOVING          ] =  7;
  bcpriority[BC_OUTLET_MOVING         ] =  6;
  bcpriority[BC_ADIABATIC_WALL_FIXED  ] =  5;
  bcpriority[BC_ISOTHERMAL_WALL_FIXED ] =  4;
  bcpriority[BC_SLIP_WALL_FIXED       ] =  3;
  bcpriority[BC_INLET_FIXED           ] =  2;
  bcpriority[BC_OUTLET_FIXED          ] =  1;
  bcpriority[BC_SYMMETRY              ] =  0;
  bcpriority[BC_INTERNAL              ] = -1;


  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setFaceType(facemap);

  nodeType = new int*[numLocSub];

  CommPattern<int> ndC(subTopo, com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenNodes(1, ndC);

  ndC.finalize();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setNodeType(bcpriority, ndC);

  ndC.exchange();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    nodeType[iSub] = subDomain[iSub]->completeNodeType(bcpriority, ndC);

  /*
  DistSVec<double,1> nt(*nodeDistInfo);
#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub) {
    double (*_nt)[1] = nt.subData(iSub);
    for (int i=0; i<nt.subSize(iSub); ++i)
      _nt[i][0] = static_cast<double>(nodeType[iSub][i]);
  }
  writeVectorToFile("nodetype", 0, 0.0, nt);
  */
  // Now we find how many slip surfaces we have
  int numSlipSurfs = 0;
  map<int,SurfaceData *> &surfaceMap = ioData.surfaces.surfaceMap.dataMap;
  map<int,SurfaceData *>::iterator it = surfaceMap.begin();
  SurfaceData *slipSurface[sizeof(int)];

  while(it != surfaceMap.end()) {
     if(it->second->nx != 0.0 || it->second->ny != 0.0 || it->second->nz != 0.0) {
        com->fprintf(stderr," ... surface %3d is sliding with normal = %e %e %e\n",it->first,
                     it->second->nx,it->second->ny,it->second->nz);
        it->second->setBit(1<<numSlipSurfs);
	slipSurface[numSlipSurfs] = it->second;
	numSlipSurfs++;
	if(numSlipSurfs > sizeof(int) ) {
	   com->fprintf(stderr, " *** ERROR: There are more sliding surfaces than the maximum allowed (%d)\n",
	       sizeof(int));
	   throw "ERROR";
	}
     }
     it++;
  }

  if(numSlipSurfs)  
    com->fprintf(stderr," ... There is %d sliding surfaces.\n",numSlipSurfs); 
   
  int **slipSurfOwn = new int*[numLocSub];

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    slipSurfOwn[iSub] = subDomain[iSub]->getSlipSurfOwnership(ndC, surfaceMap);

  ndC.exchange();

  meshMotionBCs = new BCApplier(numLocSub,this,ioData);
#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->createSlipSurfProjection(slipSurfOwn[iSub], ndC, meshMotionBCs, slipSurface);

  if(slipSurfOwn){
#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; ++iSub)
      if(slipSurfOwn[iSub]) { delete [] slipSurfOwn[iSub]; slipSurfOwn[iSub] = 0; }
    delete [] slipSurfOwn; slipSurfOwn = 0;
  }
  
/* ARL: creation of nodeFaceType which 
 *      tells which kind of faces a node
 *      is connected to.
 */

  nodeFaceType = new int*[numLocSub];
  CommPattern<int> ndCFace(subTopo, com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenNodes(1, ndCFace);

  ndCFace.finalize();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setNodeFaceType(ndCFace);

  ndCFace.exchange();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    nodeFaceType[iSub] = subDomain[iSub]->completeNodeFaceType(ndCFace);

}
//------------------------------------------------------------------------------

void Domain::setInletNodes(IoData &ioData)
{
// set up the inlet nodes if extrapolation is to be used for boundary treatment
  int iSub;
  if (ioData.schemes.bc.type == BoundarySchemeData::CONSTANT_EXTRAPOLATION ||
      ioData.schemes.bc.type == BoundarySchemeData::LINEAR_EXTRAPOLATION){

    //creation of inletNodeDistInfo and inletNodes
    //inletNodeDistInfo = new DistInfo(numLocThreads, numLocSub, numGlobSub, locSubToGlobSub, com);
    //  this part done in the creation of the domain.
#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->setInletNodes(ioData);
      subDomain[iSub]->markLenInletNodes(*inletNodeDistInfo);
    }

    if (ioData.schemes.bc.type == BoundarySchemeData::LINEAR_EXTRAPOLATION){
#pragma omp parallel for
      for (iSub = 0; iSub<numLocSub; ++iSub) 
        subDomain[iSub]->setInletNodes2(ioData);
    }

    inletNodeDistInfo->finalize(false);

//creation of connectivity sharedInletNodes
#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->createSharedInletNodeConnectivity(iSub);

//creation of the Communicator patterns needed.
    inletVec3DPat = new CommPattern<double>(subTopo, com, CommPattern<double>::CopyOnSend);
    inletCountPat = new CommPattern<int>(subTopo, com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; iSub++){
      subDomain[iSub]->setComLenInletNodes(3, *inletVec3DPat);
      subDomain[iSub]->setComLenInletNodes(1, *inletCountPat);
    }
    inletVec3DPat->finalize();
    inletCountPat->finalize();
  }
}

//------------------------------------------------------------------------------

void Domain::makeRotationOwnership(IoData &ioData)  {

  // Find numbers of "rotating" surfaces & "translating" walls (only needed for output display)
  map<int,SurfaceData *> &surfaceMap = ioData.surfaces.surfaceMap.dataMap;
  map<int,RotationData*> &rotationMap= ioData.rotations.rotationMap.dataMap;
  map<int,SurfaceData *>::iterator it = surfaceMap.begin();
  int numRotSurfs  = 0;
  int numTransWalls= 0;
 
  while(it != surfaceMap.end()) {
    map<int,RotationData*>::iterator it1 = rotationMap.find(it->second->rotationID); 
    if(it1!=rotationMap.end()) {
      if(it1->second->infRadius) {
        numTransWalls++;
        com->fprintf(stderr," ... surface %2d is ``translating''\n",it->first, it1->first);
        com->fprintf(stderr,"     -> uniform velocity V = %3.2e in direction %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      } else {
        numRotSurfs++;
        com->fprintf(stderr," ... surface %2d is ``rotating'' using rotation data %2d\n",it->first, it1->first);
        com->fprintf(stderr,"     -> omega = %3.2e, rotation axis = %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      }
    } 
    it++;
  }
  if(numRotSurfs)  
    com->fprintf(stderr," ... There is %d ``rotating'' surfaces.\n",numRotSurfs); 
  if(numTransWalls)  
    com->fprintf(stderr," ... There is %d ``translating'' surfaces.\n",numTransWalls); 

  CommPattern<int> ndC(subTopo, com, CommPattern<int>::CopyOnSend);

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenNodes(1, ndC);

  ndC.finalize();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->getRotSurfaceOwnership(ndC, surfaceMap);

  ndC.exchange();

#pragma omp parallel for
  for (iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->completeRotateSurfaceOwnership(ndC);
}

//------------------------------------------------------------------------------

void Domain::setFaceToElementConnectivity()
{
  int nswap = 0;
#pragma omp parallel for reduction(+: nswap)
  for (int iSub=0; iSub<numLocSub; ++iSub)
    nswap += subDomain[iSub]->setFaceToElementConnectivity();

  com->globalSum(1, &nswap);
  if (nswap > 0)
    com->printf(1, "*** Warning: changed the orientation of %d boundary face%s\n", 
		nswap, nswap>1 ? "s":"");
}

//------------------------------------------------------------------------------

void Domain::printElementStatistics()
{
  int iSub;
  int (*num)[4] = reinterpret_cast<int (*)[4]>(alloca(sizeof(int) * numLocSub * 4));

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->getElementStatistics(num[iSub][0], num[iSub][1], num[iSub][2], num[iSub][3]);

  int minNum[4], maxNum[4], totNum[4];

  int k;
  for (k=0; k<4; ++k) {
    minNum[k] = num[0][k];
    maxNum[k] = num[0][k];
    totNum[k] = num[0][k];
  }
    
  for (iSub=1; iSub<numLocSub; ++iSub) {
    for (k=0; k<4; ++k) {
      minNum[k] = min(minNum[k], num[iSub][k]);
      maxNum[k] = max(maxNum[k], num[iSub][k]);
      totNum[k] += num[iSub][k];
    }
  }
    
  com->globalMin(4, minNum);
  com->globalMax(4, maxNum);
  com->globalSum(4, totNum);

  com->printf(2, "Node statistics: min=%d, max=%d, total=%d\n", 
	      minNum[0], maxNum[0], totNum[0]);
  com->printf(2, "Edge statistics: min=%d, max=%d, total=%d\n", 
	      minNum[1], maxNum[1], totNum[1]);
  com->printf(2, "Face statistics: min=%d, max=%d, total=%d\n", 
	      minNum[2], maxNum[2], totNum[2]);
  com->printf(2, "Elem statistics: min=%d, max=%d, total=%d\n", 
	      minNum[3], maxNum[3], totNum[3]);
}

//------------------------------------------------------------------------------

int Domain::computeControlVolumes(double lscale, DistSVec<double,3> &X, DistVec<double> &ctrlVol)
{
  int iSub, ierr = 0;

#pragma omp parallel for reduction(+: ierr)
  for (iSub=0; iSub<numLocSub; ++iSub) {
    ierr += subDomain[iSub]->computeControlVolumes(0, lscale, X(iSub), ctrlVol(iSub));
    double (*locctrlVol)[1] = reinterpret_cast<double (*)[1]>(ctrlVol.subData(iSub));
    subDomain[iSub]->sndData(*volPat, locctrlVol);
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locctrlVol)[1] = reinterpret_cast<double (*)[1]>(ctrlVol.subData(iSub));
    subDomain[iSub]->addRcvData(*volPat, locctrlVol);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    com->fprintf(stderr, "*** Error: %d negative volume%s\n", ierr, ierr>1? "s":"");

#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->computeControlVolumes(ierr, lscale, X(iSub), ctrlVol(iSub));

    exit(1);
  }

  return ierr;
}

//------------------------------------------------------------------------------

void Domain::computeFaceNormals(DistSVec<double,3> &X, DistVec<Vec3D> &faceNorm)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->computeFaceNormals(X(iSub), faceNorm(iSub));
}

//------------------------------------------------------------------------------

void Domain::computeNormalsGCL1(DistSVec<double,3> &Xn, DistSVec<double,3> &Xnp1, 
				DistSVec<double,3> &Xdot,
				DistVec<Vec3D> &edgeNorm, DistVec<double> &edgeNormVel,
				DistVec<Vec3D> &faceNorm, DistVec<double> &faceNormVel)
{
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeNormalsGCL1(Xn(iSub), Xnp1(iSub), Xdot(iSub), 
					edgeNorm(iSub), edgeNormVel(iSub), 
					faceNorm(iSub), faceNormVel(iSub));
    subDomain[iSub]->sndNormals(*edgePat, edgeNorm.subData(iSub), 
				edgeNormVel.subData(iSub));
  }

  edgePat->exchange();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->rcvNormals(*edgePat, edgeNorm.subData(iSub), 
				edgeNormVel.subData(iSub));
}

//------------------------------------------------------------------------------

void Domain::computeNormalsGCL2(TimeData &data, DistVec<Vec3D> &edgeNorm, 
				DistVec<Vec3D> &edgeNorm_old, DistVec<double> &edgeNormVel, 
				DistVec<double> &edgeNormVel_old, DistVec<Vec3D> &faceNorm, 
				DistVec<Vec3D> &faceNorm_old, DistVec<double> &faceNormVel, 
				DistVec<double> &faceNormVel_old)
{
  if (data.exist_nm1) {
    double cnp1 = data.alpha_np1;
    double cnm1 = (data.alpha_np1 + data.alpha_n) / data.tau_n;

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      Vec3D *edgenorm = edgeNorm.subData(iSub);
      Vec3D *edgenorm_old = edgeNorm_old.subData(iSub);
      double *edgenormvel = edgeNormVel.subData(iSub);
      double *edgenormvel_old = edgeNormVel_old.subData(iSub);
      Vec3D *facenorm = faceNorm.subData(iSub);
      Vec3D *facenorm_old = faceNorm_old.subData(iSub);
      double *facenormvel = faceNormVel.subData(iSub);
      double *facenormvel_old = faceNormVel_old.subData(iSub);

      double tmpvel;
      Vec3D tmpnorm;

      int i;
      for (i=0; i<edgeNorm.subSize(iSub); ++i) {
	tmpnorm = edgenorm[i];
	edgenorm[i] = cnp1 * tmpnorm + cnm1 * edgenorm_old[i];
	edgenorm_old[i] = tmpnorm;

	tmpvel = edgenormvel[i];
	edgenormvel[i] = cnp1 * tmpvel + cnm1 * edgenormvel_old[i];
	edgenormvel_old[i] = tmpvel;
      }
      for (i=0; i<faceNorm.subSize(iSub); ++i) {
	tmpnorm = facenorm[i];
	facenorm[i] = cnp1 * tmpnorm + cnm1 * facenorm_old[i];
	facenorm_old[i] = tmpnorm;

	tmpvel = facenormvel[i];
	facenormvel[i] = cnp1 * tmpvel + cnm1 * facenormvel_old[i];
	facenormvel_old[i] = tmpvel;
      }
    }
  } 
  else {
    edgeNorm_old = edgeNorm;
    edgeNormVel_old = edgeNormVel;
    faceNorm_old = faceNorm;
    faceNormVel_old = faceNormVel;
  }
}

//------------------------------------------------------------------------------

void Domain::computeNormalsEZGCL1(double oodt, DistSVec<double,3>& Xn, DistSVec<double,3>& Xnp1, 
				  DistVec<Vec3D>& edgeNorm, DistVec<double>& edgeNormVel,
				  DistVec<Vec3D>& faceNorm, DistVec<double>& faceNormVel)
{
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeNormalsEZGCL1(oodt, Xn(iSub), Xnp1(iSub),
					  edgeNorm(iSub), edgeNormVel(iSub), 
					  faceNorm(iSub), faceNormVel(iSub));
    subDomain[iSub]->sndNormals(*edgePat, edgeNorm.subData(iSub), 
				edgeNormVel.subData(iSub));
  }

  edgePat->exchange();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->rcvNormals(*edgePat, edgeNorm.subData(iSub), 
				edgeNormVel.subData(iSub));
}

//------------------------------------------------------------------------------

void Domain::computeNormalsEZGCL2(TimeData& data, DistVec<double>& edgeNormVel, 
				  DistVec<double>& edgeNormVel_old, 
				  DistVec<double>& faceNormVel, 
				  DistVec<double>& faceNormVel_old)
{
  if (data.exist_nm1) {
    double cnp1 = data.alpha_np1;
    double cnm1 = (data.alpha_np1 + data.alpha_n) / data.tau_n;

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double *edgenormvel = edgeNormVel.subData(iSub);
      double *edgenormvel_old = edgeNormVel_old.subData(iSub);
      double *facenormvel = faceNormVel.subData(iSub);
      double *facenormvel_old = faceNormVel_old.subData(iSub);

      int i;
      for (i=0; i<edgeNormVel.subSize(iSub); ++i) {
	double tmpvel = edgenormvel[i];
	edgenormvel[i] = cnp1 * tmpvel + cnm1 * edgenormvel_old[i];
	edgenormvel_old[i] = tmpvel;
      }
      for (i=0; i<faceNormVel.subSize(iSub); ++i) {
	double tmpvel = facenormvel[i];
	facenormvel[i] = cnp1 * tmpvel + cnm1 * facenormvel_old[i];
	facenormvel_old[i] = tmpvel;
      }
    }
  } 
  else {
    edgeNormVel_old = edgeNormVel;
    faceNormVel_old = faceNormVel;
  }
}

//------------------------------------------------------------------------------

void Domain::computeNormalsEZGCL3(TimeData& data, DistVec<double>& edgeNormVel, 
				  DistVec<double>& edgeNormVel_nm1, 
				  DistVec<double>& edgeNormVel_nm2,
				  DistVec<double>& faceNormVel, 
				  DistVec<double>& faceNormVel_nm1,
				  DistVec<double>& faceNormVel_nm2)
{
  if (data.exist_nm1 && data.exist_nm2) {
    double cnp1 = data.alpha_np1;
    double cn = (data.alpha_np1 + data.alpha_n) / data.tau_n;
    double cnm1 = (data.alpha_np1 + data.alpha_n + data.alpha_nm1) * data.dt_nm2 / data.dt_n;

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      double* evel = edgeNormVel.subData(iSub);
      double* evel_nm1 = edgeNormVel_nm1.subData(iSub);
      double* evel_nm2 = edgeNormVel_nm2.subData(iSub);
      double* fvel = faceNormVel.subData(iSub);
      double* fvel_nm1 = faceNormVel_nm1.subData(iSub);
      double* fvel_nm2 = faceNormVel_nm2.subData(iSub);

      int i;
      for (i=0; i<edgeNormVel.subSize(iSub); ++i) {
	double tmpvel = evel[i];
	evel[i] = cnp1*tmpvel + cn*evel_nm1[i] + cnm1*evel_nm2[i];
	evel_nm2[i] = evel_nm1[i];
	evel_nm1[i] = tmpvel;
      }
      for (i=0; i<faceNormVel.subSize(iSub); ++i) {
	double tmpvel = fvel[i];
	fvel[i] = cnp1*tmpvel + cn*fvel_nm1[i] + cnm1*fvel_nm2[i];
	fvel_nm2[i] = fvel_nm1[i];
	fvel_nm1[i] = tmpvel;
      }
    }
  } 
  else {
    if (data.exist_nm1) {
      edgeNormVel_nm2 = edgeNormVel_nm1;
      faceNormVel_nm2 = faceNormVel_nm1;
    }
    computeNormalsEZGCL2(data, edgeNormVel, edgeNormVel_nm1, faceNormVel, faceNormVel_nm1);
  }
}

//------------------------------------------------------------------------------

void Domain::testNormals(DistVec<Vec3D> &edgeNorm, DistVec<double> &edgeNormVel,
			 DistVec<Vec3D> &faceNorm, DistVec<double> &faceNormVel)
{
  int iSub = 0;
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->testNormals(edgeNorm(iSub), edgeNormVel(iSub), 
				 faceNorm(iSub), faceNormVel(iSub));

  com->barrier();

  exit(1);

}

//------------------------------------------------------------------------------
                                                                                                  
void Domain::computeInletNormals(DistVec<Vec3D>& inletNodeNorm, DistVec<Vec3D>& faceNorm, DistVec<int>& numFaceNeighb)
{
                                                                                                  
  if (inletRhsPat){     //better to put inletVec3DPat????
                                                                                                  
        //first we sum up the normals
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub){
      subDomain[iSub]->sumInletNormals(inletNodeNorm(iSub), faceNorm(iSub), numFaceNeighb(iSub));
      subDomain[iSub]->sndInletData(*inletVec3DPat,
                reinterpret_cast<double(*)[3]>(inletNodeNorm.subData(iSub)));
      subDomain[iSub]->sndInletData(*inletCountPat,
                reinterpret_cast<int(*)[1]>(numFaceNeighb.subData(iSub)));
    }
                                                                                                  
    inletVec3DPat->exchange();
    inletCountPat->exchange();
                                                                                                  
                                                                                                  
#pragma omp parallel for
    for ( int iSub = 0; iSub<numLocSub; ++iSub){
      subDomain[iSub]->addRcvInletData(*inletVec3DPat,
                reinterpret_cast<double(*)[3]>(inletNodeNorm.subData(iSub)));
      subDomain[iSub]->addRcvInletData(*inletCountPat,
                reinterpret_cast<int(*)[1]>(numFaceNeighb.subData(iSub)));
    }
        //then we divide
#pragma omp parallel for
    for( int iSub = 0; iSub<numLocSub; ++iSub)  {
      subDomain[iSub]->numDivideNormals(inletNodeNorm(iSub), numFaceNeighb(iSub));
    }
  }
}

//------------------------------------------------------------------------------

void Domain::computeVelocities(ImplicitData::Velocities typeVel, TimeData &timeData, 
			       DistSVec<double,3> &Xsdot, DistSVec<double,3> &Xnm1, 
			       DistSVec<double,3> &Xn, DistSVec<double,3> &X, 
			       DistSVec<double,3> &Xdot)
{

  timeData.computeVelocities(typeVel, Xnm1, Xn, X, Xdot);

  if (typeVel == ImplicitData::IMPOSED_VEL) {
    Xdot = Xsdot;
  }
  else if (typeVel == ImplicitData::IMPOSED_BACKWARD_EULER_VEL ||
	   typeVel == ImplicitData::IMPOSED_THREE_POINT_BDF_VEL) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      int *type = nodeType[iSub];
      double (*xsdot)[3] = Xsdot.subData(iSub);
      double (*xdot)[3] = Xdot.subData(iSub);

      for (int i=0; i<Xdot.subSize(iSub); ++i) {
	if (type[i] < BC_INTERNAL) {
	  xdot[i][0] = xsdot[i][0];
	  xdot[i][1] = xsdot[i][1];
	  xdot[i][2] = xsdot[i][2];
	}
      }
    }
  }
  Xsdot = Xdot;

}

//------------------------------------------------------------------------------

void Domain::computeWeightsLeastSquares(DistSVec<double,3> &X, DistSVec<double,6> &R)
{

  int iSub;
  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeWeightsLeastSquaresEdgePart(X(iSub), R(iSub));
    subDomain[iSub]->sndData(*weightPat, R.subData(iSub));
  }

  weightPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*weightPat, R.subData(iSub));
    subDomain[iSub]->computeWeightsLeastSquaresNodePart(R(iSub));
  }

  timer->addNodalWeightsTime(t0);

}

//------------------------------------------------------------------------------
//  least square gradient involving only nodes of same fluid (multiphase flow)
void Domain::computeWeightsLeastSquares(DistSVec<double,3> &X, DistVec<double> &Phi,
                                        DistSVec<double,6> &R)
{

  int iSub;
  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeWeightsLeastSquaresEdgePart(X(iSub), Phi(iSub), R(iSub));
    subDomain[iSub]->sndData(*weightPat, R.subData(iSub));
  }

  weightPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*weightPat, R.subData(iSub));
    subDomain[iSub]->computeWeightsLeastSquaresNodePart(R(iSub));
  }

  timer->addNodalWeightsTime(t0);

}

//------------------------------------------------------------------------------

void Domain::computeWeightsGalerkin(DistSVec<double,3> &X, DistSVec<double,3> &wii, 
				    DistSVec<double,3> &wij, DistSVec<double,3> &wji)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->computeWeightsGalerkin(X(iSub), wii(iSub), wij(iSub), wji(iSub));

  timer->addNodalWeightsTime(t0);

}

//------------------------------------------------------------------------------

void Domain::getReferenceMeshPosition(DistSVec<double,3> &x)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->getReferenceMeshPosition(x(iSub));

}

//------------------------------------------------------------------------------
//HB: modified to return only the matched nodes as interface nodes
void Domain::getNdAeroLists(int *&nInterfNd, int **&interfNd, int *&nInfNd, 
			    int **&infNd, int *&nInternalNd, int **&internalNd, MatchNodeSet** matchNodes)
{

  nInterfNd   = new int[numLocSub];
  interfNd    = new int *[numLocSub];
  nInfNd      = new int[numLocSub];
  infNd       = new int *[numLocSub];
  nInternalNd = new int[numLocSub];
  internalNd  = new int *[numLocSub];

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    MatchNodeSet* subMatchNodes = (matchNodes) ? matchNodes[iSub] : 0;
    subDomain[iSub]->getNdAeroLists(nInterfNd[iSub], interfNd[iSub], nInfNd[iSub], 
				    infNd[iSub], nInternalNd[iSub], internalNd[iSub], subMatchNodes);
  }
}

//------------------------------------------------------------------------------

void Domain::applySmoothing(DistVec<double> &ctrlVol, DistVec<double> &Q)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(Q.subData(iSub)));

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(Q.subData(iSub)));

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applySmoothing(ctrlVol(iSub), Q(iSub));

}
//------------------------------------------------------------------------------

void Domain::computeDelRatios(DistMacroCellSet *macroCells, DistVec<double> &ctrlVol, int scopeDepth,
                              double *rmax, double *rmin, double *rsum, double *rrsum, int *numNodes)
{
  macroCells->computeDelRatios(ctrlVol, scopeDepth, rmin, rmax, rsum, rrsum, numNodes);
}
//------------------------------------------------------------------------------

void Domain::computeDynamicLESTerm(DynamicLESTerm *dles, DistSVec<double,2> &CsDeltaSq,
                                   DistSVec<double,3> &X, DistVec<double> &Cs, DistVec<double> &VolSum)
{
#pragma omp parallel for
   for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->computeDynamicLESTerm(dles, CsDeltaSq(iSub), X(iSub), Cs(iSub), VolSum(iSub));
}

//------------------------------------------------------------------------------
/*
double Domain::reinitLS(DistSVec<double,3> &X, DistVec<double> &Phi, DistSVec<double,dim> &U, int iti)
{
  int iSub;
  Communicator *com = getCommunicator();
  double (*rb) = reinterpret_cast<double (*)>(alloca(sizeof(double) * numLocSub));                                                                                                   
  double rbg  = 0.0;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    rb[iSub]  = subDomain[iSub]->reinitLS(X(iSub),Phi(iSub), U(iSub), iti);
    rbg       = max(rbg, rb[iSub]);
  }
                                                                                                  
  com->globalMax(1, &rbg);
                                                                                                  
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->reinitLS(X(iSub),Phi(iSub), U(iSub), rbg, iti);
                                                                                                  
  return rbg;
}
*/
//------------------------------------------------------------------------------

void Domain::solveLSequation(DistSVec<double,3> &X, DistVec<double> &Phi, 
			DistSVec<double,6> &ddx,DistSVec<double,6> &ddy, 
			DistSVec<double,6> &ddz, DistVec<double> &PhiF)
{
  int iSub;
  Communicator *com = getCommunicator();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->SolveLSequation(X(iSub), Phi(iSub),ddx(iSub), ddy(iSub), ddz(iSub), PhiF(iSub));
}
//------------------------------------------------------------------------------
void Domain::solveLS(DistVec<double> &b, DistVec<double> &dPhi, double dt)
{
  int iSub;
  Communicator *com = getCommunicator();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->solveLS(b(iSub), dPhi(iSub), dt);
}

//------------------------------------------------------------------------------
void Domain::printPhi(DistSVec<double,3> &X, DistVec<double> &Phi, int it)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printPhi(X(iSub), Phi(iSub), numLocSub);
  com->barrier();
}


// ------------------------------------------------------------------------------------------
