/*
 * CrackingSurface.h
 *
 *  Created on: Feb 2, 2011 (the night before the Year of Rabbit)
 *  Author: Kevin Wang
 */

#ifndef CRACKINGSURFACE_H_
#define CRACKINGSURFACE_H_
#include<LOCAL_LEVELSET.h>
#include<map>
#include<set>

//------------------------------------------------------------------------------

struct PhantomElement {
  int nNodes;
  double *phi;
  int *nodes;

  // constructors
  PhantomElement(): nNodes(-1), phi(0), nodes(0) {}
  PhantomElement(int n, int* nod, double* ph);
  PhantomElement(int a, int b, int c, int d, 
                 double phia, double phib, double phic, double phid);
  // destructor
  ~PhantomElement() {if(phi) delete[] phi;  if(nodes) delete[] nodes;}
};

//------------------------------------------------------------------------------

struct LatestCracking {
  std::set<int> phantomQuads;
  std::map<int,int> phantomNodes; //Note: "phantomNodes" are NOT equivalent to "nodes of phantomQuads"!
  LatestCracking() {/*nothing*/}
};

//------------------------------------------------------------------------------

class CrackingSurface : public LocalLevelSet {
  const int elemType; //currently only support quadrangles.
  int nTotalQuads, nUsedQuads;
  int nTotalTrias, nUsedTrias;
  int nTotalNodes, nUsedNodes;

  std::map<int,PhantomElement> phantoms; //size: number of cracked (quad) elements
  LatestCracking latest;
//  std::set<int> latestCrackedQuads;
  int (*tria2quad)[2]; //size: nTotalTrias
  int (*quad2tria)[2]; //size: nTotalQuads
  bool *cracked; //size: nTotalQuads
   
public:
  CrackingSurface(int eType, int nUsed, int nTotal, int nUsedNd, int nTotNodes);
  ~CrackingSurface() {if(tria2quad) delete[] tria2quad; if(cracked) delete[] cracked;}

  //called by EmbeddedStructure only!
  int splitQuads(int* quadTopo, int nQuads, int(*triaTopo)[3]);
  int updateCracking(int nNew, int* newPhan, double* phi, int(*triaTopo)[3], int nUsedNd);

  int numCrackedElements() {return phantoms.size();}
  bool hasCracked(int trId);
  double getPhi(int trId, double xi1, double xi2, bool* hasCracked=0);

  int totNodes()  const {return nTotalNodes;}
  int usedNodes() const {return nUsedNodes;}
  int totTrias()  const {return nTotalTrias;}
  int usedTrias() const {return nUsedTrias;}
  std::set<int> getLatestPhantomQuads() const {return latest.phantomQuads;}
  std::map<int,int> getLatestPhantomNodes() const {return latest.phantomNodes;}
  void getQuad2Tria(int quad, int &trId1, int &trId2) {trId1=quad2tria[quad][0]; trId2=quad2tria[quad][1];}
};

//------------------------------------------------------------------------------
#endif /* CRACKINGSURFACE_H_ */

