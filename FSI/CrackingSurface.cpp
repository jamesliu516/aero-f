/*
 * CrackingSurface.h
 *
 *  Created on: Feb 9, 2011
 *  Author: Kevin Wang
 */

#include<FSI/CrackingSurface.h>
using std::map;
using std::set;

//------------------------------------------------------------------------------

PhantomElement::PhantomElement(int n, int* nod, double* ph): nNodes(n) {
  phi = new double[n];
  nodes = new int[n];
  for(int i=0; i<n; i++) {
    phi[i] = ph[i];  nodes[i] = nod[i];}
}

//------------------------------------------------------------------------------

PhantomElement::PhantomElement(int a, int b, int c, int d,
               double phia, double phib, double phic, double phid): nNodes(4) {
  phi = new double[4];
  nodes = new int[4];
  phi[0] = phia; phi[1] = phib; phi[2] = phic; phi[3] = phid;
  nodes[0] = a; nodes[1] = b; nodes[2] = c; nodes[3] = d;
}

//------------------------------------------------------------------------------

CrackingSurface::CrackingSurface(int eType, int nUsed, int nTotal, int nUsedNd, int nTotNodes): elemType(eType)
{
  if(eType!=4) {fprintf(stderr,"ERROR: ElemType(%d) is not supported for CrackingSurface!\n");exit(1);}
  nTotalNodes = nTotNodes; nUsedNodes = nUsedNd;
  nTotalQuads = nTotal;    nUsedQuads = nUsed;
  nTotalTrias = nTotal*2,  nUsedTrias = 0; //surface not splitted yet!

  tria2quad = new int[nTotalTrias][2];
  quad2tria = new int[nTotalQuads][2];
  cracked = new bool[nTotalQuads];

  for(int i=0; i<nTotalTrias; i++)
    tria2quad[i][0] = tria2quad[i][1] = -1;
  for(int i=0; i<nTotalQuads; i++) {
    cracked[i] = false;
    quad2tria[i][0] = quad2tria[i][1] = -1;
  }
}

//------------------------------------------------------------------------------

CrackingSurface::~CrackingSurface()
{
  for(map<int,PhantomElement*>::iterator it=phantoms.begin(); it!=phantoms.end(); it++)
    delete it->second;

  if(tria2quad) delete[] tria2quad; 
  if(cracked) delete[] cracked;
}

//------------------------------------------------------------------------------

int CrackingSurface::splitQuads(int* quadTopo, int nQuads, int(*triaTopo)[3])
{
  if(nQuads!=nUsedQuads) {fprintf(stderr,"Software bug in CrackingSurface::splitQuads!\n");exit(-1);}

  int count = 0;
  for(int i=0; i<nQuads; i++) {
    triaTopo[count][0] = quadTopo[i*4];
    triaTopo[count][1] = quadTopo[i*4+1];
    triaTopo[count][2] = quadTopo[i*4+2];
    tria2quad[count][0] = i;
    tria2quad[count][1] = 0;
    quad2tria[i][0] = count;
    count++;

    if(quadTopo[i*4+2]==quadTopo[i*4+3]) { //this is a degenerated triangle
      quad2tria[i][1] = -1;
      continue;
    }

    triaTopo[count][0] = quadTopo[i*4];
    triaTopo[count][1] = quadTopo[i*4+2];
    triaTopo[count][2] = quadTopo[i*4+3];
    tria2quad[count][0] = i;
    tria2quad[count][1] = 1;
    quad2tria[i][1] = count;
    count++;
  }

  nUsedTrias = count;

  for(int i=nQuads; i<nTotalQuads; i++) {
    tria2quad[count][0] = i; tria2quad[count][1] = 0;
    quad2tria[i][0] = count;
    count++;
    tria2quad[count][0] = i; tria2quad[count][1] = 1;
    quad2tria[i][1] = count;
    count++;
  }

  for(int i=nUsedTrias; i<nTotalTrias; i++)
    triaTopo[i][0] = triaTopo[i][1] = triaTopo[i][2] = -1;

  return nUsedTrias;
}

//------------------------------------------------------------------------------

int CrackingSurface::updateCracking(int nNew, int* newPhan, double* phi, int(*triaTopo)[3], int nUsedNd)
{
  if(nNew==0)
    return 0;

  latest.phantomQuads.clear();
  latest.phantomNodes.clear(); //flush
  int quadId,maxtrId=0,count=0;

  for(int i=0; i<2*nNew; i++) {
    quadId = newPhan[5*i];
    cracked[quadId] = true;
    latest.phantomQuads.insert(quadId);
    //insert a new phantom element
    PhantomElement *thisone = new PhantomElement(newPhan[5*i+1],newPhan[5*i+2],newPhan[5*i+3],newPhan[5*i+4],
                                                 phi[4*i],phi[4*i+1],phi[4*i+2],phi[4*i+3]);
    phantoms[quadId] = thisone; 
    //modify the triangle mesh connectivity
    int trId1, trId2;
    if(quadId>=nUsedQuads) { //this is a new quad
      if(quadId>=nUsedQuads+nNew) {
        fprintf(stderr,"ERROR: nUsed = %d, nNew = %d, currentId = %d!\n", nUsedQuads, nNew, quadId+1);exit(-1);}
      count++;
      trId1 = nUsedTrias + 2*(quadId-nUsedQuads);
      trId2 = trId1 + 1;
      if(maxtrId<trId2) maxtrId = trId2;
    } else { //this is an existing quad
      trId1 = quad2tria[quadId][0];
      trId2 = quad2tria[quadId][1];
      if(trId2<0) {fprintf(stderr,"SOFTWARE BUG: CAPTURED TRIANGLE ID %d for QUAD %D\n", trId2+1, quadId+1);exit(-1);}
    }

    triaTopo[trId1][0] = newPhan[5*i+1];
    triaTopo[trId1][1] = newPhan[5*i+2];
    triaTopo[trId1][2] = newPhan[5*i+3];
    triaTopo[trId2][0] = newPhan[5*i+1];
    triaTopo[trId2][1] = newPhan[5*i+3];
    triaTopo[trId2][2] = newPhan[5*i+4];
  }

  //construct latest.phantomNodes.
  for(int i=0; i<2*nNew; i+=2) {
    int j=i+1;
    for(int myNode=1; myNode<5; myNode++)
      if(newPhan[5*i+myNode]>=nUsedNodes && newPhan[5*j+myNode]<nUsedNodes) //newPhan[5*i+myNode] is new
        latest.phantomNodes[newPhan[5*i+myNode]] = newPhan[5*j+myNode];
      else if(newPhan[5*j+myNode]>=nUsedNodes && newPhan[5*i+myNode]<nUsedNodes) //newPhan[5*j+myNode] is new
        latest.phantomNodes[newPhan[5*j+myNode]] = newPhan[5*i+myNode];
      else {fprintf(stderr,"SOFTWARE BUG: Node numbering is wrong! Aborting.\n");exit(-1);}
  }

  if(count!=nNew) {fprintf(stderr,"SOFTWARE BUG: Inconsistent on the number of new quads(%d v.s. %d)!\n", nNew, count);exit(-1);}
  nUsedQuads += nNew;
  if(maxtrId+1!=nUsedTrias+2*nNew) {
    fprintf(stderr,"SOFTWARE BUG: Violated the ordering of new elements (%d v.s. %d)\n", maxtrId+1,nUsedTrias+2*nNew);exit(-1);}
  nUsedTrias += 2*nNew;

  nUsedNodes = nUsedNd;

  return 2*nNew;
}

//------------------------------------------------------------------------------

bool CrackingSurface::hasCracked(int trId)
{
  if(trId>=nUsedTrias) {
    fprintf(stderr,"ERROR: Unable to access Triangle %d of the embedded surface(%d)!\n", trId+1, nUsedTrias);
    exit(-1);
  }
  if(cracked[tria2quad[trId][0]])
    return true;
  else
    return false;
}

//------------------------------------------------------------------------------

double CrackingSurface::getPhi(int trId, double xi1, double xi2, bool* hasCracked)
{
  if(trId>=nUsedTrias) {
    fprintf(stderr,"ERROR: Unable to access Triangle %d of the embedded surface(%d)!\n", trId+1, nUsedTrias);
    exit(-1);
  }

  if(!cracked[tria2quad[trId][0]]) { //no cracking :)
    if(hasCracked) *hasCracked = false;
    return 1.0;
  }

  // This element really cracked.
  if(hasCracked) *hasCracked = true;
  if(phantoms.find(tria2quad[trId][0])==phantoms.end()) {
    fprintf(stderr,"ERROR:Triangle %d (in Quad %d) contains no cracking!\n",trId,tria2quad[trId][0]);exit(-1);}

  double *phi = phantoms[tria2quad[trId][0]]->phi;
  double xi3;
  switch (tria2quad[trId][1]) {
    case 0: // This triangle is ABC
      xi3 = 1.0 - xi1 - xi2;
//      fprintf(stderr,"Now in getPhi! input:(%d,%e,%e), Quad %d, phi = %e\n", trId+1, xi1, xi2, tria2quad[trId][0]+1, phi[0]*xi1*(1.0-xi3) + phi[1]*(1.0-xi1)*(1.0-xi3) + phi[2]*(1.0-xi1)*xi3 + phi[3]*xi1*xi3);
      return phi[0]*xi1*(1.0-xi3) + phi[1]*(1.0-xi1)*(1.0-xi3) + phi[2]*(1.0-xi1)*xi3 + phi[3]*xi1*xi3;
    case 1: // This triangle is ACD
//      fprintf(stderr,"Now in getPhi! input:(%d,%e,%e), Quad %d, phi = %e\n", trId+1, xi1, xi2, tria2quad[trId][0]+1, phi[0]*xi1*(1.0-xi2) + phi[1]*xi1*xi2 + phi[2]*(1.0-xi1)*xi2 + phi[3]*(1.0-xi1)*(1.0-xi2));
      return phi[0]*xi1*(1.0-xi2) + phi[1]*xi1*xi2 + phi[2]*(1.0-xi1)*xi2 + phi[3]*(1.0-xi1)*(1.0-xi2);
    default:
      fprintf(stderr,"Software bug in the cracking surface...\n");
      exit(-1);
  }
}

//----------------------------------------------------------------------------

void CrackingSurface::printInfo(char* filename)
{
  FILE *myout = fopen(filename,"w");

  fprintf(myout, "...Information about the cracking surface...\n");
  fprintf(myout, "  elemType = %d, n{Total,Used}Quads = %d/%d, n{Total,Used}Trias = %d/%d, n{Total,Used}Nodes = %d/%d.\n",
                  elemType, nTotalQuads, nUsedQuads, nTotalTrias, nUsedTrias, nTotalNodes, nUsedNodes); 

  fprintf(myout, "phantoms(%d)...\n", phantoms.size());
  for(map<int,PhantomElement*>::iterator it=phantoms.begin(); it!=phantoms.end(); it++) {
    fprintf(myout, "  %d --> nodes(%d): ", it->first+1, it->second->nNodes);
    int* nod = it->second->nodes;
    for(int i=0; i<it->second->nNodes; i++)
      fprintf(myout, "%d ", nod[i]+1);
    fprintf(myout, "\n");
    fprintf(myout, "         phi: ");
    double* ph = it->second->phi;
    for(int i=0; i<it->second->nNodes; i++)
      fprintf(myout, "%e ", ph[i]);
    fprintf(myout, "\n");
  }

  fprintf(myout, "latest...\n");
  fprintf(myout, "  phantomQuads(%d): ", latest.phantomQuads.size());
  for(set<int>::iterator it=latest.phantomQuads.begin(); it!=latest.phantomQuads.end(); it++)
    fprintf(myout,"%d ", (*it)+1);
  fprintf(myout, "\n");
  fprintf(myout, "  phantomNodes(%d): ", latest.phantomNodes.size());
  for(map<int,int>::iterator it=latest.phantomNodes.begin(); it!=latest.phantomNodes.end(); it++)
    fprintf(myout,"%d(%d)  ", it->first+1, it->second+1);
  fprintf(myout, "\n");

  fprintf(myout, "tria2quad(%d)...\n", nTotalTrias);
  for(int i=0; i<nTotalTrias; i++)
    fprintf(myout, "  Tria %d <-> Quad %d.%d\n", i+1, tria2quad[i][0]+1, tria2quad[i][1]+1);

  fprintf(myout, "quad2tria(%d)...\n", nTotalQuads);
  for(int i=0; i<nTotalQuads; i++)
    fprintf(myout,"  Quad %d <-> Tria %d and %d\n", i+1, quad2tria[i][0]+1, quad2tria[i][1]+1);

  fprintf(myout, "cracked(%d)...\n", nTotalQuads);
  for(int i=0; i<nTotalQuads; i++)
    fprintf(myout,"  Quad %d -> %d\n", i+1, cracked[i]);

  fclose(myout);
}

//----------------------------------------------------------------------------
