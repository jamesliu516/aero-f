#include <PolygonReconstructionData.h>

#include <Elem.h>
#include "LevelSet/LevelSetStructure.h"
#include <Vector3D.h>

#include <cstdlib>
#include <cstdio>
//------------------------------------------------------------------------------

int getPolygons(Elem &elem, LevelSetStructure &LSS, PolygonReconstructionData* polygons)
{
    int numberOfPolygons=0;
    int oppositeNodes[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    int edgeToNodes[6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
    int edgeToOppositeNodes[6][2] = {{2,3}, {0,3}, {3,1}, {1,2}, {0,2}, {0,1}};
    int T[4]; for (int i=0; i<4; ++i) T[i] = elem[i]; //nodes in a tet.
    bool isBlocked[4][4],exit_early(true);
    for(int i=0; i<6; ++i){int ni=edgeToNodes[i][0],nj=edgeToNodes[i][1];
        isBlocked[ni][nj] = isBlocked[nj][ni] = LSS.edgeIntersectsStructure(0,T[ni],T[nj]);
        if(isBlocked[ni][nj]) exit_early=false;}

    if(exit_early) return 0; // No intersections detected; don't bother finding a reconstructed surface that does not exist.

    bool finished[4]; for(int i=0; i<4; ++i) finished[i]=false;

    // First, identify any of the simple cases
    for(int current_node=0; current_node<4; ++current_node){
        int n0=oppositeNodes[current_node][0],n1=oppositeNodes[current_node][1],n2=oppositeNodes[current_node][2];
        if(isBlocked[current_node][n0] && isBlocked[current_node][n1] && isBlocked[current_node][n2]){
            finished[current_node]=true;
            polygons[numberOfPolygons++].AssignTriangle(T[current_node],T[n0],T[n1],T[n2]);
            if(!isBlocked[n0][n1] && !isBlocked[n0][n2] && !isBlocked[n1][n2]){
                polygons[numberOfPolygons++].AssignTriangle(T[current_node],T[n0],T[n1],T[n2],false);
                return numberOfPolygons;}}} // REALLY simple case...
    
    for(int current_node=0; current_node<4; ++current_node){if(finished[current_node]) continue;
        for(int n=0; n<3; ++n){int neighbor_node=oppositeNodes[current_node][n];
            int e0=oppositeNodes[current_node][(n+1)%3],e1=oppositeNodes[current_node][(n+2)%3];
            if(!isBlocked[current_node][neighbor_node] && isBlocked[current_node][e0] && isBlocked[current_node][e1]){
                if(isBlocked[neighbor_node][e0] && isBlocked[neighbor_node][e1]){
                    finished[current_node] = finished[neighbor_node] = true;
                    polygons[numberOfPolygons++].AssignQuadrilateral(T[current_node],T[neighbor_node],T[e0],T[e1]);}
                else if(isBlocked[neighbor_node][e0]){
                    finished[current_node] = finished[neighbor_node] = true;
                    polygons[numberOfPolygons++].AssignQuadTriangle(T[current_node],T[neighbor_node],T[e0],T[e1]);}
                else if(isBlocked[neighbor_node][e1]){
                    finished[current_node] = finished[neighbor_node] = true;
                    polygons[numberOfPolygons++].AssignQuadTriangle(T[current_node],T[neighbor_node],T[e1],T[e0]);}
                else {
                    finished[current_node] = true;
                    polygons[numberOfPolygons++].AssignTwoEdges(T[current_node],T[e1],T[e0]);}}}}

    // Finally, check for single-intersection edge cases
    for(int current_node=0; current_node<4; ++current_node) if(!finished[current_node]){
        for(int i=0; i<3; ++i) if(isBlocked[current_node][oppositeNodes[current_node][i]]){
            polygons[numberOfPolygons++].AssignSingleEdge(T[current_node],T[oppositeNodes[current_node][i]]);}}

    assert(numberOfPolygons);
    return numberOfPolygons;
}

//--------------------------------------------------------------------------

int getPolygon(Elem &elem, LevelSetStructure &LSS, int polygon[4][2])
{
  int facet[4][3] = {{0,1,3}, {1,2,3}, {0,3,2}, {0,2,1}};
  int edgeIndex[4][3] = {{0,4,3}, {1,5,4}, {3,5,2}, {2,1,0}};
  int edgeToNodes[6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
  int T[4]; //nodes in a tet.

  for (int i=0; i<4; i++)
    T[i] = elem[i];
  for (int i=0; i<4; i++)
    polygon[i][0] = polygon[i][1] = -1;

  int nextEdge[6];
  for (int i=0; i<6; i++) nextEdge[i] = -1;
  int firstEdge = -1;

  int outside_status=-5;
  for (int iFacet=0; iFacet<4; iFacet++) for (int j=0; j<3; j++) if(LSS.fluidModel(0,T[facet[iFacet][j]])==0) outside_status=0;
  if(outside_status != 0) for (int iFacet=0; iFacet<4; iFacet++) for (int j=0; j<3; j++) outside_status=std::max(outside_status,LSS.fluidModel(0,T[facet[iFacet][j]]));
  for (int iFacet=0; iFacet<4; iFacet++) {
    int status = 0;
    for (int j=0; j<3; j++)
      if(LSS.fluidModel(0,T[facet[iFacet][j]])==outside_status) // TODO(jontg): This DOES NOT WORK for thin-shells.
         status |= 1 << j;  //KW: set the j-th (counting from 0) binary digit to 1
    switch (status) {
      case 1:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][0];
        nextEdge[edgeIndex[iFacet][0]] = edgeIndex[iFacet][2];
        break;
      case 2:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][1];
        nextEdge[edgeIndex[iFacet][1]] = edgeIndex[iFacet][0];
        break;
      case 3:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][1];
        nextEdge[edgeIndex[iFacet][1]] = edgeIndex[iFacet][2];
        break;
      case 4:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][2];
        nextEdge[edgeIndex[iFacet][2]] = edgeIndex[iFacet][1];
        break;
      case 5:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][0];
        nextEdge[edgeIndex[iFacet][0]] = edgeIndex[iFacet][1];
        break;
      case 6:
        if (firstEdge<0)
          firstEdge = edgeIndex[iFacet][2];
        nextEdge[edgeIndex[iFacet][2]] = edgeIndex[iFacet][0];
        break;
      default:
        break;
    }
  }

  if (firstEdge<0) return 0;
  int curEdge = firstEdge;
  int edgeCount = 0;
  do {
    if (edgeCount>4) {
      fprintf(stderr,"Error in getPolygon. edgeCount = %d > 4, Abort...\n", edgeCount);
      exit(-1);
    }
    if (LSS.fluidModel(0,T[edgeToNodes[curEdge][0]])==outside_status) {
      polygon[edgeCount][0] = T[edgeToNodes[curEdge][0]];
      polygon[edgeCount][1] = T[edgeToNodes[curEdge][1]];
    }else {
      polygon[edgeCount][0] = T[edgeToNodes[curEdge][1]];
      polygon[edgeCount][1] = T[edgeToNodes[curEdge][0]];
    }
    edgeCount++;
    curEdge = nextEdge[curEdge];
  } while (curEdge>=0 && curEdge!=firstEdge);

  return edgeCount;
}
//--------------------------------------------------------------------------

void getPolygonNormal(SVec<double,3>& X, Vec3D &normal, LevelSetStructure &LSS, PolygonReconstructionData &polygon)
{
  int nEdges = polygon.numberOfEdges;
  LevelSetResult lsRes[nEdges];
  Vec3D Xinter[nEdges];
  Vec3D tempoLoc;
  for(int m=0; m<nEdges; m++){
    lsRes[m] = LSS.getLevelSetDataAtEdgeCenter(0,polygon.edgeWithVertex[m][0],polygon.edgeWithVertex[m][1]);
    double alpha = lsRes[m].alpha;
    if (alpha<0) {fprintf(stderr,"Unable to get intersection results at edge center! Abort...\n"); exit(-1);}
    Xinter[m][0] = alpha*X[polygon.edgeWithVertex[m][0]][0] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][0];
    Xinter[m][1] = alpha*X[polygon.edgeWithVertex[m][0]][1] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][1];
    Xinter[m][2] = alpha*X[polygon.edgeWithVertex[m][0]][2] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][2];
  }
  switch(nEdges){
    case 3: // got a triangle.
      normal = (Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
      if (normal.norm() != 0.0) {normal = 1.0/normal.norm()*normal;}
      // Then we check the orientation of the normal.
      for (int i=0; i<3; i++) {tempoLoc[i] = X[polygon.edgeWithVertex[0][0]][i];}
      if (LSS.isActive(0,polygon.edgeWithVertex[0][0])) {
        if (normal*(tempoLoc-Xinter[0]) < 0) {normal = -1.0*normal;}
      } else if (normal*(tempoLoc-Xinter[0]) > 0) {normal = -1.0*normal;}
      break;
    case 4: // got a quadrangle... We cut it into two triangles and we average the two resulting normals. 
      Vec3D normal1, normal2;
      normal1 = (Xinter[1]-Xinter[0])^(Xinter[3]-Xinter[0]);
      if (normal1.norm() != 0.0) {normal1 = 1.0/normal1.norm()*normal1;}
      normal2 = (Xinter[2]-Xinter[1])^(Xinter[3]-Xinter[1]);
      if (normal2.norm() != 0.0) {normal2 = 1.0/normal2.norm()*normal2;}
      normal = 0.5*(normal1+normal2);
      if (normal.norm() != 0.0) {normal = 1.0/normal.norm()*normal;}
      // Then we check the orientation of the normal.
      for (int i=0; i<3; i++) {tempoLoc[i] = X[polygon.edgeWithVertex[0][0]][i];}
      if (LSS.isActive(0,polygon.edgeWithVertex[0][0])) {
        if (normal*(tempoLoc-Xinter[0]) < 0) {normal = -1.0*normal;}
      } else if (normal*(tempoLoc-Xinter[0]) > 0) {normal = -1.0*normal;}
      break;
  }
}

