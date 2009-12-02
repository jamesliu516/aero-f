#include "SparseGrid.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//------------------------------------------------------------------------------

template<typename T>
void SparseGrid::tabulate(void (T::*fn)(double *, double *, double *), 
                          T &object, bool restart){

    tabulate(Functor<T>(fn,object), restart);

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::tabulate(FnType fn, bool restart){

  messages(1);

  if(!restart)
    initialize(fn);

  bool success = false;
  bool adaptivity = false;
  int currentMultiIndex, addedSubGrids;

  while(nPoints<maxPoints && !success){

    messages(2);
    
    // get next active subgrid
    currentMultiIndex = currentSubGrid(adaptivity);
    if(currentMultiIndex<0){
      fprintf(stdout, "SparseGrid: no more subgrids available\n");
      exit(1);
    }
    active[currentMultiIndex] = false;
    messages(3,currentMultiIndex);

    // resize arrays if necessary
    while(nSubGrids+dim>sizeMultiIndex) resizeMultiIndex();

    addedSubGrids = 0;
    // for each forward neighbour(one in each direction/dimension),
    //   check if it is admissible
    for(int idim=0; idim<dim; idim++){
      if(admissible(currentMultiIndex, idim)){
        messages(4,currentMultiIndex);
        // find its neighbours and index it as a neighbour of its own neighbours
        findNeighbours(currentMultiIndex, idim, addedSubGrids);
        // the forward admissible neighbour is now active
        active[neighbour[currentMultiIndex][idim]] = true;
        addedSubGrids++;
      }
    }
    messages(5,currentMultiIndex);

    // get the hierarchical surplus of the points on the admissible 
    // forward neighbouring subgrids (which become active subgrids)
    for(int forwardAdmissible=0; forwardAdmissible<addedSubGrids; forwardAdmissible++){
      int numNewPoints = integrateForwardAdmissible(nSubGrids, fn);
      nPoints += numNewPoints;
      nSubGrids++;

      if(adaptivity) nAdaptivePoints += numNewPoints;
      messages(6,forwardAdmissible);
    }

    if(nPoints>minPoints)
      success = checkAccuracy();

  }
  messages(10);

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::initialize(FnType fn){
// initialization of the data structures for the first subgrid
// in the Clenshaw-Curtis grid, which contains only one point located
// at the center of the hypercube (this is level of refinement 0 in all directions).

  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 

  double firstPoint[dim]; double scaledCoord[dim]; double res[out];
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5;
  scale(firstPoint,scaledCoord, 0);
  fn(scaledCoord,res,parameters);
  for(int iout=0; iout<out; iout++){
    surplus[0][iout] = res[iout];
    fnmin[iout]      = res[iout];
    fnmax[iout]      = res[iout];
  }
  
  activeError[0] = surplus[0][0]; //max(surplus[0]);
  for(int iout=1; iout<out; iout++)
    if(activeError[0]<surplus[0][iout]) activeError[0] = surplus[0][iout];
  activeCost[0]  = 1.0;
  activeHeapError.insert(0,activeError);
  activeHeapCost.insert(0,activeCost);
  active[0] = true;

  for(int iout=0; iout<out; iout++)
    error[0][iout] = surplus[0][iout];

  // no backward neighbour for subGrid multiIndex[0]
  for(int isize=0; isize<sizeMultiIndex; isize++)
    for(int neigh=0; neigh<2*dim; neigh++)
      neighbour[isize][neigh] = -1;
}

//------------------------------------------------------------------------------

template<typename FnType>
int SparseGrid::integrateForwardAdmissible(const int newSubGridIndex,
                            FnType fn){

  int nPointsSubGrid;
  double ** subGrid = generateSubGrid(newSubGridIndex,nPointsSubGrid);
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  evaluatePreviousInterpolation(subGrid, nPointsSubGrid);

  updateError(nPointsSubGrid);
  updateCost();
  activeHeapError.insert(newSubGridIndex,activeError);
  activeHeapCost.insert(newSubGridIndex,activeCost);

  for (int i=0; i<nPointsSubGrid; i++) delete [] subGrid[i];
  delete [] subGrid;
  return nPointsSubGrid;

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::evaluateFunctionOnGrid(double **subGrid,
                                        const int nPointsSubGrid, FnType fn){

  while(nPoints+nPointsSubGrid > sizeSurplus) resizeSurplus();

  double scaledCoord[dim];
  double res[out];
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    scale(subGrid[iPts],scaledCoord, 0);
    fn(scaledCoord,res,parameters);
    for(int iout=0; iout<out; iout++)
      surplus[nPoints+iPts][iout] = res[iout];
  }
  
  bounds(nPointsSubGrid); //find the min and max values of the target function
                          //on the grid points

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::test(FnType fn){



  // to test evaluatePreviousInterpolation and evaluateFunctionOnGrid
  //initialize
/*  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 
  Coord firstPoint; double res;
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5*(range[idim][0]+range[idim][1]);
  evaluateFunction(firstPoint[0],firstPoint[1],0.0,res,fn);
  for(int iout=0; iout<out; iout++) surplus[0][iout] = res;
  
  for(int iPts=0; iPts<nPoints; iPts++)
    for(int iout=0; iout<out; iout++)
      fprintf(stderr, "## initialization -- surplus[%d] = %e\n", iPts, surplus[iPts][iout]);

  //new subgrid to consider
  //multiIndex[1][0]=0; multiIndex[1][1]=0;
  multiIndex[1][0]=1; multiIndex[1][1]=1;
  multiIndex[2][0]=1; multiIndex[2][1]=2;

  int nPointsSubGrid;
  Coord * subGrid = generateSubGrid(1,nPointsSubGrid);
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  for(int iPts=nPoints; iPts<nPoints+nPointsSubGrid; iPts++)
    for(int iout=0; iout<out; iout++)
      fprintf(stderr, "## after function evaluation -- surplus[%d] = %e\n", iPts, surplus[iPts][iout]);
  evaluatePreviousInterpolation(subGrid, nPointsSubGrid);

  for(int iPts=nPoints; iPts<nPoints+nPointsSubGrid; iPts++)
    for(int iout=0; iout<out; iout++)
      fprintf(stderr, "## after interpolation -- surplus[%d] = %e\n", iPts, surplus[iPts][iout]);

  nPoints += nPointsSubGrid;
  nSubGrids++;

  subGrid = generateSubGrid(2,nPointsSubGrid);
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  for(int iPts=nPoints; iPts<nPoints+nPointsSubGrid; iPts++)
    for(int iout=0; iout<out; iout++)
      fprintf(stderr, "## after function evaluation -- surplus[%d] = %e\n", iPts, surplus[iPts][iout]);
  evaluatePreviousInterpolation(subGrid, nPointsSubGrid);

  for(int iPts=nPoints; iPts<nPoints+nPointsSubGrid; iPts++)
    for(int iout=0; iout<out; iout++)
      fprintf(stderr, "## after interpolation -- surplus[%d] = %e\n", iPts, surplus[iPts][iout]);
*/

  // to test generateSubGrid and tensorize
/*
  assert(dim==4 && out==1);

  multiIndex[0][0]=1; multiIndex[0][1]=1; multiIndex[0][2] = 3; multiIndex[0][3] = 2;
  multiIndex[1][0]=1; multiIndex[1][1]=2; multiIndex[1][2] = 0; multiIndex[1][3] = 2;
  nSubGrids = 2;
  int nPointsSubGrid;
  Coord * subGrid0 = generateSubGrid(0,nPointsSubGrid);
  fprintf(stderr, "number of points in subgrid0 is %d\n", nPointsSubGrid);
  for(int nPts = 0; nPts<nPointsSubGrid; nPts++)
    fprintf(stderr, "coord[%d] = (%e %e %e %e)\n",nPts, subGrid0[nPts][0], subGrid0[nPts][1], subGrid0[nPts][2], subGrid0[nPts][3]);
  Coord * subGrid1 = generateSubGrid(1,nPointsSubGrid);
  fprintf(stderr, "number of points in subgrid1 is %d\n", nPointsSubGrid);
  for(int nPts = 0; nPts<nPointsSubGrid; nPts++)
    fprintf(stderr, "coord[%d] = (%e %e %e)\n",nPts, subGrid1[nPts][0], subGrid1[nPts][1], subGrid1[nPts][2], subGrid1[nPts][3]);
*/


  // to test heap class
  /*nSubGrids = 8;
  multiIndex[0][0]=0; multiIndex[0][1]=0;
  multiIndex[1][0]=1; multiIndex[1][1]=0;
  multiIndex[2][0]=0; multiIndex[2][1]=1;
  multiIndex[3][0]=0; multiIndex[3][1]=2;
  multiIndex[4][0]=1; multiIndex[4][1]=1;
  multiIndex[5][0]=0; multiIndex[5][1]=3;
  multiIndex[6][0]=2; multiIndex[6][1]=0;
  multiIndex[7][0]=7; multiIndex[7][1]=7;
  for(int toto=0; toto<8; toto++)
    fprintf(stderr, "multiIndex[%d]= %d, %d\n", toto, multiIndex[toto][0], multiIndex[toto][1]);

  activeError[0] = 0.1;
  activeError[1] = 0.2;
  activeError[2] = 0.3;
  activeError[3] = 0.4;
  activeError[4] = 0.5;
  activeError[5] = 0.6;
  activeError[6] = 0.7;
  activeError[7] = 0.8;
  activeError[8] = 0.9;

  activeHeapError.insert(0,activeError);
  activeHeapError.insert(1,activeError);
  activeHeapError.insert(2,activeError);
  activeHeapError.insert(3,activeError);
  activeHeapError.insert(4,activeError);
  activeHeapError.insert(5,activeError);
  activeHeapError.insert(6,activeError);
  activeHeapError.insert(7,activeError);

  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  
  int popped = activeHeapError.pop(activeError);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  popped = activeHeapError.pop(activeError);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  popped = activeHeapError.pop(activeError);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);

  activeHeapError.insert(1,activeError);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  */
}
