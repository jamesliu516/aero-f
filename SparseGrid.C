#include <SparseGrid.h>



//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::SparseGrid(){

  nPoints     = 0;
  sizeSurplus = 0;
  surplus     = 0;

  nSubGrids      = 0;
  sizeMultiIndex = 0;
  multiIndex     = 0;

  maxPoints = 100;
  minPoints = 100;
  absAccuracy = 1.0;
  relAccuracy = 0.1;
  for(int idim=0; idim<dim; idim++){
    range[idim][0] = 0.0;
    range[idim][1] = 0.0;
  }
  dimAdaptDegree = 0.0;

  nAdaptivePoints = 0;
  nInactiveSet = 0;
  inactiveSet = 0;
  activeSetAll = 0;

  errorActiveSet = 0;
  costActiveSet  = 0;

  neighbour = 0;
  error = 0;

}

//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::~SparseGrid(){

  delete [] surplus;
  delete [] multiIndex;
  delete [] inactiveSet;
  delete [] activeSetAll;
  delete [] errorActiveSet;
  delete [] costActiveSet;
  delete [] neighbour;
  delete [] error;

}

//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::SparseGrid(int maxPts_, int minPts_, double absAcc_,
          double relAcc_, double range_[dim][2], double dimAdaptDegree_){

  maxPoints = maxPts_;
  minPoints = minPts_;
  absAccuracy = absAcc_;
  relAccuracy = relAcc_;
  for(int idim=0; idim<dim; idim++){
    range[idim][0] = range_[idim][0];
    range[idim][1] = range_[idim][1];
  }
  dimAdaptDegree = dimAdaptDegree_; //min(max(dimAdaptDegree_,0),1);

  nPoints         = 0;
  sizeSurplus     = maxPoints;
  surplus         = new Output[sizeSurplus];

  nSubGrids       = 0;
  sizeMultiIndex  = 10;
  multiIndex      = new MultiIndex[sizeMultiIndex];

  nAdaptivePoints = 0;
  nInactiveSet    = 0;
  inactiveSet     = new int[sizeMultiIndex];
  activeSetAll    = new bool[sizeMultiIndex];

  errorActiveSet  = new double[sizeMultiIndex];
  costActiveSet   = new double[sizeMultiIndex];

  neighbour       = new Neighbour[sizeMultiIndex];
  error           = new Output[sizeMultiIndex];
  
}

//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::SparseGrid(SparseGridData &data){

  maxPoints = data.maxPoints;
  minPoints = data.minPoints;
  absAccuracy = data.absAccuracy;
  relAccuracy = data.relAccuracy;
  range[0][0] = data.range1min;
  range[0][1] = data.range1max;
  range[1][0] = data.range2min;
  range[1][1] = data.range2max;
  range[2][0] = data.range3min;
  range[2][1] = data.range3max;
  dimAdaptDegree = data.dimAdaptDegree; //min(max(dimAdaptDegree_,0),1);

  nPoints         = 0;
  sizeSurplus     = maxPoints;
  surplus         = new Output[sizeSurplus];

  nSubGrids       = 0;
  sizeMultiIndex  = 10;
  multiIndex      = new MultiIndex[sizeMultiIndex];

  nAdaptivePoints = 0;
  nInactiveSet    = 0;
  inactiveSet     = new int[sizeMultiIndex];
  activeSetAll    = new bool[sizeMultiIndex];
  errorActiveSet  = new double[sizeMultiIndex];
  costActiveSet   = new double[sizeMultiIndex];

  neighbour       = new Neighbour[sizeMultiIndex];
  error           = new Output[sizeMultiIndex];
  
}

//------------------------------------------------------------------------------

template<int dim, int out>
void SparseGrid<dim,out>::tabulate(void (*fn)(double, double, double, double &)){

  fprintf(stdout, "###########################################\n");
  fprintf(stdout, "#     INITIALIZATION\n");
  fprintf(stdout, "###########################################\n");
  initialize(fn);

  bool success = false;
  bool adaptivity = false;
  int currentMultiIndex, addedSubGrids;

  while(nPoints<maxPoints && !success){
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "#     BEGIN ITERATION\n");
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "#     STATE AT BEGINNING OF ITERATION\n");
    fprintf(stdout, "# activeSetAll is:");
    for(int kk = 0; kk<nSubGrids; kk++)
      fprintf(stdout, "  %d", activeSetAll[kk]);
    fprintf(stdout, "\n");
    fprintf(stdout, "# multiIndex is:");
    for(int idim = 0; idim<dim; idim++){
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %d", multiIndex[kk][idim]);
      fprintf(stdout, "\n#               ");
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "# errorActiveSet is:");
    for(int kk = 0; kk<nSubGrids; kk++)
      fprintf(stdout, "  %e", errorActiveSet[kk]);
    fprintf(stdout, "\n");
    fprintf(stdout, "# costActiveSet is:");
    for(int kk = 0; kk<nSubGrids; kk++)
      fprintf(stdout, "  %e", costActiveSet[kk]);
    fprintf(stdout, "\n");
    fprintf(stdout, "# errors are:");
    for(int iout = 0; iout<out; iout++){
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %e", error[kk][iout]);
      fprintf(stdout, "\n#               ");
    }
    fprintf(stdout, "\n");
    // get next active subgrid
    currentMultiIndex = currentSubGrid(adaptivity);
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "# adding grid number %d(%d) to inactive grid set\n", currentMultiIndex, adaptivity);
    fprintf(stdout, "###########################################\n");
    // get the backward neighbours of that subgrid (why?)

    // resize arrays if potential number of grids exceeds sizeMultiIndex at end of iteration
    while(nSubGrids+dim>sizeMultiIndex) resizeMultiIndex();

    addedSubGrids = 0;
    // for each forward neighbour(one in each direction),
    //   check if it is admissible
    for(int idim=0; idim<dim; idim++){
      if(admissible(currentMultiIndex, idim)){
        fprintf(stdout, "###########################################\n");
        fprintf(stdout, "# neighbour of %d in direction %d is admissible\n", currentMultiIndex, idim);
        fprintf(stdout, "###########################################\n");
        // find its own neighbours
        findNeighbours(currentMultiIndex, idim, addedSubGrids);
        activeSetAll[neighbour[currentMultiIndex][idim]] = true;
        addedSubGrids++;
      }
    }
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "# neighbours of %d are:", currentMultiIndex);
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      fprintf(stdout, "  %d", neighbour[currentMultiIndex][ineigh]);
    fprintf(stdout, "\n");
    fprintf(stdout, "###########################################\n");
      

    // get the hierarchical surplus of the points on the admissible 
    // forward neighbouring subgrids (which become active subgrids)
    for(int forwardAdmissible=0; forwardAdmissible<addedSubGrids; forwardAdmissible++){
      int numNewPoints = integrateForwardAdmissible(nSubGrids, fn);
      nPoints += numNewPoints;
      nSubGrids++;
      if(adaptivity) nAdaptivePoints += numNewPoints;
      fprintf(stdout, "###########################################\n");
      fprintf(stdout, "# forward admissible grid %d leads to %d points and %d subgrids\n", forwardAdmissible, nPoints, nSubGrids);
      fprintf(stdout, "###########################################\n");
    }
    activeSetAll[currentMultiIndex] = false;

    if(nPoints>minPoints)
      success = checkAccuracy();

    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "###########################################\n");
  fprintf(stdout, "#     SPARSE GRID HAS BEEN CREATED\n");
  fprintf(stdout, "###########################################\n");

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::initialize(void (*fn)(double,double,double,double &)){
// initialization of the data structures for the first subgrid
// in the Clenshaw-Curtis grid, which contains only one point located
// at the center of the hypercube (this is level of refinement 0 in all directions).
  double temp = 0.0;

  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 
  Coord firstPoint; double res;
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5*(range[idim][0]+range[idim][1]);
  if(out==1) fn(firstPoint[0],firstPoint[1],temp,res);
  for(int iout=0; iout<out; iout++) surplus[0][iout] = res;
  
  
  errorActiveSet[0] = surplus[0][0]; //max(surplus[0]);
  for(int iout=1; iout<out; iout++)
    if(errorActiveSet[0]<surplus[0][iout]) errorActiveSet[0] = surplus[0][iout];
  costActiveSet[0]  = 1.0;
  activeHeapError.insert(0,errorActiveSet);
  activeHeapCost.insert(0,costActiveSet);
  activeSetAll[0] = true;

  for(int iout=0; iout<out; iout++)
    error[0][iout] = surplus[0][iout];
  // no backward neighbour for subGrid multiIndex[0]
  for(int isize=0; isize<sizeMultiIndex; isize++)
    for(int neigh=0; neigh<2*dim; neigh++)
      neighbour[isize][neigh] = -1;
}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
int SparseGrid<dim,out>::currentSubGrid(bool &adaptivity){

  bool done = false;
  int current  = -1;

  adaptivity = false;

  while(!done){
    current = -1;
    if(nAdaptivePoints>dimAdaptDegree*nPoints) // non-adaptive
      current = activeHeapCost.pop(costActiveSet);

    if(current<0){ // adaptive
      current = activeHeapError.pop(errorActiveSet);
      if(current>=0) adaptivity = true;
    }

    if(current<0)
      current = activeHeapCost.pop(costActiveSet);

    if(current<0) return current;

    // make usre it is really active
    if(activeSetAll[current]) done = true;

  }

  return current;

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
bool SparseGrid<dim,out>::admissible(const int currentMultiIndex,
                                     const int forwardDir){

  // currentMultiIndex has multi-index multiIndex[currentMultiIndex]
  // the forward neighbour considered here has multi-index:
  //    multiIndex[currentMultiIndex] + canonical_vector(forwardDir)
  // To check its admissibility, all its backward neighbours must be inactive
  //   subgrids. In particular, its backward neighbour currentMultiIndex
  //   refers to an inactive subgrid. 

// NOT VALID FOR FIRST POINT of FIRST GRID!
  int backwardOfCurrent, backwardOfForward;
  for(int idim=0; idim<dim; idim++){
    if(forwardDir == idim) continue;
    backwardOfCurrent = neighbour[currentMultiIndex][dim+idim];
    if(backwardOfCurrent<0){ // already at lowest level in direction idim
      if(multiIndex[currentMultiIndex][idim]!=0){
        fprintf(stdout, "*** Error: contradiction detected (no backward neighbours but not at lowest level...\n");
        exit(1);
      }
      continue;
    }
    backwardOfForward = neighbour[backwardOfCurrent][forwardDir];
    if(backwardOfForward<0) return false; // not part of sparse grid 
    if(activeSetAll[backwardOfForward]) return false; // part of sparse grid, but active
    
  }

  return true;

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::findNeighbours(const int currentMultiIndex,
                                         const int forwardDir, const int addedSubGrids){
// NOT VALID FOR FIRST POINT of FIRST GRID!
  int backwardOfCurrent, backwardOfForward;

  for(int idim=0; idim<dim; idim++){
    if(forwardDir == idim){
      neighbour[currentMultiIndex      ][forwardDir    ] = nSubGrids+addedSubGrids;
      neighbour[nSubGrids+addedSubGrids][dim+forwardDir] = currentMultiIndex;
    }else{
      backwardOfCurrent = neighbour[currentMultiIndex][dim+idim];
      if(!(backwardOfCurrent<0)){ 
        backwardOfForward = neighbour[backwardOfCurrent][forwardDir];
        neighbour[backwardOfForward      ][idim    ] = nSubGrids+addedSubGrids;
        neighbour[nSubGrids+addedSubGrids][dim+idim] = backwardOfForward;
      }
    }
  }

  for(int idim=0; idim<dim; idim++){
//    fprintf(stderr, "### findNeighbours -- %d is a backward neighbour of %d in direction %d\n", neighbour[nSubGrids+addedSubGrids][idim+dim], nSubGrids+addedSubGrids, idim);
//    fprintf(stderr, "### findNeighbours -- it has level %d\n", multiIndex[neighbour[nSubGrids+addedSubGrids][idim+dim]][idim]);
    if(neighbour[nSubGrids+addedSubGrids][idim+dim]>=0)
      multiIndex[nSubGrids+addedSubGrids][idim] = multiIndex[neighbour[nSubGrids+addedSubGrids][idim+dim]][idim]+1;
    else multiIndex[nSubGrids+addedSubGrids][idim] = 0;
  }

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
int SparseGrid<dim,out>::integrateForwardAdmissible(const int newSubGridIndex,
                            void (*fn)(double,double,double,double &)){

  int nPointsSubGrid;
  //fprintf(stderr, "### integrateForwardAdmissible -- generating SubGrid\n");
  Coord * subGrid = generateSubGrid(newSubGridIndex,nPointsSubGrid);
  //fprintf(stderr, "### integrateForwardAdmissible -- evaluating Function on SubGrid\n");
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  //fprintf(stderr, "### integrateForwardAdmissible -- evaluating Interpolation on Grid\n");
  evaluatePreviousInterpolation(subGrid, nPointsSubGrid);
  //fprintf(stderr, "### integrateForwardAdmissible -- updating Errors\n");
  updateError(nPointsSubGrid);
  //fprintf(stderr, "### integrateForwardAdmissible -- updating Costs\n");
  updateCost();
  //fprintf(stderr, "### integrateForwardAdmissible -- inserting new grid in active heap for error\n");
  activeHeapError.insert(newSubGridIndex,errorActiveSet);
  //fprintf(stderr, "### integrateForwardAdmissible -- inserting new grid in active heap for cost\n");
  activeHeapCost.insert(newSubGridIndex,costActiveSet);
  //fprintf(stderr, "### integrateForwardAdmissible -- about to return\n");

  delete [] subGrid;
  return nPointsSubGrid;

}

//------------------------------------------------------------------------------

template<int dim, int out>
typename SparseGrid<dim,out>::Coord *
SparseGrid<dim,out>::generateSubGrid(const int newSubGrid,
                                     int &nPointsSubGrid){

  // find number of points in subgrid (Clenshaw-Curtis type)
  //   and coordinates of the points in each direction
  // tensorization is done later.
  fprintf(stderr, "the new subgrid (%d) has multi indices:\n", newSubGrid);
  for(int idim=0; idim<dim; idim++)
    fprintf(stderr, "%d  ", multiIndex[newSubGrid][idim]);
  nPointsSubGrid = 1;
  int nPointsDim[dim];
  int nCumulatedPointsDim[dim];
  double **coordDim = new double*[dim];
  for(int i=0; i<dim; i++){
    if(multiIndex[newSubGrid][i] == 0){
      nPointsDim[i] = 1;
      coordDim[i] = new double[1];
      coordDim[i][0] = 0.5;
    }
    else if(multiIndex[newSubGrid][i] < 3){
      nPointsDim[i] = 2;
      coordDim[i] = new double[2];
      if(multiIndex[newSubGrid][i] == 1){
        coordDim[i][0] = 0.0;
        coordDim[i][1] = 1.0;
      }else{
        coordDim[i][0] = 0.25;
        coordDim[i][1] = 0.75;
      }
    }
    else{
      nPointsDim[i] = static_cast<int>(pow(2,multiIndex[newSubGrid][i]-1));
      coordDim[i] = new double[nPointsDim[i]];
      for(int j=0; j<nPointsDim[i]; j++)
        coordDim[i][j] = (2.0*j+1.0)/(2.0*nPointsDim[i]);
    }

    nCumulatedPointsDim[i] = nPointsDim[i];
    nPointsSubGrid *= nPointsDim[i];

    if(i>0) nCumulatedPointsDim[i] *= nCumulatedPointsDim[i-1];
    if(i==dim-1 && nCumulatedPointsDim[i]!=nPointsSubGrid)
      fprintf(stdout, "*** Error: the number of points in this subgrid is not correct\n");
  }

  // tensorization of the coordinates is returned in res
  Coord *res = new Coord[nPointsSubGrid];
  tensorize(res,coordDim,nPointsDim,nCumulatedPointsDim);

  for(int i=0; i<dim; i++) delete [] coordDim[i];
  delete [] coordDim;

  return res;

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::tensorize(Coord *res, double ** coordDim,
                                    int *nPtsDim, int *nCumulatedPtsDim){

  for(int newdim=0; newdim<dim; newdim++){

    if(newdim==0){
      for(int iPts=0; iPts<nPtsDim[newdim]; iPts++)
        res[iPts][newdim] = coordDim[newdim][iPts];

    }else{
      for(int iPts=0; iPts<nPtsDim[newdim]; iPts++)
        for(int iPtsRepeat=0; iPtsRepeat<nCumulatedPtsDim[newdim-1]; iPtsRepeat++)
          res[nCumulatedPtsDim[newdim-1]*iPts+iPtsRepeat][newdim] = coordDim[newdim][iPts];
  
      for(int iPts=0; iPts<nCumulatedPtsDim[newdim-1]; iPts++)
        for(int iPtsRepeat=0; iPtsRepeat<nPtsDim[newdim]; iPtsRepeat++)
          for(int otherdim=newdim-1; otherdim>=0; otherdim--)
            res[iPts+iPtsRepeat*nCumulatedPtsDim[newdim-1]][otherdim] = res[iPts][otherdim];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::evaluateFunctionOnGrid(Coord *subGrid,
                                                 const int nPointsSubGrid,
                                                 void (*fn)(double,double,double,double &)){

  while(nPoints+nPointsSubGrid > sizeSurplus) resizeSurplus();

  Coord scaledCoord;
  double res, temp = 1.0;
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    scale(subGrid[iPts],scaledCoord, 0);
    if(dim==2 && out==1) fn(scaledCoord[0],scaledCoord[1],temp,res);
    //else if(dim==3) res = fn(scaledCoord[0],scaledCoord[1],scaledCoord[2]);
    else{
      fprintf(stdout, "*** Error: (Sparse Grid) wrong function call\n");
      exit(1);
    }
    surplus[nPoints+iPts][0] = res;
  }

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::resizeSurplus(){

  Output *larger = new Output[2*sizeSurplus];
  for(int i=0; i<sizeSurplus; i++)
    for(int iout=0; iout<out; iout++)
      larger[i][iout] = surplus[i][iout];

  sizeSurplus *= 2;

  delete [] surplus;
  surplus = larger;

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::resizeMultiIndex(){

  MultiIndex *largerMultiIndex = new MultiIndex[2*sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int idim=0; idim<dim; idim++)
      largerMultiIndex[i][idim] = multiIndex[i][idim];

  delete [] multiIndex;
  multiIndex = largerMultiIndex;

  double *largerErrorActiveSet = new double[2*sizeMultiIndex];
  double *largerCostActiveSet = new double[2*sizeMultiIndex];
  bool *largerActiveSetAll = new bool[2*sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++){
    largerErrorActiveSet[i]= errorActiveSet[i];
    largerCostActiveSet[i]= costActiveSet[i];
    largerActiveSetAll[i]= activeSetAll[i];
  }

  delete [] errorActiveSet; delete [] costActiveSet; delete [] activeSetAll;
  errorActiveSet = largerErrorActiveSet;
  costActiveSet = largerCostActiveSet;
  activeSetAll = largerActiveSetAll;

  Output *largerError = new Output[2*sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int iout=0; iout<out; iout++)
      largerError[i][iout] = error[i][iout];

  delete [] error;
  error = largerError;

  Neighbour *largerNeighbour = new Neighbour[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = -1;
  for(int i=0; i<nSubGrids; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = neighbour[i][ineigh];

  delete [] neighbour;
  neighbour = largerNeighbour;

  sizeMultiIndex *= 2;

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::scale(Coord subGrid, Coord &scaledCoord, int op){

  if(op==0)
    for(int i=0; i<dim; i++)
      scaledCoord[i] = subGrid[i]*(range[i][1]-range[i][0])+range[i][0];
  else if(op==1)
    for(int i=0; i<dim; i++)
      scaledCoord[i] = (subGrid[i]-range[i][0])/(range[i][1]-range[i][0]);

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::evaluatePreviousInterpolation(Coord *subGrid,
                                                 const int nPointsSubGrid){

  Output temp;
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    singleInterpolation(subGrid[iPts], temp);
    for(int iout=0; iout<out; iout++) 
      surplus[nPoints+iPts][iout] -= temp[iout];
  }

  /*for(int iPts=0; iPts<nPointsSubGrid; iPts++)
    for(int iout=0; iout<out; iout++) 
      fprintf(stderr, "### surplus[%d][%d] = %e\n", nPoints+iPts, iout, surplus[nPoints+iPts][iout]); 
*/
}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::updateError(const int nPointsSubGrid){

  Output indicator; double temp;

  // computes the "errors" given by the surpluses for a subgrid in each output
  for(int iout=0; iout<out; iout++) indicator[iout] = 0.0;
  for(int iout=0; iout<out; iout++)
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      indicator[iout] += ((temp>0) ? temp : -temp);
    }

  errorActiveSet[nSubGrids] = 0.0;
  for(int iout=0; iout<out; iout++)
    if(indicator[iout] > errorActiveSet[nSubGrids])
      errorActiveSet[nSubGrids] = indicator[iout];

  //fprintf(stderr, "### updateError -- errorActiveSet[%d] = %e\n", nSubGrids, errorActiveSet[nSubGrids]);

  // computes the "errors" for the 
  for(int iout=0; iout<out; iout++){
    error[nSubGrids][iout] = 0.0;
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      temp = (temp>0) ? temp : -temp;
      if(temp>error[nSubGrids][iout]) 
        error[nSubGrids][iout] = temp;
    }
  }
      
}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::updateCost(){

  costActiveSet[nSubGrids] = 0.0;
  for(int idim=0; idim<dim; idim++){
    costActiveSet[nSubGrids] += static_cast<double>(multiIndex[nSubGrids][idim]);
    //fprintf(stderr, "### updateCost -- costActiveSet[%d] = %e\n", nSubGrids, costActiveSet[nSubGrids]);
  }
  costActiveSet[nSubGrids] = 1.0/costActiveSet[nSubGrids];

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
bool SparseGrid<dim,out>::checkAccuracy(){
// absolute accuracy only for now
  double temp;
  double maxsurplus[out];
  for(int iout=0; iout<out; iout++){
    for(int iactiveErr=0; iactiveErr<activeHeapError.size(); iactiveErr++){
      if(activeSetAll[activeHeapError[iactiveErr]]){
        temp = error[activeHeapError[iactiveErr]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }

    for(int iactiveCost=0; iactiveCost<activeHeapCost.size(); iactiveCost++){
      if(activeSetAll[activeHeapCost[iactiveCost]]){
        temp = error[activeHeapCost[iactiveCost]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }
  }

  // for now assume that out = 1;
  double absAcc = maxsurplus[0];
  for(int iout=1; iout<out; iout++)
    if(absAcc<maxsurplus[iout]) absAcc = maxsurplus[iout];
  if(absAcc<absAccuracy){
    fprintf(stdout, "current absolute accuracy is %e and less than desired absolute accuracy %e\n", absAcc, absAccuracy);
    return true;
  }
  else return false;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//assumes that each coord is between 0 and 1.
template<int dim, int out>
inline
void SparseGrid<dim,out>::singleInterpolation(Coord coord, Output &output){

  int index = 0;
  int nPointsSubGrid;
  int index2[dim], index3;
  double scale; int xp;
  for(int iout=0; iout<out; iout++) output[iout] = 0.0;

  //fprintf(stderr, "### singleInterpolation -- nSubGrids = %d\n", nSubGrids);
  for(int subGrid=0; subGrid<nSubGrids; subGrid++){
    //for(int debugDim=0; debugDim<dim; debugDim++)
    //  fprintf(stderr, "###      for multiIndex[%d][%d] = %d\n", subGrid, debugDim, multiIndex[subGrid][debugDim]);
    // compute subGrid data structures
    nPointsSubGrid = 1;
    int nPointsDim[dim];
    int nCumulatedPointsDim[dim];
    for(int i=0; i<dim; i++){
      if(multiIndex[subGrid][i] == 0)
        nPointsDim[i] = 1;
      else if(multiIndex[subGrid][i] < 3)
        nPointsDim[i] = 2;
      else
        nPointsDim[i] = static_cast<int>(pow(2,multiIndex[subGrid][i]-1));
  
      nCumulatedPointsDim[i] = nPointsDim[i];
      nPointsSubGrid *= nPointsDim[i];
    }

    // interpolation
    double basisFnVal = 1.0;
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 1){
        if(coord[idim] == 1.0) index2[idim] = 1;
        else{
          xp = static_cast<int>(floor(coord[idim]*2));
          if(xp == 0) basisFnVal *= 2.0*(0.5 - coord[idim]);
          else        basisFnVal *= 2.0*(coord[idim] - 0.5);
          index2[idim] = xp;
        }
      }else if(multiIndex[subGrid][idim] == 0) index2[idim] = 0;
      else if(coord[idim] == 0.0){
        basisFnVal = 0.0;
        break;
      }else{
        scale = pow(2.0,multiIndex[subGrid][idim]);
        xp = static_cast<int>(floor(coord[idim] * scale / 2.0));
        basisFnVal *= (1.0 - scale*fabs(coord[idim]-(2.0*static_cast<double>(xp)+1.0)/scale));
      }

      if(basisFnVal == 0.0) break;

    }

    //fprintf(stderr, "### singleInterpolation -- basisFnVal = %e\n", basisFnVal);
    if(basisFnVal > 0.0){
      index3 = index + index2[0];
      for(int idim=1; idim<dim; idim++) index3 += nCumulatedPointsDim[idim-1]*index2[idim];
      for(int iout=0; iout<out; iout++) output[iout] += basisFnVal*surplus[index3][iout];
      //fprintf(stderr, "### singleInterpolation -- surplus[%d][0] = %e\n", index3, surplus[index3][0]);
    }
    index += nPointsSubGrid;
    //fprintf(stderr, "### singleInterpolation -- index = %d\n", index);

  }

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::interpolate(int numRes, Coord *coord,
                                      Output *res){

  Coord scaledCoord;
  for(int iPts=0; iPts<numRes; iPts++){
    scale(coord[iPts],scaledCoord, 1);
    //if(outOfRange(scaledCoord)) closestPointInRange(scaledCoord);
    singleInterpolation(scaledCoord, res[iPts]);
  }

}

//------------------------------------------------------------------------------
// SparseGrid::print and read sparse grid in a file
//------------------------------------------------------------------------------

template<int dim, int out>
void SparseGrid<dim,out>::printToFile(){

  // print to file the number of subgrids and the corresponding multiIndex array
  //           and the number of points and the surpluses.
  FILE *fpSPARSEGRID = fopen("SparseGrid", "w");
  if (!fpSPARSEGRID){
    fprintf(stderr, "*** Error: could not open 'SparseGrid' file\n");
    exit(1);
  }

  // header
  fprintf(fpSPARSEGRID, "# sparse grid: hierarchical subgrids and the surpluses\n");

  // characteristics used to create the sparse grid
  fprintf(fpSPARSEGRID, "%d %d\n", dim, out);
  fprintf(fpSPARSEGRID, "%d %d\n", minPoints, maxPoints);
  fprintf(fpSPARSEGRID, "%.12e %.12e %.12e\n", absAccuracy, relAccuracy, dimAdaptDegree);
  for(int idim=0; idim<dim; idim++)
    fprintf(fpSPARSEGRID, "%.12e %.12e ", range[idim][0], range[idim][1]);
  fprintf(fpSPARSEGRID, "\n");

  // the sparse grid data
  fprintf(fpSPARSEGRID, "%d\n", nSubGrids);
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
    for(int idim=0; idim<dim; idim++)
      fprintf(fpSPARSEGRID, "%d ", multiIndex[isubgrid][idim]);
    fprintf(fpSPARSEGRID, "\n"); 
  }
  
  fprintf(fpSPARSEGRID, "%d\n", nPoints);
  for(int ipts=0; ipts<nPoints; ipts++){
    for(int iout=0; iout<out; iout++)
      fprintf(fpSPARSEGRID, "%.12e ", surplus[ipts][iout]);
    fprintf(fpSPARSEGRID, "\n"); 
  }
  
  fflush(fpSPARSEGRID);
  fclose(fpSPARSEGRID);

}

//------------------------------------------------------------------------------

template<int dim, int out>
void SparseGrid<dim,out>::readFromFile(){

  char mystring [100];
  FILE *fpSPARSEGRID = fopen("SparseGrid", "r");
  if (!fpSPARSEGRID){
    fprintf(stderr, "*** Error: could not open 'SparseGrid' file to read\n");
    exit(1);
  }

  // header
  fgets(mystring, 100, fpSPARSEGRID);

  int readdim, readout;
  // characteristics used to create the sparse grid
  fscanf(fpSPARSEGRID, "%d %d\n", &readdim, &readout);
  assert(dim==readdim && out==readout);

  fscanf(fpSPARSEGRID, "%d %d", &minPoints, &maxPoints);
  fscanf(fpSPARSEGRID, "%lf %lf %lf", &absAccuracy, &relAccuracy, &dimAdaptDegree);
  for(int idim=0; idim<dim; idim++)
    fscanf(fpSPARSEGRID, "%lf %lf ", &(range[idim][0]), &(range[idim][1]));

  // the sparse grid data
  fscanf(fpSPARSEGRID, "%d", &nSubGrids);
  multiIndex = new MultiIndex[nSubGrids];
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
    for(int idim=0; idim<dim; idim++)
      fscanf(fpSPARSEGRID, "%d ", &(multiIndex[isubgrid][idim]));
  }

  fscanf(fpSPARSEGRID, "%d", &nPoints);
  surplus = new Output[nPoints];
  for(int ipts=0; ipts<nPoints; ipts++){
    for(int iout=0; iout<out; iout++)
      fscanf(fpSPARSEGRID, "%lf ", &(surplus[ipts][iout]));
  }

  fclose(fpSPARSEGRID);

  //fill in the rest
  sizeSurplus = nPoints;
  sizeMultiIndex = nSubGrids;

  fprintf(stdout, "????????????????????????????????????\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "min/maxPoints = %d %d\n", minPoints, maxPoints);
  fprintf(stdout, "abs/relAccuracy = %e %e\n", absAccuracy, relAccuracy);
  fprintf(stdout, "degree of dimensional adaptivity = %e\n", dimAdaptDegree);
  fprintf(stdout, "range1 is [%e %e]\n", range[0][0], range[0][1]);
  fprintf(stdout, "range2 is [%e %e]\n", range[1][0], range[1][1]);

}

//------------------------------------------------------------------------------
// SparseGrid::Heap implementation
//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::Heap::Heap(){

  size_    = 0;
  numElem_ = 0;
  elem_    = 0;

}

//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::Heap::~Heap(){

  delete [] elem_;

}

//------------------------------------------------------------------------------

template<int dim, int out>
SparseGrid<dim,out>::Heap::Heap(const Heap &heap){

  size_    = heap.size_;
  numElem_ = heap.numElem_;
  elem_    = new int[size_];
  for(int i=0; i<size_; i++)
    elem_[i] = heap.elem_[i];

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::Heap::sort(int index, double *value){

  if(index>numElem_-1 || index<0) return;

  int parent, temp;
  while(index>0){
    parent = static_cast<int>(floor((index-1)/2));
    if(value[elem_[index]]>value[elem_[parent]]){
      temp = elem_[parent];
      elem_[parent] = elem_[index];
      elem_[index] = temp;
    }else 
      break;
    index = parent;
  }

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
void SparseGrid<dim,out>::Heap::insert(int newElem, double *value){

  if(numElem_+1>size_){//resize
    if(size_>0) size_ *= 2;
    else size_++;
    int *temp = new int[size_];
    for(int i=0; i<numElem_; i++) temp[i] = elem_[i];
    delete [] elem_;
    elem_ = temp;
  }

  elem_[numElem_] = newElem;
  numElem_++;
  sort(numElem_-1,value);

}

//------------------------------------------------------------------------------

template<int dim, int out>
inline
int SparseGrid<dim,out>::Heap::pop(double *value){

  if(numElem_==0) return -1;

  int res = elem_[0];

  // reorder the heap 
  // first, fillup empty spot by children
  int i = 0;
  int firstChild, secondChild;
  while(2*i+1 < numElem_-1){ // condition for i to have to 2 children
    firstChild = 2*i+1;
    secondChild = firstChild + 1;
    if(value[elem_[firstChild]]>value[elem_[secondChild]]){
      elem_[i] = elem_[firstChild];
      i = firstChild;
    }else{
      elem_[i] = elem_[secondChild];
      i = secondChild;
    }
  }
  // second, put last element in last empty spot
  elem_[i] = elem_[numElem_-1];
  // third, make sure elem at position i is well positioned
  sort(i,value);

  numElem_--;

  return res;


}

//------------------------------------------------------------------------------
// TEST for DEBUG PURPOSES
//------------------------------------------------------------------------------

template<int dim, int out>
void SparseGrid<dim,out>::test(void (*fn)(double, double, double, double &)){


  assert(dim==2 && out==1);
/*
  // to test evaluatePreviousInterpolation and evaluateFunctionOnGrid
  //initialize
  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 
  Coord firstPoint; double res;
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5*(range[idim][0]+range[idim][1]);
  fn(firstPoint[0],firstPoint[1],0.0,res);
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

  errorActiveSet[0] = 0.1;
  errorActiveSet[1] = 0.2;
  errorActiveSet[2] = 0.3;
  errorActiveSet[3] = 0.4;
  errorActiveSet[4] = 0.5;
  errorActiveSet[5] = 0.6;
  errorActiveSet[6] = 0.7;
  errorActiveSet[7] = 0.8;
  errorActiveSet[8] = 0.9;

  activeHeapError.insert(0,errorActiveSet);
  activeHeapError.insert(1,errorActiveSet);
  activeHeapError.insert(2,errorActiveSet);
  activeHeapError.insert(3,errorActiveSet);
  activeHeapError.insert(4,errorActiveSet);
  activeHeapError.insert(5,errorActiveSet);
  activeHeapError.insert(6,errorActiveSet);
  activeHeapError.insert(7,errorActiveSet);

  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  
  int popped = activeHeapError.pop(errorActiveSet);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  popped = activeHeapError.pop(errorActiveSet);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  popped = activeHeapError.pop(errorActiveSet);
  fprintf(stderr, "retrieved element %d with multiIndex %d %d\n", popped, multiIndex[popped][0],multiIndex[popped][1]);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);

  activeHeapError.insert(1,errorActiveSet);
  fprintf(stderr, "\n");
  for(int toto=0; toto<activeHeapError.size(); toto++)
    fprintf(stderr, "element %d of the heap points to %d-th entry of the multiIndex\n", toto, activeHeapError[toto]);
  */
}
