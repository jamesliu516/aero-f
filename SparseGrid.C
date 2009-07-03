#include <SparseGrid.h>
#include <math.h>



//------------------------------------------------------------------------------

SparseGrid::SparseGrid(){

  verbose = 0;

  dim = 0;
  out = 0;

  nPoints     = 0;
  sizeSurplus = 0;
  surplus     = 0;

  nSubGrids      = 0;
  sizeMultiIndex = 0;
  multiIndex     = 0;

  maxPoints = 100;
  minPoints = 100;
  absAccuracy = -1.0;
  relAccuracy = -1.0;
  range = 0;
  dimAdaptDegree = 0.0;

  nAdaptivePoints = 0;
  active = 0;

  activeError = 0;
  activeCost  = 0;

  neighbour = 0;
  error = 0;
  fnmin = 0;
  fnmax = 0;

}

//------------------------------------------------------------------------------

SparseGrid::~SparseGrid(){

  for(int i=0; i<sizeSurplus; i++) delete [] surplus[i];
  delete [] surplus;
  for(int i=0; i<sizeMultiIndex; i++) delete [] multiIndex[i];
  delete [] multiIndex;
  delete [] range;
  delete [] active;
  delete [] activeError;
  delete [] activeCost;
  if(neighbour){ // not initialized if sparse grid read from a file
    for(int i=0; i<sizeMultiIndex; i++) delete [] neighbour[i];
    delete [] neighbour;
  }
  if(error){ // not initialized if sparse grid read from a file
    for(int i=0; i<sizeMultiIndex; i++) delete [] error[i];
    delete [] error;
  }
  delete [] fnmin;
  delete [] fnmax;

}

//------------------------------------------------------------------------------

SparseGrid::SparseGrid(SparseGridData &data){

  dim = data.numInputs;
  out = data.numOutputs;
  verbose = data.verbose;

  maxPoints = data.maxPoints;
  minPoints = data.minPoints;
  absAccuracy = data.absAccuracy;
  relAccuracy = data.relAccuracy;
  range = new Range[dim];
  for(int idim=0; idim<dim; idim++){
    range[idim][0] = data.range[idim][0];
    range[idim][1] = data.range[idim][1];
  }
  dimAdaptDegree = data.dimAdaptDegree;

  nPoints         = 0;
  sizeSurplus     = maxPoints;
  surplus         = new double *[sizeSurplus];
  for(int i=0; i<sizeSurplus; i++) surplus[i] = new double[out];

  nSubGrids       = 0;
  sizeMultiIndex  = 10;
  multiIndex      = new int *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) multiIndex[i] = new int[dim];

  nAdaptivePoints = 0;
  active          = new bool[sizeMultiIndex];
  activeError  = new double[sizeMultiIndex];
  activeCost   = new double[sizeMultiIndex];

  neighbour       = new int *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) neighbour[i] = new int[2*dim];
  
  error           = new double *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) error[i] = new double[out];
  
  fnmin = new double[out];
  fnmax = new double[out];
}

//------------------------------------------------------------------------------

void SparseGrid::tabulate(void (*fn)(double * , double * )){

  messages(1);
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

}

//------------------------------------------------------------------------------

inline
void SparseGrid::initialize(void (*fn)(double *, double * )){
// initialization of the data structures for the first subgrid
// in the Clenshaw-Curtis grid, which contains only one point located
// at the center of the hypercube (this is level of refinement 0 in all directions).
  double temp = 0.0;

  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 

  double firstPoint[dim]; double res[out];
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5*(range[idim][0]+range[idim][1]);
  fn(firstPoint,res);
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

inline
int SparseGrid::currentSubGrid(bool &adaptivity){
// pick next subgrid: choose from the two available heaps

  bool done = false;
  int current  = -1;

  adaptivity = false;

  while(!done){
    current = -1;
    if(nAdaptivePoints>dimAdaptDegree*nPoints) // non-adaptive
      current = activeHeapCost.pop(activeCost);

    if(current<0){ // adaptive
      current = activeHeapError.pop(activeError);
      if(current>=0) adaptivity = true;
    }

    if(current<0)
      current = activeHeapCost.pop(activeCost);

    if(current<0) return current;

    // make sure it is really active
    if(active[current]) done = true;

  }

  return current;

}

//------------------------------------------------------------------------------

inline
bool SparseGrid::admissible(const int currentMultiIndex,
                                     const int forwardDir){

  // currentMultiIndex has multi-index multiIndex[currentMultiIndex]
  // To check the admissibility of the newMultiIndex 
  // (defined by the forward neighbour of currentMultiIndex in the 
  //  direction forwardDir, see sketch SparseGrid::findNeighbours),
  // all its backward neighbours must be inactive (integrated) subgrids.
  // In particular, its backward neighbour currentMultiIndex
  //   refers to an inactive subgrid. 

  int backwardOfCurrent, backwardOfNew;

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
    backwardOfNew = neighbour[backwardOfCurrent][forwardDir];
    if(backwardOfNew<0)       return false; // not part of sparse grid 
    if(active[backwardOfNew]) return false; // part of sparse grid, but active
    
  }

  return true;

}

//------------------------------------------------------------------------------

inline
void SparseGrid::findNeighbours(const int currentMultiIndex,
                                         const int forwardDir, const int addedSubGrids){
  int backwardOfCurrent, backwardOfNew;
  const int newMultiIndex = nSubGrids+addedSubGrids;

//    ____________________
//   |         |          |
//   | Current |   New    |  idim
//   |_________|__________|   ^
//   |         |          |   |
//   |  bOfC   |  bOfN    |   |
//   |_________|__________|   ----> forwardDir
//

  // find the neighbours
  for(int idim=0; idim<dim; idim++){
    if(forwardDir == idim){
      neighbour[currentMultiIndex][forwardDir    ] = newMultiIndex;
      neighbour[newMultiIndex    ][forwardDir+dim] = currentMultiIndex;
    }else{
      backwardOfCurrent = neighbour[currentMultiIndex][dim+idim];
      if(!(backwardOfCurrent<0)){ 
        backwardOfNew = neighbour[backwardOfCurrent][forwardDir];
        neighbour[backwardOfNew][idim    ] = newMultiIndex;
        neighbour[newMultiIndex][idim+dim] = backwardOfNew;
      }
    }
  }

  // get the levels of refinement of the new active subgrids
  for(int idim=0; idim<dim; idim++){
    backwardOfNew = neighbour[newMultiIndex][idim+dim];
    if(backwardOfNew>=0)
      multiIndex[newMultiIndex][idim] = multiIndex[backwardOfNew][idim]+1;
    else multiIndex[newMultiIndex][idim] = 0;
  }

}

//------------------------------------------------------------------------------

inline
int SparseGrid::integrateForwardAdmissible(const int newSubGridIndex,
                            void (*fn)(double * , double * )){

  int nPointsSubGrid;
  double ** subGrid = generateSubGrid(newSubGridIndex,nPointsSubGrid);
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  bounds(nPointsSubGrid);
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
inline
double **SparseGrid::generateSubGrid(const int newSubGrid,
                                     int &nPointsSubGrid){

  // find number of points in subgrid (Clenshaw-Curtis type)
  // and coordinates of the points in each direction.
  nPointsSubGrid = 1;
  int nPointsDim[dim]; // number of points in each dimension
  int nCumulatedPointsDim[dim]; // number of cumulated points in each dimension
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
  double **res = new double *[nPointsSubGrid];
  for(int i=0; i<nPointsSubGrid; i++)
    res[i] = new double[dim];
  tensorize(res,coordDim,nPointsDim,nCumulatedPointsDim);

  for(int i=0; i<dim; i++) delete [] coordDim[i];
  delete [] coordDim;

  return res;

}

//------------------------------------------------------------------------------

inline
void SparseGrid::tensorize(double **res, double ** coordDim,
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

inline
void SparseGrid::evaluateFunctionOnGrid(double **subGrid,
                                        const int nPointsSubGrid,
                                        void (*fn)(double * , double * )){

  while(nPoints+nPointsSubGrid > sizeSurplus) resizeSurplus();

  double scaledCoord[dim];
  double res[out];
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    scale(subGrid[iPts],scaledCoord, 0);
    fn(scaledCoord,res);
    for(int iout=0; iout<out; iout++)
      surplus[nPoints+iPts][iout] = res[iout];
  }

}

//------------------------------------------------------------------------------

inline
void SparseGrid::scale(double *subGrid, double *scaledCoord, int op){

  if(op==0)
    for(int i=0; i<dim; i++)
      scaledCoord[i] = subGrid[i]*(range[i][1]-range[i][0])+range[i][0];
  else if(op==1)
    for(int i=0; i<dim; i++)
      scaledCoord[i] = (subGrid[i]-range[i][0])/(range[i][1]-range[i][0]);

}

//------------------------------------------------------------------------------

inline
void SparseGrid::bounds(const int nPointsSubGrid){

  double temp;

  // computes the "errors" given by the surpluses for a subgrid in each output
  for(int iout=0; iout<out; iout++){
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      if(temp<fnmin[iout]) fnmin[iout] = temp;
      if(temp>fnmax[iout]) fnmax[iout] = temp;
    }
  }

}

//------------------------------------------------------------------------------

inline
void SparseGrid::evaluatePreviousInterpolation(double **subGrid,
                                                 const int nPointsSubGrid){

  double temp[out];
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    singleInterpolation(subGrid[iPts], temp);
    for(int iout=0; iout<out; iout++) 
      surplus[nPoints+iPts][iout] -= temp[iout];
  }

}

//------------------------------------------------------------------------------

inline
void SparseGrid::updateError(const int nPointsSubGrid){

  double indicator[out]; double temp;

  // computes the "errors" given by the surpluses for a subgrid in each output
  for(int iout=0; iout<out; iout++){
    indicator[iout] = 0.0;
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      indicator[iout] += ((temp>0) ? temp : -temp);
    }
  }

  activeError[nSubGrids] = 0.0; // = max(indicator)
  for(int iout=0; iout<out; iout++)
    if(indicator[iout] > activeError[nSubGrids])
      activeError[nSubGrids] = indicator[iout];

  // computes the "errors" for the accuracy check
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

inline
void SparseGrid::updateCost(){

  activeCost[nSubGrids] = 0.0;
  for(int idim=0; idim<dim; idim++)
    activeCost[nSubGrids] += static_cast<double>(multiIndex[nSubGrids][idim]);

  activeCost[nSubGrids] = 1.0/activeCost[nSubGrids];

}

//------------------------------------------------------------------------------

inline
bool SparseGrid::checkAccuracy(){

  double temp;
  double maxsurplus[out];
  for(int iout=0; iout<out; iout++){
    for(int iactiveErr=0; iactiveErr<activeHeapError.size(); iactiveErr++){
      if(active[activeHeapError[iactiveErr]]){
        temp = error[activeHeapError[iactiveErr]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }

    for(int iactiveCost=0; iactiveCost<activeHeapCost.size(); iactiveCost++){
      if(active[activeHeapCost[iactiveCost]]){
        temp = error[activeHeapCost[iactiveCost]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }
  }

  double absAcc = 0.0;
  double relAcc = 0.0;
  for(int iout=0; iout<out; iout++){
    if(absAcc<maxsurplus[iout]) absAcc = maxsurplus[iout];
    if(fnmax[iout]>fnmin[iout])
      if(relAcc<maxsurplus[iout]/(fnmax[iout]-fnmin[iout])) 
        relAcc = maxsurplus[iout]/(fnmax[iout]-fnmin[iout]);
  }


  if(absAcc<absAccuracy || relAcc<relAccuracy){
    messages(7);
    return true;
  }
  else return false;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//assumes that each coord is between 0 and 1.
inline
void SparseGrid::singleInterpolation(double *coord, double *output){

  int firstSurplus = 0; // points to first surplus of the considered subgrid
  int nPointsSubGrid;
  int surplusLocalCoord[dim], contributingSurplus;
  for(int iout=0; iout<out; iout++) output[iout] = 0.0;

  messages(8);

  for(int subGrid=0; subGrid<nSubGrids; subGrid++){

    // compute subGrid data structures
    nPointsSubGrid = 1;
    int nPointsDim[dim];
    int nCumulatedPointsDim[dim];
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 0)
        nPointsDim[idim] = 1;
      else if(multiIndex[subGrid][idim] < 3)
        nPointsDim[idim] = 2;
      else
        nPointsDim[idim] = static_cast<int>(pow(2,multiIndex[subGrid][idim]-1));
  
      nCumulatedPointsDim[idim] = nPointsDim[idim];
      nPointsSubGrid *= nPointsDim[idim];
    }

    // interpolation: value of the basis function at the considered subgrid (basisFnVal)
    //              : detection of the corresponding surplus in numbering of the
    //                  considered subgrid (surplusLocalCoord).
    double basisFnVal = 1.0;
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 1){
        if(coord[idim] == 1.0) surplusLocalCoord[idim] = 1;
        else{
          int xp = static_cast<int>(floor(coord[idim]*2));
          if(xp == 0) basisFnVal *= 2.0*(0.5 - coord[idim]);
          else        basisFnVal *= 2.0*(coord[idim] - 0.5);
          surplusLocalCoord[idim] = xp;
        }
      }else if(multiIndex[subGrid][idim] == 0) surplusLocalCoord[idim] = 0;
      else if(coord[idim] == 0.0){
        basisFnVal = 0.0;
        break;
      }else{
        double scale = pow(2.0,multiIndex[subGrid][idim]);
        int xp = static_cast<int>(floor(coord[idim] * scale / 2.0));
        basisFnVal *= (1.0 - scale*fabs(coord[idim]-(2.0*static_cast<double>(xp)+1.0)/scale));
        surplusLocalCoord[idim] = xp;
      }

      if(basisFnVal == 0.0) break;

    }

    // contributing surplus in the considered subgrid is added to result
    if(basisFnVal > 0.0){
      // position of the contributing surplus in the array of surpluses
      contributingSurplus = firstSurplus + surplusLocalCoord[0];
      for(int idim=1; idim<dim; idim++)
        contributingSurplus += nCumulatedPointsDim[idim-1]*surplusLocalCoord[idim];
      // contribution added
      for(int iout=0; iout<out; iout++)
        output[iout] += basisFnVal*surplus[contributingSurplus][iout];
    }
    firstSurplus += nPointsSubGrid;

  }

}

//------------------------------------------------------------------------------

void SparseGrid::interpolate(int numRes, double **coord, double **res){

  double scaledCoord[dim];
  for(int iPts=0; iPts<numRes; iPts++){
    scale(coord[iPts],scaledCoord, 1);
    if(outOfRange(scaledCoord)) closestPointInRange(scaledCoord);
    singleInterpolation(scaledCoord, res[iPts]);
  }

}

//------------------------------------------------------------------------------

inline
bool SparseGrid::outOfRange(double *coord){

  for(int idim=0; idim<dim; idim++)
    if(coord[idim]<0 || coord[idim]>1){
      fprintf(stdout, "*** Warning: coordinates(%e) in dimension %d is out of bounds [0,1]\n", coord[idim], idim);
      return true;
    }

  return false;

}

//------------------------------------------------------------------------------

inline
void SparseGrid::closestPointInRange(double *coord){

  for(int idim=0; idim<dim; idim++){
    if(coord[idim]<0) coord[idim] = 0.0;
    if(coord[idim]>1) coord[idim] = 1.0;
  }

}

//------------------------------------------------------------------------------
// SparseGrid::print and read sparse grid in a file
//------------------------------------------------------------------------------

void SparseGrid::printToFile(){

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

void SparseGrid::readFromFile(){

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
  fscanf(fpSPARSEGRID, "%d %d\n", &dim, &out);

  fscanf(fpSPARSEGRID, "%d %d", &minPoints, &maxPoints);
  fscanf(fpSPARSEGRID, "%lf %lf %lf", &absAccuracy, &relAccuracy, &dimAdaptDegree);
  range = new Range[dim];
  for(int idim=0; idim<dim; idim++)
    fscanf(fpSPARSEGRID, "%lf %lf ", &(range[idim][0]), &(range[idim][1]));

  // the sparse grid data
  fscanf(fpSPARSEGRID, "%d", &nSubGrids);
  multiIndex = new int *[nSubGrids];
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
    multiIndex[isubgrid] = new int[dim];
    for(int idim=0; idim<dim; idim++)
      fscanf(fpSPARSEGRID, "%d ", &(multiIndex[isubgrid][idim]));
  }

  fscanf(fpSPARSEGRID, "%d", &nPoints);
  surplus = new double *[nPoints];
  for(int ipts=0; ipts<nPoints; ipts++){
    surplus[ipts] = new double[out];
    for(int iout=0; iout<out; iout++)
      fscanf(fpSPARSEGRID, "%lf ", &(surplus[ipts][iout]));
  }

  fclose(fpSPARSEGRID);

  //fill in the rest
  sizeSurplus = nPoints;
  sizeMultiIndex = nSubGrids;

}

//------------------------------------------------------------------------------
// Auxiliary functions
//------------------------------------------------------------------------------

inline
void SparseGrid::resizeSurplus(){

  double **larger = new double *[2*sizeSurplus];
  for(int i=0; i<sizeSurplus; i++)
    for(int iout=0; iout<out; iout++)
      larger[i][iout] = surplus[i][iout];

  for(int i=0; i<sizeSurplus; i++) delete [] surplus[i];
  delete [] surplus;

  sizeSurplus *= 2;
  surplus = larger;

}

//------------------------------------------------------------------------------

inline
void SparseGrid::resizeMultiIndex(){

  int **largerMultiIndex = new int*[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerMultiIndex[i] = new int[dim];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int idim=0; idim<dim; idim++)
      largerMultiIndex[i][idim] = multiIndex[i][idim];

  for(int i=0; i<sizeMultiIndex; i++) delete [] multiIndex[i];
  delete [] multiIndex;
  multiIndex = largerMultiIndex;

  double *largerErrorActiveSet = new double[2*sizeMultiIndex];
  double *largerCostActiveSet = new double[2*sizeMultiIndex];
  bool *largerActive= new bool[2*sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++){
    largerErrorActiveSet[i]= activeError[i];
    largerCostActiveSet[i]= activeCost[i];
    largerActive[i]= active[i];
  }

  delete [] activeError; delete [] activeCost; delete [] active;
  activeError = largerErrorActiveSet;
  activeCost = largerCostActiveSet;
  active= largerActive;

  double **largerError = new double *[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerError[i] = new double[out];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int iout=0; iout<out; iout++)
      largerError[i][iout] = error[i][iout];

  for(int i=0; i<sizeMultiIndex; i++) delete [] error[i];
  delete [] error;
  error = largerError;

  int **largerNeighbour = new int *[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerNeighbour[i] = new int[2*dim];
  for(int i=0; i<2*sizeMultiIndex; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = -1;
  for(int i=0; i<nSubGrids; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = neighbour[i][ineigh];

  for(int i=0; i<sizeMultiIndex; i++) delete [] neighbour[i];
  delete [] neighbour;
  neighbour = largerNeighbour;

  sizeMultiIndex *= 2;

}

//------------------------------------------------------------------------------
// SparseGrid::Heap implementation
//------------------------------------------------------------------------------

SparseGrid::Heap::Heap(){

  size_    = 0;
  numElem_ = 0;
  elem_    = 0;

}

//------------------------------------------------------------------------------

SparseGrid::Heap::~Heap(){

  delete [] elem_;

}

//------------------------------------------------------------------------------

SparseGrid::Heap::Heap(const Heap &heap){

  size_    = heap.size_;
  numElem_ = heap.numElem_;
  elem_    = new int[size_];
  for(int i=0; i<size_; i++)
    elem_[i] = heap.elem_[i];

}

//------------------------------------------------------------------------------

inline
void SparseGrid::Heap::sort(int index, double *value){

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

inline
void SparseGrid::Heap::insert(int newElem, double *value){

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

inline
int SparseGrid::Heap::pop(double *value){

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
// functions for DEBUG PURPOSES
//------------------------------------------------------------------------------

void SparseGrid::messages(const int flag, const int arg){

  if(verbose==0) return;

  if(flag==1){
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "#     INITIALIZATION\n");
    fprintf(stdout, "###########################################\n");
  }

  if(flag==2){
    fprintf(stdout, "#     BEGIN ITERATION\n");
    fprintf(stdout, "###########################################\n");
    if(verbose>1){
      fprintf(stdout, "#     STATE AT BEGINNING OF ITERATION\n");
      fprintf(stdout, "# activity of grids is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %d", active[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# multiIndex is:");
      for(int idim = 0; idim<dim; idim++){
        for(int kk = 0; kk<nSubGrids; kk++)
          fprintf(stdout, "  %d", multiIndex[kk][idim]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "# activeError is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %e", activeError[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# activeCost is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %e", activeCost[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# errors are:");
      for(int iout = 0; iout<out; iout++){
        for(int kk = 0; kk<nSubGrids; kk++)
          fprintf(stdout, "  %e", error[kk][iout]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "\n");
      fprintf(stdout, "# surpluses are:");
      for(int iout = 0; iout<out; iout++){
        for(int kk = 0; kk<nSubGrids; kk++)
          fprintf(stdout, "  %e", surplus[kk][iout]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "\n###########################################\n");
    }
  }

  if(flag==3 && verbose>3){
    fprintf(stdout, "# adding grid %d to inactive grid set\n", arg);
    fprintf(stdout, "###########################################\n");
  }
  
  if(flag==4 && verbose>3){
    fprintf(stdout, "# neighbour of %d is admissible\n", arg);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==5 && verbose>3){
    fprintf(stdout, "# neighbours of %d are:", arg);
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      fprintf(stdout, "  %d", neighbour[arg][ineigh]);
    fprintf(stdout, "\n");
    fprintf(stdout, "###########################################\n");
  }

  if(flag==6 && verbose>3){
    fprintf(stdout, "# forward admissible grid %d leads to %d points and %d subgrids\n", arg, nPoints, nSubGrids);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==7 && verbose>3){
    fprintf(stdout, "# desired absolute accuracy %e has been reached\n", absAccuracy);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==8 && verbose>2){
    fprintf(stderr, "### singleInterpolation -- nSubGrids = %d\n", nSubGrids);
    for(int subGrid=0; subGrid<nSubGrids; subGrid++){
      for(int debugDim=0; debugDim<dim; debugDim++)
      fprintf(stderr, "###      for multiIndex[%d][%d] = %d\n", subGrid, debugDim, multiIndex[subGrid][debugDim]);
    }
  }

  if(flag==9){
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "# After reading sparse grid file,");
    fprintf(stdout, "# min/maxPoints = %d / %d\n", minPoints, maxPoints);
    fprintf(stdout, "# abs/relAccuracy = %e / %e\n", absAccuracy, relAccuracy);
    fprintf(stdout, "# degree of dimensional adaptivity = %e\n", dimAdaptDegree);
    for(int idim=0; idim<dim; idim++)
      fprintf(stdout, "# range[%d] is [%e %e]\n", idim, range[idim][0], range[idim][1]);
    fprintf(stdout, "###########################################\n");
  }

}

//------------------------------------------------------------------------------

void SparseGrid::test(void (*fn)(double * , double * )){


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
