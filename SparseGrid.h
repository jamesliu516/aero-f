#ifndef _SPARSE_GRID_H_
#define _SPARSE_GRID_H_

/*------------------------------------------------------------------------------
COMMENTS on Sparse Grids.
Clenshaw-Curtis grid
multilinear functions
dimensional adaptivity

Here, sparse grids are used to tabulate some target functions
with multilinear basis functions. Dimension adaptivity is
considered to select the subgrids (subspaces) that increase the most
the accuracy of the tabulation.

The following class has two purposes:
1 -- create the tabulation and write it in a file
2 -- read a tabulation from a file and use it for interpolation (instead of
     computing the target function)


Bibliography:
Gerstner T. and Griebel M.
Dimension-Adaptive Tensor-Product Quadrature
Computing, 2003

Klimke A.
Sparse Grid Interpolation Toolbox User's Guide (Matlab)

Klimke A. and Wohlmuth B. 
Spinter:piecewise multilinear hierarchical sparse
grid interpolation in Matlab. 
ACM Transactions on Mathematical Software, 2005.

Garcke J.
Sparse Grid Tutorial
2006
online


------------------------------------------------------------------------------*/


#include<IoData.h>




//------------------------------------------------------------------------------

template<int dim, int out>
class SparseGrid {

  int verbose;

  typedef double Output[out];
  typedef int    MultiIndex[dim];
  typedef double Coord[dim];
  typedef int    Neighbour[2*dim];

  int nPoints;            // number of stored points (on the sparse grid)
  Output *surplus;        // hierarchical surplus for each point for each output
  int sizeSurplus;

  int nSubGrids;          // number of stored subgrids
  MultiIndex *multiIndex; // multi-index for each stored subgrid
  int sizeMultiIndex;

// user-specified values for the creation of the sparse grid
  int maxPoints;          // max number of points in the grid
  int minPoints;          // min number of points in the grid
  double absAccuracy;     // required absolute accuracy
  double relAccuracy;     // required relative accuracy
  double range[dim][2];   // range of the tabulation (min and max in each dir)
  double dimAdaptDegree;  

// data structures necessary to construct a sparse grid

  double fnmin[out];      // minimum values of the interpolation
  double fnmax[out];      // maximum values of the interpolation

  Output *error;          // errors for each subgrid in each dimension
                          // used for accuracy check
                          // max_{pts of subgrid}(abs(surpluses[idim] of a subgrid))

// data structures necessary to construct a sparse grid using
// dimensional adaptivity

  int nAdaptivePoints;    // tracks the number of adaptive points during 
                          // creation of the sparse grid

  bool *active;           // for each subgrid, checks if it is active
                          // (same indexation as multiIndex)

  class Heap{
    int size_;
    int numElem_;
    int *elem_;

  public:
    Heap();
    ~Heap();
    Heap(const Heap &heap);
    Heap &operator=(const Heap &heap);
    int size(){ return numElem_; }
    int &operator[](int i){return elem_[i];}
    void sort(int index, double *value);     // sort after addition of a new element
                                             // placed at position 'index' and assuming
                                             // the rest was already sorted
    void insert(int newElem, double *value);
    int  pop(double *value);
  };

  Heap activeHeapError;   // heap that contains the indices of the active subgrids
                          // stored in multiIndex. The indices are ordered
                          // by the error of their corresponding subgrid.
  double *activeError;    // error of each subgrid (same indexation as multiIndex)
                          // max_{iout}(sum(abs(surpluses of a subgrid))

  Heap activeHeapCost;    // heap that contains the indices of the active subgrids
                          // stored in multiIndex. The indices are ordered
                          // by the cost of their corresponding subgrid.
  double *activeCost;     // cost of each subgrid (same indexation as multiIndex)

  Neighbour *neighbour;   // list of neighbours of a stored subgrid
                          // (same indexation as multiIndex and
                          // points to multiIndex)
                          // first the forward indices and then the backward ones

public:
  SparseGrid();
  ~SparseGrid();
  SparseGrid(SparseGridData &data);

  // functions to create the sparse grid with function fn to tabulate
  void tabulate(void (*fn)(Coord , Output ));
  void printToFile();


  // functions to perform interpolation on sparse grid
  void readFromFile();
  void interpolate(int numRes, Coord *coord, Output *res);

  // test function for debugging
  void test(void (*fn)(Coord , Output ));

private:
  SparseGrid(const SparseGrid<dim, out> &sparseGrid); // to prevent copying such objects
  SparseGrid& operator=(const SparseGrid<dim,out> &sparseGrid);

  void initialize(void (*fn)(Coord , Output ));
  int  currentSubGrid(bool &adaptivity);
  bool admissible(const int currentMultiIndex, const int forwardDir);
  void findNeighbours(const int currentMultiIndex, const int forwardDir, 
                      const int addedSubGrids);
  int integrateForwardAdmissible(const int newSubGrid, 
                      void (*fn)(Coord , Output ));
  Coord *generateSubGrid(const int newSubGrid, int &nPointsSubGrid);
  void evaluateFunctionOnGrid(Coord *subGrid, const int nPointsSubGrid, 
                              void (*fn)(Coord , Output ));
  void resizeSurplus();
  void resizeMultiIndex();
  void scale(Coord subGrid, Coord &scaledCoord, int op);
  void bounds(const int nPointsSubGrid);
  bool outOfRange(Coord coord);
  void closestPointInRange(Coord &coord);
  void evaluatePreviousInterpolation(Coord *subGrid, const int nPointsSubGrid);
  void updateError(const int nPointsSubGrid);
  void updateCost();
  void tensorize(Coord *res, double ** coordDim, int *nPtsDim, 
                 int *nCumulatedPtsDim);
  bool checkAccuracy();

  void singleInterpolation(Coord coord, Output &output);

  void messages(const int flag, const int arg=0);


};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SparseGrid.C>
#endif

#endif
