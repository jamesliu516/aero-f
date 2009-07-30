#ifndef SPARSEGRID_HPP_
#define SPARSEGRID_HPP_

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

#include "IoData.hpp"

//------------------------------------------------------------------------------

class SparseGrid {
  double *parameters;

  int dim, out;           // number of inputs and outputs
  
  int verbose;

  typedef double Range[2];

  int nPoints;            // number of stored points (on the sparse grid)
  double **surplus;       // hierarchical surplus for each point for each output
  int sizeSurplus;

  int nSubGrids;          // number of stored subgrids
  int **multiIndex;       // multi-index for each stored subgrid
  int sizeMultiIndex;

// user-specified values for the creation of the sparse grid
  int maxPoints;          // max number of points in the grid
  int minPoints;          // min number of points in the grid
  double absAccuracy;     // required absolute accuracy
  double relAccuracy;     // required relative accuracy
  Range *range;           // range of the tabulation (min and max in each dir)
  double dimAdaptDegree;  // degree of dimensional adaptivity

// data structures necessary to construct a sparse grid

  double *fnmin;          // minimum values of the interpolation
  double *fnmax;          // maximum values of the interpolation

  double **error;         // errors for each subgrid in each dimension
                          // used for accuracy check
                          // max_{pts of subgrid}(abs(surpluses[idim] of a subgrid))

// data structures necessary to construct a sparse grid using
// dimensional adaptivity

  int nAdaptivePoints;    // tracks the number of adaptive points during 
                          // creation of the sparse grid

  bool *active;           // for each subgrid, checks if it is active
                          // (same indexation as multiIndex)

  class Heap{             // customed heap class
    int size_;
    int numElem_;
    int *elem_;

  public:
    Heap();
    ~Heap();
    Heap(const Heap &heap);
    Heap &operator=(const Heap &heap);
    int size() const{ return numElem_; }
    int &operator[](int i) const {return elem_[i];}
    void sort(int index, const double *value);     // sort after addition of a new element
                                             // placed at position 'index' and assuming
                                             // the rest was already sorted
    void insert(const int newElem, const double *value);
    int  pop(const double *value);
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

  int **neighbour;        // list of neighbours of a stored subgrid
                          // (same indexation as multiIndex and
                          // points to multiIndex)
                          // first the forward indices and then the backward ones

protected:
  //Adapter function object in order to tabulate member functions of a class
  template<typename T>
  class Functor{
  public:
    typedef void (T::*MemberFunctionType)(double *, double *, double *);

    Functor(MemberFunctionType memberfn, T &object) : 
      memberfn_(memberfn), object_(object){}

    void operator()(double *argin, double *argout, double *argparameters){
      (object_.*memberfn_)(argin, argout, argparameters);
    }

  private:
    MemberFunctionType memberfn_;
    T &object_;
  };


public:
  SparseGrid();
  ~SparseGrid();
  SparseGrid(SparseGridData &data, double *param);

  // one of the two following functions in order to create a tabulation
  // choice depends if it is a member function of a class or not
  template<typename T>
  void tabulate(void (T::*fn)(double *, double *, double *), T &object);
  
  template<typename FnType>
  void tabulate(FnType fn);
  
  // prints the tabulation in an ASCII file
  // which can be later read by readFromFile()
  void printToFile() const;

  // functions to perform interpolation on sparse grid
  void readFromFile();
  void interpolate(const int numRes, double **coord, double **res);

  // test function for debugging
  template<typename FnType>
  void test(FnType fn);

private:
  SparseGrid(const SparseGrid &sparseGrid); // to prevent copying such objects
  SparseGrid& operator=(const SparseGrid &sparseGrid);

  template<typename FnType>
  void initialize(FnType fn);
  int  currentSubGrid(bool &adaptivity);
  bool admissible(const int currentMultiIndex, const int forwardDir);
  void findNeighbours(const int currentMultiIndex, const int forwardDir, 
                      const int addedSubGrids);
  template<typename FnType>
  int integrateForwardAdmissible(const int newSubGrid, FnType fn);
  double **generateSubGrid(const int newSubGrid, int &nPointsSubGrid) const;
  void tensorize(double **res, double ** coordDim, const int *nPtsDim, 
                 const int *nCumulatedPtsDim) const;
  template<typename FnType>
  void evaluateFunctionOnGrid(double **subGrid, const int nPointsSubGrid, FnType fn);
  void evaluatePreviousInterpolation(double **subGrid, const int nPointsSubGrid);
  void updateError(const int nPointsSubGrid);
  void updateCost();
  bool checkAccuracy();

  void singleInterpolation(const double * coord, double * output) const;

  // depending on verbose, outputs messages on standard output
  void messages(const int flag, const int arg=0) const;

  void scale(const double *subGrid, double *scaledCoord, const int op) const;
  void bounds(const int nPointsSubGrid);
  bool outOfRange(double *coord) const;
  void closestPointInRange(double *coord);

  // resizing data structures if necessary  
  void resizeSurplus();
  void resizeMultiIndex();
  
};

//------------------------------------------------------------------------------
#ifdef TEMPLATE_FIX
#include "SparseGrid.cpp"
#endif

#endif /*SPARSEGRID_HPP_*/
