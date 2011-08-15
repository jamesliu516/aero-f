#ifndef _MVP_MATRIX_H_
#define _MVP_MATRIX_H_

#include <Vector.h>
#include <GenMatrix.h>
#include <cstdio>

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class MvpMat : public GenMat<Scalar,dim> {

  int n;
  int numEdges;

  SVec<Scalar,dim*dim> a;

  typedef std::map< std::pair<int,int>, Scalar* > AuxilliaryRows;
  AuxilliaryRows realAuxilliaryRows;
  AuxilliaryRows ghostAuxilliaryRows;
  AuxilliaryRows ghostGhostAuxilliaryRows;

public:

  MvpMat(int nn, int ne, int nBC = 0) : a(nn + 2*(ne+nBC)) { n = nn; numEdges = ne; }
  ~MvpMat() {}

  MvpMat<Scalar,dim> &operator= (const Scalar x) { a = x; return *this; }
  MvpMat<Scalar,dim> &operator*= (const Scalar x) { a *= x; return *this; }

  int numNonZeroBlocks() const { return a.size(); }
  Scalar (*data())[dim*dim] { return a.data(); }

  Scalar *getElem_ii(int i) { return *(a.v + i); }
  Scalar *getElem_ij(int l) { return *(a.v + n + 2 * l); }
  Scalar *getElem_ji(int l) { return *(a.v + n + 2 * l + 1); }

  Scalar *getBcElem_ij(int l) { return *(a.v + n + 2 * numEdges + 2 * l); }
  Scalar *getBcElem_ji(int l) { return *(a.v + n + 2 * numEdges + 2 * l + 1); }
  void addContrib(int nnd, int *nd, double *K) {}

  // -------------------------------------------------------------------------
  // Auxilliary terms (for ghost points)
 
  // Return the edge data corresponding to real node i and ghost node j
  Scalar* getRealNodeElem_ij(int i,int j) {
    return getAuxilliaryRow(realAuxilliaryRows,i,j);
  }

  // Return the edge data correponding to ghost node i and real node j
  Scalar* getGhostNodeElem_ij(int i,int j) {
    return getAuxilliaryRow(ghostAuxilliaryRows,i,j);
  }
  
  // Return the edge data correponding to ghost node i and ghost node j
  Scalar* getGhostGhostElem_ij(int i,int j) {
    return getAuxilliaryRow(ghostGhostAuxilliaryRows,i,j);
  }


  struct MvpAuxilliaryIterator : public GenMat<Scalar,dim>::AuxilliaryIterator {
    typename AuxilliaryRows::iterator itr;
    AuxilliaryRows* map_ptr;
  };

  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_realNodes() {
    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = realAuxilliaryRows.begin();
    mvpItr->map_ptr = &realAuxilliaryRows;
    updateIterator(mvpItr);
    return mvpItr;
  }

  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_ghostNodes() {
    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = ghostAuxilliaryRows.begin();
    mvpItr->map_ptr = &ghostAuxilliaryRows;
    updateIterator(mvpItr); 
    return mvpItr;
  }
  
  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_ghostGhostNodes() {
    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = ghostGhostAuxilliaryRows.begin();
    mvpItr->map_ptr = &ghostGhostAuxilliaryRows;
    updateIterator(mvpItr); 
    return mvpItr;
  }
  
  bool next(typename GenMat<Scalar,dim>::AuxilliaryIterator* genItr) { 
    MvpAuxilliaryIterator* mvpItr = static_cast< MvpAuxilliaryIterator*>(genItr);
    ++mvpItr->itr;
    if (mvpItr->itr == mvpItr->map_ptr->end())
      return false;
    else {
      updateIterator(mvpItr);
      return true;
    }
  }

  void free(typename GenMat<Scalar,dim>::AuxilliaryIterator* genItr) { 
    delete genItr;
  }
  
 protected:

  Scalar* getAuxilliaryRow(AuxilliaryRows& A, int i, int j) {
    std::pair<int,int> ij(i,j);
    typename AuxilliaryRows::iterator itr = A.find( ij );
    if (itr != A.end())
      return (A[ij] = new Scalar[dim*dim]);
    else
      return itr->second;
  }

  void updateIterator(MvpAuxilliaryIterator* mvpItr) {

    mvpItr->row = mvpItr->itr->first.first;
    mvpItr->col = mvpItr->itr->first.second;
    mvpItr->pData = mvpItr->itr->second;
  }
};

//------------------------------------------------------------------------------

#endif
