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

public:

  MvpMat(int nn, int ne, int nBC = 0) : a(nn + 2*(ne+nBC)) { n = nn; numEdges = ne; }
  ~MvpMat() {}

  MvpMat<Scalar,dim> &operator= (const Scalar x) { a = x; return *this; }
  MvpMat<Scalar,dim> &operator*= (const Scalar x) { a *= x; return *this; }

  double norm() {return a.norm();}

  int numNonZeroBlocks() const { return a.size(); }
  Scalar (*data())[dim*dim] { return a.data(); }

  Scalar *getElem_ii(int i) { return *(a.v + i); }
  Scalar *getElem_ij(int l) { return *(a.v + n + 2 * l); }
  Scalar *getElem_ji(int l) { return *(a.v + n + 2 * l + 1); }

  Scalar *getBcElem_ij(int l) { return *(a.v + n + 2 * numEdges + 2 * l); }
  Scalar *getBcElem_ji(int l) { return *(a.v + n + 2 * numEdges + 2 * l + 1); }
  void addContrib(int nnd, int *nd, double *K) {}

};

//------------------------------------------------------------------------------

#endif
