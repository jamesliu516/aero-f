#ifndef _GEN_MATRIX_H_
#define _GEN_MATRIX_H_

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class GenMat {

public:

  GenMat() {}
  ~GenMat() {}

  virtual GenMat<Scalar,dim> &operator= (const Scalar) = 0;
  virtual GenMat<Scalar,dim> &operator*= (const Scalar) = 0;

  virtual Scalar (*data())[dim*dim] = 0;

  virtual Scalar *getElem_ii(int) = 0;
  virtual Scalar *getElem_ij(int) = 0;
  virtual Scalar *getElem_ji(int) = 0;

  virtual Scalar *getBcElem_ij(int l) {fprintf(stderr, "No Implementation\n");}
  virtual Scalar *getBcElem_ji(int l) {fprintf(stderr, "No Implementation\n");}
  virtual void addContrib(int, int *, double *) = 0;

};

//------------------------------------------------------------------------------

#endif
