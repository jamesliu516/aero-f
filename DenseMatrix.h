#ifndef _DENSE_MATRIX_H_
#define _DENSE_MATRIX_H_

// GenFullM = Full Matrix class
//         stores an mxn matrix
//         certain member functions only work for
//         the square matrix case (nxn)

#include <Vector.h>

typedef Vec<double> Vector;

template<class Scalar>
class GenFullM {
 protected:
   int nrow;	// number of rows
   int ncolumn; // number of columns
   Scalar *v;   // pointer to matrix data

 public:

   // constructors
   GenFullM<Scalar>(); // Creates an empty matrix
   GenFullM(int _nr);
   GenFullM(int _nr, int _nc);
   GenFullM(const GenFullM &,int _nr, int sr, int _nc, int sc);
   GenFullM(const GenFullM &);

   // destructor
   ~GenFullM();

   void setNewSize(int _nr, int _nc, double d=0.0);
   void setNewSize(int _nr, double d=0.0);

   // OPERATORS
   void  operator = (const GenFullM &);
   void  operator = (const Scalar c);
   void  operator *= (const Scalar c);
   GenFullM<Scalar> operator *(GenFullM<Scalar>&);
   GenFullM<Scalar> operator ^(GenFullM<Scalar>&); // product A^T*B
   GenFullM<Scalar> operator %(GenFullM<Scalar>&); // product A*B^T

//   Vector operator *(const Vector &v);

//   GenFullM invert();
//   GenFullM transpose();

   int dim()    { return nrow;    }
   int numRow() { return nrow;    }
   int numCol() { return ncolumn; }

   Scalar *operator[](int i) const;
   Scalar* data() const { return v; }

//   double max();
   void print(char *msg = "");
   void factor();
   void reSolve(double *d);
   void zero();
   void add(GenFullM&, int, int);

//   void transposeAssign(GenFullM&);

//   void transposeMult(GenFullM&, GenFullM&);
};

template<class Scalar>
inline
Scalar *
GenFullM<Scalar>::operator[](int i) const
 { return v+i*ncolumn; }

/*
class StackFullM : public GenFullM {
 public:
   StackFullM(int nr, int nc, double *data);
   ~StackFullM() { v = 0; }
};

template<class Scalar>
inline
StackFullM::StackFullM(int nr, int nc, double *data)
{
 nrow    = nr;
 ncolumn = nc;
 v       = data;
}

*/
typedef GenFullM<double> FullM;
#ifdef TEMPLATE_FIX
#include <DenseMatrix.C>
#endif
#endif
