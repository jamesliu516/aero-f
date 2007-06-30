#ifndef _BC_FCN_H_
#define _BC_FCN_H_

#include<complex.h>
typedef complex<double> bcomp;
class IoData;

//------------------------------------------------------------------------------

class BcFcn {

public:

  BcFcn() {}
  ~BcFcn() {}

  virtual void applyToSolutionVector(int, double *, double *);
  virtual void applyToResidualTerm(int, double *, double *, double *);
  virtual void applyToDiagonalTerm(int, double *, double *, float *);
  virtual void applyToDiagonalTerm(int, double *, double *, double *);
  virtual void applyToDiagonalTerm(int, double *, double *, bcomp *);
  virtual void applyToOffDiagonalTerm(int, float *);
  virtual void applyToOffDiagonalTerm(int, double *);
  virtual void applyToOffDiagonalTerm(int, bcomp *);

  virtual void zeroDiagonalTerm(int, float *);
  virtual void zeroDiagonalTerm(int, double *);
  virtual void zeroDiagonalTerm(int, bcomp *);

};

//------------------------------------------------------------------------------

class BcFcnNS : public BcFcn {

public:

  BcFcnNS() {}
  ~BcFcnNS() {}

  static void template_applyToSolutionVectorTerm(int, double *, double *);

  static void template_applyToResidualTerm(int, double *, double *, double *);

  template<class Scalar, int neq>
  static void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
  
  template<class Scalar, int neq>
  static void template_applyToOffDiagonalTerm(int, Scalar *);

  void applyToSolutionVector(int, double *, double *);
  void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void zeroDiagonalTerm(int, float *);
  void zeroDiagonalTerm(int, double *);
  void zeroDiagonalTerm(int, bcomp *);
  template<class Scalar, int neq>
  static void template_zeroDiagonalTerm(int, Scalar *);

};

//------------------------------------------------------------------------------
/*
@BOOK{white-74,
  author = "White, F. M.",
  title = "Viscous fluid flow",
  publisher = "McGraw-Hill",
  year = 1974,
} 
*/

class BcFcnSA : public BcFcn {

  bool wallFcn;

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

public:

  BcFcnSA(IoData&);
  ~BcFcnSA() {}

  void applyToSolutionVector(int, double *, double *);
  void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

};

//------------------------------------------------------------------------------

class BcFcnSAturb : public BcFcn {

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

public:

  BcFcnSAturb() {}
  ~BcFcnSAturb() {}

  void applyToSolutionVector(int, double *, double *);
  //void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

};

//------------------------------------------------------------------------------
/*
@ARTICLE{jaeger-dhatt-92,
  author = "Jaeger, M. and Dhatt, G.",
  title = "An extended k--$\epsilon$ finite element model",
  journal = ijnmf,
  year = 1992,
  volume = 14,
  pages = "1325--1345",
} 
*/

class BcFcnKE : public BcFcn {

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

public:

  BcFcnKE() {}
  ~BcFcnKE() {}

  void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);
};

//------------------------------------------------------------------------------

class BcFcnKEturb : public BcFcn {

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

public:

  BcFcnKEturb() {}
  ~BcFcnKEturb() {}

  //void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <BcFcn.C>
#endif

#endif
