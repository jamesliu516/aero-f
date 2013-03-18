#ifndef _KSP_BINARY_OUTPUT_H_
#define _KSP_BINARY_OUTPUT_H_

#include <cstdio>
#include <Vector.h>
#include <DenseMatrix.h>
#include <VectorSet.h>
#include <DistVector.h>
#include <DistEmbeddedVector.h>
#include <string.h>

template <class VecType>
class KspBinaryOutput {

//This class is used to output the krylov vectors from the linear solver (GMRES only) 

  protected:
  Domain* domain;
  Communicator* com; 
  IoData* ioData;

  int krylovFreqTime;
  int krylovFreqNewton;
  double krylovEnergy;
  int* timeIt;
  int* newtonIt;
  char* fileName;

  public:

  KspBinaryOutput(Communicator *, IoData *, Domain *);
  ~KspBinaryOutput();

  template <int dim>
  void writeKrylovVectors(VecSet<DistSVec<double, dim> >&, Vec<double>, int);

  template <int dim>
  void writeKrylovVectors(VecSet<DistSVec<bcomp, dim> >&, Vec<bcomp>, int);

  template <int dim, class Scalar>
  void writeKrylovVectors(VecSet<DistEmbeddedVec<Scalar, dim> >&, Vec<Scalar>, int);
};

#include "KspBinaryOutput.C"
#endif
