#ifndef _DIST_MATRIX_H_
#define _DIST_MATRIX_H_

#include <Domain.h>
#include <SubDomain.h>
#include <GenMatrix.h>
#include <Communicator.h>

#include <complex.h>
typedef complex<double> bcomp;

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class DistMat {

protected:

  int numLocSub;

  CommPattern<double> *vecPat;
  CommPattern<bcomp> *compVecPat;

  CommPattern<Scalar> *diagMatPat;

  CommPattern<Scalar> *offDiagMatPat;

  SubDomain **subDomain;

  Communicator *com;

public:

  DistMat(Domain *);
  ~DistMat();

  virtual DistMat<Scalar,dim> &operator= (const Scalar) = 0;

  virtual GenMat<Scalar,dim> &operator() (int) = 0;

  CommPattern<Scalar> *getDiagMatPat() const { return diagMatPat; }

  CommPattern<Scalar> *getOffDiagMatPat() const { return offDiagMatPat; }

  CommPattern<double> *getCommPat(DistSVec<double, dim> &)  { return vecPat; }
  CommPattern<bcomp> *getCommPat(DistSVec<bcomp, dim> &)  { return compVecPat; }
};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistMat<Scalar,dim>::DistMat(Domain *domain)
{

  numLocSub = domain->getNumLocSub();

  subDomain = domain->getSubDomain();

  com = domain->getCommunicator();

  vecPat = new CommPattern<double>(domain->getSubTopo(), com, CommPattern<double>::CopyOnSend);
  compVecPat = new CommPattern<bcomp>(domain->getSubTopo(), com, CommPattern<bcomp>::CopyOnSend);

  diagMatPat = new CommPattern<Scalar>(domain->getSubTopo(), com, 
				       CommPattern<Scalar>::CopyOnSend);

  offDiagMatPat = new CommPattern<Scalar>(domain->getSubTopo(), com, 
					  CommPattern<Scalar>::CopyOnSend);

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->setComLenNodes(dim, *vecPat);
    subDomain[iSub]->setComLenNodes(dim, *compVecPat);

    subDomain[iSub]->setComLenNodes(dim*dim, *diagMatPat);
    subDomain[iSub]->setComLenEdges(2*dim*dim, *offDiagMatPat);
  }

  vecPat->finalize();
  compVecPat->finalize();
  diagMatPat->finalize();
  offDiagMatPat->finalize();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistMat<Scalar,dim>::~DistMat()
{

  if (vecPat) delete vecPat;
  if (compVecPat) delete compVecPat;
  if (diagMatPat) delete diagMatPat;
  if (offDiagMatPat) delete offDiagMatPat;

}

//------------------------------------------------------------------------------

#endif
