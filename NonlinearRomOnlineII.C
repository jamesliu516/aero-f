#include <NonlinearRomOnlineII.h>
#include <Modal.h>
#include <TsInput.h>
#include <cmath>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

extern "C"      {
   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}

template<int dim> 
NonlinearRomOnlineII<dim>::NonlinearRomOnlineII(Communicator* _com, IoData& _ioData, Domain& _domain)  : 
  NonlinearRom<dim>(_com, _ioData, _domain)
{ 

  // this->ioData->example, this->com->example, this->domain.example

  if (this->ioData->problem.alltype != ProblemData::_NONLINEAR_ROM_OFFLINE_) { //projection error
    if (this->nClusters>1) readClosestCenterInfoModelII();
   
    if (this->ioData->romOnline.storeAllClusters==NonlinearRomOnlineData::STORE_ALL_CLUSTERS_TRUE)
      this->readAllClusteredOnlineQuantities();
  }


  if (this->ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_FAST_APPROX)
    this->readApproxMetricLowRankFactor("full"); 
}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRomOnlineII<dim>::~NonlinearRomOnlineII() 
{
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::readClosestCenterInfoModelII() {

  if (this->ioData->romOnline.distanceComparisons) {
    switch (this->ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        this->readDistanceComparisonInfo("noUpdates");
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        this->com->fprintf(stderr, "*** Error: fast distance comparisons are incompatible with simple ROB updates (use Exact)\n");
        exit(-1);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        this->readDistanceComparisonInfo("exactUpdates");
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        this->readDistanceComparisonInfo("approxUpdates");
        break;
      default:
        this->com->fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
    }
  } else {
    this->readClusterCenters("centers");
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::readClusteredOnlineQuantities(int iCluster) {

  // only need to read the state basis and update quantities for model II
  this->readClusteredBasis(iCluster, "state");

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::updateBasis(int iCluster, DistSVec<double, dim> &U) {

 switch (this->ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        this->com->fprintf(stdout, " ... Applying simple rank one basis update\n");
        updateBasisSimple(iCluster, U);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        this->com->fprintf(stderr, "*** Warning: Fast exact updates not yet implemented for Model II (using simple, which is also exact)\n");
        updateBasisSimple(iCluster, U);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        this->com->fprintf(stdout, " ... Applying rank one basis update using approximate metric\n");
        updateBasisFastApprox(iCluster, U);
        break;
      default:
        this->com->fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
   }


}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::updateBasisSimple(int iCluster, DistSVec<double, dim> &U) {

/* 
  When updateBasis is called the following quantities are available:
    - nClusters: number of clusters
    - snapsInCluster:  a vector containing the number of snapshots in each cluster, size nClusters
    - columnSumsV: a vector containing the sums of the columns of V, size ny + buffer (same buffer as basis and sVals?)
    - basis: the local ROB, ny + buffer
    - sVals: the singular values, size ny + buffer
  (As well as all of the ioData values)

*/
  this->readClusteredUpdateInfo(iCluster, "state");

  this->readClusteredReferenceState(iCluster, "state"); // reads Uref

  int robSize = this->basis->numVectors();
  int kSize = robSize+1;

  DistSVec<double, dim> a(this->domain.getNodeDistInfo());
  a = *(this->Uref) - U;
  
  if (a.norm() >= 1e-6) {  // only update if Uref is different than U (this handles the case of time=0) 

    double m[robSize];
    for (int iVec=0; iVec<robSize; ++iVec) {
      m[iVec] = (*(this->basis))[iVec] * a;
    }

    DistSVec<double, dim> p(this->domain.getNodeDistInfo());
    p = a;

    for (int iVec=0; iVec<robSize; ++iVec) {
      p -= (*(this->basis))[iVec] * m[iVec];
    }

    double Ra = p.norm();
    double RaInv = 1/Ra; 
    p *= RaInv;

    double *K = new double[(kSize)*(kSize)];

    for (int iCol = 0; iCol < (kSize); ++iCol){
      for (int iRow = 0; iRow < (kSize); ++iRow) {
        if ((iCol == iRow) && (iCol < kSize-1)) {
          K[iCol*(kSize) + iRow] = (*this->sVals)[iCol];
        } else {
          K[iCol*(kSize) + iRow] = 0.0;
        }
      }
    }

    double q = 0;
    for (int iVec=robSize; iVec<(this->columnSumsV->size()); ++iVec) {
      q += pow((*(this->columnSumsV))[iVec], 2);
    }
    q = pow(q, 0.5);

    for (int iRow = 0; iRow < (kSize-1); ++iRow) {
      for (int iCol = 0; iCol < (kSize-1); ++iCol){
        K[iCol*(kSize) + iRow] += m[iRow] * ((*(this->columnSumsV))[iCol]);
      }
      K[(kSize-1)*(kSize) + iRow] = m[iRow] * q;
    }

    for (int iCol = 0; iCol < kSize-1; ++iCol){
      K[iCol*(kSize) + (kSize-1)] += Ra * ((*(this->columnSumsV))[iCol]);
    }
    K[(kSize-1)*(kSize) + kSize-1] += Ra * q;

    double *sigma = new double[kSize];
    double *error = new double[kSize];
    double *work = new double[kSize];
    int info;
    double *zVec = new double[kSize*kSize]; // right singular vectors
    double *yVec = new double[kSize*kSize]; // left singular vectors

    this->com->fprintf(stdout, " ... computing rank one update to basis using current state\n");
    F77NAME(dsvdc)(K, kSize, kSize, kSize, sigma, error, yVec, kSize, zVec, kSize, work, 11, info);

    VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

    for (int iVec=0; iVec<robSize; ++iVec)
        basisOld[iVec] = (*(this->basis))[iVec];

    for (int iVec=0; iVec<robSize; ++iVec) {
      (*(this->basis))[iVec] = p * yVec[(iVec*kSize) + kSize-1];
      for (int jVec=0; jVec<robSize; ++jVec) {
        (*(this->basis))[iVec] += yVec[(iVec*kSize) + jVec] * basisOld[jVec];
      }
    }
  }
  delete (this->Uref);
  this->Uref = NULL;

}
//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::updateBasisFastApprox(int iCluster, DistSVec<double, dim> &U) {

  this->readClusteredUpdateInfo(iCluster, "state");

  this->readClusteredReferenceState(iCluster,"state"); // reads Uref

  int robSize = this->basis->numVectors();
  int kSize = robSize+1;

  DistSVec<double, dim> a(this->domain.getNodeDistInfo());
  a = *(this->Uref) - U;

  if (a.norm() >= 1e-6*U.norm()) {  // only update if Uref is different than U (this handles the case of time=0) 

    double m[robSize];
    double temp1[this->nLowRankFactors];
    double temp2[this->nLowRankFactors];
    for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank)
      temp1[iRank] = (*this->lowRankFactor)[iRank] * a;

    for (int iVec=0; iVec<robSize; ++iVec) {
      m[iVec]  = 0.0;
      for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank) {
        temp2[iRank] = (*this->lowRankFactor)[iRank] * (*(this->basis))[iVec];
        m[iVec] += (temp2[iRank]*temp1[iRank]);
      }
    }

    DistSVec<double, dim> p(this->domain.getNodeDistInfo());
    p = a;

    for (int iVec=0; iVec<robSize; ++iVec)
      p -= (*(this->basis))[iVec] * m[iVec];


    double Ra = 0.0;
    for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank) {
      temp1[iRank] = (*this->lowRankFactor)[iRank] * p;
      Ra += pow(temp1[iRank],2);
    }
    Ra = pow(Ra,0.5);
    double RaInv = 1/Ra;
    p *= RaInv;

    double *K = new double[(kSize)*(kSize)];
    for (int iCol = 0; iCol < (kSize); ++iCol){
      for (int iRow = 0; iRow < (kSize); ++iRow) {
        if ((iCol == iRow) && (iCol < kSize-1)) {
          K[iCol*(kSize) + iRow] = (*this->sVals)[iCol];
        } else {
          K[iCol*(kSize) + iRow] = 0.0;
        }
      }
    }

    double q = 0;
    for (int iVec=robSize; iVec<(this->columnSumsV->size()); ++iVec) {
      q += pow((*(this->columnSumsV))[iVec], 2);
    }
    q = pow(q, 0.5);

    for (int iRow = 0; iRow < (kSize-1); ++iRow) {
      for (int iCol = 0; iCol < (kSize-1); ++iCol){
        K[iCol*(kSize) + iRow] += m[iRow] * ((*(this->columnSumsV))[iCol]);
      }
      K[(kSize-1)*(kSize) + iRow] = m[iRow] * q;
    }

    for (int iCol = 0; iCol < kSize-1; ++iCol){
      K[iCol*(kSize) + (kSize-1)] += Ra * ((*(this->columnSumsV))[iCol]);
    }
    K[(kSize-1)*(kSize) + kSize-1] += Ra * q;

    double *sigma = new double[kSize];
    double *error = new double[kSize];
    double *work = new double[kSize];
    int info;
    double *zVec = new double[kSize*kSize]; // right singular vectors
    double *yVec = new double[kSize*kSize]; // left singular vectors

    this->com->fprintf(stdout, " ... computing rank one update to basis using current state\n");
    F77NAME(dsvdc)(K, kSize, kSize, kSize, sigma, error, yVec, kSize, zVec, kSize, work, 11, info);

    VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

    for (int iVec=0; iVec<robSize; ++iVec)
        basisOld[iVec] = (*(this->basis))[iVec];

    for (int iVec=0; iVec<robSize; ++iVec) {
      (*(this->basis))[iVec] = p * yVec[(iVec*kSize) + kSize-1];
      for (int jVec=0; jVec<robSize; ++jVec) {
        (*(this->basis))[iVec] += yVec[(iVec*kSize) + jVec] * basisOld[jVec];
      }
    }
  }

  delete (this->Uref);
  this->Uref = NULL;


  if (this->ioData->romOnline.distanceComparisons)
    this->resetDistanceComparisonQuantitiesApproxUpdates();   

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::appendNonStateDataToBasis(int cluster, char* basisType, bool relProjError) {

  int robSize = this->basis->numVectors();
  VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

  for (int iVec=0; iVec<robSize; ++iVec)
    basisOld[iVec] = (*(this->basis))[iVec];
 
  this->readClusteredBasis(cluster, basisType, relProjError);
  int nonStateSize = this->basis->numVectors();
  VecSet< DistSVec<double, dim> > nonStateBasis(nonStateSize, this->domain.getNodeDistInfo());

  for (int iVec=0; iVec<nonStateSize; ++iVec)
    nonStateBasis[iVec] = (*(this->basis))[iVec];

  this->basis->resize(robSize + nonStateSize);

  for (int iVec=0; iVec<robSize; ++iVec)
    (*(this->basis))[iVec] = basisOld[iVec];

  for (int iVec=robSize; iVec<(robSize+nonStateSize); iVec++)
    (*(this->basis))[iVec] = nonStateBasis[iVec-robSize];

  bool gramSchmidt;
  if (relProjError) {
    if (strcmp(basisType, "krylov")==0) {
      gramSchmidt = (this->ioData->romOffline.rob.relativeProjectionError.krylov.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
    } else if (strcmp(basisType, "sensitivity")==0) {
      gramSchmidt = (this->ioData->romOffline.rob.relativeProjectionError.sensitivity.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
    } else {
      this->com->fprintf(stderr, "*** Error: unexpected basis type passed to appendNonStateDataToBasis (%s)\n",basisType);
      exit(-1);
    }
  } else {
    if (strcmp(basisType, "krylov")==0) {
      gramSchmidt = (this->ioData->romOnline.krylov.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
    } else if (strcmp(basisType, "sensitivity")==0) {
      gramSchmidt = (this->ioData->romOnline.sensitivity.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
    } else {
      this->com->fprintf(stderr, "*** Error: unexpected basis type passed to appendNonStateDataToBasis (%s)\n",basisType);
      exit(-1);
    }
  }

  if (gramSchmidt) { 
    int uniqueVecs = robSize;
    for (int iVec = robSize; iVec<(robSize+nonStateSize); iVec++) {
      for (int jVec = 0; jVec<iVec; jVec++) {
        (*(this->basis))[iVec] -= (*(this->basis))[jVec] * ((*(this->basis))[iVec] * (*(this->basis))[jVec]);
        double norm = (*(this->basis))[iVec].norm();
        if (norm>=1e-14) {
          (*(this->basis))[iVec] *= 1/norm;
          ++uniqueVecs;
        } else {
          this->com->fprintf(stderr, "*** Warning: removing linearly dependent vector (norm=%e)\n",norm);
          (*(this->basis))[iVec] *= 0;    
        }
      }
    }
    if (uniqueVecs < (robSize+nonStateSize)) {
      basisOld.resize(robSize+nonStateSize);
      for (int iVec=0; iVec<(robSize+nonStateSize); ++iVec)
        basisOld[iVec] = (*(this->basis))[iVec];
      this->basis->resize(uniqueVecs);
      int added = 0;
      for (int iVec=0; iVec<(robSize+nonStateSize); ++iVec) {
        if (basisOld[iVec].norm()>=1e-14) {
          (*(this->basis))[added] = basisOld[iVec];
          ++added;
        }
      }
    }
  }
}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineII<dim>::appendVectorToBasis(DistSVec<double,dim> &vec, int numKeep) {

  this->com->fprintf(stdout, " ... appending Residual to ROB (absolute norm = %e)\n", vec.norm());
  
  int robSize = (numKeep == 0) ? this->basis->numVectors() : numKeep;

  VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

  for (int iVec=0; iVec<robSize; ++iVec)
    basisOld[iVec] = (*(this->basis))[iVec];

  int newRobSize = robSize + 1;

  this->basis->resize(newRobSize);

  for (int iVec=0; iVec<robSize; ++iVec)
    (*(this->basis))[iVec] = basisOld[iVec];

  (*(this->basis))[robSize] = vec;

  if (this->ioData->romOnline.onlineResiduals.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) {
    for (int iVec = 0; iVec<robSize; iVec++) {
      (*(this->basis))[robSize] -= (*(this->basis))[iVec] * ((*(this->basis))[iVec] * (*(this->basis))[robSize]);
    }  
    double norm = (*(this->basis))[robSize].norm();
    if (norm>=1e-10) {
      //this->com->fprintf(stdout, " ... norm = %e\n", norm);
      (*(this->basis))[robSize] *= 1/norm;
    } else {
      this->com->fprintf(stderr, "*** Warning: vector is already contained in basis (norm=%e)\n",norm);
      newRobSize = robSize;
      this->basis->resize(newRobSize);
      for (int iVec=0; iVec<newRobSize; ++iVec)
        (*(this->basis))[iVec] = basisOld[iVec];
    }
  }

  this->com->fprintf(stdout, " ... using ROB with %d vectors\n", newRobSize);

}


