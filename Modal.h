#ifndef _MODAL_H_
#define _MODAL_H_

#ifdef TYPE_PREC
#define PreScalar TYPE_PREC
#else
#define PreScalar double
#endif

class IoData;
class Domain;
template<int dim, class Scalar, int neq> class MatVecProdH2;
template<int dim, int neq> class MatVecProd;
template<int dim> class PostOperator;

#include <DistVector.h>
#include <VectorSet.h>
#include <SpaceOperator.h>
#include <Communicator.h>
#include <KspSolver.h>
#include <KspPrec.h>
#include <TsInput.h>
#include <TsOutput.h>
#include <TsRestart.h>
#include <GnatPreprocessing.h>
#include <GnatPreprocessingStep1.h>
#include <GnatPreprocessingStep2.h>
#include <SurfMeshGen.h>
#include <ReducedMeshShapeChanger.h>
#include <ParallelRom.h> 
#include <time.h> 

#ifdef DO_MODAL
  #include <arpack++/include/ardnsmat.h>
  #include <arpack++/include/ardssym.h>
#endif

#include <complex>
typedef std::complex<double> bcomp;

//-----------------------------------------------------------------------------------

template <int dim>
class ModalSolver {

    Communicator *com;
    Domain &domain;
    int nStrMode;
    int nPadeDeriv;
    IoData *ioData;
    DistBcData<dim> *bcData;
    DistGeoState *geoState;
    DistTimeState<dim> *tState;
    TsInput *tInput;
    TsOutput<dim> *tOutput;
    TsRestart *tRestart;

    SpaceOperator<dim> *spaceOp;
    PostOperator<dim> *postOp;
    //MatVecProd<dim, dim> *HOp;
    MatVecProdH2<dim, double, dim> *HOp;
    MatVecProdH2<dim, double, dim> *HOp2step1;
    MatVecProdH2<dim, double, dim> *HOp2;
    MatVecProdH2<dim, double, dim> *HOpstep2;
    MatVecProdH2<dim, double, dim> *HOpstep3;
    MatVecProdH2<dim, bcomp, dim> *HOpC;

    /*GmresSolver<DistSVec<double,dim>, MatVecProd<dim, 5>, KspPrec<dim, double>, Communicator> *ksp;
    GmresSolver<DistSVec<bcomp,dim>, MatVecProd<dim, 5>, KspPrec<dim, bcomp>, Communicator, bcomp> *kspComp;*/
    KspSolver<DistSVec<double,dim>, MatVecProd<dim, dim>, KspPrec<dim, double>, Communicator> *ksp;
    KspSolver<DistSVec<double,dim>, MatVecProd<dim, dim>, KspPrec<dim, double>, Communicator> *ksp2;
    KspSolver<DistSVec<double,dim>, MatVecProd<dim, dim>, KspPrec<dim, double>, Communicator> *ksp3;


    KspSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>, KspPrec<dim, bcomp>, Communicator, bcomp> *kspComp;
    GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>, KspPrec<dim, bcomp>, Communicator, bcomp> *kspCompGcr;

    KspPrec<dim, double> *pc;
    KspPrec<dim, bcomp> *pcComplex;
    KspData kspData;

    VecSet< DistSVec<double,3> > mX;
    DistSVec<double,dim> Uref;
    
    DistSVec<double,3> Xref;
    double dt;
    double dt0;
    double *K;
    int podMethod;

    VecSet< DistSVec<double,dim> > DX;
    VecSet< DistSVec<double,dim> > DE;
    DistVec<double> controlVol;
    DistVec<bcomp> controlVolComp;  // failing in this destructor
		double totalEnergy;
 
  public:
    ModalSolver(Communicator *, IoData &, Domain &);
    void preProcess();
    void createPODInTime();
    void constructPOD();
    void solve();
    void solveInTimeDomain();
    void solveInFreqDomain();
    void timeIntegrate(VecSet<DistSVec<double, dim> > &, int nSteps, int freq,
                       double *delU, double *delY, int &iSnap, double sdt, char *snapFile = 0);
    void freqIntegrate(VecSet<DistSVec<double, dim> > &, double, int &, VecSet<DistSVec<bcomp, dim > > &);
    void freqIntegrateMultipleRhs(VecSet<DistSVec<double, dim> > &, double, int &, VecSet<DistSVec<bcomp, dim > > &);
   void constructROM2(double *romOp, VecSet<Vec<double> > &romOp0, double *romOp1, double *romOp2, VecSet<Vec<double> > &ec,
                      VecSet<Vec<double> > &g, VecSet<DistSVec<double, dim> > &, int);

    void evalAeroSys(VecSet<Vec<double> > &, VecSet<DistSVec<double, dim> > &pVec, int);
    void evalFluidSys(VecSet<DistSVec<double, dim> > &pVec, int);
    void formOutputRom(VecSet<Vec<double> > &, VecSet<DistSVec<double, dim> > &, int nPodVecs);
    void timeIntegrateROM(double *romOp, VecSet<Vec<double> > &romOp0, double *romOp1, double *romOp2, VecSet<Vec<double> > &ecMat,
                          VecSet<Vec<double> > &gMat, VecSet<DistSVec<double, dim> > &podVecs,
                          int nSteps, int nPodVecs, double *delU, double *delY, double sdt);

void updateModalValues(double sdt, double *delU, double *delY, Vec<double> &modalF, int timeIt);

void computeModalDisp(double sdt, DistSVec<double, 3> &xPos, DistSVec<double, dim> &delW, double *delU, double *delY, Vec<double> &refModalF, int timeIt);

void computeModalDisp(double sdt, Vec<double> &delWRom, double *delU, double *delY, Vec<double> &refModalF, VecSet<Vec<double> > &PtimesPhi, int nPodVecs, int timeIt);

   // void computeModalDisp(double, DistSVec<double, 3> &, DistSVec<double, dim> &, double *, double *, Vec<double> &, VecSet<Vec<double> > &, Vec<double> &);
   // void computeModalDispStep1(double, DistSVec<double, 3> &, DistSVec<double, dim> &, double *, double *, Vec<double> &, VecSet<Vec<double> > &, Vec<double> &);
   //void computeModalDisp(double, DistSVec<double, 3> &, DistSVec<double, dim> &, double *, double *, Vec<double> &);
   // void computeModalDispStep1(double, DistSVec<double, 3> &, DistSVec<double, dim> &, double *, double *, Vec<double> &);
    void outputModalDisp(double *, double *, double, int, int, FILE *);
    void makeFreqPOD(VecSet<DistSVec<double, dim> > &, int, int, bool);
    void buildGlobalPOD();
    void wait(const int seconds);
    void normalizeSnap(DistSVec<double, dim>  &, const int, const int);
    //void projectFullSoltn();
    void interpolatePOD();
    template<class Scalar>
    void readPodVecs(VecSet<DistSVec<Scalar, dim> > &, int &);

    void checkROBType(VecSet<DistSVec<double, dim> > &, int );
#ifdef DO_MODAL
    void outputPODVectors(ARluSymStdEig<double> &podEigProb, VecSet<DistSVec<double, dim> > &, int nPod, int numSnaps);
#endif
    void outputPODVectors(VecSet<DistSVec<double, dim> > &U, Vec<double> &, int nPod);
    void computeRelativeEnergy(FILE *sValsFile, const Vec<double> &sVals, const int nPod);
};

#include "Modal.C"
#endif
