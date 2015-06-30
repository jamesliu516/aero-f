#ifndef _EMB_FLUID_SENSITIVITY_ANALYSIS_HANDLER_H_
#define _EMB_FLUID_SENSITIVITY_ANALYSIS_HANDLER_H_

#include <ImplicitEmbeddedCoupledTsDesc.h>
#include <string>

class IoData;
class Domain;
class GeoSource;

/*
template<int dim, int neq> class MatVecProd;
template<int dim, class Scalar2> class KspPrec; //???
template<class Scalar, int dim> class DistSVec;
template<int dim> class ImplicitEmbeddedCoupledTsDesc;
*/

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------
template<int dim>
class EmbeddedFluidSensitivityAnalysisHandler : public ImplicitEmbeddedCoupledTsDesc<dim> {

private:

  Domain *domain;

  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;

private:

  int step;
  int actvar;
  int numLocSub;

  double reynolds0;
  double kenergy0; 

  double length;
  double surface;
  double xmach;
  double alprad;
  double teta;
  double DFSPAR[3];

  DistVec<double> Pin;
  //DistVec<double> dAdS;
  DistVec<double> *Ap;
  DistVec<double> *Am;
 
  DistSVec<double,3> p;
  DistSVec<double,3> dPdS;
  DistSVec<double,3> dXdS;
  DistSVec<double,3> dXdSb;
  DistSVec<double,3> Xc;
  DistSVec<double,3> *Xp;
  DistSVec<double,3> *Xm;

  DistSVec<double,3> *Lp;
  DistSVec<double,3> *Lm;

  DistSVec<double,dim> Flux;
  DistSVec<double,dim> FluxFD;
  DistSVec<double,dim> *Fp;
  DistSVec<double,dim> *Fm;
  DistSVec<double,dim> dFdS;
  DistSVec<double,dim> *Up;
  DistSVec<double,dim> *Um;
  DistSVec<double,dim> dUdS;
  DistSVec<double,dim> Uc;

  DistSVec<double,3> Xplus;
  DistSVec<double,3> Xminus;
  DistSVec<double,3> dX;

  FILE* outFile;
  
public:

  EmbeddedFluidSensitivityAnalysisHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom
  );

  ~EmbeddedFluidSensitivityAnalysisHandler();


  void fsaPrintTextOnScreen(const char *);

  void fsaRestartBcFluxs(IoData &);

  /*
  void fsaGetEfforts(IoData &, DistSVec<double,3> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsaGetDerivativeOfEffortsFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double>&, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsaGetDerivativeOfEffortsAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);
  */

  void fsaSemiAnalytical
  (
   IoData &, 
   DistSVec<double,3> &, 
   DistVec<double> &, 
   DistSVec<double,dim> &, 
   DistSVec<double,dim> &
  );

  /*
  void fsaAnalytical
  (
    IoData &,
    DistSVec<double,3> &,
    DistVec<double> &,
    DistSVec<double,dim> &,
    DistSVec<double,dim> &
  );
  */

  void fsaSetUpLinearSolver
  (
   IoData &, 
   DistSVec<double,3> &, 
   DistVec<double> &, 
   DistSVec<double,dim> &, 
   DistSVec<double,dim> &
  );

  void fsaLinearSolver
  (
   IoData &, 
   DistSVec<double,dim> &, 
   DistSVec<double,dim> &
  );

  int fsaHandler(IoData &, DistSVec<double,dim> &);

  void fsaComputeDerivativesOfFluxAndSolution
  (
   IoData &, 
   DistSVec<double,3> &,
   DistVec<double> &, 
   DistSVec<double,dim> &
  );

  void fsaComputeSensitivities
  (
   IoData &, 
   const char *, 
   const char *, 
   DistSVec<double,3> &, 
   DistSVec<double,dim> &
  );

};


//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <EmbeddedFluidSensitivityAnalysisHandler.C>
#endif

#endif
