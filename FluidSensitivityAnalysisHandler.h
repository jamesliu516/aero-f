#ifndef _FLUID_SENSITIVITY_ANALYSIS_HANDLER_H_
#define _FLUID_SENSITIVITY_ANALYSIS_HANDLER_H_

#include <ImplicitCoupledTsDesc.h>

#include <string>

class IoData;
class Domain;
class GeoSource;
class MeshMotionSolver;
//class StructExc;
class RigidMeshMotionHandler;

template<int dim, int neq> class MatVecProd;
template<int dim, class Scalar2> class KspPrec;
template<class Scalar, int dim> class DistSVec;
template<int dim> class ImplicitCoupledTsDesc;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------
template<int dim>
class FluidSensitivityAnalysisHandler : public ImplicitCoupledTsDesc<dim> {

private:

  Domain *domain;
  MeshMotionSolver *mms;

  // UH (08/10) This pointer is never used.
  //StructExc *strExc;

  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;

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
  DistVec<double> dAdS;
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
  DistSVec<double,3> *Z;

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

  /// Constructor
  /// \param[in] ioData  Reference to an 'IoData' object.
  /// \param[in] geoSource  Reference to a 'GeoSource' object.
  /// \param[in] dom  Pointer to the domain.
  /// \note The pointer 'dom' is passed to ImplicitCoupledTsDesc
  /// and stored in the member variable domain.
  /// It seems redundant (UH - 08/10)
  FluidSensitivityAnalysisHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom
  );

  ~FluidSensitivityAnalysisHandler();

  /// \note This function is implemented but never called.
  void fsaOutput1D(const char *, DistVec<double> &);

  /// \note This function is implemented but never called.
  void fsaOutput3D(const char *, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsaOutputDimD(const char *, DistSVec<double,dim> &);

  void fsaPrintTextOnScreen(const char *);

  void fsaRestartBcFluxs(IoData &);

  void fsaGetEfforts(IoData &, DistSVec<double,3> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsaGetDerivativeOfEffortsFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double>&, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsaGetDerivativeOfEffortsAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  /// \note This function is implemented but never called.
  void fsaGetDerivativeOfLoadFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsaGetDerivativeOfLoadAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  void fsaSemiAnalytical
  (
    IoData &,
    DistSVec<double,3> &,
    DistVec<double> &,
    DistSVec<double,dim> &,
    DistSVec<double,dim> &
  );

  void fsaAnalytical
  (
    IoData &,
    DistSVec<double,3> &,
    DistVec<double> &,
    DistSVec<double,dim> &,
    DistSVec<double,dim> &
  );

  void fsaSetUpLinearSolver(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  void fsaLinearSolver
  (
    IoData &, DistSVec<double,dim> &, DistSVec<double,dim> &
  );

  int fsaHandler
  (
    IoData &,
    DistSVec<double,dim> &
  );

  void fsaComputeDerivativesOfFluxAndSolution
  (
    IoData &,
    DistSVec<double,3> &,
    DistVec<double> &,
    DistSVec<double,dim> &
  );

  void fsaComputeSensitivities(IoData &, const char *, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);

};


//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <FluidSensitivityAnalysisHandler.C>
#endif

#endif

