#ifndef _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_
#define _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_

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
class FluidShapeOptimizationHandler : public ImplicitCoupledTsDesc<dim> {

private:

  Domain *domain;
//	TsSolver<ImplicitCoupledTsDesc<dim> > *tsSolver;
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
	DistSVec<double,3> dXb;
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
  FluidShapeOptimizationHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom//,
//		TsSolver<ImplicitCoupledTsDesc<dim> > *_tsSolver
  );

  ~FluidShapeOptimizationHandler();

  /// \note This function is implemented but never called.
  void fsoOutput1D(const char *, DistVec<double> &);

  /// \note This function is implemented but never called.
  void fsoOutput3D(const char *, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoOutputDimD(const char *, DistSVec<double,dim> &);

  void fsoPrintTextOnScreen(const char *);

  void fsoRestartBcFluxs(IoData &);

  void fsoGetEfforts(IoData &, DistSVec<double,3> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double>&, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  void fsoSemiAnalytical(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoAnalytical(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoSetUpLinearSolver(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoLinearSolver(IoData &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  int fsoHandler(IoData &, DistSVec<double,dim> &);
  void fsoComputeDerivativesOfFluxAndSolution(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);
  void fsoComputeSensitivities(IoData &, const char *, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);
	void fsoMoveMesh(IoData &ioData, DistSVec<double,dim> &U);

};


//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <FluidShapeOptimizationHandler.C>
#endif

#endif


