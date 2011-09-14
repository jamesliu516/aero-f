#ifndef _TS_OUTPUT_H_
#define _TS_OUTPUT_H_

#include <PostFcn.h>
#include <Vector3D.h>

#include <cstdio>

class IoData;
class RefVal;
class Domain;
class MeshMotionHandler;
class RigidMeshMotionHandler;
class HeavingMeshMotionHandler;
class PitchingMeshMotionHandler;
class DeformingMeshMotionHandler;
class AccMeshMotionHandler;
class Communicator;

template<int dimLS> class LevelSet;
template<int dim> class PostOperator;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistTimeState;
template<int dim> class DistExactRiemannSolver;
//------------------------------------------------------------------------------

template<int dim>
class TsOutput {

private:

  RefVal *refVal;
  PostOperator<dim> *postOp;
  RigidMeshMotionHandler *rmmh;
  HeavingMeshMotionHandler *hmmh;
  PitchingMeshMotionHandler *pmmh;
  DeformingMeshMotionHandler *dmmh;
  Domain *domain;
  Communicator *com;

private:

  bool steady;
  int it0;
  int frequency;
  double frequency_dt, prtout;
  int numFluidPhases; //excludes "ghost" solids
  double length;
  double surface;
  static int counter;
  Vec3D x0;
  int *output_newton_step;	// points to domain's

  double sscale[PostFcn::SSIZE];
  double vscale[PostFcn::SSIZE];
  double avsscale[PostFcn::AVSSIZE];
  double avvscale[PostFcn::AVVSIZE];

  char *solutions;
  char *scalars[PostFcn::SSIZE];
  char *vectors[PostFcn::VSIZE];
  char *avscalars[PostFcn::AVSSIZE];
  char *avvectors[PostFcn::AVVSIZE];
  char *forces;
  char *tavforces;
  char *hydrostaticforces;
  char *hydrodynamicforces;
  char *lift;
  char *tavlift;
  char *hydrostaticlift;
  char *hydrodynamiclift;
  char *generalizedforces;
  char *residuals;
  char *material_volumes;
  char *conservation;
  char *modeFile;
  Vec3D *TavF, *TavM; 
  Vec3D *TavL;
  VecSet< DistSVec<double,3> > *mX;

	// Gappy POD
  char *newtonresiduals;
	char *jacobiandeltastate;
	char *reducedjac;
	char *staterom;
	char *error;

  double tprevf, tprevl, tinit;
  double tener,tenerold;

  FILE **fpForces;
  FILE **fpLift;
  FILE **fpTavForces;
  FILE **fpTavLift;
  FILE **fpHydroStaticForces;
  FILE **fpHydroDynamicForces;
  FILE **fpHydroStaticLift;
  FILE **fpHydroDynamicLift;
  FILE *fpResiduals;
  FILE *fpMatVolumes;
  FILE *fpConservationErr;
  FILE *fpGnForces;
  FILE *fpStateRom;
  FILE *fpError;


  DistVec<double>    *Qs;
  DistSVec<double,3> *Qv;
  
  DistVec<double>    *AvQs[PostFcn::AVSSIZE];
  DistSVec<double,3> *AvQv[PostFcn::AVVSIZE];

// Included (MB)
  bool switchOpt;

  double dSscale[PostFcn::DSSIZE];
  double dVscale[PostFcn::DVSIZE];

  char *dScalars[PostFcn::DSSIZE];
  char *dVectors[PostFcn::DVSIZE];
  char *dSolutions;
  char *dForces;

  FILE *fpdForces;

  char *heatfluxes;
  FILE **fpHeatFluxes;

  struct {

    double* results;
    int numNodes;
    int* subId;
    int* locNodeId;
    int* last;
    int step;
    std::vector<Vec3D> locations;
  } nodal_output;

  char *nodal_scalars[PostFcn::SSIZE];
  char *nodal_vectors[PostFcn::VSIZE];
  
  
private:
  bool toWrite(int it, bool lastIt, double t);
  int getStep(int it, bool lastIt, double t);

public:

  TsOutput(IoData &, RefVal *, Domain *, PostOperator<dim> *);
  ~TsOutput();

  void updatePrtout(double t);
  void setMeshMotionHandler(IoData&, MeshMotionHandler*);
  FILE *backupAsciiFile(char *);
  void openAsciiFiles();
  void closeAsciiFiles();
  void writeForcesToDisk(bool, int, int, int, double, double, double*, 
			 DistSVec<double,3> &, DistSVec<double,dim> &,
                         DistVec<int> * = 0); 
  void writeForcesToDisk(DistExactRiemannSolver<dim>&, bool, int, int, int, double, double, double*,
                         DistSVec<double,3> &, DistSVec<double,dim> &,
                         DistVec<int> * = 0);
  void writeHydroForcesToDisk(bool, int, int, int, double, double, double*, 
			 DistSVec<double,3> &, DistSVec<double,dim> &,
                         DistVec<int> * = 0);
  void writeLiftsToDisk(IoData &, bool, int, int, int, double, double, double*,
                         DistSVec<double,3> &, DistSVec<double,dim> &,
                         DistVec<int> * = 0);
  void writeHydroLiftsToDisk(IoData &, bool, int, int, int, double, double, double*,
                         DistSVec<double,3> &, DistSVec<double,dim> &,
                         DistVec<int> * = 0);
  void writeHeatFluxesToDisk(bool, int, int, int, double, double,
                             double* , DistSVec<double,3> &, DistSVec<double,dim> &,
                             DistVec<int> * = 0);
  void writeResidualsToDisk(int, double, double, double);
  void writeMaterialVolumesToDisk(int, double, DistVec<double>&, DistVec<int>* = 0);
  void writeStateRomToDisk(int, double, int, const Vec<double> &);
  void writeErrorToDisk(const int, const double, const int, const double *);
  void writeConservationErrors(IoData &iod, int it, double t, int numPhases,
                               double **expected, double **computed);
  void writeDisplacementVectorToDisk(int step, double tag, DistSVec<double,3> &X,
                                     DistSVec<double,dim> &U);
	void writeBinaryVectorsToDiskRom(bool, int, double, DistSVec<double,dim> *,
			DistSVec<double,dim> *, VecSet<DistSVec<double,dim> > *);
  void writeBinaryVectorsToDisk(bool, int, double, DistSVec<double,3> &, 
				DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *);

  void cleanProbesFile();
  
  void writeProbesToDisk(bool, int, double, DistSVec<double,3> &, 
			 DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *);
  
  template<int dimLS>
  void writeBinaryVectorsToDisk(bool, int, double, DistSVec<double,3> &,
                                DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *,
                                DistVec<int> &,DistSVec<double,dimLS>* = NULL);
  
  template<int dimLS>
    void writeProbesToDisk(bool, int, double, DistSVec<double,3> &,
			   DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *,
			   DistVec<int> &,DistSVec<double,dimLS>* = NULL);
  
  void writeBinaryVectorsToDisk(bool, int, double, DistSVec<double,3> &,
                                DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *,
                                DistVec<int> &);

  void writeProbesToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
			 DistVec<double> &A, DistSVec<double,dim> &U, 
			 DistTimeState<dim> *timeState,
			 DistVec<int> &fluidId);
  
  void writeAvgVectorsToDisk(bool,int,double,DistSVec<double,3> &,
                             DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *);

// Included (MB)
  void rstVar(IoData &);
  void writeDerivativeOfForcesToDisk(int, int, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double &, double &);
  void writeBinaryDerivativeOfVectorsToDisk(int, int, double [3], DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> *);

};


template<int dim>
int TsOutput<dim>::counter = 0;

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TsOutput.C>
#endif

#endif
