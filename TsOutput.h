#ifndef _TS_OUTPUT_H_
#define _TS_OUTPUT_H_

#include <PostFcn.h>
#include <Vector3D.h>

#include <stdio.h>

class IoData;
class RefVal;
class Domain;
class RigidMeshMotionHandler;
class AccMeshMotionHandler;
class LandingMeshMotionHandler;
class Communicator;
class LevelSet;

template<int dim> class PostOperator;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistTimeState;

//------------------------------------------------------------------------------

template<int dim>
class TsOutput {

private:

  RefVal *refVal;
  PostOperator<dim> *postOp;
  RigidMeshMotionHandler *rmmh;
  Domain *domain;
  Communicator *com;

private:

  bool steady;
  int it0;
  int frequency;
  double length;
  double surface;
  static int counter;
  Vec3D x0;

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
  char *modeFile;
  Vec3D *TavF, *TavM; 
  Vec3D *TavL;
  VecSet< DistSVec<double,3> > *mX;

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
  FILE *fpGnForces;


  DistVec<double> *Qs;
  DistSVec<double,3> *Qv;
  
// Included (MB)
  bool switchOpt;

  double dSscale[PostFcn::DSSIZE];
  double dVscale[PostFcn::DVSIZE];

  char *dScalars[PostFcn::DSSIZE];
  char *dVectors[PostFcn::DVSIZE];
  char *dSolutions;
  char *dForces;

  FILE *fpdForces;

public:

  TsOutput(IoData &, RefVal *, Domain *, PostOperator<dim> *);
  ~TsOutput();

  void setMeshMotionHandler(RigidMeshMotionHandler *);
  FILE *backupAsciiFile(char *);
  void openAsciiFiles();
  void closeAsciiFiles();
  void writeForcesToDisk(bool, int, int, int, double, double, double*, 
			 DistSVec<double,3> &, DistSVec<double,dim> &);
  void writeHydroForcesToDisk(bool, int, int, int, double, double, double*, 
			 DistSVec<double,3> &, DistSVec<double,dim> &);
  void writeLiftsToDisk(IoData &, bool, int, int, int, double, double, double*,
                         DistSVec<double,3> &, DistSVec<double,dim> &);
  void writeHydroLiftsToDisk(IoData &, bool, int, int, int, double, double, double*,
                         DistSVec<double,3> &, DistSVec<double,dim> &);
  void writeResidualsToDisk(int, double, double, double);
  void writeBinaryVectorsToDisk(bool, int, double, DistSVec<double,3> &, 
				DistVec<double> &, DistSVec<double,dim> &, DistTimeState<dim> *);
  void writeBinaryVectorsToDisk(bool, int, double, DistSVec<double,3> &,
                                DistVec<double> &, DistSVec<double,dim> &,
                                DistVec<double> &);
  void writeBinaryVectorsToDiskLS(bool, int, double, DistSVec<double,3> &,
                                DistVec<double> &, DistSVec<double,dim> &, LevelSet *);
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
