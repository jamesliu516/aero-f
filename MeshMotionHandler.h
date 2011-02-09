#ifndef _MESH_MOTION_HANDLER_H_
#define _MESH_MOTION_HANDLER_H_

#include <IoData.h>
#include <DistVector.h>
#include <Vector3D.h>
#include <LevelSet/LevelSetStructure.h>
#include <FSI/DynamicNodalTransfer.h>

class VarFcn;
class MatchNodeSet;
class Domain;
class StructExc;
class CorotSolver;
class MeshMotionSolver;
class MemoryPool;

template<int dim> class PostOperator;

//------------------------------------------------------------------------------

class MeshMotionHandler {

protected:

  double tscale;
  double oolscale;
  double Wn;

  DistSVec<double,3> F;
  DistSVec<double,3> F0;
  DistVec<double> Pin;
  DistSVec<double,3> X0;
  DistSVec<double,3> dX;

  Domain *domain;
  Communicator *com;

public:

  MeshMotionHandler(IoData &, Domain *);
  virtual ~MeshMotionHandler() {}

  virtual double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * =0) {return 0.0;} 
  virtual double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &) {return 0.0;}

  template<int dim>
  void computeInterfaceWork(double, PostOperator<dim>*, DistSVec<double,3>&, 
			    DistSVec<double,dim>&, DistSVec<double,3>&,
			    DistSVec<double,dim>&, double*,
			    DistVec<int> * = 0, DistVec<int> * = 0);

  virtual int getAlgNum()  { return 0; }

};

//------------------------------------------------------------------------------

class RigidMeshMotionHandler {

protected: 
  RigidMeshMotionData::Tag typeTag;
  RigidMeshMotionData::LawType typeLaw;
  double tscale;

  Vec3D vin;
  double cin;

  Vec3D v;
  Vec3D a;
  Vec3D  deltaX;

  int nvpts;
  double *timeVpts;
  Vec3D *velVpts;

  DistSVec<double,3> Xr;

public:

  RigidMeshMotionHandler(IoData &, VarFcn *, double *, Domain *);
  ~RigidMeshMotionHandler();

  void setupVelocityPoints(IoData &);
  void addRigidMotion(double, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &, Communicator *com = 0);

  DistSVec<double,3> &getRelativePositionVector(double, DistSVec<double,3> &);
  DistSVec<double,3> &getFlightPositionVector(double, DistSVec<double,3> &);

  const char* getTagName();
  double getTagValue(double);
  void updateMomentArm(Vec3D &x0);

};

//------------------------------------------------------------------------------

class AccMeshMotionHandler : public MeshMotionHandler, 
			     public RigidMeshMotionHandler {

  double dt;

public:

  AccMeshMotionHandler(IoData &, VarFcn *, double *, Domain *);
  ~AccMeshMotionHandler() {}

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * =0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

class AeroMeshMotionHandler : public MeshMotionHandler {

protected:


  AeroelasticData::Force forceComputation;
  enum TimeIntegrator {IMPLICIT_FIRST_ORDER, IMPLICIT_SECOND_ORDER} timeIntegrator;
  bool steady;
  int it0;
  double mppFactor;

  DistSVec<double,3>* Fn;
  DistSVec<double,3>* Fnp1;
  DistSVec<double,3>* Favg;

  StructExc *strExc;

  MeshMotionSolver *mms;
  char *posFile;

public:

  AeroMeshMotionHandler(IoData &, VarFcn *, double *, MatchNodeSet **, 
			Domain *, MemoryPool *);
  ~AeroMeshMotionHandler();

  template<int dim>
  void setup(int *, double *, PostOperator<dim>*, DistSVec<double,3>&, 
             DistSVec<double,dim>&, DistVec<int> * = 0);

  template<int dim>
  void resetOutputToStructure(PostOperator<dim>*, DistSVec<double,3>&,
                              DistSVec<double,dim>&, DistVec<int> * = 0);

  template<int dim>
  void updateOutputToStructure(double, double, PostOperator<dim>*,
                               DistSVec<double,3>&, DistSVec<double,dim>&,
                               DistVec<int> * = 0);

  int getModalMotion(DistSVec<double,3> &);

  virtual double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * =0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  int getAlgNum(); 

};

//------------------------------------------------------------------------------

class AccAeroMeshMotionHandler : public AeroMeshMotionHandler, 
				 public RigidMeshMotionHandler {

public:

  AccAeroMeshMotionHandler(IoData &, VarFcn *, double *, MatchNodeSet **, 
			   Domain *, MemoryPool *);
  ~AccAeroMeshMotionHandler() {}

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * =0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

class DeformingMeshMotionHandler : public MeshMotionHandler {

  double dt;
  double omega;

  DistSVec<double,3> *dXmax;

  MeshMotionSolver *mms;

public:

  DeformingMeshMotionHandler(IoData &, Domain *);
  ~DeformingMeshMotionHandler();

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * = 0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

  DistSVec<double,3> getModes();

  void setup(DistSVec<double, 3> &X);

};

//------------------------------------------------------------------------------

class PitchingMeshMotionHandler : public MeshMotionHandler {

  double alpha_in;
  double alpha_max;

  double dt;
  double omega;

  double x1[3];
  double x2[3];

  double u, v, w;
  double ix, iy, iz;

  MeshMotionSolver *mms;

public:

  PitchingMeshMotionHandler(IoData &, Domain *);
  ~PitchingMeshMotionHandler();

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * = 0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

  DistSVec<double,3> getModes();

  void setup(DistSVec<double, 3> &X);

};

//------------------------------------------------------------------------------

class HeavingMeshMotionHandler : public MeshMotionHandler {

  double dt;
  double omega;

  double delta[3];

  MeshMotionSolver *mms;

public:

  HeavingMeshMotionHandler(IoData &, Domain *);
  ~HeavingMeshMotionHandler();

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * = 0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

  DistSVec<double,3> getModes();

  void setup(DistSVec<double, 3> &X);

};

//------------------------------------------------------------------------------

class AccForcedMeshMotionHandler : public DeformingMeshMotionHandler, 
				   public RigidMeshMotionHandler {

public:

  AccForcedMeshMotionHandler(IoData &, VarFcn *, double *, Domain *);
  ~AccForcedMeshMotionHandler() {}

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * = 0);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

class RigidRollMeshMotionHandler : public MeshMotionHandler {

  double alpha_in;
  double beta_in;

  double dt;
  double maxtime;

  Vec3D xo;
  Vec3D u;

  MeshMotionSolver *mms;
  
private:

  double computeRotationAngle(double);

public:

  RigidRollMeshMotionHandler(IoData &, double *, Domain *);
  ~RigidRollMeshMotionHandler() {}

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

class EmbeddedMeshMotionHandler : public MeshMotionHandler {  //<! For embedded fluid-structure interactions
 
protected:

  double dts;            //<! structure time-step.
  //Vec3D *structX0;       //<! initial position of structure nodes.
  //Vec3D *structXn;       //<! position of struct nodes at t^n.
  //Vec3D *structXnPlus1;  //<! position of struct nodes at t^{n+1}.
  //Vec3D *structVel;      //<! velocity of struct nodes.

  int it0;            //<! restart timestep
  //int structVelocity; //<! 0: use the velocity received by structure;  1: compute by finite-difference.
  
  DynamicNodalTransfer *dynNodalTransfer; 
  DistLevelSetStructure *distLSS; //<! interface finder (not necessarily a levelset solver).

public:

  EmbeddedMeshMotionHandler(IoData &, Domain *, DynamicNodalTransfer *, DistLevelSetStructure *);
  ~EmbeddedMeshMotionHandler();

  void setup(double *);
  double update(bool *lastIt, int it, double t, DistSVec<double,3> &Xdot, DistSVec<double,3> &X) {return dts;}
  double updateStep1(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &, double * = 0); 
  void step1ForA6(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &); 
  void step1ForC0FEM(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  void step1ForC0XFEM(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  void step1ForC0XFEM3D(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &); 
  void step2ForA6(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &); 
  void step2ForC0(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &); 
  void step2ForC0XFEM3D(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &); 
  /** get displacement and pass it to distLSS.  send force to structure. */

  int getAlgNum()  { return 0; }

};

//------------------------------------------------------------------------------

class RbmExtractor : public MeshMotionHandler {

  char *name_in;
  char *name_out1;
  char *name_out2;

  CorotSolver *cs;

public:

  RbmExtractor(IoData &, Domain *);
  ~RbmExtractor(){}

  double update(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);
  double updateStep2(bool *, int, double, DistSVec<double,3> &, DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MeshMotionHandler.C>
#endif

#endif
