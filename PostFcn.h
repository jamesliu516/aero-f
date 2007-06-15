#ifndef _POST_FCN_H_
#define _POST_FCN_H_

#include <NavierStokesTerm.h>
#include <SpalartAllmarasTerm.h>
#include <DESTerm.h>
#include <KEpsilonTerm.h>

class VarFcn;
class WallFcn;

struct Vec3D;

//------------------------------------------------------------------------------

class PostFcn {

public:

  enum ScalarType {DENSITY = 0, MACH = 1, PRESSURE = 2, TEMPERATURE = 3, TOTPRESSURE = 4,
		   VORTICITY = 5, NUT_TURB = 6, K_TURB = 7, EPS_TURB = 8, EDDY_VISCOSITY = 9, 
		   DELTA_PLUS = 10, PSENSOR = 11, CSDLES = 12, CSDVMS = 13, MUT_OVER_MU = 14,
                   PHILEVEL = 15, DIFFPRESSURE = 16, SPEED = 17, HYDROSTATICPRESSURE = 18,
                   HYDRODYNAMICPRESSURE = 19, WTMACH = 20, WTSPEED = 21, SSIZE = 22};
  enum VectorType {VELOCITY = 0, DISPLACEMENT = 1, FLIGHTDISPLACEMENT = 2, LOCALFLIGHTDISPLACEMENT = 3, VSIZE = 4};
  enum ScalarAvgType {DENSITYAVG = 0, MACHAVG = 1, PRESSUREAVG = 2, TEMPERATUREAVG = 3,
                     TOTPRESSUREAVG = 4, VORTICITYAVG = 5, AVSSIZE = 6};
  enum VectorAvgType {VELOCITYAVG = 0, DISPLACEMENTAVG = 1, AVVSIZE = 2};

protected:

  VarFcn *varFcn;

public:

  PostFcn(VarFcn *);
  ~PostFcn() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, double = 0);
  virtual double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
					   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
		double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,double *nodalForceWeight, int = 0) = 0;
  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]) = 0;
  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double) = 0;
  virtual bool doesFaceNeedGradientP1Function() { return false; }

  virtual double* getMeshVel()  { return varFcn->getMeshVel(); }
  
};

//------------------------------------------------------------------------------

class PostFcnEuler : public PostFcn {

protected:

  static const double third;

  double mach;
  double pinfty;
  double depth;
  double gravity;
  double alpha;
  double beta;
  double nGravity[3];
  double refRho;
  double refPressure;
  double refVel;

public:

  PostFcnEuler(IoData &, VarFcn *);
  ~PostFcnEuler() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, double = 0);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double *nodalForceWeight, int = 0);
  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]);
  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double);


};

//------------------------------------------------------------------------------

class PostFcnNS : public PostFcnEuler, public NavierStokesTerm {

protected:

  WallFcn* wallFcn;

private:

  Vec3D computeViscousForce(double [4][3], Vec3D&, double [3], double*, double* [3], double* [4]);

public:

  PostFcnNS(IoData &, VarFcn *);
  ~PostFcnNS();

  double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
				   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double *nodalForceWeight, int = 0);
  double computeHeatPower(double [4][3], Vec3D&, double [3],
			  double*, double* [3], double* [4]);
  double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
			      double* [3], double* [4], double);
  bool doesFaceNeedGradientP1Function() { return ((wallFcn) ? false : true); }
  
};

//------------------------------------------------------------------------------

class PostFcnSA : public PostFcnNS, public SATerm {

public:

  PostFcnSA(IoData &, VarFcn *);
  ~PostFcnSA() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
};

//------------------------------------------------------------------------------

class PostFcnDES : public PostFcnNS, public DESTerm {

public:

  PostFcnDES(IoData &, VarFcn *);
  ~PostFcnDES() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
};

//------------------------------------------------------------------------------

class PostFcnKE : public PostFcnNS, public KEpsilonTerm {

public:

  PostFcnKE(IoData &, VarFcn *);
  ~PostFcnKE() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
};

//------------------------------------------------------------------------------

#endif
