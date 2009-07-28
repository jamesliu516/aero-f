#ifndef _POST_FCN_H_
#define _POST_FCN_H_

#include <NavierStokesTerm.h>
#include <SpalartAllmarasTerm.h>
#include <DESTerm.h>
#include <KEpsilonTerm.h>

class VarFcn;
class WallFcn;

// Included (MB)
class Communicator;

struct Vec3D;

//------------------------------------------------------------------------------

class PostFcn {

public:

// Included (MB)
  enum ScalarType {DENSITY = 0, MACH = 1, PRESSURE = 2, TEMPERATURE = 3, TOTPRESSURE = 4,
		   VORTICITY = 5, CSDLES = 6, CSDVMS = 7, SKIN_FRICTION = 8, NUT_TURB = 9, 
                   K_TURB = 10, EPS_TURB = 11, EDDY_VISCOSITY = 12, DELTA_PLUS = 13, 
                   PSENSOR = 14, MUT_OVER_MU = 15, PHILEVEL = 16, DIFFPRESSURE = 17, 
                   SPEED = 18, HYDROSTATICPRESSURE = 19, HYDRODYNAMICPRESSURE = 20, 
                   WTMACH = 21, WTSPEED = 22, VELOCITY_NORM = 23, TEMPERATURE_NORMAL_DERIVATIVE = 24, 
                   SURFACE_HEAT_FLUX = 25, PRESSURECOEFFICIENT = 26, PHILEVEL_STRUCTURE = 27,
                   SSIZE = 28};


// Original
/*
  enum ScalarType {DENSITY = 0, MACH = 1, PRESSURE = 2, TEMPERATURE = 3, TOTPRESSURE = 4,
		   VORTICITY = 5, NUT_TURB = 6, K_TURB = 7, EPS_TURB = 8, EDDY_VISCOSITY = 9, 
		   DELTA_PLUS = 10, PSENSOR = 11, CSDLES = 12, CSDVMS = 13, MUT_OVER_MU = 14,
                   PHILEVEL = 15, DIFFPRESSURE = 16, SPEED = 17, HYDROSTATICPRESSURE = 18,
                   HYDRODYNAMICPRESSURE = 19, WTMACH = 20, WTSPEED = 21, SSIZE = 22};
*/
  enum VectorType {VELOCITY = 0, DISPLACEMENT = 1, FLIGHTDISPLACEMENT = 2, LOCALFLIGHTDISPLACEMENT = 3, VSIZE = 4};
  enum ScalarAvgType {DENSITYAVG = 0, MACHAVG = 1, PRESSUREAVG = 2, TEMPERATUREAVG = 3,
                      TOTPRESSUREAVG = 4, VORTICITYAVG = 5, CSDLESAVG = 6, CSDVMSAVG = 7, 
                      SKIN_FRICTIONAVG =8, AVSSIZE = 9};
  enum VectorAvgType {VELOCITYAVG = 0, DISPLACEMENTAVG = 1, AVVSIZE = 2};

// Included (MB)
  enum ScalarDerivativeType {DERIVATIVE_DENSITY = 0, DERIVATIVE_MACH = 1, DERIVATIVE_PRESSURE = 2, 
                             DERIVATIVE_TEMPERATURE = 3, DERIVATIVE_TOTPRESSURE = 4, DERIVATIVE_NUT_TURB = 5, 
                             DERIVATIVE_EDDY_VISCOSITY = 6, DERIVATIVE_VELOCITY_SCALAR = 7, DSSIZE = 8};
  enum VectorDerivativeType {DERIVATIVE_VELOCITY_VECTOR = 0, DERIVATIVE_DISPLACEMENT = 1, DVSIZE = 2};

protected:

  VarFcn *varFcn;

public:

  PostFcn(VarFcn *);
  ~PostFcn() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, double = 0);
  virtual double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
					   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
		double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0) = 0;
  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3], 
	double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0) = 0;

  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]) = 0;

  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool) =0;

  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double) = 0;
  virtual bool doesFaceNeedGradientP1Function() { return false; }

  virtual double* getMeshVel()  { return varFcn->getMeshVel(); }
  
// Included (MB)
  virtual double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);
  virtual void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0) = 0;
  virtual void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0) = 0;
  virtual void rstVar(IoData &, Communicator*) = 0;
  virtual double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]) = 0;

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

// Included (MB)
  bool dimFlag;
  double dpinfty;
  double dPin;

public:

  PostFcnEuler(IoData &, VarFcn *);
  ~PostFcnEuler() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, double = 0);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0);
  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0);
  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]);
  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool);
  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double);


// Included (MB)
  virtual double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);
  virtual void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  virtual void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  void rstVar(IoData &, Communicator*);
  virtual double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);

};

//------------------------------------------------------------------------------

class PostFcnNS : public PostFcnEuler, public NavierStokesTerm {

protected:

  WallFcn* wallFcn;

private:

  Vec3D computeViscousForce(double [4][3], Vec3D&, double [3], double*, double* [3], double* [4]);

// Included (MB)
  Vec3D computeDerivativeOfViscousForce(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);

public:

  PostFcnNS(IoData &, VarFcn *);
  ~PostFcnNS();

  double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
				   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0);
  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0);
  double computeHeatPower(double [4][3], Vec3D&, double [3],
			  double*, double* [3], double* [4]);
  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool);
  double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
			      double* [3], double* [4], double);
  bool doesFaceNeedGradientP1Function() { return ((wallFcn) ? false : true); }
  
// Included (MB)
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

  void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);

};

//------------------------------------------------------------------------------

class PostFcnSA : public PostFcnNS, public SATerm {

public:

  PostFcnSA(IoData &, VarFcn *);
  ~PostFcnSA() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

class PostFcnDES : public PostFcnNS, public DESTerm {

public:

  PostFcnDES(IoData &, VarFcn *);
  ~PostFcnDES() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

class PostFcnKE : public PostFcnNS, public KEpsilonTerm {

public:

  PostFcnKE(IoData &, VarFcn *);
  ~PostFcnKE() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, double);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

#endif
