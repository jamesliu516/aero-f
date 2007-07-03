#ifndef _WALL_FCN_H_
#define _WALL_FCN_H_

class IoData;
class VarFcn;
class ViscoFcn;
struct Vec3D;

//------------------------------------------------------------------------------

class WallFcn {

  const static double third;
  const static double eleventh;

  double gam;
  double prandtl;

  VarFcn *varFcn;
  ViscoFcn *viscoFcn;

protected:

  double vkcst;
  double reynolds;

  virtual void computeWallValues(double utau, double delta, double rho,
				 double ut, double mu, double *V) {}

private:

  Vec3D computeTangentVector(Vec3D &, Vec3D &);
  void computeFaceValues(double [3], double *, double *[3], double &, Vec3D &, 
			 double &, Vec3D &, double &, double &);
  double computeFrictionVelocity(Vec3D &, double, double, Vec3D &, double);
  double computeFrictionTemperature(double, double, double, double, double);

public:

  WallFcn(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcn() {}

  void computeSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  Vec3D computeForce(Vec3D &, double [3], double *, double *[3]);
  double computeInterfaceWork(Vec3D&, double [3], double*, double* [3]);
  double computeHeatPower(Vec3D &, double [3], double *, double *[3]);
  double computeDeltaPlus(Vec3D &, double [3], double *, double *[3]);

  template<int neq>
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, 
				  double *[3], double (*)[neq][neq]);

};

//------------------------------------------------------------------------------
/*
@ARTICLE{reichardt-42,
  author = "Reichardt, H.",
  title = "Gesetzm{\"a}ssigkeiten der freien {T}urbulenz",
  journal = "VDI-Forschungsheft 414, 1st edition, Berlin",
  year = 1942,
} 
*/

class WallFcnSA : public WallFcn {

  double cv1_pow3;

public:

  WallFcnSA(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcnSA() {}

  void computeWallValues(double, double, double, double, double, double *);

};

//------------------------------------------------------------------------------

class WallFcnKE : public WallFcn {

  double orcmu;

public:

  WallFcnKE(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcnKE() {}

  void computeWallValues(double, double, double, double, double, double *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <WallFcn.C>
#endif

#endif
