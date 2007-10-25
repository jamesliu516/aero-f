#ifndef _FEM_EQUATION_TERM_DESC_H_
#define _FEM_EQUATION_TERM_DESC_H_

#include <FemEquationTerm.h>
#include <NavierStokesTerm.h>
#include <SpalartAllmarasTerm.h>
#include <KEpsilonTerm.h>
#include <DESTerm.h>
#include <Vector3D.h>

#include <stdlib.h>
#include <stdio.h>

class IoData;

//------------------------------------------------------------------------------

class FemEquationTermNS : public FemEquationTerm, public NavierStokesTerm {

public:

  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermNS(IoData &, VarFcn *);
  ~FemEquationTermNS() {}

  double computeViscousTimeStep(double *, double *);

  bool computeVolumeTerm(double [4][3], double [4], double *[4],
			 double *, double *, double *, double, 
                         SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double,
				 SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
			  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
			  double *, double *[4], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
				  double *, double *[4], double *);

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double [4][3], double [4][3], double [4], double *[4], double *[4],
			 double, double *, double *, double *, double, SVec<double,3> &, int [4], int);

  void computeDerivativeOfSurfaceTerm(int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double, double *);

  void computeDerivativeOfSurfaceTerm(double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
			  double *, double *, double *[4], double *[4], double, double *);

  void rstVar(IoData &ioData, Communicator *com) { rstVarNS(ioData, com); if (wallFcn) wallFcn->rstVar(ioData, com);}

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {}

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermSA : public FemEquationTerm, public NavierStokesTerm, 
			  public SATerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermSA(IoData &, VarFcn *);
  ~FemEquationTermSA() {}

  double computeViscousTimeStep(double *, double *);

  bool computeVolumeTerm(double [4][3], double [4], double *[4],
			 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *, double *, double,
                                 SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
			  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
			  double *, double *[4], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
				  double *, double *[4], double *);
  bool doesSourceTermExist() { return true; }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double [4][3], double [4][3], double [4], double *[4], double *[4],
			 double, double *, double *, double *, double, SVec<double,3> &, int [4], int);

  void computeDerivativeOfSurfaceTerm(int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double, double *);

  void computeDerivativeOfSurfaceTerm(double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
			  double *, double *, double *[4], double *[4], double, double *);

  void rstVar(IoData &ioData, Communicator *com) { rstVarNS(ioData, com); if (wallFcn) wallFcn->rstVar(ioData, com);}

  void computeBCsJacobianWallValues(int, Vec3D &, double [3], double *, double *, double *[3]);

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermDES : public FemEquationTerm, public NavierStokesTerm, 
			  public DESTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  double cdes;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermDES(IoData &, VarFcn *);
  ~FemEquationTermDES() {}

  double computeViscousTimeStep(double *, double *);

  double max(double a, double b) { return (a>b) ? a : b; }
  double min(double a, double b) { return (a<b) ? a : b; }

  bool computeVolumeTerm(double [4][3], double [4], double *[4],
			 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *, double *, double,
                                 SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
			  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
			  double *, double *[4], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
				  double *, double *[4], double *);
  bool doesSourceTermExist() { return true; }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double [4][3], double [4][3], double [4], double *[4], double *[4],
			 double, double *, double *, double *, double, SVec<double,3> &, int [4], int);

  void computeDerivativeOfSurfaceTerm(int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double, double *);

  void computeDerivativeOfSurfaceTerm(double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
			  double *, double *, double *[4], double *[4], double, double *);

  void rstVar(IoData &ioData, Communicator *com) { rstVarNS(ioData, com); if (wallFcn) wallFcn->rstVar(ioData, com);}

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermSAmean : public FemEquationTerm, public NavierStokesTerm, 
			      public SATerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermSAmean(IoData &, VarFcn *);
  ~FemEquationTermSAmean() {}

  double computeViscousTimeStep(double *, double *);

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
				  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
				  double *, double *[4], double *);

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------
class FemEquationTermDESmean : public FemEquationTerm, public NavierStokesTerm, 
			      public DESTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  ////map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermDESmean(IoData &, VarFcn *);
  ~FemEquationTermDESmean() {}

  double computeViscousTimeStep(double *, double *);

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
				  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
				  double *, double *[4], double *);

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermSAturb : public FemEquationTerm, public NavierStokesTerm,
			      public SATerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;

  FemEquationTermSAturb(IoData &, VarFcn *);
  ~FemEquationTermSAturb() {}

  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemSAturb\n");
    exit(1);
  }

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
				  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double){
    fprintf(stderr, "*** Error: computeViscousDerivativeOfTimeStep should not be called in FemSAturb\n");
    exit(1);
  }

};

//------------------------------------------------------------------------------
class FemEquationTermDESturb : public FemEquationTerm, public NavierStokesTerm,
			      public DESTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;

  FemEquationTermDESturb(IoData &, VarFcn *);
  ~FemEquationTermDESturb() {}

  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemDESturb\n");
    exit(1);
  }

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, 
                         double, SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
				  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double){
    fprintf(stderr, "*** Error: computeDerivativeOfViscousTimeStep should not be called in FemDESturb\n");
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FemEquationTermKE : public FemEquationTerm, public NavierStokesTerm, 
			  public KEpsilonTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermKE(IoData &, VarFcn *);
  ~FemEquationTermKE() {}

  double computeViscousTimeStep(double *, double *);

  bool computeVolumeTerm(double [4][3], double [4], double *[4],
			 double *, double *, double *, double, 
                         SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *,
                                 double *, double, SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
			  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  bool doesSourceTermExist() { return true; }

  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double [4][3], double [4][3], double [4], double *[4], double *[4],
			 double, double *, double *, double *, double, SVec<double,3> &, int [4], int);

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &, Vec3D &, double [3],
			  double *, double *, double *[3], double *[3], double, double *);

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) { rstVarNS(ioData, com); if (wallFcn) wallFcn->rstVar(ioData, com);}

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermKEmean : public FemEquationTerm, public NavierStokesTerm, 
			      public KEpsilonTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;

  FemEquationTermKEmean(IoData &, VarFcn *);
  ~FemEquationTermKEmean() {}

  double computeViscousTimeStep(double *, double *);

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
				  double *, double *[3], double *);

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
					  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }

  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double);

};

//------------------------------------------------------------------------------

class FemEquationTermKEturb : public FemEquationTerm, public NavierStokesTerm,
			      public KEpsilonTerm {

public:

  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;

  FemEquationTermKEturb(IoData &, VarFcn *);
  ~FemEquationTermKEturb() {}

  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemKEturb\n");
    exit(1);
  }

  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
				 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }

  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
			 double *r, double *s, double *, double , 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
			  double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
				  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
			  double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }

// Included (MB)
  bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4], double *dv[4],
			 double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid) {
    fprintf(stderr, "*** Error: computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
			  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
			  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) {
    fprintf(stderr, "*** Error: computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }

  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: rstVar should not be called\n");
    exit(1);
  }

  void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double){
    fprintf(stderr, "*** Error: computeDerivativeOfViscousTimeStep should not be called in FemKEturb\n");
    exit(1);
  }

};

//------------------------------------------------------------------------------

#endif
