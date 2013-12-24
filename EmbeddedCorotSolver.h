#ifndef _EMBCOROT_SOLVER_H_
#define _EMBCOROT_SOLVER_H_

#include <IoData.h>
#include <DistVector.h>

class MatchNodeSet;
class Domain;
class Communicator;
class BCApplier;

//------------------------------------------------------------------------------

class EmbeddedCorotSolver {

  int numStNodes;

  int numLocSub;

  double *Xs0;

  DistSVec<double,3> X0;

  double cg0[3], cgN[3];

  double n[3];

  double R[3][3];

  Domain *domain;
  Communicator *com;

  BCApplier* meshMotionBCs;

  enum SymmetryAxis {NONE, AXIS_X, AXIS_Y, AXIS_Z} SymAxis;
  enum Type {BASIC, COROTATIONAL} type;

private:

  void computeCG(double *, double [3]);
  void computeRotGradAndJac(double *, double [3][3], 
			    double [3], double [3], double [3][3]);
  void computeRotMat(double *, double [3][3]);
  void rotLocVec(double [3][3], double [3]);
  void solveDeltaRot(double *, double[3][3], double [3]);
  void solveRotMat(double [3][3], double [3]);
  void computeNodeRot(double [3][3], DistSVec<double,3> &, double [3], double [3]);
  void computeDeltaNodeRot(double [3][3], DistSVec<double,3> &, DistSVec<double,3> &, double [3], double [3]);

  void printRotMat(double mat[3][3]);

public:

  EmbeddedCorotSolver(DefoMeshMotionData &, Domain *, double*, int);
  ~EmbeddedCorotSolver() {};

  void applyProjector(DistSVec<double,3> &Xdot);
  void findProjection(DistSVec<double,3> &dX);
  void computeMeanDXForSlidingPlane(DistSVec<double,3> &dX, double meandX[3]);

  void setup(double *Xtilde, int nNodes);
  void solve(double *Xtilde, int nNodes, DistSVec<double,3> &X, DistSVec<double,3> &dX);
  // setup computes the rotation and CG that will fit best for restarting
		
};

//------------------------------------------------------------------------------

#endif
