#ifndef _MESH_MOTION_SOLVER_H_
#define _MESH_MOTION_SOLVER_H_

#include <IoData.h>
#include <Domain.h>
#include <DistVector.h>
#include <KspPrec.h>
#include <BCApplier.h>
#include <TsParameters.h>
#include <ErrorHandler.h>

class MatchNodeSet;
class CorotSolver;
class Communicator;
class MemoryPool;
class Timer;

template<class Scalar, int dim> class StiffMat;
template<class ProbDesc> class NewtonSolver;
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------------

class MeshMotionSolver {

public:

  MeshMotionSolver() {};
  virtual ~MeshMotionSolver() {}

  virtual int solve(DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual int solveSensitivity(DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual void setup(DistSVec<double,3> &) = 0;
  virtual void applyProjector(DistSVec<double,3> &dX) {};

};

//------------------------------------------------------------------------------

class TetMeshMotionSolver : public MeshMotionSolver {

public:

  typedef DistSVec<double,3> SolVecType;
  typedef DistSVec<double,1> PhiVecType;
  typedef DistVec<double> VolVecType;

protected:

  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  DefoMeshMotionData::Element typeElement;

  DistSVec<double,3> *F0;
  DistSVec<double,3> *dX0;

  CorotSolver *cs;
  StiffMat<double,3> *mvp;
  KspPrec<3> *pc;
  KspSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator> *ksp;
  NewtonSolver<TetMeshMotionSolver> *ns;

  Domain *domain;

  Communicator *com;
  Timer *timer;

  double volStiff;

  BCApplier* meshMotionBCs; //HB

public:

  TetMeshMotionSolver(DefoMeshMotionData &, MatchNodeSet **, Domain *, MemoryPool *);
  TetMeshMotionSolver(Domain *dom) : domain(dom) {};
  virtual ~TetMeshMotionSolver();

  virtual int solve(DistSVec<double,3> &, DistSVec<double,3> &);
  int solveSensitivity(DistSVec<double,3> &, DistSVec<double,3> &) { return 0.0; };

  void applyProjector(DistSVec<double,3> &X);
 
  void printf(int, const char *, ...);
  virtual void computeFunction(int, DistSVec<double,3> &, DistSVec<double,3> &);
  void recomputeFunction(DistSVec<double,3> &, DistSVec<double,3> &) {}
  double recomputeResidual(DistSVec<double,3> &, DistSVec<double,3> &) { return 0.0; }
  void computeJacobian(int, DistSVec<double,3> &, DistSVec<double,3> &);
  void setOperators(DistSVec<double,3> &);
  int solveLinearSystem(int, DistSVec<double,3> &, DistSVec<double,3> &);

  int checkSolution(DistSVec<double,3> &X) { return 0; }
  int checkFailSafe(DistSVec<double,3>& X) { return 0; }
  void resetFixesTag() { return;}
  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }
  DistInfo &getVecInfo() const { return domain->getNodeDistInfo(); }
  
  virtual void setup(DistSVec<double,3> &X);
  TsParameters* getTsParams() { return NULL; }
  ErrorHandler* getErrorHandler() {return NULL; }

  // Included (MB)
  int fixSolution(DistSVec<double,3> &X, DistSVec<double,3> &dX) { return 0; }
	void writeBinaryVectorsToDiskRom(bool lastIt, int it, double t,
			DistSVec<double,3> *F1 = NULL, DistSVec<double,3> *F2 = NULL,
			VecSet< DistSVec<double,3> > *F3 = NULL) {};

  void printNodalDebug(int globNodeId, int identifier, DistSVec<double,3> *U, DistVec<int> *Id=0, DistVec<int> *Id0=0) {}

};

//------------------------------------------------------------------------------
//
class EmbeddedALETetMeshMotionSolver : public TetMeshMotionSolver {

public:

  EmbeddedALETetMeshMotionSolver(DefoMeshMotionData &, MatchNodeSet **, Domain *, MemoryPool *);
  ~EmbeddedALETetMeshMotionSolver(){};
  int solve(DistSVec<double,3> &, DistSVec<double,3> &);
};

//------------------------------------------------------------------------------

class TetMeshSensitivitySolver : public TetMeshMotionSolver {

protected:

 DistSVec<double,3> *currentPosition;
 NewtonSolver<TetMeshSensitivitySolver> *nss;

public:

  TetMeshSensitivitySolver(DefoMeshMotionData &, MatchNodeSet **, Domain *, MemoryPool *);
  ~TetMeshSensitivitySolver();
  int solveSensitivity(DistSVec<double,3> &, DistSVec<double,3> &);

  void computeFunction(int, DistSVec<double,3> &, DistSVec<double,3> &);
  void setup(DistSVec<double,3> &X);
};

//------------------------------------------------------------------------------

#endif
