#include <MeshMotionSolver.h>

#include <MatchNode.h>
#include <StiffMatrix.h>
#include <CorotSolver.h>
#include <KspSolver.h>
#include <NewtonSolver.h>
#include <MeshMotionHandler.h>
#include <Communicator.h>
#include <MemoryPool.h>
#include <Timer.h>
#include <BCApplier.h> 

#include <cstdio>

#ifdef TYPE_PREC_MESH
#define PrecScalar TYPE_PREC_MESH
#else
#define PrecScalar double
#endif

//#define HB_MESHMOTION_DEBUG

//------------------------------------------------------------------------------

TetMeshMotionSolver::TetMeshMotionSolver
(
  DefoMeshMotionData &data, MatchNodeSet **matchNodes, 
  Domain *dom, MemoryPool *mp
) 
: domain(dom)
{

  com = domain->getCommunicator();

  KspData &kspData = data.newton.ksp;
  PcData &pcData = data.newton.ksp.pc;

  typeElement = data.element;
  maxItsNewton = data.newton.maxIts;
  if (data.element == DefoMeshMotionData::TORSIONAL_SPRINGS || data.element == DefoMeshMotionData::BALL_VERTEX)
    maxItsNewton = 1;
  epsNewton = data.newton.eps;
  epsAbsResNewton = data.newton.epsAbsRes;
  epsAbsIncNewton = data.newton.epsAbsInc;

  timer = domain->getTimer();

  F0 = new DistSVec<double,3>(domain->getNodeDistInfo());

  if (data.type == DefoMeshMotionData::COROTATIONAL)
    cs = new CorotSolver(data, matchNodes, domain);
  else
    cs = 0;

  //int **ndType = domain->getNodeType();
  int **ndType = 0;

  meshMotionBCs = domain->getMeshMotionBCs(); //HB

  if (meshMotionBCs)   {
    meshMotionBCs->setDofType(matchNodes);

  }

  mvp = new StiffMat<double,3>(domain, ndType, mp, meshMotionBCs);

  if (pcData.type == PcData::IDENTITY)
    pc = new IdentityPrec<3>(meshMotionBCs);
  else if (pcData.type == PcData::JACOBI)
    //pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DIAGONAL, domain, ndType, meshMotionBCs);
    pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DENSE, domain, ndType, meshMotionBCs);

  else if (pcData.type == PcData::AS || pcData.type == PcData::RAS || pcData.type == PcData::ASH || pcData.type == PcData::AAS)
    pc = new IluPrec<PrecScalar,3>(pcData, domain, ndType);

  if (kspData.type == KspData::RICHARDSON)
    ksp = new RichardsonSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::CG)
    ksp = new CgSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::GMRES)
    ksp = new GmresSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);

  ns = new NewtonSolver<TetMeshMotionSolver>(this);

  volStiff = data.volStiff;

}  

//------------------------------------------------------------------------------

TetMeshMotionSolver::~TetMeshMotionSolver()
{
  if (F0) delete F0;

  if (cs) delete cs;
  if (ns) delete ns;
  if (mvp) delete mvp;
  if (pc) delete pc;
}  

//------------------------------------------------------------------------------
//HB: X <- P.X where P is the projector onto the sliding type of constraints
void
TetMeshMotionSolver::applyProjector(DistSVec<double,3> &X)
{
   if(meshMotionBCs) meshMotionBCs->applyP(X);
}

//------------------------------------------------------------------------------

/*
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/
int TetMeshMotionSolver::solve(DistSVec<double,3> &dX, DistSVec<double,3> &X)  {

  // HB: dX <- P.dX where P is the projector onto  the sliding type of constraints
  // WARNING: assume the homogeneous Dirichlet BCs have already been applied to dX
  //com->fprintf(stdout, "Received 'unprojected' incr disp = %e\n",dX.norm());

  applyProjector(dX); 

  //com->fprintf(stdout, "Received 'projected' incr disp = %e\n",dX.norm());

  dX0 = &dX;

  if (cs) cs->solve(dX, X, meshMotionBCs); //HB: we may also want to apply the rotation to 
  //the normal directions used in the projections in BCApplier ...Not currently done, because 
  //it is assumed that in most practical cases (to date) there would be only one sliding plane,
  //and this sliding plane would also be the symmetry plane so that the plane normal wouldn't be
  // affected by the rotation of the fluid mesh around the symmetry plane normal ...

  ns->solve(X);

  return 0;
}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::setup(DistSVec<double,3> &X)
{
  if(cs) cs->setup(X);

}
//------------------------------------------------------------------------------

void TetMeshMotionSolver::printf(int verbose, const char *format, ...)
{
  if (com->cpuNum() == 0 && verbose <= com->getMaxVerbose()) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }
}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::computeFunction(int it, DistSVec<double,3> &X, 
					  DistSVec<double,3> &F) 
{

  DistMat<PrecScalar,3> *_pc = dynamic_cast<DistMat<PrecScalar,3> *>(pc);

  // PJSA FIX
  if(it == 0 && (typeElement == DefoMeshMotionData::NON_LINEAR_FE 
     || typeElement == DefoMeshMotionData::NL_BALL_VERTEX)) {
    X += *dX0; 
  }

  domain->computeStiffAndForce(typeElement, X, F, *mvp, _pc, volStiff);

  // PJSA FIX 
  if (it == 0 && !(typeElement == DefoMeshMotionData::NON_LINEAR_FE 
      || typeElement == DefoMeshMotionData::NL_BALL_VERTEX)) { // compute F0 <- F0 + [Kib*dXb,0] & X <- X + [0,dXb]
      mvp->BCs = 0;
      mvp->apply(*dX0, *F0);
      mvp->BCs = meshMotionBCs;
      F += *F0;
      X += *dX0;
    }

  // PJSA FIX
  if(meshMotionBCs) meshMotionBCs->applyPD(F);

}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::computeJacobian(int it, DistSVec<double,3> &X, 
					  DistSVec<double,3> &F) 
{

}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::setOperators(DistSVec<double,3> &X)
{

  double t0 = timer->getTime();

  pc->setup();
  
  double t = timer->addMeshPrecSetupTime(t0);

  com->printf(6, "Mesh preconditioner computation: %f s\n", t);

}

//------------------------------------------------------------------------------

int TetMeshMotionSolver::solveLinearSystem(int it, DistSVec<double,3> &rhs, 
					   DistSVec<double,3> &dX) 
{

  double t0 = timer->getTime();

  dX = 0.0;

  ksp->setup(it, maxItsNewton, rhs);

  int lits = ksp->solve(rhs, dX);

  // PJSA FIX (note rhs has already been projected in computeFunction)
  if(meshMotionBCs) meshMotionBCs->applyPD(dX);

  timer->addMeshKspTime(t0);
  
  return lits;

}

//--------------------------------------------------------------------------------------------------

