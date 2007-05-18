#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include <stdlib.h>
#include <math.h>

/*@ARTICLE{brown-saad-90,
  author = "Brown, P. N. and Saad, Y.",
  title = "Hybrid {K}rylov Methods for
           Nonlinear Systems of Equations",
  journal = siamjscistat,
  year = 1990,
  volume = 11,
  number = 3,
  pages = "450--481",
}
@ARTICLE{keyes-venkatakrishnan-96,
  author = "Keyes, D. E. and Venkatakrishnan, V.",
  title = "{N}ewton-{K}rylov-{S}chwarz Methods: Interfacing Sparse Linear
  Solvers with Nonlinear Applications",
  journal = zamm,
  year = 1996,
  volume = 76,
  pages = "147--150",
} 
*/
//------------------------------------------------------------------------------

template<class ProblemDescriptor>
class NewtonSolver {

  ProblemDescriptor *probDesc;

  typedef DistVec<double> VolVecType;
  typename ProblemDescriptor::SolVecType F;  // nonlinear function
  typename ProblemDescriptor::SolVecType Finlet;  // nonlinear function at inlet nodes
  typename ProblemDescriptor::SolVecType dQ; // gradient of F
  typename ProblemDescriptor::SolVecType rhs; // right hand side
  typename ProblemDescriptor::VolVecType dPhi; // 
  typename ProblemDescriptor::VolVecType PhiF; // 
  typename ProblemDescriptor::VolVecType rhsPhi; // 

public:

  NewtonSolver(ProblemDescriptor *);
  ~NewtonSolver() {}
  int solve(typename ProblemDescriptor::SolVecType &);
  int solve(typename ProblemDescriptor::SolVecType &, typename ProblemDescriptor::VolVecType &);
  int solveLS(typename ProblemDescriptor::VolVecType &,typename ProblemDescriptor::VolVecType &,typename ProblemDescriptor::VolVecType &, typename ProblemDescriptor::SolVecType &);

};

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
NewtonSolver<ProblemDescriptor>::NewtonSolver(ProblemDescriptor *prbd) : 
  F(prbd->getVecInfo()), Finlet(prbd->getVecInfo()), dQ(prbd->getVecInfo()), rhs(prbd->getVecInfo()), dPhi(prbd->getVecInfo()), PhiF(prbd->getVecInfo()), rhsPhi(prbd->getVecInfo())
{

  probDesc = prbd;

}

//------------------------------------------------------------------------------------------------------

template<class ProblemDescriptor>
int
NewtonSolver<ProblemDescriptor>::solve(typename ProblemDescriptor::SolVecType &Q )
{
  double res, target;

  int fsIt = 0;
  int maxIts = probDesc->getMaxItsNewton();
  double eps = probDesc->getEpsNewton();

  double res0, res2=0.0;
  int it;
  for (it=0; it<maxIts; ++it) {

    // compute the nonlinear function value
    probDesc->computeFunction(it, Q, F);
    res2 = probDesc->recomputeResidual(F, Finlet);
    res = F*F-res2;
    res = sqrt(F*F-res2);

    if (it == 0) {
      target = eps*res; 
      res0 = res;
    }

    if (res == 0.0 || res <= target) break;

    rhs = -1.0 * F;

    probDesc->recomputeFunction(Q, rhs);
    probDesc->computeJacobian(it, Q, F);

    // apply preconditioner if available
    probDesc->setOperators(Q);

    probDesc->solveLinearSystem(it, rhs, dQ);

    rhs = Q;
    Q += dQ;

    // verify that the solution is physical
    if (probDesc->checkSolution(Q)) {
      if (probDesc->checkFailSafe(Q) && fsIt < 5) {
	probDesc->printf(1, "*** Warning: Newton solver redoing iteration %d\n", it+1);
	Q = rhs;
	--it;
	++fsIt;
      }
      else
	exit(1);
    }

  }
    if (fsIt > 0 && probDesc->checkFailSafe(Q) == 1)
      probDesc->resetFixesTag();

  if (it == maxIts && maxIts != 1) {
    probDesc->printf(1, "*** Warning: Newton solver reached %d its", maxIts);
    probDesc->printf(1, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);    
  }

  return it;

}

//------------------------------------------------------------------------------
template<class ProblemDescriptor>
int
NewtonSolver<ProblemDescriptor>::solveLS(typename ProblemDescriptor::VolVecType &Phi, typename ProblemDescriptor::VolVecType &Phi1, typename ProblemDescriptor::VolVecType &Phi2, typename ProblemDescriptor::SolVecType &U )
{

  double res, target, res1;

  int fsIt = 0;
  //int maxIts = probDesc->getMaxItsNewton();
  //double eps = probDesc->getEpsNewton();
  int maxIts  = 2;
  double eps  = 0.01;
  int it;
  for (it=0; it<maxIts; ++it) {

    // compute the nonlinear function value for the Level Set Equation
    probDesc->computeFunctionLS(it, Phi, Phi1, Phi2, U, PhiF);

    res = sqrt(PhiF*PhiF);

    if (it == 0) target = eps*res;
                         
    if (it ==0) res1  = res;
                                                                                           
    if (res == 0.0 || res <= target) break;
                                                                                                                    
    rhsPhi = -1.0* PhiF;
    
    // compute the Jacobian for the Level Set Equation 
    probDesc->computeJacobianLS(it, Phi, Phi1, Phi2, U, PhiF);
 
    // Solve the linearized system of equations
    probDesc->solveLinearSystemLS(it, rhsPhi, dPhi);

    rhsPhi  = Phi;
    Phi += dPhi;

  }
  if (it == maxIts && maxIts != 1) {
    //probDesc->printf(1, "*** Warning: Newton solver for LS reached %d its", maxIts);
    //probDesc->printf(1, " (initial=%.2e, res=%.2e, target=%.2e)\n", res1, res, target);
  }

  int iti = probDesc->getInitialIteration();
  //if (iti % 100  == 0 )  probDesc->printf(1,"Radius of Bubble %f\n",probDesc->reinitLS(Phi, U, iti));

  return it;
}
//------------------------------------------------------------------------------

#endif
