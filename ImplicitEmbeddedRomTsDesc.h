//
// Created by lei on 5/16/16.
// for running with _EMBEDDED_ALS_ROM_ONLINE_

#ifndef PROJECT_IMPLICITEMBEDDEDROMTSDESC_H
#define PROJECT_IMPLICITEMBEDDEDROMTSDESC_H

#include <MatVecProd.h>                             // for Jacobian
#include <ImplicitEmbeddedTsDesc.h>                 // base calss, time integrator
#include <VectorSet.h>                              // vector and matrix class
#include <ParallelRom.h>                            // for scalapack least square
#include <EmbeddedAlternatingLeastSquare.h>         // for reading basis

// TODO: inherit from ImplicitEmbeddedTsDesc instead
// TODO: inherit from ImplictRomTsDesc instead
// see ImplicitEmbeddedCoupledTsDesc.C
template<int dim>
class ImplicitEmbeddedRomTsDesc : public ImplicitEmbeddedTsDesc<dim> {
private:
    typedef ImplicitEmbeddedTsDesc<dim> super;
    ParallelRom<dim> *LeastSquareSolver;
    MatVecProd<dim, dim> *Jacobian;
    VecSet<DistSVec<double, dim> > reducedJacobian;
    Vec<double> reducedNewtonDirection;
    //TODO: two methods: project U into Qy and lift y to U

    int reducedDimension;
    VecSet<DistSVec<double, dim> > reducedBasis;
    double **result; //<! scrachpad for scalapack to store result, fixed size from initialization

    // internal methods
    void expandVector(Vec<double>& p, DistSVec<double, dim>& dQ);
    void projectVector(VecSet<DistSVec<double, dim> > &mat, DistSVec<double,dim> &vec, Vec<double> &buffer);

    // misc variables to pass c++ compiler hoops.
    EmbeddedAlternatingLeastSquare<dim> embeddedALS;

public:
    /** @name interface to NewtonSolver
     * REQUIRED functions to call NewtonSolver (backtracking linear search)
     * listed here explicitly for clarity
     */
    ///@{
    int getMaxItsNewton() { return super::getMaxItsNewtion(); };
    double getEpsNewton() { return super::getMaxItsNewtion(); };
    double getEpsAbsResNewton() { return super::getEpsAbsResNewton(); };
    double getEpsAbsIncNewton() {return supper::getAbsIncNewton(); };
    FILE* getOutputNewton() { return super::getOutputNewton(); };
    bool getLineSearch() { return super::getLineSearch(); };
    double getContractionLinSearch() {return super::getContractionLinSearch(); };
    double getSufficientDecreaseLinSearch() {return super::getSufficientDecreaseLinSearch();};
    double getMaxItsLineSearch() {return super::getMaxItsLineSearch(); };
    double getNewtonIt(){return super::getNewtonIt(); };
    double getNumResidualOutputCurrentNewtonIt() {return super::getNumResidualOutputCurrentNewtonIt(); };
    // equivalent to computeFullResidual(), result in F ?
    void computeFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F) { super::computeFunction(it, Q, F); };
    void recomputeResidual(DistSVec<double, dim> &F, DistSVec<double, dim> &Finlet) { super::recomputeResidual(F, Finlet); };
    // add inlet contribution, result in rhs ?
    void recomputeFunction(DistSVec<double, dim> &Q, DistSVec<double, dim> &rhs){ super::recomputeFunction(Q, rhs); };
    void setOperators(DistSVec<double, dim> &) {};
    // TODO: emulate coupledTsDesc (compute Hessian of function)
    void computeJacobian(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F);
    // TODO: emulate solveNewtonSystem() (get search direction in dQ)
    int solveLinearSystem(int it, DistSVec<double, dim> &rhs, DistSVec<double, dim> &dQ);
    ///@}

    /** @name interface to NewtonSolver
     * OPTIONAL functions to call NewtonSolver (backtracking linear search)
     */
    ///@{
    //void fprintf(FILE *, const char *, ...); // already in TsDesc.h
    void writeBinaryVectorsToDiskRom(bool, int, DistSVec<double, dim> &, DistSVec<double, dim> &) {};
    void setCurrentStateForKspBinaryOutput(DistSVec<double, dim> &Q) { super::setCurrentStateForKspBinaryOutput(Q); };
    int checkSolution(DistSVec<double, dim> &Q) {return super::checkSolution(Q); };
    void fixSolution(DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ) { super::fixSolution(Q, dQ); };
    bool checkFailSafe(DistSVec<double, dim> &Q) {return super::checkFailSafe(Q); };
    void resetFixesTag() { super::resetFixesTag(); };
    ///@}

    /** @name interface to TsSolver
     * REQUIRED functions to run within TsSolver
     * override base methods
     */
    ///@{
    //TODO: overwrite this to use reduced coordinate ?
    int solveNonlinearSystem(DistSVec<double, dim> &Q, const int totalTimeSteps) { return super::solveNonLinearSystem(Q, totalTimeSteps); };
    ///@}
    /*
    /* must implement, or compilation error
    void solveNewtonSystem(const int &it, double &res, bool &breakloop,
                           DistSVec<double, dim> &U,
                           const int &totalTimeSteps = 0);
    */
    ImplicitEmbeddedRomTsDesc(IoData &, GeoSource &, Domain *);
    ~ImplicitEmbeddedRomTsDesc();
};

#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedRomTsDesc.cpp>
#endif

#endif //PROJECT_IMPLICITEMBEDDEDROMTSDESC_H
