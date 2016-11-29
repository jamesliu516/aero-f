//
// Created by lei on 5/16/16.
//

#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedRomTsDesc.h>
#endif

#include <VecSetOp.h>
#include <RefVector.h>

#ifndef DEBUG
#define DEBUG 10
#endif

template<int dim>
ImplicitEmbeddedRomTsDesc<dim>::ImplicitEmbeddedRomTsDesc(IoData &_ioData,
                                                          GeoSource &geoSource,
                                                          Domain *dom) :
        ImplicitEmbeddedCoupledTsDesc<dim>(_ioData, geoSource, dom),
        reducedJacobian(0, dom->getNodeDistInfo()),
        reducedBasis(0, dom->getNodeDistInfo()),
        reducedNewtonDirection(0),
        embeddedALS(dom->getCommunicator(), _ioData, *dom),
        referenceState(dom->getNodeDistInfo())/*,
        residualRef(this->F) */{
    // load reducedBasis from file; this needs to be done before MatrixVectorProduct class is set
    embeddedALS.readBasisFiles(this->reducedBasis);
    int n = this->reducedBasis.numVectors();
    this->com->barrier();
    this->printf(DEBUG, " ... basis read into memory\n");
    // load reference state from file;
    embeddedALS.readReferenceStateFiles(this->referenceState);
    this->printf(DEBUG, " ... reference state read into memory\n");
    // initialize reduced Dimension
    this->reducedDimension = n;


    // initialize MatVecProd
    this->printf(DEBUG, " ... initialize ImplicitEmbeddedRomTsDesc\n");
    this->printf(DEBUG, " ... ioData.forced.type is %d ( 0 == heaving ) \n", _ioData.forced.type);
    ImplicitData &implicitData = _ioData.ts.implicit;
    this->printf(DEBUG, " ... implicitData.mvp is %d\n", implicitData.mvp);
    switch(implicitData.mvp) {
        case ImplicitData::FD: // finite difference
            this->Jacobian = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp,this->domain, _ioData);
            this->printf(DEBUG, " ... FD matrix vector product set\n");
            break;
        case ImplicitData::H1: // approximate
            this->Jacobian = new MatVecProdH1<dim, double, dim>(this->timeState, this->spaceOp, this->domain, _ioData);
           // this->Jacobian = new MatVecProdRomH1<dim,double,dim>(this->timeState, this->spaceOp, this->domain, _ioData, this->reducedBasis);
            this->printf(DEBUG, " ... H1 matrix vector product set, only works for ROM\n");
            break;
        case ImplicitData::H2: // exact
            this->Jacobian = new MatVecProdH2<dim,double,dim>(_ioData, this->varFcn, this->timeState, this->spaceOp, this->domain, this->geoState);
            this->printf(DEBUG, " ... H2 matrix vector product set\n");
            break;
    }
    /*
    if (implicitData.mvp == ImplicitData::FD){

        mvp = new MatVecProdFD<dim,dim>(implicitData,this->timeState, this->geoState,
                                        this->spaceOp,this->domain,ioData);

    } else if (implicitData.mvp == ImplicitData::H1){

        mvp = new MatVecProdH1<dim,double,dim>(this->timeState, this->spaceOp, this->domain, ioData);

    } else if (implicitData.mvp == ImplicitData::H2){

        mvp = new MatVecProdH2<dim,double,dim>(ioData, this->varFcn, this->timeState,
                                               this->spaceOp, this->domain, this->geoState);

    }
     */

    typename MatVecProd<dim,dim>::_fsi fsi = {
            this->distLSS,
            &this->nodeTag,
            this->riemann,
            this->linRecAtInterface,
            this->viscSecOrder,
            &this->Wtemp,
            this->riemannNormal,
            this->ghostPoints,
    };

    Jacobian->AttachStructure(fsi);
    //test_Jacobian->AttachStructure(fsi);

    if (this->modifiedGhidaglia){
        Jacobian->attachHH(this->embeddedU);
        //test_Jacobian->attachHH(this->embeddedU);
    }

    // create two krylov subspace solver
    /* from parent class
    pc = ImplicitEmbeddedTsDesc<dim>::template
    createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);

    ksp = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);
     */
    rom_pc = ImplicitEmbeddedTsDesc<dim>::template createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);

    rom_ksp = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, Jacobian, rom_pc, this->com);

    // initialize reducedJacobian, = copy assignment not defined
    //reducedJacobian = VecSet<DistSVec<double, dim> >(0, dom->getNodeDistInfo());
    reducedJacobian.resize(reducedDimension);

    // initialize reducedNewtonDirection
    reducedNewtonDirection.resize(this->reducedDimension);

    // initialize parallelRom and scratch pad
    this->LeastSquareSolver = new ParallelRom<dim>(*dom, this->com, dom->getNodeDistInfo());
    RefVec<DistSVec<double, dim> > F_temp(*(this->F));
    LeastSquareSolver->parallelLSMultiRHSInit(this->reducedJacobian, F_temp);
    result = new double*[1];
    result[0] = new double[this->reducedDimension];

    this->printf(DEBUG, " ... initialization completed\n");
}

template<int dim>
ImplicitEmbeddedRomTsDesc<dim>::~ImplicitEmbeddedRomTsDesc() {
    if(Jacobian)            delete Jacobian;
    if(LeastSquareSolver)   delete LeastSquareSolver;
    if(result)              delete [] result;
    // TODO: clear up points to various data structure
}

/* See ImplicitEmbeddedCoupledTsDesc and ImplicitRomTsDesc
 * super::computeJacobian(it, Q, F), result is this->mvp;
 * this->computeJacobian(it, Q, F), retsult in this->Jacobian;
 * */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::computeJacobian(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F) {
    this->printf(DEBUG, " ... entering parent class probDesc->computeJacobian()\n");
    super::computeJacobian(it, Q, F);
    this->printf(DEBUG, " ... leaving parent class  probDesc->computeJacobian()\n");
    this->printf(DEBUG, " ... entering function probDesc->computeJacobian()\n");
    MatVecProdH1<dim, double, dim> *approximateJacobian = dynamic_cast<MatVecProdH1<dim, double, dim> *>(this->Jacobian);
    if (approximateJacobian)
        approximateJacobian->clearGhost();

    if (this->modifiedGhidaglia)
        Jacobian->evaluateHH(*this->hhResidual, *this->bcData->getBoundaryStateHH());

    Jacobian->evaluate(it, *(this->X), *(this->A), Q, F);

    if (approximateJacobian)
        this->domain->setExactBoundaryJacobian(Q, *this->X, this->ioData,
                                               this->currentTime + this->currentTimeStep,
                                               this->spaceOp->getVarFcn(), *approximateJacobian);
    // todo: modify this for embedded case as well? use EmbeddedDistSVec
    for(int i = 0; i < this->reducedDimension; i++){
        Jacobian->apply(reducedBasis[i], reducedJacobian[i]);
    }

    this->printf(DEBUG, " ... leaving function probDesc->computejacobian()\n");
}

/*
template<int dim>
int ImplicitEmbeddedRomTsDesc<dim>::solveLinearSystem(int it, DistSVec<double, dim> &b,
                                                       DistSVec<double, dim> &dQ) {
    this->printf(DEBUG, " ... entering parent class probDesc->solveLinearSystem()\n");
    super::computeJacobian(it, b, dQ);
    this->printf(DEBUG, " ... leaving parent class  probDesc->solveLinearSystem()\n");
    this->printf(DEBUG, " ... entering child class probDesc->solveLinearSystem()\n");

    double t0 = this->timer->getTime();

    this->embeddeddQ = 0.0;
    //this->embeddedB.ghost() = 0.0;
    this->embeddedB = 0.0;
    this->embeddedB.real() = b;

    if (this->modifiedGhidaglia)
        this->embeddedB.hh() = -1.0*(*this->hhResidual);

    Vec<double> rom_b(this->reducedDimension, NULL);
    Vec<double> rom_dQ(this->reducedDimension, NULL);
    rom_dQ = 0.0;
    projectVector(this->reducedBasis, b, rom_b);
    rom_ksp->setup(it, this->maxItsNewton, rom_b);
    int lits = rom_ksp->solve(rom_b, rom_dQ);

    if(lits == rom_ksp->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;

    DistSVec<double, dim> result(this->domain->getNodeDistInfo());
    expandVector(rom_dQ, result);
    DistSVec<double, dim> difference(this->domain->getNodeDistInfo());
    difference = result - dQ;
    this->printf(DEBUG, " ... fom result is %e, rom result is %e and difference is %e\n", dQ.norm(), result.norm(), difference.norm());
    this->embeddedU.ghost() += this->embeddeddQ.ghost();
    if (this->modifiedGhidaglia) {
        this->embeddedU.hh() += this->embeddeddQ.hh();

        *this->bcData->getBoundaryStateHH() = this->embeddedU.hh();
    }

    this->timer->addKspTime(t0);
    this->printf(DEBUG, " ... leaving child class probDesc->solveLinearSystem()\n");

    return lits;

}
*/

/**
 * solve Jx = b for x, where J is this->mvp or this->jacobian
 * parent class result in DistSVec dQ_prime;
 * child class result in dQ;
 */

template<int dim>
int ImplicitEmbeddedRomTsDesc<dim>::solveLinearSystem(int it, DistSVec<double, dim> &rhs,
                                                       DistSVec<double, dim> &dQ) {
    // step 1: compute using parent class
    this->printf(DEBUG, " debugging ImplicitEmbeddedRomTsDesc::solveLinearSystem: calling parent method\n");
    DistSVec<double, dim> dQ_prime(this->domain->getNodeDistInfo());
    super::solveLinearSystem(it, rhs, dQ_prime);
    this->printf(DEBUG, " debugging ImplicitEmbeddedRomTsDesc::solveLinearSystem: done calling parent method\n");

    // step 2: compute using scalapack
    int k = this->reducedDimension;
    this->printf(DEBUG, " ... calling parallelLSMultiRHS\n");
    RefVec<DistSVec<double, dim> > rhs_temp(rhs);
    LeastSquareSolver->parallelLSMultiRHS(this->reducedJacobian, rhs_temp, this->reducedDimension, 1, this->result);
    //this->printf(DEBUG, " ... least square solver done, reduced newton direction is %p\n", (void *) &this->reducedNewtonDirection);
    for(int i = 0; i < k; i++)
        this->reducedNewtonDirection[i] = -(this->result)[0][i];

    // step 3: compute using normal equation
    this->printf(DEBUG, " ... calling GenFullM normal equation solver\n");
    Vec<double> newdirection(k);
    this->_solveLinearSystem(it, rhs, newdirection);
    Vec<double> difference = this->reducedNewtonDirection - newdirection;
    this->printf(DEBUG, " ... comparing normal equation and parallelLS difference: answer is %f\n", difference.norm());
    this->printf(DEBUG, " ... reduced newton search direction expanded into dQ\n");
    expandVector(newdirection, dQ);

    // step 4: compare their values:
    DistSVec<double, dim> error1(this->domain->getNodeDistInfo());
    DistSVec<double, dim> error2(this->domain->getNodeDistInfo());
    error1 = dQ - this
    error2 = dQ - dQ_prime;
    this->printf(DEBUG, " ... comparing HFM and ROM error: answer is %f\n", error.norm());
    this->printf(DEBUG, " leaving probDesc->solveLinearSystem(), answer set to HFM\n");
    dQ = dQ_prime;
    return it;
}


/*
 * return dQ = U * p, U is reduced order bases
 */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::expandVector(Vec<double>& p, DistSVec<double, dim>& dQ){
    dQ = 0.0;
    for (int i = 0; i < this->reducedDimension; i++)
        dQ += this->reducedBasis[i] * p[i];
}

/*
 * return buffer = transpose(mat) * vec
 */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::projectVector(VecSet<DistSVec<double, dim> > &mat, DistSVec<double,dim> &vec, Vec<double> &buffer) {
    Vec<double> temp(this->reducedDimension, NULL);
    transMatVecProd(mat, vec, temp); // temp = transpose(mat) * vec;
    for (int i = 0; i < this->reducedDimension; i++)
        buffer[i] = temp[i];
}

/*
 * see NonlinearRomOnlineII::projectSwitchStateOntoAffineSubspace()
 */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::projectStateOntoROB(DistSVec<double, dim> &U) {
    DistSVec<double, dim> V(this->domain->getNodeDistInfo());

    V = U - this->referenceState;
    int n = this->reducedDimension;
    Vec<double> q(n);
    for(int i = 0; i < n; i++)
        q[i] = (this->reducedBasis)[i] * V;

    DistSVec<double, dim> Vdiff(this->domain->getNodeDistInfo());
    Vdiff = 0;

    for(int i = 0; i < n; i++)
        Vdiff += q[i] * (this->reducedBasis)[i];

    V = V - Vdiff;
    this->com->printf(DEBUG, " ... || original - projected ||_2 / || original ||_2 = %e\n", V.norm() / U.norm());
    U = this->referenceState + Vdiff;

    /**
     * additional diagonastics: checking ghost node values
     */
    this->com->printf(DEBUG, "testing print ghost nodes\n");
    DistVec<bool> is_swept = this->distLSS->getIsActive();
    this->domain->printDistVecBool(is_swept, false);
}

/*
 * print ghost point from DistEmbeddedVec or distLSS ?
 */
template <int dim>
void ImplicitEmbeddedRomTsDesc<dim>::printGhostPoint(){

}


/*
 * various diagnostic tests:
 * 1. testing difference between mvpop->apply()
 * 2.
 */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::test(){
    // quick hacks to fix dt == 0 bug
    this->timeState->setDt(0.05);
    // load initial U
    DistSVec<double, dim> U(this->domain->getNodeDistInfo());
    // initialize solutions and geometry
    this->setupTimeStepping(&U, this->ioData);
    DistSVec<double, dim> F(this->domain->getNodeDistInfo());
    this->computeFunction(0, U, F);
    this->computeJacobian(0, U, F);
    this->com->printf(DEBUG, " ... setting up U (%e),  F (%e) and jacobian\n", U.norm(), F.norm());

    for(int i = 0; i < this->reducedDimension; i++) {
        // start comparing the results
        DistEmbeddedVec<double, dim> embedded_x(this->domain->getNodeDistInfo());
        DistEmbeddedVec<double, dim> embedded_b(this->domain->getNodeDistInfo());
        embedded_b = 0.0;
        embedded_x = 0.0;
        embedded_x.real() = reducedBasis[i];

        // result from apply on distEmbeddedSVec type
        Jacobian->apply(embedded_x, embedded_b);
        DistSVec<double, dim> result_1(this->domain->getNodeDistInfo());
        DistSVec<double, dim> result_2(this->domain->getNodeDistInfo());
        Jacobian->apply(embedded_x, embedded_b);
        result_1 = embedded_b.real();

        // result from apply on distSVec type: same. maybe setting ghost to be different values would be different?
        Jacobian->apply(reducedBasis[i], result_2);
        DistSVec<double, dim> difference(this->domain->getNodeDistInfo());
        difference = result_1 - result_2;
        this->com->printf(DEBUG, " ... compare result of apply on different types of DistSVec, %e, %e, %e\n",
                          result_1.norm(), result_1.norm(), difference.norm());
    }
}


/**
 * internal method of solving linear system
 */
template<int dim>
int ImplicitEmbeddedRomTsDesc<dim>::_solveLinearSystem(int it , DistSVec<double, dim> &rhs, Vec<double> &sol)
{
    int k = this->reducedDimension;
    // initialize x and jacobian
    Vec<double> x(k);
    GenFullM<double> jacobian(k);
    // use transMatMatProd and transMatVecProd function
    double buffer[k * k];
    transMatMatProd(this->reducedJacobian, this->reducedJacobian, buffer);
    for(int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            jacobian[i][j] = buffer[i + j * k];
    double vector[k];
    transMatVecProd(this->reducedJacobian, rhs, vector);
    for(int i = 0; i < k; i++)
        x[i] = - 1.0 * vector[i]; // not sure if it is correct
    /*
    for(int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            jacobian[i][j] = this->reducedJacobian[i] * this->reducedJacobian[j];
    for(int i = 0; i < k; i++)
        x[i] = this->reducedJacobian[i] * rhs;
    */
    // solve the normal equation
    jacobian.factor();
    jacobian.reSolve(x.data()); // solution is stored in x
    sol = x;

    return 0;

}