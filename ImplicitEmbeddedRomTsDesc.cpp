//
// Created by lei on 5/16/16.
//

#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedRomTsDesc.h>
#endif

#include <VecSetOp.h>
#include <RefVector.h>

#ifndef DEBUG
#define DEBUG 9
#endif

template<int dim>
ImplicitEmbeddedRomTsDesc<dim>::ImplicitEmbeddedRomTsDesc(IoData &_ioData,
                                                          GeoSource &geoSource,
                                                          Domain *dom) :
        ImplicitEmbeddedTsDesc<dim>(_ioData, geoSource, dom),
        reducedJacobian(0, dom->getNodeDistInfo()),
        reducedBasis(0, dom->getNodeDistInfo()),
        reducedNewtonDirection(0),
        embeddedALS(dom->getCommunicator(), _ioData, *dom)/*,
        residualRef(this->F) */{
    // initialize MatVecProd
    this->com->fprintf(stdout, "initialize ImplicitEmbeddedRomTsDesc\n");
    ImplicitData &implicitData = _ioData.ts.implicit;
    switch(implicitData.mvp) {
        case ImplicitData::FD: // finite difference
            this->Jacobian = new MatVecProdFD<dim,dim>(implicitData,this->timeState, this->geoState, this->spaceOp,this->domain, _ioData);
            break;
        case ImplicitData::H1: // approximate
            this->Jacobian = new MatVecProdH1<dim,double,dim>(this->timeState, this->spaceOp, this->domain, _ioData);
            break;
        case ImplicitData::H2: // exact
            this->Jacobian = new MatVecProdH2<dim,double,dim>(_ioData, this->varFcn, this->timeState, this->spaceOp, this->domain, this->geoState);
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

    if (this->modifiedGhidaglia)
        Jacobian->attachHH(this->embeddedU);

    // load reducedBasis from file
    embeddedALS.readBasisFiles(this->reducedBasis);
    int n = this->reducedBasis.numVectors();
    this->com->barrier();
    this->com->fprintf(stdout, "basis read into memory\n");
    // initialize reduced Dimension
    this->reducedDimension = n;

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

    this->com->fprintf(stdout, "initialization completed\n");
}

template<int dim>
ImplicitEmbeddedRomTsDesc<dim>::~ImplicitEmbeddedRomTsDesc() {
    if(Jacobian)            delete Jacobian;
    if(LeastSquareSolver)   delete LeastSquareSolver;
    if(result)              delete [] result;
    // TODO: clear up points to various data structure
}

/* See ImplicitEmbeddedCoupledTsDesc and ImplicitRomTsDesc */
template<int dim>
void ImplicitEmbeddedRomTsDesc<dim>::computeJacobian(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F) {
    this->printf(DEBUG, " debugging: entering function probDesc->computeJacobian()\n");
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

    for(int i = 0; i < this->reducedDimension; i++){
        Jacobian->apply(reducedBasis[i], reducedJacobian[i]);
    }

    this->printf(DEBUG, " debugging: leaving function probDesc->computejacobian()\n");
}

template<int dim>
int ImplicitEmbeddedRomTsDesc<dim>::solveLinearSystem(int it, DistSVec<double, dim> &rhs,
                                                       DistSVec<double, dim> &dQ) {
    this->printf(DEBUG, " ... solve linear system\n");
    RefVec<DistSVec<double, dim> > rhs_temp(rhs);
    // need this to circumvent function only taking reference
    this->printf(Debug, "parameters are (%p, %p, %i, 1, %p)\n", (void *) &reducedJacobian, (void *) &rhs_temp, this->reducedDimension, (void*) &this->result);
    LeastSquareSolver->parallelLSMultiRHS(this->reducedJacobian, rhs_temp, this->reducedDimension, 1, this->result);
    // TODO: segfault here!!!
    this->printf(DEBUG, " ... least square solver done, %p\n", (void *) &this->reducedNewtonDirection);
    for(int i = 0; i < this->reducedDimension; i++){
        this->com->fprintf(stdout, " ... %i iteration\n", i);
        this->reducedNewtonDirection[i] = -(this->result)[0][i];
    }
    this->printf(DEBUG, " ... reduced direction computed\n");
    expandVector(this->reducedNewtonDirection, dQ);
    this->printf(DEBUG, " ... linear system solved \n");
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
    Vec<double> temp;
    transMatVecProd(mat, vec, temp); // temp = transpose(mat) * vec;
    for (int i = 0; i < this->reducedDimension; i++)
        buffer[i] = temp;
}