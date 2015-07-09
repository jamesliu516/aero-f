#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitCollocationTsDesc<dim>::ImplicitCollocationTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitGappyTsDesc<dim>(ioData, geoSource, dom)
	//, leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{

  this->projVectorTmp = NULL;
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitCollocationTsDesc<dim>::~ImplicitCollocationTsDesc() {}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCollocationTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  // Form the normal equations
  double t0 = this->timer->getTime();

  this->projectVector(*this->AJRestrict, *this->ResRestrict, From);     // different from PG
  rhs = -1.0 * From;
  
  res = rhs*rhs;
  
  if (res < 0.0){
    fprintf(stderr, "*** negative residual: %e\n", res);
    exit(1);
  }
  res = sqrt(res);
 
  transMatMatProd(*this->AJRestrict,*this->AJRestrict,this->jactmp);
  for (int iRow = 0; iRow < this->nPod; ++iRow) {
    for (int iCol = 0; iCol < this->nPod; ++iCol) { // different from PG
      this->jac[iRow][iCol] = this->jactmp[iRow + iCol * this->pod.numVectors()];
    }
  }

  this->timer->addLinearSystemFormTime(t0);

  // Solve the normal equations
  t0 = this->timer->getTime();
  this->solveLinearSystem(it, rhs, this->dUromNewtonIt);
  this->timer->addLinearSystemSolveTime(t0);

  // Convergence criterion
  t0 = this->timer->getTime();
  res = this->dUromNewtonIt.norm();
  if (it == 0) {
    this->res0 = res;
    this->target = this->epsNewton * this->res0;
  }
  breakloop = (res == 0.0) || (res <= this->target);
  this->timer->addCheckConvergenceTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCollocationTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {

  this->jac.setNewSize(this->nPod,this->nPod);
  From.resize(this->nPod);
  rhs.resize(this->nPod);

  if (this->projVectorTmp) delete [] (this->projVectorTmp);
  this->projVectorTmp = new double [this->nPod];

  //nSampleNodes = this->rom->sampleNodes.size();
  //leastSquaresSolver.problemSizeIs(nSampleNodes*dim, this->nPod);

  if (this->jactmp) delete [] this->jactmp;
  if (this->column) delete [] this->column;
  this->jactmp = new double [this->nPod * this->nPod];
  this->column = new double [this->nPod];

  this->AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->rom->getRestrictedDistInfo()));
  this->ResRestrict.reset(new DistSVec<double, dim> (this->rom->getRestrictedDistInfo())); 

}

//------------------------------------------------------------------------------

