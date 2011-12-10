#include <VecSetOp.h>
//------------------------------------------------------------------------------

template<int dim>
ImplicitGalerkinTsDesc<dim>::ImplicitGalerkinTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {
  
  this->jac.setNewSize(this->nPod,this->nPod);
  jactmp = new double [this->nPod * this->nPod];
  this->projVectorTmp = new double [this->nPod];
}

template<int dim>
ImplicitGalerkinTsDesc<dim>::~ImplicitGalerkinTsDesc()
{

  if (jactmp) delete [] jactmp;
  if (this->projVectorTmp) delete [] this->projVectorTmp;
}
//------------------------------------------------------------------------------


template<int dim>
void ImplicitGalerkinTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		projectVector(this->pod, this->F, From);	// different from PG
		rhs = -1.0 * From;

		// KTC FIX!
		// saving residual vectors (for GappyPOD)
		//writeBinaryVectorsToDisk1(false, it, 0.0, this->F, Dummy);

    res = rhs*rhs;

    if (res < 0.0){
      fprintf(stderr, "*** negative residual: %e\n", res);
      exit(1);
    }
    res = sqrt(res);

    if (it == 0) {
      this->target = this->epsNewton*res;
      this->res0 = res;
    }

    if (res == 0.0 || res <= this->target) {
			breakloop = true;
			return;	// do not solve the system
		}

		// form reduced Jacobian

		transMatMatProd(this->pod,this->AJ,jactmp);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->pod.numVectors()];
			}
		} 

    solveLinearSystem(it, rhs, this->dUrom);
}
