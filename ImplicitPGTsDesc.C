#include <VecSetOp.h>
//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {

		this->jac.setNewSize(this->nPod,this->nPod);
		jactmp = new double [this->nPod * this->nPod];
  
}

template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc()
{

		if (jactmp) delete [] jactmp;
  
}
//------------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

		projectVector(this->AJ, this->F, From);
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

		// faster way
		transMatMatProd(this->AJ,this->AJ,jactmp);
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol < this->nPod; ++iCol) {	// different from PG
				this->jac[iRow][iCol] = jactmp[iRow + iCol * this->AJ.numVectors()];
			}
		} 

		/*
		// old way
		for (int iRow = 0; iRow < this->nPod; ++iRow) {
			for (int iCol = 0; iCol <= iRow; ++iCol) {
				this->jac[iRow][iCol] = this->AJ[iRow] * this->AJ[iCol];
				if (iCol < iRow)
					this->jac[iCol][iRow] = this->jac[iRow][iCol];
			}
		} 
		*/

    solveLinearSystem(it, rhs, this->dUrom);

}
