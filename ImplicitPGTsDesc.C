//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  residualRef(this->F),
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {
  parallelRom = new ParallelRom<dim>(*dom,this->com);
  parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);  
  lsCoeff = new double*[1];
  lsCoeff[0] = new double[this->nPod];
}

//------------------------------------------------------------------------------
template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc(){

    delete [] lsCoeff[0];
    delete [] lsCoeff;  

}
//-----------------------------------------------------------------------------
template<int dim>
void ImplicitPGTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

  projectVector(this->AJ, this->F, From);
  Vec<double> rhs(this->nPod);
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

/*
    // Option 1: Solve the normal equations
    // form reduced Jacobian
    this->jac.setNewSize(this->nPod,this->nPod);
    for (int iRow = 0; iRow < this->nPod; ++iRow) {
      for (int iCol = 0; iCol <= iRow; ++iCol) {
	this->jac[iRow][iCol] = this->AJ[iRow]*this->AJ[iCol];
	if (iRow > iCol)
	  this->jac[iCol][iRow] = this->jac[iRow][iCol];
      }
    } 

    solveLinearSystem(it, rhs, this->dUrom);
*/
    //Option 2: solve a Least Squares Problem
    RefVec<DistSVec<double, dim> > residualRef2(this->F);
    parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoeff);
    for (int iPod=0; iPod<this->nPod; ++iPod)
      this->dUrom[iPod] = -lsCoeff[0][iPod];



}
