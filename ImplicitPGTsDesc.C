//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  rhsVS(1, (*this->domain).getNodeDistInfo()),
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), parallelRom(*dom,this->com), From(this->nPod), rhs(this->nPod) {
  parallelRom.parallelLSMultiRHSInit(this->AJ,rhsVS);  
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

//    res=this->F*this->F;
  
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
    //double **lsCoeff;
    //lsCoeff = new double*[1];
    //lsCoeff[0] = new double[this->nPod];
// copy this->F in rhsVS
  //parallelRom.parallelLSMultiRHSInit(this->AJ,rhsVS);


    rhsVS[0] = (this->F);
     parallelRom.parallelLSMultiRHS(this->AJ,rhsVS,this->nPod,1,lsCoeff);
    for (int iPod=0; iPod<this->nPod; ++iPod)
      this->dUrom[iPod] = -lsCoeff[0][iPod];

    //delete [] lsCoeff[0];
    //delete [] lsCoeff;


}
