//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom) {

	// read in sample nodes
	
	dom->readPodBasis(ioData.input.AFile, nJ,A);
	dom->readPodBasis(ioData.input.BFile, nJ,B);
	// read in Afull, Bfull (temporary) (binary files)
	// restrict Afull and Bfull to be A, B
		// read in A, B
  
}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitGappyTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

}
