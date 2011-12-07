#include <Communicator.h>

#include <cmath>
#include <iostream>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomPostproTsDesc<dim>::ImplicitRomPostproTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), Uinitial(dom->getNodeDistInfo())
{
	this->maxItsNewton = 0;	// never do iterations

 char *redCoordsFile = this->input->staterom;
 readRedCoords = fopen(redCoordsFile, "r");
 if (!readRedCoords)  {
   this->com->fprintf(stderr, "*** Warning: cannot open %s\n", redCoordsFile);
   exit (-1);
 }

 while (fgetc(readRedCoords) != '\n') ;	// ignore first line

 int tmp, _n;
 double tmp2;
 _n = fscanf(readRedCoords, "%d %lf", &tmp, &tmp2);	// ignore first time step (initial condition)
 for (int iPod = 0; iPod < this->nPod; ++iPod) {
	 _n = fscanf(readRedCoords, "%lf", &tmp2);
 }
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeFullResidual(int , DistSVec<double, dim> &) {


	// do nothing
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeAJ(int, DistSVec<double, dim> &)  {

	// do nothing
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop)  {

  breakloop = true;	// after loop, exit will occur because maxItsNewton = 1
}

template<int dim>
void ImplicitRomPostproTsDesc<dim>::postProStep(DistSVec<double, dim> &U, int totalTimeSteps)  {
	
	if (totalTimeSteps == 1) {	// first time step
		Uinitial = U;
	}

	int tmp, _n;
	double tmp2;
	_n = fscanf(readRedCoords, "%d %lf", &tmp, &tmp2);
	for (int iPod = 0; iPod < this->nPod; ++iPod) {
		_n = fscanf(readRedCoords, "%lf", &tmp2);
		this->UromTotal[iPod] = tmp2;	// set dUrom = UromTotal (subtract next)
	}

	expandVector(this->UromTotal, U); // solution increment in full coordinates
	U += Uinitial;
}
