// Exact Solution.h

#pragma once

#include <DistVector.h>
#include <VarFcn.h>

class ExactSolution {

 public:

  static void AcousticBeam(IoData&, double x, double y, double z,
			   double t, double* V);

  static void 
  AcousticBeamStructure(IoData& iod,double x, double y, double z,
  	              double t, double& uy, double& vy);
  
  template <void (*F)(IoData&, double,double,double,
				     double,double*), int dim >
    static void Fill(DistSVec<double,dim>& U, DistSVec<double,3>& X,
		     IoData& iod, double t, VarFcn* vf) {

#pragma omp parallel for
    for (int iSub = 0; iSub < U.numLocSub(); ++iSub) {

      double v[dim];
      SVec<double,3>& x(X(iSub));
      for (int i = 0; i < U(iSub).size(); ++i) {

	F(iod, x[i][0], x[i][1], x[i][2],
	  t, v);
	vf->primitiveToConservative(v, U(iSub)[i], 0);
      }
    }
    
  }
};
