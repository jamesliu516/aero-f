#ifndef _LOCAL_RIEMANN_H
#define _LOCAL_RIEMANN_H

#include <LinkF77.h>

class VarFcn;

extern "C" {
	void F77NAME(eriemanngw) (const double&, const double&, const double&,
	                          const double&, const double&, const double&,
	                          const double&, const double&, const double &,
	                          const double &, const double&, const double&,
	                          const double &, const double&);
  void F77NAME(eriemanngg) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double &,
                            const double &, const double&, const double&,
 	                          const double &, const double&);
	void F77NAME(eriemannww) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&);
};
//------------------------------------------------------------------------------

class LocalRiemann {
//this is a virtual class.
//only subclasses will be created:
//    -LocalRiemannGfmpGasGas
//    -LocalRiemannGfmpTaitTait
//    -LocalRiemannGfmparGasGas
//    -LocalRiemannGfmparGasTait
//    -LocalRiemannGfmparTaitTait
protected:

public:
	LocalRiemann() {}
	~LocalRiemann() {}

	virtual void computeRiemannSolution(double *Vi, double *Vj,
	    double Phii, double Phij, double *nphi, VarFcn *vf,
	    int &epsi, int &epsj, double *Wi, double *Wj,
			double *rupdatei, double *rupdatej,
			double &weighti, double &weightj, int it){}

};

//------------------------------------------------------------------------------

#endif

