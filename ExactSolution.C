// Exact Solution.cpp

#include <IoData.h>

#include "ExactSolution.h"

#include <math.h>

void ExactSolution::
AcousticBeam(IoData& iod,double x, double y, double z,
	     double t, double* V) {

 /* 
  double alpha = 3.810627283;//6.28267340353564;
  double omega0 = 2697.56348148;//26.975634814780758;
  double omegatilde = 0.60774519311;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;
  */

  double alpha = 3.8116859813491843;//6.28267340353564;
  double omega0 = 2696.889343398955;//26.975634814780758;
  double omegatilde = 0.6077988242743707;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;
 
  for (int i = 0; i < 5; ++i)
    V[i] = 0.0;
  
  double u = omega*k*what* 
    (cosh(alpha*y) /(alpha*sinh(alpha*H)))*(cos(k*x-omega*t*iod.ref.rv.time)+
					       cos(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
  double v = omega*what*(sinh(alpha*y)/sinh(alpha*H))*
    (sin(k*x-omega*t*iod.ref.rv.time)-
     sin(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
  double p = omega*omega*rhof*cosh(alpha*y)*what/(alpha*sinh(alpha*H))*
    (cos(k*x-omega*t*iod.ref.rv.time)-
     cos(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.pressure + iod.bc.inlet.pressure;

  V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.3 / iod.ref.rv.density);
  V[1] = u;
  V[2] = v;
  V[3] = 0.0;
  V[4] = p;

}

void ExactSolution::
AcousticBeamStructure(IoData& iod,double x, double y, double z,
  	              double t, double& uy, double& vy) {

  
  double alpha = 3.8116859813491843;//6.28267340353564;
  double omega0 = 2696.889343398955;//26.975634814780758;
  double omegatilde = 0.6077988242743707;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;

  uy = what* (cos(k*x-omega*t*iod.ref.rv.time)-cos(-k*x-omega*t*iod.ref.rv.time));
  vy = omega*what* (sin(k*x-omega*t*iod.ref.rv.time)-
	             sin(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
}


