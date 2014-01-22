// Exact Solution.cpp

#include <IoData.h>

#include "ExactSolution.h"

#include <math.h>

#ifdef AEROACOUSTIC

#include "gsl/gsl_sf.h"

#endif // AEROACOUSTIC


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

#ifdef AEROACOUSTIC
double j0prime(double r) {

  return -gsl_sf_bessel_J1(r);
}

double y0prime(double r) {

  return -gsl_sf_bessel_Y1(r);
}

#endif

void ExactSolution::
CylindricalBubble(IoData& iod,double x, double y, double z,
		  double t, double* V, double* phi, int& fid) {

#ifdef AEROACOUSTIC
  double omega = 733.2541686104948;
  double Bn = 10.0;//1.0e-1;
  double r = sqrt(x*x+y*y);
  double rhoo = 10.0;
  double co = sqrt(4.4*1e5/rhoo);
  double rhoi = 1.0;
  double ci = sqrt(1.4*1e5/rhoi);
  double a = 0.5;
  double R = 1.0;
  double pinf = 1e5;
  
  t *= iod.ref.rv.time;

  double theta = atan2(y,x);

  double Cn = -Bn*j0prime(omega*R / co) / y0prime(omega*R / co);

  double An = (Bn*gsl_sf_bessel_J0(omega*a / co) + Cn*gsl_sf_bessel_Y0(omega*a / co))/gsl_sf_bessel_J0(omega*a / ci);

  double ur, utheta = 0;
  double p;


  if (r < a) {
    
    ur = -An/(rhoi*ci)*j0prime(omega*r/ci)*cos(omega*t);
    p = -An*gsl_sf_bessel_J0(omega*r / ci) *sin(omega*t);
  } else {
    ur = -1.0/(rhoo*co)*(Bn*j0prime(omega*r/co)+Cn*y0prime(omega*r/co))*cos(omega*t);
    p = -(Bn*gsl_sf_bessel_J0(omega*r / co) + Cn*gsl_sf_bessel_Y0(omega*r / co))*sin(omega*t);
  }

  double ux = ur*cos(theta), uy = ur*sin(theta);

  p += pinf;

  p /= iod.ref.rv.pressure;

  if (r < a)
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.0 / iod.ref.rv.density);
  else
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/4.4))*(10.0 / iod.ref.rv.density);

  V[1] = ux/ iod.ref.rv.velocity;
  V[2] = uy/ iod.ref.rv.velocity;
  
  V[3] = 0.0;
  
  V[4] = p;  

  fid = (r <= a);

  *phi = a-r;
  
#else
  std::cout << "Error: Cylindrical bubble exact solution requires the code to be compiled with the AEROACOUSTIC option" << std::endl;
  exit(-1);
#endif
}
