#include <WallFcn.h>
#include <VarFcn.h>
#include <Vector3D.h>

#include <math.h>

//------------------------------------------------------------------------------

template<int neq>
void WallFcn::computeJacobianSurfaceTerm(int code, Vec3D &normal, double d2wall[3],
					 double *Vwall, double *V[3], 
					 double (*dRdU)[neq][neq])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double oomuw = 1.0 / muw;

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);
  
  double dplus = reynolds * utau * delta * rhow * oomuw;

  double f = 2.5 * log(1.0 + vkcst*dplus) 
    + 7.8 * (1.0 - exp(-eleventh*dplus) - eleventh*dplus*exp(-0.33*dplus));

  double df = 2.5*vkcst/(1.0 + vkcst*dplus) + 7.8*eleventh* 
    (exp(-eleventh*dplus) - exp(-0.33*dplus) + 0.33*dplus*exp(-0.33*dplus));

  double dfdutau = df * reynolds * delta * rhow * oomuw;

  double dfdrho = df * third * reynolds * utau * delta * oomuw;

  double oodFdutau = 1.0 / (f + utau * dfdutau); 

  for (int k=0; k<3; ++k) {
    double oorho3 = third / V[k][0];
    Vec3D u = this->varFcn->getVelocity(V[k]);

    double dFdU0 = oorho3 * (u*t) + utau * dfdrho;
    double dFdU1 = - oorho3 * t[0];
    double dFdU2 = - oorho3 * t[1];
    double dFdU3 = - oorho3 * t[2];

    double dutaudU0 = - dFdU0 * oodFdutau;
    double dutaudU1 = - dFdU1 * oodFdutau;
    double dutaudU2 = - dFdU2 * oodFdutau;
    double dutaudU3 = - dFdU3 * oodFdutau;
    
    double trut = 2.0 * rhow * utau;

    double drut2du0 = third * utau*utau + trut * dutaudU0;
    double drut2du1 = trut * dutaudU1;
    double drut2du2 = trut * dutaudU2;
    double drut2du3 = trut * dutaudU3;

//    t *= norm;

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = - drut2du0 * t[0] * norm;
    dRdU[k][1][1] = - drut2du1 * t[0] * norm;
    dRdU[k][1][2] = - drut2du2 * t[0] * norm;
    dRdU[k][1][3] = - drut2du3 * t[0] * norm;
    dRdU[k][1][4] = 0.0;

    dRdU[k][2][0] = - drut2du0 * t[1] * norm;
    dRdU[k][2][1] = - drut2du1 * t[1] * norm;
    dRdU[k][2][2] = - drut2du2 * t[1] * norm;
    dRdU[k][2][3] = - drut2du3 * t[1] * norm;
    dRdU[k][2][4] = 0.0;

    dRdU[k][3][0] = - drut2du0 * t[2] * norm;
    dRdU[k][3][1] = - drut2du1 * t[2] * norm;
    dRdU[k][3][2] = - drut2du2 * t[2] * norm;
    dRdU[k][3][3] = - drut2du3 * t[2] * norm;
    dRdU[k][3][4] = 0.0;

    double utw = uw * t * norm;

    dRdU[k][4][0] = - drut2du0 * utw;
    dRdU[k][4][1] = - drut2du1 * utw;
    dRdU[k][4][2] = - drut2du2 * utw;
    dRdU[k][4][3] = - drut2du3 * utw;
    dRdU[k][4][4] = 0.0;
  }

}

//------------------------------------------------------------------------------
