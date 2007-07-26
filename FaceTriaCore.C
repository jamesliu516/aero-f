#include <FaceTria.h>

#include <RefVal.h>
#include <BcDef.h>
#include <Edge.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
using std::min;
using std::max;
using std::swap;
#endif

const double FaceTria::third = 1.0/3.0;
const int FaceTria::edgeEndT[Face::MaxNumNd][2] = { {0,1}, {1,2}, {2,0}, {-1,-1} };


//------------------------------------------------------------------------------
// Computation of the OUTWARD subface normals (only one, they are all equal)
void FaceTria::computeNormal(SVec<double,3> &X, Vec<Vec3D> &faceNorm)
{
  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  faceNorm[normNum] = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));
}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormal(SVec<double,3> &X, Vec3D &faceNorm)
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  faceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormalGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				 SVec<double,3> &Xdot, Vec<Vec3D> &faceNorm, 
				 Vec<double> &faceNormVel)
{

  static double twelfth = 1.0/12.0;

  Vec3D x_n[3] = {Xn[ nodeNum(0) ], Xn[ nodeNum(1) ], Xn[ nodeNum(2) ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum(0) ], Xnp1[ nodeNum(1) ], Xnp1[ nodeNum(2) ]};
  Vec3D xdot[3] = {Xdot[ nodeNum(0) ], Xdot[ nodeNum(1) ], Xdot[ nodeNum(2) ]};

  Vec3D x01_n = x_n[1] - x_n[0];
  Vec3D x02_n = x_n[2] - x_n[0];

  Vec3D x01_np1 = x_np1[1] - x_np1[0];
  Vec3D x02_np1 = x_np1[2] - x_np1[0];

  faceNorm[normNum] = twelfth * (((2.0*x02_np1 + x02_n) ^ x01_np1) + 
				 ((2.0*x02_n + x02_np1) ^ x01_n));

  faceNormVel[normNum] = third * (xdot[0] + xdot[1] + xdot[2]) * faceNorm[normNum];

}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormalEZGCL1(double oodt, SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				   Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{

  Vec3D x_n[3] = {Xn[ nodeNum(0) ], Xn[ nodeNum(1) ], Xn[ nodeNum(2) ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum(0) ], Xnp1[ nodeNum(1) ], Xnp1[ nodeNum(2) ]};

  faceNorm[normNum] = 0.5 * ((x_np1[2] - x_np1[0]) ^ (x_np1[1] - x_np1[0]));
  //EZ1 faceNorm = 0.5 * ((x_n[2] - x_n[0]) ^ (x_n[1] - x_n[0]));

  double vol = Face::computeVolume(x_n[0], x_n[2], x_n[1], x_np1[0], x_np1[2], x_np1[1]);

  faceNormVel[normNum] = oodt * vol;

}

//------------------------------------------------------------------------------
// Get OUTWARD face normal
Vec3D FaceTria::getNormal(Vec<Vec3D> &faceNorm) {

  return faceNorm[normNum];

}

double FaceTria::getNormalVel(Vec<double> &faceNormVel) {
  return faceNormVel[normNum];
}

//------------------------------------------------------------------------------
// Get i-th OUTWARD subface normal
Vec3D FaceTria::getNormal(Vec<Vec3D> &faceNorm, int i) {

  return third*faceNorm[normNum];

}

double FaceTria::getNormalVel(Vec<double> &faceNormVel, int i) {
  return third*faceNormVel[normNum];
}
