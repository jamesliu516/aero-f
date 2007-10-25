#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <stdio.h>
#include <math.h>

//------------------------------------------------------------------------------

struct Vec3D {

  double v[3];

  Vec3D() {v[0] = v[1] = v[2] = 0.0; }
  Vec3D(double x[3]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; }
  Vec3D(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
  Vec3D(const Vec3D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; }
  ~Vec3D() {}

  Vec3D &operator=(const double);
  Vec3D &operator=(const Vec3D &);
  Vec3D &operator+=(const Vec3D &);
  Vec3D &operator+=(const double &);
  Vec3D &operator-=(const Vec3D &);
  Vec3D &operator*=(double);
  Vec3D &operator/=(double);
  Vec3D operator/(double) const;

  Vec3D operator+(const Vec3D &);
  Vec3D operator-(const Vec3D &);
  Vec3D operator^(const Vec3D &);

  double operator*(const Vec3D &);

  operator double*() { return v; }

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stderr, "%s(%e %e %e)\n", msg, v[0], v[1], v[2]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])); }
};

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;
  v[2] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const Vec3D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];
  v[2] = v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator+=(const Vec3D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];
  v[2] += v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;
  v[2] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator-=(const Vec3D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];
  v[2] -= v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec3D Vec3D::operator/(double cst)  const
{
  Vec3D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;
  vnew[2] = v[2] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator+(const Vec3D &v2)
{

  Vec3D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];
  res.v[2] = v[2]+v2.v[2];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator-(const Vec3D &v2)
{

  Vec3D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];
  res.v[2] = v[2]-v2.v[2];

  return res;

}

//------------------------------------------------------------------------------
// define vector cross product

inline 
Vec3D Vec3D::operator^(const Vec3D &v2)
{

  Vec3D res;

  res.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  res.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  res.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];

  return res;

}

//------------------------------------------------------------------------------

inline 
double Vec3D::operator*(const Vec3D &v2)
{

  return v[0]*v2.v[0] + v[1]*v2.v[1] + v[2]*v2.v[2];

}

//------------------------------------------------------------------------------

inline 
Vec3D operator*(double c, const Vec3D &v)
{

  Vec3D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];

  return res;

}

//------------------------------------------------------------------------------

#endif
