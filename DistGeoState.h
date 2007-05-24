#ifndef _DIST_GEO_STATE_H_
#define _DIST_GEO_STATE_H_

#include <GeoData.h>

class TimeData;
class Domain;
class GeoState;
class Communicator;

struct Vec3D;

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

class DistGeoState {

  int numLocSub;
  double lscale;
  double oolscale;

  GeoData data;

  DistSVec<double,3> *Xn;
  DistSVec<double,3> *Xnm1;
  DistSVec<double,3> *Xnm2;
  DistSVec<double,3> *Xdot;

  DistVec<double> *ctrlVol_n;
  DistVec<double> *ctrlVol_nm1;
  DistVec<double> *ctrlVol_nm2;

  DistVec<double> *d2wall;

  DistVec<Vec3D> *edgeNorm;
  DistVec<Vec3D> *faceNorm;
  DistVec<double> *edgeNormVel;
  DistVec<double> *faceNormVel;
  DistVec<Vec3D> *edgeNorm_nm1;
  DistVec<Vec3D> *faceNorm_nm1;
  DistVec<double> *edgeNormVel_nm1;
  DistVec<double> *faceNormVel_nm1;
  DistVec<Vec3D> *edgeNorm_nm2;
  DistVec<Vec3D> *faceNorm_nm2;
  DistVec<double> *edgeNormVel_nm2;
  DistVec<double> *faceNormVel_nm2;
  
  DistVec<Vec3D> *inletNodeNorm;
  DistVec<int> *numFaceNeighb; 		//number of faces connected to an inletnode, 
  									//independantly of subdomains
 
  
  Domain *domain;

  GeoState **subGeoState;

  Communicator *com;

public:

  DistGeoState(IoData &, Domain *);
  ~DistGeoState();

  GeoState &operator() (int i) const { return *subGeoState[i]; }

  void setup(char *, TimeData &, DistSVec<double,3> *, DistVec<double> *);
  void setup1(char *, DistSVec<double,3> *, DistVec<double> *);
  void setup2(TimeData &);
  void compute(TimeData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &);
  void interpolate(double, double, DistSVec<double,3> &, DistSVec<double,3> &);
  void update(DistSVec<double,3> &, DistVec<double> &);
  void writeToDisk(char *);

  int getConfig() const { return data.config; }
  DistSVec<double,3> &getXn() const { return *Xn; }
  DistVec<Vec3D> &getInletNodeNorm() const { return *inletNodeNorm; }

};

//------------------------------------------------------------------------------

#endif
