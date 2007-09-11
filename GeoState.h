#ifndef _GEO_STATE_H_
#define _GEO_STATE_H_

class GeoData;

struct Vec3D;

template<class Scalar> class Vec;
template<class Scalar, int dim> class SVec;

//------------------------------------------------------------------------------

class GeoState {

  GeoData &data;

  Vec<double> &ctrlVol_n;
  Vec<double> &ctrlVol_nm1;
  Vec<double> &ctrlVol_nm2;

  Vec<double> &d2wall;

  Vec<Vec3D> &edgeNorm;
  Vec<Vec3D> &faceNorm;
  Vec<double> &edgeNormVel;
  Vec<double> &faceNormVel;
  Vec<Vec3D> &inletNodeNorm;
  Vec<int> &numFaceNeighb;


public:

  GeoState(GeoData &, Vec<double> &, Vec<double> &, Vec<double> &, Vec<double> &, 
	   Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &, Vec<Vec3D> &, Vec<int>&);
  ~GeoState() {}

  Vec<double> &getCtrlVol_n() const { return ctrlVol_n; }
  Vec<double> &getCtrlVol_nm1() const { return ctrlVol_nm1; }
  Vec<double> &getCtrlVol_nm2() const { return ctrlVol_nm2; }
  Vec<double> &getDistanceToWall() const { return d2wall; }
  Vec<Vec3D> &getEdgeNormal() const { return edgeNorm; }
  Vec<Vec3D> &getFaceNormal() const { return faceNorm; }
  Vec<double> &getEdgeNormalVel() const { return edgeNormVel; }
  Vec<double> &getFaceNormalVel() const { return faceNormVel; }
  Vec<Vec3D> &getInletNodeNormal() const { return inletNodeNorm; }
  Vec<int> &getNumFaceNeighb() const {return numFaceNeighb; }

};

//------------------------------------------------------------------------------

#endif
