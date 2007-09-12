#ifndef _GEO_DATA_H_
#define _GEO_DATA_H_

#include <IoData.h>

//------------------------------------------------------------------------------

class GeoData {

public:

  ImplicitData::Normals typeNormals;
  ImplicitData::Velocities typeVelocities;

  int config;

  bool use_n;
  bool use_nm1;
  bool use_nm2;

public:

  GeoData(IoData &);
  ~GeoData() {}

};

//------------------------------------------------------------------------------

#endif
