#ifndef _TIME_DATA_H_
#define _TIME_DATA_H_

#include <IoData.h>

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

class TimeData {

public:

  ImplicitData::Type typeIntegrator;
  ImplicitData::Startup typeStartup;
  TsData::TypeTimeStep typeTimeStep;

  double dt_imposed;
  double dt_n;
  double dt_nm1;
  double dt_nm2;

  double tau_n;
  double tau_nm1;
  double alpha_np1;
  double alpha_n;
  double alpha_nm1;
  double alpha_nm2;

  bool exist_nm1;
  bool exist_nm2;
  bool use_nm1;
  bool use_nm2;

  bool use_freq;
  bool use_modal;

public:

  TimeData(IoData &);
  ~TimeData() {}

  void update();
  void computeCoefficients(DistVec<double> &, double);
  void computeVelocities(DGCLData::Velocities, DistSVec<double,3> &,
			 DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &);

// Included
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
