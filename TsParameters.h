#ifndef _TS_PARAMETERS_H_
#define _TS_PARAMETERS_H_

class IoData;

//------------------------------------------------------------------------------

class TsParameters {

  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflMax;
  double cflMin;
  double ser;

public:

  int maxIts;
  int resType;
  double eps;
  double maxTime;
  double cfl;
  double residual;

  char *output;

public:

  TsParameters(IoData &);
  ~TsParameters();

  void computeCflNumber(int, double);
  double getCflMinOverCfl0(){return (cflMin/cfl0);}


// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
