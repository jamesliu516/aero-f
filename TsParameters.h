#ifndef _TS_PARAMETERS_H_
#define _TS_PARAMETERS_H_

#include <ErrorHandler.h>

class IoData;

//------------------------------------------------------------------------------

class TsParameters {

  int cfllaw;

  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflCoef3;
  double cflMax;
  double cflMin;
  double ser;
 
  double angle_growth;
  double angle_zero;
  int dft_history;
  int dft_freqcutoff;
  double dft_growth;
  int fixedunsteady_counter;
 
  double* reshistory;
  complex<double>* dft;

  ErrorHandler* errorHandler;

public:

  int maxIts;
  int resType;
  double eps;
  double maxTime;
  double cfl;
  double residual;
  double dualtimecfl;

  char *output;

  int checksol;
  int checklinsolve;
  bool unphysical;
  bool badlinsolve;
  bool allowstop;
  int forbidreduce;

public:

  TsParameters(IoData &);
  ~TsParameters();

  void computeCflNumber(int, double, double);
  void resolveErrors();
  double getCflMinOverCfl0(){return (cflMin/cfl0);}
  void assignErrorHandler(ErrorHandler* in){errorHandler = in; }

// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
