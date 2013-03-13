#include <cstring>
#include <TsParameters.h>
#include <IoData.h>

#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif
#include <cstring>

//------------------------------------------------------------------------------

TsParameters::TsParameters(IoData &ioData)
{

  resType = ioData.ts.residual;
  cfl0 = ioData.ts.cfl0;
  cflCoef1 = ioData.ts.cflCoef1;
  cflCoef2 = ioData.ts.cflCoef2;
  cflCoef3 = ioData.ts.cflCoef3;  // undocumented -- coefficient for quadratic CFL law
  cflMax = ioData.ts.cflMax;
  cflMaxSlope = ioData.ts.cflMaxSlope;  //undocumented -- slope of cflMax curve
  cflMin = ioData.ts.cflMin;
  ser = ioData.ts.ser;
  dualtimecfl = ioData.ts.dualtimecfl;

  maxIts = ioData.ts.maxIts;
  eps = ioData.ts.eps;
  maxTime = ioData.ts.maxTime;

  output = new char[strlen(ioData.ts.output) + 1];
  sprintf(output, "%s", ioData.ts.output);

  cfl = cfl0;
  residual = 1.0;

}

//------------------------------------------------------------------------------

TsParameters::~TsParameters()
{

  if (output) delete [] output;

}

//------------------------------------------------------------------------------

void TsParameters::computeCflNumber(int its, double res)
{

  cfl = min( max( max(cflCoef1, max(cflCoef2*its, cflCoef3*its*its)), cfl0/pow(res,ser) ), (cflMax + cflMaxSlope*its) );

}

//------------------------------------------------------------------------------

// Included (MB)
void TsParameters::rstVar(IoData &ioData) {

  maxTime = ioData.ts.maxTime;

}

//------------------------------------------------------------------------------

