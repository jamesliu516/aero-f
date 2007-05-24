#include <TsParameters.h>
#include <IoData.h>

#include <math.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif

//------------------------------------------------------------------------------

TsParameters::TsParameters(IoData &ioData)
{

  resType = ioData.ts.residual;
  cfl0 = ioData.ts.cfl0;
  cflCoef1 = ioData.ts.cflCoef1;
  cflCoef2 = ioData.ts.cflCoef2;
  cflMax = ioData.ts.cflMax;
  ser = ioData.ts.ser;

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

  cfl = min( max( max(cflCoef1, cflCoef2*its), cfl0/pow(res,ser) ), cflMax );

}

//------------------------------------------------------------------------------

