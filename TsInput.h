#ifndef TS_INPUT_H
#define TS_INPUT_H

#include <string>

class IoData;

//------------------------------------------------------------------------------

struct TsInput {
  char *solutions;
  char *positions;
  char *levelsets;
  char *podFile;
  char *snapFile;
  char *snapRefSolutionFile;

// Gappy offline
  char *podFileResJac;
  char *podFileResJacHat;

// Gappy online
  char *sampleNodes;
  char *aMatrix;
  char *bMatrix;
  char *staterom;
  char *reducedfullnodemap;
  char *mesh;

  char *shapederivatives;

  TsInput(IoData &);
  ~TsInput();

private:
  static char * absolutePath(const std::string & rawPath, const std::string & prefix);

  // Disallow copy and assignment
  TsInput(const TsInput &);
  TsInput & operator=(const TsInput &);
};

//------------------------------------------------------------------------------

#endif
