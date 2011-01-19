#include <TsInput.h>

#include <IoData.h>

#include <string>
#include <cstring>

//------------------------------------------------------------------------------

TsInput::TsInput(IoData &iod) {
  const std::string prefix(iod.input.prefix); 

  solutions = absolutePath(iod.input.solutions, prefix);
  positions = absolutePath(iod.input.positions, prefix);
  levelsets = absolutePath(iod.input.levelsets, prefix);
  podFile   = absolutePath(iod.input.podFile,   prefix);
  snapFile  = absolutePath(iod.input.snapFile,  prefix);
  snapRefSolutionFile = absolutePath(iod.input.snapRefSolutionFile, prefix);
  podFileResJac = absolutePath(iod.input.podFileResJac, prefix);
  sampleNodes = absolutePath(iod.input.sampleNodes, prefix);
  aMatrix = absolutePath(iod.input.aMatrix, prefix);
  bMatrix = absolutePath(iod.input.bMatrix, prefix);
  shapederivatives = absolutePath(iod.input.shapederivatives, prefix);
}

//------------------------------------------------------------------------------

char *
TsInput::absolutePath(const std::string & rawPath, const std::string & prefix) {
  const bool isRelativePath = rawPath.size() > 0 && rawPath[0] != '/';
  const std::string finalPath(isRelativePath ? prefix + rawPath : rawPath);

  char * result = new char[finalPath.size() + 1];
  std::strcpy(result, finalPath.c_str());

  return result;
}

//------------------------------------------------------------------------------

TsInput::~TsInput() {
  delete[] solutions;
  delete[] positions;
  delete[] levelsets;
  delete[] podFile;
  delete[] snapFile;
  delete[] snapRefSolutionFile;
  delete[] podFileResJac;
  delete[] sampleNodes;
  delete[] aMatrix;
  delete[] bMatrix;
  delete[] shapederivatives; 
}

//------------------------------------------------------------------------------
