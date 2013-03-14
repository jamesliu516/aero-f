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
  //snapFile  = absolutePath(iod.input.snapFile,  prefix);
  //snapRefSolutionFile = absolutePath(iod.input.snapRefSolutionFile, prefix);
  //podFileState = absolutePath(iod.input.podFileState, prefix);
  //podFileRes = absolutePath(iod.input.podFileRes, prefix);
  //podFileJac = absolutePath(iod.input.podFileJac, prefix);
  //podFileResHat = absolutePath(iod.input.podFileResHat, prefix);
  //podFileJacHat = absolutePath(iod.input.podFileJacHat, prefix);
  //sampleNodes = absolutePath(iod.input.sampleNodes, prefix);
  //jacMatrix = absolutePath(iod.input.jacMatrix, prefix);
  //resMatrix = absolutePath(iod.input.resMatrix, prefix);
	wallsurfacedisplac = absolutePath(iod.input.wallsurfacedisplac, prefix);
  shapederivatives = absolutePath(iod.input.shapederivatives, prefix);
  //staterom = absolutePath(iod.input.staterom, prefix);
  //reducedfullnodemap = absolutePath(iod.input.reducedfullnodemap , prefix);
  //mesh = absolutePath(iod.input.mesh , prefix);
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
 // delete[] snapFile;
 // delete[] snapRefSolutionFile;
 // delete[] podFileRes;
 // delete[] podFileJac;
 // delete[] podFileResHat;
 // delete[] podFileJacHat;
 // delete[] sampleNodes;
 // delete[] jacMatrix;
 // delete[] resMatrix;
	delete[] wallsurfacedisplac;
  delete[] shapederivatives; 
 // delete[] staterom;
 // delete[] reducedfullnodemap ;
 // delete[] mesh ;
}

//------------------------------------------------------------------------------
