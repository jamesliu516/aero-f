#include <string.h>

#include <TsInput.h>

#include <IoData.h>

//------------------------------------------------------------------------------

TsInput::TsInput(IoData &iod)
{

  int sp = strlen(iod.input.prefix) + 1;

  if (iod.input.solutions[0] != 0) {
    solutions = new char[sp + strlen(iod.input.solutions)];
    if (strncmp(iod.input.solutions, "/", 1) == 0)
      sprintf(solutions, "%s", iod.input.solutions);
    else
      sprintf(solutions, "%s%s", iod.input.prefix, iod.input.solutions);
  }
  else {
    solutions = new char[1];
    sprintf(solutions, "");
  }

  if (iod.input.positions[0] != 0) {
    positions = new char[sp + strlen(iod.input.positions)];
    if (strncmp(iod.input.positions, "/", 1) == 0)
      sprintf(positions, "%s", iod.input.positions);
    else
      sprintf(positions, "%s%s", iod.input.prefix, iod.input.positions);
  }
  else {
    positions = new char[1];
    sprintf(positions, "");
  }

  if (iod.input.levelsets[0] != 0) {
    levelsets = new char[sp + strlen(iod.input.levelsets)];
    if (strncmp(iod.input.levelsets, "/", 1) == 0)
      sprintf(levelsets, "%s", iod.input.levelsets);
    else
      sprintf(levelsets, "%s%s", iod.input.prefix, iod.input.levelsets);
  }
  else {
    levelsets = new char[1];
    sprintf(levelsets, "");
  }

  if (iod.input.podFile[0] != 0) {
    podFile = new char[sp + strlen(iod.input.podFile)];
    if (strncmp(iod.input.podFile, "/", 1) == 0)
      sprintf(podFile, "%s", iod.input.podFile);
    else
      sprintf(podFile, "%s%s", iod.input.prefix, iod.input.podFile);
  }
  else{
    podFile = new char[1];
    sprintf(podFile, "");
  }

  if (iod.input.snapFile[0] != 0) {
    snapFile = new char[sp + strlen(iod.input.snapFile)];
    if (strncmp(iod.input.snapFile, "/", 1) == 0)
      sprintf(snapFile, "%s", iod.input.snapFile);
    else
      sprintf(snapFile, "%s%s", iod.input.prefix, iod.input.snapFile);
  }
  else{
    snapFile = new char[1];
    sprintf(snapFile, "");
  }

  if (iod.input.sampleNodes[0] != 0) {
    sampleNodes = new char[sp + strlen(iod.input.sampleNodes)];
    if (strncmp(iod.input.sampleNodes, "/", 1) == 0)
      sprintf(sampleNodes, "%s", iod.input.sampleNodes);
    else
      sprintf(sampleNodes, "%s%s", iod.input.prefix, iod.input.sampleNodes);
  }
  else{
    sampleNodes = new char[1];
    sprintf(sampleNodes, "");
  }

  if (iod.input.podFileResJac[0] != 0) {
    podFileResJac = new char[sp + strlen(iod.input.podFileResJac)];
    if (strncmp(iod.input.podFileResJac, "/", 1) == 0)
      sprintf(podFileResJac, "%s", iod.input.podFileResJac);
    else
      sprintf(podFileResJac, "%s%s", iod.input.prefix, iod.input.podFileResJac);
  }
  else{
    podFileResJac = new char[1];
    sprintf(podFileResJac, "");
  }

  if (iod.input.aMatrix[0] != 0) {
    aMatrix = new char[sp + strlen(iod.input.aMatrix)];
    if (strncmp(iod.input.aMatrix, "/", 1) == 0)
      sprintf(aMatrix, "%s", iod.input.aMatrix);
    else
      sprintf(aMatrix, "%s%s", iod.input.prefix, iod.input.aMatrix);
  }
  else{
    aMatrix = new char[1];
    sprintf(aMatrix, "");
  }
  if (iod.input.bMatrix[0] != 0) {
    bMatrix = new char[sp + strlen(iod.input.bMatrix)];
    if (strncmp(iod.input.bMatrix, "/", 1) == 0)
      sprintf(bMatrix, "%s", iod.input.bMatrix);
    else
      sprintf(bMatrix, "%s%s", iod.input.prefix, iod.input.bMatrix);
  }
  else{
    bMatrix = new char[1];
    sprintf(bMatrix, "");
  }

// Included
  if (iod.input.shapederivatives[0] != 0) {
    shapederivatives = new char[sp + strlen(iod.input.shapederivatives)];
    if (strncmp(iod.input.shapederivatives, "/", 1) == 0)
      sprintf(shapederivatives, "%s", iod.input.shapederivatives);
    else
      sprintf(shapederivatives, "%s%s", iod.input.prefix, iod.input.shapederivatives);
  }
  else {
    shapederivatives = new char[1];
    sprintf(shapederivatives, "");
  }

}

//------------------------------------------------------------------------------
TsInput::~TsInput()
{

  if (solutions) delete [] solutions;
  if (positions) delete [] positions;
  if (levelsets) delete [] levelsets;
  if (podFile)   delete [] podFile;

// Included
  if (shapederivatives)   delete [] shapederivatives;

}
