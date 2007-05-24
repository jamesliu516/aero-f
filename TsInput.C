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
    sprintf(podFile, "%s%s", iod.input.prefix, iod.input.podFile);
  }
  else{
    podFile = new char[1];
    sprintf(podFile, "");
  }

}

//------------------------------------------------------------------------------
TsInput::~TsInput()
{

  if (solutions) delete [] solutions;
  if (positions) delete [] positions;
  if (levelsets) delete [] levelsets;
  if (podFile)   delete [] podFile;

}
