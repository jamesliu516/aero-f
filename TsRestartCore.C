#include <cstring>

#include <TsRestart.h>

#include <IoData.h>
#include <RefVal.h>

//------------------------------------------------------------------------------

TsRestart::TsRestart(IoData &iod, RefVal *rv) : refVal(rv)
  , writeKPtraces(false)
{

  int sp = strlen(iod.output.restart.prefix) + 1;

  solutions[0] = new char[sp + strlen(iod.output.restart.solutions)];
  if (iod.output.restart.solutions[0] != 0)
    sprintf(solutions[0], "%s%s", iod.output.restart.prefix, iod.output.restart.solutions);
  else
    sprintf(solutions[0], "");

  positions[0] = new char[sp + strlen(iod.output.restart.positions)];
  if (iod.output.restart.positions[0] != 0 && (iod.problem.framework == ProblemData::BODYFITTED || iod.problem.framework == ProblemData::EMBEDDEDALE))
    sprintf(positions[0], "%s%s", iod.output.restart.prefix, iod.output.restart.positions);
  else
    sprintf(positions[0], "");

  levelsets[0] = new char[sp + strlen(iod.output.restart.levelsets)];
  if (iod.output.restart.levelsets[0] != 0)
    sprintf(levelsets[0], "%s%s", iod.output.restart.prefix, iod.output.restart.levelsets);
  else
    sprintf(levelsets[0], "");

  cracking[0] = new char[sp + strlen(iod.output.restart.cracking)];
  if (iod.output.restart.cracking[0] != 0)
    sprintf(cracking[0], "%s%s", iod.output.restart.prefix, iod.output.restart.cracking);
  else
    sprintf(cracking[0], "");


  data[0] = new char[sp + strlen(iod.output.restart.data)];
  if (iod.output.restart.data[0] != 0)
    sprintf(data[0], "%s%s", iod.output.restart.prefix, iod.output.restart.data);
  else
    sprintf(data[0], "");

  structPos = new char[sp + strlen(iod.output.restart.embeddedpositions)];
  if (iod.output.restart.embeddedpositions[0] != 0 && (iod.problem.framework == ProblemData::EMBEDDED || iod.problem.framework == ProblemData::EMBEDDEDALE) )
    sprintf(structPos, "%s%s", iod.output.restart.prefix, iod.output.restart.embeddedpositions);
  else
    sprintf(structPos, "");

  if (iod.output.restart.type == RestartData::SINGLE) {
    deleteCharStar = false;
    for (int i=1; i<3; ++i) {
      solutions[i] = solutions[0];
      positions[i] = positions[0];
      levelsets[i] = levelsets[0];
      cracking[i] = cracking[0];
      data[i] = data[0];

    }

    index = 0;
  }
  else {
    deleteCharStar = true;
    for (int i=1; i<3; ++i) {
      solutions[i] = new char[strlen(solutions[0]) + 5 + 1];
      if (solutions[0][0] != 0)
	sprintf(solutions[i], "%s.%drst", solutions[0], i);
      else
	sprintf(solutions[i], "");

      positions[i] = new char[strlen(positions[0]) + 5 + 1];
      if (positions[0][0] != 0)
	sprintf(positions[i], "%s.%drst", positions[0], i);
      else
	sprintf(positions[i], "");

      levelsets[i] = new char[strlen(levelsets[0]) + 5 + 1];
      if (levelsets[0][0] != 0)
	sprintf(levelsets[i], "%s.%drst", levelsets[0], i);
      else
	sprintf(levelsets[i], "");

      cracking[i] = new char[strlen(cracking[0]) + 5 + 1];
      if (cracking[0][0] != 0)
	sprintf(cracking[i], "%s.%drst", cracking[0], i);
      else
	sprintf(cracking[i], "");


      data[i] = new char[strlen(data[0]) + 5 + 1];
      if (data[0][0] != 0)
	sprintf(data[i], "%s.%drst", data[0], i);
      else
	sprintf(data[i], "");

    }

    index = 1;
  }

  iteration = iod.restart.iteration;
  etime = iod.restart.etime;
  residual = iod.restart.residual;
  energy[0] = iod.restart.energy;
  energy[1] = iod.restart.energy;
  frequency = iod.output.restart.frequency;
  frequency_dt = iod.output.restart.frequency_dt;
  prtout = 0.0;

  //
  // Check whether pressure snapshots are needed.
  // UH (07/2012)
  //
  if (strlen(iod.output.restart.strKPtraces) > 0)
  {
    writeKPtraces = true;
  }

}

//------------------------------------------------------------------------------

bool TsRestart::toWrite(int it, bool lastIt, double t)
{
  if(frequency_dt<=0.0)
    return (((frequency > 0) && (it % frequency == 0)) || lastIt);

  return (t>=prtout || lastIt);
}

//------------------------------------------------------------------------------

void TsRestart::updatePrtout(double t)
{
  if(frequency_dt<=0.0)
    return;
  if(t>=prtout)
    prtout += frequency_dt;
}

//------------------------------------------------------------------------------

void TsRestart::writeStructPosToDisk(int cpuNum, bool lastIt, Vec<Vec3D>& Xs)
{
  
  if(cpuNum>0) return; //only Proc.#1 will work.

  if (structPos[0]!=0 && toWrite(iteration, lastIt, etime)) {
  //if ((lastIt || (frequency > 0 && iteration % frequency == 0)) && structPos[0]!=0) {
    FILE *fp = fopen(structPos,"w");
    if (!fp) {
      fprintf(stderr, "*** Error: could not open \'%s\'\n", structPos);
      exit(1);
    }

    for (int i=0; i<Xs.size(); i++)
      fprintf(fp,"%d %e %e %e\n", i+1, Xs[i][0], Xs[i][1], Xs[i][2]);

   fclose(fp);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void TsRestart::rstVar(IoData &ioData) {

  etime = ioData.restart.etime;
  energy[0] = ioData.restart.energy;
  energy[1] = ioData.restart.energy;

}

//------------------------------------------------------------------------------
