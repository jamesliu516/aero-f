#include <cstring>
#include <complex>
using std::complex;
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

#ifndef PI
#define PI 3.14159
#endif

//------------------------------------------------------------------------------

TsParameters::TsParameters(IoData &ioData)
{

  resType = ioData.ts.residual;

  cfllaw = ioData.ts.cfl.strategy;

  cfl0 = ioData.ts.cfl.cfl0;
  cflCoef1 = ioData.ts.cfl.cflCoef1;
  cflCoef2 = ioData.ts.cfl.cflCoef2;
  cflMax = ioData.ts.cfl.cflMax;
  cflMin = ioData.ts.cfl.cflMin;
  dualtimecfl = ioData.ts.cfl.dualtimecfl;

  checksol = !(!ioData.ts.adaptivetime.checksol || !ioData.ts.cfl.checksol || !ioData.ts.checksol);
  checklinsolve = !(!ioData.ts.cfl.checklinsolve || !ioData.ts.adaptivetime.checksol);
  checkriemann = checksol;
  checklargevelocity = ioData.ts.checkvelocity;
  checkpclipping = ioData.ts.checkpressure;
  rapidpchangecutoff = max(0,ioData.ts.deltapressurethreshold);

  ser = ioData.ts.cfl.ser;
  angle_growth = ioData.ts.cfl.angle_growth;
  angle_zero = ioData.ts.cfl.angle_zero;
  dft_history = ioData.ts.cfl.dft_history;
  dft_freqcutoff = ioData.ts.cfl.dft_freqcutoff;
  dft_growth = ioData.ts.cfl.dft_growth;
  fixedunsteady_counter = 0;

  maxIts = ioData.ts.maxIts;
  eps = ioData.ts.eps;
  maxTime = ioData.ts.maxTime;

  forbidreduce = ioData.ts.cfl.forbidreduce;

  output = new char[strlen(ioData.ts.output) + 1];
  sprintf(output, "%s", ioData.ts.output);

  cfl = cfl0;
  residual = 1.0;

  reshistory = new double[dft_history];
  for (int i=0; i<dft_history; i++) reshistory[i] = 1.0;
  dft = new complex<double>[dft_history];
  for (int i=0; i<dft_history; i++) dft[i]=0.0;

  unphysical = false;
  badlinsolve = false;
  allowstop = true;

  errorHandler = NULL;

}

//------------------------------------------------------------------------------

TsParameters::~TsParameters()
{

  if (output) delete [] output;

}

//------------------------------------------------------------------------------

//Figure out what to do with errors (related to time step)
void TsParameters::resolveErrors(){

  if (checklinsolve && errorHandler->globalErrors[ErrorHandler::SATURATED_LS]){
    errorHandler->com -> printf(1,"Detected saturated linear solver. Reducing time step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
  }

  if (checksol && errorHandler->globalErrors[ErrorHandler::UNPHYSICAL]){
    errorHandler->com -> printf(1,"Detected unphysical solution. Reducing time step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
  }
  
  if (checkriemann && errorHandler->globalErrors[ErrorHandler::BAD_RIEMANN]){
    errorHandler->com -> printf(1,"Detected error in Riemann solver. Reducing time step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
  }

  if (checkpclipping && errorHandler->globalErrors[ErrorHandler::PRESSURE_CLIPPING]){
    errorHandler->com -> printf(1,"Detected pressure clipping. Reducing error. If the problem is expected to produce cavitation, set CheckSolution to Off.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
  }

  if (checklargevelocity && errorHandler->globalErrors[ErrorHandler::LARGE_VELOCITY]){
    errorHandler->com -> printf(1,"Detected abnormally large velocities. Reducing time step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
  }

  if (rapidpchangecutoff && errorHandler->globalErrors[ErrorHandler::RAPIDLY_CHANGING_PRESSURE] >= rapidpchangecutoff){
    errorHandler->com -> printf(1,"Detected multiple rapidly changing pressures. Reducing time step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    //errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
  }

  errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP_TIME] = errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP];

}

//------------------------------------------------------------------------------

void TsParameters::computeCflNumber(int its, double res, double angle)
{

  // for backwards compatibility
  if (cfllaw == CFLData::OLD){
    cfl = min( max( max(cflCoef1, cflCoef2*its), cfl0/pow(res,ser) ), cflMax );
    return;
  }

  // First run automatic CFL checks
/*
  if (unphysical){
    unphysical=false;
    badlinsolve=false;
    cfl *= 0.5;
    fixedunsteady_counter = 1;
    //std::printf("Reduction params: cfl0=%f, cfl=%f\n",cfl0,cfl);
    //std::printf("Unphysicality detected. Reducing CFL number to %f.\n",cfl);
    if (cfl < cfl0/10000. && allowstop ) {
      std::printf("Could not resolve unphysicality by reducing CFL number. Aborting.\n"); 
      std::printf("Params: cfl0=%f, cfl=%f\n",cfl0,cfl);
      exit(-1);
    }
    return;
  }
  if (badlinsolve){
    unphysical=false;
    badlinsolve=false;
    cfl *= 0.5;
    fixedunsteady_counter = 1;
    //std::printf("Saturated linear solver detected. Reducing CFL number to %f.\n",cfl);
    if (cfl < cfl0/10000. && allowstop) {std::printf("Linear solver does not converge for any feasible CFL number. Aborting.\n"); exit(-1);}
    return;
  }
*/

  if (errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP]){
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP]=0;
    double cflold=cfl;
    cfl *= 0.5;
    fixedunsteady_counter = 1;
    errorHandler->com->printf(1,"Reducing CFL. Previous cfl=%e, new cfl=%e.\n",cflold,cfl);
    //std::printf("Saturated linear solver detected. Reducing CFL number to %f.\n",cfl);
    if (cfl < cfl0/10000. && allowstop) {std::printf("Cannot further reduce CFL number. Aborting.\n"); exit(-1);}
    return;
  }

  // Now run main code
  double cfl_res, cfl_dir, cfl_dft, hf_ratio;

  double cfl_prev = cfl;

  for (int i=dft_history-1; i>0; i--) reshistory[i]=reshistory[i-1];
  reshistory[0]=res;

  if (cfllaw == CFLData::RESIDUAL || cfllaw == CFLData::HYBRID){
    // compute residual strategy proposal
    cfl = cfl_prev;
    cfl *= pow(reshistory[1]/reshistory[0],ser);
    cfl_res = cfl;
    //std::printf("CFL residual strategy: old residual: %e, new residual: %e, CFL proposal: %e\n",reshistory[1],reshistory[0],cfl_res);
  }
  if (cfllaw == CFLData::DIRECTION || cfllaw == CFLData::HYBRID){
    // compute direction strategy proposal
    cfl = cfl_prev;
    if (angle != -2.0) cfl *= pow(angle_growth,angle-angle_zero);
    cfl_dir = cfl;
    //std::printf("CFL direction strategy: angle: %e, CFL proposal: %e\n",angle,cfl_dir);
  }
  if (cfllaw == CFLData::DFT || cfllaw == CFLData::HYBRID){
    // compute dft strategy proposal
    cfl = cfl_prev;

    if (reshistory[0] != 0.0){
      double e_hf, e_ac, e_dc, e_total;
      int cutofflow = dft_history/2 - (dft_freqcutoff-1)/2; 
      int cutoffhigh = (dft_history+1)/2 + (dft_freqcutoff-1)/2;

      for (int i=cutofflow; i<cutoffhigh+1; i++){
	dft[i]=0.0;
	for (int j=0; j<dft_history; j++) 
	  dft[i] += 1/sqrt((double)dft_history)*complex<double>(reshistory[j],0.0)*exp(-2.0*complex<double>(0.0,1.0)*PI*(double)i*(double)j/(double)dft_history);
      }
      e_hf = 0.0;
      for (int i=cutofflow; i<cutoffhigh+1; i++)
	e_hf += norm(dft[i]);

      e_total = 0.0;
      for (int i=0; i<dft_history; i++) e_total += pow(reshistory[i],2);
      e_dc = 0.0;
      for (int i=0; i<dft_history; i++) e_dc += reshistory[i]/sqrt((double)dft_history);
      e_dc *= e_dc;
      e_ac = e_total-e_dc;

      if (e_ac<e_dc*1e-5) hf_ratio = 1.0;
      else hf_ratio = e_hf/e_ac;

      if (hf_ratio<0 || hf_ratio>1) {
        std::printf("Found invalid hf_ratio: %e, check for bugs in CFL Law\n",hf_ratio); 
        char buffer[250];
        sprintf(buffer,"Res history:");
        for (int i=0; i<dft_history; i++) sprintf(buffer,"%s %e",buffer,reshistory[i]);
        sprintf(buffer,"%s \n",buffer);
        sprintf(buffer,"%s DFT:",buffer);
        for (int i=0; i<dft_history; i++) sprintf(buffer,"%s %e+%ei", buffer,dft[i].real(),dft[i].imag());
        sprintf(buffer,"%s \n",buffer);
        std::printf("%s e_total=%e, e_dc=%e, e_ac=%e, e_hf=%e, hf_ratio=%e\n",buffer,e_total,e_dc,e_ac,e_hf,hf_ratio);
        exit(-1);
      }

      cfl *= pow(dft_growth, 1-2*hf_ratio);
    }
    cfl_dft = cfl;
    //std::printf("CFL DFT strategy: e_ac: %e, e_hf: %e, CFL proposal: %e\n",e_ac,e_hf,cfl_dft);
  }
  if (cfllaw == CFLData::HYBRID){
    // compute hybrid strategy cfl
    if (hf_ratio > 0.66) cfl = cfl_dft;
    else if (angle < angle_zero) cfl = cfl_dir;
    else cfl = max(cfl_dir, cfl_res);
    //std::printf("CFL Hybrid strategy: dft_proposal: %e, direction proposal: %e, residual proposal: %e, chosen cfl: %e\n",cfl_dft,cfl_dir,cfl_res,cfl);
  }
  if (cfllaw == CFLData::FIXEDUNSTEADY){
    // compute fixed unsteady cfl law
    // attempt to keep cfl fixed, if automatic reductions lower the CFL, re-increase it after a few iterations
    cfl = cfl_prev;
    if(fixedunsteady_counter >= 4){
      fixedunsteady_counter = 2;
      cfl *= 2.0;
      cfl = min(cfl,cfl0);
    }
    fixedunsteady_counter++;
  }

  cfl = (min(max(cfl, cflCoef1),cflMax));
  if (forbidreduce && cfl<cfl_prev) cfl = cfl_prev;
  //std::printf("CFL number chosen: %e. Forbidreduce: %i. cfl_prev=%e\n",cfl,forbidreduce,cfl_prev);
  return;
}

//------------------------------------------------------------------------------

// Included (MB)
void TsParameters::rstVar(IoData &ioData) {

  maxTime = ioData.ts.maxTime;

}

//------------------------------------------------------------------------------

