#include <stdlib.h>
#include <string.h>

#include <TsOutput.h>

#include <IoData.h>
#include <RefVal.h>
#include <Domain.h>
#include <PostOperator.h>
#include <MeshMotionHandler.h>
#include <DistVector.h>
#include <BinFileHandler.h>

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::TsOutput(IoData &iod, RefVal *rv, Domain *dom, PostOperator<dim> *po) : 
  refVal(rv), domain(dom), postOp(po), rmmh(0)
{

  steady = !iod.problem.type[ProblemData::UNSTEADY];
  com = domain->getCommunicator();

  int sp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.solutions[0] != 0) {
    solutions = new char[sp + strlen(iod.output.transient.solutions)];
    sprintf(solutions, "%s%s", iod.output.transient.prefix, iod.output.transient.solutions);
  }
  else
    solutions = 0;

  int i;
  for (i=0; i<PostFcn::SSIZE; ++i) {
    sscale[i] = 1.0;
    scalars[i] = 0;
  }
  for (i=0; i<PostFcn::AVSSIZE; ++i) {
    avsscale[i] = 1.0;
    avscalars[i] = 0;
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    vscale[i] = 1.0;
    vectors[i] = 0;
  }
  for (i=0; i<PostFcn::AVVSIZE; ++i) {
    avvscale[i] = 1.0;
    avvectors[i] = 0;
  }

  if (iod.output.transient.density[0] != 0) {
    sscale[PostFcn::DENSITY] = iod.ref.rv.density;
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.density);
  }
  if (iod.output.transient.tavdensity[0] != 0) {
    avsscale[PostFcn::DENSITYAVG] = iod.ref.rv.density;
    avscalars[PostFcn::DENSITYAVG] = new char[sp + strlen(iod.output.transient.tavdensity)];
    sprintf(avscalars[PostFcn::DENSITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavdensity);
  }
  if (iod.output.transient.mach[0] != 0) {
    scalars[PostFcn::MACH] = new char[sp + strlen(iod.output.transient.mach)];
    sprintf(scalars[PostFcn::MACH], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.mach);
  }
  if (iod.output.transient.wtmach[0] != 0) {
    scalars[PostFcn::WTMACH] = new char[sp + strlen(iod.output.transient.wtmach)];
    sprintf(scalars[PostFcn::WTMACH], "%s%s", iod.output.transient.prefix, iod.output.transient.wtmach);
  }
  if (iod.output.transient.wtspeed[0] != 0) {
    scalars[PostFcn::WTSPEED] = new char[sp + strlen(iod.output.transient.wtmach)];
    sprintf(scalars[PostFcn::WTSPEED], "%s%s", iod.output.transient.prefix, iod.output.transient.wtspeed);
    sscale[PostFcn::WTSPEED] = iod.ref.rv.velocity;
  }
  if (iod.output.transient.speed[0] != 0) {
    sscale[PostFcn::SPEED] = iod.ref.rv.velocity;
    scalars[PostFcn::SPEED] = new char[sp + strlen(iod.output.transient.speed)];
    sprintf(scalars[PostFcn::SPEED], "%s%s",
            iod.output.transient.prefix, iod.output.transient.speed);
  }
  if (iod.output.transient.tavmach[0] != 0) {
    avscalars[PostFcn::MACHAVG] = new char[sp + strlen(iod.output.transient.tavmach)];
    sprintf(avscalars[PostFcn::MACHAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavmach);
  }
  if (iod.output.transient.pressure[0] != 0) {
    sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.pressure);
  }
  if (iod.output.transient.diffpressure[0] != 0) {
    sscale[PostFcn::DIFFPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::DIFFPRESSURE] = new char[sp + strlen(iod.output.transient.diffpressure)];
    sprintf(scalars[PostFcn::DIFFPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.diffpressure);
  }
  if (iod.output.transient.tavpressure[0] != 0) {
    avsscale[PostFcn::PRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::PRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavpressure)];
    sprintf(avscalars[PostFcn::PRESSUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavpressure);
  }
  if (iod.output.transient.hydrostaticpressure[0] != 0) {
    sscale[PostFcn::HYDROSTATICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDROSTATICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrostaticpressure)];
    sprintf(scalars[PostFcn::HYDROSTATICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrostaticpressure);
  }
  if (iod.output.transient.hydrodynamicpressure[0] != 0) {
    sscale[PostFcn::HYDRODYNAMICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDRODYNAMICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrodynamicpressure)];
    sprintf(scalars[PostFcn::HYDRODYNAMICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrodynamicpressure);
  }
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.temperature);
  }
  if (iod.output.transient.tavtemperature[0] != 0) {
    avsscale[PostFcn::TEMPERATUREAVG] = iod.ref.rv.temperature;
    avscalars[PostFcn::TEMPERATUREAVG] = new char[sp + strlen(iod.output.transient.tavtemperature)];
    sprintf(avscalars[PostFcn::TEMPERATUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavtemperature);
  }
  if (iod.output.transient.totalpressure[0] != 0) {
    sscale[PostFcn::TOTPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::TOTPRESSURE] = new char[sp + strlen(iod.output.transient.totalpressure)];
    sprintf(scalars[PostFcn::TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.totalpressure);
  }
  if (iod.output.transient.tavtotalpressure[0] != 0) {
    avsscale[PostFcn::TOTPRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::TOTPRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavtotalpressure)];
    sprintf(avscalars[PostFcn::TOTPRESSUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavtotalpressure);
  }
  if (iod.output.transient.vorticity[0] != 0) {
    sscale[PostFcn::VORTICITY] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    scalars[PostFcn::VORTICITY] = new char[sp + strlen(iod.output.transient.vorticity)];
    sprintf(scalars[PostFcn::VORTICITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.vorticity);
  }
  if (iod.output.transient.tavvorticity[0] != 0) {
    avsscale[PostFcn::VORTICITYAVG] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    avscalars[PostFcn::VORTICITYAVG] = new char[sp + strlen(iod.output.transient.tavvorticity)];
    sprintf(avscalars[PostFcn::VORTICITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavvorticity);
  }
  if (iod.output.transient.nutturb[0] != 0) {
    sscale[PostFcn::NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    scalars[PostFcn::NUT_TURB] = new char[sp + strlen(iod.output.transient.nutturb)];
    sprintf(scalars[PostFcn::NUT_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.nutturb);
  }
  if (iod.output.transient.kturb[0] != 0) {
    sscale[PostFcn::K_TURB] = iod.ref.rv.kenergy;
    scalars[PostFcn::K_TURB] = new char[sp + strlen(iod.output.transient.kturb)];
    sprintf(scalars[PostFcn::K_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.kturb);
  }
  if (iod.output.transient.epsturb[0] != 0) {
    sscale[PostFcn::EPS_TURB] = iod.ref.rv.epsilon;
    scalars[PostFcn::EPS_TURB] = new char[sp + strlen(iod.output.transient.epsturb)];
    sprintf(scalars[PostFcn::EPS_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.epsturb);
  }
  if (iod.output.transient.eddyvis[0] != 0) {
    sscale[PostFcn::EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    scalars[PostFcn::EDDY_VISCOSITY] = new char[sp + strlen(iod.output.transient.eddyvis)];
    sprintf(scalars[PostFcn::EDDY_VISCOSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.eddyvis);
  }
  if (iod.output.transient.dplus[0] != 0) {
#if defined(HEAT_FLUX)
    /*
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    double dT = iod.bc.wall.temperature - 1.0 / (gam*(gam-1.0)*iod.bc.inlet.mach*iod.bc.inlet.mach);
    sscale[PostFcn::DELTA_PLUS] = iod.ref.reynolds * iod.eqs.thermalCondModel.prandtl / (gam * dT);
    */
    sscale[PostFcn::DELTA_PLUS] = iod.ref.rv.tpower / (iod.ref.length*iod.ref.length);
#endif
    scalars[PostFcn::DELTA_PLUS] = new char[sp + strlen(iod.output.transient.dplus)];
    sprintf(scalars[PostFcn::DELTA_PLUS], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dplus);
  }
  if (iod.output.transient.psensor[0] != 0) {
    scalars[PostFcn::PSENSOR] = new char[sp + strlen(iod.output.transient.psensor)];
    sprintf(scalars[PostFcn::PSENSOR], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.psensor);
  }
  if (iod.output.transient.csdles[0] != 0) {
    scalars[PostFcn::CSDLES] = new char[sp + strlen(iod.output.transient.csdles)];
    sprintf(scalars[PostFcn::CSDLES], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.csdles);
  }
  if (iod.output.transient.csdvms[0] != 0) {
    scalars[PostFcn::CSDVMS] = new char[sp + strlen(iod.output.transient.csdvms)];
    sprintf(scalars[PostFcn::CSDVMS], "%s%s",
            iod.output.transient.prefix, iod.output.transient.csdvms);
  }
  if (iod.output.transient.mutOmu[0] != 0) {
    scalars[PostFcn::MUT_OVER_MU] = new char[sp + strlen(iod.output.transient.mutOmu)];
    sprintf(scalars[PostFcn::MUT_OVER_MU], "%s%s",
            iod.output.transient.prefix, iod.output.transient.mutOmu);
  }
  if (iod.output.transient.philevel[0] != 0) {
    sscale[PostFcn::PHILEVEL] = 1.0;
    scalars[PostFcn::PHILEVEL] = new char[sp + strlen(iod.output.transient.philevel)];
    sprintf(scalars[PostFcn::PHILEVEL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel);
  }
  if (iod.output.transient.velocity[0] != 0) {
    vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocity);
  }
  if (iod.output.transient.tavvelocity[0] != 0) {
    avvscale[PostFcn::VELOCITYAVG] = iod.ref.rv.velocity;
    avvectors[PostFcn::VELOCITYAVG] = new char[sp + strlen(iod.output.transient.tavvelocity)];
    sprintf(avvectors[PostFcn::VELOCITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavvelocity);
  }
  if (iod.output.transient.displacement[0] != 0) {
    vscale[PostFcn::DISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::DISPLACEMENT] = new char[sp + strlen(iod.output.transient.displacement)];
    sprintf(vectors[PostFcn::DISPLACEMENT], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.displacement);
  }
  if (iod.output.transient.tavdisplacement[0] != 0) {
    avvscale[PostFcn::DISPLACEMENTAVG] = iod.ref.rv.tlength;
    avvectors[PostFcn::DISPLACEMENTAVG] = new char[sp + strlen(iod.output.transient.tavdisplacement)];
    sprintf(avvectors[PostFcn::DISPLACEMENTAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavdisplacement);
  }

  if (iod.output.transient.flightDisplacement[0] != 0) {
    vscale[PostFcn::FLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::FLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.flightDisplacement)];
    sprintf(vectors[PostFcn::FLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.flightDisplacement);
  }

  if (iod.output.transient.localFlightDisplacement[0] != 0) {
    vscale[PostFcn::LOCALFLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::LOCALFLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.localFlightDisplacement)];
    sprintf(vectors[PostFcn::LOCALFLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.localFlightDisplacement);
  }

  if (iod.output.transient.forces[0] != 0) {
    forces = new char[sp + strlen(iod.output.transient.forces)];
    sprintf(forces, "%s%s", iod.output.transient.prefix, iod.output.transient.forces);
  }
  else
    forces = 0;

  if (iod.output.transient.tavforces[0] != 0) {
    tavforces = new char[sp + strlen(iod.output.transient.tavforces)];
    sprintf(tavforces, "%s%s", iod.output.transient.prefix, iod.output.transient.tavforces);
  }
  else
    tavforces = 0;

  if (iod.output.transient.hydrostaticforces[0] != 0) {
    hydrostaticforces = new char[sp + strlen(iod.output.transient.hydrostaticforces)];
    sprintf(hydrostaticforces, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrostaticforces);
  }
  else
    hydrostaticforces = 0;

  if (iod.output.transient.hydrodynamicforces[0] != 0) {
    hydrodynamicforces = new char[sp + strlen(iod.output.transient.hydrodynamicforces)];
    sprintf(hydrodynamicforces, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrodynamicforces);
  }
  else
    hydrodynamicforces = 0;


  if (iod.output.transient.lift[0] != 0){
    lift = new char[sp + strlen(iod.output.transient.lift)];
    sprintf(lift, "%s%s", iod.output.transient.prefix, iod.output.transient.lift);
  }
  else
    lift = 0;

  if (iod.output.transient.tavlift[0] != 0){
    tavlift = new char[sp + strlen(iod.output.transient.tavlift)];
    sprintf(tavlift, "%s%s", iod.output.transient.prefix, iod.output.transient.tavlift);
  }
  else
    tavlift = 0;

  if (iod.output.transient.hydrostaticlift[0] != 0) {
    hydrostaticlift = new char[sp + strlen(iod.output.transient.hydrostaticlift)];
    sprintf(hydrostaticlift, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrostaticlift);
  }
  else
    hydrostaticlift = 0;
                                                                                                                                                                                                     
  if (iod.output.transient.hydrodynamiclift[0] != 0) {
    hydrodynamiclift = new char[sp + strlen(iod.output.transient.hydrodynamiclift)];
    sprintf(hydrodynamiclift, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrodynamiclift);
  }
  else
    hydrodynamiclift = 0;

  if (iod.output.transient.residuals[0] != 0) {
    residuals = new char[sp + strlen(iod.output.transient.residuals)];
    sprintf(residuals, "%s%s", iod.output.transient.prefix, iod.output.transient.residuals);
  }
  else
    residuals = 0;

  it0 = iod.restart.iteration;
  frequency = iod.output.transient.frequency;
  length = iod.output.transient.length;
  surface = iod.output.transient.surface;
  x0[0] = iod.output.transient.x0;
  x0[1] = iod.output.transient.y0;
  x0[2] = iod.output.transient.z0;

  Qs = 0;
  Qv = 0;
  fpResiduals = 0;

  int nSurf = postOp->getNumSurf();
  fpForces             = new FILE *[nSurf];
  fpLift               = new FILE *[nSurf];
  fpTavForces          = new FILE *[nSurf];
  fpTavLift            = new FILE *[nSurf];
  fpHydroStaticForces  = new FILE *[nSurf];
  fpHydroDynamicForces = new FILE *[nSurf];
  fpHydroStaticLift    = new FILE *[nSurf];
  fpHydroDynamicLift   = new FILE *[nSurf];
  for (int iSurf = 0; iSurf < nSurf; iSurf++)  {
    fpForces[iSurf]          = 0;
    fpLift[iSurf]            = 0;
    fpTavForces[iSurf]       = 0;
    fpTavLift[iSurf]         = 0;
    fpHydroStaticForces[iSurf]  = 0;
    fpHydroDynamicForces[iSurf] = 0;
    fpHydroStaticLift[iSurf]    = 0;
    fpHydroDynamicLift[iSurf]   = 0;
  }

  TavF = 0;
  TavM = 0;
  TavL = 0;

// Included (MB)
  if (iod.output.transient.velocitynorm[0] != 0) {
    sscale[PostFcn::VELOCITY_NORM] = iod.ref.rv.velocity;
    scalars[PostFcn::VELOCITY_NORM] = new char[sp + strlen(iod.output.transient.velocitynorm)];
    sprintf(scalars[PostFcn::VELOCITY_NORM], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocitynorm);
  }

  if (iod.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {

  int dsp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.dSolutions[0] != 0) {
    dSolutions = new char[dsp + strlen(iod.output.transient.dSolutions)];
    sprintf(dSolutions, "%s%s", iod.output.transient.prefix, iod.output.transient.dSolutions);
  }
  else
    dSolutions = 0;

  int i;
  for (i=0; i<PostFcn::DSSIZE; ++i) {
    dSscale[i] = 1.0;
    dScalars[i] = 0;
  }
  for (i=0; i<PostFcn::DVSIZE; ++i) {
    dVscale[i] = 1.0;
    dVectors[i] = 0;
  }

  if (iod.output.transient.dDensity[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_DENSITY] = iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_DENSITY] = new char[dsp + strlen(iod.output.transient.dDensity)];
    sprintf(dScalars[PostFcn::DERIVATIVE_DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dDensity);
  }
  if (iod.output.transient.dMach[0] != 0) {
    dScalars[PostFcn::DERIVATIVE_MACH] = new char[dsp + strlen(iod.output.transient.dMach)];
    sprintf(dScalars[PostFcn::DERIVATIVE_MACH], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dMach);
  }
  if (iod.output.transient.dPressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_PRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_PRESSURE] = new char[dsp + strlen(iod.output.transient.dPressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dPressure);
  }
  if (iod.output.transient.dTemperature[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TEMPERATURE] = iod.ref.rv.temperature;
    dScalars[PostFcn::DERIVATIVE_TEMPERATURE] = new char[dsp + strlen(iod.output.transient.dTemperature)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTemperature);
  }
  if (iod.output.transient.dTotalpressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TOTPRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_TOTPRESSURE] = new char[dsp + strlen(iod.output.transient.dTotalpressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTotalpressure);
  }
  if (iod.output.transient.dNutturb[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_NUT_TURB] = new char[dsp + strlen(iod.output.transient.dNutturb)];
    sprintf(dScalars[PostFcn::DERIVATIVE_NUT_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dNutturb);
  }
  if (iod.output.transient.dEddyvis[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = new char[dsp + strlen(iod.output.transient.dEddyvis)];
    sprintf(dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dEddyvis);
  }
  if (iod.output.transient.dVelocityScalar[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = iod.ref.rv.velocity;
    dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = new char[dsp + strlen(iod.output.transient.dVelocityScalar)];
    sprintf(dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dVelocityScalar);
  }
  if (iod.output.transient.dVelocityVector[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = iod.ref.rv.velocity;
    dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = new char[dsp + strlen(iod.output.transient.dVelocityVector)];
    sprintf(dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dVelocityVector);
  }
  if (iod.output.transient.dDisplacement[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_DISPLACEMENT] = iod.ref.rv.tlength;
    dVectors[PostFcn::DERIVATIVE_DISPLACEMENT] = new char[dsp + strlen(iod.output.transient.dDisplacement)];
    sprintf(dVectors[PostFcn::DERIVATIVE_DISPLACEMENT], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dDisplacement);
  }

  if (iod.output.transient.dForces[0] != 0) {
    dForces = new char[dsp + strlen(iod.output.transient.dForces)];
    sprintf(dForces, "%s%s", iod.output.transient.prefix, iod.output.transient.dForces);
  }
  else
    dForces = 0;

  fpdForces = 0;

  switchOpt = true;
  }
  else if (iod.problem.alltype != ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {
    switchOpt = false;
  }

}

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::~TsOutput()
{

  if (Qs) delete Qs;
  if (Qv) delete Qv;
  if(TavF) delete [] TavF;
  if(TavM) delete [] TavM;
  if(TavL) delete [] TavL;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::setMeshMotionHandler(RigidMeshMotionHandler *mmh)
{

  rmmh = mmh;

}

//------------------------------------------------------------------------------

template<int dim>
FILE *TsOutput<dim>::backupAsciiFile(char *name)
{

  FILE *fp = fopen(name, "r");
  if (!fp) 
    return 0;

  char nameback[MAXLINE], line[MAXLINE];
  sprintf(nameback, "%s.back", name);
  FILE *fpback = fopen(nameback, "w");
  if (!fpback) {
    fprintf(stderr, "*** Error: could not open backup file \'%s\'\n", nameback);
    return 0;
  }
  while (fgets(line, MAXLINE, fp) != 0)
    fprintf(fpback, "%s", line);
  fclose(fp);
  fclose(fpback);

  fpback = fopen(nameback, "r");
  fp = fopen(name, "w");
  fgets(line, MAXLINE, fpback);
  fprintf(fp, "%s", line);
  while (fgets(line, MAXLINE, fpback) != 0) {
    int iter;
    sscanf(line, "%d", &iter);
    if (iter <= it0) 
      fprintf(fp, "%s", line);
  }
  fclose(fpback);

  return fp;

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::openAsciiFiles()
{

  if (com->cpuNum() != 0) return;

  int nSurf = postOp->getNumSurf();
  int *surfNums = 0;
  int iSurf;
  if (nSurf > 0)  {
    surfNums = new int[nSurf];
    map<int, int> surfMap = postOp->getSurfMap();
    map<int, int>::iterator it;
    iSurf = 1;
    surfNums[0] = 0;
    for (it = surfMap.begin(); it != surfMap.end(); it++)
      if (it->second > 0)
        surfNums[iSurf++] = it->first;
  }
      
  if (forces) {
    if (it0 != 0) 
      fpForces[0] = backupAsciiFile(forces);
    if (it0 == 0 || fpForces[0] == 0) {
      fpForces[0] = fopen(forces, "w");
      if (!fpForces[0]) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", forces);
	exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpForces[0], "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
	fprintf(fpForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
	fprintf(fpForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpForces[0]);
  } 

  if (forces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", forces, surfNums[iSurf]);
      if (it0 != 0)
        fpForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpForces[iSurf] == 0) {
        fpForces[iSurf] = fopen(filename, "w");
        if (!fpForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpForces[iSurf], "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }   
      fflush(fpForces[iSurf]);
    }
  }

// Included (MB)
  if (switchOpt) {
    if (dForces) {
      fpdForces = fopen(dForces, "w");
      if (!fpdForces) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dForces);
        exit(1);
      }

      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdForces, "Step Variable Cfx Cfy Cfz Cmx Cmy Cmz sboom dCfx dCfy dCfz dCmx dCmy dCmz dSboom \n");
      else
        fprintf(fpdForces, "Step Variable Fx Fy Fz Mx My Mz sboom dFx dFy dFz dMx dMy dMz dSboom \n");

      fflush(fpdForces);
    }
  }

  if (hydrostaticforces) {
    if (it0 != 0){
      fpHydroStaticForces[0] = backupAsciiFile(hydrostaticforces);
    }
    if (it0 == 0 || fpHydroStaticForces[0] == 0) {
      fpHydroStaticForces[0] = fopen(hydrostaticforces, "w");
      if (!fpHydroStaticForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrostaticforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHydroStaticForces[0], "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHydroStaticForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpHydroStaticForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpHydroStaticForces[0]);
  }

  if (hydrostaticforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", hydrostaticforces, surfNums[iSurf]);
      if (it0 != 0)
        fpHydroStaticForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroStaticForces[iSurf] == 0) {
        fpHydroStaticForces[iSurf] = fopen(filename, "w");
        if (!fpHydroStaticForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroStaticForces[iSurf], "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHydroStaticForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpHydroStaticForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpHydroStaticForces[iSurf]);
    }
  }

  if (hydrodynamicforces) {
    if (it0 != 0)
      fpHydroDynamicForces[0] = backupAsciiFile(hydrodynamicforces);
    if (it0 == 0 || fpHydroDynamicForces[0] == 0) {
      fpHydroDynamicForces[0] = fopen(hydrodynamicforces, "w");
      if (!fpHydroDynamicForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrodynamicforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHydroDynamicForces[0], "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHydroDynamicForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpHydroDynamicForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpHydroDynamicForces[0]);
  }
                                                                                                                                                                  
  if (hydrodynamicforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", hydrodynamicforces, surfNums[iSurf]);
      if (it0 != 0)
        fpHydroDynamicForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroDynamicForces[iSurf] == 0) {
        fpHydroDynamicForces[iSurf] = fopen(filename, "w");
        if (!fpHydroDynamicForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroDynamicForces[iSurf], "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHydroDynamicForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpHydroDynamicForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpHydroDynamicForces[iSurf]);
    }
  }

  if (tavforces) {
    if (it0 != 0) 
      fpTavForces[0] = backupAsciiFile(tavforces);
    if (it0 == 0 || fpTavForces[0] == 0) {
      fpTavForces[0] = fopen(tavforces, "w");
      if (!fpTavForces[0]) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", tavforces);
	exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpTavForces[0], "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
	fprintf(fpTavForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
	fprintf(fpTavForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpTavForces[0]);
  }
      
  if (tavforces) { 
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", tavforces, surfNums[iSurf]);
      if (it0 != 0)
        fpTavForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpTavForces[iSurf] == 0) {
      	fpTavForces[iSurf] = fopen(filename, "w");
        if (!fpTavForces[iSurf]) {
	  fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
	  exit(1);
        }
	const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpTavForces[iSurf], "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpTavForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpTavForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpTavForces[iSurf]);
    }
  }

  if (lift) {
   if(it0 !=0)
      fpLift[0] = backupAsciiFile(lift);
   if (it0 == 0 || fpLift[0] == 0){
     fpLift[0] = fopen(lift,"w");
     if(!fpLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", lift);
       exit(1);
     }    
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpLift[0]);
  }

  if (lift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", lift, surfNums[iSurf]);
      if(it0 !=0)
        fpLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpLift[iSurf] == 0){
        fpLift[iSurf] = fopen(filename,"w");
        if (!fpLift[iSurf]) {
	  fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
	  exit(1);
        }             
        const char *addvar = "";
	if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpLift[iSurf]);
    }   
  }

  if (hydrostaticlift) {
   if(it0 !=0)
      fpHydroStaticLift[0] = backupAsciiFile(hydrostaticlift);
   if (it0 == 0 || fpHydroStaticLift[0] == 0){
     fpHydroStaticLift[0] = fopen(hydrostaticlift,"w");
     if(!fpHydroStaticLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrostaticlift);
       exit(1);
     }    
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpHydroStaticLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpHydroStaticLift[0]);
  }

  if (hydrostaticlift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", hydrostaticlift, surfNums[iSurf]);
      if(it0 !=0)
        fpHydroStaticLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroStaticLift[iSurf] == 0){
        fpHydroStaticLift[iSurf] = fopen(filename,"w");
        if (!fpHydroStaticLift[iSurf]) {
	  fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
	  exit(1);
        }             
        const char *addvar = "";
	if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroStaticLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpHydroStaticLift[iSurf]);
    }   
  }

  if (hydrodynamiclift) {
   if(it0 !=0)
      fpHydroDynamicLift[0] = backupAsciiFile(hydrodynamiclift);
   if (it0 == 0 || fpHydroDynamicLift[0] == 0){
     fpHydroDynamicLift[0] = fopen(hydrodynamiclift,"w");
     if(!fpHydroDynamicLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrodynamiclift);
       exit(1);
     }
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpHydroDynamicLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpHydroDynamicLift[0]);
  }
                                                                                                                                                                  
  if (hydrodynamiclift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", hydrodynamiclift, surfNums[iSurf]);
      if(it0 !=0)
        fpHydroDynamicLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroDynamicLift[iSurf] == 0){
        fpHydroDynamicLift[iSurf] = fopen(filename,"w");
        if (!fpHydroDynamicLift[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroDynamicLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpHydroDynamicLift[iSurf]);
    }
  }

  if(tavlift) {
    if (it0 != 0) 
      fpTavLift[0] = backupAsciiFile(tavlift);
   if (it0 == 0 || fpTavLift[0] == 0){
     fpTavLift[0] = fopen(tavlift,"w");
     if(!fpTavLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", tavlift);
       exit(1);
     }    
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpTavLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
   }
   fflush(fpTavLift[0]);
  }
  
  if(tavlift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", tavlift, surfNums[iSurf]);
      if(it0 !=0)
        fpTavLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpTavLift[iSurf] == 0){
        fpTavLift[iSurf] = fopen(filename,"w");
        if (!fpTavLift[iSurf]) {
	  fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
	  exit(1);
        }      
        const char *addvar = "";
	if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpTavLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }       
      fflush(fpTavLift[iSurf]);
    }
  }

  if (residuals) {
    if (it0 != 0) 
      fpResiduals = backupAsciiFile(residuals);
    if (it0 == 0 || fpResiduals == 0) {
      fpResiduals = fopen(residuals, "w");
      if (!fpResiduals) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", residuals);
	exit(1);
      }
      fprintf(fpResiduals, "# TimeIteration ElapsedTime RelativeResidual CflNumber\n");
    }
    fflush(fpResiduals);
  }

 delete [] surfNums;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::closeAsciiFiles()
{

  for (int iSurf = 0; iSurf < postOp->getNumSurf(); iSurf++)  {
    if (fpForces[iSurf]) fclose(fpForces[iSurf]);
    if (fpHydroDynamicForces[iSurf]) fclose(fpHydroDynamicForces[iSurf]);
    if (fpHydroStaticForces[iSurf]) fclose(fpHydroStaticForces[iSurf]);
    if (fpLift[iSurf]) fclose(fpLift[iSurf]);
    if (fpTavForces[iSurf]) fclose(fpTavForces[iSurf]);
    if (fpTavLift[iSurf]) fclose(fpTavLift[iSurf]);
  }
  if (fpResiduals) fclose(fpResiduals);

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeForcesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  double del_t;

  if (TavF ==0 && TavM == 0) {
    TavF = new Vec3D[nSurfs];
    TavM = new Vec3D[nSurfs];
    for(int i = 0; i < nSurfs; ++i) {
      TavF[i] = 0.0;
      TavM[i] = 0.0;
    }
  }
 

  if(counter == 0) {
    tinit = refVal->time*t;
    tprevf = refVal->time*t;
    tener = 0.0; tenerold = 0.0;
    del_t = 0.0;
  }

  double time = refVal->time * t;
  if (forces || tavforces)  {
    Vec3D rVec = x0;
    postOp->computeForceAndMoment(rVec, X, U, Fi, Mi, Fv, Mv);    
    Vec3D F = Fi[0] + Fv[0];
    Vec3D moment = Mi[0] + Mv[0];
    F *= refVal->force;
    moment *= refVal->energy;
  }

  int iSurf;
  if (fpForces[0] || fpTavForces[0]) {
 
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D F = Fi[iSurf] + Fv[iSurf];
      Vec3D M = Mi[iSurf] + Mv[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        F *= 2.0 * refVal->length*refVal->length / surface;
        M *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        F *= refVal->force;
        M *= refVal->energy;
      }
      double energy = refVal->energy * e[0];


      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e %e \n", 
	        it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy, tag);
      }
      else
        fprintf(fpForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e \n", 
       	        it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);

      fflush(fpForces[iSurf]);

      if(fpTavForces[iSurf] && counter == 0){
        fprintf(fpTavForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e \n", 
          it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      del_t = time - tprevf; 
      F *= del_t; TavF[iSurf] += F;
      M *= del_t; TavM[iSurf] += M;
      tener += energy*del_t; 

      if(fpTavForces[iSurf] && counter > 0){
        F = TavF[iSurf]; M = TavM[iSurf]; energy = tener; 
        F /= (time - tinit); M /= (time - tinit);
        energy /= (time - tinit); 
        fprintf(fpTavForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e \n", 
           it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      fflush(fpTavForces[iSurf]);
    }
  }

  if (!steady) {
    if (rmmh) {
      double tag = rmmh->getTagValue(t);
      com->printf(0, "It %d (%d,%d): Time=%g, Mach=%g, Elapsed Time=%.2e s\n", 
		  it, itSc, itNl, t*refVal->time, tag, cpu);
    }
    else
      com->printf(0, "It %d (%d,%d): Time=%g, Elapsed Time=%.2e s\n", 
		  it, itSc, itNl, t*refVal->time, cpu);
  }
  tprevf = time;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::writeDerivativeOfForcesToDisk(int it, int actvar, Vec3D & F, Vec3D & dF, Vec3D & M, Vec3D & dM, double &sboom, double &dSboom)
{

  if (fpdForces) {
    fprintf(fpdForces, "%d %d %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e \n", it, actvar, F[0], F[1], F[2], M[0], M[1], M[2], sboom, dF[0], dF[1], dF[2], dM[0], dM[1], dM[2], dSboom);
    fflush(fpdForces);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeHydroForcesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *FiS = new Vec3D[nSurfs];
  Vec3D *MiS = new Vec3D[nSurfs];
  Vec3D *FvS = new Vec3D[nSurfs];
  Vec3D *MvS = new Vec3D[nSurfs];
  Vec3D *FiD = new Vec3D[nSurfs];
  Vec3D *MiD = new Vec3D[nSurfs];
  Vec3D *FvD = new Vec3D[nSurfs];
  Vec3D *MvD = new Vec3D[nSurfs];

  double time = refVal->time * t;

  if (hydrostaticforces)
    postOp->computeForceAndMoment(x0, X, U, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamicforces)
    postOp->computeForceAndMoment(x0, X, U, FiD, MiD, FvD, MvD, 2);

  int iSurf;
  if (fpHydroStaticForces[0]) {
 
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D FS = FiS[iSurf] + FvS[iSurf];
      Vec3D MS = MiS[iSurf] + MvS[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        FS *= 2.0 * refVal->length*refVal->length / surface;
        MS *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        FS *= refVal->force;
        MS *= refVal->energy;
      }
      double energy = refVal->energy * e[0];

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroStaticForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e %e \n", 
	        it, time, itSc, itNl, FS[0], FS[1], FS[2], MS[0], MS[1], MS[2], energy, tag);
      }
      else
        fprintf(fpHydroStaticForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e \n", 
       	        it, time, itSc, itNl, FS[0], FS[1], FS[2], MS[0], MS[1], MS[2], energy);

      fflush(fpHydroStaticForces[iSurf]);
    }
  }

  if (fpHydroDynamicForces[0]) {
                                                                                                                                                                                                     
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D FD = FiD[iSurf] + FvD[iSurf];
      Vec3D MD = MiD[iSurf] + MvD[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        FD *= 2.0 * refVal->length*refVal->length / surface;
        MD *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        FD *= refVal->force;
        MD *= refVal->energy;
      }
      double energy = refVal->energy * e[0];
                                                                                                                                                                                                     
      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroDynamicForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e %e\n",
                it, time, itSc, itNl, FD[0], FD[1], FD[2], MD[0], MD[1], MD[2], energy, tag);
      }
      else
        fprintf(fpHydroDynamicForces[iSurf], "%d %e %d %d %e %e %e %e %e %e %e \n",
                it, time, itSc, itNl, FD[0], FD[1], FD[2], MD[0], MD[1], MD[2], energy);
                                                                                                                                                                                                     
      fflush(fpHydroDynamicForces[iSurf]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

// This routine outputs both the non-averaged and time-averaged values of the lift and drag 
// in ascii format 

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  double del_t;

  if (TavL == 0) {
  TavL = new Vec3D[nSurfs];
  for(int i = 0; i < nSurfs; ++i) {
      TavL[i] = 0.0;
  }
  }

  int iSurf;
  if(counter == 0) {
    tinit = refVal->time*t;
    tprevl = refVal->time*t;
    del_t = 0.0;
  }

  if (lift || tavlift)
    postOp->computeForceAndMoment(x0, X, U, Fi, Mi, Fv, Mv);    

  double time = refVal->time * t;

  if (fpLift[0] || fpTavLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D L;
      Vec3D F = Fi[iSurf] + Fv[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        F *= 2.0 * refVal->length*refVal->length / surface;
      else
       F *= refVal->force;

      L[0] = F[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
             F[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
             F[2]*sin(iod.bc.inlet.alpha);

      L[1] = -F[0]*sin(iod.bc.inlet.beta) + F[1]*cos(iod.bc.inlet.beta);

      L[2] = -F[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
              F[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              F[2]*cos(iod.bc.inlet.alpha); 

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpLift[iSurf], "%d %e %d %d %e %e %e %e \n", 
  	        it, time, itSc, itNl, L[0], L[1], L[2], tag);
      }
      else
        fprintf(fpLift[iSurf], "%d %e %d %d %e %e %e \n", 
      	        it, time, itSc, itNl, L[0], L[1], L[2]);
      fflush(fpLift[iSurf]);

      if(fpTavLift[iSurf] && counter == 0){
        fprintf(fpTavLift[iSurf], "%d %e %d %d %e %e %e \n", 
             it, time, itSc, itNl, L[0], L[1], L[2]);
      }

      del_t = time - tprevl; 
      L *= del_t;  TavL[iSurf] += L;

      if(fpTavLift[iSurf] && counter > 0){
        L = TavL[iSurf];
        L /= (time - tinit); 
        fprintf(fpTavLift[iSurf], "%d %e %d %d %e %e %e \n", 
	        it, time, itSc, itNl, L[0], L[1], L[2]);
      }
      fflush(fpTavLift[iSurf]);
    }
  }

  tprevl = time;

}
//------------------------------------------------------------------------------
                                                                                                                                                                                                     
template<int dim>
void TsOutput<dim>::writeHydroLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *FiS = new Vec3D[nSurfs];
  Vec3D *MiS = new Vec3D[nSurfs];
  Vec3D *FvS = new Vec3D[nSurfs];
  Vec3D *MvS = new Vec3D[nSurfs];
  Vec3D *FiD = new Vec3D[nSurfs];
  Vec3D *MiD = new Vec3D[nSurfs];
  Vec3D *FvD = new Vec3D[nSurfs];
  Vec3D *MvD = new Vec3D[nSurfs];

  double time = refVal->time * t;
                                                                                                                                                                                                     
  if (hydrostaticlift)
    postOp->computeForceAndMoment(x0, X, U, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamiclift)
    postOp->computeForceAndMoment(x0, X, U, FiD, MiD, FvD, MvD, 2);

  int iSurf;
  if (fpHydroStaticLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D LS;
      Vec3D FS = FiS[iSurf] + FvS[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        FS *= 2.0 * refVal->length*refVal->length / surface;
      else
       FS *= refVal->force;

      LS[0] = FS[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
              FS[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              FS[2]*sin(iod.bc.inlet.alpha);
                                                                                                                                                                                                     
      LS[1] = -FS[0]*sin(iod.bc.inlet.beta) + FS[1]*cos(iod.bc.inlet.beta);
                                                                                                                                                                                                     
      LS[2] = -FS[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
               FS[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
               FS[2]*cos(iod.bc.inlet.alpha);

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroStaticLift[iSurf], "%d %e %d %d %e %e %e %e \n",
                it, time, itSc, itNl, LS[0], LS[1], LS[2], tag);
      }
      else
        fprintf(fpHydroStaticLift[iSurf], "%d %e %d %d %e %e %e \n",
                it, time, itSc, itNl, LS[0], LS[1], LS[2]);
      fflush(fpHydroStaticLift[iSurf]);
    }
  }

  if (fpHydroDynamicLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D LD;
      Vec3D FD = FiD[iSurf] + FvD[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        FD *= 2.0 * refVal->length*refVal->length / surface;
      else
        FD *= refVal->force;
                                                                                                                                                                                                     
      LD[0] = FD[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
              FD[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              FD[2]*sin(iod.bc.inlet.alpha);
                                                                                                                                                                                                     
      LD[1] = -FD[0]*sin(iod.bc.inlet.beta) + FD[1]*cos(iod.bc.inlet.beta);
                                                                                                                                                                                                     
      LD[2] = -FD[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
               FD[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
               FD[2]*cos(iod.bc.inlet.alpha);
                                                                                                                                                                                                     
      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroDynamicLift[iSurf], "%d %e %d %d %e %e %e %e \n",
                it, time, itSc, itNl, LD[0], LD[1], LD[2], tag);
      }
      else
        fprintf(fpHydroDynamicLift[iSurf], "%d %e %d %d %e %e %e \n",
                it, time, itSc, itNl, LD[0], LD[1], LD[2]);
      fflush(fpHydroDynamicLift[iSurf]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeResidualsToDisk(int it, double cpu, double res, double cfl)
{

  if (com->cpuNum() != 0) return;

  if (fpResiduals) {
    fprintf(fpResiduals, "%d %e %e %e\n", it, cpu, res, cfl);
    fflush(fpResiduals);
  }

  if (steady)
    com->printf(0, "It %5d: Res = %e, Cfl = %e, Elapsed Time = %.2e s\n", it, res, cfl, cpu);

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X, 
					     DistVec<double> &A, DistSVec<double,dim> &U, DistTimeState<dim> *timeState)
{

  if (((frequency > 0) && (it % frequency == 0)) || lastIt) {
    int step = 0;
    if (frequency > 0) {
      step = it / frequency;
      if (lastIt && (it % frequency != 0))
	++step;
    }
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;

    if (solutions)
      domain->writeVectorToFile(solutions, step, tag, U);

    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {
	if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
	postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState);
	DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
	domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));
      }
    }
    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());

        if (static_cast<PostFcn::VectorType>(i) == PostFcn::FLIGHTDISPLACEMENT)  {

          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getFlightPositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Flight Displacement Output not available\n");
        }
        else if (static_cast<PostFcn::VectorType>(i) == PostFcn::LOCALFLIGHTDISPLACEMENT)  {
          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Local Flight Displacement Output not available\n");

        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::writeBinaryDerivativeOfVectorsToDisk(int it, int actvar, double dS[3], DistSVec<double,3> &X, DistSVec<double,3> &dX, DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistTimeState<dim> *timeState)
{

  int    step = it-1;
  double tag  = (double)actvar;

  if (dSolutions)
    domain->writeVectorToFile(dSolutions, step, tag, dU);
 
  int i;
  for (i=0; i<PostFcn::DSSIZE; ++i) {
    if (dScalars[i]) {
      if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
      postOp->computeDerivativeOfScalarQuantity(static_cast<PostFcn::ScalarDerivativeType>(i), dS, X, dX, U, dU, *Qs, timeState);
      DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
      domain->writeVectorToFile(dScalars[i], step, tag, Qs1, &(dSscale[i]));
    }
  }
  for (i=0; i<PostFcn::DVSIZE; ++i) {
    if (dVectors[i]) {        
      if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
      if (static_cast<PostFcn::VectorType>(i) != PostFcn::FLIGHTDISPLACEMENT && static_cast<PostFcn::VectorType>(i) != PostFcn::LOCALFLIGHTDISPLACEMENT)
        postOp->computeDerivativeOfVectorQuantity(static_cast<PostFcn::VectorDerivativeType>(i), X, dX, U, dU, *Qv);
      domain->writeVectorToFile(dVectors[i], step, tag, *Qv, &(dVscale[i]));
    } 
  }

}

//------------------------------------------------------------------------------
                                                                                                      
template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, 
                                             DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U,
                                             DistVec<double> &Phi)
{
                                                                                                      
  if (((frequency > 0) && (it % frequency == 0)) || lastIt) {
    int step = 0;
    if (frequency > 0) {
      step = it / frequency;
      if (lastIt && (it % frequency != 0))
        ++step;
    }
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;
                                                                                                      
    if (solutions)
      domain->writeVectorToFile(solutions, step, tag, U);
                                                                                                      
    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, *Qs, Phi);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));
      }
    }
    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        if (rmmh) {
          DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, Phi);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }
                                                                                                      
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeAvgVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U, DistTimeState<dim> *timeState)
{

// This routine outputs the time averaged values of the scalar and vector output files
// in binary format

  int i;
  static DistVec<double> *AvQs[PostFcn::AVSSIZE];
  static DistSVec<double,3> *AvQv[PostFcn::AVVSIZE];
  static double tprev,tinit;
  double time = refVal->time*t;
  double del_t;

  int step = 0;
  if (frequency > 0) {
    step = it / frequency;
    if (lastIt && (it % frequency != 0))
    ++step;
  }

  double tag;
  if (rmmh)
    tag = rmmh->getTagValue(t);
  else
    tag = t * refVal->time;


  if (counter == 0){
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(!AvQs[i]) AvQs[i] = new DistVec<double>(domain->getNodeDistInfo());
      *AvQs[i] = 0.0;
    }
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A,*Qs, timeState);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(avscalars[i], step, tag, Qs1, &(avsscale[i]));
      }
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(!AvQv[i]) AvQv[i] = new DistSVec<double,3>(domain->getNodeDistInfo());
      *AvQv[i] = 0.0;
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(avvectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        if (rmmh) {
          DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
        }
        else
        postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        domain->writeVectorToFile(avvectors[i], step, tag, *Qv, &(avvscale[i]));
      }
    }

    tinit = time;
    tprev = time;
  }


  if (counter > 0){
    del_t = time - tprev; 
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if (avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState);
        *Qs *= del_t;
        *AvQs[i] += *Qs; 
      }
    }
    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if (avvectors[i]) {
	if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
	if (rmmh) {
	  DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
	  postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
	}
	else
	  postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        *Qv *= del_t;
        *AvQv[i] += *Qv;
      }
    }
    tprev = time;
  }

  
  if (((frequency > 0) && (it % frequency == 0) && (counter > 0)) || lastIt) {
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        *Qs = *AvQs[i];
	*Qs *= (1.0/(time-tinit));
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(avscalars[i], step, tag, Qs1, &(avsscale[i]));
      }
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(avvectors[i]) {
	if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
	*Qv = *AvQv[i];
	*Qv *= (1.0/(time-tinit));
	domain->writeVectorToFile(avvectors[i], step, tag, *Qv, &(avvscale[i])); 
      }
    } 
  }

  if (lastIt) {
    for (i=0; i<PostFcn::AVSSIZE; ++i) delete AvQs[i];
    for (i=0; i<PostFcn::AVVSIZE; ++i) delete AvQv[i];
  }

  counter += 1; // increment the counter for keeping track of the averaging
 
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::rstVar(IoData &iod) {

  int sp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.density[0] != 0) {
    sscale[PostFcn::DENSITY] = iod.ref.rv.density;
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.density);
  }
  if (iod.output.transient.tavdensity[0] != 0) {
    avsscale[PostFcn::DENSITYAVG] = iod.ref.rv.density;
    avscalars[PostFcn::DENSITYAVG] = new char[sp + strlen(iod.output.transient.tavdensity)];
    sprintf(avscalars[PostFcn::DENSITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavdensity);
  }
  if (iod.output.transient.speed[0] != 0) {
    sscale[PostFcn::SPEED] = iod.ref.rv.velocity;
    scalars[PostFcn::SPEED] = new char[sp + strlen(iod.output.transient.speed)];
    sprintf(scalars[PostFcn::SPEED], "%s%s",
            iod.output.transient.prefix, iod.output.transient.speed);
  }
  if (iod.output.transient.pressure[0] != 0) {
    sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.pressure);
  }
  if (iod.output.transient.diffpressure[0] != 0) {
    sscale[PostFcn::DIFFPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::DIFFPRESSURE] = new char[sp + strlen(iod.output.transient.diffpressure)];
    sprintf(scalars[PostFcn::DIFFPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.diffpressure);
  }
  if (iod.output.transient.tavpressure[0] != 0) {
    avsscale[PostFcn::PRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::PRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavpressure)];
    sprintf(avscalars[PostFcn::PRESSUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavpressure);
  }
  if (iod.output.transient.hydrostaticpressure[0] != 0) {
    sscale[PostFcn::HYDROSTATICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDROSTATICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrostaticpressure)];
    sprintf(scalars[PostFcn::HYDROSTATICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrostaticpressure);
  }
  if (iod.output.transient.hydrodynamicpressure[0] != 0) {
    sscale[PostFcn::HYDRODYNAMICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDRODYNAMICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrodynamicpressure)];
    sprintf(scalars[PostFcn::HYDRODYNAMICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrodynamicpressure);
  }
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.temperature);
  }
  if (iod.output.transient.tavtemperature[0] != 0) {
    avsscale[PostFcn::TEMPERATUREAVG] = iod.ref.rv.temperature;
    avscalars[PostFcn::TEMPERATUREAVG] = new char[sp + strlen(iod.output.transient.tavtemperature)];
    sprintf(avscalars[PostFcn::TEMPERATUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavtemperature);
  }
  if (iod.output.transient.totalpressure[0] != 0) {
    sscale[PostFcn::TOTPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::TOTPRESSURE] = new char[sp + strlen(iod.output.transient.totalpressure)];
    sprintf(scalars[PostFcn::TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.totalpressure);
  }
  if (iod.output.transient.tavtotalpressure[0] != 0) {
    avsscale[PostFcn::TOTPRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::TOTPRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavtotalpressure)];
    sprintf(avscalars[PostFcn::TOTPRESSUREAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavtotalpressure);
  }
  if (iod.output.transient.vorticity[0] != 0) {
    sscale[PostFcn::VORTICITY] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    scalars[PostFcn::VORTICITY] = new char[sp + strlen(iod.output.transient.vorticity)];
    sprintf(scalars[PostFcn::VORTICITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.vorticity);
  }
  if (iod.output.transient.tavvorticity[0] != 0) {
    avsscale[PostFcn::VORTICITYAVG] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    avscalars[PostFcn::VORTICITYAVG] = new char[sp + strlen(iod.output.transient.tavvorticity)];
    sprintf(avscalars[PostFcn::VORTICITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavvorticity);
  }
  if (iod.output.transient.nutturb[0] != 0) {
    sscale[PostFcn::NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    scalars[PostFcn::NUT_TURB] = new char[sp + strlen(iod.output.transient.nutturb)];
    sprintf(scalars[PostFcn::NUT_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.nutturb);
  }
  if (iod.output.transient.eddyvis[0] != 0) {
    sscale[PostFcn::EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    scalars[PostFcn::EDDY_VISCOSITY] = new char[sp + strlen(iod.output.transient.eddyvis)];
    sprintf(scalars[PostFcn::EDDY_VISCOSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.eddyvis);
  }
  if (iod.output.transient.dplus[0] != 0) {
#if defined(HEAT_FLUX)
    /*
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    double dT = iod.bc.wall.temperature - 1.0 / (gam*(gam-1.0)*iod.bc.inlet.mach*iod.bc.inlet.mach);
    sscale[PostFcn::DELTA_PLUS] = iod.ref.reynolds * iod.eqs.thermalCondModel.prandtl / (gam * dT);
    */
    sscale[PostFcn::DELTA_PLUS] = iod.ref.rv.tpower / (iod.ref.length*iod.ref.length);
#endif
    scalars[PostFcn::DELTA_PLUS] = new char[sp + strlen(iod.output.transient.dplus)];
    sprintf(scalars[PostFcn::DELTA_PLUS], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dplus);
  }
  if (iod.output.transient.philevel[0] != 0) {
    sscale[PostFcn::PHILEVEL] = 1.0;
    scalars[PostFcn::PHILEVEL] = new char[sp + strlen(iod.output.transient.philevel)];
    sprintf(scalars[PostFcn::PHILEVEL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel);
  }
  if (iod.output.transient.velocity[0] != 0) {
    vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocity);
  }
  if (iod.output.transient.tavvelocity[0] != 0) {
    avvscale[PostFcn::VELOCITYAVG] = iod.ref.rv.velocity;
    avvectors[PostFcn::VELOCITYAVG] = new char[sp + strlen(iod.output.transient.tavvelocity)];
    sprintf(avvectors[PostFcn::VELOCITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavvelocity);
  }
  if (iod.output.transient.displacement[0] != 0) {
    vscale[PostFcn::DISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::DISPLACEMENT] = new char[sp + strlen(iod.output.transient.displacement)];
    sprintf(vectors[PostFcn::DISPLACEMENT], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.displacement);
  }
  if (iod.output.transient.tavdisplacement[0] != 0) {
    avvscale[PostFcn::DISPLACEMENTAVG] = iod.ref.rv.tlength;
    avvectors[PostFcn::DISPLACEMENTAVG] = new char[sp + strlen(iod.output.transient.tavdisplacement)];
    sprintf(avvectors[PostFcn::DISPLACEMENTAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavdisplacement);
  }

  if (iod.output.transient.flightDisplacement[0] != 0) {
    vscale[PostFcn::FLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::FLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.flightDisplacement)];
    sprintf(vectors[PostFcn::FLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.flightDisplacement);
  }

  if (iod.output.transient.localFlightDisplacement[0] != 0) {
    vscale[PostFcn::LOCALFLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::LOCALFLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.localFlightDisplacement)];
    sprintf(vectors[PostFcn::LOCALFLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.localFlightDisplacement);
  }

  if (iod.output.transient.velocitynorm[0] != 0) {
    sscale[PostFcn::VELOCITY_NORM] = iod.ref.rv.velocity;
    scalars[PostFcn::VELOCITY_NORM] = new char[sp + strlen(iod.output.transient.velocitynorm)];
    sprintf(scalars[PostFcn::VELOCITY_NORM], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocitynorm);
  }

  int dsp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.dDensity[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_DENSITY] = iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_DENSITY] = new char[dsp + strlen(iod.output.transient.dDensity)];
    sprintf(dScalars[PostFcn::DERIVATIVE_DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dDensity);
  }

  if (iod.output.transient.dPressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_PRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_PRESSURE] = new char[dsp + strlen(iod.output.transient.dPressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dPressure);
  }

  if (iod.output.transient.dTemperature[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TEMPERATURE] = iod.ref.rv.temperature;
    dScalars[PostFcn::DERIVATIVE_TEMPERATURE] = new char[dsp + strlen(iod.output.transient.dTemperature)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTemperature);
  }

  if (iod.output.transient.dTotalpressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TOTPRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_TOTPRESSURE] = new char[dsp + strlen(iod.output.transient.dTotalpressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTotalpressure);
  }

  if (iod.output.transient.dNutturb[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_NUT_TURB] = new char[dsp + strlen(iod.output.transient.dNutturb)];
    sprintf(dScalars[PostFcn::DERIVATIVE_NUT_TURB], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dNutturb);
  }

  if (iod.output.transient.dEddyvis[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = new char[dsp + strlen(iod.output.transient.dEddyvis)];
    sprintf(dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dEddyvis);
  }

  if (iod.output.transient.dVelocityScalar[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = iod.ref.rv.velocity;
    dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = new char[dsp + strlen(iod.output.transient.dVelocityScalar)];
    sprintf(dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dVelocityScalar);
  }

  if (iod.output.transient.dVelocityVector[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = iod.ref.rv.velocity;
    dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = new char[dsp + strlen(iod.output.transient.dVelocityVector)];
    sprintf(dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dVelocityVector);
  }

  if (iod.output.transient.dDisplacement[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_DISPLACEMENT] = iod.ref.rv.tlength;
    dVectors[PostFcn::DERIVATIVE_DISPLACEMENT] = new char[dsp + strlen(iod.output.transient.dDisplacement)];
    sprintf(dVectors[PostFcn::DERIVATIVE_DISPLACEMENT], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dDisplacement);
  }

}

//------------------------------------------------------------------------------
