#include <cstdlib>
#include <cstring>

#include <TsOutput.h>

#include <IoData.h>
#include <RefVal.h>
#include <Domain.h>
#include <PostOperator.h>
#include <MeshMotionHandler.h>
#include <DistVector.h>
#include <DistExactRiemannSolver.h>
#include <BinFileHandler.h>
#include <VectorSet.h>
#include <GhostPoint.h>

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::TsOutput(IoData &iod, RefVal *rv, Domain *dom, PostOperator<dim> *po) : 
  refVal(rv), domain(dom), postOp(po), rmmh(0)
{
  int i;

  modeFile = 0;
  TavF = 0;
  TavM = 0;
  TavL = 0;
  Qs = 0;
  Qv = 0;
  output_newton_step = domain->getOutputNewtonStep();
  for (i=0; i<PostFcn::AVSSIZE; ++i) 
    {
      AvQs[i] = 0;
      AvQv[i] = 0;
    }

  steady = !iod.problem.type[ProblemData::UNSTEADY];
  com = domain->getCommunicator();

  int sp = strlen(iod.output.transient.prefix) + 1;
  int spn = strlen(iod.output.transient.probes.prefix) + 1;
  int sprom = strlen(iod.output.rom.prefix) + 1;

  if (iod.output.transient.solutions[0] != 0) {
    solutions = new char[sp + strlen(iod.output.transient.solutions)];
    sprintf(solutions, "%s%s", iod.output.transient.prefix, iod.output.transient.solutions);
  }
  else
    solutions = 0;

  // GAPPY POD STUFF (CBM+KTC)
  if (iod.output.rom.newtonresiduals[0] != 0) {
    newtonresiduals = new char[sprom + strlen(iod.output.rom.newtonresiduals)];
    sprintf(newtonresiduals, "%s%s", iod.output.rom.prefix, iod.output.rom.newtonresiduals);
  }
  else
    newtonresiduals = 0;

  if (iod.output.rom.jacobiandeltastate[0] != 0) {
    jacobiandeltastate = new char[sprom + strlen(iod.output.rom.jacobiandeltastate)];
    sprintf(jacobiandeltastate, "%s%s", iod.output.rom.prefix, iod.output.rom.jacobiandeltastate);
  }
  else
    jacobiandeltastate = 0;

  if (iod.output.rom.reducedjac[0] != 0) {
    reducedjac = new char[sprom + strlen(iod.output.rom.reducedjac)];
    sprintf(reducedjac, "%s%s", iod.output.rom.prefix, iod.output.rom.reducedjac);
  }
  else
    reducedjac = 0;

  for (i=0; i<PostFcn::SSIZE; ++i) {
    sscale[i] = 1.0;
    scalars[i] = 0;
    nodal_scalars[i] = 0;
  }
  for (i=0; i<PostFcn::AVSSIZE; ++i) {
    avsscale[i] = 1.0;
    avscalars[i] = 0;
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    vscale[i] = 1.0;
    vectors[i] = 0;
    nodal_vectors[i] = 0;
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
  if (iod.output.transient.pressurecoefficient[0] != 0) {
    sscale[PostFcn::PRESSURECOEFFICIENT] = 1.0;
    scalars[PostFcn::PRESSURECOEFFICIENT] = new char[sp + strlen(iod.output.transient.pressurecoefficient)];
    sprintf(scalars[PostFcn::PRESSURECOEFFICIENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressurecoefficient);
  }
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
//    sscale[PostFcn::TEMPERATURE] = 1;
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
  if (iod.output.transient.surfaceheatflux[0] != 0) {
    sscale[PostFcn::SURFACE_HEAT_FLUX] = iod.ref.rv.power /(iod.ref.rv.length * iod.ref.rv.length);
    scalars[PostFcn::SURFACE_HEAT_FLUX] = new char[sp + strlen(iod.output.transient.surfaceheatflux)];
    sprintf(scalars[PostFcn::SURFACE_HEAT_FLUX], "%s%s",
            iod.output.transient.prefix, iod.output.transient.surfaceheatflux);
  }
  if (iod.output.transient.tempnormalderivative[0] != 0) {
    sscale[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = iod.ref.rv.temperature/iod.ref.rv.length;
    scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = new char[sp + strlen(iod.output.transient.tempnormalderivative)];
    sprintf(scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tempnormalderivative);
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
  if (iod.output.transient.sfric[0] != 0) {
    scalars[PostFcn::SKIN_FRICTION] = new char[sp + strlen(iod.output.transient.sfric)];
    sprintf(scalars[PostFcn::SKIN_FRICTION], "%s%s",
            iod.output.transient.prefix, iod.output.transient.sfric);
  }
  if (iod.output.transient.tavsfric[0] != 0) {
    avscalars[PostFcn::SKIN_FRICTIONAVG] = new char[sp + strlen(iod.output.transient.tavsfric)];
    sprintf(avscalars[PostFcn::SKIN_FRICTIONAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavsfric);
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
  if (iod.output.transient.tavcsdles[0] != 0) {
    avscalars[PostFcn::CSDLESAVG] = new char[sp + strlen(iod.output.transient.tavcsdles)];
    sprintf(avscalars[PostFcn::CSDLESAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavcsdles);
  }
  if (iod.output.transient.csdvms[0] != 0) {
    scalars[PostFcn::CSDVMS] = new char[sp + strlen(iod.output.transient.csdvms)];
    sprintf(scalars[PostFcn::CSDVMS], "%s%s",
            iod.output.transient.prefix, iod.output.transient.csdvms);
  }
  if (iod.output.transient.tavcsdvms[0] != 0) {
    avscalars[PostFcn::CSDVMSAVG] = new char[sp + strlen(iod.output.transient.tavcsdvms)];
    sprintf(avscalars[PostFcn::CSDVMSAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavcsdvms);
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
  if (iod.output.transient.fluidid[0] != 0) {
    sscale[PostFcn::FLUIDID] = 1.0;
    scalars[PostFcn::FLUIDID] = new char[sp + strlen(iod.output.transient.fluidid)];
    sprintf(scalars[PostFcn::FLUIDID], "%s%s",
            iod.output.transient.prefix, iod.output.transient.fluidid);
  }
  if (iod.output.transient.controlvolume[0] != 0) {
    sscale[PostFcn::CONTROL_VOLUME] = iod.ref.rv.length * iod.ref.rv.length * iod.ref.rv.length;
    scalars[PostFcn::CONTROL_VOLUME] = new char[sp + strlen(iod.output.transient.controlvolume)];
    sprintf(scalars[PostFcn::CONTROL_VOLUME], "%s%s",
            iod.output.transient.prefix, iod.output.transient.controlvolume);
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
    
  mX = 0;
  if (iod.output.transient.generalizedforces[0] != 0 && strcmp(iod.input.strModesFile, "") != 0) {
    generalizedforces = new char[sp + strlen(iod.output.transient.generalizedforces)];
    sprintf(generalizedforces, "%s%s", iod.output.transient.prefix, iod.output.transient.generalizedforces);
    modeFile = new char [MAXLINE];
    sprintf(modeFile, "%s%s", iod.input.prefix, iod.input.strModesFile);
    DistSVec<double, dim> tmpVec(domain->getNodeDistInfo());
    DistSVec<double, 3> Xtmp(domain->getNodeDistInfo());
    double f;

    if (strcmp(modeFile, "") != 0)  {
      com->fprintf(stderr, " ... Reading Modefile %s\n", modeFile);
      domain->readVectorFromFile(modeFile, 0, &f, Xtmp);
      int nStrMode = int(f);

      // We read the modal deformations
      mX = new VecSet<DistSVec<double,3> >(nStrMode, domain->getNodeDistInfo());

      for(int iMode = 0; iMode < nStrMode; ++iMode) 
        domain->readVectorFromFile(modeFile, iMode+1, &f, (*mX)[iMode]);
     }
     else 
      mX = 0;
  }
  else 
    {
      generalizedforces = 0;
      modeFile =0;
    }
  
  if (iod.output.transient.generalizedforces[0] != 0 && !(iod.problem.type[ProblemData::FORCED]) && strcmp(iod.input.strModesFile, "") == 0)
    fprintf(stderr, "Error : StrModes file for Generalized Forces not given.. Aborting !! \n");
    
  
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

  if (iod.output.transient.heatfluxes[0] != 0) {
    heatfluxes = new char[sp + strlen(iod.output.transient.heatfluxes)];
    sprintf(heatfluxes, "%s%s", iod.output.transient.prefix, iod.output.transient.heatfluxes);
  }
  else
    heatfluxes = 0;

  if (iod.output.transient.residuals[0] != 0) {
    residuals = new char[sp + strlen(iod.output.transient.residuals)];
    sprintf(residuals, "%s%s", iod.output.transient.prefix, iod.output.transient.residuals);
  }
  else
    residuals = 0;

      if (iod.output.rom.staterom[0] != 0) {
    staterom = new char[sprom + strlen(iod.output.rom.staterom)];
    sprintf(staterom, "%s%s", iod.output.rom.prefix, iod.output.rom.staterom);
  }
  else
    staterom = 0;

  if (iod.output.rom.error[0] != 0) {
    error = new char[sprom + strlen(iod.output.rom.error)];
    sprintf(error, "%s%s", iod.output.rom.prefix, iod.output.rom.error);
  }
  else
    error = 0;
    
  if (iod.output.transient.materialVolumes[0] != 0) {
    material_volumes = new char[sp + strlen(iod.output.transient.materialVolumes)];
    sprintf(material_volumes, "%s%s", iod.output.transient.prefix, iod.output.transient.materialVolumes);
  }
  else
    material_volumes = 0;

  if (iod.output.transient.embeddedsurface[0] != 0) {
    embeddedsurface = new char[sp + strlen(iod.output.transient.embeddedsurface)];
    sprintf(embeddedsurface, "%s%s", iod.output.transient.prefix, iod.output.transient.embeddedsurface); 
  }
  else
    embeddedsurface = 0;

  if (iod.output.transient.cputiming[0] != 0) {
    cputiming = new char[sp + strlen(iod.output.transient.cputiming)];
    sprintf(cputiming, "%s%s", iod.output.transient.prefix, iod.output.transient.cputiming); 
  }
  else
    cputiming = 0;

  if (iod.output.transient.conservation[0] != 0) {
    conservation = new char[sp + strlen(iod.output.transient.conservation)];
    sprintf(conservation, "%s%s", iod.output.transient.prefix, iod.output.transient.conservation);
  }
  else
    conservation = 0;

  it0 = iod.restart.iteration;
  //std::cout << "it0 = " << it0 << std::endl;
  numFluidPhases = iod.eqs.numPhase;
  frequency = iod.output.transient.frequency;
  frequency_dt = iod.output.transient.frequency_dt;
  prtout = iod.restart.etime;
  length = iod.output.transient.length;
  surface = iod.output.transient.surface;
  x0[0] = iod.output.transient.x0;
  x0[1] = iod.output.transient.y0;
  x0[2] = iod.output.transient.z0;

  fpResiduals = 0;
  fpStateRom = 0;
  fpMatVolumes = 0;
  fpConservationErr = 0;
  fpGnForces  = 0;
  fpError = 0;

  int nSurf = postOp->getNumSurf();
  int nSurfHF = postOp->getNumSurfHF();
  fpForces             = new FILE *[nSurf];
  fpLift               = new FILE *[nSurf];
  fpTavForces          = new FILE *[nSurf];
  fpTavLift            = new FILE *[nSurf];
  fpHydroStaticForces  = new FILE *[nSurf];
  fpHydroDynamicForces = new FILE *[nSurf];
  fpHydroStaticLift    = new FILE *[nSurf];
  fpHydroDynamicLift   = new FILE *[nSurf];
  fpHeatFluxes         = new FILE *[nSurfHF]; 

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

  for (int iSurf = 0; iSurf < nSurfHF; iSurf++)  {
    fpHeatFluxes[iSurf]         = 0; 
  }


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

  // Initialize nodal output structures
  Probes& myProbes = iod.output.transient.probes;
  nodal_output.step = 0;
  nodal_output.results = new double[Probes::MAXNODES*3];
  nodal_output.subId = new int[Probes::MAXNODES];
  nodal_output.locNodeId = new int[Probes::MAXNODES];
  nodal_output.last = new int[Probes::MAXNODES];
  nodal_output.locations.resize(Probes::MAXNODES);

  nodal_output.step = 0;
					      
  for (i = 0; i < Probes::MAXNODES; ++i) {
    nodal_output.locations[i] = Vec3D(myProbes.myNodes[i].locationX/iod.ref.rv.length,
                                      myProbes.myNodes[i].locationY/iod.ref.rv.length,
                                      myProbes.myNodes[i].locationZ/iod.ref.rv.length);
     
    nodal_output.last[i] = 0;
    if (myProbes.myNodes[i].id >= 0) {
      com->fprintf(stdout,"[Probe] Node %d: NodeId = %d.\n", i+1, myProbes.myNodes[i].id);

      int flag = -1;
      int locid = -1;
      int lis = -1;
#pragma omp parallel for
      for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) {
	locid = dom->getSubDomain()[iSub]->getLocalNodeNum( myProbes.myNodes[i].id );
//	fprintf(stdout,"locid = %i\n",locid);
	if (locid >= 0) {
	  lis = iSub;
	  flag = com->cpuNum();
	  break;
	}
      }

      com->globalMax(1,&flag);

      if ( flag == com->cpuNum() ) {
	nodal_output.locNodeId[i] = locid;
	nodal_output.subId[i] = lis;
      } else {
	nodal_output.locNodeId[i] = -1;
	nodal_output.subId[i] = -1;
      }
    } else if (myProbes.myNodes[i].locationX < -1.0e10)
      break;
    else
      com->fprintf(stdout,"[Probe] Node %d: Coords = (%e, %e, %e).\n", i+1, 
                   myProbes.myNodes[i].locationX, myProbes.myNodes[i].locationY, myProbes.myNodes[i].locationZ);
  }

  com->fprintf(stdout,"[Probe] Number of probing nodes is %d\n",i);

  nodal_output.numNodes = i;

  if (iod.output.transient.probes.density[0] != 0) {
    nodal_scalars[PostFcn::DENSITY] = new char[spn + strlen(iod.output.transient.probes.density)];
    sprintf(nodal_scalars[PostFcn::DENSITY], "%s%s", 
	    iod.output.transient.probes.prefix, iod.output.transient.probes.density);
  }
  if (iod.output.transient.probes.pressure[0] != 0) {
    nodal_scalars[PostFcn::PRESSURE] = new char[spn + strlen(iod.output.transient.probes.pressure)];
    sprintf(nodal_scalars[PostFcn::PRESSURE], "%s%s", 
	    iod.output.transient.probes.prefix, iod.output.transient.probes.pressure);
  }
  if (iod.output.transient.probes.temperature[0] != 0) {
    nodal_scalars[PostFcn::TEMPERATURE] = new char[spn + strlen(iod.output.transient.probes.temperature)];
    sprintf(nodal_scalars[PostFcn::TEMPERATURE], "%s%s", 
	    iod.output.transient.probes.prefix, iod.output.transient.probes.temperature);
  }
  if (iod.output.transient.probes.velocity[0] != 0) {
    nodal_vectors[PostFcn::VELOCITY] = new char[spn + strlen(iod.output.transient.probes.velocity)];
    sprintf(nodal_vectors[PostFcn::VELOCITY], "%s%s", 
	    iod.output.transient.probes.prefix, iod.output.transient.probes.velocity);
  }
  if (iod.output.transient.probes.displacement[0] != 0) {
    nodal_vectors[PostFcn::DISPLACEMENT] = new char[spn + strlen(iod.output.transient.probes.displacement)];
    sprintf(nodal_vectors[PostFcn::DISPLACEMENT], "%s%s", 
	    iod.output.transient.probes.prefix, iod.output.transient.probes.displacement);
  }

  tscale = iod.ref.rv.time;
  xscale = iod.ref.rv.length;
}

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::~TsOutput()
{

  for (int i=0; i<PostFcn::AVSSIZE; ++i) {
    delete AvQs[i];
    delete AvQv[i];
  }
  if (Qs) delete Qs;
  if (Qv) delete Qv;
  if(TavF) delete [] TavF;
  if(TavM) delete [] TavM;
  if(TavL) delete [] TavL;
  if(modeFile) delete [] modeFile;
  if(mX) delete mX;

  if (switchOpt) //STEADY_SENSITIVITY_ANALYSIS
    {
      delete[] dForces;
      
      int i;
      for (i=0; i<PostFcn::DSSIZE; ++i) {
	delete[]  dScalars[i];
      }
      for (i=0; i<PostFcn::DVSIZE; ++i) {
	delete[] dVectors[i];
      }
      delete[] dSolutions;
    }

  delete[]  fpForces            ;
  delete[]  fpLift              ;
  delete[]  fpTavForces         ;
  delete[]  fpTavLift           ;
  delete[]  fpHydroStaticForces ;
  delete[]  fpHydroDynamicForces;
  delete[]  fpHydroStaticLift   ;
  delete[]  fpHydroDynamicLift  ;
  delete[]  fpHeatFluxes        ; 

  delete[] heatfluxes;
  delete[] residuals;
  delete[] material_volumes;
  delete[] embeddedsurface;
  delete[] cputiming;
  delete[] staterom;
  delete[] error;
  delete[] conservation;

  delete[] lift;
  delete[] tavlift;
  delete[] hydrostaticlift;
  delete[] hydrodynamiclift;

  delete mX;

  delete[] generalizedforces;
  delete[] modeFile;

  delete[] forces;
  delete[] tavforces;
  delete[] hydrostaticforces;
  delete[] hydrodynamicforces;

  int i;
  for (i=0; i<PostFcn::SSIZE; ++i) {
    delete[] scalars[i];
  }
  for (i=0; i<PostFcn::AVSSIZE; ++i) {
     delete[]  avscalars[i];
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
     delete[]  vectors[i];
  }
  for (i=0; i<PostFcn::AVVSIZE; ++i) {
     delete[]  avvectors[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
bool TsOutput<dim>::toWrite(int it, bool lastIt, double t)
{
  if(frequency_dt<=0.0)
    return (((frequency > 0) && (it % frequency == 0)) || lastIt);

  return (t>=prtout || lastIt);
}

//------------------------------------------------------------------------------

template<int dim>
int TsOutput<dim>::getStep(int it, bool lastIt, double t)
{
  int step = 0;
  if(frequency_dt<=0.0) {
    if (frequency > 0) {
      step = it / frequency;
      if (lastIt && (it % frequency != 0))
        ++step;
    }
  } else {
    step = (int)(t / frequency_dt);
    if (lastIt && (t<prtout))
      ++step;
  }

  return step;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::updatePrtout(double t)
{
  if(frequency_dt<=0.0)
    return;
  if(t>=prtout)
    prtout += frequency_dt;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::setMeshMotionHandler(IoData &ioData, MeshMotionHandler *mmh)
{

  rmmh = dynamic_cast<RigidMeshMotionHandler *>(mmh);

  if (ioData.problem.type[ProblemData::FORCED]) {
     if (ioData.forced.type == ForcedData::HEAVING)
       hmmh = dynamic_cast<HeavingMeshMotionHandler *>(mmh);
     else if (ioData.forced.type  == ForcedData::PITCHING)
       pmmh = dynamic_cast<PitchingMeshMotionHandler *>(mmh);
     else if (ioData.forced.type  == ForcedData::DEFORMING)
       dmmh = dynamic_cast<DeformingMeshMotionHandler *>(mmh);
  }

  if (ioData.output.transient.generalizedforces[0] != 0 && ioData.problem.type[ProblemData::FORCED] && strcmp(ioData.input.strModesFile, "") == 0) {
    int sp = strlen(ioData.output.transient.prefix) + 1;
    generalizedforces = new char[sp + strlen(ioData.output.transient.generalizedforces)];
    sprintf(generalizedforces, "%s%s", ioData.output.transient.prefix, ioData.output.transient.generalizedforces);

    mX = new VecSet<DistSVec<double,3> >(1, domain->getNodeDistInfo());

    if (hmmh)
       (*mX)[0] = hmmh->getModes(); 
    else if(pmmh)
       (*mX)[0] = pmmh->getModes(); 
    else if(dmmh)
       (*mX)[0] = dmmh->getModes();

  }
//  else
//    generalizedforces = 0;

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
  char* toto = fgets(line, MAXLINE, fpback);
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


  if (generalizedforces) {
    if (it0 != 0)
      fpGnForces = backupAsciiFile(generalizedforces);
    if (it0 == 0 || fpGnForces == 0) {
      fpGnForces = fopen(generalizedforces, "w");
      if (!fpGnForces) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", generalizedforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpGnForces, "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpGnForces, "Coefficient of Generalized Forces %s\n", addvar);
      else
        fprintf(fpGnForces, "Generalized Forces %s\n", addvar);
    }
    fflush(fpGnForces);
  }

  if (generalizedforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", generalizedforces, surfNums[iSurf]);
      if (it0 != 0)
        fpGnForces = backupAsciiFile(filename);
      if (it0 == 0 || fpGnForces == 0) {
        fpGnForces = fopen(filename, "w");
        if (!fpGnForces) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpGnForces, "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpGnForces, "Coefficient of Generalized Forces %s\n", addvar);
        else
          fprintf(fpGnForces, "Generalized Forces %s\n", addvar);
      }
      fflush(fpGnForces);
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

//For Heat Fluxes the flags are different

  int nSurfHF = postOp->getNumSurfHF();
   int *surfNumsHF = 0;
    if (nSurfHF > 0)  {
    surfNumsHF = new int[nSurfHF];
    map<int,int> surfMapHF = postOp->getSurfMapHF();
    map<int, int>::iterator it;
    iSurf = 1;
    surfNumsHF[0] = 0;
    for (it = surfMapHF.begin(); it != surfMapHF.end(); it++)
      if (it->second > 0)
        surfNumsHF[iSurf++] = it->first;
  }

  if (heatfluxes) {
    if (it0 != 0)
      fpHeatFluxes[0] = backupAsciiFile(heatfluxes);
    if (it0 == 0 || fpHeatFluxes[0] == 0) {
      fpHeatFluxes[0] = fopen(heatfluxes, "w");
      if (!fpHeatFluxes[0]) {
        fprintf(stderr, "*** HF Error: could not open \'%s\'\n", heatfluxes);
        exit(1);
      }
      const char *addvar = "";
 //     if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHeatFluxes[0], "# TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHeatFluxes[0], "Nondimensional HeatFlux %s\n", addvar);
      else
        fprintf(fpHeatFluxes[0], "HeatFlux %s\n", addvar);
    }
    fflush(fpHeatFluxes[0]);
  }

  if (heatfluxes) {
    for (iSurf = 1; iSurf < nSurfHF; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", heatfluxes, surfNumsHF[iSurf]);
      if (it0 != 0)
        fpHeatFluxes[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHeatFluxes[iSurf] == 0) {
        fpHeatFluxes[iSurf] = fopen(filename, "w");
        if (!fpHeatFluxes[iSurf]) {
           fprintf(stderr, "*** HF Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHeatFluxes[iSurf], "# TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHeatFluxes[iSurf], "Nondimensional HeatFlux  %s\n", addvar);
        else
          fprintf(fpHeatFluxes[iSurf], "HeatFlux %s\n", addvar);
      }
      fflush(fpHeatFluxes[iSurf]);
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
  
  if (staterom) {
    if (it0 != 0) 
      fpStateRom = backupAsciiFile(staterom);
    if (it0 == 0 || fpStateRom == 0) {
      fpStateRom = fopen(staterom, "w");
      if (!fpStateRom) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", staterom);
	exit(1);
      }
      fprintf(fpStateRom, "# TimeIteration ElapsedTime StatePodCoords \n");
    }
    fflush(fpStateRom);
  }

  if (error) {
    if (it0 != 0) 
      fpError = backupAsciiFile(error);
    if (it0 == 0 || fpError == 0) {
      fpError = fopen(error, "w");
      if (!fpError) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", error);
	exit(1);
      }
      fprintf(fpError, "# TimeIteration ElapsedTime RelativeError AbsoluteError \n");
    }
    fflush(fpError);
  }

  if (material_volumes) {
    if (it0 != 0)
      fpMatVolumes = backupAsciiFile(material_volumes);
    if (it0 == 0 || fpMatVolumes == 0) {
      fpMatVolumes = fopen(material_volumes, "w");
      if (!fpMatVolumes) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", material_volumes);
        exit(1);
      }
      fprintf(fpMatVolumes, "# TimeIteration ElapsedTime ");
      for(int i=0; i<numFluidPhases; i++)
        fprintf(fpMatVolumes, "Volume[FluidID==%d] ", i);
      fprintf(fpMatVolumes, "Volume[FluidID==%d(GhostSolid)] TotalVolume\n", numFluidPhases);
    }
    fflush(fpMatVolumes);
  }

  if (embeddedsurface) {
    if (it0 != 0)
      fpEmbeddedSurface = backupAsciiFile(embeddedsurface);
    if (it0 == 0 || fpEmbeddedSurface == 0) {
      fpEmbeddedSurface = fopen(embeddedsurface, "w");
      if (!fpEmbeddedSurface) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", embeddedsurface);
        exit(1);
      }
      fprintf(fpEmbeddedSurface, "Vector DISP under NLDynamic for EmbeddedNodes\n");
    }
    fflush(fpEmbeddedSurface);
  }

  if (cputiming) {
    if (it0 != 0)
      fpCpuTiming = backupAsciiFile(cputiming);
    if (it0 == 0 || fpCpuTiming== 0) {
      fpCpuTiming= fopen(cputiming, "w");
      if (!fpCpuTiming) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", cputiming);
        exit(1);
      }
    }
    fflush(fpCpuTiming);
  }

  if (conservation) {
    if (it0 != 0) 
      fpConservationErr = backupAsciiFile(conservation);
    if (it0 == 0 || fpConservationErr == 0) {
      fpConservationErr = fopen(conservation, "w");
      if (!fpConservationErr) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", conservation);
	exit(1);
      }
      fprintf(fpConservationErr, "# 1TimeIteration 2Time ");
      fprintf(fpConservationErr, "3TotalExpectedMass 4TotalExpectedMomentumx 5TotalExpectedMomentumy 6TotalExpectedMomentumz 7TotalExpectedEnergy ");
      fprintf(fpConservationErr, "8Fluid1ExpectedMass 9Fluid1ExpectedMomentumx 10Fluid1ExpectedMomentumx 11Fluid1ExpectedMomentumz 12Fluid1ExpectedEnergy ");
      fprintf(fpConservationErr, "13Fluid2ExpectedMass 14Fluid2ExpectedMomentumx 15Fluid2ExpectedMomentumx 16Fluid2ExpectedMomentumz 17Fluid2ExpectedEnergy ");
      fprintf(fpConservationErr, "18TotalComputedMass 19TotalComputedMomentumx 20TotalComputedMomentumy 21TotalComputedMomentumz 22TotalComputedEnergy ");
      fprintf(fpConservationErr, "23Fluid1ComputedMass 24Fluid1ComputedMomentumx 25Fluid1ComputedMomentumx 26Fluid1ComputedMomentumz 27Fluid1ComputedEnergy ");
      fprintf(fpConservationErr, "28Fluid2ComputedMass 29Fluid2ComputedMomentumx 30Fluid2ComputedMomentumx 31Fluid2ComputedMomentumz 32Fluid2ComputedEnergy\n");
    }
    fflush(fpConservationErr);
 }

 delete [] surfNums;
 delete [] surfNumsHF;
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
  for (int iSurf = 0; iSurf < postOp->getNumSurfHF(); iSurf++)  {
     if (fpHeatFluxes[iSurf]) fclose(fpHeatFluxes[iSurf]);
  }
  if (fpResiduals) fclose(fpResiduals);
  if (fpMatVolumes) fclose(fpMatVolumes);
  if (fpEmbeddedSurface) fclose(fpEmbeddedSurface);
  if (fpCpuTiming) fclose(fpCpuTiming);
  if (fpStateRom) fclose(fpStateRom);
  if (fpError) fclose(fpError);
  if (fpGnForces) fclose(fpGnForces);
  if (fpConservationErr) fclose(fpConservationErr);
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeForcesToDisk(DistExactRiemannSolver<dim> &riemann,
                                      bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec<double> GF( mX ? mX->numVectors(): 0);
  if(mX) GF = 0.0;

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
    // if single-phase flow -- fluidId is a null pointer
    // if multi-phase flow  -- fluidId points to a DistVec<int>
    postOp->computeForceAndMoment(riemann, rVec, X, U, fluidId, Fi, Mi, Fv, Mv, 0, mX, mX ? &GF: 0);    
    if(mX != 0) 
      com->globalSum(GF.size(), GF.data());

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
  if (fpGnForces) {

     fprintf(fpGnForces,"%d %e %d %d", it, time, itSc, itNl);
     for (int i=0; i < mX->numVectors(); i++) {
       if (refVal->mode == RefVal::NON_DIMENSIONAL) 
         fprintf(fpGnForces," %e", 2.0*refVal->length*refVal->length*GF[i]/surface);
       else 
         fprintf(fpGnForces," %e", refVal->force*GF[i]);
     } 

      fprintf(fpGnForces, "\n");

      fflush(fpGnForces);
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

  delete[] Fi;
  delete[] Fv;
  delete[] Mi;
  delete[] Mv;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeForcesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec<double> GF( mX ? mX->numVectors(): 0);
  if(mX) GF = 0.0;

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
    // if single-phase flow -- fluidId is a null pointer
    // if multi-phase flow  -- fluidId points to a DistVec<int>
    postOp->computeForceAndMoment(rVec, X, U, fluidId, Fi, Mi, Fv, Mv, 0, mX, mX ? &GF: 0);    
    if(mX != 0) 
      com->globalSum(GF.size(), GF.data());

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
  if (fpGnForces) {

     fprintf(fpGnForces,"%d %e %d %d", it, time, itSc, itNl);
     for (int i=0; i < mX->numVectors(); i++) {
       if (refVal->mode == RefVal::NON_DIMENSIONAL) 
         fprintf(fpGnForces," %e", 2.0*refVal->length*refVal->length*GF[i]/surface);
       else 
         fprintf(fpGnForces," %e", refVal->force*GF[i]);
     } 

      fprintf(fpGnForces, "\n");

      fflush(fpGnForces);
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

  delete[] Fi;
  delete[] Fv;
  delete[] Mi;
  delete[] Mv;
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
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
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

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId points to a DistVec<int>
  if (hydrostaticforces)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamicforces)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiD, MiD, FvD, MvD, 2);

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

  delete []FiS;
  delete []MiS;
  delete []FvS;
  delete []MvS;
  delete []FiD;
  delete []MiD;
  delete []FvD;
  delete []MvD;

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu, 
				      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
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

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId points to a DistVec<int>
  if (lift || tavlift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, Fi, Mi, Fv, Mv);    

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

  delete [] Fi;
  delete [] Fv;
  delete [] Mi;
  delete [] Mv;

}
//------------------------------------------------------------------------------
                                                                                                                                                                                                     
template<int dim>
void TsOutput<dim>::writeHydroLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
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

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId  points to a DistVec<int>
  if (hydrostaticlift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamiclift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiD, MiD, FvD, MvD, 2);

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

  delete [] FiS;
  delete [] FvS;
  delete [] MiS;
  delete [] MvS;
  delete [] FiD;
  delete [] FvD;
  delete [] MiD;
  delete [] MvD;
}

//------------------------------------------------------------------------------
template<int dim>
  void TsOutput<dim>::writeHeatFluxesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{
  int nSurfs = postOp->getNumSurfHF();

  double *HF = new double[nSurfs];
  for(int index =0; index < nSurfs; index++){
    HF[index] = 0.0;
  }

  double time = refVal->time * t;

  if (heatfluxes)
    postOp->computeHeatFluxes(X,U,HF);

  int iSurf;
  if (fpHeatFluxes[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        HF[iSurf] *= 2.0 * refVal->length*refVal->length / surface; 
      else
       HF[iSurf] *= refVal->power; 
        fprintf(fpHeatFluxes[iSurf], "%d %e %d %d %e \n",
                it, time, itSc, itNl, HF[iSurf]);
      fflush(fpHeatFluxes[iSurf]);
    }
  }
  delete[] HF;
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
void TsOutput<dim>::writeMaterialVolumesToDisk(int it, double t, DistVec<double> &A, DistVec<int> *fluidId)
{
  if(!material_volumes)
    return;

  int myLength = numFluidPhases + 1/*ghost*/;
  double Vol[myLength];
  for(int i=0; i<myLength; i++)
    Vol[i] = 0.0;

  domain->computeMaterialVolumes(Vol,myLength,A,fluidId); //computes Vol

  if (com->cpuNum() !=0 ) return;

  double length3 = length*length*length;
  for(int i=0; i<myLength; i++)
    Vol[i] *= length3; //dimensionalize

  fprintf(fpMatVolumes, "%d %e ", it, (refVal->time)*t);
  for(int i=0; i<numFluidPhases+1; i++)
    fprintf(fpMatVolumes, "%e ", Vol[i]);
  
  double totVol = 0.0;
  for(int i=0; i<myLength; i++)
    totVol += Vol[i];

  fprintf(fpMatVolumes, "%e\n", totVol);

  fflush(fpMatVolumes);
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeEmbeddedSurfaceToDisk(bool lastIt, int it, double t, Vec<Vec3D>& solidX, Vec<Vec3D>& solidX0)
{
  if(!embeddedsurface)
    return;

  if(toWrite(it,lastIt,t)) {
    if(com->cpuNum() != 0) return;
    fprintf(fpEmbeddedSurface, "%e\n", t*tscale);
    for(int i=0; i<solidX.size(); i++) 
      fprintf(fpEmbeddedSurface, "%e %e %e\n", xscale*(solidX[i][0]-solidX0[i][0]), xscale*(solidX[i][1]-solidX0[i][1]), 
              xscale*(solidX[i][2]-solidX0[i][2]));
    fflush(fpEmbeddedSurface); 
    fprintf(stdout, "Wrote solution %d to \'%s\'\n", getStep(it, lastIt, t), embeddedsurface);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeCPUTimingToDisk(bool lastIt, int it, double t, Timer *timer)
{
  if(!cputiming) 
    return;

  if(toWrite(it,lastIt,t)) {
    com->fprintf(fpCpuTiming, "It %d, Time = %e.", it, t*tscale);
    timer->setRunTime();
    timer->print(NULL, fpCpuTiming);
    if(com->cpuNum()==0)
      fflush(fpCpuTiming);
    com->fprintf(stdout, "Wrote CPUTiming %d to \'%s\'\n", getStep(it, lastIt, t), cputiming);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeConservationErrors(IoData &iod, int it, double t,
                    int numPhases, double **expected, double **computed)
{

  if (com->cpuNum() != 0) return;

  if (!steady)
    if (fpConservationErr) {
      fprintf(fpConservationErr, "%d %e ", it, t);
      for(int iPhase=0; iPhase<numPhases; iPhase++)
        for(int i=0; i<dim; i++) 
          fprintf(fpConservationErr, "%e ", expected[iPhase][i]);
      for(int iPhase=0; iPhase<numPhases; iPhase++)
        for(int i=0; i<dim; i++) 
          fprintf(fpConservationErr, "%e ", computed[iPhase][i]);
      fprintf(fpConservationErr, "\n");
      fflush(fpConservationErr);
      
    }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDiskRom(bool lastIt, int it, double t,
		DistSVec<double,dim> *U1 = NULL, DistSVec<double,dim> *U2 = NULL, VecSet< DistSVec<double,dim> > *U3 = NULL)
{

	double tag = 0.0;	// do not tag for now

	if (newtonresiduals && U1)  {	// for both pg residuals and fom residuals
		DistSVec<double,dim> soltn(*U1);
		if (refVal->mode == RefVal::DIMENSIONAL)
			domain->scaleSolution(soltn, refVal);
		domain->writeVectorToFile(newtonresiduals, *output_newton_step, tag, *U1);	//TODO: output_newton_step should accumulate over restarts
	}

	if (jacobiandeltastate && U2)  {
		DistSVec<double,dim> soltn(*U2);
		if (refVal->mode == RefVal::DIMENSIONAL)
			domain->scaleSolution(soltn, refVal);
		domain->writeVectorToFile(jacobiandeltastate, *output_newton_step, tag, *U2);	//TODO: output_newton_step should accumulate over restarts
	}

	if (reducedjac && U3)  {	// entire JPhi
		for (int i = 0; i < U3->numVectors(); ++i) {
			DistSVec<double,dim> soltn((*U3)[i]);
			if (refVal->mode == RefVal::DIMENSIONAL)
				domain->scaleSolution(soltn, refVal);

			domain->writeVectorToFile(reducedjac, *output_newton_step, tag, soltn);	//TODO: output_newton_step should accumulate over restarts
		}
	}

	++(*output_newton_step);
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeDisplacementVectorToDisk(int step, double tag, 
                      DistSVec<double,3> &X, DistSVec<double,dim> &U){

  if(vectors[PostFcn:: DISPLACEMENT]){
    if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
    postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(PostFcn::DISPLACEMENT), X, U, *Qv);
    domain->writeVectorToFile(vectors[PostFcn::DISPLACEMENT], 1, 1.0, *Qv, &(vscale[PostFcn::DISPLACEMENT]));
  }


}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X, 
					     DistVec<double> &A, DistSVec<double,dim> &U, DistTimeState<dim> *timeState)
{

  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
		double tag;
		if (rmmh)
			tag = rmmh->getTagValue(t);
		else
			tag = t * refVal->time;

		if (solutions)  {
			DistSVec<double,dim> soltn(U);
			if (refVal->mode == RefVal::DIMENSIONAL)
				domain->scaleSolution(soltn, refVal);
			domain->writeVectorToFile(solutions, step, tag, soltn);
		}

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

template<int dim>
template<int dimLS>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U, 
                                             DistTimeState<dim> *timeState,
                                             DistVec<int> &fluidId,DistSVec<double,dimLS>* Phi)
{
  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
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
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState,fluidId,Phi);
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
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, fluidId);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}

static void copyFile(const char* fname) {

  FILE* f = fopen(fname,"rb");
  fseek (f , 0 , SEEK_END);
  size_t lSize = ftell (f);
  rewind (f);
  char* buffer = new char[lSize];
  size_t err = fread (buffer,1,lSize,f);
  fclose(f);

  char nn[256];
  sprintf(nn,"%s.back",fname);
  f = fopen(nn,"wb");
  fwrite(buffer,1,lSize,f);
  fclose(f);

  delete [] buffer;
}

template<int dim>
void TsOutput<dim>::cleanProbesFile() {

  char nn[256];
  int iter,i,n;
  double time,res;
  if (it0 == 0) 
    return;

  for (i=0; i<PostFcn::SSIZE; ++i) {
    if (nodal_scalars[i]) {
      if (com->cpuNum() == 0) {
        copyFile(nodal_scalars[i]);
        sprintf(nn,"%s.back",nodal_scalars[i]);
        FILE* scalar_file = fopen(nodal_scalars[i],"w");
        FILE* scalar_file_old = fopen(nn,"r");
        while (!feof(scalar_file_old)) {
          n = fscanf(scalar_file_old,"%d",&iter);
          n = fscanf(scalar_file_old,"%lf",&time);
          if (iter > it0)
            break;          
          fprintf(scalar_file,"%d %e ",iter,time);
	  for (int k =0 ; k < nodal_output.numNodes; ++k) {
	    n = fscanf(scalar_file_old,"%lf",&res);
            fprintf(scalar_file,"%e ",res);
          }
          fprintf(scalar_file,"\n");
	  
	}
        fclose(scalar_file);
        fclose(scalar_file_old);
      }
    }
  }
  
  for (i=0; i<PostFcn::VSIZE; ++i) {
    if (nodal_vectors[i]) {
	
      if (com->cpuNum() == 0) {
        copyFile(nodal_vectors[i]);
        sprintf(nn,"%s.back",nodal_vectors[i]);
        FILE* vector_file = fopen(nodal_vectors[i],"w");
        FILE* vector_file_old = fopen(nn,"r");
        while (!feof(vector_file_old)) {
          n = fscanf(vector_file_old,"%d",&iter);
          n = fscanf(vector_file_old,"%lf",&time);
          if (iter > it0)
            break;          
          fprintf(vector_file,"%d %e ",iter,time);
          for (int k =0 ; k < nodal_output.numNodes; ++k) {
            for (int l = 0; l < 3; ++l) {
	      n = fscanf(vector_file_old,"%lf",&res);
              fprintf(vector_file,"%e ",res);
            }
          }
          fprintf(vector_file,"\n");
	}
        fclose(vector_file);
        fclose(vector_file_old);
      }
    }
  }
}

template<int dim>
template<int dimLS>
void TsOutput<dim>::writeProbesToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
				      DistVec<double> &A, DistSVec<double,dim> &U, 
				      DistTimeState<dim> *timeState, DistVec<int> &fluidId,
                                      DistSVec<double,dimLS>* Phi, DistLevelSetStructure *distLSS,
                                      DistVec<GhostPoint<dim>*> *ghostPoints)
{
  //if (toWrite(it,lastIt,t)) {
  if (nodal_output.numNodes == 0)
    return;
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;
    
    // if (solutions)
    // domain->writeVectorToFile(solutions, step, tag, U);
    
    int i;
    const char* mode = nodal_output.step ? "a" : "w";
    if (it0 > 0)
      mode = "a";

    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (nodal_scalars[i]) {
	
	postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A,
				      timeState,fluidId,
				      nodal_output.subId, nodal_output.locNodeId,
				      nodal_output.last,nodal_output.numNodes,nodal_output.results,
                                      nodal_output.locations,
                                      Phi, distLSS, ghostPoints);
	if (com->cpuNum() == 0) {
	  FILE* scalar_file = fopen(nodal_scalars[i],mode);
          if (scalar_file != 0) {
	    fprintf(scalar_file,"%d %e ",nodal_output.step+it0, tag);
	    for (int k =0 ; k < nodal_output.numNodes; ++k)
	      fprintf(scalar_file,"%e ",nodal_output.results[k]*sscale[i]);
            fprintf(scalar_file,"\n");
	    fclose(scalar_file);
          } else {

            this->com->fprintf(stderr,"Warning: Cannot open probe file %s",nodal_scalars[i]);
          }
	}
      }
    }

    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (nodal_vectors[i]) {
	
	postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U,
				      nodal_output.subId, nodal_output.locNodeId,
				      nodal_output.last,nodal_output.numNodes,nodal_output.results,
                                      nodal_output.locations,
				      fluidId, distLSS, ghostPoints);

	if (com->cpuNum() == 0) {
	  FILE* vector_file = fopen(nodal_vectors[i],mode);
          if (vector_file != 0) {
	    fprintf(vector_file,"%d %e ",nodal_output.step+it0, tag);
	    for (int k =0 ; k < nodal_output.numNodes; ++k)
	      fprintf(vector_file,"%e %e %e ",
	  	      nodal_output.results[k*3]*vscale[i],
		      nodal_output.results[k*3+1]*vscale[i],
		      nodal_output.results[k*3+2]*vscale[i]);
            fprintf(vector_file,"\n");
	    fclose(vector_file);
          } else {

            this->com->fprintf(stderr,"Warning: Cannot open probe file %s",nodal_vectors[i]);
          }
	}
      }
    }
    // }
    ++nodal_output.step;
}

//----------------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U, 
                                             DistTimeState<dim> *timeState,
                                             DistVec<int> &fluidId)
{
  writeBinaryVectorsToDisk(lastIt,it,t,X,A,U,timeState,fluidId, (DistSVec<double,1>*)0);
}

template<int dim>
void TsOutput<dim>::writeProbesToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
				      DistVec<double> &A, DistSVec<double,dim> &U, 
				      DistTimeState<dim> *timeState,
				      DistVec<int> &fluidId, DistLevelSetStructure *distLSS,
                                      DistVec<GhostPoint<dim>*> *ghostPoints)
{
  writeProbesToDisk(lastIt,it,t,X,A,U,timeState,fluidId, (DistSVec<double,1>*)0, distLSS, ghostPoints);
}

//----------------------------------------------------------------------------------------


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
/*
template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, 
                                             DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U,
                                             DistSVec<double,1> &Phi, DistVec<int> &fluidId)
{
 
  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
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
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, *Qs, Phi, fluidId);
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
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, fluidId);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}
*/
//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeAvgVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U, DistTimeState<dim> *timeState)
{

// This routine outputs the time averaged values of the scalar and vector output files
// in binary format
  int i;
  static double tprev,tinit;
  double time = refVal->time*t;
  double del_t;

  int step = getStep(it,lastIt,t);

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

  if (toWrite(it,lastIt,t) && counter>0) {
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
    // Before deletion, check that pointers have been allocated
    for (i=0; i<PostFcn::AVSSIZE; ++i) 
      if ((avscalars[i]) && (AvQs[i]))
        delete AvQs[i];
    for (i=0; i<PostFcn::AVVSIZE; ++i) 
      if ((avvectors[i]) && (AvQv[i]))
        delete AvQv[i];
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
  if (iod.output.transient.pressurecoefficient[0] != 0) {
    sscale[PostFcn::PRESSURECOEFFICIENT] = 1.0;
    scalars[PostFcn::PRESSURECOEFFICIENT] = new char[sp + strlen(iod.output.transient.pressurecoefficient)];
    sprintf(scalars[PostFcn::PRESSURECOEFFICIENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressurecoefficient);
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
  if (iod.output.transient.surfaceheatflux[0] != 0) {
    scalars[PostFcn::SURFACE_HEAT_FLUX] = new char[sp + strlen(iod.output.transient.surfaceheatflux)];
    sprintf(scalars[PostFcn::SURFACE_HEAT_FLUX], "%s%s",
            iod.output.transient.prefix, iod.output.transient.surfaceheatflux);
  }
  if (iod.output.transient.tempnormalderivative[0] != 0) {
    scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = new char[sp + strlen(iod.output.transient.tempnormalderivative)];
    sprintf(scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE], "%s%s",
             iod.output.transient.prefix, iod.output.transient.tempnormalderivative);
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
  if (iod.output.transient.fluidid[0] != 0) {
    sscale[PostFcn::FLUIDID] = 1.0;
    scalars[PostFcn::FLUIDID] = new char[sp + strlen(iod.output.transient.fluidid)];
    sprintf(scalars[PostFcn::FLUIDID], "%s%s",
            iod.output.transient.prefix, iod.output.transient.fluidid);
  }
  if (iod.output.transient.controlvolume[0] != 0) {
    sscale[PostFcn::CONTROL_VOLUME] = iod.ref.rv.length * iod.ref.rv.length * iod.ref.rv.length;
    scalars[PostFcn::CONTROL_VOLUME] = new char[sp + strlen(iod.output.transient.controlvolume)];
    sprintf(scalars[PostFcn::CONTROL_VOLUME], "%s%s",
            iod.output.transient.prefix, iod.output.transient.controlvolume);
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
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTotalpressure);
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.dTotalpressure);
  }

  int dsp = strlen(iod.output.transient.prefix) + 1;

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

template<int dim>
void TsOutput<dim>::writeStateRomToDisk(int it, double cpu, int nPod, const Vec<double> &UromTotal)
{

  if (com->cpuNum() != 0) return;

  if (fpStateRom) {
    fprintf(fpStateRom, "%d %e", it, cpu);
		for (int iPod = 0; iPod < nPod; ++iPod) {
			fprintf(fpStateRom, " %23.15e", UromTotal[iPod]);	// write out with high precision
		}
    fprintf(fpStateRom, "\n");
    fflush(fpStateRom);
  }

}

template<int dim>
void TsOutput<dim>::writeErrorToDisk(const int it, const double cpu, const int nErr, const double *error)
{

  if (com->cpuNum() != 0) return;

  if (fpError) {
    fprintf(fpError, "%d %e", it, cpu);
		for (int iErr = 0; iErr < nErr; ++iErr) {
			fprintf(fpError, " %23.15e", error[iErr]);	// write out with high precision
		}
    fprintf(fpError, "\n");
    fflush(fpError);
  }

}
