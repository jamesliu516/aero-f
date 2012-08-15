#include <IoData.h>

#include <Communicator.h>
#include <parser/Assigner.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <dlfcn.h>
#include <set>
#include <ProgrammedBurn.h>
using std::set;

#ifdef COUGAR
extern int optind;
extern "C" int getopt(int, char **, char *);
#endif

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------
FluidRemapData::FluidRemapData() {

  oldID = newID = -1;
}

void FluidRemapData::setup(const char * name, ClassAssigner * father) {
  
  ClassAssigner *ca = new ClassAssigner(name, 1, father);
  //new ClassInt<FluidRemapData>(ca, "OldID", this, &FluidRemapData::oldID);
  new ClassInt<FluidRemapData>(ca, "FluidIDReceptor", this, &FluidRemapData::newID);
}

Assigner* FluidRemapData::getAssigner() {
  
  ClassAssigner *ca = new ClassAssigner("normal", 1, nullAssigner);
  //new ClassInt<FluidRemapData>(ca, "OldID", this, &FluidRemapData::oldID);
  new ClassInt<FluidRemapData>(ca, "FluidIDReceptor", this, &FluidRemapData::newID);
  return ca;
}

OneDimensionalInputData::OneDimensionalInputData() {

  file = "";
  x0 = y0 = z0 = 0.0;
}

void OneDimensionalInputData::setup(const char * name, ClassAssigner * father) {

  ClassAssigner *ca = new ClassAssigner(name, 5, father);
  new ClassStr<OneDimensionalInputData>(ca, "File", this, &OneDimensionalInputData::file);
  new ClassDouble<OneDimensionalInputData>(ca, "X0", this, &OneDimensionalInputData::x0);
  new ClassDouble<OneDimensionalInputData>(ca, "Y0", this, &OneDimensionalInputData::y0);
  new ClassDouble<OneDimensionalInputData>(ca, "Z0", this, &OneDimensionalInputData::z0);
  
  fluidRemap.setup("FluidIDMap",ca);
}

Assigner* OneDimensionalInputData::getAssigner() {

  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);
  new ClassStr<OneDimensionalInputData>(ca, "File", this, &OneDimensionalInputData::file);
  new ClassDouble<OneDimensionalInputData>(ca, "X0", this, &OneDimensionalInputData::x0);
  new ClassDouble<OneDimensionalInputData>(ca, "Y0", this, &OneDimensionalInputData::y0);
  new ClassDouble<OneDimensionalInputData>(ca, "Z0", this, &OneDimensionalInputData::z0);
  
  fluidRemap.setup("FluidIDMap",ca);
  return ca;
}

InputData::InputData()
{

  prefix = "";
  geometryprefix = "";
  connectivity = "";
  geometry = "";
  decomposition = "";
  cpumap = "";
  match = "";
  d2wall = "";
  perturbed = "";
  solutions = "";
  positions = "";
  levelsets = "";
  rstdata = "";
	podFile = "";
	snapRefSolutionFile = "";
	staterom = "";
	snapFile = "";

	// gnat
	gnatPrefix = "";
	sampleNodes = "";
	resMatrix = "";
	jacMatrix = "";
  podFileState = "";
  podFileRes = "";
  podFileJac = "";
  podFileResHat = "";
  podFileJacHat = "";
  mesh = "";
  reducedfullnodemap = "";

  convergence_file = "";

// Included (MB)
  shapederivatives = "";
  strModesFile = "";
  embeddedSurface= "";
  oneDimensionalSolution = "";

  exactInterfaceLocation = "";
}

//------------------------------------------------------------------------------

void InputData::setup(const char *name, ClassAssigner *father)
{

// Modified (MB)
  ClassAssigner *ca = new ClassAssigner(name, 32, father);
  new ClassStr<InputData>(ca, "Prefix", this, &InputData::prefix);
  new ClassStr<InputData>(ca, "GeometryPrefix", this, &InputData::geometryprefix);
  new ClassStr<InputData>(ca, "Connectivity", this, &InputData::connectivity);
  new ClassStr<InputData>(ca, "Geometry", this, &InputData::geometry);
  new ClassStr<InputData>(ca, "Decomposition", this, &InputData::decomposition);
  new ClassStr<InputData>(ca, "CpuMap", this, &InputData::cpumap);
  new ClassStr<InputData>(ca, "Matcher", this, &InputData::match);
  new ClassStr<InputData>(ca, "WallDistance", this, &InputData::d2wall);
  new ClassStr<InputData>(ca, "Perturbed", this, &InputData::perturbed);
  new ClassStr<InputData>(ca, "Solution", this, &InputData::solutions);
  new ClassStr<InputData>(ca, "Position", this, &InputData::positions);
  new ClassStr<InputData>(ca, "LevelSet", this, &InputData::levelsets);
  new ClassStr<InputData>(ca, "RestartData", this, &InputData::rstdata);
  new ClassStr<InputData>(ca, "PODData", this, &InputData::podFile);
  new ClassStr<InputData>(ca, "SnapshotData", this, &InputData::snapFile);
  new ClassStr<InputData>(ca, "SnapshotsReferenceSolution", this, &InputData::snapRefSolutionFile);
  new ClassStr<InputData>(ca, "ROBStateCoordinates", this, &InputData::staterom);

	// sample mesh
  new ClassStr<InputData>(ca, "GNATPrefix", this, &InputData::gnatPrefix);
  new ClassStr<InputData>(ca, "SampledNodes", this, &InputData::sampleNodes);
  new ClassStr<InputData>(ca, "GappyJacMat", this, &InputData::jacMatrix);
  new ClassStr<InputData>(ca, "GappyResMat", this, &InputData::resMatrix);
  new ClassStr<InputData>(ca, "ROBState", this, &InputData::podFileState);
  new ClassStr<InputData>(ca, "ROBResidual", this, &InputData::podFileRes);
  new ClassStr<InputData>(ca, "ROBJacobian", this, &InputData::podFileJac);
  new ClassStr<InputData>(ca, "SampledROBResidual", this, &InputData::podFileResHat);
  new ClassStr<InputData>(ca, "SampledROBJacobian", this, &InputData::podFileJacHat);
	// ???
  new ClassStr<InputData>(ca, "ReducedMesh", this, &InputData::mesh);
  new ClassStr<InputData>(ca, "ReducedFullNodeMap", this, &InputData::reducedfullnodemap);

// Included (MB)
  new ClassStr<InputData>(ca, "ShapeDerivative", this, &InputData::shapederivatives);
  new ClassStr<InputData>(ca, "StrModes", this, &InputData::strModesFile);

  new ClassStr<InputData>(ca, "EmbeddedSurface", this, &InputData::embeddedSurface);
  new ClassStr<InputData>(ca, "ConvergenceFile", this, &InputData::convergence_file);
  new ClassStr<InputData>(ca, "ExactInterfaceLocation", this, &InputData::exactInterfaceLocation);

  oneDimensionalInput.setup("1DRestartData",ca);

}

//------------------------------------------------------------------------------

PreconditionData::PreconditionData()
{
  mach = 1.0;
  k = 1.0;
  cmach = 1.0;
  betav = 0.0;
  shockreducer = 0.0;
}

//------------------------------------------------------------------------------

void PreconditionData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name,5,father);

  new ClassDouble<PreconditionData>(ca,"Mach", this, &PreconditionData::mach);
  new ClassDouble<PreconditionData>(ca,"CutOffMach", this, &PreconditionData::cmach);
  new ClassDouble<PreconditionData>(ca,"k", this, &PreconditionData::k);
  new ClassDouble<PreconditionData>(ca,"Betav", this, &PreconditionData::betav);
  new ClassDouble<PreconditionData>(ca,"ShockReducer", this, &PreconditionData::shockreducer);

}

//------------------------------------------------------------------------------
OutputData::OutputData()
{

}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  transient.setup("Postpro", ca);
  restart.setup("Restart", ca);
  transient.probes.setup("Probes", ca);
  rom.setup("NonlinearROM", ca);
}

//------------------------------------------------------------------------------

void Probes::Node::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassInt<Probes::Node>(ca, "ID", this, &Probes::Node::id);
  new ClassDouble<Probes::Node>(ca, "LocationX",this,&Probes::Node::locationX);
  new ClassDouble<Probes::Node>(ca, "LocationY",this,&Probes::Node::locationY);
  new ClassDouble<Probes::Node>(ca, "LocationZ",this,&Probes::Node::locationZ);
}

Probes::Probes() {

  prefix = "";
  density = "";
  pressure = "";
  temperature = "";
  velocity = "";
  displacement = "";
}

void Probes::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 56, father);
  new ClassStr<Probes>(ca, "Prefix", this, &Probes::prefix);
  new ClassStr<Probes>(ca, "Density", this, &Probes::density);
  new ClassStr<Probes>(ca, "Pressure", this, &Probes::pressure);
  new ClassStr<Probes>(ca, "Temperature", this, &Probes::temperature);
  new ClassStr<Probes>(ca, "Velocity", this, &Probes::velocity);
  new ClassStr<Probes>(ca, "Displacement", this, &Probes::displacement);

  char nodename[12];
  for (int i = 0; i < MAXNODES; ++i) {
    sprintf(nodename,"Node%d",i+1);
    
    myNodes[i].setup(nodename, ca);
  }

}

TransientData::TransientData()
{

  prefix = "";
  solutions = "";
  density = "";
  tavdensity = "";
  mach = "";
  speed = "";
  wtmach = "";
  wtspeed = "";
  absvelocity = "";
  tavmach = "";
  pressure = "";
  diffpressure = "";
  tavpressure = "";
  hydrostaticpressure = "";
  hydrodynamicpressure = "";
  pressurecoefficient = "";
  temperature = "";
  tavtemperature = "";
  totalpressure = "";
  tavtotalpressure = "";
  vorticity = "";
  tavvorticity = "";
  nutturb = "";
  kturb = "";
  epsturb = "";
  eddyvis = "";
  dplus = "";
  sfric = "";
  tavsfric = "";
  psensor = "";
  csdles = "";
  tavcsdles = "";
  csdvms = "";
  tavcsdvms = "";
  mutOmu = "";
  velocity = "";
  tavvelocity = "";
  displacement = "";
  tavdisplacement = "";
  flightDisplacement = "";
  localFlightDisplacement = "";
  forces = "";
  tavforces = "";
  hydrostaticforces = "";
  hydrodynamicforces = "";
  generalizedforces = "";
  lift = "";
  tavlift = "";
  hydrostaticlift = "";
  hydrodynamiclift = "";
  residuals = "";
  materialVolumes = "";
  conservation = "";
  podFile = "";
  romFile = "";
  robProductFile = "";
  gendispFile = "";
  philevel = "";
  controlvolume = "";
  fluidid="";
  embeddedsurface = "";
  cputiming = "";

// Included (MB)
  velocitynorm = "";
  dSolutions = "";
  dDensity = "";
  dMach = "";
  dPressure = "";
  dTotalpressure = "";
  dTemperature = "";
  dNutturb = "";
  dVelocityScalar = "";
  dVelocityVector = "";
  dDisplacement = "";
  dForces = "";
  dEddyvis = "";

  tempnormalderivative = "";
  surfaceheatflux = "";
  heatfluxes = "";
  sparseGrid = "SparseGrid";

  bubbleRadius = "";

  frequency = 0;
  frequency_dt = -1.0;
  length = 1.0;
  surface = 1.0;
  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;

}

//------------------------------------------------------------------------------

void TransientData::setup(const char *name, ClassAssigner *father)
{

// Modified (MB)
  ClassAssigner *ca = new ClassAssigner(name, 85, father); 

  new ClassStr<TransientData>(ca, "Prefix", this, &TransientData::prefix);
  new ClassStr<TransientData>(ca, "StateVector", this, &TransientData::solutions);
  new ClassStr<TransientData>(ca, "Density", this, &TransientData::density);
  new ClassStr<TransientData>(ca, "TavDensity", this, &TransientData::density);
  new ClassStr<TransientData>(ca, "Mach", this, &TransientData::mach);
  new ClassStr<TransientData>(ca, "VelocityMagnitude", this, &TransientData::speed);
  new ClassStr<TransientData>(ca, "HWTMach", this, &TransientData::wtmach);
  new ClassStr<TransientData>(ca, "HWTVelocityMagnitude", this, &TransientData::wtspeed);
  new ClassStr<TransientData>(ca, "AbsVelocity", this, &TransientData::absvelocity);
  new ClassStr<TransientData>(ca, "TavMach", this, &TransientData::tavmach);

  new ClassStr<TransientData>(ca, "Pressure", this, &TransientData::pressure);
  new ClassStr<TransientData>(ca, "DeltaPressure", this, &TransientData::diffpressure);
  new ClassStr<TransientData>(ca, "HydroStaticPressure", this, &TransientData::hydrostaticpressure);
  new ClassStr<TransientData>(ca, "HydroDynamicPressure", this, &TransientData::hydrodynamicpressure);
  new ClassStr<TransientData>(ca, "PressureCoefficient", this, &TransientData::pressurecoefficient);
  new ClassStr<TransientData>(ca, "TavPressure", this, &TransientData::tavpressure);
  new ClassStr<TransientData>(ca, "Temperature", this, &TransientData::temperature);
  new ClassStr<TransientData>(ca, "TavTemperature", this, &TransientData::tavtemperature);
  new ClassStr<TransientData>(ca, "TotalPressure", this, &TransientData::totalpressure);
  new ClassStr<TransientData>(ca, "TavTotalPressure", this, &TransientData::tavtotalpressure);
  new ClassStr<TransientData>(ca, "Vorticity", this, &TransientData::vorticity);

  new ClassStr<TransientData>(ca, "TavVorticity", this, &TransientData::tavvorticity);
  new ClassStr<TransientData>(ca, "NuTilde", this, &TransientData::nutturb);
  new ClassStr<TransientData>(ca, "K", this, &TransientData::kturb);
  new ClassStr<TransientData>(ca, "Eps", this, &TransientData::epsturb);
  new ClassStr<TransientData>(ca, "EddyViscosity", this, &TransientData::eddyvis);
  new ClassStr<TransientData>(ca, "DeltaPlus", this, &TransientData::dplus);
  new ClassStr<TransientData>(ca, "SkinFriction", this, &TransientData::sfric);
  new ClassStr<TransientData>(ca, "TavSkinFriction", this, &TransientData::tavsfric);
  new ClassStr<TransientData>(ca, "PressureSensor", this, &TransientData::psensor);
  new ClassStr<TransientData>(ca, "CsDLES", this, &TransientData::csdles);
  new ClassStr<TransientData>(ca, "TavCsDLES", this, &TransientData::tavcsdles);
  new ClassStr<TransientData>(ca, "CsDVMS", this, &TransientData::csdvms);
  new ClassStr<TransientData>(ca, "TavCsDVMS", this, &TransientData::tavcsdvms);
  new ClassStr<TransientData>(ca, "MutOverMu", this, &TransientData::mutOmu);

  new ClassStr<TransientData>(ca, "Velocity", this, &TransientData::velocity);
  new ClassStr<TransientData>(ca, "TavVelocity", this, &TransientData::tavvelocity);
  new ClassStr<TransientData>(ca, "Displacement", this, &TransientData::displacement);
  new ClassStr<TransientData>(ca, "TavDisplacement", this, &TransientData::tavdisplacement);
  new ClassStr<TransientData>(ca, "FlightDisplacement", this, &TransientData::flightDisplacement);
  new ClassStr<TransientData>(ca, "LocalFlightDisplacement", this, &TransientData::localFlightDisplacement);
  new ClassStr<TransientData>(ca, "Force", this, &TransientData::forces);
  new ClassStr<TransientData>(ca, "TavForce", this, &TransientData::tavforces);
  new ClassStr<TransientData>(ca, "HydroStaticForce", this, &TransientData::hydrostaticforces);
  new ClassStr<TransientData>(ca, "HydroDynamicForce", this, &TransientData::hydrodynamicforces);

  new ClassStr<TransientData>(ca, "GeneralizedForce", this, &TransientData::generalizedforces);
  new ClassStr<TransientData>(ca, "LiftandDrag", this, &TransientData::lift);
  new ClassStr<TransientData>(ca, "HydroStaticLiftandDrag", this, &TransientData::hydrostaticlift);
  new ClassStr<TransientData>(ca, "HydroDynamicLiftandDrag", this, &TransientData::hydrodynamiclift);
  new ClassStr<TransientData>(ca, "TavLiftandDrag", this, &TransientData::tavlift);
  new ClassStr<TransientData>(ca, "Residual", this, &TransientData::residuals);
  new ClassStr<TransientData>(ca, "MaterialVolumes", this, &TransientData::materialVolumes);
  new ClassInt<TransientData>(ca, "Frequency", this, &TransientData::frequency);
  new ClassDouble<TransientData>(ca, "TimeInterval", this, &TransientData::frequency_dt);
  new ClassDouble<TransientData>(ca, "Length", this, &TransientData::length);
  new ClassDouble<TransientData>(ca, "Surface", this, &TransientData::surface);
  new ClassDouble<TransientData>(ca, "XM", this, &TransientData::x0);

  new ClassDouble<TransientData>(ca, "YM", this, &TransientData::y0);
  new ClassDouble<TransientData>(ca, "ZM", this, &TransientData::z0);
  new ClassStr<TransientData>(ca, "PODData", this, &TransientData::podFile);
  new ClassStr<TransientData>(ca, "ROM", this, &TransientData::romFile);
  new ClassStr<TransientData>(ca, "ROBInnerProducts", this, &TransientData::robProductFile);
  new ClassStr<TransientData>(ca, "GeneralizedDisplacement", this, &TransientData::gendispFile);
  new ClassStr<TransientData>(ca, "Philevel", this, &TransientData::philevel);
  new ClassStr<TransientData>(ca, "ConservationErrors", this, &TransientData::conservation);
  new ClassStr<TransientData>(ca, "FluidID", this, &TransientData::fluidid);
  new ClassStr<TransientData>(ca, "ControlVolume", this, &TransientData::controlvolume);
  new ClassStr<TransientData>(ca, "EmbeddedSurfaceDisplacement", this, &TransientData::embeddedsurface);
  new ClassStr<TransientData>(ca, "CPUTiming", this, &TransientData::cputiming);

	// Gappy POD offline
	// Gappy POD snapshots
// Included (MB)
  new ClassStr<TransientData>(ca, "VelocityNorm", this, &TransientData::velocitynorm);
  new ClassStr<TransientData>(ca, "StateVectorSensitivity", this, &TransientData::dSolutions); //KW(Aug.17,2010): used to be SolutionSensitivity
  new ClassStr<TransientData>(ca, "DensitySensitivity", this, &TransientData::dDensity);
  new ClassStr<TransientData>(ca, "MachSensitivity", this, &TransientData::dMach);
  new ClassStr<TransientData>(ca, "PressureSensitivity", this, &TransientData::dPressure);

  new ClassStr<TransientData>(ca, "TemperatureSensitivity", this, &TransientData::dTemperature);
  new ClassStr<TransientData>(ca, "TotalPressureSensitivity", this, &TransientData::dTotalpressure);
  new ClassStr<TransientData>(ca, "NuTildeSensitivity", this, &TransientData::dNutturb);
  new ClassStr<TransientData>(ca, "EddyViscositySensitivity", this, &TransientData::dEddyvis);
  new ClassStr<TransientData>(ca, "VelocityNormSensitivity", this, &TransientData::dVelocityScalar);
  new ClassStr<TransientData>(ca, "VelocitySensitivity", this, &TransientData::dVelocityVector);
  new ClassStr<TransientData>(ca, "DisplacementSensitivity", this, &TransientData::dDisplacement);
  new ClassStr<TransientData>(ca, "ForceSensitivity", this, &TransientData::dForces);

  new ClassStr<TransientData>(ca, "TemperatureNormalDerivative", this, &TransientData::tempnormalderivative);
  new ClassStr<TransientData>(ca, "HeatFluxPerUnitSurface", this, &TransientData::surfaceheatflux); 
  new ClassStr<TransientData>(ca, "HeatFlux", this, &TransientData::heatfluxes);
  new ClassStr<TransientData>(ca, "SparseGrid", this, &TransientData::sparseGrid);

  new ClassStr<TransientData>(ca, "BubbleRadius", this, &TransientData::bubbleRadius);

	//do defaults

}

ROMOutputData::ROMOutputData()
{
  prefix = "";

  newtonresiduals = "";
  jacobiandeltastate = "";
  reducedjac = "";
  stateRom = "";

  gnatPrefix = "";

  sampleNodes = "";
  onlineMatrix = "";
  podStateRed = "";
  podNonlinRed = "";
  solution = "";
  wallDistanceRed = "";
  staterom = "";
  error = "";
  dUnormAccum = "";

  sampleNodesFull = "";
  onlineMatrixFull = "";
  mesh = "";
  reducedfullnodemap = "";

  frequency = 0;
  frequency_dt = -1.0;
// Gappy POD

  //mesh = "";
  //reducedfullnodemap = "";
  //sampleNodes = "";
  //sampleNodesFull = "";
  //onlineMatrix = "";
  //onlineMatrixFull = "";
  //podStateRed = "";
  //wallDistanceRed = "";
  //newtonresiduals = "";
  //jacobiandeltastate = "";
  //reducedjac = "";
  //staterom = "";
  //error = "";
  //dUnormAccum = "";
}

void ROMOutputData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 18, father); 
  new ClassStr<ROMOutputData>(ca, "Prefix", this, &ROMOutputData::prefix);
	
  new ClassStr<ROMOutputData>(ca, "Residual", this, &ROMOutputData::newtonresiduals);
  new ClassStr<ROMOutputData>(ca, "JacobianDeltaState", this, &ROMOutputData::jacobiandeltastate);
  new ClassStr<ROMOutputData>(ca, "ReducedJac", this, &ROMOutputData::reducedjac);

  new ClassStr<ROMOutputData>(ca, "GNATPrefix", this, &ROMOutputData::gnatPrefix);
  new ClassStr<ROMOutputData>(ca, "StateReducedCoordinates", this, &ROMOutputData::staterom);
  new ClassStr<ROMOutputData>(ca, "SampledNodes", this, &ROMOutputData::sampleNodes);
  new ClassStr<ROMOutputData>(ca, "GappyPODMatrix", this, &ROMOutputData::onlineMatrix);
  new ClassStr<ROMOutputData>(ca, "SampledROBState", this, &ROMOutputData::podStateRed);
  new ClassStr<ROMOutputData>(ca, "SampledROBNonlinear", this, &ROMOutputData::podNonlinRed);
  new ClassStr<ROMOutputData>(ca, "SampledSolution", this, &ROMOutputData::solution);
  new ClassStr<ROMOutputData>(ca, "SampledWallDistance", this, &ROMOutputData::wallDistanceRed);
  new ClassStr<ROMOutputData>(ca, "Error", this, &ROMOutputData::error);
  new ClassStr<ROMOutputData>(ca, "NetStateReducedCoordinates", this, &ROMOutputData::dUnormAccum);
  new ClassStr<ROMOutputData>(ca, "SampledNodesFullMesh", this, &ROMOutputData::sampleNodesFull);
  new ClassStr<ROMOutputData>(ca, "GappyPODMatrixFullMesh", this, &ROMOutputData::onlineMatrixFull);
  new ClassStr<ROMOutputData>(ca, "SampledMesh", this, &ROMOutputData::mesh);
  new ClassStr<ROMOutputData>(ca, "SampledFullNodeMap", this, &ROMOutputData::reducedfullnodemap);
}

//------------------------------------------------------------------------------

RestartData::RestartData()
{

  type = SINGLE;
  prefix = "";
  solutions = "DEFAULT.SOL";
  positions = "DEFAULT.POS";
  levelsets= "DEFAULT.LEV";
  data = "DEFAULT.RST";


  frequency = 0;
  frequency_dt = -1.0;

}

//------------------------------------------------------------------------------

void RestartData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);
  new ClassToken<RestartData>(ca, "Type", this,
			      reinterpret_cast<int RestartData::*>(&RestartData::type), 2,
			      "Single", 0, "Double", 1);

  new ClassStr<RestartData>(ca, "Prefix", this, &RestartData::prefix);
  new ClassStr<RestartData>(ca, "Solution", this, &RestartData::solutions);
  new ClassStr<RestartData>(ca, "Position", this, &RestartData::positions);
  new ClassStr<RestartData>(ca, "LevelSet", this, &RestartData::levelsets);
  new ClassStr<RestartData>(ca, "RestartData", this, &RestartData::data);
  new ClassInt<RestartData>(ca, "Frequency", this, &RestartData::frequency);
  new ClassDouble<RestartData>(ca, "TimeInterval", this, &RestartData::frequency_dt);

}

//------------------------------------------------------------------------------

RestartParametersData::RestartParametersData()
{

  iteration = 0;

  etime = 0.0;
  dt_nm1 = 1.0;
  dt_nm2 = 1.0;
  residual = 1.0;
  energy = 0.0;
  output_newton_step = 0;

}

//------------------------------------------------------------------------------

void RestartParametersData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassInt<RestartParametersData>(ca, "Iteration", this, &RestartParametersData::iteration);
  new ClassDouble<RestartParametersData>(ca, "Time", this, &RestartParametersData::etime);
  new ClassDouble<RestartParametersData>(ca, "TimeStep1", this, &RestartParametersData::dt_nm1);
  new ClassDouble<RestartParametersData>(ca, "TimeStep2", this, &RestartParametersData::dt_nm2);
  new ClassDouble<RestartParametersData>(ca, "Residual", this, &RestartParametersData::residual);
  new ClassDouble<RestartParametersData>(ca, "Energy", this, &RestartParametersData::energy);
  new ClassInt<RestartParametersData>(ca, "NewtonOutputStep", this, &RestartParametersData::output_newton_step);

}

//------------------------------------------------------------------------------
/*
  0: only time iterations
  1: warnings
  2: initialization step + cpu report
  3: file reading
  4: file writing
  5: status of linear solver + global dt + rotation angle
  6: timings for jacobians and preconditioners
  7: send + receive
  8: #cpus
*/

ProblemData::ProblemData()
{

  alltype = _STEADY_;
  mode = NON_DIMENSIONAL;
  prec = NON_PRECONDITIONED;
  framework = BODYFITTED;
  solvefluid = ON;

  test = REGULAR;
  verbose = 4;

}

//------------------------------------------------------------------------------

void ProblemData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);
  new ClassToken<ProblemData>
    (ca, "Type", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::alltype), 33,
     "Steady", 0, "Unsteady", 1, "AcceleratedUnsteady", 2, "SteadyAeroelastic", 3,
     "UnsteadyAeroelastic", 4, "AcceleratedUnsteadyAeroelastic", 5,
     "SteadyAeroThermal", 6, "UnsteadyAeroThermal", 7, "SteadyAeroThermoElastic", 8,
		 "UnsteadyAeroThermoElastic", 9, "Forced", 10, "AcceleratedForced", 11,
		 "RigidRoll", 12, "RbmExtractor", 13, "UnsteadyLinearizedAeroelastic", 14,
		 "UnsteadyLinearized", 15, "ROBConstruction", 16, "ROMAeroelastic", 17,
		 "ROM", 18, "ForcedLinearized", 19, "PODInterpolation", 20,
		 "SteadySensitivityAnalysis", 21, "SparseGridGeneration", 22,
		 "1D", 23, "NonlinearROM", 24, "NonlinearROMPreprocessing", 25,
		 "NonlinearROMSurfaceMeshConstruction",26, "SampledMeshShapeChange", 27,
		 "NonlinearROMPreprocessingStep1", 28, "NonlinearROMPreprocessingStep2", 29,
		 "NonlinearROMPostprocessing", 30, "PODConstruction", 31, "ROBInnerProduct", 32);

  new ClassToken<ProblemData>
    (ca, "Mode", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::mode), 2,
     "NonDimensional", 0, "Dimensional", 1);

  new ClassToken<ProblemData>
    (ca, "Prec", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::prec), 2,
     "NonPreconditioned", 0, "LowMach", 1);

  new ClassToken<ProblemData>
    (ca, "Framework", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::framework), 2,
     "BodyFitted", 0, "Embedded", 1);

  new ClassToken<ProblemData>
    (ca, "SolveFluid", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::solvefluid), 2,
     "Off", 0, "On", 1);

  new ClassToken<ProblemData>
    (ca, "Test", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::test), 1,
     "Regular", 0);

}

//------------------------------------------------------------------------------

ReferenceStateData::ReferenceStateData()
{

  mach = -1.0;
  velocity = -1.0;
  density = -1.0;
  pressure = -1.0;
  temperature = -1.0;
  reynolds_mu = -1.0;
  length = 1.0;

}

//------------------------------------------------------------------------------

void ReferenceStateData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassDouble<ReferenceStateData>(ca, "Mach", this, &ReferenceStateData::mach);
  new ClassDouble<ReferenceStateData>(ca, "Velocity", this, &ReferenceStateData::velocity);
  new ClassDouble<ReferenceStateData>(ca, "Density", this, &ReferenceStateData::density);
  new ClassDouble<ReferenceStateData>(ca, "Pressure", this, &ReferenceStateData::pressure);
  new ClassDouble<ReferenceStateData>(ca, "Temperature", this, &ReferenceStateData::temperature);
  new ClassDouble<ReferenceStateData>(ca, "Reynolds", this, &ReferenceStateData::reynolds_mu);
  new ClassDouble<ReferenceStateData>(ca, "Length", this, &ReferenceStateData::length);

}

//------------------------------------------------------------------------------

BcsFreeStreamData::BcsFreeStreamData()
{

  type = EXTERNAL;
  mach = -1.0;
  velocity = -1.0;
  alpha = 400.0;
  beta = 400.0;
  density = -1.0;
  pressure = -1.0;
  temperature = -1.0;
  nutilde = -1.0;
  kenergy = -1.0;
  eps = -1.0;

}

//------------------------------------------------------------------------------

void BcsFreeStreamData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 11, father);

  new ClassToken<BcsFreeStreamData>
    (ca, "Type", this, reinterpret_cast<int BcsFreeStreamData::*>(&BcsFreeStreamData::type), 2,
     "External", 0, "Internal", 1);
  new ClassDouble<BcsFreeStreamData>(ca, "Mach", this, &BcsFreeStreamData::mach);
  new ClassDouble<BcsFreeStreamData>(ca, "Velocity", this, &BcsFreeStreamData::velocity);
  new ClassDouble<BcsFreeStreamData>(ca, "Alpha", this, &BcsFreeStreamData::alpha);
  new ClassDouble<BcsFreeStreamData>(ca, "Beta", this, &BcsFreeStreamData::beta);
  new ClassDouble<BcsFreeStreamData>(ca, "Density", this, &BcsFreeStreamData::density);
  new ClassDouble<BcsFreeStreamData>(ca, "Pressure", this, &BcsFreeStreamData::pressure);
  new ClassDouble<BcsFreeStreamData>(ca, "Temperature", this, &BcsFreeStreamData::temperature);
  new ClassDouble<BcsFreeStreamData>(ca, "NuTilde", this, &BcsFreeStreamData::nutilde);
  new ClassDouble<BcsFreeStreamData>(ca, "K", this, &BcsFreeStreamData::kenergy);
  new ClassDouble<BcsFreeStreamData>(ca, "Eps", this, &BcsFreeStreamData::eps);

}

//------------------------------------------------------------------------------

BcsWallData::BcsWallData()
{

  type = ADIABATIC;
  integration = AUTO;
  reconstruction = CONSTANT;
  temperature = -1.0;
  delta = -1.0;

}

//------------------------------------------------------------------------------

void BcsWallData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<BcsWallData>(ca, "Type", this,
			      reinterpret_cast<int BcsWallData::*>(&BcsWallData::type), 2,
			      "Isothermal", 0, "Adiabatic", 1);
  new ClassToken<BcsWallData>(ca, "Integration", this,
			      reinterpret_cast<int BcsWallData::*>(&BcsWallData::integration), 2,
			      "WallFunction", 1, "Full", 2);
  new ClassToken<BcsWallData>(ca, "Reconstruction", this,
			      reinterpret_cast<int BcsWallData::*>(&BcsWallData::reconstruction), 2,
			      "Constant", 0, "ExactRiemann", 1);
  new ClassDouble<BcsWallData>(ca, "Temperature", this, &BcsWallData::temperature);
  new ClassDouble<BcsWallData>(ca, "Delta", this, &BcsWallData::delta);

}

//------------------------------------------------------------------------------

BcsHydroData::BcsHydroData()
{

  depth   = 0.0;

}

//------------------------------------------------------------------------------

void BcsHydroData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassDouble<BcsHydroData>(ca, "Depth", this, &BcsHydroData::depth);

}

//------------------------------------------------------------------------------

BcsData::BcsData()
{

}

//------------------------------------------------------------------------------

void BcsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  inlet.setup("Inlet", ca);
  outlet.setup("Outlet", ca);
  wall.setup("Wall", ca);
  hydro.setup("Hydro", ca);

}

//------------------------------------------------------------------------------

GasModelData::GasModelData()
{

  type = IDEAL;
  specificHeatRatio = 1.4;
  idealGasConstant = 287.1;
  specificHeatPressure = -1.0;
  pressureConstant = 0.0;

}

//------------------------------------------------------------------------------

void GasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<GasModelData>(ca, "Type", this,
			       reinterpret_cast<int GasModelData::*>(&GasModelData::type), 2,
			       "Ideal", 0, "Stiffened", 1);
  new ClassDouble<GasModelData>(ca, "SpecificHeatRatio", this,
                                &GasModelData::specificHeatRatio);
  new ClassDouble<GasModelData>(ca, "IdealGasConstant", this,
                                &GasModelData::idealGasConstant);
  new ClassDouble<GasModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                &GasModelData::specificHeatPressure);
  new ClassDouble<GasModelData>(ca, "PressureConstant", this,
                                &GasModelData::pressureConstant);

}

//------------------------------------------------------------------------------

JWLModelData::JWLModelData()
{

  type = JWL;
  omega = 0.4;
  idealGasConstant = 287.1;
  A1 = 0.0;
  A2 = 0.0;
  R1 = 1.0;
  R2 = 1.0;
  rhoref = 1.0;

}

//------------------------------------------------------------------------------

void JWLModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<JWLModelData>(ca, "Type", this,
			       reinterpret_cast<int JWLModelData::*>(&JWLModelData::type), 2,
			       "Ideal", 0, "JWL", 1);
  new ClassDouble<JWLModelData>(ca, "Omega", this,
                                &JWLModelData::omega);
  new ClassDouble<JWLModelData>(ca, "IdealGasConstant", this,
                                &JWLModelData::idealGasConstant);
  new ClassDouble<JWLModelData>(ca, "A1", this, &JWLModelData::A1);
  new ClassDouble<JWLModelData>(ca, "A2", this, &JWLModelData::A2);
  new ClassDouble<JWLModelData>(ca, "R1", this, &JWLModelData::R1);
  new ClassDouble<JWLModelData>(ca, "R2", this, &JWLModelData::R2);
  new ClassDouble<JWLModelData>(ca, "ReferenceDensity", this, &JWLModelData::rhoref);

}

//------------------------------------------------------------------------------

LiquidModelData::LiquidModelData()
{
  //data is available in 'Properties of Water and Steam', Wagner & Kruse, Springer, 1997
  type = COMPRESSIBLE;
  check = YES;
  specificHeat = -1.0;
  k1water     = 2.07e9;
  k2water     = 7.15;
  Prefwater   = -1.0;
  RHOrefwater = -1.0;
  burnable = NO;

  //these parameters are the adimensionalized parameters used by the VarFcn where the
  // liquid state equation is implemented. It it were dimensionalized, that is what their
  //value would be.
//  Pref  = -k1water/k2water;
//  alpha = (Prefwater + k1water/k2water)/pow(RHOrefwater, k2water);
//  beta  = k2water;

}

//------------------------------------------------------------------------------

void LiquidModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassToken<LiquidModelData>(ca, "Type", this,
            reinterpret_cast<int LiquidModelData::*>(&LiquidModelData::type), 1,
            "Compressible", 0);
  new ClassToken<LiquidModelData>(ca, "Check", this,
            reinterpret_cast<int LiquidModelData::*>(&LiquidModelData::check), 2,
            "Yes", 0, "No", 1);
  new ClassToken<LiquidModelData>(ca, "Burnable", this,
            reinterpret_cast<int LiquidModelData::*>(&LiquidModelData::burnable), 2,
            "Yes", 0, "No", 1);
  new ClassDouble<LiquidModelData>(ca, "k1", this, &LiquidModelData::k1water);
  new ClassDouble<LiquidModelData>(ca, "k2", this, &LiquidModelData::k2water);
  new ClassDouble<LiquidModelData>(ca, "Pressure", this, &LiquidModelData::Prefwater);
  new ClassDouble<LiquidModelData>(ca, "Density", this, &LiquidModelData::RHOrefwater);
  new ClassDouble<LiquidModelData>(ca, "SpecificHeat", this, &LiquidModelData::specificHeat);

}

//------------------------------------------------------------------------------

FluidModelData::FluidModelData()
{

  fluid = PERFECT_GAS;
  rhomin = -1.e9;
  pmin = -1.e9;

}

//------------------------------------------------------------------------------

Assigner *FluidModelData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 6, nullAssigner);

  new ClassToken<FluidModelData>(ca, "Fluid", this,
                                 reinterpret_cast<int FluidModelData::*>(&FluidModelData::fluid), 4,
                                 "PerfectGas",   FluidModelData::PERFECT_GAS,
				 "Liquid",       FluidModelData::LIQUID,
				 "StiffenedGas", FluidModelData::STIFFENED_GAS,
				 "JWL",          FluidModelData::JWL);
  new ClassDouble<FluidModelData>(ca, "DensityCutOff", this, &FluidModelData::rhomin);
  new ClassDouble<FluidModelData>(ca, "PressureCutOff", this, &FluidModelData::pmin);

  gasModel.setup("GasModel", ca);
  jwlModel.setup("JWLModel", ca);
  liquidModel.setup("LiquidModel", ca);

  return ca;

};

//------------------------------------------------------------------------------

void FluidModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<FluidModelData>(ca, "Fluid", this,
                                 reinterpret_cast<int FluidModelData::*>(&FluidModelData::fluid), 4,
                                 "PerfectGas",   FluidModelData::PERFECT_GAS,
				 "Liquid",       FluidModelData::LIQUID,
				 "StiffenedGas", FluidModelData::STIFFENED_GAS,
				 "JWL",          FluidModelData::JWL);
  new ClassDouble<FluidModelData>(ca, "DensityCutOff", this, &FluidModelData::rhomin);
  new ClassDouble<FluidModelData>(ca, "PressureCutOff", this, &FluidModelData::pmin);

  gasModel.setup("GasModel", ca);
  jwlModel.setup("JWLModel", ca);
  liquidModel.setup("LiquidModel", ca);

};

//------------------------------------------------------------------------------

/*
FluidModelData & FluidModelData::operator=(const FluidModelData &fm){

  this->fluid = fm.fluid;
  this->rhomin  = fm.rhomin;
  this->pmin  = fm.pmin;

  this->gasModel.type              = fm.gasModel.type;
  this->gasModel.specificHeatRatio = fm.gasModel.specificHeatRatio;
  this->gasModel.idealGasConstant  = fm.gasModel.idealGasConstant;
  this->gasModel.pressureConstant  = fm.gasModel.pressureConstant;

  this->jwlModel.type              = fm.jwlModel.type;
  this->jwlModel.omega             = fm.jwlModel.omega;
  this->jwlModel.idealGasConstant  = fm.jwlModel.idealGasConstant;
  this->jwlModel.A1                = fm.jwlModel.A1;
  this->jwlModel.A2                = fm.jwlModel.A2;
  this->jwlModel.R1                = fm.jwlModel.R1;
  this->jwlModel.R2                = fm.jwlModel.R2;
  this->jwlModel.rhoref            = fm.jwlModel.rhoref;

  this->liquidModel.type           = fm.liquidModel.type;
  this->liquidModel.check          = fm.liquidModel.check;
  this->liquidModel.specificHeatRatio = fm.liquidModel.specificHeatRatio;
  this->liquidModel.Cv             = fm.liquidModel.Cv;
  this->liquidModel.k1water        = fm.liquidModel.k1water;
  this->liquidModel.k2water        = fm.liquidModel.k2water;
  this->liquidModel.Prefwater      = fm.liquidModel.Prefwater;
  this->liquidModel.RHOrefwater    = fm.liquidModel.RHOrefwater;
  this->liquidModel.Pref           = fm.liquidModel.Pref;
  this->liquidModel.alpha          = fm.liquidModel.alpha;
  this->liquidModel.beta           = fm.liquidModel.beta;

  return *this;

}
*/
//------------------------------------------------------------------------------

ViscosityModelData::ViscosityModelData()
{

  type = SUTHERLAND;
  sutherlandReferenceTemperature = 110.6;
  sutherlandConstant = 1.458e-6;
  dynamicViscosity   = -1.0;
  bulkViscosity      = 0.0;

}

//------------------------------------------------------------------------------

void ViscosityModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<ViscosityModelData>
    (ca, "Type", this, reinterpret_cast<int ViscosityModelData::*>
                       (&ViscosityModelData::type), 3,
                       "Constant",   ViscosityModelData::CONSTANT, 
		       "Sutherland", ViscosityModelData::SUTHERLAND, 
		       "Prandtl",    ViscosityModelData::PRANDTL);
  new ClassDouble<ViscosityModelData>(ca, "SutherlandReferenceTemperature", this,
				      &ViscosityModelData::sutherlandReferenceTemperature);
  new ClassDouble<ViscosityModelData>(ca, "SutherlandConstant", this,
				      &ViscosityModelData::sutherlandConstant);
  new ClassDouble<ViscosityModelData>(ca, "DynamicViscosity", this,
				      &ViscosityModelData::dynamicViscosity);
  new ClassDouble<ViscosityModelData>(ca, "BulkViscosity", this,
				      &ViscosityModelData::bulkViscosity);

}

//------------------------------------------------------------------------------

ThermalCondModelData::ThermalCondModelData()
{

  type = CONSTANT_PRANDTL;
  prandtl = 0.72;
  conductivity = 0.0;

}

//------------------------------------------------------------------------------

void ThermalCondModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<ThermalCondModelData>
     (ca, "Type", this, reinterpret_cast<int ThermalCondModelData::*>
     (&ThermalCondModelData::type), 2,
     "ConstantPrandtl",      ThermalCondModelData::CONSTANT_PRANDTL,
     "ConstantConductivity", ThermalCondModelData::CONSTANT);

  new ClassDouble<ThermalCondModelData>
     (ca, "Prandtl", this, &ThermalCondModelData::prandtl);
  new ClassDouble<ThermalCondModelData>
     (ca, "HeatConductivity", this, &ThermalCondModelData::conductivity);

}

//------------------------------------------------------------------------------

SAModelData::SAModelData()
{

  cb1 = 0.1355;
  cb2 = 0.622;
  cw2 = 0.3;
  cw3 = 2.0;
  cv1 = 7.1;
  cv2 = 5.0;
  sigma = 2.0/3.0;
  vkcst = 0.41;
  form = ORIGINAL;

}

//------------------------------------------------------------------------------

void SAModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassDouble<SAModelData>(ca, "Cb1", this, &SAModelData::cb1);
  new ClassDouble<SAModelData>(ca, "Cb2", this, &SAModelData::cb2);
  new ClassDouble<SAModelData>(ca, "Cw2", this, &SAModelData::cw2);
  new ClassDouble<SAModelData>(ca, "Cw3", this, &SAModelData::cw3);
  new ClassDouble<SAModelData>(ca, "Cv1", this, &SAModelData::cv1);
  new ClassDouble<SAModelData>(ca, "Cv2", this, &SAModelData::cv2);
  new ClassDouble<SAModelData>(ca, "Sigma", this, &SAModelData::sigma);
  new ClassDouble<SAModelData>(ca, "Kappa", this, &SAModelData::vkcst);
  new ClassToken<SAModelData>
    (ca, "Form", this,
     reinterpret_cast<int SAModelData::*>(&SAModelData::form), 2,
     "Original", 0, "Fv3", 1);

}

//------------------------------------------------------------------------------

DESModelData::DESModelData()
{

  cb1 = 0.1355;
  cb2 = 0.622;
  cw2 = 0.3;
  cw3 = 2.0;
  cv1 = 7.1;
  cv2 = 5.0;
  cdes = 0.65;
  sigma = 2.0/3.0;
  vkcst = 0.41;
  form = ORIGINAL;

}

//------------------------------------------------------------------------------

void DESModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassDouble<DESModelData>(ca, "Cb1", this, &DESModelData::cb1);
  new ClassDouble<DESModelData>(ca, "Cb2", this, &DESModelData::cb2);
  new ClassDouble<DESModelData>(ca, "Cw2", this, &DESModelData::cw2);
  new ClassDouble<DESModelData>(ca, "Cw3", this, &DESModelData::cw3);
  new ClassDouble<DESModelData>(ca, "Cv1", this, &DESModelData::cv1);
  new ClassDouble<DESModelData>(ca, "Cv2", this, &DESModelData::cv2);
  new ClassDouble<DESModelData>(ca, "CDes", this, &DESModelData::cdes);
  new ClassDouble<DESModelData>(ca, "Sigma", this, &DESModelData::sigma);
  new ClassDouble<DESModelData>(ca, "Kappa", this, &DESModelData::vkcst);
  new ClassToken<DESModelData>
    (ca, "Form", this,
     reinterpret_cast<int DESModelData::*>(&DESModelData::form), 2,
     "Original", 0, "Fv3", 1);

}

//------------------------------------------------------------------------------

KEModelData::KEModelData()
{

  sigma_k = 1.0;
  sigma_eps = 1.0/1.4245; // 1.0/1.3
  sigma_eps1 = 1.44;
  sigma_eps2 = 11.0/6.0; // 1.92
  c_mu = 0.09;

}

//------------------------------------------------------------------------------

void KEModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassDouble<KEModelData>(ca, "SigmaK", this, &KEModelData::sigma_k);
  new ClassDouble<KEModelData>(ca, "SigmaEps", this, &KEModelData::sigma_eps);
  new ClassDouble<KEModelData>(ca, "SigmaEps1", this, &KEModelData::sigma_eps1);
  new ClassDouble<KEModelData>(ca, "SigmaEps2", this, &KEModelData::sigma_eps2);
  new ClassDouble<KEModelData>(ca, "Cmu", this, &KEModelData::c_mu);

}

//------------------------------------------------------------------------------

TurbulenceModelData::TurbulenceModelData()
{

  type = ONE_EQUATION_SPALART_ALLMARAS;

}

//------------------------------------------------------------------------------

void TurbulenceModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassToken<TurbulenceModelData>
    (ca, "Type", this, reinterpret_cast<int TurbulenceModelData::*>
     (&TurbulenceModelData::type), 3,
     "SpalartAllmaras", 0, "DES", 1, "KEpsilon", 2);

  sa.setup("SpalartAllmaras", ca);
  des.setup("DES", ca);
  ke.setup("KEpsilon", ca);

}

//------------------------------------------------------------------------------

SmagorinskyLESData::SmagorinskyLESData()
{

  c_s = 0.1;

}

//------------------------------------------------------------------------------

void SmagorinskyLESData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassDouble<SmagorinskyLESData>(ca, "Cs", this, &SmagorinskyLESData::c_s);

}

//------------------------------------------------------------------------------

WaleLESData::WaleLESData()
{

  c_w = 0.325;

}

//------------------------------------------------------------------------------

void WaleLESData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassDouble<WaleLESData>(ca, "Cw", this, &WaleLESData::c_w);

}

//------------------------------------------------------------------------------

DynamicLESData::DynamicLESData()
{

}

//------------------------------------------------------------------------------

void DynamicLESData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  clip.setup("Clipping", ca);

}

//------------------------------------------------------------------------------

ClippingData::ClippingData()
{

  cs_max = 0.4;
  pt_min = 0.05;
  pt_max = 1.6;

}

//------------------------------------------------------------------------------

void ClippingData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

 new ClassDouble<ClippingData>(ca, "CsMax", this, &ClippingData::cs_max);
 new ClassDouble<ClippingData>(ca, "PtMin", this, &ClippingData::pt_min);
 new ClassDouble<ClippingData>(ca, "PtMax", this, &ClippingData::pt_max);

}


//------------------------------------------------------------------------------

DynamicVMSData::DynamicVMSData()
{

  type = D2VMSLES;
  c_s_prime = 0.1;
  agglomeration_width = 1;
  agglomeration_depth1 = 1;
  agglomeration_depth2 = agglomeration_depth1 + 1;

}

//------------------------------------------------------------------------------

void DynamicVMSData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<DynamicVMSData>
    (ca, "Type", this, reinterpret_cast<int DynamicVMSData::*>
     (&DynamicVMSData::type), 3, "D1VMSLES", 0, "D2VMSLES", 1, "D3VMSLES", 2);
  new ClassDouble<DynamicVMSData>(ca, "Csprime", this, &DynamicVMSData::c_s_prime);
  new ClassInt<DynamicVMSData>(ca, "AgglomerationLayer", this, &DynamicVMSData::agglomeration_width);
  new ClassInt<DynamicVMSData>(ca, "AgglomerationDepth1", this, &DynamicVMSData::agglomeration_depth1);
  new ClassInt<DynamicVMSData>(ca, "AgglomerationDepth2", this, &DynamicVMSData::agglomeration_depth2);
  clip.setup("Clipping", ca);

}

//------------------------------------------------------------------------------

VMSLESData::VMSLESData()
{

  c_s_prime = 0.1;
  agglomeration_width = 1;
  agglomeration_depth = 1;

}

//------------------------------------------------------------------------------

void VMSLESData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassDouble<VMSLESData>(ca, "Csprime", this, &VMSLESData::c_s_prime);
  new ClassInt<VMSLESData>(ca, "AgglomerationLayer", this, &VMSLESData::agglomeration_width);
  new ClassInt<VMSLESData>(ca, "AgglomerationDepth", this, &VMSLESData::agglomeration_width);

}


//------------------------------------------------------------------------------

LESModelData::LESModelData()
{

  type = VMS;
  delta = VOLUME;

}

//------------------------------------------------------------------------------

void LESModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassToken<LESModelData>
    (ca, "Type", this, reinterpret_cast<int LESModelData::*>
     (&LESModelData::type), 5, "Smagorinsky", 0, "Dynamic", 1, "VMS", 2, "DynamicVMS", 3, "WALE", 4);
  new ClassToken<LESModelData>
	(ca, "Delta", this,
	 reinterpret_cast<int LESModelData::*>(&LESModelData::delta), 2,
         "Volume", 0, "Side", 1);

  sma.setup("Smagorinsky", ca);
  dles.setup("Dynamic", ca);
  vms.setup("VMS", ca);
  dvms.setup("DynamicVMS", ca);
  wale.setup("WALE", ca);

}

//------------------------------------------------------------------------------

CFixData::CFixData()
{

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  x1 = 0.0;
  y1 = 0.0;
  z1 = 0.0;

  r0 = -1.0;
  r1 = -1.0;

  failsafeN = -1;
  failsafe = ALWAYSON;
}

//------------------------------------------------------------------------------

void CFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassDouble<CFixData>(ca, "X0", this, &CFixData::x0);
  new ClassDouble<CFixData>(ca, "Y0", this, &CFixData::y0);
  new ClassDouble<CFixData>(ca, "Z0", this, &CFixData::z0);
  new ClassDouble<CFixData>(ca, "Radius0", this, &CFixData::r0);

  new ClassDouble<CFixData>(ca, "X1", this, &CFixData::x1);
  new ClassDouble<CFixData>(ca, "Y1", this, &CFixData::y1);
  new ClassDouble<CFixData>(ca, "Z1", this, &CFixData::z1);
  new ClassDouble<CFixData>(ca, "Radius1", this, &CFixData::r1);

  ClassToken<CFixData>* sf = new ClassToken<CFixData>
        (ca, "FailSafe", this,
         reinterpret_cast<int CFixData::*>(&CFixData::failsafe), 2,
         "Off", 0, "On",1, "AlwayOn", 2);

  sf->allowIntPair(&CFixData::failsafeN);
}
//------------------------------------------------------------------------------

TBFixData::TBFixData()
{

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  x1 = -1.0;
  y1 = -1.0;
  z1 = -1.0;

}

//------------------------------------------------------------------------------

void TBFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<TBFixData>(ca, "X0", this, &TBFixData::x0);
  new ClassDouble<TBFixData>(ca, "Y0", this, &TBFixData::y0);
  new ClassDouble<TBFixData>(ca, "Z0", this, &TBFixData::z0);
  new ClassDouble<TBFixData>(ca, "X1", this, &TBFixData::x1);
  new ClassDouble<TBFixData>(ca, "Y1", this, &TBFixData::y1);
  new ClassDouble<TBFixData>(ca, "Z1", this, &TBFixData::z1);


}

//------------------------------------------------------------------------------

TripDomainData::TripDomainData()
{

}

//------------------------------------------------------------------------------

void TripDomainData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  bfix.setup("Box1", ca);

}

//------------------------------------------------------------------------------

TurbulenceClosureData::TurbulenceClosureData()
{

  type = NONE;
  prandtlTurbulent = 0.9;

}

//------------------------------------------------------------------------------

void TurbulenceClosureData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<TurbulenceClosureData>
    (ca, "Type", this, reinterpret_cast<int TurbulenceClosureData::*>
     (&TurbulenceClosureData::type), 3,
     "None", 0, "TurbulenceModel", 1, "LESModel", 2);

  new ClassDouble<TurbulenceClosureData>
    (ca, "PrandtlTurbulent", this, &TurbulenceClosureData::prandtlTurbulent);

  tm.setup("TurbulenceModel", ca);
  les.setup("LESModel", ca);
  tr.setup("Tripping", ca);

}

//------------------------------------------------------------------------------

EquationsData::EquationsData()
{

  dimension = 3;
  type = EULER;
  numPhase = 1;
  gravity_x = 0.0;
  gravity_y = 0.0;
  gravity_z = 0.0;

}

//------------------------------------------------------------------------------

void EquationsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 10, father);
  new ClassInt<EquationsData>(ca, "Dimension", this, &EquationsData::dimension);

  new ClassToken<EquationsData>(ca, "Type", this,
                                reinterpret_cast<int EquationsData::*>(&EquationsData::type), 2,
                                "Euler", 0, "NavierStokes", 1);

  new ClassInt<EquationsData>(ca, "NumPhases", this,
                                &EquationsData::numPhase);

  new ClassDouble<EquationsData>(ca, "GravityX", this, &EquationsData::gravity_x);
  new ClassDouble<EquationsData>(ca, "GravityY", this, &EquationsData::gravity_y);
  new ClassDouble<EquationsData>(ca, "GravityZ", this, &EquationsData::gravity_z);

  fluidModelMap.setup("FluidModel", 0);
  //fluidModelMap.setup("FluidModel", ca); //which one to take? ca or 0?
  viscosityModel.setup("ViscosityModel", ca);
  thermalCondModel.setup("ThermalConductivityModel", ca);
  tc.setup("TurbulenceClosure", ca);

}

//------------------------------------------------------------------------------

SphereData::SphereData()
{

  fluidModelID = -1;

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  radius = -1.0;

}

//------------------------------------------------------------------------------

Assigner *SphereData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 7, nullAssigner);

  new ClassInt<SphereData> (ca, "FluidID", this, &SphereData::fluidModelID);

  new ClassDouble<SphereData> (ca, "Center_x", this, &SphereData::cen_x);
  new ClassDouble<SphereData> (ca, "Center_y", this, &SphereData::cen_y);
  new ClassDouble<SphereData> (ca, "Center_z", this, &SphereData::cen_z);
  new ClassDouble<SphereData> (ca, "Radius", this, &SphereData::radius);

  initialConditions.setup("InitialState", ca);

  programmedBurn.setup("ProgrammedBurn", ca);

  return ca;
}

//------------------------------------------------------------------------------

PrismData::PrismData()
{

  fluidModelID = -1;

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  w_x = w_y = w_z = -1.0;

}

//------------------------------------------------------------------------------

Assigner *PrismData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 9, nullAssigner);

  new ClassInt<PrismData> (ca, "FluidID", this, &PrismData::fluidModelID);

  //new ClassDouble<PrismData> (ca, "Center_x", this, &PrismData::cen_x);
  //new ClassDouble<PrismData> (ca, "Center_y", this, &PrismData::cen_y);
  //new ClassDouble<PrismData> (ca, "Center_z", this, &PrismData::cen_z);
  //new ClassDouble<PrismData> (ca, "Width_x", this, &PrismData::w_x);
  //new ClassDouble<PrismData> (ca, "Width_y", this, &PrismData::w_y);
  //new ClassDouble<PrismData> (ca, "Width_z", this, &PrismData::w_z);

  new ClassDouble<PrismData> (ca, "X0", this, &PrismData::X0);
  new ClassDouble<PrismData> (ca, "Y0", this, &PrismData::Y0);
  new ClassDouble<PrismData> (ca, "Z0", this, &PrismData::Z0);
  new ClassDouble<PrismData> (ca, "X1", this, &PrismData::X1);
  new ClassDouble<PrismData> (ca, "Y1", this, &PrismData::Y1);
  new ClassDouble<PrismData> (ca, "Z1", this, &PrismData::Z1);

  initialConditions.setup("InitialState", ca);

  programmedBurn.setup("ProgrammedBurn", ca);

  return ca;
}

bool PrismData::inside(double x,double y, double z) const {

  return (x <= cen_x+w_x*0.5 && x >=cen_x-w_x*0.5 &&
	  y <= cen_y+w_y*0.5 && y >=cen_y-w_y*0.5 &&
	  z <= cen_z+w_z*0.5 && z >=cen_z-w_z*0.5);
}

//------------------------------------------------------------------------------

PlaneData::PlaneData()
{

  fluidModelID = -1;
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 0.0;

}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);

  new ClassInt<PlaneData> (ca, "FluidID", this, &PlaneData::fluidModelID);

  new ClassDouble<PlaneData> (ca, "Point_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData> (ca, "Point_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData> (ca, "Point_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData> (ca, "Normal_x", this, &PlaneData::nx);
  new ClassDouble<PlaneData> (ca, "Normal_y", this, &PlaneData::ny);
  new ClassDouble<PlaneData> (ca, "Normal_z", this, &PlaneData::nz);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

PointData::PointData()
{
  fluidModelID = -1;
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
}

//------------------------------------------------------------------------------

//void PointData::setup(const char *name, ClassAssigner *father)
Assigner *PointData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);

  new ClassInt<PointData>
    (ca, "FluidID", this, &PointData::fluidModelID);
  new ClassDouble<PointData>
    (ca, "X", this, &PointData::x);
  new ClassDouble<PointData>
    (ca, "Y", this, &PointData::y);
  new ClassDouble<PointData>
    (ca, "Z", this, &PointData::z);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

void MultiInitialConditionsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  sphereMap.setup("Sphere", ca);
  prismMap.setup("Box", ca);
  planeMap.setup("Plane", ca);
  pointMap.setup("Point", ca);

}

//------------------------------------------------------------------------------

SparseGridData::SparseGridData()
{
  tabulationFileName = "";
  numberOfTabulations = 0;

  verbose = 0;

  minPoints = 100;
  maxPoints = 100;
  relAccuracy = 1.e-3;
  absAccuracy = 1.e-1;

  dimAdaptDegree = 0.75;

  range1min = 0.0; range1max = 1.0; mapBaseValue1 = 0.0; numDomainDim1 = 1;
  range2min = 0.0; range2max = 1.0; mapBaseValue2 = 0.0; numDomainDim2 = 1;
  range3min = 0.0; range3max = 1.0; mapBaseValue3 = 0.0; numDomainDim3 = 1;
  range4min = 0.0; range4max = 1.0; mapBaseValue4 = 0.0; numDomainDim4 = 1;
  range5min = 0.0; range5max = 1.0; mapBaseValue5 = 0.0; numDomainDim5 = 1;
  range6min = 0.0; range6max = 1.0; mapBaseValue6 = 0.0; numDomainDim6 = 1;

  range = 0;
  mapBaseValue = 0;
  numDomainDim = 0;

  numOutputs = 2;
  numInputs  = 2;

}

//------------------------------------------------------------------------------

void SparseGridData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 34, father);

  new ClassStr<SparseGridData>
    (ca, "FileName", this, &SparseGridData::tabulationFileName);
  new ClassInt<SparseGridData>
    (ca, "NumberOfFiles", this, &SparseGridData::numberOfTabulations);
  new ClassInt<SparseGridData>
    (ca, "Verbose", this, &SparseGridData::verbose);
  new ClassInt<SparseGridData>
    (ca, "MinimumNumberOfPoints", this, &SparseGridData::minPoints);
  new ClassInt<SparseGridData>
    (ca, "MaximumNumberOfPoints", this, &SparseGridData::maxPoints);
  new ClassDouble<SparseGridData>
    (ca, "RelativeAccuracy", this, &SparseGridData::relAccuracy);
  new ClassDouble<SparseGridData>
    (ca, "AbsoluteAccuracy", this, &SparseGridData::absAccuracy);
  new ClassDouble<SparseGridData>
    (ca, "Input1Minimum", this, &SparseGridData::range1min);
  new ClassDouble<SparseGridData>
    (ca, "Input1Maximum", this, &SparseGridData::range1max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase1", this, &SparseGridData::mapBaseValue1);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains1", this, &SparseGridData::numDomainDim1);
  new ClassDouble<SparseGridData>
    (ca, "Input2Minimum", this, &SparseGridData::range2min);
  new ClassDouble<SparseGridData>
    (ca, "Input2Maximum", this, &SparseGridData::range2max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase2", this, &SparseGridData::mapBaseValue2);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains2", this, &SparseGridData::numDomainDim2);
  new ClassDouble<SparseGridData>
    (ca, "Input3Minimum", this, &SparseGridData::range3min);
  new ClassDouble<SparseGridData>
    (ca, "Input3Maximum", this, &SparseGridData::range3max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase3", this, &SparseGridData::mapBaseValue3);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains3", this, &SparseGridData::numDomainDim3);
  new ClassDouble<SparseGridData>
    (ca, "Input4Minimum", this, &SparseGridData::range4min);
  new ClassDouble<SparseGridData>
    (ca, "Input4Maximum", this, &SparseGridData::range4max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase4", this, &SparseGridData::mapBaseValue4);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains4", this, &SparseGridData::numDomainDim4);
  new ClassDouble<SparseGridData>
    (ca, "Input5Minimum", this, &SparseGridData::range5min);
  new ClassDouble<SparseGridData>
    (ca, "Input5Maximum", this, &SparseGridData::range5max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase5", this, &SparseGridData::mapBaseValue5);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains5", this, &SparseGridData::numDomainDim5);
  new ClassDouble<SparseGridData>
    (ca, "Input6Minimum", this, &SparseGridData::range6min);
  new ClassDouble<SparseGridData>
    (ca, "Input6Maximum", this, &SparseGridData::range6max);
  new ClassDouble<SparseGridData>
    (ca, "LogarithmMappingBase6", this, &SparseGridData::mapBaseValue6);
  new ClassInt<SparseGridData>
    (ca, "NumberOfDomains6", this, &SparseGridData::numDomainDim6);
  new ClassInt<SparseGridData>
    (ca, "NumberOfOutputs", this, &SparseGridData::numOutputs);
  new ClassInt<SparseGridData>
    (ca, "NumberOfInputs",  this, &SparseGridData::numInputs);
  new ClassDouble<SparseGridData>
    (ca, "DegreeDimAdapt", this, &SparseGridData::dimAdaptDegree);

}

//struct ProgrammedBurnData {

//  int unburnedEOS,burnedEOS;
//  double ignitionX0,ignitionY0,ignitionZ0;
//  double cjEnergy;
//  double cjDetonationVelocity;
  
ProgrammedBurnData::ProgrammedBurnData() {

  unburnedEOS = burnedEOS = -1;
  ignitionX0 = ignitionY0 = ignitionZ0 = 0.0;
  e0 = 0.0;
  cjDetonationVelocity = -1.0;
  cjPressure = 0.0;
  cjDensity = 0.0;
  cjEnergy = 0.0;
  ignitionTime = 0.0;
  ignited = 0;
  factorB = 1.0;
  factorS = 1.0;
  limitPeak = 0;
}

ProgrammedBurnData::~ProgrammedBurnData() { 

}

void ProgrammedBurnData::setup(const char* name, ClassAssigner* father) {

  ClassAssigner *ca = new ClassAssigner(name, 15, father);

  new ClassInt<ProgrammedBurnData>
    (ca, "UnburnedEOS", this, &ProgrammedBurnData::unburnedEOS);
  new ClassInt<ProgrammedBurnData>
    (ca, "BurnedEOS", this, &ProgrammedBurnData::burnedEOS);
  new ClassDouble<ProgrammedBurnData>
    (ca, "IgnitionX0", this, &ProgrammedBurnData::ignitionX0);
  new ClassDouble<ProgrammedBurnData>
    (ca, "IgnitionY0", this, &ProgrammedBurnData::ignitionY0);
  new ClassDouble<ProgrammedBurnData>
    (ca, "IgnitionZ0", this, &ProgrammedBurnData::ignitionZ0);
  new ClassDouble<ProgrammedBurnData>
    (ca, "IgnitionTime", this, &ProgrammedBurnData::ignitionTime);
  new ClassToken<ProgrammedBurnData>
    (ca, "Ignite", this, &ProgrammedBurnData::ignited, 2, "False", 1, "True", 0);

  new ClassDouble<ProgrammedBurnData>
    (ca, "E0", this, &ProgrammedBurnData::e0);
  new ClassDouble<ProgrammedBurnData>
    (ca, "ChapmanJouguetDetonationVelocity", this, &ProgrammedBurnData::cjDetonationVelocity);
  new ClassDouble<ProgrammedBurnData>
    (ca, "ChapmanJouguetDensity", this, &ProgrammedBurnData::cjDensity);
  new ClassDouble<ProgrammedBurnData>
    (ca, "ChapmanJouguetPressure", this, &ProgrammedBurnData::cjPressure);
  new ClassDouble<ProgrammedBurnData>
    (ca, "ChapmanJouguetEnergy", this, &ProgrammedBurnData::cjEnergy);
  new ClassDouble<ProgrammedBurnData>
    (ca, "FactorB", this, &ProgrammedBurnData::factorB);
  new ClassDouble<ProgrammedBurnData>
    (ca, "FactorS", this, &ProgrammedBurnData::factorS);
  
  new ClassToken<ProgrammedBurnData>(ca, "LimitPeak", this,
             reinterpret_cast<int ProgrammedBurnData::*>(&ProgrammedBurnData::limitPeak), 2,
             "False", 0, "True", 1);
}


//------------------------------------------------------------------------------

MultiFluidData::MultiFluidData()
{

  method = GHOSTFLUID_FOR_POOR;
  problem = BUBBLE; //hidden
  typePhaseChange = RIEMANN_SOLUTION; //hidden
  riemannComputation = RK2;
  localtime  = GLOBAL; //hidden
  typeTracking = LINEAR; //hidden
  bandlevel = 5;
  subIt = 10; //hidden
  cfl = 0.7; //hidden
  frequency = 0;
  eps = 1.e-6; //hidden
  outputdiff = 0; //hidden
  copy = TRUE; //hidden

  testCase = 0; // hidden

  lsInit = VOLUMES; //hidden
  interfaceType = FSF; //hidden
  jwlRelaxationFactor = 1.0;

  interfaceTreatment = FIRSTORDER;
  interfaceExtrapolation = EXTRAPOLATIONFIRSTORDER;
  levelSetMethod = CONSERVATIVE;
  interfaceOmitCells = 0;
}

//------------------------------------------------------------------------------

void MultiFluidData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 18, father);
  new ClassToken<MultiFluidData>(ca, "Method", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::method), 3,
             "None", 0, "GhostFluidForThePoor", 1, "FiniteVolumeWithExactTwoPhaseRiemann", 2);
  new ClassToken<MultiFluidData>(ca, "Problem", this,
				 reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::problem), 2,
             "Bubble", 0, "ShockTube", 1);
  new ClassToken<MultiFluidData>(ca, "PhaseChange", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::typePhaseChange), 3,
             "None", 0, "RiemannSolution", 1, "Extrapolation", 2);
  new ClassToken<MultiFluidData>(ca, "RiemannComputation", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::riemannComputation), 4,
             "FirstOrder", 0, "SecondOrder", 1, "TabulationRiemannInvariant", 2,
             "TabulationRiemannProblem", 3);
  new ClassToken<MultiFluidData>(ca, "FictitiousTimeStepping", this,
		         reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::localtime),2,
             "Global", 0, "Local", 1);
  new ClassToken<MultiFluidData>(ca, "InterfaceTracking", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::typeTracking),3,
             "Linear", 0, "Gradient", 1, "Hermite", 2);
  new ClassInt<MultiFluidData>(ca, "BandLevel", this,
             &MultiFluidData::bandlevel);
  new ClassInt<MultiFluidData>(ca, "SubIt", this,
             &MultiFluidData::subIt);
  new ClassDouble<MultiFluidData>(ca, "Cfl", this,
             &MultiFluidData::cfl);
  new ClassInt<MultiFluidData>(ca, "LevelSetReinitializationFrequency", this,
             &MultiFluidData::frequency);
  new ClassDouble<MultiFluidData>(ca, "Epsilon", this,
             &MultiFluidData::eps);
  new ClassInt<MultiFluidData>(ca, "OutputDiff", this,
             &MultiFluidData::outputdiff);
  new ClassToken<MultiFluidData>(ca, "CopyCloseNodes", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::copy),2,
             "False", 0, "True", 1);
  new ClassToken<MultiFluidData>(ca, "LSInit", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::lsInit),3,
             "Old", 0, "Volumes", 1, "Geometric", 2);
  new ClassToken<MultiFluidData>(ca, "InterfaceType", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::interfaceType),3,
             "FluidStructureFluid", 0, "FluidFluid", 1, "BOTH", 2);

  new ClassToken<MultiFluidData>(ca, "InterfaceTreatment", this,
				 reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::interfaceTreatment),2,
				 "FirstOrder", 0, "SecondOrder", 1);
  new ClassToken<MultiFluidData>(ca, "InterfaceExtrapolation", this,
				 reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::interfaceExtrapolation),2,
				 "FirstOrder", 0, "SecondOrder", 1);

  new ClassToken<MultiFluidData>(ca, "LevelSetMethod", this,
				 reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::levelSetMethod),4,
				 "Conservative", 0, "HJWENO", 1,"Scalar", 2, "Primitive",3);

  new ClassDouble<MultiFluidData>(ca, "JwlRelaxationFactor", this,
				  &MultiFluidData::jwlRelaxationFactor);

  new ClassInt<MultiFluidData>(ca, "TestCase", this,
			       &MultiFluidData::testCase);
  
  new ClassInt<MultiFluidData>(ca, "OmitCells", this,
			       &MultiFluidData::interfaceOmitCells);

  multiInitialConditions.setup("InitialConditions", ca);
  sparseGrid.setup("SparseGrid",ca);

}

//------------------------------------------------------------------------------

SchemeData::SchemeData(int af) : allowsFlux(af)
{

  advectiveOperator = FINITE_VOLUME;
  flux = ROE;
  reconstruction = LINEAR;
  limiter = NONE;
  gradient = LEAST_SQUARES;
  dissipation = SECOND_ORDER;

  beta = 1.0/3.0;
  gamma = 1.0;
  xiu  = -2.0/15.0;
  xic = -1.0/30.0;
  eps = 0.1;

}

//------------------------------------------------------------------------------

void SchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner* ca;
  if (allowsFlux)
    ca = new ClassAssigner(name, 11, father);
  else
    ca = new ClassAssigner(name, 10, father);

  new ClassToken<SchemeData>
    (ca, "AdvectiveOperator", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::advectiveOperator), 2,
     "FiniteVolume", 0, "Galerkin", 1);

  if (allowsFlux) {
    new ClassToken<SchemeData>
      (ca, "Flux", this,
       reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 4,
       "Roe", 0, "VanLeer", 1, "HLLE", 2, "HLLC", 3);
  }

  new ClassToken<SchemeData>
    (ca, "Reconstruction", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::reconstruction), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<SchemeData>
    (ca, "Limiter", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::limiter), 5,
     "None", 0, "VanAlbada", 1, "Barth", 2, "Venkatakrishnan", 3, "PressureSensor", 4);

  new ClassToken<SchemeData>
    (ca, "Gradient", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::gradient), 3,
     "LeastSquares", 0, "Galerkin", 1, "NonNodal", 2);

  new ClassToken<SchemeData>
    (ca, "Dissipation", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::dissipation), 2,
     "SecondOrder", 0, "SixthOrder", 1);

  new ClassDouble<SchemeData>(ca, "Beta", this, &SchemeData::beta);
  new ClassDouble<SchemeData>(ca, "Gamma", this, &SchemeData::gamma);
  new ClassDouble<SchemeData>(ca, "XiU", this, &SchemeData::xiu);
  new ClassDouble<SchemeData>(ca, "XiC", this, &SchemeData::xic);
  new ClassDouble<SchemeData>(ca, "Eps", this, &SchemeData::eps);

}

//------------------------------------------------------------------------------

SFixData::SFixData()
{

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  r = -1.0;

  failsafe = ALWAYSON;
  failsafeN = -1;
}
//------------------------------------------------------------------------------

void SFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassDouble<SFixData>(ca, "X0", this, &SFixData::x0);
  new ClassDouble<SFixData>(ca, "Y0", this, &SFixData::y0);
  new ClassDouble<SFixData>(ca, "Z0", this, &SFixData::z0);
  new ClassDouble<SFixData>(ca, "Radius", this, &SFixData::r);

  ClassToken<SFixData>* sf = new ClassToken<SFixData>
	(ca, "FailSafe", this,
	 reinterpret_cast<int SFixData::*>(&SFixData::failsafe), 2,
         "Off", 0, "On",1, "AlwayOn", 2);
  
  sf->allowIntPair(&SFixData::failsafeN);
}

//------------------------------------------------------------------------------

BFixData::BFixData()
{

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  x1 = -1.0;
  y1 = -1.0;
  z1 = -1.0;

  failsafe = ALWAYSON;
  failsafeN = -1;
}

//------------------------------------------------------------------------------

void BFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<BFixData>(ca, "X0", this, &BFixData::x0);
  new ClassDouble<BFixData>(ca, "Y0", this, &BFixData::y0);
  new ClassDouble<BFixData>(ca, "Z0", this, &BFixData::z0);
  new ClassDouble<BFixData>(ca, "X1", this, &BFixData::x1);
  new ClassDouble<BFixData>(ca, "Y1", this, &BFixData::y1);
  new ClassDouble<BFixData>(ca, "Z1", this, &BFixData::z1);

  ClassToken<BFixData>* sf = new ClassToken<BFixData>
	(ca, "FailSafe", this,
	 reinterpret_cast<int BFixData::*>(&BFixData::failsafe), 2,
         "Off", 0, "On",1, "AlwayOn", 2);
  
  sf->allowIntPair(&BFixData::failsafeN);

}

//------------------------------------------------------------------------------

SchemeFixData::SchemeFixData()
{

  symmetry = NONE;
  dihedralAngle = -1.0;

  spheres[0] = &sfix1;
  spheres[1] = &sfix2;
  spheres[2] = &sfix3;
  spheres[3] = &sfix4;
  spheres[4] = &sfix5;
  spheres[5] = &sfix6;
  spheres[6] = &sfix7;
  spheres[7] = &sfix8;
  spheres[8] = &sfix9;
  spheres[9] = &sfix10;

  boxes[0] = &bfix1;
  boxes[1] = &bfix2;
  boxes[2] = &bfix3;
  boxes[3] = &bfix4;
  boxes[4] = &bfix5;
  boxes[5] = &bfix6;
  boxes[6] = &bfix7;
  boxes[7] = &bfix8;
  boxes[8] = &bfix9;
  boxes[9] = &bfix10;

  cones[0] = &cfix1;
  cones[1] = &cfix2;
  cones[2] = &cfix3;
  cones[3] = &cfix4;
  cones[4] = &cfix5;
  cones[5] = &cfix6;
  cones[6] = &cfix7;
  cones[7] = &cfix8;
  cones[8] = &cfix9;
  cones[9] = &cfix10;


}

//------------------------------------------------------------------------------

void SchemeFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 32, father);

  new ClassToken<SchemeFixData>
    (ca, "Symmetry", this,
     reinterpret_cast<int SchemeFixData::*>(&SchemeFixData::symmetry), 4,
     "None", 0, "X", 1, "Y", 2, "Z", 3);
  new ClassDouble<SchemeFixData>(ca, "DihedralAngle", this, &SchemeFixData::dihedralAngle);

  sfix1.setup("Sphere1", ca);
  sfix2.setup("Sphere2", ca);
  sfix3.setup("Sphere3", ca);
  sfix4.setup("Sphere4", ca);
  sfix5.setup("Sphere5", ca);
  sfix6.setup("Sphere6", ca);
  sfix7.setup("Sphere7", ca);
  sfix8.setup("Sphere8", ca);
  sfix9.setup("Sphere9", ca);
  sfix10.setup("Sphere10", ca);

  bfix1.setup("Box1", ca);
  bfix2.setup("Box2", ca);
  bfix3.setup("Box3", ca);
  bfix4.setup("Box4", ca);
  bfix5.setup("Box5", ca);
  bfix6.setup("Box6", ca);
  bfix7.setup("Box7", ca);
  bfix8.setup("Box8", ca);
  bfix9.setup("Box9", ca);
  bfix10.setup("Box10", ca);

  cfix1.setup("Cone1", ca);
  cfix2.setup("Cone2", ca);
  cfix3.setup("Cone3", ca);
  cfix4.setup("Cone4", ca);
  cfix5.setup("Cone5", ca);
  cfix6.setup("Cone6", ca);
  cfix7.setup("Cone7", ca);
  cfix8.setup("Cone8", ca);
  cfix9.setup("Cone9", ca);
  cfix10.setup("Cone10", ca);
}

//------------------------------------------------------------------------------

BoundarySchemeData::BoundarySchemeData()
{

  type = STEGER_WARMING;

}

//------------------------------------------------------------------------------

void BoundarySchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<BoundarySchemeData>(ca, "Type", this,
         reinterpret_cast<int BoundarySchemeData::*>(&BoundarySchemeData::type), 5,
         "StegerWarming", 0, "ConstantExtrapolation", 1, "LinearExtrapolation", 2,
	 "Ghidaglia", 3, "ModifiedGhidaglia", 4);

}

//------------------------------------------------------------------------------

SchemesData::SchemesData() : ls(0), tm(0)
{

}

//------------------------------------------------------------------------------

void SchemesData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  ns.setup("NavierStokes", ca);
  tm.setup("TurbulenceModel", ca);
  ls.setup("LevelSet",ca);
  fixes.setup("Fixes", ca);
  bc.setup("Boundaries", ca);

  tm.reconstruction = SchemeData::CONSTANT;
}

//------------------------------------------------------------------------------

ExplicitData::ExplicitData()
{

  type = RUNGE_KUTTA_4;

}

//------------------------------------------------------------------------------

void ExplicitData::setup(const char *name, ClassAssigner *father)
{

 ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ExplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 5,
     "RungeKutta4", 0, "RungeKutta2", 1, "ForwardEuler", 2, 
     "OneBlockRK2", 3, "OneBlockRK2bis", 4);

}

//------------------------------------------------------------------------------

PcData::PcData()
{

  type = RAS;
  renumbering = RCM;

  fill = 0;

}

//------------------------------------------------------------------------------

void PcData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<PcData>(ca, "Type", this,
			 reinterpret_cast<int PcData::*>(&PcData::type), 7,
			 "Identity", 0, "Jacobi", 1, "As", 2, "Ras", 3, "Has", 4, "Aas", 5, "MultiGrid", 6);

  new ClassToken<PcData>(ca, "Renumbering", this,
			 reinterpret_cast<int PcData::*>(&PcData::renumbering), 2,
			 "Natural", 0, "Rcm", 1);

  new ClassInt<PcData>(ca, "Fill", this, &PcData::fill);

}

//------------------------------------------------------------------------------

KspData::KspData()
{

  type = GMRES;
  epsFormula = CONSTANT;
  checkFinalRes = NO;

  maxIts = 30;
  numVectors = 30;
  eps = 1.e-2;

  output = "";

}

//------------------------------------------------------------------------------

void KspData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<KspData>(ca, "Type", this,
			   reinterpret_cast<int KspData::*>(&KspData::type), 4,
			   "Richardson", 0, "Cg", 1, "Gmres", 2, "Gcr", 3);

  new ClassToken<KspData>(ca, "EpsFormula", this,
			   reinterpret_cast<int KspData::*>(&KspData::epsFormula), 2,
			   "Constant", 0, "Eisenstadt", 1);

  new ClassToken<KspData>(ca, "CheckFinalRes", this,
			   reinterpret_cast<int KspData::*>(&KspData::checkFinalRes), 2,
			   "No", 0, "Yes", 1);

  new ClassInt<KspData>(ca, "MaxIts", this, &KspData::maxIts);

  new ClassInt<KspData>(ca, "KrylovVectors", this, &KspData::numVectors);

  new ClassDouble<KspData>(ca, "Eps", this, &KspData::eps);

  new ClassStr<KspData>(ca, "Output", this, &KspData::output);

  pc.setup("Preconditioner", ca);

}

//------------------------------------------------------------------------------

KspFluidData::KspFluidData()
{

}

//------------------------------------------------------------------------------

void KspFluidData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  ns.setup("NavierStokes", ca);
  tm.setup("TurbulenceModel", ca);
  lsi.setup("LevelSet", ca);

}

//------------------------------------------------------------------------------

ImplicitData::ImplicitData()
{

  type = BACKWARD_EULER;
  startup = REGULAR;
  tmcoupling = WEAK;
  mvp = H1;
  fdOrder = FIRST_ORDER;
  fvmers_3pbdf = BDF_SCHEME2;
  //normals = AUTO;
  //velocities = AUTO_VEL;
 
  // (Slave) Flag for the Jacobian of the flux function
  ffjacobian = APPROXIMATE;

}

//------------------------------------------------------------------------------

void ImplicitData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<ImplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::type), 4,
     "BackwardEuler", 0, "CrankNicolson", 1, "ThreePointBackwardDifference", 2,
     "FourPointBackwardDifference", 3);

  new ClassToken<ImplicitData>
    (ca, "Startup", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::startup), 2,
     "Regular", 0, "Modified", 1);

  new ClassToken<ImplicitData>
    (ca, "TurbulenceModelCoupling", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::tmcoupling), 2,
     "Weak", 0, "Strong", 1);

  new ClassToken<ImplicitData>
    (ca, "MatrixVectorProduct", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::mvp), 4,
     "FiniteDifference", 0, "Approximate", 1, "Exact", 2,
     "ApproximateFiniteDifference", 3);

  new ClassToken<ImplicitData>
    (ca, "FiniteDifferenceOrder", this,
      reinterpret_cast<int ImplicitData::*>(&ImplicitData::fdOrder), 2,
      "FirstOrder", FIRST_ORDER, "SecondOrder", SECOND_ORDER);


  new ClassToken<ImplicitData>
    (ca, "FVMERSBDFScheme", this,
      reinterpret_cast<int ImplicitData::*>(&ImplicitData::fvmers_3pbdf), 2,
      "Scheme1", BDF_SCHEME1, "Scheme2", BDF_SCHEME2);
 

  /*new ClassToken<ImplicitData>
    (ca, "Normals", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::normals), 7,
     "FirstOrderGcl", 1, "SecondOrderGcl", 2, "FirstOrderEZGcl", 3,
     "SecondOrderEZGcl", 4, "ThirdOrderEZGcl", 5, "CurrentConfiguration", 6,
     "LatestConfiguration", 7);

  new ClassToken<ImplicitData>
    (ca, "Velocities", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::velocities), 5,
     "BackwardEuler", 1, "ThreePointBackwardDifference", 2, "Imposed", 3,
     "ImposedBackwardEuler", 4, "ImposedThreePointBackwardDifference", 5);
  */
  newton.setup("Newton", ca);

}

//------------------------------------------------------------------------------

TsData::TsData()
{

  type = IMPLICIT;
  typeTimeStep = AUTO;
  typeClipping = FREESTREAM;
  timeStepCalculation = CFL;
  dualtimestepping = OFF;

  prec = NO_PREC;
  viscousCst = 0.0;

  maxIts = 100;
  eps = 1.e-6;
  timestep = -1.0;
  timestepinitial = -1.0;
  maxTime = 1.e99;

  residual = -1;
  cfl0 = 5.0;
  cflCoef1 = 0.0;
  cflCoef2 = 0.0;
  cflMax = 1000.0;
  cflMin = 1.0;
  ser = 0.7;
  errorTol = 1.e-10;
  form = NONDESCRIPTOR;
  dualtimecfl = 100.0;

  output = "";

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 23, father);

  new ClassToken<TsData>(ca, "Type", this,
			 reinterpret_cast<int TsData::*>(&TsData::type), 2,
			 "Explicit", 0, "Implicit", 1);
  new ClassToken<TsData>(ca, "TypeTimeStep", this,
			 reinterpret_cast<int TsData::*>(&TsData::typeTimeStep), 2,
			 "Local", 1, "Global", 2);
  new ClassToken<TsData>(ca, "DualTimeStepping", this,
                         reinterpret_cast<int TsData::*>(&TsData::dualtimestepping), 2,
                         "Off", 0, "On", 1);
  new ClassToken<TsData>(ca, "Clipping", this,
			 reinterpret_cast<int TsData::*>(&TsData::typeClipping), 3,
			 "None", 0, "AbsoluteValue", 1, "Freestream", 2);
  new ClassToken<TsData>(ca, "TimeStepAdaptation", this,
			 reinterpret_cast<int TsData::*>(&TsData::timeStepCalculation), 2,
			 "Cfl", 0, "ErrorEstimation", 1);
  new ClassToken<TsData>(ca, "Prec", this,
                         reinterpret_cast<int TsData::*>(&TsData::prec), 2,
                         "NonPreconditioned", 0, "LowMach", 1);
  new ClassDouble<TsData>(ca, "ViscosityParameter", this, &TsData::viscousCst);

  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "Eps", this, &TsData::eps);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "TimeStepInitial", this, &TsData::timestepinitial);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);
  new ClassInt<TsData>(ca, "Residual", this, &TsData::residual);
  new ClassDouble<TsData>(ca, "Cfl0", this, &TsData::cfl0);
  new ClassDouble<TsData>(ca, "Cfl1", this, &TsData::cflCoef1);
  new ClassDouble<TsData>(ca, "Cfl2", this, &TsData::cflCoef2);
  new ClassDouble<TsData>(ca, "CflMax", this, &TsData::cflMax);
  new ClassDouble<TsData>(ca, "CflMin", this, &TsData::cflMin);
  new ClassDouble<TsData>(ca, "Ser", this, &TsData::ser);
  new ClassDouble<TsData>(ca, "ErrorTol", this, &TsData::errorTol);
  new ClassDouble<TsData>(ca, "DualTimeCfl", this, &TsData::dualtimecfl);
  new ClassStr<TsData>(ca, "Output", this, &TsData::output);
  new ClassToken<TsData> (ca, "Form", this, reinterpret_cast<int TsData::*>(&TsData::form), 3, "NonDescriptor", 0, "Descriptor", 1, "Hybrid", 2);  

  expl.setup("Explicit", ca);
  implicit.setup("Implicit", ca);

}

//------------------------------------------------------------------------------

DGCLData::DGCLData()
{

  normals = AUTO;
  velocities = AUTO_VEL;

}

//------------------------------------------------------------------------------

void DGCLData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassToken<DGCLData>
    (ca, "Normals", this,
     reinterpret_cast<int DGCLData::*>(&DGCLData::normals), 8,
     "FirstOrderGcl", 1, "SecondOrderGcl", 2, "FirstOrderEZGcl", 3,
     "SecondOrderEZGcl", 4, "ThirdOrderEZGcl", 5, "CurrentConfiguration", 6,
     "LatestConfiguration", 7, "RK2Gcl", 8);

  new ClassToken<DGCLData>
    (ca, "Velocities", this,
     reinterpret_cast<int DGCLData::*>(&DGCLData::velocities), 6,
     "BackwardEuler", 1, "ThreePointBackwardDifference", 2, "Imposed", 3,
     "ImposedBackwardEuler", 4, "ImposedThreePointBackwardDifference", 5,
     "RK2Gcl", 7);

}

//------------------------------------------------------------------------------


// Included (MB)
SensitivityAnalysis::SensitivityAnalysis()
{
  method  = DIRECT;
  scFlag = ANALYTICAL;
  eps = 0.00001;
  sensMesh = OFF_SENSITIVITYMESH;
  sensMach = OFF_SENSITIVITYMACH;
  sensAlpha = OFF_SENSITIVITYALPHA;
  sensBeta = OFF_SENSITIVITYBETA;
  si = 0;
  sf = -1;

  // For debugging purposes
  excsol = OFF_EXACTSOLUTION;
  homotopy = OFF_HOMOTOPY;
  comp3d = ON_COMPATIBLE3D;
  angleRad = OFF_ANGLERAD;
  machref = -1.0;
  alpharef = 400.0;
  betaref = 400.0;
  sensoutput = "";
  fres = 0.0;
  fixsol = NONEFIX;
  avgsIt = 0;

}

//------------------------------------------------------------------------------

// Included (MB)
void SensitivityAnalysis::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 21, father);
  new ClassToken<SensitivityAnalysis>(ca, "Method", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::method), 2, "Direct", 0, "Adjoint", 1);
  new ClassToken<SensitivityAnalysis>(ca, "SensitivityComputation", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::scFlag), 3, "Analytical", 0, "SemiAnalytical", 1, "FiniteDifference", 2);
  new ClassDouble<SensitivityAnalysis>(ca, "EpsFD", this, &SensitivityAnalysis::eps);
  new ClassToken<SensitivityAnalysis>(ca, "SensitivityMesh", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::sensMesh), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "SensitivityMach", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::sensMach), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "SensitivityAlpha", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::sensAlpha), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "SensitivityBeta", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::sensBeta), 2, "Off", 0, "On", 1);
  new ClassInt<SensitivityAnalysis>(ca, "ShapeVariableInitial", this, &SensitivityAnalysis::si);
  new ClassInt<SensitivityAnalysis>(ca, "ShapeVariableFinal", this, &SensitivityAnalysis::sf);

// For debugging purposes
  new ClassToken<SensitivityAnalysis>(ca, "ExactSolution", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::excsol), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "HomotopyComputation", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::homotopy), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "Compatible3D", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::comp3d), 2, "Off", 0, "On", 1);
  new ClassToken<SensitivityAnalysis>(ca, "AngleRadians", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::angleRad), 2, "Off", 0, "On", 1);

  new ClassDouble<SensitivityAnalysis>(ca, "MachReference", this, &SensitivityAnalysis::machref);
  new ClassDouble<SensitivityAnalysis>(ca, "AlphaReference", this, &SensitivityAnalysis::alpharef);
  new ClassDouble<SensitivityAnalysis>(ca, "BetaReference", this, &SensitivityAnalysis::betaref);
  new ClassStr<SensitivityAnalysis>(ca, "SensitivityOutput", this, &SensitivityAnalysis::sensoutput);
  new ClassDouble<SensitivityAnalysis>(ca, "ForceResidual", this, &SensitivityAnalysis::fres);
  new ClassToken<SensitivityAnalysis>(ca, "FixSolution", this, reinterpret_cast<int SensitivityAnalysis::*>(&SensitivityAnalysis::fixsol), 2, "None", 0, "PreviousValues", 1);
  new ClassInt<SensitivityAnalysis>(ca, "AverageStateIterations", this, &SensitivityAnalysis::avgsIt);

  ksp.setup("LinearSolver", ca);

}

//------------------------------------------------------------------------------

SymmetryData::SymmetryData()
{
  nx = 0.0;
  ny = 0.0;
  nz = 0.0;
}

//------------------------------------------------------------------------------

void SymmetryData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassDouble<SymmetryData>(ca, "Nx", this, &SymmetryData::nx);
  new ClassDouble<SymmetryData>(ca, "Ny", this, &SymmetryData::ny);
  new ClassDouble<SymmetryData>(ca, "Nz", this, &SymmetryData::nz);

  //HB: ToDo: check if it is a cannonical plane ...
  //double nrm = sqrt(nx*nx+ny*ny+nz*nz);
  //if(nrm==0.0)
  //if( ((nx==1.0) & (ny!=0.0) & (nz!=0.0)) || ((nx!=0.0) & (ny==1.0) & (nz!=0.0)) || ((nx!=0.0) & (ny!=0.0) & (nz==1.0)) )
}

//------------------------------------------------------------------------------

DefoMeshMotionData::DefoMeshMotionData()
{

  type = BASIC;
  element = BALL_VERTEX;
  volStiff = 0.0;

  mode = NonRecursive;
  numIncrements = 1;

  newton.ksp.type = KspData::CG;
  newton.ksp.epsFormula = KspData::CONSTANT;
  newton.ksp.checkFinalRes = KspData::NO;
  newton.ksp.maxIts = 20;
  newton.ksp.eps = 1.e-3;
  newton.ksp.pc.type = PcData::JACOBI;

}

//------------------------------------------------------------------------------

void DefoMeshMotionData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<DefoMeshMotionData>
    (ca, "Type", this,
     reinterpret_cast<int DefoMeshMotionData::*>(&DefoMeshMotionData::type), 2,
     "Basic", 0, "Corotational", 1);

  new ClassToken<DefoMeshMotionData>
    (ca, "Element", this,
     reinterpret_cast<int DefoMeshMotionData::*>(&DefoMeshMotionData::element), 5,
     "LinearFiniteElement", 0, "NonLinearFiniteElement", 1, "TorsionalSprings", 2, "BallVertexSprings", 3,
         "NonLinearBallVertex", 4);

  new ClassDouble<DefoMeshMotionData>(ca, "VolumeStiffness", this, &DefoMeshMotionData::volStiff);
  new ClassToken<DefoMeshMotionData>
    (ca, "Mode", this,
     reinterpret_cast<int DefoMeshMotionData::*>(&DefoMeshMotionData::mode), 2,
     "Recursive", 1, "NonRecursive", 2);
  new ClassInt<DefoMeshMotionData>(ca, "NumIncrements", this, &DefoMeshMotionData::numIncrements);

  symmetry.setup("Symmetry", ca);
  newton.setup("Newton", ca);
  blmeshmotion.setup("BoundaryLayer" , ca);

}

//------------------------------------------------------------------------------

BLMeshMotionData::BLMeshMotionData()
{

  type = PSEUDOSTRUCTURAL;
  numLayers = 0;
  power = 4.0;
  neiSelectionDistFactor = 1.0;
  fractionalStrategy  = Distance; // Based on Dist
  bestNeiStrategy = Fixed; // is On
  feedbackFrequency = 0; // 0 will turn off the feedback
  numIncrements = 1;
}

//------------------------------------------------------------------------------
void BLMeshMotionData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 8, father);

  new ClassToken<BLMeshMotionData>(ca, "Type", this,reinterpret_cast<int BLMeshMotionData::*>(&BLMeshMotionData::type), 2 ,"PseudoStructural", 0, "Algebraic", 1);
  new ClassInt<BLMeshMotionData>(ca, "NumLayers", this, &BLMeshMotionData::numLayers);
  new ClassDouble<BLMeshMotionData>(ca, "Power", this, &BLMeshMotionData::power);
  new ClassDouble<BLMeshMotionData>(ca, "DistanceFactor", this, &BLMeshMotionData::neiSelectionDistFactor);

  new ClassToken<BLMeshMotionData>(ca, "FractionalStrategy", this,reinterpret_cast<int BLMeshMotionData::*>(&BLMeshMotionData::fractionalStrategy), 2 ,"Distance", 1, "DotProduct", 2);
  new ClassToken<BLMeshMotionData>(ca, "ClosestNeighbor", this,reinterpret_cast<int BLMeshMotionData::*>(&BLMeshMotionData::bestNeiStrategy), 2 ,"Fixed", 1, "Variable", 0 );
  new ClassInt<BLMeshMotionData>(ca, "FeedBack", this,&BLMeshMotionData::feedbackFrequency);
  new ClassInt<BLMeshMotionData>(ca, "NumIncrements", this, &BLMeshMotionData::numIncrements);

}
//------------------------------------------------------------------------------

VelocityPoints::VelocityPoints()
{

  time      = -1.0;
  velocityX =  0.0;
  velocityY =  0.0;
  velocityZ =  0.0;

}

//------------------------------------------------------------------------------

void VelocityPoints::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassDouble<VelocityPoints>(ca, "Time", this, &VelocityPoints::time);
  new ClassDouble<VelocityPoints>(ca, "VelocityX", this, &VelocityPoints::velocityX);
  new ClassDouble<VelocityPoints>(ca, "VelocityY", this, &VelocityPoints::velocityY);
  new ClassDouble<VelocityPoints>(ca, "VelocityZ", this, &VelocityPoints::velocityZ);
}

//------------------------------------------------------------------------------

ForcePoints::ForcePoints()
{

  time  = -1.0;
  force =  0.0;

}

//------------------------------------------------------------------------------

void ForcePoints::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassDouble<ForcePoints>(ca, "Time", this, &ForcePoints::time);
  new ClassDouble<ForcePoints>(ca, "Force", this, &ForcePoints::force);

}

//------------------------------------------------------------------------------

RigidMeshMotionData::RigidMeshMotionData()
{

  tag = MACH;
  lawtype = CONSTANTACCELERATION;

  vx = 0.0;
  vy = 0.0;
  vz = 0.0;

  ax = 0.0;
  ay = 0.0;
  az = 0.0;

  vpts[0] = &vpts1;
  vpts[1] = &vpts2;
  vpts[2] = &vpts3;
  vpts[3] = &vpts4;
  vpts[4] = &vpts5;
  vpts[5] = &vpts6;
  vpts[6] = &vpts7;
  vpts[7] = &vpts8;
  vpts[8] = &vpts9;
  vpts[9] = &vpts10;

  timestep = -1.0;

}

//------------------------------------------------------------------------------

void RigidMeshMotionData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 19, father);

  new ClassToken<RigidMeshMotionData>
    (ca, "Tag", this,
     reinterpret_cast<int RigidMeshMotionData::*>(&RigidMeshMotionData::tag), 3,
     "Mach", 0, "Time", 1, "Velocity", 2);
  new ClassToken<RigidMeshMotionData>
    (ca, "LawType", this,
     reinterpret_cast<int RigidMeshMotionData::*>(&RigidMeshMotionData::lawtype), 2,
     "VelocityLaw", 0, "ConstantAcceleration", 1);

  new ClassDouble<RigidMeshMotionData>(ca, "VelocityX", this, &RigidMeshMotionData::vx);
  new ClassDouble<RigidMeshMotionData>(ca, "VelocityY", this, &RigidMeshMotionData::vy);
  new ClassDouble<RigidMeshMotionData>(ca, "VelocityZ", this, &RigidMeshMotionData::vz);

  new ClassDouble<RigidMeshMotionData>(ca, "AccelerationX", this, &RigidMeshMotionData::ax);
  new ClassDouble<RigidMeshMotionData>(ca, "AccelerationY", this, &RigidMeshMotionData::ay);
  new ClassDouble<RigidMeshMotionData>(ca, "AccelerationZ", this, &RigidMeshMotionData::az);

  vpts1.setup("TimeVelocity1", ca);
  vpts2.setup("TimeVelocity2", ca);
  vpts3.setup("TimeVelocity3", ca);
  vpts4.setup("TimeVelocity4", ca);
  vpts5.setup("TimeVelocity5", ca);
  vpts6.setup("TimeVelocity6", ca);
  vpts7.setup("TimeVelocity7", ca);
  vpts8.setup("TimeVelocity8", ca);
  vpts9.setup("TimeVelocity9", ca);
  vpts10.setup("TimeVelocity10", ca);

  new ClassDouble<RigidMeshMotionData>(ca, "TimeStep", this, &RigidMeshMotionData::timestep);

}

//------------------------------------------------------------------------------

AeroelasticData::AeroelasticData()
{

  force = LAST;
  pressure = -1.0;
  displacementScaling = 1.0;
  forceScaling = 1.0;
  powerScaling = 1.0;

}

//------------------------------------------------------------------------------

void AeroelasticData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<AeroelasticData>
    (ca, "Force", this,
     reinterpret_cast<int AeroelasticData::*>(&AeroelasticData::force), 3,
     "Last", 0, "Averaged", 1, "LastKris", 2);

  new ClassDouble<AeroelasticData>(ca, "InternalPressure", this, &AeroelasticData::pressure);
  new ClassDouble<AeroelasticData>(ca, "DisplacementScaling", this, &AeroelasticData::displacementScaling);
  new ClassDouble<AeroelasticData>(ca, "ForceScaling", this, &AeroelasticData::forceScaling);
  new ClassDouble<AeroelasticData>(ca, "PowerScaling", this, &AeroelasticData::powerScaling);

}

//------------------------------------------------------------------------------

ForcedData::ForcedData()
{

  type = HEAVING;
  frequency = -1.0;
  timestep = -1.0;

}

//------------------------------------------------------------------------------

void ForcedData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<ForcedData>
    (ca, "Type", this,
     reinterpret_cast<int ForcedData::*>(&ForcedData::type), 3,
     "Heaving", 0, "Pitching", 1, "Deforming", 2);

  new ClassDouble<ForcedData>(ca, "Frequency", this, &ForcedData::frequency);
  new ClassDouble<ForcedData>(ca, "TimeStep", this, &ForcedData::timestep);

  hv.setup("Heaving", ca);
  pt.setup("Pitching", ca);
  df.setup("Deforming", ca);

}

//------------------------------------------------------------------------------

HeavingData::HeavingData()
{

  domain = VOLUME;
  ax = 0.0;
  ay = 0.0;
  az = 0.0;

}

//------------------------------------------------------------------------------

void HeavingData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassToken<HeavingData>
    (ca, "Domain", this,
     reinterpret_cast<int HeavingData::*>(&HeavingData::domain), 2,
     "Volume", 0, "Surface", 1);

  new ClassDouble<HeavingData>(ca, "AX", this, &HeavingData::ax);
  new ClassDouble<HeavingData>(ca, "AY", this, &HeavingData::ay);
  new ClassDouble<HeavingData>(ca, "AZ", this, &HeavingData::az);

}

//------------------------------------------------------------------------------

PitchingData::PitchingData()
{

  domain = VOLUME;
  alpha_in = 0.0;
  alpha_max = 0.0;
  x1 =  0.0;
  y1 = -1.0;
  z1 =  0.0;
  x2 =  0.0;
  y2 =  1.0;
  z2 =  0.0;

}

//------------------------------------------------------------------------------

void PitchingData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassToken<PitchingData>
    (ca, "Domain", this,
     reinterpret_cast<int PitchingData::*>(&PitchingData::domain), 2,
     "Volume", 0, "Surface", 1);

  new ClassDouble<PitchingData>(ca, "Alpha0", this, &PitchingData::alpha_in);
  new ClassDouble<PitchingData>(ca, "AlphaMax", this, &PitchingData::alpha_max);
  new ClassDouble<PitchingData>(ca, "X1", this, &PitchingData::x1);
  new ClassDouble<PitchingData>(ca, "Y1", this, &PitchingData::y1);
  new ClassDouble<PitchingData>(ca, "Z1", this, &PitchingData::z1);
  new ClassDouble<PitchingData>(ca, "X2", this, &PitchingData::x2);
  new ClassDouble<PitchingData>(ca, "Y2", this, &PitchingData::y2);
  new ClassDouble<PitchingData>(ca, "Z2", this, &PitchingData::z2);

}

//------------------------------------------------------------------------------

DeformingData::DeformingData()
{

  domain = VOLUME;
  positions = "";
  amplification = 1.0;

}

//------------------------------------------------------------------------------

void DeformingData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<DeformingData>
    (ca, "Domain", this,
     reinterpret_cast<int DeformingData::*>(&DeformingData::domain), 2,
     "Volume", 0, "Surface", 1);

  new ClassStr<DeformingData>(ca, "Position", this, &DeformingData::positions);
  new ClassDouble<DeformingData>(ca, "Amplification", this, &DeformingData::amplification);

}

ModelReductionData::ModelReductionData()
{
	projection = PETROV_GALERKIN;
	systemApproximation = SYSTEM_APPROXIMATION_NONE;
	basisType = POD;
	lsSolver = QR;
	dimension = 0;
}

SnapshotsData::SnapshotsData()
{
	normalizeSnaps = NORMALIZE_FALSE;
	incrementalSnaps = INCREMENTAL_FALSE;
	subtractIC = SUBTRACT_IC_FALSE;
	relProjError = REL_PROJ_ERROR_OFF;
//	sampleFreq = 1; // this is now included in the snapshot ascii file	
	snapshotWeights = UNIFORM;
  
}

DataCompressionData::DataCompressionData()
{
	type = POD;
	podMethod = SVD;
	maxVecStorage = 0;
	energyOnly = ENERGY_ONLY_FALSE;	// if ROB computation should only compute total energy of snapshots
  tolerance = 1e-8;
}

GNATData::GNATData()
{

	nRobState = -1;
	robNonlinear = UNSPECIFIED_NONLIN;
	nRobNonlin = -1;
	nRobRes = -1;
	nRobJac = -1;

	robGreedy = UNSPECIFIED_GREEDY;
	nRobGreedy = 0;
  sampleNodeFactor = -1;	// specify default of 2.0 elsewhere
  nSampleNodes = 0;
	layers = 2;	// by default, include two layers of nodes (2nd-order flux)

	includeLiftFaces = NONE_LIFTFACE;

	computeGappyRes = YES_GAPPYRES;
	// if NO, only output things corresponding to the jacobian (first pod basis),
	// and assume podTpod = I. Useful when you want to use another basis for
	// determining sample nodes

	sampleMeshUsed = SAMPLE_MESH_USED;
	pseudoInverseNodes = 20;	// if pod computation should only compute total energy of snapshots
}

//------------------------------------------------------------------------------

void ModelReductionData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassInt<ModelReductionData>(ca, "Dimension", this, &ModelReductionData::dimension);
	new ClassToken<ModelReductionData> (ca, "Projection", this, reinterpret_cast<int
			ModelReductionData::*>(&ModelReductionData::projection), 3, "PetrovGalerkin", 0, "Galerkin", 1,
			"ProjError", 2);	// ProjErrorcomputes projection error onto a basis
	new ClassToken<ModelReductionData> (ca, "SystemApproximation", this, reinterpret_cast<int
			ModelReductionData::*>(&ModelReductionData::systemApproximation), 4,
			"None", 0, "GNAT", 1, "Collocation", 2, "Broyden", 3);
	new ClassToken<ModelReductionData> (ca, "BasisType", this, reinterpret_cast<int
			ModelReductionData::*>(&ModelReductionData::basisType), 3, "Snaps", 0, "POD", 1,
			"None", 2);
	new ClassToken<ModelReductionData> (ca, "LeastSquaresSolver", this, reinterpret_cast<int
			ModelReductionData::*>(&ModelReductionData::lsSolver), 2, "QR", 0, "NormalEquations", 1);
}

//------------------------------------------------------------------------------

void SnapshotsData::setup(const char *name, ClassAssigner *father)
{

	ClassAssigner *ca = new ClassAssigner(name, 6, father);
	new ClassToken<SnapshotsData> (ca, "NormalizeSnaps", this, reinterpret_cast<int
			SnapshotsData::*>(&SnapshotsData::normalizeSnaps), 2, "False", 0, "True", 1);
	new ClassToken<SnapshotsData> (ca, "IncrementalSnaps", this, reinterpret_cast<int
			SnapshotsData::*>(&SnapshotsData::incrementalSnaps), 2, "False", 0, "True", 1);
	new ClassToken<SnapshotsData> (ca, "SubtractIC", this, reinterpret_cast<int
			SnapshotsData::*>(&SnapshotsData::subtractIC), 2, "False", 0, "True", 1);
 	new ClassToken<SnapshotsData> (ca, "RelProjError", this, reinterpret_cast<int
			SnapshotsData::*>(&SnapshotsData::relProjError), 2, "Off", 0, "On", 1);
 // new ClassInt<SnapshotsData>(ca, "Frequency", this, &SnapshotsData::sampleFreq);

	new ClassToken<SnapshotsData> (ca, "SnapshotWeights", this, reinterpret_cast<int
			SnapshotsData::*>(&SnapshotsData::snapshotWeights), 2, "Uniform", 0, "RBF", 1);

	dataCompression.setup("DataCompression",ca);

}

void DataCompressionData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
	new ClassToken<DataCompressionData> (ca, "Type", this, reinterpret_cast<int DataCompressionData::*>(&DataCompressionData::type), 2, "POD", 0, "Balanced POD", 1);
	new ClassToken<DataCompressionData> (ca, "PODMethod", this, reinterpret_cast<int
			DataCompressionData::*>(&DataCompressionData::podMethod), 2, "SVD", 0, "Eig", 1);
  new ClassInt<DataCompressionData>(ca, "MaxNumStoredVectors", this, &DataCompressionData::maxVecStorage);
	new ClassToken<DataCompressionData> (ca, "EnergyOnly", this, reinterpret_cast<int
			DataCompressionData::*>(&DataCompressionData::energyOnly), 2, "False", 0, "True", 1);
  new ClassDouble<DataCompressionData>(ca, "Tolerance", this, &DataCompressionData::tolerance);

}

void GNATData::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 14, father);
	// optional: document
  new ClassInt<GNATData>(ca, "DimensionROBState", this, &GNATData::nRobState);	// default: full size

	new ClassToken<GNATData>(ca, "IncludeLiftDragFaces", this, reinterpret_cast<int GNATData::*>(&GNATData::includeLiftFaces), 3, "None", 0, "Specified", 1, "All", 2);	// ProjErrorcomputes projection error onto a basis


	new ClassToken<GNATData>(ca, "ROBNonlinear", this, reinterpret_cast<int GNATData::*>(&GNATData::robNonlinear), 4, "Unspecified", -1, "Residual", 0, "Jacobian", 1, "Both", 2);

  new ClassInt<GNATData>(ca, "DimensionROBNonlinear", this, &GNATData::nRobNonlin);	// default: 0 (make sure C file reads in everything)
  new ClassInt<GNATData>(ca, "DimensionROBResidual", this, &GNATData::nRobRes);	// default: nRobNonlin
  new ClassInt<GNATData>(ca, "DimensionROBJacobian", this, &GNATData::nRobJac);	// default: nRobNonlin


	new ClassToken<GNATData>(ca, "ROBGreedy", this, reinterpret_cast<int GNATData::*>(&GNATData::robGreedy), 4, "Unspecified", -1, "Residual", 0, "Jacobian", 1, "Both", 2);	// ProjErrorcomputes projection error onto a basis
  new ClassInt<GNATData>(ca, "DimensionROBGreedy", this, &GNATData::nRobGreedy);


  new ClassDouble<GNATData>(ca, "SampledNodeFactor", this, &GNATData::sampleNodeFactor); // default: 2
  new ClassInt<GNATData>(ca, "NumSampledNodes", this, &GNATData::nSampleNodes); // overrides SampleNodeFactor to determine the number of sample nodes
  new ClassInt<GNATData>(ca, "NumSampledMeshLayers", this, &GNATData::layers);	// default: 2

	// optional: undocumented
	new ClassToken<GNATData> (ca, "ComputeGappyRes", this, reinterpret_cast<int
			GNATData::*>(&GNATData::computeGappyRes), 2, "False", 0, "True", 1);	

	new ClassToken<GNATData> (ca, "SampledMeshUsed", this, reinterpret_cast<int
			GNATData::*>(&GNATData::sampleMeshUsed), 2, "False", 0, "True", 1);	

  new ClassInt<GNATData>(ca, "NumPseudoInvNodesAtATime", this, &GNATData::pseudoInverseNodes);	// how many nodes of the pseudo inverse are calculated at a time. If this is too high, memory problems may ensue.

}
//------------------------------------------------------------------------------

LinearizedData::LinearizedData()
{

  type = DEFAULT;
  padeReconst = FALSE;
  domain = TIME;
  initCond = DISPLACEMENT;
  amplification = 1.0;
  frequency = 10.0;
  stepsize = -1.0;
  stepsizeinitial = -1.0;
  eps = 1e-4;
  eps2 = 5.0;
  tolerance = 1e-8;
  strModesFile = "";
  modeNumber = 1;
  numSteps = 0;
  numPOD = 0;
  numStrModes = 0;
  refLength = 1;
  freqStep = 0;

}

//------------------------------------------------------------------------------

void LinearizedData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 16, father);

  new ClassToken<LinearizedData> (ca, "Type", this, reinterpret_cast<int LinearizedData::*>(&LinearizedData::type), 3, "Default", 0, "Rom", 1, "Forced", 2);
  new ClassToken<LinearizedData> (ca, "Domain", this, reinterpret_cast<int LinearizedData::*>(&LinearizedData::domain), 2, "Time", 0, "Frequency", 1);
  new ClassToken<LinearizedData> (ca, "InitialCondition", this, reinterpret_cast<int LinearizedData::*>(&LinearizedData::initCond), 2, "Displacement", 0, "Velocity", 1);

  new ClassDouble<LinearizedData>(ca, "Amplification", this, &LinearizedData::amplification);
  new ClassDouble<LinearizedData>(ca, "Frequency", this, &LinearizedData::frequency);
  new ClassDouble<LinearizedData>(ca, "FreqStep", this, &LinearizedData::freqStep);
  new ClassDouble<LinearizedData>(ca, "Eps", this, &LinearizedData::eps);
  new ClassDouble<LinearizedData>(ca, "Eps2", this, &LinearizedData::eps2);
  new ClassDouble<LinearizedData>(ca, "Tolerance", this, &LinearizedData::tolerance);
  new ClassStr<LinearizedData>(ca, "StrModes", this, &LinearizedData::strModesFile);
  new ClassInt<LinearizedData>(ca, "ExcMode", this, &LinearizedData::modeNumber);
  new ClassInt<LinearizedData>(ca, "NumSteps", this, &LinearizedData::numSteps);
  new ClassInt<LinearizedData>(ca, "NumPOD", this, &LinearizedData::numPOD);
  new ClassInt<LinearizedData>(ca, "NumStrModes", this, &LinearizedData::numStrModes);

  pade.setup("Pade", ca);
	dataCompression.setup("DataCompression", ca);

}

//------------------------------------------------------------------------------

PadeData::PadeData()
{

  freq1 = -1;
  freq2 = -1;
  freq3 = -1;
  freq4 = -1;
  freq5 = -1;
  freq6 = -1;
  freq7 = -1;
  freq8 = -1;
  freq9 = -1;
  freq10 = -1;
  freq11 = -1;

  freq[0] = -1;
  freq[1] = -1;
  freq[2] = -1;
  freq[3] = -1;
  freq[4] = -1;
  freq[5] = -1;
  freq[6] = -1;
  freq[7] = -1;
  freq[8] = -1;
  freq[9] = -1;
  freq[10] = -1;
  nPoints = 0;
  degNum = 3;
  degDen = 4;


}


//------------------------------------------------------------------------------

void PadeData::setup(const char *name, ClassAssigner *father)
{


 ClassAssigner *ca = new ClassAssigner(name, 14, father);
  new ClassDouble<PadeData>(ca, "Freq1", this, &PadeData::freq1);
  new ClassDouble<PadeData>(ca, "Freq2", this, &PadeData::freq2);
  new ClassDouble<PadeData>(ca, "Freq3", this, &PadeData::freq3);
  new ClassDouble<PadeData>(ca, "Freq4", this, &PadeData::freq4);
  new ClassDouble<PadeData>(ca, "Freq5", this, &PadeData::freq5);
  new ClassDouble<PadeData>(ca, "Freq6", this, &PadeData::freq6);
  new ClassDouble<PadeData>(ca, "Freq7", this, &PadeData::freq7);
  new ClassDouble<PadeData>(ca, "Freq8", this, &PadeData::freq8);
  new ClassDouble<PadeData>(ca, "Freq9", this, &PadeData::freq9);
  new ClassDouble<PadeData>(ca, "Freq10", this, &PadeData::freq10);
  new ClassDouble<PadeData>(ca, "Freq11", this, &PadeData::freq11);
  new ClassInt<PadeData>(ca, "NumPoints", this, &PadeData::nPoints);
  new ClassInt<PadeData>(ca, "L", this, &PadeData::degNum);
  new ClassInt<PadeData>(ca, "M", this, &PadeData::degDen);


}

//------------------------------------------------------------------------------

SurfaceData::SurfaceData()  {

  nx = 0.0;
  ny = 0.0;
  nz = 0.0;

  sBit = 0;

  computeForces = (ComputeForces) UNSPECIFIED;
  forceResults = NO;

  rotationID = -1;
  velocity = 0.0;

  type = (Type) UNSPECIFIED; 
  temp = -1.0;
  computeHeatFluxes = (ComputeHeatPower) UNSPECIFIED_HF;
  heatFluxResults = NO_HF; 
}

//------------------------------------------------------------------------------
Assigner *SurfaceData::getAssigner()  {

  ClassAssigner *ca = new ClassAssigner("normal", 12, nullAssigner);

  new ClassDouble<SurfaceData>(ca, "Nx", this, &SurfaceData::nx);
  new ClassDouble<SurfaceData>(ca, "Ny", this, &SurfaceData::ny);
  new ClassDouble<SurfaceData>(ca, "Nz", this, &SurfaceData::nz);
  new ClassToken<SurfaceData> (ca, "ComputeForces", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::computeForces), 2, "False", 0, "True", 1);
  new ClassToken<SurfaceData> (ca, "SeparateForces", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::forceResults), 2, "False", 0, "True", 1);
  new ClassToken<SurfaceData> (ca, "SeparateFile", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::forceResults), 2, "False", 0, "True", 1); //I think this variable is never used
  
  new ClassToken<SurfaceData> (ca, "ComputeHeatFlux", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::computeHeatFluxes), 2, "False", 0, "True", 1);
  new ClassToken<SurfaceData> (ca, "SeparateHeatFlux", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::heatFluxResults), 2, "False", 0, "True", 1);

  new ClassInt<SurfaceData>(ca, "VelocityID", this, &SurfaceData::rotationID);
  new ClassDouble<SurfaceData>(ca, "Velocity", this, &SurfaceData::velocity);

  new ClassToken<SurfaceData>(ca, "Type", this,
                              (int SurfaceData::*)(&SurfaceData::type), 2,
                              "Isothermal", ISOTHERMAL, "Adiabatic", ADIABATIC);
  new ClassDouble<SurfaceData>(ca, "Temperature", this, &SurfaceData::temp);

  return ca;
}

//------------------------------------------------------------------------------

void Surfaces::setup(const char *name)  {

  ClassAssigner *ca = new ClassAssigner(name, 0, 0);
  surfaceMap.setup("SurfaceData", 0);
}

//------------------------------------------------------------------------------

void Velocity::setup(const char *name)  {

  ClassAssigner *ca = new ClassAssigner(name, 0, 0);
  rotationMap.setup("RotationAxis", 0);
}

//------------------------------------------------------------------------------

void Volumes::setup(const char *name, ClassAssigner *father)  {
  ClassAssigner *ca = new ClassAssigner(name, 0, father);
  volumeMap.setup("VolumeData", 0);
}


//------------------------------------------------------------------------------

RotationData::RotationData()  {

  nx = 0.0;
  ny = 0.0;
  nz = 0.0;

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;

  omega = 0.0;
  infRadius = FALSE;
}

//------------------------------------------------------------------------------


Assigner *RotationData::getAssigner()  {

  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);

  new ClassDouble<RotationData>(ca, "Nx", this, &RotationData::nx);
  new ClassDouble<RotationData>(ca, "Ny", this, &RotationData::ny);
  new ClassDouble<RotationData>(ca, "Nz", this, &RotationData::nz);

  new ClassDouble<RotationData>(ca, "X0", this, &RotationData::x0);
  new ClassDouble<RotationData>(ca, "Y0", this, &RotationData::y0);
  new ClassDouble<RotationData>(ca, "Z0", this, &RotationData::z0);

  new ClassDouble<RotationData>(ca, "Omega", this, &RotationData::omega);
  new ClassToken<RotationData> (ca, "InfiniteRadius", this, reinterpret_cast<int RotationData::*>(&RotationData::infRadius), 2, "False" , 0, "True", 1);

  return ca;
}

//------------------------------------------------------------------------------

VolumeData::VolumeData()  {

  type = FLUID;
  fluidModelID = -1;

}

//------------------------------------------------------------------------------

Assigner *VolumeData::getAssigner()  {

  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassToken<VolumeData> (ca, "Type", this, reinterpret_cast<int VolumeData::*>(&VolumeData::type), 2,
                              "Fluid", 0, "Porous", 1);
  new ClassInt<VolumeData> (ca, "FluidID", this, &VolumeData::fluidModelID);

  porousMedia.setup("PorousMedium", ca);
  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

PorousMedia::PorousMedia()  {

  iprimex = 1.0;
  iprimey = 0.0;
  iprimez = 0.0;

  jprimex = 0.0;
  jprimey = 1.0;
  jprimez = 0.0;

  kprimex = 0.0;
  kprimey = 0.0;
  kprimez = 1.0;

  alphax = 0.0;
  alphay = 0.0;
  alphaz = 0.0;

  betax = 0.0;
  betay = 0.0;
  betaz = 0.0;

  idr = 0.01;
  ldr = 1.0;

}

//------------------------------------------------------------------------------

//Assigner *PorousMedia::getAssigner()  {
void PorousMedia::setup(const char *name, ClassAssigner *father)  {

  //ClassAssigner *ca = new ClassAssigner("normal", 17, nullAssigner);
  ClassAssigner *ca = new ClassAssigner(name, 17, father);

  new ClassDouble<PorousMedia>(ca, "Ix", this, &PorousMedia::iprimex);
  new ClassDouble<PorousMedia>(ca, "Iy", this, &PorousMedia::iprimey);
  new ClassDouble<PorousMedia>(ca, "Iz", this, &PorousMedia::iprimez);

  new ClassDouble<PorousMedia>(ca, "Jx", this, &PorousMedia::jprimex);
  new ClassDouble<PorousMedia>(ca, "Jy", this, &PorousMedia::jprimey);
  new ClassDouble<PorousMedia>(ca, "Jz", this, &PorousMedia::jprimez);

  new ClassDouble<PorousMedia>(ca, "Kx", this, &PorousMedia::kprimex);
  new ClassDouble<PorousMedia>(ca, "Ky", this, &PorousMedia::kprimey);
  new ClassDouble<PorousMedia>(ca, "Kz", this, &PorousMedia::kprimez);

  new ClassDouble<PorousMedia>(ca, "Alphax", this, &PorousMedia::alphax);
  new ClassDouble<PorousMedia>(ca, "Alphay", this, &PorousMedia::alphay);
  new ClassDouble<PorousMedia>(ca, "Alphaz", this, &PorousMedia::alphaz);

  new ClassDouble<PorousMedia>(ca, "Betax", this, &PorousMedia::betax);
  new ClassDouble<PorousMedia>(ca, "Betay", this, &PorousMedia::betay);
  new ClassDouble<PorousMedia>(ca, "Betaz", this, &PorousMedia::betaz);

  new ClassDouble<PorousMedia>(ca, "Idr", this, &PorousMedia::idr);

  new ClassDouble<PorousMedia>(ca, "Ldr", this, &PorousMedia::ldr);

  //return ca;

}

//------------------------------------------------------------------------------

InitialConditions::InitialConditions() {

  mach = -1.0;
  velocity = -1.0;
  alpha = 400.0;
  beta = 400.0;
  density = -1.0;
  pressure = -1.0;
  temperature = -1.0;

}

//------------------------------------------------------------------------------

void InitialConditions::setup(const char *name, ClassAssigner *father) {

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassDouble<InitialConditions>(ca, "Mach",        this, &InitialConditions::mach);
  new ClassDouble<InitialConditions>(ca, "Velocity",    this, &InitialConditions::velocity);
  new ClassDouble<InitialConditions>(ca, "Alpha",       this, &InitialConditions::alpha);
  new ClassDouble<InitialConditions>(ca, "Beta",        this, &InitialConditions::beta);
  new ClassDouble<InitialConditions>(ca, "Density",     this, &InitialConditions::density);
  new ClassDouble<InitialConditions>(ca, "Pressure",    this, &InitialConditions::pressure);
  new ClassDouble<InitialConditions>(ca, "Temperature", this, &InitialConditions::temperature);

}

//------------------------------------------------------------------------------

EmbeddedFramework::EmbeddedFramework() {

  intersectorName = PHYSBAM;
  structNormal = ELEMENT_BASED;
  eosChange = NODAL_STATE;
  forceAlg = RECONSTRUCTED_SURFACE;
  riemannNormal = STRUCTURE;
  phaseChangeAlg = AVERAGE;
  interfaceAlg = MID_EDGE;
  alpha = 0.1;

  nLevelset = 0;

  //debug variables
  crackingWithLevelset = OFF;
  coupling = TWOWAY;
  dim2Treatment = NO;    
  reconstruct = CONSTANT;

}

//------------------------------------------------------------------------------

void EmbeddedFramework::setup(const char *name) {

  ClassAssigner *ca = new ClassAssigner(name, 13, 0); //father);
  new ClassToken<EmbeddedFramework> (ca, "Intersector", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::intersectorName), 2,
                                      "PhysBAM", 0, "FRG", 1);
  new ClassToken<EmbeddedFramework> (ca, "StructureNormal", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::structNormal), 2,
                                      "ElementBased", 0, "NodeBased", 1);
  new ClassToken<EmbeddedFramework> (ca, "EOSChange", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::eosChange), 2,
                                      "NodalState", 0, "RiemannSolution", 1);
  new ClassToken<EmbeddedFramework> (ca, "SurrogateSurface", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::forceAlg), 2,
                                      "Reconstructed", 0, "ControlVolumeFace", 1);
  new ClassToken<EmbeddedFramework> (ca, "RiemannNormal", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::riemannNormal), 3,
                                      "Structure", 0, "Fluid", 1, "AveragedStructure", 2);
  new ClassToken<EmbeddedFramework> (ca, "PhaseChangeAlgorithm", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::phaseChangeAlg), 2, "Average", 0, "LeastSquares", 1);
  new ClassToken<EmbeddedFramework> (ca, "InterfaceAlgorithm", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::interfaceAlg), 2, "MidEdge", 0, "Intersection", 1);
  embedIC.setup("InitialConditions", ca); 

  new ClassDouble<EmbeddedFramework>(ca, "Alpha", this, &EmbeddedFramework::alpha);

  //debug variables
  new ClassToken<EmbeddedFramework> (ca, "CrackingWithLevelSet", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::crackingWithLevelset), 2,
                                      "Off", 0, "On", 1);
  new ClassToken<EmbeddedFramework> (ca, "Coupling", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::coupling), 2,
                                      "TwoWay", 0, "OneWay", 1);
  new ClassToken<EmbeddedFramework> (ca, "TwoDimension", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::dim2Treatment), 2,
                                      "No", 0, "Yes", 1);
  new ClassToken<EmbeddedFramework> (ca, "Reconstruction", this, reinterpret_cast<int EmbeddedFramework::*>(&EmbeddedFramework::reconstruct), 2,
                                      "Constant", 0, "Linear", 1);
}

//------------------------------------------------------------------------------
OneDimensionalInfo::OneDimensionalInfo(){

  coordType  = SPHERICAL;//CARTESIAN;
  volumeType = CONSTANT_VOLUME;//REAL_VOLUME;

  mode = NORMAL;

  maxDistance = 0.0;
  numPoints = 101;
  interfacePosition = 0.5;

  density1 = 1.0; velocity1 = 0.0; pressure1 = 1.0;
  density2 = 1.0; velocity2 = 0.0; pressure2 = 1.0;
  temperature1 = 1.0; temperature2 = 1.0;

  sourceTermOrder = 1;
}
//------------------------------------------------------------------------------
void OneDimensionalInfo::setup(const char *name){

  ClassAssigner *ca = new ClassAssigner(name, 12, 0);

  new ClassToken<OneDimensionalInfo>(ca, "Coordinates", this, reinterpret_cast<int OneDimensionalInfo::*>(&OneDimensionalInfo::coordType), 3, "Cartesian", 0, "Cylindrical", 1, "Spherical", 2);
  new ClassToken<OneDimensionalInfo>(ca, "Volumes", this, reinterpret_cast<int OneDimensionalInfo::*>(&OneDimensionalInfo::volumeType), 2, "Constant", 0, "Real", 1);
  new ClassDouble<OneDimensionalInfo>(ca, "Radius", this, &OneDimensionalInfo::maxDistance);
  new ClassInt<OneDimensionalInfo>(ca, "NumberOfPoints", this, &OneDimensionalInfo::numPoints);
  new ClassInt<OneDimensionalInfo>(ca, "SourceTermOrder", this, &OneDimensionalInfo::sourceTermOrder);
  new ClassDouble<OneDimensionalInfo>(ca, "InterfacePosition", this, &OneDimensionalInfo::interfacePosition);

  new ClassDouble<OneDimensionalInfo>(ca, "Density1", this, &OneDimensionalInfo::density1);
  new ClassDouble<OneDimensionalInfo>(ca, "Velocity1", this, &OneDimensionalInfo::velocity1);
  new ClassDouble<OneDimensionalInfo>(ca, "Pressure1", this, &OneDimensionalInfo::pressure1);
  new ClassDouble<OneDimensionalInfo>(ca, "Density0", this, &OneDimensionalInfo::density2);
  new ClassDouble<OneDimensionalInfo>(ca, "Velocity0", this, &OneDimensionalInfo::velocity2);
  new ClassDouble<OneDimensionalInfo>(ca, "Pressure0", this, &OneDimensionalInfo::pressure2);
  new ClassToken<OneDimensionalInfo>(ca, "Mode", this, reinterpret_cast<int OneDimensionalInfo::*>(&OneDimensionalInfo::mode), 2, "Normal", 0, "ConvergenceTest1", 1);

  programmedBurn.setup("ProgrammedBurn",ca);
  
}

//-----------------------------------------------------------------------------

ImplosionSetup::ImplosionSetup() {
  // for buckling of cylinder
  type = LINEAR;
  Prate = -1.0;
  Pinit = -1.0;
  intersector_freq = 1;
  tmax = -1.0;
}

//-----------------------------------------------------------------------------

void ImplosionSetup::setup(const char *name) {
  ClassAssigner *ca = new ClassAssigner(name, 5, 0);
  new ClassToken<ImplosionSetup>(ca, "Type", this, reinterpret_cast<int ImplosionSetup::*>(&ImplosionSetup::type), 2, "Linear", 0, "SmoothStep", 1);
  new ClassDouble<ImplosionSetup>(ca, "RampupRate", this, &ImplosionSetup::Prate);
  new ClassDouble<ImplosionSetup>(ca, "InitialPressure", this, &ImplosionSetup::Pinit);
  new ClassDouble<ImplosionSetup>(ca, "Tmax", this, &ImplosionSetup::tmax);
  new ClassInt<ImplosionSetup>(ca, "InterfaceTrackingFrequency", this, &ImplosionSetup::intersector_freq);
}

MultigridInfo::MultigridInfo() {

  fineMesh = "";
  coarseMesh = "";
  fineDec = "";
  coarseDec = "";
  packageFile = "";
  collectionFile = "";

  radius0 = 1.0;
  radiusf = 2.0;
  threshold = 10;
}

void MultigridInfo::setup(const char * name) {
  
  ClassAssigner *ca = new ClassAssigner(name, 9, 0);
  new ClassDouble<MultigridInfo>(ca, "Radius0", this, &MultigridInfo::radius0);
  new ClassDouble<MultigridInfo>(ca, "RadiusF", this, &MultigridInfo::radiusf);
  new ClassInt<MultigridInfo>(ca, "Threshold", this, &MultigridInfo::threshold);

  new ClassStr<MultigridInfo>(ca, "FineMesh", this, &MultigridInfo::fineMesh);
  new ClassStr<MultigridInfo>(ca, "CoarseMesh", this, &MultigridInfo::coarseMesh);

  new ClassStr<MultigridInfo>(ca, "FineDec", this, &MultigridInfo::fineDec);
  new ClassStr<MultigridInfo>(ca, "CoarseDec", this, &MultigridInfo::coarseDec);
  
  new ClassStr<MultigridInfo>(ca, "PackageFile", this, &MultigridInfo::packageFile);
  new ClassStr<MultigridInfo>(ca, "CollectionFile", this, &MultigridInfo::collectionFile);

}

//-----------------------------------------------------------------------------

IoData::IoData(Communicator *communicator)
{

  com = communicator;

}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{

  int c;

  while ((c = getopt(argc, argv, "v:")) != -1)
    switch (c) {
    case 'v':
      problem.verbose = atoi(optarg);
      break;
    }

  if (optind < argc - 1) {
    com->fprintf(stderr, "*** Error: options must come before file name\n");
    exit(-1);
  }

  if (optind == argc-1 && argc > 1)
    cmdFileName = argv[argc-1];
  else {
    com->fprintf(stderr, "*** Error: no command input file given\n");
    exit(-1);
  }

}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{

  input.setup("Input");
  output.setup("Output");
  prec.setup("Preconditioner");
  restart.setup("RestartParameters");
  problem.setup("Problem");
  ref.setup("ReferenceState");
  bc.setup("BoundaryConditions");
  eqs.setup("Equations");
  mf.setup("MultiPhase");
  schemes.setup("Space");
  ts.setup("Time");
  dgcl.setup("DGCL");

// Included (MB)
  sa.setup("SensitivityAnalysis");

  dmesh.setup("MeshMotion");
  rmesh.setup("Accelerated");
  aero.setup("Aeroelastic");
  forced.setup("Forced");
	rom.setup("NonlinearModelReduction");
	gnat.setup("GNAT");
	snapshots.setup("Snapshots");
  linearizedData.setup("Linearized");
  surfaces.setup("Surfaces");
  rotations.setup("Velocity");
  volumes.setup("Volumes");
  embed.setup("EmbeddedFramework");
  oneDimensionalInfo.setup("1DGrid");
  implosion.setup("ImplosionSetup");
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{

  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    com->fprintf(stderr,"*** Error: could not open \'%s\'\n", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    com->fprintf(stderr,"*** Error: command file contained parsing errors\n");
    exit(error);
  }
  fclose(cmdFilePtr);

  if (input.rstdata[0] != 0) {
    char *name = new char[strlen(input.prefix) + strlen(input.rstdata) + 1];
    if (strncmp(input.rstdata, "/", 1) == 0)
      sprintf(name, "%s", input.rstdata);
    else
      sprintf(name, "%s%s", input.prefix, input.rstdata);
//    FILE *fp = freopen(name, "r", stdin);
    FILE *fp = yyCmdfin = fopen(name, "r");
    if (!fp) {
      com->fprintf(stderr, "*** Error: could not open \'%s\'\n", name);
      exit(-1);
    }
    error = yyCmdfparse();
    if (error) {
      com->fprintf(stderr, "*** Error: parameter file contained parsing errors\n");
      exit(error);
    }
    fclose(fp);
  }

  resetInputValues();
  error = checkFileNames();
  error += checkInputValues();
  error += checkSolverValues(surfaces.surfaceMap.dataMap);
  if (error) {
    com->fprintf(stderr, "*** Error: command file contained %d error%s\n",
		 error, error>1? "s":"");
    exit(-1);
  }
//  printDebug();

  com->setMaxVerbose(problem.verbose);

}

//------------------------------------------------------------------------------

void IoData::resetInputValues()
{

  // part 0

   if (strcmp(input.connectivity, "") == 0 && strcmp(input.geometryprefix, "") != 0 ) {
     char* name = new char[strlen(input.geometryprefix) + strlen(".con") + 1];
     sprintf(name, "%s%s", input.geometryprefix, ".con");
     input.connectivity = new char[strlen(name) + 1];
     strcpy( (char*)input.connectivity , name );
     delete [] name;
   }
   
   if (strcmp(input.geometry, "") == 0 && strcmp(input.geometryprefix, "") != 0 ) {
     char* name = new char[strlen(input.geometryprefix) + strlen(".msh") + 1];
     sprintf(name, "%s%s", input.geometryprefix, ".msh");
     input.geometry = new char[strlen(name) + 1];
     strcpy( (char*)input.geometry , name );
     delete [] name;
   }

   if (strcmp(input.decomposition, "") == 0 && strcmp(input.geometryprefix, "") != 0 ) {
     char* name = new char[strlen(input.geometryprefix) + strlen(".dec") + 1];
     sprintf(name, "%s%s", input.geometryprefix, ".dec");
     input.decomposition = new char[strlen(name) + 1];
     strcpy( (char*)input.decomposition , name );
     delete [] name;
   }

   if (strcmp(input.cpumap, "") == 0 && strcmp(input.geometryprefix, "") != 0 ) {
     char* numcpu = new char[6];
     sprintf(numcpu, "%i", com->size());
     char* name = new char[strlen(input.geometryprefix) + strlen(".cpu") + strlen(numcpu) + 1];
     sprintf(name, "%s%s%s%s", input.geometryprefix,".",numcpu,"cpu");
     input.cpumap = new char[strlen(name) + 1];
     strcpy( (char*)input.cpumap , name );
     delete [] numcpu;
     delete [] name;
   }

   if (eqs.type == EquationsData::NAVIER_STOKES && 
       eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY &&
       ( eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
         eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES ) )  {
      if (strcmp(input.d2wall, "") == 0 && strcmp(input.geometryprefix, "") != 0 ) {
       char* name = new char[strlen(input.geometryprefix) + strlen(".dwall") + 1];
       sprintf(name, "%s%s", input.geometryprefix,".dwall");
       input.d2wall = new char[strlen(name) + 1];
       strcpy( (char*)input.d2wall , name );
       delete [] name;
      }
   }

  // part 1

  for (int i=0; i<ProblemData::SIZE; ++i)
    problem.type[i] = false;

  if (problem.alltype == ProblemData::_UNSTEADY_ ||
      problem.alltype == ProblemData::_NONLINEAR_ROM_ ||
      problem.alltype == ProblemData::_ACC_UNSTEADY_ ||
      problem.alltype == ProblemData::_UNSTEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ACC_UNSTEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_UNSTEADY_THERMO_ ||
      problem.alltype == ProblemData::_UNSTEADY_AEROTHERMOELASTIC_ ||
      problem.alltype == ProblemData::_FORCED_ ||
      problem.alltype == ProblemData::_ACC_FORCED_ ||
      problem.alltype == ProblemData::_ROLL_)
    problem.type[ProblemData::UNSTEADY] = true;

  if (problem.alltype == ProblemData::_ACC_UNSTEADY_ ||
      problem.alltype == ProblemData::_ACC_UNSTEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ACC_FORCED_)
    problem.type[ProblemData::ACCELERATED] = true;

  if (problem.alltype == ProblemData::_STEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_UNSTEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ACC_UNSTEADY_AEROELASTIC_ ||
      problem.alltype == ProblemData::_STEADY_AEROTHERMOELASTIC_ ||
      problem.alltype == ProblemData::_UNSTEADY_AEROTHERMOELASTIC_)
    problem.type[ProblemData::AERO] = true;

  if (problem.alltype == ProblemData::_STEADY_THERMO_ ||
      problem.alltype == ProblemData::_UNSTEADY_THERMO_ ||
      problem.alltype == ProblemData::_STEADY_AEROTHERMOELASTIC_ ||
      problem.alltype == ProblemData::_UNSTEADY_AEROTHERMOELASTIC_)
    problem.type[ProblemData::THERMO] = true;

  if (problem.alltype == ProblemData::_FORCED_ ||
      problem.alltype == ProblemData::_ACC_FORCED_)
    problem.type[ProblemData::FORCED] = true;

  if (problem.alltype == ProblemData::_ROLL_)
    problem.type[ProblemData::ROLL] = true;

  if (problem.alltype == ProblemData::_RBM_)
    problem.type[ProblemData::RBM] = true;

  if (problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_ ||
      problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_ ||
      problem.alltype == ProblemData::_ROB_CONSTRUCTION_ ||
      problem.alltype == ProblemData::_POD_CONSTRUCTION_ ||
      problem.alltype == ProblemData::_ROM_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ROM_ ||
      problem.alltype == ProblemData::_INTERPOLATION_ ||
      problem.alltype == ProblemData::_ROB_INNER_PRODUCT_ ||
			problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_ ||
			problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_1_ ||
			problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_STEP_2_ ||
			problem.alltype == ProblemData::_SURFACE_MESH_CONSTRUCTION_ || 
			problem.alltype == ProblemData::_SAMPLE_MESH_SHAPE_CHANGE_) 
    problem.type[ProblemData::LINEARIZED] = true;

  // part 2

  // Included (MB)
  if (problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) 
  {

    //
    // Check that the code is running within the "correct" limits
    //

    if (sa.method == SensitivityAnalysis::ADJOINT)
    {
      com->fprintf(stderr, " ----- SA >> SensitivityAnalysis.Method has to be set to Direct -----\n");
      exit(1);
    }

    if (dmesh.type != DefoMeshMotionData::BASIC) 
    {
      com->fprintf(stderr, " ----- SA >> MeshMotion.Type has to be set to Basic -----\n");
      exit(1);
    }

    if (schemes.bc.type != BoundarySchemeData::STEGER_WARMING) 
    {
      com->fprintf(stderr, " ----- SA >> Boundaries.Type has to be set to StegerWarming -----\n");
      exit(1);
    }

    if (eqs.fluidModel.fluid != FluidModelData::PERFECT_GAS) 
    {
      com->fprintf(stderr, " ----- SA >> Equations.FluidModel.Type has to be set to Perfect Gas -----\n");
      exit(1);
    }

    if (problem.mode == ProblemData::NON_DIMENSIONAL) 
    {
      com->fprintf(stderr, " ----- SA >> Problem.Mode has to be set to Dimensional -----\n");
      exit(1);
    }

    //
    // Overwite some parameters if needed
    //


    if (problem.prec == ProblemData::PRECONDITIONED)
    {
      if (ts.implicit.mvp != ImplicitData::FD)
        com->fprintf(stderr, " ----- SA >> Time.Implicit.Mvp set to FiniteDifference of Order 2 for Problems with Low-Mach Preconditioning -----\n");
      ts.implicit.mvp = ImplicitData::FD;
      ts.implicit.fdOrder = ImplicitData::SECOND_ORDER;
    }

    if ((eqs.type == EquationsData::NAVIER_STOKES) &&
        (eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
    {
      //---------------
      if (ts.implicit.tmcoupling == ImplicitData::WEAK)
      {
        com->fprintf(stderr, " ----- SA >> Time.Implicit.TurbulenceModelCoupling set to Strong -----\n");
        ts.implicit.tmcoupling = ImplicitData::STRONG;
      }
      //---------------
      if (ts.implicit.mvp != ImplicitData::FD)
      {
        com->fprintf(stderr, " ----- SA >> Time.Implicit.Mvp set to FiniteDifference for Navier-Stokes Problems with Turbulence -----\n");
      }
      //---------------
      ts.implicit.mvp = ImplicitData::FD;
      //---------------
      if (ts.implicit.fdOrder == ImplicitData::FIRST_ORDER)
      {
        com->fprintf(stderr, " ----- SA >> Second-Order Finite Differencing is Recommended -----\n");
      }
    }


    if ((ts.implicit.mvp == ImplicitData::H1) || (ts.implicit.mvp == ImplicitData::H1FD))
    {
      com->fprintf(stderr, " ----- SA >> Time.Implicit.MatrixVectorProduct set to FiniteDifference -----\n");
      ts.implicit.mvp = ImplicitData::FD;
    }


    if (ts.implicit.ffjacobian != ImplicitData::EXACT) {
      // The overwriting is silent because the feature is not documented.
      //com->fprintf(stderr, " ----- SA >> Time.Implicit.FluxJacobian set to Exact -----\n");
      ts.implicit.ffjacobian = ImplicitData::EXACT;
    }


    int trip;
    if ( eqs.tc.tr.bfix.x0 > eqs.tc.tr.bfix.x1 ||
         eqs.tc.tr.bfix.y0 > eqs.tc.tr.bfix.y1 ||
         eqs.tc.tr.bfix.z0 > eqs.tc.tr.bfix.z1 )
      trip = 0;
    else
      trip = 1;


    if (ts.implicit.mvp == ImplicitData::H2 && trip)
    {
      com->fprintf(stderr,
                   " ----- SA >> MatrixVectorProduct set to"
                   " FiniteDifference to account for tripping -----\n");
      ts.implicit.mvp = ImplicitData::FD;
    }


  } // END if (problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_)

  //
  // Check parameters for the matrix-vector product in implicit simulations.
  //

//  if ((problem.prec == ProblemData::PRECONDITIONED) ||
//      (ts.prec == TsData::PREC))
//  {
//    if (ts.implicit.mvp == ImplicitData::H2)
//    {
//      com->fprintf(stderr, "*** Warning: Exact Matrix-Vector Product not supported with Low-Mach Preconditioning.\n");
//      com->fprintf(stderr, "             Second Order Finite Difference will be used.\n");
//      ts.implicit.mvp = ImplicitData::FD;
//      ts.implicit.fdOrder = ImplicitData::SECOND_ORDER;
//    }
//  } // END of if ((problem.prec == ProblemData::PRECONDITIONED) || ...

  if (problem.prec == ProblemData::PRECONDITIONED && 
      ts.prec == TsData::PREC && problem.type[ProblemData::UNSTEADY])
  {
    if (ts.dualtimestepping == TsData::OFF) {
      com->fprintf(stderr, "*** Warning: Dual Time-stepping required for unsteady Low-Mach Preconditioning.\n");
      com->fprintf(stderr, "             Turning on Dual Time-stepping.\n");
      ts.dualtimestepping = TsData::ON;
    }
  } // END of if (problem.prec == ProblemData::PRECONDITIONED && ...

  if (schemes.ns.flux != SchemeData::ROE)
  {
    if (ts.implicit.mvp == ImplicitData::H2)
    {
      com->fprintf(stderr, "*** Warning: Exact Matrix-Vector Product only supported with Roe flux.\n");
      com->fprintf(stderr, "             Second Order Finite Difference will be used.\n");
      ts.implicit.mvp = ImplicitData::FD;
      ts.implicit.fdOrder = ImplicitData::SECOND_ORDER;
    }
    if (eqs.fluidModel.fluid != FluidModelData::PERFECT_GAS && 
        eqs.fluidModel.fluid != FluidModelData::STIFFENED_GAS)
    {
      com->fprintf(stderr, "*** Warning: Roe flux has to be used for Tait or JWL simulations.\n");
      schemes.ns.flux = SchemeData::ROE;
    }
    if(!eqs.fluidModelMap.dataMap.empty())
    {
      map<int, FluidModelData *>::iterator it;
      for (it=eqs.fluidModelMap.dataMap.begin(); it!=eqs.fluidModelMap.dataMap.end(); it++){
	if (it->second->fluid != FluidModelData::PERFECT_GAS &&
            it->second->fluid != FluidModelData::STIFFENED_GAS)
	{
	  com->fprintf(stderr, "*** Warning: Roe flux has to be used for Tait or JWL simulations.\n");
	  schemes.ns.flux = SchemeData::ROE;
	}    
      }
    }
  } // END of if (schemes.ns.flux != SchemeData::ROE)

  if (ts.implicit.mvp == ImplicitData::H2)
  {
    if (problem.prec != ProblemData::PRECONDITIONED)
    // The overwriting is silent because ffjacobian is a "slave" flag.
    ts.implicit.ffjacobian = ImplicitData::EXACT;
  }

  //
  // Part 3
  //

  if (problem.type[ProblemData::AERO] || problem.type[ProblemData::THERMO] ||
      problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ROM_AEROELASTIC_)
    problem.mode = ProblemData::DIMENSIONAL;

  if (!problem.type[ProblemData::UNSTEADY]) {
    ts.implicit.type = ImplicitData::BACKWARD_EULER;
  }

  if (ts.type == TsData::IMPLICIT &&
      (ts.implicit.newton.failsafe == NewtonData<KspFluidData>::YES ||
        ts.implicit.newton.failsafe == NewtonData<KspFluidData>::ALWAYS)) {
    if (schemes.ns.reconstruction == SchemeData::CONSTANT)
      ts.implicit.newton.failsafe = NewtonData<KspFluidData>::NO;
    //else
      //schemes.ns.advectiveOperator = SchemeData::FE_GALERKIN;
  }

  if (problem.type[ProblemData::AERO] && !problem.type[ProblemData::UNSTEADY]) {
    //ts.implicit.normals = ImplicitData::LATEST_CFG;
    //ts.implicit.velocities = ImplicitData::ZERO;
    dgcl.normals = DGCLData::IMPLICIT_LATEST_CFG;
    dgcl.velocities = DGCLData::IMPLICIT_ZERO;
  }

  if (bc.wall.integration == BcsWallData::AUTO) {
    if (eqs.type == EquationsData::NAVIER_STOKES &&
	eqs.tc.type == TurbulenceClosureData::NONE)
      bc.wall.integration = BcsWallData::FULL;
    else
      bc.wall.integration = BcsWallData::WALL_FUNCTION;
  }

  if (problem.type[ProblemData::THERMO])
    bc.wall.type = BcsWallData::ISOTHERMAL;

  // to avoid having inlet nodes when computing Internal BCs
  if(bc.inlet.type == BcsFreeStreamData::INTERNAL)
    schemes.bc.type = BoundarySchemeData::STEGER_WARMING;

  int nCoarseFreq = 0;
  if (problem.alltype == ProblemData::_POD_CONSTRUCTION_) {

    // Assign the values for the coarse freq to an array
    linearizedData.pade.freq[0] = linearizedData.pade.freq1;
    linearizedData.pade.freq[1] = linearizedData.pade.freq2;
    linearizedData.pade.freq[2] = linearizedData.pade.freq3;
    linearizedData.pade.freq[3] = linearizedData.pade.freq4;
    linearizedData.pade.freq[4] = linearizedData.pade.freq5;
    linearizedData.pade.freq[5] = linearizedData.pade.freq6;
    linearizedData.pade.freq[6] = linearizedData.pade.freq7;
    linearizedData.pade.freq[7] = linearizedData.pade.freq8;
    linearizedData.pade.freq[8] = linearizedData.pade.freq9;
    linearizedData.pade.freq[9] = linearizedData.pade.freq10;
    linearizedData.pade.freq[10] = linearizedData.pade.freq11;

    for (int i=0; i<linearizedData.pade.num; i++) {
      if (linearizedData.pade.freq[i] >= 0)
        nCoarseFreq++;

    }
    if (nCoarseFreq > 0)
      linearizedData.padeReconst = LinearizedData::TRUE;

  }


  if (problem.type[ProblemData::LINEARIZED] && ts.form == TsData::HYBRID) {
    com->fprintf(stderr, "*** Error: the hybrid form is not supported in the linearized module \n");
      exit(1);
  }

  
  if (linearizedData.padeReconst == LinearizedData::TRUE) {
    if (ts.implicit.newton.ksp.ns.type != KspData::GCR) {
      com->fprintf(stderr, "*** Error: for a Pade Reconstruction, a GCR solver needs to be used \n");
      exit(1);
    }

    if (nCoarseFreq < linearizedData.pade.nPoints) {

      com->fprintf(stderr, "*** Error: the number of specified coarse-grid frequencies is lower than the number of points needed for each Pade Reconstruction\n");
      exit(1);
    }

    if ((linearizedData.pade.degNum + linearizedData.pade.degDen + 1)%linearizedData.pade.nPoints != 0 ) {
      com->fprintf(stderr, "*** Error: In the Pade reconstruction, L+M+1 has to be a multiple of the number of points used \n");
      exit(1);
    }

    for (int i=0; i < nCoarseFreq-1; i++) {

      if (linearizedData.pade.freq[i] >= linearizedData.pade.freq[i+1]) {


        com->fprintf(stderr, "*** Error: the coarse frequencies specified for the Pade Reconstruction are not sorted in an ascending order\n");
        exit(1);
      }
    }

    if (ts.implicit.newton.ksp.ns.type == KspData::CG) {
      com->fprintf(stderr, "*** Error: CG solver cannot be used for unsymmetric problems\n");
      exit(1);
    }
  }
   
  // Sparse Grid Generation does not call the flow solver
  if(problem.alltype == ProblemData::_SPARSEGRIDGEN_)
    problem.mode = ProblemData::DIMENSIONAL;

}

//------------------------------------------------------------------------------

int IoData::checkFileNames()
{

  int error = 0;

  // no input files for Sparse Grid generation, hence no check
  // or if we are doing a one dimensional problem!
  if(problem.alltype == ProblemData::_SPARSEGRIDGEN_ || problem.alltype == ProblemData::_ONE_DIMENSIONAL_)
    return 0;
  
  // flow solver requires input files, hence check
  if (strcmp(input.connectivity, "") == 0) {
    com->fprintf(stderr, "*** Error: no global file given\n");
    ++error;
  }
  if (strcmp(input.geometry, "") == 0) {
    com->fprintf(stderr, "*** Error: no geometry file given\n");
    ++error;
  }
  if (strcmp(input.decomposition, "") == 0) {
    com->fprintf(stderr, "*** Error: no decomposition file given\n");
    ++error;
  }
  if (strcmp(input.cpumap, "") == 0) {
    com->fprintf(stderr, "*** Error: no CPU map file given\n");
    ++error;
  }
  if ((problem.type[ProblemData::AERO] || problem.type[ProblemData::THERMO]) &&
       problem.framework!=ProblemData::EMBEDDED && strcmp(input.match, "") == 0) {
    com->fprintf(stderr, "*** Error: no matcher file given\n");
    ++error;
  }
  if (eqs.type == EquationsData::NAVIER_STOKES &&
      eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    if (eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
      bc.wall.integration = BcsWallData::WALL_FUNCTION;
    if (eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS &&
	strcmp(input.d2wall, "") == 0) {
      com->fprintf(stderr, "*** Error: no distance to wall file given\n");
      ++error;
    }
   if (eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES &&
	strcmp(input.d2wall, "") == 0) {
      com->fprintf(stderr, "*** Error: no distance to wall file given\n");
    ++error;
   }
    if (bc.wall.integration == BcsWallData::WALL_FUNCTION) {
      if (bc.wall.delta < 0.0) {
	com->fprintf(stderr, "*** Error: no delta value given\n");
	++error;
      }
    }
    else
      bc.wall.delta = 0.0;
  }
  else {
    output.transient.nutturb = "";
    output.transient.kturb = "";
    output.transient.epsturb = "";
    output.transient.eddyvis = "";
  }

  return error;

}

//------------------------------------------------------------------------------
// alpha_bar * ref.rv.alpha = alpha_SI
// alpha_nonSI * scaling = alpha_SI
// alpha_bar * ref.rv.talpha = alpha_nonSI
// -> ref.rv.talpha = ref.rv.alpha / scaling

int IoData::checkInputValues()
{

  int error = 0;

  // input values for flow solver
  error += checkInputValuesAllEquationsOfState();

  // no need for all input values for Sparse Grid generation
  if(problem.alltype == ProblemData::_SPARSEGRIDGEN_){
    //eqs.fluidModel  =  mf.fluidModel;
    //eqs.fluidModel2 =  mf.fluidModel2;
    return checkInputValuesSparseGrid(mf.sparseGrid);
  }

  // Check the input values for programmed burn first,
  // because if we are doing programmed burn we can set some
  // initial conditions that do not need to be specified in the
  // input file
  checkInputValuesProgrammedBurn();
  if (problem.alltype == ProblemData::_ONE_DIMENSIONAL_) {

    bc.inlet.pressure = bc.outlet.pressure;
    bc.inlet.density = bc.outlet.density;
    bc.inlet.alpha = bc.inlet.beta = 0.0;
    bc.inlet.temperature = bc.outlet.temperature;
    setupOneDimensional();
  }
    

  error += checkInputValuesAllInitialConditions();

  error += checkInputValuesEssentialBC();
  
  // check input values for the Embedded Framework.
  if(problem.framework == ProblemData::EMBEDDED) 
    error += checkInputValuesEmbeddedFramework();

  error += checkInputValuesNonDimensional();
  error += checkInputValuesDimensional(surfaces.surfaceMap.dataMap);

  checkInputValuesTurbulence();
  checkInputValuesDefaultOutlet();
  nonDimensionalizeAllEquationsOfState(); 
  nonDimensionalizeViscosityModel(eqs.viscosityModel);
  nonDimensionalizeThermalCondModel(eqs.thermalCondModel);
  nonDimensionalizeForcedMotion();

  if(problem.alltype == ProblemData::_ONE_DIMENSIONAL_) {
    nonDimensionalizeOneDimensionalProblem();
    return 0;
  }
  nonDimensionalizeAllInitialConditions();

  bc.inlet.alpha *= acos(-1.0) / 180.0;
  bc.inlet.beta *= acos(-1.0) / 180.0;
  bc.outlet.alpha *= acos(-1.0) / 180.0;
  bc.outlet.beta *= acos(-1.0) / 180.0;

  if (aero.pressure < 0.0)
    aero.pressure = bc.inlet.pressure;


  eqs.tc.tr.bfix.x0 /= ref.rv.tlength;
  eqs.tc.tr.bfix.x1 /= ref.rv.tlength;
  eqs.tc.tr.bfix.y0 /= ref.rv.tlength;
  eqs.tc.tr.bfix.y1 /= ref.rv.tlength;
  eqs.tc.tr.bfix.z0 /= ref.rv.tlength;
  eqs.tc.tr.bfix.z1 /= ref.rv.tlength;


  return error;

}
//------------------------------------------------------------------------------

int IoData::checkInputValuesAllEquationsOfState(){

  int error = 0;

  // count number of phases and check input values of each EOS
  // copy FluidModel[0] in FluidModel, if the first exists (always true unless a one-phase simulation is done using old inputfile type)
  eqs.numPhase = 0;
  bool fluidModelZeroExists = false;

  if(!eqs.fluidModelMap.dataMap.empty()){
    map<int, FluidModelData *>::iterator it;
    for (it=eqs.fluidModelMap.dataMap.begin(); it!=eqs.fluidModelMap.dataMap.end(); it++){
      eqs.numPhase++;
      error += checkInputValuesEquationOfState(*(it->second), it->first);
      // copy FluidModel[0] in FluidModel
      if(it->first == 0){
        fluidModelZeroExists = true;
        eqs.fluidModel = *(it->second);
      }
    }
    // if FluidModel[0] was not specified for a multiphase flow simulation
    if(!fluidModelZeroExists && eqs.numPhase > 0){
      error++;
      com->fprintf(stderr, "*** Error: FluidModel[0] must be specified for a multiphase flow simulation\n");
    }
  }
  // if no map, then take the default FluidModel (cannot be user-specified!)
  if(eqs.numPhase == 0){
    error += checkInputValuesEquationOfState(eqs.fluidModel, -1);
    eqs.numPhase = 1;
  }

  // check if fluid-fluid interfaces are expected
  if(eqs.numPhase > 1 && 
     (!mf.multiInitialConditions.sphereMap.dataMap.empty() || 
      !mf.multiInitialConditions.planeMap.dataMap.empty()  ||
      !mf.multiInitialConditions.prismMap.dataMap.empty() || 
      !input.oneDimensionalInput.dataMap.empty()))
    mf.interfaceType = MultiFluidData::FF;
  else
    mf.interfaceType = MultiFluidData::FSF;

  return error;

}

//------------------------------------------------------------------------------

int IoData::checkInputValuesAllInitialConditions(){

  int error = 0;

  // check input values of initial conditions for volumeData
  if(!volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator it;
    for (it=volumes.volumeMap.dataMap.begin(); it!=volumes.volumeMap.dataMap.end();it++)
      if(it->second->type==VolumeData::FLUID)
        error += checkInputValuesInitialConditions(it->second->initialConditions, it->second->fluidModelID);
  }

  // check input values of initial conditions for multiFluidData
  if(!mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator it;
    for (it=mf.multiInitialConditions.sphereMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.sphereMap.dataMap.end();
         it++){
      error += checkInputValuesInitialConditions(it->second->initialConditions, it->second->fluidModelID);
      
      if(it->second->radius <= 0.0){
        error++;
        com->fprintf(stderr, "*** Error: a positive radius must be specified for the sphere[%d] initial conditions\n", it->first);
      }
    }
  }
  if(!mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator it;
    for (it=mf.multiInitialConditions.prismMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.prismMap.dataMap.end();
         it++){
      error += checkInputValuesInitialConditions(it->second->initialConditions, it->second->fluidModelID);

      it->second->w_x = it->second->X1-it->second->X0;
      it->second->w_y = it->second->Y1-it->second->Y0;
      it->second->w_z = it->second->Z1-it->second->Z0;

      it->second->cen_x = (it->second->X1+it->second->X0)*0.5;
      it->second->cen_y = (it->second->Y1+it->second->Y0)*0.5;
      it->second->cen_z = (it->second->Z1+it->second->Z0)*0.5;

      
      if(it->second->w_x <= 0.0 || 
	 it->second->w_y <= 0.0 || 
	 it->second->w_z <= 0.0 ){
        error++;
        com->fprintf(stderr, "*** Error: a size must be specified for box[%d] initial conditions\n", it->first);
      }
    }
  }
  if(!mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator it;
    for (it=mf.multiInitialConditions.planeMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.planeMap.dataMap.end();
         it++){
      error += checkInputValuesInitialConditions(it->second->initialConditions, it->second->fluidModelID);
      double norm = it->second->cen_x*it->second->cen_x+it->second->cen_y*it->second->cen_y+it->second->cen_z*it->second->cen_z;
      if(norm <= 0.0){
        error++;
        com->fprintf(stderr, "*** Error: a positive vector norm must be specified for the plane[%d] initial conditions\n", it->first);
      }
    }
  }

  // check input values of initial conditions for EmbeddedStructure
  if(!embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator it;
    for (it=embed.embedIC.pointMap.dataMap.begin();
         it!=embed.embedIC.pointMap.dataMap.end();
         it++) {
//      if(it->second->fluidModelID==0)
        //com->fprintf(stderr,"*** WARNING: FluidID on a user-specified point (for initial condition setup) \n");
      error += checkInputValuesInitialConditions(it->second->initialConditions, it->second->fluidModelID);
    }
  }

  embed.nLevelset = (embed.crackingWithLevelset==EmbeddedFramework::ON) ? 1 : 0;

  // count number of levelsets (consider only bubbles!) for the Embedded Framework.
  set<int> usedModels; 
  for (map<int, SphereData *>::iterator it=mf.multiInitialConditions.sphereMap.dataMap.begin();
       it!=mf.multiInitialConditions.sphereMap.dataMap.end();
       it++)
    usedModels.insert(it->second->fluidModelID);

  for (map<int, PrismData *>::iterator it=mf.multiInitialConditions.prismMap.dataMap.begin();
       it!=mf.multiInitialConditions.prismMap.dataMap.end();
       it++)
    usedModels.insert(it->second->fluidModelID);
 
  if (!input.oneDimensionalInput.dataMap.empty()) {
    if (input.oneDimensionalInput.dataMap.size() > 1) 
      std::cout << "Warning: having more than one 1D->3D remap has not been considered" << std::endl;

    for (map<int, OneDimensionalInputData *>::iterator it=input.oneDimensionalInput.dataMap.begin();
         it != input.oneDimensionalInput.dataMap.end(); it++) {

      for (map<int,FluidRemapData*>::iterator it2 = (it->second)->fluidRemap.dataMap.begin();
           it2 != (it->second)->fluidRemap.dataMap.end(); it2++) {

        if (it2->second->newID > 0)
          usedModels.insert(it2->second->newID);
      }
    }      
 
    //embed.nLevelset++;
  }

  int nModels = usedModels.size();
  if (nModels > 0) {
    int minModel = *(usedModels.begin());
    int maxModel = *(--(set<int>::iterator)(usedModels.end()));
    if(minModel!=1 || nModels != maxModel - minModel + 1) {
      com->fprintf(stderr,"*** Error: FluidID(s) specified under 'volume(s)' cannot be accepted!\n");
      com->fprintf(stderr,"***        Currently, when fluid-fluid interfaces are present in an embedded FSI simulation, the \n");
      com->fprintf(stderr,"***        FluidID(s) (integer) specified under MultiPhase for spheres or planes must be consecutive \n");
      com->fprintf(stderr,"***        starting at 1!\n");
      com->fprintf(stderr,"*** Re-order the fluid models to satisfy this constraint, then re-run this simulation!\n");
      error++;
    } else 
      embed.nLevelset += nModels;
  }

  return error;

}

//------------------------------------------------------------------------------

void IoData::nonDimensionalizeAllEquationsOfState(){

  if (problem.mode == ProblemData::NON_DIMENSIONAL) return;

  //non-dimensionalize the time interval for output
  output.transient.frequency_dt /= ref.rv.time;
  output.restart.frequency_dt /= ref.rv.time;

  if(!eqs.fluidModelMap.dataMap.empty()){
    map<int, FluidModelData *>::iterator it;
    for (it=eqs.fluidModelMap.dataMap.begin(); it!=eqs.fluidModelMap.dataMap.end();it++)
      nonDimensionalizeFluidModel(*(it->second));
  }
  nonDimensionalizeFluidModel(eqs.fluidModel);

}

//------------------------------------------------------------------------------

void IoData::nonDimensionalizeAllInitialConditions(){

  if (problem.mode == ProblemData::NON_DIMENSIONAL) return;

  // non-dimensionalize initial conditions for volumeData
  if(!volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator it;
    for (it=volumes.volumeMap.dataMap.begin(); it!=volumes.volumeMap.dataMap.end();it++)
      if(it->second->type==VolumeData::FLUID)
        nonDimensionalizeInitialConditions(it->second->initialConditions);
  }

  // non-dimensionalize initial conditions for multiFluidData
  if(!mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator it;
    for (it=mf.multiInitialConditions.sphereMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.sphereMap.dataMap.end();
         it++){
      nonDimensionalizeInitialConditions(it->second->initialConditions);
      it->second->cen_x  /= ref.rv.length;
      it->second->cen_y  /= ref.rv.length;
      it->second->cen_z  /= ref.rv.length;
      it->second->radius /= ref.rv.length;

      // Nondimensionalize  programmed burn
      it->second->programmedBurn.ignitionX0 /= ref.rv.length;
      it->second->programmedBurn.ignitionY0 /= ref.rv.length;
      it->second->programmedBurn.ignitionZ0 /= ref.rv.length;
      it->second->programmedBurn.ignitionTime /= ref.rv.time;
      it->second->programmedBurn.e0 /= ref.rv.energy;
      it->second->programmedBurn.cjDetonationVelocity /= ref.rv.velocity;
      it->second->programmedBurn.cjPressure /= ref.rv.pressure;
      it->second->programmedBurn.cjDensity /= ref.rv.density;
      it->second->programmedBurn.cjEnergy /= ref.rv.energy;
    }
  }

  // non-dimensionalize initial conditions for multiFluidData
  if(!mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator it;
    for (it=mf.multiInitialConditions.prismMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.prismMap.dataMap.end();
         it++){
      nonDimensionalizeInitialConditions(it->second->initialConditions);
      it->second->cen_x  /= ref.rv.length;
      it->second->cen_y  /= ref.rv.length;
      it->second->cen_z  /= ref.rv.length;
      it->second->w_x /= ref.rv.length;
      it->second->w_y /= ref.rv.length;
      it->second->w_z /= ref.rv.length;

      // Nondimensionalize  programmed burn
      it->second->programmedBurn.ignitionX0 /= ref.rv.length;
      it->second->programmedBurn.ignitionY0 /= ref.rv.length;
      it->second->programmedBurn.ignitionZ0 /= ref.rv.length;
      it->second->programmedBurn.ignitionTime /= ref.rv.time;
      it->second->programmedBurn.e0 /= ref.rv.energy;
      it->second->programmedBurn.cjDetonationVelocity /= ref.rv.velocity;
      it->second->programmedBurn.cjPressure /= ref.rv.pressure;
      it->second->programmedBurn.cjDensity /= ref.rv.density;
      it->second->programmedBurn.cjEnergy /= ref.rv.energy;
    }
  }

  if(!mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator it;
    for (it=mf.multiInitialConditions.planeMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.planeMap.dataMap.end();
         it++){
      nonDimensionalizeInitialConditions(it->second->initialConditions);
      it->second->cen_x /= ref.rv.length;
      it->second->cen_y  /= ref.rv.length;
      it->second->cen_z  /= ref.rv.length;
    }
  }

  implosion.Pinit /= ref.rv.pressure;
  implosion.Prate /= ref.rv.pressure/ref.rv.time;
  implosion.tmax  /= ref.rv.time;

  // non-dimensionalize initial conditions for EmbeddedStructure
  if(!embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator it;
    for (it=embed.embedIC.pointMap.dataMap.begin();
         it!=embed.embedIC.pointMap.dataMap.end();
         it++) {
      nonDimensionalizeInitialConditions(it->second->initialConditions);
      it->second->x /= ref.rv.length;
      it->second->y /= ref.rv.length;
      it->second->z /= ref.rv.length;
    }
  }
}

//------------------------------------------------------------------------------

void IoData::setupOneDimensional() {

  if(!mf.multiInitialConditions.sphereMap.dataMap.empty()){
    if (mf.multiInitialConditions.sphereMap.dataMap.size() > 1) {
      fprintf(stderr,"*** Error: more than one sphere specified for a one dimensional problem!\n");
      exit(1);
    }
      
    map<int, SphereData *>::iterator it;
    for (it=mf.multiInitialConditions.sphereMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.sphereMap.dataMap.end();
         it++){
      if (it->second->cen_x != 0.0 ||
	  it->second->cen_x != 0.0 ||
	  it->second->cen_x != 0.0) {
	fprintf(stderr,"*** Error: non zero center specified for 1d spherical problem!\n");
	exit(1);
      }

      oneDimensionalInfo.interfacePosition = it->second->radius;
      memcpy(&oneDimensionalInfo.programmedBurn, &it->second->programmedBurn, sizeof(it->second->programmedBurn));
      
      oneDimensionalInfo.fluidId2 = it->second->fluidModelID;
      oneDimensionalInfo.density1 = it->second->initialConditions.density;
      oneDimensionalInfo.velocity1 = it->second->initialConditions.velocity;
      oneDimensionalInfo.temperature1 = it->second->initialConditions.temperature;
      oneDimensionalInfo.pressure1 = it->second->initialConditions.pressure;

      oneDimensionalInfo.density2 = bc.outlet.density;
      oneDimensionalInfo.velocity2 = 0.0;//bc.outlet.velocity;
      oneDimensionalInfo.temperature2 = bc.outlet.temperature;
      oneDimensionalInfo.pressure2 = bc.outlet.pressure;
   
      
    }
  }
}

//------------------------------------------------------------------------------

void IoData:: nonDimensionalizeOneDimensionalProblem(){

  if (problem.mode == ProblemData::NON_DIMENSIONAL) return;

  oneDimensionalInfo.maxDistance /= ref.rv.length;
  oneDimensionalInfo.interfacePosition /= ref.rv.length;

  oneDimensionalInfo.density1  /= ref.rv.density;
  oneDimensionalInfo.velocity1 /= ref.rv.velocity;
  oneDimensionalInfo.pressure1 /= ref.rv.pressure;
  oneDimensionalInfo.temperature1  /= ref.rv.temperature;
  oneDimensionalInfo.density2  /= ref.rv.density;
  oneDimensionalInfo.velocity2 /= ref.rv.velocity;
  oneDimensionalInfo.pressure2 /= ref.rv.pressure;
  oneDimensionalInfo.temperature2 /= ref.rv.temperature;


  oneDimensionalInfo.programmedBurn.ignitionX0 /= ref.rv.length;
  oneDimensionalInfo.programmedBurn.ignitionY0 /= ref.rv.length;
  oneDimensionalInfo.programmedBurn.ignitionZ0 /= ref.rv.length;
  oneDimensionalInfo.programmedBurn.ignitionTime /= ref.rv.time;
  oneDimensionalInfo.programmedBurn.e0 /= ref.rv.energy;
  oneDimensionalInfo.programmedBurn.cjDetonationVelocity /= ref.rv.velocity;
  oneDimensionalInfo.programmedBurn.cjPressure /= ref.rv.pressure;
  oneDimensionalInfo.programmedBurn.cjDensity /= ref.rv.density;
  oneDimensionalInfo.programmedBurn.cjEnergy /= ref.rv.energy;
}

//------------------------------------------------------------------------------

void IoData:: nonDimensionalizeForcedMotion(){

  if (problem.mode == ProblemData::NON_DIMENSIONAL) return;

  //heaving
  forced.hv.ax /= ref.rv.length;
  forced.hv.ay /= ref.rv.length;
  forced.hv.az /= ref.rv.length;

  //pitching
  forced.pt.x1 /= ref.rv.length;
  forced.pt.y1 /= ref.rv.length;
  forced.pt.z1 /= ref.rv.length;
  forced.pt.x2 /= ref.rv.length;
  forced.pt.y2 /= ref.rv.length;
  forced.pt.z2 /= ref.rv.length;
}

//------------------------------------------------------------------------------
int IoData::checkInputValuesNonDimensional()
{
  int error = 0;
  if (problem.mode == ProblemData::NON_DIMENSIONAL) {

    // no multiphase flow in non-dimensional
    if(eqs.numPhase > 1 && problem.framework != ProblemData::EMBEDDED){ 
      com->fprintf(stderr, "*** Error: multiphase flow are possible only in Dimensional Mode \n");
      ++error;
      return error;
    }

    if (eqs.type != EquationsData::EULER) {
      if (ref.reynolds_mu < 0.0) {
        com->fprintf(stderr, "*** Error: no valid Reynolds number (%d) given\n", ref.reynolds_mu);
        ++error;
      }

      if (eqs.viscosityModel.type == ViscosityModelData::SUTHERLAND && ref.temperature < 0.0) {
        com->fprintf(stderr, "*** Error: no valid reference temperature (%d) given\n", ref.temperature);
        ++error;
      }
    }

// Included (MB)
    ref.dRe_mudMach = 0.0;

    if (bc.inlet.density < 0.0)
      sa.densFlag = false;
    else
      sa.densFlag = true;

    if (bc.inlet.pressure < 0.0)
      sa.pressFlag = false;
    else
      sa.pressFlag = true;


    double gamma = eqs.fluidModel.gasModel.specificHeatRatio;

    double omega = eqs.fluidModel.jwlModel.omega;
    double A1    = eqs.fluidModel.jwlModel.A1;
    double A2    = eqs.fluidModel.jwlModel.A2;
    double R1    = eqs.fluidModel.jwlModel.R1;
    double R2    = eqs.fluidModel.jwlModel.R2;
    double rhoref= eqs.fluidModel.jwlModel.rhoref;

    double Prefwater = eqs.fluidModel.liquidModel.Prefwater;
    double k1water = eqs.fluidModel.liquidModel.k1water;
    double k2water = eqs.fluidModel.liquidModel.k2water;
    double RHOrefwater = eqs.fluidModel.liquidModel.RHOrefwater;
    double Pref, awater, bwater;
    if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
      Pref = -k1water/k2water;
      awater = (Prefwater + k1water/k2water)/pow(RHOrefwater, k2water);
      bwater = k2water;
    }

    // set up density
    if (bc.inlet.density < 0.0)
      bc.inlet.density = 1.0;

    // set up pressure
    if (bc.inlet.pressure < 0.0)
      if (eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
          eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS)
        if(ref.mach>0.0) {
//          bc.inlet.pressure = bc.inlet.pressure / (gamma * ref.mach * ref.mach * (bc.inlet.pressure + eqs.fluidModel.gasModel.pressureConstant));
          bc.inlet.pressure = bc.inlet.density / (gamma * ref.mach * ref.mach * (1.0 + eqs.fluidModel.gasModel.pressureConstant));
          eqs.fluidModel.gasModel.pressureConstant *= bc.inlet.pressure;
        }
        else
          com->fprintf(stderr, "*** Error: no valid Mach number for non-dimensional simulation\n");
      else if (eqs.fluidModel.fluid == FluidModelData::JWL)
        if(ref.mach>0.0){
          double frho  = A1*(1-omega*bc.inlet.density/(R1*rhoref))*exp(-R1*rhoref/bc.inlet.density) +
                         A2*(1-omega*bc.inlet.density/(R2*rhoref))*exp(-R2*rhoref/bc.inlet.density);
          double frhop = A1*(-omega/(R1*rhoref) + (1-omega*bc.inlet.density/(R1*rhoref))*R1*rhoref/(bc.inlet.density*bc.inlet.density))
                           *exp(-R1*rhoref/bc.inlet.density)
                       + A2*(-omega/(R2*rhoref) + (1-omega*bc.inlet.density/(R2*rhoref))*R2*rhoref/(bc.inlet.density*bc.inlet.density))
                           *exp(-R2*rhoref/bc.inlet.density);
          bc.inlet.pressure = bc.inlet.pressure / (ref.mach*ref.mach * ((omega+1.0)*bc.inlet.pressure -frho+bc.inlet.density*frhop));
        }
        else
          com->fprintf(stderr, "*** Error: no valid Mach number for non-dimensional simulation\n");
      else if(eqs.fluidModel.fluid == FluidModelData::LIQUID)
        bc.inlet.pressure = Prefwater/((Prefwater+k1water/k2water)*k2water*ref.mach*ref.mach);

    // set up temperature (for Tait)
    if (bc.inlet.temperature < 0.0 && eqs.fluidModel.fluid == FluidModelData::LIQUID){
      com->fprintf(stderr, "*** Error: no valid non-dimensionalized temperature (%f) given\n", bc.inlet.temperature);
      error ++;
    }
    if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
      eqs.fluidModel.liquidModel.Pref = Pref / (awater * bwater * pow(RHOrefwater, bwater) * pow(ref.mach, 2.0));
      eqs.fluidModel.liquidModel.alpha = 1.0/(bwater * ref.mach*ref.mach);
      eqs.fluidModel.liquidModel.beta = bwater;
    }
  }

  return error;

}

//------------------------------------------------------------------------------------
int IoData::checkInputValuesDimensional(map<int,SurfaceData*>& surfaceMap)
{

/* Non dimensionalization of two phase flows is somewhat tricky,
 * in the sense that we don't get the same state equation for
 * both phases, even if the two fluids are both perfect gases
 * or are both barotropic liquids.
 * For instance, the perfect gas laws (P = \rho R T) in the non dimensionalized
 * form does not write the same for both gases. One could think that
 * the first law being written with \gamma1, the second should write
 * exactly the same way with \gamma2 instead of \gamma1. What actually
 * must be done is replace \gamma1 by R2*(\gamma1-1)/R1 = Cv2*(\gamma2-1)/Cv1
 * in order to get the correct non dimensional perfect gas law for
 * the second gas!
 * However note that you don t need the perfect gas law written in
 * this form to do the computation. Using just P = (\gamma - 1.0) rho e
 * you recover the right physics, except everything where the temperature
 * is involved. (Watch out for Navier-Stokes equations, but Euler are fine)
 */

  int error = 0;

  double R = eqs.fluidModel.gasModel.idealGasConstant;
  double gamma = eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = eqs.fluidModel.gasModel.pressureConstant;
  double cpgas  = eqs.fluidModel.gasModel.specificHeatPressure;

  double Rjwl   = eqs.fluidModel.jwlModel.idealGasConstant;
  double omegajwl  = eqs.fluidModel.jwlModel.omega;
  double A1jwl     = eqs.fluidModel.jwlModel.A1;
  double A2jwl     = eqs.fluidModel.jwlModel.A2;
  double R1jwl     = eqs.fluidModel.jwlModel.R1;
  double R2jwl     = eqs.fluidModel.jwlModel.R2;
  double rhorefjwl = eqs.fluidModel.jwlModel.rhoref;

  double Cwater = eqs.fluidModel.liquidModel.specificHeat;
  double k1water = eqs.fluidModel.liquidModel.k1water;
  double k2water = eqs.fluidModel.liquidModel.k2water;
  double Prefwater = eqs.fluidModel.liquidModel.Prefwater;
  double RHOrefwater = eqs.fluidModel.liquidModel.RHOrefwater;
  double Pref, awater, bwater;
  if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
    Pref = -k1water/k2water;
    awater = (Prefwater + k1water/k2water)/pow(RHOrefwater, k2water);
    bwater = k2water;
  }

  if (problem.mode == ProblemData::DIMENSIONAL) {
    if (bc.inlet.pressure < 0.0)  {
      if(eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
         eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS || 
         eqs.fluidModel.fluid == FluidModelData::JWL ){
        com->fprintf(stderr, "*** Error: no valid inlet pressure (%f) given\n", bc.inlet.pressure);
        ++error;
      }else if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
        bc.inlet.pressure = Prefwater;
        com->fprintf(stderr, "*** Warning: inlet pressure set to reference pressure for Tait's EOS\n");
      }
    }
    if (bc.inlet.density < 0.0)
      if(eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
         eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS || 
         eqs.fluidModel.fluid == FluidModelData::JWL ){
        com->fprintf(stderr, "*** Error: no valid inlet density (%f) given\n", bc.inlet.density);
        ++error;
      }

    if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
      if (mf.problem == MultiFluidData::BUBBLE)
        bc.inlet.density = pow( (bc.inlet.pressure - Pref)/awater, 1.0/bwater);
    }

    if(bc.inlet.temperature < 0.0) {
      if(eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
         eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS)
        //bc.inlet.temperature = (bc.inlet.pressure + gamma*Pstiff)/(R*bc.inlet.density);
        bc.inlet.temperature = gamma/(gamma-1.0) * (bc.inlet.pressure + Pstiff)/(cpgas*bc.inlet.density);
      else if(eqs.fluidModel.fluid == FluidModelData::JWL){
        double frhoref = A1jwl*(1-omegajwl*bc.inlet.density/(R1jwl*rhorefjwl))*exp(-R1jwl*rhorefjwl/bc.inlet.density) +
                         A2jwl*(1-omegajwl*bc.inlet.density/(R2jwl*rhorefjwl))*exp(-R2jwl*rhorefjwl/bc.inlet.density);
        bc.inlet.temperature = (bc.inlet.pressure - frhoref)/(bc.inlet.density * Rjwl);
      }
      else if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
        com->fprintf(stderr, "*** Error: no valid inlet temperature (%f) given\n", bc.inlet.temperature);
        ++error;
      }
    }

    if(eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
       eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS){
      if (ref.density < 0.0)
        ref.density = bc.inlet.density;
      if (ref.pressure < 0.0)
        ref.pressure = bc.inlet.pressure;
      if (ref.mach <= 0.0)
        ref.mach = ref.velocity /sqrt(gamma * (ref.pressure+Pstiff) / ref.density);
      double velocity = ref.mach * sqrt(gamma * (ref.pressure+Pstiff) / ref.density);
      //ref.temperature = (ref.pressure + gamma*Pstiff)/ (ref.density * R);
      ref.temperature = gamma/(gamma-1.0) * (ref.pressure + Pstiff)/ (ref.density * cpgas);
      // if the gas is perfect, Pstiff = 0 and Cp = gamma*R/(gamma-1);
      ref.energy = (ref.pressure + gamma*Pstiff) / (ref.density * (gamma-1));
      double Cv = cpgas/gamma;
      double viscosity = 1.0;
      if (eqs.type == EquationsData::NAVIER_STOKES){
        viscosity = eqs.viscosityModel.sutherlandConstant * sqrt(ref.temperature) /
        (1.0 + eqs.viscosityModel.sutherlandReferenceTemperature/ref.temperature);
        if(eqs.viscosityModel.type == ViscosityModelData::CONSTANT)
	  viscosity = eqs.viscosityModel.dynamicViscosity;
        ref.reynolds_mu = velocity * ref.length * ref.density / viscosity;
      }

// Included (MB)
      double dvelocitydMach = sqrt(gamma * ref.pressure / ref.density);
      //if (eqs.type == EquationsData::NAVIER_STOKES)
        //com->fprintf(stderr, "\n\n Reynolds = %e \n\n",ref.reynolds_mu);
      ref.dRe_mudMach = dvelocitydMach * ref.length * ref.density / viscosity;
            
      ref.rv.mode = RefVal::DIMENSIONAL;
      ref.rv.density = ref.density;
      ref.rv.velocity = velocity;
      ref.rv.pressure = ref.density * velocity*velocity;
      ref.rv.temperature = velocity*velocity/Cv;
//      ref.rv.temperature = gamma*(gamma - 1.0) * ref.mach*ref.mach * ref.temperature;
      ref.rv.viscosity_mu = viscosity;
      ref.rv.nutilde = viscosity / ref.density;
      ref.rv.kenergy = velocity*velocity;
      ref.rv.epsilon = velocity*velocity*velocity / ref.length;
      ref.rv.time = ref.length / velocity; 
      ref.rv.force = ref.density * velocity*velocity * ref.length*ref.length;
      ref.rv.energy = ref.density * velocity*velocity * ref.length*ref.length*ref.length;
      ref.rv.power = ref.density * velocity*velocity*velocity * ref.length*ref.length;
      ref.rv.tvelocity = velocity / aero.displacementScaling;
      ref.rv.tforce = ref.rv.force / aero.forceScaling;
      ref.rv.tpower = ref.rv.power / aero.powerScaling;
      ref.rv.entropy = pow(ref.rv.density,1.0-gamma)*velocity*velocity;

// Included (MB)
      ref.rv.dvelocitydMach = dvelocitydMach;
      ref.rv.dtimedMach = - ref.length / (velocity * velocity) * dvelocitydMach;
      /*
      fprintf(stderr,"Reference State Calculated in IoDataCore.C\n");
      fprintf(stderr,"You are asking for a Dimensional Simulation for a Gas\n");
      fprintf(stderr,"Constants      : gamma   = %f, R        = %f, Pstiff   = %f, Length  = %f\n",gamma, R, Pstiff,ref.length);
      fprintf(stderr,"State          : density = %f, velocity = %f, pressure = %f, temp    = %f, viscosity = %f\n",
	      ref.density, velocity, ref.pressure, ref.temperature, viscosity);
      fprintf(stderr,"Flow Parameters: mach    = %f, reynolds = %f\n", ref.mach,ref.reynolds_mu);
      fprintf(stderr,"Ref Values     : length  = %f, density  = %f, velocity = %f, pressure = %f, temp = %f, viscosity = %f\n",
	      ref.length,ref.rv.density, ref.rv.velocity, ref.rv.pressure, ref.rv.temperature, ref.rv.viscosity_mu);
      fprintf(stderr,"\n");
      */
    }
    else if(eqs.fluidModel.fluid == FluidModelData::JWL){
      if (ref.density < 0.0)
        ref.density = bc.inlet.density;
      double frhoref = A1jwl*(1-omegajwl*ref.density/(R1jwl*rhorefjwl))*exp(-R1jwl*rhorefjwl/ref.density) +
                       A2jwl*(1-omegajwl*ref.density/(R2jwl*rhorefjwl))*exp(-R2jwl*rhorefjwl/ref.density);
      double frhorefp = A1jwl*(-omegajwl/(R1jwl*rhorefjwl) + (1.0-omegajwl*ref.density/(R1jwl*rhorefjwl))*R1jwl*rhorefjwl/(ref.density*ref.density))
                          *exp(-R1jwl*rhorefjwl/ref.density)
                      + A2jwl*(-omegajwl/(R2jwl*rhorefjwl) + (1.0-omegajwl*ref.density/(R2jwl*rhorefjwl))*R2jwl*rhorefjwl/(ref.density*ref.density))
                          *exp(-R2jwl*rhorefjwl/ref.density);
      if (ref.pressure < 0.0)
        ref.pressure = bc.inlet.pressure;
      if (ref.mach <= 0.0)
        ref.mach = ref.velocity /sqrt(((omegajwl+1.0)*ref.pressure - frhoref + ref.density*frhorefp)/ref.density);
      double velocity = ref.mach * sqrt(((omegajwl+1.0)*ref.pressure - frhoref + ref.density*frhorefp)/ref.density);
      ref.temperature = (ref.pressure - frhoref)/(ref.density * Rjwl);
      double viscosity = 0.0000000001;
      ref.reynolds_mu = velocity * ref.length * ref.density / viscosity;

      ref.rv.mode = RefVal::DIMENSIONAL;
      ref.rv.density = ref.density;
      ref.rv.velocity = velocity;
      ref.rv.pressure = ref.density * velocity*velocity;
      ref.rv.temperature = omegajwl*(omegajwl + 1.0) * ref.mach*ref.mach * ((omegajwl+1.0)*ref.pressure - frhoref + ref.density*frhorefp);
      ref.rv.viscosity_mu = viscosity;
      ref.rv.nutilde = viscosity / ref.density;
      ref.rv.kenergy = velocity*velocity;
      ref.rv.epsilon = velocity*velocity*velocity / ref.length;
      ref.rv.time = ref.length / velocity;
      ref.rv.force = ref.density * velocity*velocity * ref.length*ref.length;
      ref.rv.energy = ref.density * velocity*velocity * ref.length*ref.length*ref.length;
      ref.rv.power = ref.density * velocity*velocity*velocity * ref.length*ref.length;
      ref.rv.tvelocity = velocity / aero.displacementScaling;
      ref.rv.tforce = ref.rv.force / aero.forceScaling;
      ref.rv.tpower = ref.rv.power / aero.powerScaling;
      ref.rv.entropy = pow(ref.rv.density,-omegajwl)*velocity*velocity;

    }
    else if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
      if (ref.density < 0.0)
        ref.density = bc.inlet.density;
      if (ref.pressure < 0.0)
        ref.pressure = bc.inlet.pressure;
      if (ref.temperature < 0.0)
        ref.temperature = bc.inlet.temperature;
      if (ref.mach <= 0.0)
        ref.mach = ref.velocity / sqrt(bwater*awater*pow(ref.density, bwater - 1.0));
      double velocity = ref.mach * sqrt(bwater*awater*pow(ref.density, bwater - 1.0));
      double soundvelocity = sqrt(bwater*awater*pow(ref.density, bwater - 1.0));
//      Cwater has not been specified by user (not possible), we set it automatically using
//      EOS compatibility condition and inlet BCs.
      Cwater = awater*bwater*pow(bc.inlet.density, bwater-1.0)/((bwater-1.0)*bc.inlet.temperature);
      eqs.fluidModel.liquidModel.specificHeat = Cwater;
      ref.energy = Cwater * ref.temperature; // this is actually enthalpy, but for the sake of non-dimensionalization
      double Cv = Cwater;

      double viscosity = eqs.viscosityModel.dynamicViscosity;
      ref.reynolds_mu = velocity * ref.length * ref.density / viscosity;
      ref.rv.mode = RefVal::DIMENSIONAL;

      ref.rv.density = ref.density;
      ref.rv.velocity = velocity;
      ref.rv.pressure = ref.density * velocity*velocity;
      ref.rv.temperature = velocity * velocity / Cv;


      ref.rv.viscosity_mu = viscosity;
      ref.rv.nutilde = viscosity / ref.density;
      ref.rv.kenergy = velocity*velocity;
      ref.rv.epsilon = velocity*velocity*velocity / ref.length;
      ref.rv.time = ref.length / velocity;
      ref.rv.force = ref.density * velocity*velocity * ref.length*ref.length;
      ref.rv.energy = ref.density * velocity*velocity * ref.length*ref.length*ref.length;
      ref.rv.power = ref.density * velocity*velocity*velocity * ref.length*ref.length;
      ref.rv.entropy = pow(ref.rv.density,1.0-bwater)*velocity*velocity;

      ref.rv.tvelocity = velocity / aero.displacementScaling;

      ref.rv.tforce = ref.rv.force / aero.forceScaling;
      ref.rv.tpower = ref.rv.power / aero.powerScaling;

    }

    // non-dimensionalize boundary conditions
    bc.inlet.density /= ref.rv.density;
    bc.inlet.pressure /= ref.rv.pressure;
    bc.inlet.velocity /= ref.rv.velocity;
    bc.inlet.temperature /= ref.rv.temperature;
    bc.inlet.nutilde /= ref.rv.nutilde;
    bc.inlet.kenergy /= ref.rv.kenergy;
    bc.inlet.eps /= ref.rv.epsilon;
    bc.outlet.density /= ref.rv.density;
    bc.outlet.pressure /= ref.rv.pressure;
    bc.outlet.velocity /= ref.rv.velocity;
    bc.outlet.temperature /= ref.rv.temperature;
    bc.outlet.nutilde /= ref.rv.nutilde;
    bc.outlet.kenergy /= ref.rv.kenergy;
    bc.outlet.eps /= ref.rv.epsilon;

    restart.etime /= ref.rv.time;
    restart.dt_nm1 /= ref.rv.time;
    restart.dt_nm2 /= ref.rv.time;
    restart.energy /= ref.rv.energy;
    bc.wall.temperature /= ref.rv.temperature;

       for(int j = 1; j < 8*sizeof(int); j++) {
            map<int,SurfaceData*>::iterator it = surfaceMap.find(j);
             if(it == surfaceMap.end())
               continue;
             if(it->second->type == SurfaceData::ISOTHERMAL) {
               it->second->temp /= ref.rv.temperature;
            }
       }

    linearizedData.stepsize = ts.timestep;
    linearizedData.stepsizeinitial = ts.timestepinitial;
    ts.timestep /= ref.rv.time;
    ts.maxTime /= ref.rv.time;
    rmesh.vx /= ref.rv.velocity;
    rmesh.vy /= ref.rv.velocity;
    rmesh.vz /= ref.rv.velocity;
    rmesh.ax /= ref.rv.velocity / ref.rv.time;
    rmesh.ay /= ref.rv.velocity / ref.rv.time;
    rmesh.az /= ref.rv.velocity / ref.rv.time;
    rmesh.timestep /= ref.rv.time;

    for (int j=0; j<rmesh.num; j++){
//      com->fprintf(stderr,"old time: %f, new time: %f.\n", rmesh.vpts[j]->time, rmesh.vpts[j]->time/ref.rv.time);
      rmesh.vpts[j]->time     /= ref.rv.time;
      rmesh.vpts[j]->velocityX /= ref.rv.velocity;
      rmesh.vpts[j]->velocityY /= ref.rv.velocity;
      rmesh.vpts[j]->velocityZ /= ref.rv.velocity;
    }
    aero.pressure /= ref.rv.pressure;
    forced.timestep /= ref.rv.time;
    forced.frequency *= ref.rv.time;

    bc.hydro.depth /= ref.length;
    eqs.gravity_x /= ref.rv.velocity / ref.rv.time;
    eqs.gravity_y /= ref.rv.velocity / ref.rv.time;
    eqs.gravity_z /= ref.rv.velocity / ref.rv.time;
  }

  ref.rv.length = ref.length;
  ref.rv.tlength = ref.length / aero.displacementScaling;

  bc.wall.delta /= ref.rv.tlength;
  schemes.fixes.dihedralAngle *= acos(-1.0) / 180.0;
  for (int j=0; j<schemes.fixes.num; ++j) {
    schemes.fixes.spheres[j]->x0 /= ref.rv.tlength;
    schemes.fixes.spheres[j]->y0 /= ref.rv.tlength;
    schemes.fixes.spheres[j]->z0 /= ref.rv.tlength;
    schemes.fixes.spheres[j]->r /= ref.rv.tlength;
    schemes.fixes.boxes[j]->x0 /= ref.rv.tlength;
    schemes.fixes.boxes[j]->y0 /= ref.rv.tlength;
    schemes.fixes.boxes[j]->z0 /= ref.rv.tlength;
    schemes.fixes.boxes[j]->x1 /= ref.rv.tlength;
    schemes.fixes.boxes[j]->y1 /= ref.rv.tlength;
    schemes.fixes.boxes[j]->z1 /= ref.rv.tlength;

    schemes.fixes.cones[j]->x0 /= ref.rv.tlength;
    schemes.fixes.cones[j]->y0 /= ref.rv.tlength;
    schemes.fixes.cones[j]->z0 /= ref.rv.tlength;
    schemes.fixes.cones[j]->r0 /= ref.rv.tlength;
    schemes.fixes.cones[j]->x1 /= ref.rv.tlength;
    schemes.fixes.cones[j]->y1 /= ref.rv.tlength;
    schemes.fixes.cones[j]->z1 /= ref.rv.tlength;
    schemes.fixes.cones[j]->r1 /= ref.rv.tlength;
  }
  
  return error;
}
//------------------------------------------------------------------------------------
int IoData::checkInputValuesEssentialBC()
{
  int error = 0;

  if (ref.mach < 0.0)
    ref.mach = bc.inlet.mach;

  if (bc.inlet.mach < 0.0 && bc.inlet.velocity < 0.0)
    bc.inlet.mach = ref.mach;

  if (ref.velocity < 0.0)
    ref.velocity = bc.inlet.velocity;

  if (ref.mach <= 0.0) {
    com->fprintf(stderr, "*** Warning: no valid Mach number (%e) given\n", ref.mach);
    if (ref.velocity<0.0){
      com->fprintf(stderr, "*** Error: no valid Mach number and no valid velocity given\n");
      ++error;
    }
    else
      com->fprintf(stderr, "*** Warning: velocity used instead of mach number\n");
  }
  if (bc.inlet.alpha > 360.0) {
    com->fprintf(stderr, "*** Error: no valid angle of attack (%e) given\n", bc.inlet.alpha);
    ++error;
  }
  if (bc.inlet.beta > 360.0) {
    com->fprintf(stderr, "*** Error: no valid yaw angle (%e) given\n", bc.inlet.beta);
    ++error;
  }

  if(schemes.bc.type==BoundarySchemeData::MODIFIED_GHIDAGLIA && ts.type!=TsData::EXPLICIT) {
    com->fprintf(stderr, "*** Warning: The Modified Ghidaglia scheme is only supported by explicit time-integrators.\n");
    com->fprintf(stderr, "             Reset to the standard Ghidaglia scheme.\n");
    schemes.bc.type = BoundarySchemeData::GHIDAGLIA;
  }

// Included (MB)
  if (sa.machref < 0.0)
    sa.machref = bc.inlet.mach;
  if (sa.alpharef > 360.0)
    sa.alpharef = bc.inlet.alpha;
  if (sa.betaref > 360.0)
    sa.betaref = bc.inlet.beta;
  if (!sa.angleRad) {
    sa.alpharef *= acos(-1.0) / 180.0;
    sa.betaref *= acos(-1.0) / 180.0;
  }

  return error;

}
//------------------------------------------------------------------------------
// must be done after the non-dimensionalization of input values!
void IoData::checkInputValuesTurbulence()
{
   if (bc.inlet.nutilde < 0.0)
    bc.inlet.nutilde = 0.1;
  double theta_k = 1.0;
  double theta_w = 10.0;
  if (bc.inlet.kenergy < 0.0)
    bc.inlet.kenergy = pow(10.0, -theta_k) * theta_w / ref.reynolds_mu;
  if (bc.inlet.eps < 0.0)
    bc.inlet.eps = eqs.tc.tm.ke.c_mu * bc.inlet.kenergy * theta_w;

}

//------------------------------------------------------------------------------

void IoData::checkInputValuesDefaultOutlet()
{
  if (bc.outlet.mach < 0.0)
    bc.outlet.mach = bc.inlet.mach;
  if (bc.outlet.velocity < 0.0)
    bc.outlet.velocity = bc.inlet.velocity;
  if (bc.outlet.density < 0.0)
    bc.outlet.density = bc.inlet.density;
  if (bc.outlet.pressure < 0.0)
    bc.outlet.pressure = bc.inlet.pressure;
  if (bc.outlet.temperature < 0.0)
    bc.outlet.temperature = bc.inlet.temperature;
  if (bc.outlet.nutilde < 0.0)
    bc.outlet.nutilde = bc.inlet.nutilde;
  if (bc.outlet.kenergy < 0.0)
    bc.outlet.kenergy = bc.inlet.kenergy;
  if (bc.outlet.eps < 0.0)
    bc.outlet.eps = bc.inlet.eps;
  if (bc.outlet.alpha > 360.0)
    bc.outlet.alpha = bc.inlet.alpha;
  if (bc.outlet.beta > 360.0)
    bc.outlet.beta = bc.inlet.beta;

/*
// TDL: This is a bug... already done in :checkInputValues()
// Modified (MB)
  if (!sa.angleRad) {
    bc.inlet.alpha *= acos(-1.0) / 180.0;
    bc.inlet.beta *= acos(-1.0) / 180.0;
    bc.outlet.alpha *= acos(-1.0) / 180.0;
    bc.outlet.beta *= acos(-1.0) / 180.0;
  }
*/

// Included (MB)
  if (aero.pressure < 0.0)
    sa.apressFlag = false;
  else
    sa.apressFlag = true;

}

//------------------------------------------------------------------------------

int IoData::checkSolverValues(map<int,SurfaceData*>& surfaceMap)
{

  int error = 0;

  if (problem.type[ProblemData::ACCELERATED] && !problem.type[ProblemData::AERO] &&
      rmesh.timestep < 0.0) {
    com->fprintf(stderr, "*** Error: no valid timestep (%d) given\n", rmesh.timestep);
    ++error;
  }
  if (problem.type[ProblemData::FORCED]) {
    if (forced.timestep < 0.0) {
      com->fprintf(stderr, "*** Error: no valid timestep (%d) given\n", forced.timestep);
      ++error;
    }
    if (forced.frequency < 0.0) {
      com->fprintf(stderr, "*** Error: no valid frequency (%d) given\n", forced.frequency);
      ++error;
    }
  }

  if (eqs.type != EquationsData::EULER && bc.wall.type == BcsWallData::ISOTHERMAL &&
      !problem.type[ProblemData::THERMO] && bc.wall.temperature < 0.0) {
    com->fprintf(stderr, "*** Error: no valid wall temperature (%f) given\n", bc.wall.temperature);
    ++error;
  }


for(int j = 1; j < 8*sizeof(int); j++) {
            map<int,SurfaceData*>::iterator it = surfaceMap.find(j);
             if(it == surfaceMap.end())
               continue;
             if(it->second->type == SurfaceData::ISOTHERMAL && it->second->temp < 0 && eqs.type != EquationsData::EULER && bc.wall.type == BcsWallData::ADIABATIC &&
                 !problem.type[ProblemData::THERMO] && bc.wall.temperature < 0.0) {
               com->fprintf(stderr, "*** Error: no valid wall temperature (%f) given and (%f) given for surface (%d) \n", bc.wall.temperature, it->second->temp, j);
               error++;
            }
       }


// for Multiphase flow using levelset
  if(eqs.numPhase > 1 && schemes.ls.reconstruction == SchemeData::CONSTANT
                      && mf.interfaceType != MultiFluidData::FSF){
    com->fprintf(stderr, "*** Error: Linear reconstruction of the levelset is needed!\n");
    ++error;
  }

  return error;

}

//------------------------------------------------------------------------------

int IoData::checkInputValuesInitialConditions(InitialConditions &initialConditions,
                                              int fluidModelID)
{

  int error = 0;

  // first check that there is a corresponding fluidModelID in eqs.fluidModelMap
  // and determine the type of that fluidModel
  // note that eqs.fluidModel does not need to be considered because it can only be used for single phase flow
  // and this routine is for multiphase flow.
  int fluidModelIDCount = 0;
  FluidModelData::Fluid fluidType = FluidModelData::UNDEFINED;
  if(eqs.fluidModelMap.dataMap.empty()){
    error++;
    com->fprintf(stderr, "*** Error: no valid FluidModel were specified for Initial Conditions\n");
    com->fprintf(stderr, "           FluidModelData[0] must at least be specified\n");
  }else{
    map<int, FluidModelData *>::iterator it;
    for (it=eqs.fluidModelMap.dataMap.begin(); it!=eqs.fluidModelMap.dataMap.end();it++)
      if(it->first == fluidModelID){
        fluidModelIDCount++;
        fluidType = it->second->fluid;
      }
    if(fluidModelIDCount != 1){
      error++;
      com->fprintf(stderr, "*** Error: there are %d fluidModel(s) %d\n", fluidModelIDCount, fluidModelID);
    }
  }
      

  // then check that the initial conditions for that fluidModel are adequate.
  if(initialConditions.mach < 0 && initialConditions.velocity < 0 ){
    error++;
    com->fprintf(stderr, "*** Error : an initial velocity norm or an initial mach number must be specified\n");
  }

  if(fluidType == FluidModelData::PERFECT_GAS || 
     fluidType == FluidModelData::STIFFENED_GAS ||
     fluidType == FluidModelData::JWL){
    if(initialConditions.density < 0){
      error++;
      com->fprintf(stderr, "*** Error : an initial density must be specified\n");
    }
    if(initialConditions.pressure < 0){
      error++;
      com->fprintf(stderr, "*** Error : an initial pressure must be specified\n");
    }
  }
  else if(fluidType == FluidModelData::LIQUID){
    if(initialConditions.pressure < 0 && initialConditions.density < 0 ){
      error++;
      com->fprintf(stderr, "*** Error : either initial pressure or density must be specified\n");
    }
    if(initialConditions.temperature < 0){
      error++;
      com->fprintf(stderr, "*** Error : an initial temperature must be specified\n");
    }
  }

  return error;

}

//------------------------------------------------------------------------------

int IoData::checkInputValuesEquationOfState(FluidModelData &fluidModel, int fluidModelID)
{
  int error = 0;

  // ******* TAIT ******** //
  if (fluidModel.fluid == FluidModelData::LIQUID){
    if (fluidModel.liquidModel.k2water <= 0.0){
      com->fprintf(stderr, "*** Error: k2 in the barotropic liquid EOS cannot be zero or negative\n");
      ++error;
    }
    if (fluidModel.liquidModel.Prefwater < 0.0){
      com->fprintf(stderr, "*** Error: no valid reference pressure (%e) given for Tait's EOS (fluidModelID = %d)\n", fluidModel.liquidModel.Prefwater, fluidModelID);
      ++error;
    }
    if (fluidModel.liquidModel.RHOrefwater < 0.0){
      com->fprintf(stderr, "*** Error: no valid reference density (%e) given for Tait's EOS (fluidModelID = %d)\n", fluidModel.liquidModel.RHOrefwater, fluidModelID);
      ++error;
    }
    // eqs.fluidModel.liquidModel.specificHeat is determined from the inlet temperature and density and the other
    // constants of the EOS.
  }

  // ******* PERFECT GAS ******** //
  else if (fluidModel.fluid         == FluidModelData::PERFECT_GAS){
    fluidModel.gasModel.pressureConstant = 0.0;
    double gamma = fluidModel.gasModel.specificHeatRatio;
    fluidModel.gasModel.specificHeatPressure = gamma*fluidModel.gasModel.idealGasConstant/(gamma-1.0);
  }

  // ******* STIFFENED GAS ******** //
  else if (fluidModel.fluid         == FluidModelData::STIFFENED_GAS){
    if (fluidModel.gasModel.pressureConstant < 0.0) {
      com->fprintf(stderr, "*** Error: no valid reference pressure constant (%e) given for Stiffened Gas EOS (fluidModelID = %d)\n", fluidModel.gasModel.pressureConstant, fluidModelID);
      ++error;
    }
    if(fluidModel.gasModel.specificHeatPressure < 0.0){
      if(eqs.type == EquationsData::NAVIER_STOKES) {
        com->fprintf(stderr, "*** Error: a specific heat at constant pressure must be specified for a stiffened gas in a viscous simulation.");
        error++;
      }
      else if(eqs.type == EquationsData::EULER)
        com->fprintf(stderr, "*** Note: a specific heat at constant pressure must be specified for a stiffened gas for post-processing temperature.\n");
    }
  }

  // ******* JWL GAS ******** //
  else if (fluidModel.fluid         == FluidModelData::JWL &&
	   fluidModel.jwlModel.type == JWLModelData::IDEAL){
    fluidModel.jwlModel.A1 = 0.0;
    fluidModel.jwlModel.A2 = 0.0;
  }

  // ******* JWL GAS ******** //
  else if (fluidModel.fluid == FluidModelData::JWL){
    if (fluidModel.jwlModel.R1*fluidModel.jwlModel.rhoref < 0.0){
      com->fprintf(stderr, "*** Error: negative value of R1*rhoref (%e) given for JWL EOS (fluidModelID = %d)\n", fluidModel.jwlModel.R1*fluidModel.jwlModel.rhoref, fluidModelID);
      ++error;
    }
    if (fluidModel.jwlModel.R2*fluidModel.jwlModel.rhoref < 0.0){
      com->fprintf(stderr, "*** Error: negative value of R2*rhoref (%e) given for JWL EOS (fluidModelID = %d)\n", fluidModel.jwlModel.R2*fluidModel.jwlModel.rhoref, fluidModelID);
      ++error;
    }
  }

  return error;
}

//------------------------------------------------------------------------------

void IoData::nonDimensionalizeInitialConditions(InitialConditions &initialConditions){

  initialConditions.velocity    /= ref.rv.velocity;
  initialConditions.pressure    /= ref.rv.pressure;
  initialConditions.density     /= ref.rv.density;
  initialConditions.temperature /= ref.rv.temperature;

  initialConditions.alpha       *= acos(-1.0)/180.0;
  initialConditions.beta        *= acos(-1.0)/180.0;

}

//------------------------------------------------------------------------------

void IoData::nonDimensionalizeFluidModel(FluidModelData &fluidModel){

  fluidModel.rhomin /= ref.rv.density;
  fluidModel.pmin /= ref.rv.pressure;

  if(fluidModel.fluid == FluidModelData::PERFECT_GAS ||
     fluidModel.fluid == FluidModelData::STIFFENED_GAS){
    fluidModel.gasModel.pressureConstant /= ref.rv.pressure;
    fluidModel.gasModel.specificHeatPressure /= ref.rv.velocity*ref.rv.velocity/ref.rv.temperature;
  }

  else if(fluidModel.fluid == FluidModelData::JWL){
    fluidModel.jwlModel.A1     /= ref.rv.pressure;
    fluidModel.jwlModel.A2     /= ref.rv.pressure;
    fluidModel.jwlModel.rhoref /= ref.rv.density;
  }

  else if(fluidModel.fluid == FluidModelData::LIQUID){
    double Pref = -fluidModel.liquidModel.k1water/fluidModel.liquidModel.k2water;
    double awater = (fluidModel.liquidModel.Prefwater - Pref)/pow(fluidModel.liquidModel.RHOrefwater, fluidModel.liquidModel.k2water);
    double bwater = fluidModel.liquidModel.k2water;

    fluidModel.liquidModel.specificHeat /= ref.rv.velocity*ref.rv.velocity/ref.rv.temperature;
    fluidModel.liquidModel.k1water     /= ref.rv.pressure;
    fluidModel.liquidModel.RHOrefwater /= ref.rv.density;
    fluidModel.liquidModel.Prefwater   /= ref.rv.pressure;

    fluidModel.liquidModel.Pref  = Pref / ref.rv.pressure;
    fluidModel.liquidModel.alpha = awater * pow(ref.rv.density, bwater - 1.0)/(ref.rv.velocity *ref.rv.velocity);
    fluidModel.liquidModel.beta  = bwater;
  }

}

//------------------------------------------------------------------------------
 
void IoData::nonDimensionalizeViscosityModel(ViscosityModelData &vm){

  vm.bulkViscosity /= ref.rv.viscosity_mu;

  //others do not need and must not be non-dimensionalized

}

//------------------------------------------------------------------------------

void IoData::nonDimensionalizeThermalCondModel(ThermalCondModelData &tm){

  double Cv = ref.rv.velocity*ref.rv.velocity/ref.rv.temperature;
  tm.conductivity /= (ref.rv.viscosity_mu*Cv);
}

//------------------------------------------------------------------------------

int IoData::checkInputValuesSparseGrid(SparseGridData &sparseGrid){

  int error = 0;

  if(sparseGrid.minPoints > sparseGrid.maxPoints){
    com->fprintf(stderr, "*** Error: sparse grid has incorrect number of points\n");
    error++;
  }

  if(sparseGrid.relAccuracy <= 0.0 || sparseGrid.absAccuracy <= 0.0){
    com->fprintf(stderr, "*** Error: sparse grid has incorrect accuracy specification(s)\n");
    error++;
  }

  if(sparseGrid.dimAdaptDegree<0.0 || sparseGrid.dimAdaptDegree>1.0){
    com->fprintf(stderr, "*** Error: sparse grid must have a dimension adaptivity degree between 0 and 1\n");
    error++;
  }

  if(sparseGrid.numOutputs <= 0 || sparseGrid.numInputs <= 0){
    com->fprintf(stderr, "*** Error: sparse grid must have positive numbers of inputs and outputs\n");
    error++;
  }

  typedef double Range[2];
  sparseGrid.range = new Range[sparseGrid.numInputs];
  sparseGrid.mapBaseValue = new double[sparseGrid.numInputs];
  sparseGrid.numDomainDim = new int[sparseGrid.numInputs];
  if(sparseGrid.numInputs > 5){
    sparseGrid.range[5][0] = sparseGrid.range6min;
    sparseGrid.range[5][1] = sparseGrid.range6max;
    sparseGrid.mapBaseValue[5] = sparseGrid.mapBaseValue6;
    sparseGrid.numDomainDim[5] = sparseGrid.numDomainDim6;
  }
  if(sparseGrid.numInputs > 4){
    sparseGrid.range[4][0] = sparseGrid.range5min;
    sparseGrid.range[4][1] = sparseGrid.range5max;
    sparseGrid.mapBaseValue[4] = sparseGrid.mapBaseValue5;
    sparseGrid.numDomainDim[4] = sparseGrid.numDomainDim5;
  }
  if(sparseGrid.numInputs > 3){
    sparseGrid.range[3][0] = sparseGrid.range4min;
    sparseGrid.range[3][1] = sparseGrid.range4max;
    sparseGrid.mapBaseValue[3] = sparseGrid.mapBaseValue4;
    sparseGrid.numDomainDim[3] = sparseGrid.numDomainDim4;
  }
  if(sparseGrid.numInputs > 2){
    sparseGrid.range[2][0] = sparseGrid.range3min;
    sparseGrid.range[2][1] = sparseGrid.range3max;
    sparseGrid.mapBaseValue[2] = sparseGrid.mapBaseValue3;
    sparseGrid.numDomainDim[2] = sparseGrid.numDomainDim3;
  }
  if(sparseGrid.numInputs > 1){
    sparseGrid.range[1][0] = sparseGrid.range2min;
    sparseGrid.range[1][1] = sparseGrid.range2max;
    sparseGrid.mapBaseValue[1] = sparseGrid.mapBaseValue2;
    sparseGrid.numDomainDim[1] = sparseGrid.numDomainDim2;
  }
  if(sparseGrid.numInputs > 0){
    sparseGrid.range[0][0] = sparseGrid.range1min;
    sparseGrid.range[0][1] = sparseGrid.range1max;
    sparseGrid.mapBaseValue[0] = sparseGrid.mapBaseValue1;
    sparseGrid.numDomainDim[0] = sparseGrid.numDomainDim1;
  }

  for(int i=0; i<sparseGrid.numInputs; i++){
    if(sparseGrid.range[i][0] > sparseGrid.range[i][1]){
      com->fprintf(stderr, "*** Error: sparse grid must have increasing range for input %d\n", i);
      error++;
    }
  }

  return error;

}

 int IoData::checkProgrammedBurnLocal(ProgrammedBurnData& programmedBurn,
				      InitialConditions& IC) {

  int error = 0;
  if (programmedBurn.unburnedEOS < 0)
    return 0;
  
  if (eqs.fluidModelMap.dataMap.find(programmedBurn.burnedEOS) == eqs.fluidModelMap.dataMap.end()) {
    com->fprintf(stderr, "*** Error: Cannot find burned EOS %d in fluid models\n",programmedBurn.burnedEOS);
    ++error;
  }
  
  if (programmedBurn.e0 <= 0.0) {
    com->fprintf(stderr, "*** Error: Champan-Jouguet Energy must be greater than zero");
    ++error;
  }
  
  // Check the detonation velocity/detonation pressure
  if (programmedBurn.cjDetonationVelocity <= 0.0 ||
      programmedBurn.cjPressure <= 0.0 ||
      programmedBurn.cjDensity <= 0.0 ||
      programmedBurn.cjEnergy <= 0.0) {
    
    com->fprintf(stderr, "*** At least one of the Chapman-Jouguet parameters is invalid\nRecalcuating...\n");
    const FluidModelData& burnedData = *eqs.fluidModelMap.dataMap.find(programmedBurn.burnedEOS)->second;
    if (burnedData.fluid == FluidModelData::JWL) {
      ProgrammedBurn::computeChapmanJouguetStateJWL(burnedData.jwlModel.A1, burnedData.jwlModel.A2,
						    burnedData.jwlModel.R1*burnedData.jwlModel.rhoref,
						    burnedData.jwlModel.R2*burnedData.jwlModel.rhoref,
						    burnedData.jwlModel.omega,IC.pressure,
						    IC.density,
						    programmedBurn.e0,
						    programmedBurn.cjDensity, programmedBurn.cjPressure,
						    programmedBurn.cjEnergy,programmedBurn.cjDetonationVelocity);
      
      com->fprintf(stderr,"Computed CJ state for JWL gas as: p_cj = %e rho_cj = %e e_cj = %e; Detonation Velocity = %e\n",
		   programmedBurn.burnedEOS, programmedBurn.cjPressure, programmedBurn.cjDensity, programmedBurn.cjEnergy,
		   programmedBurn.cjDetonationVelocity);
      
    } else if (burnedData.fluid == FluidModelData::PERFECT_GAS ||
               burnedData.fluid == FluidModelData::STIFFENED_GAS){
      if (burnedData.gasModel.pressureConstant != 0.0) {
	com->fprintf(stderr,"*** Warning: Correct CJ state only computed with perfect gas.  You have specified a stiffened gas for the burned state\n");
      }
      
      ProgrammedBurn::computeChapmanJouguetStatePG(burnedData.gasModel.specificHeatRatio,IC.pressure,
						   IC.density,
						   programmedBurn.e0,
						   programmedBurn.cjDensity, programmedBurn.cjPressure,
						   programmedBurn.cjEnergy,programmedBurn.cjDetonationVelocity);
      
      com->fprintf(stderr,"Computed CJ state for PG gas as: p_cj = %e rho_cj = %e e_cj = %e; Detonation Velocity = %e\n",
		   programmedBurn.burnedEOS, programmedBurn.cjPressure, programmedBurn.cjDensity,programmedBurn.cjEnergy,
		   programmedBurn.cjDetonationVelocity);
      
    }	
    
  }
  
  if (eqs.fluidModelMap.dataMap.find(programmedBurn.unburnedEOS) == eqs.fluidModelMap.dataMap.end()) {
    com->fprintf(stderr, "*** Error: Cannot find unburned EOS %d in fluid models\n",programmedBurn.unburnedEOS);
    ++error;
  } else {
    
    LiquidModelData& liq = eqs.fluidModelMap.dataMap.find(programmedBurn.unburnedEOS)->second->liquidModel;
    double B = (programmedBurn.cjPressure - liq.Prefwater) / ( pow(programmedBurn.cjDensity/liq.RHOrefwater, liq.k2water) - 1.0);
    B *= programmedBurn.factorB;
    //std::cout << "Old k1 = " << liq.k1water  << std::endl;
    liq.k1water = liq.k2water*(B - liq.Prefwater);       
    //std::cout << "New k1 = " << liq.k1water  << std::endl;
    
  }

  programmedBurn.cjDetonationVelocity *= programmedBurn.factorS;
  
  LiquidModelData& liq = eqs.fluidModelMap.dataMap.find(programmedBurn.unburnedEOS)->second->liquidModel;
  // Set the initial energy of the Tait EOS appropriately
  if (liq.burnable == LiquidModelData::YES)
    IC.temperature = programmedBurn.e0 / liq.specificHeat;//programmedBurn.cjEnergy / liq.Cv;
  else
    IC.temperature = (programmedBurn.e0 + IC.pressure/IC.density) / liq.specificHeat;
  //IC.temperature = programmedBurn.cjEnergy / liq.Cv;//programmedBurn.cjEnergy / liq.Cv;
  //std::cout << "T = " << IC.temperature << std::endl;
  return error;
}

//------------------------------------------------------------------------------
int IoData::checkInputValuesProgrammedBurn() {
  int error = 0;
  
  if(!mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator it;
    for (it=mf.multiInitialConditions.sphereMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.sphereMap.dataMap.end();
         it++){

      ProgrammedBurnData& programmedBurn = it->second->programmedBurn;
      InitialConditions& IC = it->second->initialConditions;
      error += checkProgrammedBurnLocal(programmedBurn,
					IC);
    }
  }

  if(!mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator itp;
    for (itp=mf.multiInitialConditions.prismMap.dataMap.begin(); 
         itp!=mf.multiInitialConditions.prismMap.dataMap.end();
         itp++){

      ProgrammedBurnData& programmedBurn = itp->second->programmedBurn;
      InitialConditions& IC = itp->second->initialConditions;
      error += checkProgrammedBurnLocal(programmedBurn,
					IC);
    }
  }

  

  return error;
}
//------------------------------------------------------------------------------

int IoData::checkInputValuesEmbeddedFramework() {
  int error = 0;

  if(!mf.multiInitialConditions.planeMap.dataMap.empty() ||
     !volumes.volumeMap.dataMap.empty()) {
    com->fprintf(stderr,"ERROR: Currently specifying initial conditions using 'planes' or 'volumes' are not supported by the Embedded Framework!\n");
    error ++;
  }
  return error;
}

//------------------------------------------------------------------------------

void IoData::printDebug(){

  com->fprintf(stderr, "to non-dimensionalize, use %e %e %e\n\n", ref.rv.density, ref.rv.velocity, ref.rv.pressure);

  // fluidModels
  if(!eqs.fluidModelMap.dataMap.empty()){
    map<int, FluidModelData *>::iterator it;
    for (it=eqs.fluidModelMap.dataMap.begin(); it!=eqs.fluidModelMap.dataMap.end();it++){
      com->fprintf(stderr, "FluidModelData::tag          = %d\n", it->first);
      com->fprintf(stderr, "FluidModelData::fluid        = %d\n", it->second->fluid);
      com->fprintf(stderr, "FluidModelData::rhomin         = %e\n", it->second->rhomin);
      com->fprintf(stderr, "FluidModelData::pmin         = %e\n", it->second->pmin);
      if(it->second->fluid == FluidModelData::PERFECT_GAS ||
         it->second->fluid == FluidModelData::STIFFENED_GAS){
        com->fprintf(stderr, "GasModelData::type              = %d\n", it->second->gasModel.type);
        com->fprintf(stderr, "GasModelData::specificHeatRatio = %e\n", it->second->gasModel.specificHeatRatio);
        com->fprintf(stderr, "GasModelData::idealGasConstant  = %e\n", it->second->gasModel.idealGasConstant);
        com->fprintf(stderr, "GasModelData::pressureConstant  = %e\n", it->second->gasModel.pressureConstant);
        com->fprintf(stderr, "GasModelData::specificHeatPressure  = %e\n", it->second->gasModel.specificHeatPressure);
      }
      else if(it->second->fluid == FluidModelData::JWL){
        com->fprintf(stderr, "JwlModelData::type              = %d\n", it->second->jwlModel.type);
        com->fprintf(stderr, "JwlModelData::omega             = %e\n", it->second->jwlModel.omega);
        com->fprintf(stderr, "JwlModelData::idealGasConstant  = %e\n", it->second->jwlModel.idealGasConstant);
        com->fprintf(stderr, "JwlModelData::A1                = %e\n", it->second->jwlModel.A1);
        com->fprintf(stderr, "JwlModelData::A2                = %e\n", it->second->jwlModel.A2);
        com->fprintf(stderr, "JwlModelData::R1                = %e\n", it->second->jwlModel.R1);
        com->fprintf(stderr, "JwlModelData::R2                = %e\n", it->second->jwlModel.R2);
        com->fprintf(stderr, "JwlModelData::rhoref            = %e\n", it->second->jwlModel.rhoref);
      }
      else if(it->second->fluid == FluidModelData::LIQUID){
        com->fprintf(stderr, "LiquidModelData::type           = %d\n", it->second->liquidModel.type);
        com->fprintf(stderr, "LiquidModelData::check          = %d\n", it->second->liquidModel.check);
        com->fprintf(stderr, "LiquidModelData::specificHeat   = %e\n", it->second->liquidModel.specificHeat);
        com->fprintf(stderr, "LiquidModelData::k1water        = %e\n", it->second->liquidModel.k1water);
        com->fprintf(stderr, "LiquidModelData::k2water        = %e\n", it->second->liquidModel.k2water);
        com->fprintf(stderr, "LiquidModelData::Prefwater      = %e\n", it->second->liquidModel.Prefwater);
        com->fprintf(stderr, "LiquidModelData::RHOrefwater    = %e\n", it->second->liquidModel.RHOrefwater);
        com->fprintf(stderr, "LiquidModelData::Pref           = %e\n", it->second->liquidModel.Pref);
        com->fprintf(stderr, "LiquidModelData::alpha          = %e\n", it->second->liquidModel.alpha);
        com->fprintf(stderr, "LiquidModelData::beta           = %e\n", it->second->liquidModel.beta);
      }
      com->fprintf(stderr, "\n");
    }
  }

  // initial conditions for volumeData
  if(!volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator it;
    for (it=volumes.volumeMap.dataMap.begin(); it!=volumes.volumeMap.dataMap.end();it++)
      if(it->second->type==VolumeData::FLUID){
        com->fprintf(stderr, "VolumeData::tag          = %d\n", it->first);
        com->fprintf(stderr, "VolumeData::fluidModelID = %d\n", it->second->fluidModelID);
        com->fprintf(stderr, "VolumeData::density      = %e\n", it->second->initialConditions.density);
        com->fprintf(stderr, "VolumeData::pressure     = %e\n", it->second->initialConditions.pressure);
        com->fprintf(stderr, "VolumeData::temperature  = %e\n", it->second->initialConditions.temperature);
        com->fprintf(stderr, "VolumeData::mach         = %e\n", it->second->initialConditions.mach);
        com->fprintf(stderr, "VolumeData::velocity     = %e\n", it->second->initialConditions.velocity);
        com->fprintf(stderr, "VolumeData::alpha        = %e\n", it->second->initialConditions.alpha);
        com->fprintf(stderr, "VolumeData::beta         = %e\n", it->second->initialConditions.beta);
        com->fprintf(stderr, "\n");
      }
  }

  // initial conditions for multiFluidData
  if(!mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator it;
    for (it=mf.multiInitialConditions.sphereMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.sphereMap.dataMap.end();
         it++){
      com->fprintf(stderr, "SphereData::tag          = %d\n", it->first);
      com->fprintf(stderr, "SphereData::fluidModelID = %d\n", it->second->fluidModelID);
      com->fprintf(stderr, "SphereData::density      = %e\n", it->second->initialConditions.density);
      com->fprintf(stderr, "SphereData::pressure     = %e\n", it->second->initialConditions.pressure);
      com->fprintf(stderr, "SphereData::temperature  = %e\n", it->second->initialConditions.temperature);
      com->fprintf(stderr, "SphereData::mach         = %e\n", it->second->initialConditions.mach);
      com->fprintf(stderr, "SphereData::velocity     = %e\n", it->second->initialConditions.velocity);
      com->fprintf(stderr, "SphereData::alpha        = %e\n", it->second->initialConditions.alpha);
      com->fprintf(stderr, "SphereData::beta         = %e\n", it->second->initialConditions.beta);
      com->fprintf(stderr, "SphereData::x0           = %e\n", it->second->cen_x);
      com->fprintf(stderr, "SphereData::y0           = %e\n", it->second->cen_y);
      com->fprintf(stderr, "SphereData::z0           = %e\n", it->second->cen_z);
      com->fprintf(stderr, "SphereData::R0           = %e\n", it->second->radius);
      com->fprintf(stderr, "\n");
    }
  }  
  if(!mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator it;
    for (it=mf.multiInitialConditions.prismMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.prismMap.dataMap.end();
         it++){
      com->fprintf(stderr, "PrismData::tag          = %d\n", it->first);
      com->fprintf(stderr, "PrismData::fluidModelID = %d\n", it->second->fluidModelID);
      com->fprintf(stderr, "PrismData::density      = %e\n", it->second->initialConditions.density);
      com->fprintf(stderr, "PrismData::pressure     = %e\n", it->second->initialConditions.pressure);
      com->fprintf(stderr, "PrismData::temperature  = %e\n", it->second->initialConditions.temperature);
      com->fprintf(stderr, "PrismData::mach         = %e\n", it->second->initialConditions.mach);
      com->fprintf(stderr, "PrismData::velocity     = %e\n", it->second->initialConditions.velocity);
      com->fprintf(stderr, "PrismData::alpha        = %e\n", it->second->initialConditions.alpha);
      com->fprintf(stderr, "PrismData::beta         = %e\n", it->second->initialConditions.beta);
      com->fprintf(stderr, "PrismData::x0           = %e\n", it->second->cen_x);
      com->fprintf(stderr, "PrismData::y0           = %e\n", it->second->cen_y);
      com->fprintf(stderr, "PrismData::z0           = %e\n", it->second->cen_z);
      com->fprintf(stderr, "PrismData::wx           = %e\n", it->second->w_x);
      com->fprintf(stderr, "PrismData::wy           = %e\n", it->second->w_y);
      com->fprintf(stderr, "PrismData::wz           = %e\n", it->second->w_z);
      com->fprintf(stderr, "\n");
    }
  }
  if(!mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator it;
    for (it=mf.multiInitialConditions.planeMap.dataMap.begin(); 
         it!=mf.multiInitialConditions.planeMap.dataMap.end();
         it++){
      com->fprintf(stderr, "PlaneData::tag          = %d\n", it->first);
      com->fprintf(stderr, "PlaneData::fluidModelID = %d\n", it->second->fluidModelID);
      com->fprintf(stderr, "PlaneData::density      = %e\n", it->second->initialConditions.density);
      com->fprintf(stderr, "PlaneData::pressure     = %e\n", it->second->initialConditions.pressure);
      com->fprintf(stderr, "PlaneData::temperature  = %e\n", it->second->initialConditions.temperature);
      com->fprintf(stderr, "PlaneData::mach         = %e\n", it->second->initialConditions.mach);
      com->fprintf(stderr, "PlaneData::velocity     = %e\n", it->second->initialConditions.velocity);
      com->fprintf(stderr, "PlaneData::alpha        = %e\n", it->second->initialConditions.alpha);
      com->fprintf(stderr, "PlaneData::beta         = %e\n", it->second->initialConditions.beta);
      com->fprintf(stderr, "PlaneData::x0           = %e\n", it->second->cen_x);
      com->fprintf(stderr, "PlaneData::y0           = %e\n", it->second->cen_y);
      com->fprintf(stderr, "PlaneData::z0           = %e\n", it->second->cen_z);
      com->fprintf(stderr, "PlaneData::nx           = %e\n", it->second->nx);
      com->fprintf(stderr, "PlaneData::ny           = %e\n", it->second->ny);
      com->fprintf(stderr, "PlaneData::nz           = %e\n", it->second->nz);
      com->fprintf(stderr, "\n");
    }
  }

  // initial conditions for EmbeddedStructure
  if(!embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator it;
    for (it =embed.embedIC.pointMap.dataMap.begin(); 
         it!=embed.embedIC.pointMap.dataMap.end();
         it++){
      com->fprintf(stderr, "PointData::tag          = %d\n", it->first);
      com->fprintf(stderr, "PointData::fluidModelID = %d\n", it->second->fluidModelID);
      com->fprintf(stderr, "PointData::density      = %e\n", it->second->initialConditions.density);
      com->fprintf(stderr, "PointData::pressure     = %e\n", it->second->initialConditions.pressure);
      com->fprintf(stderr, "PointData::temperature  = %e\n", it->second->initialConditions.temperature);
      com->fprintf(stderr, "PointData::mach         = %e\n", it->second->initialConditions.mach);
      com->fprintf(stderr, "PointData::velocity     = %e\n", it->second->initialConditions.velocity);
      com->fprintf(stderr, "PointData::alpha        = %e\n", it->second->initialConditions.alpha);
      com->fprintf(stderr, "PointData::beta         = %e\n", it->second->initialConditions.beta);
      com->fprintf(stderr, "PointData::x            = %e\n", it->second->x);
      com->fprintf(stderr, "PointData::y            = %e\n", it->second->y);
      com->fprintf(stderr, "PointData::z            = %e\n", it->second->z);
      com->fprintf(stderr, "\n");
    }
  }
}
