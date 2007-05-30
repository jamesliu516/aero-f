#include <IoData.h>

#include <Communicator.h>
#include <parser/Assigner.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#ifdef COUGAR
extern int optind;
extern "C" int getopt(int, char **, char *);
#endif

//------------------------------------------------------------------------------

InputData::InputData() 
{ 

  prefix = "";
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
  podFile2 = "";
}

//------------------------------------------------------------------------------

void InputData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 14, father);

  new ClassStr<InputData>(ca, "Prefix", this, &InputData::prefix);
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
  new ClassStr<InputData>(ca, "PODData2", this, &InputData::podFile2);

}

//------------------------------------------------------------------------------

PreconditionData::PreconditionData()
{
  mach = 1.0;
  k = 1.0;
  betav = 0.0;
}

//------------------------------------------------------------------------------

void PreconditionData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name,3,father);

  new ClassDouble<PreconditionData>(ca,"Mach", this, &PreconditionData::mach);
  new ClassDouble<PreconditionData>(ca,"k", this, &PreconditionData::k);
  new ClassDouble<PreconditionData>(ca,"Betav", this, &PreconditionData::betav);

}

//------------------------------------------------------------------------------
OutputData::OutputData()
{

}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  transient.setup("Postpro", ca);
  restart.setup("Restart", ca);

}

//------------------------------------------------------------------------------

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
  psensor = "";
  csdles = "";
  csdvms = "";
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
  lift = "";
  tavlift = "";
  hydrostaticlift = "";
  hydrodynamiclift = "";
  residuals = "";
  podFile = "";
  romFile = "";
  philevel = "";

  frequency = 0;
  length = 1.0;
  surface = 1.0;
  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;

}

//------------------------------------------------------------------------------

void TransientData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 54, father);

  new ClassStr<TransientData>(ca, "Prefix", this, &TransientData::prefix);
  new ClassStr<TransientData>(ca, "Solution", this, &TransientData::solutions);
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
  new ClassStr<TransientData>(ca, "PressureSensor", this, &TransientData::psensor);
  new ClassStr<TransientData>(ca, "CsDLES", this, &TransientData::csdles);
  new ClassStr<TransientData>(ca, "CsDVMS", this, &TransientData::csdvms);
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
  new ClassStr<TransientData>(ca, "LiftandDrag", this, &TransientData::lift);
  new ClassStr<TransientData>(ca, "HydroStaticLiftandDrag", this, &TransientData::hydrostaticlift);
  new ClassStr<TransientData>(ca, "HydroDynamicLiftandDrag", this, &TransientData::hydrodynamiclift);
  new ClassStr<TransientData>(ca, "TavLiftandDrag", this, &TransientData::tavlift);
  new ClassStr<TransientData>(ca, "Residual", this, &TransientData::residuals);
  new ClassInt<TransientData>(ca, "Frequency", this, &TransientData::frequency);
  new ClassDouble<TransientData>(ca, "Length", this, &TransientData::length);
  new ClassDouble<TransientData>(ca, "Surface", this, &TransientData::surface);
  new ClassDouble<TransientData>(ca, "XM", this, &TransientData::x0);
  new ClassDouble<TransientData>(ca, "YM", this, &TransientData::y0);
  new ClassDouble<TransientData>(ca, "ZM", this, &TransientData::z0);
  new ClassStr<TransientData>(ca, "PODData", this, &TransientData::podFile);
  new ClassStr<TransientData>(ca, "ROM", this, &TransientData::romFile);
  new ClassStr<TransientData>(ca, "Philevel", this, &TransientData::philevel);

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

}

//------------------------------------------------------------------------------

void RestartData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassToken<RestartData>(ca, "Type", this, 
			      reinterpret_cast<int RestartData::*>(&RestartData::type), 2,
			      "Single", 0, "Double", 1);

  new ClassStr<RestartData>(ca, "Prefix", this, &RestartData::prefix);
  new ClassStr<RestartData>(ca, "Solution", this, &RestartData::solutions);
  new ClassStr<RestartData>(ca, "Position", this, &RestartData::positions);
  new ClassStr<RestartData>(ca, "LevelSet", this, &RestartData::levelsets);
  new ClassStr<RestartData>(ca, "RestartData", this, &RestartData::data);
  new ClassInt<RestartData>(ca, "Frequency", this, &RestartData::frequency);

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

}

//------------------------------------------------------------------------------

void RestartParametersData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassInt<RestartParametersData>(ca, "Iteration", this, &RestartParametersData::iteration);
  new ClassDouble<RestartParametersData>(ca, "Time", this, &RestartParametersData::etime);
  new ClassDouble<RestartParametersData>(ca, "TimeStep1", this, &RestartParametersData::dt_nm1);
  new ClassDouble<RestartParametersData>(ca, "TimeStep2", this, &RestartParametersData::dt_nm2);
  new ClassDouble<RestartParametersData>(ca, "Residual", this, &RestartParametersData::residual);
  new ClassDouble<RestartParametersData>(ca, "Energy", this, &RestartParametersData::energy);

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
  test = REGULAR;
  verbose = 4;

}

//------------------------------------------------------------------------------

void ProblemData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  
  new ClassToken<ProblemData>
    (ca, "Type", this, 
     reinterpret_cast<int ProblemData::*>(&ProblemData::alltype), 21, 
     "Steady", 0, "Unsteady", 1, "AcceleratedUnsteady", 2, "SteadyAeroelastic", 3, 
     "UnsteadyAeroelastic", 4, "AcceleratedUnsteadyAeroelastic", 5,
     "SteadyThermal", 6, "UnsteadyThermal", 7, "SteadyAeroThermoElastic", 8, 
     "UnsteadyAeroThermoElastic", 9, "Forced", 10, "AcceleratedForced", 11, 
     "RigidRoll", 12, "RbmExtractor", 13, "UnsteadyLinearizedAeroelastic", 14,
     "UnsteadyLinearized", 15, "PODConstruction", 16, "ROMAeroelastic", 17,
     "ROM", 18, "ForcedLinearized", 19, "PODInterpolation", 20);

  new ClassToken<ProblemData>
    (ca, "Mode", this, 
     reinterpret_cast<int ProblemData::*>(&ProblemData::mode), 2,
     "NonDimensional", 0, "Dimensional", 1);

  new ClassToken<ProblemData>
    (ca, "Prec", this,
     reinterpret_cast<int ProblemData::*>(&ProblemData::prec), 2,
     "NonPreconditioned", 0, "LowMach", 1);

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
  reynolds_lambda = -1.0;
  length = 1.0;

}

//------------------------------------------------------------------------------

void ReferenceStateData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 8, father);
  
  new ClassDouble<ReferenceStateData>(ca, "Mach", this, &ReferenceStateData::mach);
  new ClassDouble<ReferenceStateData>(ca, "Velocity", this, &ReferenceStateData::velocity);
  new ClassDouble<ReferenceStateData>(ca, "Density", this, &ReferenceStateData::density);
  new ClassDouble<ReferenceStateData>(ca, "Pressure", this, &ReferenceStateData::pressure);
  new ClassDouble<ReferenceStateData>(ca, "Temperature", this, &ReferenceStateData::temperature);
  new ClassDouble<ReferenceStateData>(ca, "Reynolds", this, &ReferenceStateData::reynolds_mu);
  new ClassDouble<ReferenceStateData>(ca, "ReynoldsL", this, &ReferenceStateData::reynolds_lambda);
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
  temperature = -1.0;
  delta = -1.0;

}

//------------------------------------------------------------------------------

void BcsWallData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  
  new ClassToken<BcsWallData>(ca, "Type", this, 
			      reinterpret_cast<int BcsWallData::*>(&BcsWallData::type), 2, 
			      "Isothermal", 0, "Adiabatic", 1);
  new ClassToken<BcsWallData>(ca, "Integration", this, 
			      reinterpret_cast<int BcsWallData::*>(&BcsWallData::integration), 2, 
			      "WallFunction", 1, "Full", 2);
  new ClassDouble<BcsWallData>(ca, "Temperature", this, &BcsWallData::temperature);
  new ClassDouble<BcsWallData>(ca, "Delta", this, &BcsWallData::delta);

}

//------------------------------------------------------------------------------

BcsHydroData::BcsHydroData()
{

  type = NONE;
  gravity = 0.0;
  depth   = 0.0;
  alpha   = 0.0;
  beta    = 0.0;

}

//------------------------------------------------------------------------------

void BcsHydroData::setup(const char *name, ClassAssigner *father)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<BcsHydroData>(ca, "Type", this, 
                               reinterpret_cast<int BcsHydroData::*>(&BcsHydroData::type),
                               2, "None", 0, "Gravity", 1);
  new ClassDouble<BcsHydroData>(ca, "Gravity", this, &BcsHydroData::gravity);
  new ClassDouble<BcsHydroData>(ca, "Depth", this, &BcsHydroData::depth);
  new ClassDouble<BcsHydroData>(ca, "Alpha", this, &BcsHydroData::alpha);
  new ClassDouble<BcsHydroData>(ca, "Beta", this, &BcsHydroData::beta);

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
  pressureConstant = 0.0;

}

//------------------------------------------------------------------------------

void GasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  
  new ClassToken<GasModelData>(ca, "Type", this, 
			       reinterpret_cast<int GasModelData::*>(&GasModelData::type), 2,
			       "Ideal", 0, "Stiffened", 1);
  new ClassDouble<GasModelData>(ca, "SpecificHeatRatio", this, 
                                &GasModelData::specificHeatRatio);
  new ClassDouble<GasModelData>(ca, "IdealGasConstant", this, 
                                &GasModelData::idealGasConstant);
  new ClassDouble<GasModelData>(ca, "PressureConstant", this,
                                &GasModelData::pressureConstant);

}

//------------------------------------------------------------------------------
                                                                                                  
LiquidModelData::LiquidModelData()
{
  //data is available in 'Properties of Water and Steam', Wagner & Kruse, Springer, 1997
  type = COMPRESSIBLE;
  check = YES;
  specificHeatRatio = 1.0;
  Cv          = -1.0;
  k1water     = 2.07e9;
  k2water     = 7.15;
  Prefwater   = -1.0;
  RHOrefwater = -1.0;
                                                                                                  
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
                                                                                                  
  ClassAssigner *ca = new ClassAssigner(name, 11, father);
                                                                                                  
  new ClassToken<LiquidModelData>(ca, "Type", this,
            reinterpret_cast<int LiquidModelData::*>(&LiquidModelData::type), 1,
            "Compressible", 0);
  new ClassToken<LiquidModelData>(ca, "Check", this,
            reinterpret_cast<int LiquidModelData::*>(&LiquidModelData::check), 2,
            "Yes", 0, "No", 1);
  new ClassDouble<LiquidModelData>(ca, "SpecificHeatRatio", this,
                                &LiquidModelData::specificHeatRatio);
  new ClassDouble<LiquidModelData>(ca, "Cv", this, &LiquidModelData::Cv);
  new ClassDouble<LiquidModelData>(ca, "k1", this, &LiquidModelData::k1water);
  new ClassDouble<LiquidModelData>(ca, "k2", this, &LiquidModelData::k2water);
  new ClassDouble<LiquidModelData>(ca, "Pressure", this, &LiquidModelData::Prefwater);
  new ClassDouble<LiquidModelData>(ca, "Density", this, &LiquidModelData::RHOrefwater);
  new ClassDouble<LiquidModelData>(ca, "Pref", this, &LiquidModelData::Pref);
  new ClassDouble<LiquidModelData>(ca, "alpha", this, &LiquidModelData::alpha);
  new ClassDouble<LiquidModelData>(ca, "beta", this, &LiquidModelData::beta);
                                                                                                  
}

//------------------------------------------------------------------------------
                                                                                                  
FluidModelData::FluidModelData()
{

  fluid = GAS;
  pmin = -1.e6;
                                                                                                  
}
                                                                                                  
//------------------------------------------------------------------------------
                                                                                                  
void FluidModelData::setup(const char *name, ClassAssigner *father)
{
                                                                                                  
  ClassAssigner *ca = new ClassAssigner(name, 4, father);
                                                                                                  
  new ClassToken<FluidModelData>(ca, "Fluid", this,
                                 reinterpret_cast<int FluidModelData::*>(&FluidModelData::fluid), 3,
                                 "PerfectGas", 0, "Liquid", 1, "StiffenedGas", 0);
  new ClassDouble<FluidModelData>(ca, "Pmin", this, &FluidModelData::pmin);
                                                                                                  
  gasModel.setup("GasModel", ca);
  liquidModel.setup("LiquidModel", ca);
                                                                                                  
};

//------------------------------------------------------------------------------

ViscosityModelData::ViscosityModelData()
{

  type = SUTHERLAND;
  sutherlandReferenceTemperature = 110.6;
  sutherlandConstant = 1.458e-6;

}

//------------------------------------------------------------------------------

void ViscosityModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<ViscosityModelData>
    (ca, "Type", this, reinterpret_cast<int ViscosityModelData::*>
                       (&ViscosityModelData::type), 4,
                       "Constant", 0, "Sutherland", 1, "Prandtl", 2, "Water", 3);
  new ClassDouble<ViscosityModelData>(ca, "SutherlandReferenceTemperature", this, 
				      &ViscosityModelData::sutherlandReferenceTemperature);
  new ClassDouble<ViscosityModelData>(ca, "SutherlandConstant", this, 
				      &ViscosityModelData::sutherlandConstant);
  
}

//------------------------------------------------------------------------------

ThermalCondModelData::ThermalCondModelData()
{
  
  type = CONSTANT_PRANDTL;
  prandtl = 0.72;

}

//------------------------------------------------------------------------------

void ThermalCondModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);
 
  new ClassToken<ThermalCondModelData>
     (ca, "Type", this, reinterpret_cast<int ThermalCondModelData::*>
     (&ThermalCondModelData::type), 2,
     "ConstantPrandtl", 0, "Water", 1);

  new ClassDouble<ThermalCondModelData>
     (ca, "Prandtl", this, &ThermalCondModelData::prandtl);

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

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<LESModelData>
    (ca, "Type", this, reinterpret_cast<int LESModelData::*>
     (&LESModelData::type), 4, "Smagorinsky", 0, "Dynamic", 1, "VMS", 2, "DynamicVMS", 3);
  new ClassToken<LESModelData>
	(ca, "Delta", this,
	 reinterpret_cast<int LESModelData::*>(&LESModelData::delta), 2,
         "Volume", 0, "Side", 1);

  sma.setup("Smagorinsky", ca);
  dles.setup("Dynamic", ca);
  vms.setup("VMS", ca);
  dvms.setup("DynamicVMS", ca);

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
                                                                                                        
SphereData::SphereData()
{
                                                                                                        
  type = Fluid1;
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  r      = -1.0;
  p      = -1.0;
  rho    = -1.0;          
  t      = -1.0;
  mach   = -1.0;
	vel    = -1.0;
}
                                                                                                        
//------------------------------------------------------------------------------
                                                                                                        
void SphereData::setup(const char *name, ClassAssigner *father)
{
                                                                                                        
  ClassAssigner *ca = new ClassAssigner(name, 10, father);
                                                                                                        
  new ClassToken<SphereData>
    (ca, "Type", this, reinterpret_cast<int SphereData::*>
     (&SphereData::type), 2,
     "Fluid1", 0, "Fluid2", 1);
                                                                                                        
  new ClassDouble<SphereData>
    (ca, "Center_x", this, &SphereData::cen_x);
  new ClassDouble<SphereData>
    (ca, "Center_y", this, &SphereData::cen_y);
  new ClassDouble<SphereData>
    (ca, "Center_z", this, &SphereData::cen_z);
  new ClassDouble<SphereData>
    (ca, "Radius", this, &SphereData::r);
  new ClassDouble<SphereData>
    (ca, "Pressure", this, &SphereData::p);
  new ClassDouble<SphereData>
    (ca, "Density", this, &SphereData::rho);
  new ClassDouble<SphereData>
    (ca, "Temperature", this, &SphereData::t);
  new ClassDouble<SphereData>
    (ca, "Mach", this, &SphereData::mach);
  new ClassDouble<SphereData>
    (ca, "Velocity", this, &SphereData::vel);
}

//------------------------------------------------------------------------------
                                                                                                        
PlaneData::PlaneData()
{
                                                                                                        
  type = Fluid1;
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;   
  normal_x  = 0.0;
  normal_y  = 0.0;
  normal_z  = 0.0;   
  p      = 0.0;   
  rho    = 0.0;
  t      = 0.0;
  mach   = 0.0;
}
                                                                                                        
//------------------------------------------------------------------------------
                                                                                                        
void PlaneData::setup(const char *name, ClassAssigner *father)
{
                                                                                                        
  ClassAssigner *ca = new ClassAssigner(name, 11, father);
                                                                                                        
  new ClassToken<PlaneData>
    (ca, "Type", this, reinterpret_cast<int PlaneData::*>
     (&PlaneData::type), 2,
     "Fluid1", 0, "Fluid2", 1);
                                                                                                        
  new ClassDouble<PlaneData>
    (ca, "Center_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData>
    (ca, "Center_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData>
    (ca, "Center_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData>
    (ca, "Normal_x", this, &PlaneData::normal_x);
  new ClassDouble<PlaneData>
    (ca, "Normal_y", this, &PlaneData::normal_y);
  new ClassDouble<PlaneData>
    (ca, "Normal_z", this, &PlaneData::normal_z);
  new ClassDouble<PlaneData>
    (ca, "Pressure", this, &PlaneData::p);
  new ClassDouble<PlaneData>
    (ca, "Density", this, &PlaneData::rho);
  new ClassDouble<PlaneData>
    (ca, "Temperature", this, &PlaneData::t);
  new ClassDouble<PlaneData>
    (ca, "Mach", this, &PlaneData::mach);
}
//------------------------------------------------------------------------------
                                                                                                        
ICData::ICData()
{
                                                                                                        
  nspheres  = 0;

  sphere[0] = &s1;
  sphere[1] = &s2;
  sphere[2] = &s3;
  sphere[3] = &s4;
  sphere[4] = &s5;
  sphere[5] = &s6;
  sphere[6] = &s7;
  sphere[7] = &s8;
  sphere[8] = &s9;
  sphere[9] = &s10;
                                                                                                        
  nplanes  = 0;
  plane[0] = &p1;
  plane[1] = &p2;
  plane[2] = &p3;
  plane[3] = &p4;
  plane[4] = &p5;
  plane[5] = &p6;
  plane[6] = &p7;
  plane[7] = &p8;
  plane[8] = &p9;
  plane[9] = &p10;
                                                                                                        
}
                                                                                                        
//------------------------------------------------------------------------------
                                                                                                        
void ICData::setup(const char *name, ClassAssigner *father)
{
                                                                                                        
  ClassAssigner *ca = new ClassAssigner(name, 20, father);

  s1.setup("Sphere1", ca);
  s2.setup("Sphere2", ca);
  s3.setup("Sphere3", ca);
  s4.setup("Sphere4", ca);
  s5.setup("Sphere5", ca);
  s6.setup("Sphere6", ca);
  s7.setup("Sphere7", ca);
  s8.setup("Sphere8", ca);
  s9.setup("Sphere9", ca);
  s10.setup("Sphere10", ca);

  p1.setup("Plane1", ca);
  p2.setup("Plane2", ca);
  p3.setup("Plane3", ca);
  p4.setup("Plane4", ca);
  p5.setup("Plane5", ca);
  p6.setup("Plane6", ca);
  p7.setup("Plane7", ca);
  p8.setup("Plane8", ca);
  p9.setup("Plane9", ca);
  p10.setup("Plane10", ca);
                                                                                                        
}

//------------------------------------------------------------------------------
                                                                                                        
MultiFluidData::MultiFluidData()
{

  method = GHOSTFLUID_FOR_POOR;
	problem = BUBBLE;
	typePhaseChange = RIEMANN_SOLUTION;
  localtime  = GLOBAL;
  typeTracking = LINEAR;
  bandlevel = 3;
  subIt = 10;
	frequency = 0;
}
                                                                                                        
//------------------------------------------------------------------------------
                                                                                                        
void MultiFluidData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 9, father);

  new ClassToken<MultiFluidData>(ca, "Method", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::method), 4,
             "None", 0, "GhostFluidForThePoor", 1, "GhostFluidWithRiemann", 2,
             "RealFluidMethod", 3);
  new ClassToken<MultiFluidData>(ca, "Problem", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::problem), 2,
             "Bubble", 0, "ShockTube", 1);
	new ClassToken<MultiFluidData>(ca, "PhaseChange", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::typePhaseChange), 3,
             "None", 0, "RiemannSolution", 1, "Extrapolation", 2);
  new ClassToken<MultiFluidData>(ca, "FictitiousTimeStepping", this,
		         reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::localtime),2,
             "Global", 0, "Local", 1);
  new ClassToken<MultiFluidData>(ca, "InterfaceTracking", this,
             reinterpret_cast<int MultiFluidData::*>(&MultiFluidData::typeTracking),2,
             "Linear", 0, "Gradient", 1);
  new ClassInt<MultiFluidData>(ca, "BandLevel", this,
             &MultiFluidData::bandlevel);
  new ClassInt<MultiFluidData>(ca, "SubIt", this,
             &MultiFluidData::subIt);
  new ClassInt<MultiFluidData>(ca, "Frequency", this,
             &MultiFluidData::frequency);
                                                                                                        
  icd.setup("InitialConditions", ca);
                                                                                                        
}

//------------------------------------------------------------------------------
                                                                                                  
EquationsData::EquationsData()
{
 
  dimension = 3;
  type = EULER;
  numPhase = 1;
                                                                                                  
}

//------------------------------------------------------------------------------

void EquationsData::setup(const char *name, ClassAssigner *father)
{
                                                                                                  
  ClassAssigner *ca = new ClassAssigner(name, 8, father);
                                                                                                  
  new ClassInt<EquationsData>(ca, "Dimension", this, &EquationsData::dimension);
                                                                                                  
  new ClassToken<EquationsData>(ca, "Type", this,
                                reinterpret_cast<int EquationsData::*>(&EquationsData::type), 2,
                                "Euler", 0, "NavierStokes", 1);
                                                                                                  
  new ClassInt<EquationsData>(ca, "NumPhases", this,
                                &EquationsData::numPhase);
                                                                                                  
  fluidModel.setup("FluidModel", ca);
  fluidModel2.setup("FluidModel2", ca);
  viscosityModel.setup("ViscosityModel", ca);
  thermalCondModel.setup("ThermalConductivityModel", ca);
  tc.setup("TurbulenceClosure", ca);
 
}

//------------------------------------------------------------------------------

SchemeData::SchemeData()
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

  ClassAssigner *ca = new ClassAssigner(name, 11, father);
  
  new ClassToken<SchemeData>
    (ca, "AdvectiveOperator", this, 
     reinterpret_cast<int SchemeData::*>(&SchemeData::advectiveOperator), 2,
     "FiniteVolume", 0, "Galerkin", 1);

  new ClassToken<SchemeData>
    (ca, "Flux", this, 
     reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 2,
     "Roe", 0, "VanLeer", 1);

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

}

//------------------------------------------------------------------------------

void SFixData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassDouble<SFixData>(ca, "X0", this, &SFixData::x0);
  new ClassDouble<SFixData>(ca, "Y0", this, &SFixData::y0);
  new ClassDouble<SFixData>(ca, "Z0", this, &SFixData::z0);
  new ClassDouble<SFixData>(ca, "Radius", this, &SFixData::r);

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
         reinterpret_cast<int BoundarySchemeData::*>(&BoundarySchemeData::type), 4,
         "StegerWarming", 0, "ConstantExtrapolation", 1, "LinearExtrapolation", 2, 
	 "Ghidaglia", 3);
                                                                                                  
}

//------------------------------------------------------------------------------

SchemesData::SchemesData()
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
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 2,
     "RungeKutta4", 0, "RungeKutta2", 1);
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
			 reinterpret_cast<int PcData::*>(&PcData::type), 6, 
			 "Identity", 0, "Jacobi", 1, "As", 2, "Ras", 3, "Has", 4, "Aas", 5);

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
  coupling = WEAK;
  mvp = H1;
  jacobian = APPROXIMATE;
  normals = AUTO;
  velocities = AUTO_VEL;

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
    (ca, "Coupling", this, 
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::coupling), 2,
     "Weak", 0, "Strong", 1);

  new ClassToken<ImplicitData>
    (ca, "MatrixVectorProduct", this,
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::mvp), 4,
     "FiniteDifference", 0, "Approximate", 1, "Exact", 2, 
     "ApproximateFiniteDifference", 3);

  new ClassToken<ImplicitData>
    (ca, "FluxJacobian", this, 
     reinterpret_cast<int ImplicitData::*>(&ImplicitData::jacobian), 3,
     "FiniteDifference", 0, "Approximate", 1, "Exact", 2);

  new ClassToken<ImplicitData>
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
  
  newton.setup("Newton", ca);

}

//------------------------------------------------------------------------------

TsData::TsData()
{

  type = IMPLICIT;
  typeTimeStep = AUTO;
  typeClipping = FREESTREAM;

  prec = NO_PREC;
  viscousCst = 0.0;

  maxIts = 100;
  eps = 1.e-6;
  timestep = -1.0;
  maxTime = 1.e99;

  residual = -1;
  cfl0 = 5.0;
  cflCoef1 = 0.0;
  cflCoef2 = 0.0;
  cflMax = 1000.0;
  ser = 0.7;

  output = "";

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 18, father);
  
  new ClassToken<TsData>(ca, "Type", this, 
			 reinterpret_cast<int TsData::*>(&TsData::type), 2,
			 "Explicit", 0, "Implicit", 1);
  new ClassToken<TsData>(ca, "TypeTimeStep", this, 
			 reinterpret_cast<int TsData::*>(&TsData::typeTimeStep), 2,
			 "Local", 1, "Global", 2);
  new ClassToken<TsData>(ca, "Clipping", this, 
			 reinterpret_cast<int TsData::*>(&TsData::typeClipping), 3,
			 "None", 0, "AbsoluteValue", 1, "Freestream", 2);
  new ClassToken<TsData>(ca, "Prec", this,
                         reinterpret_cast<int TsData::*>(&TsData::prec), 2,
                         "NonPreconditioned", 0, "LowMach", 1);
  new ClassDouble<TsData>(ca, "ViscosityParameter", this, &TsData::viscousCst);

  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "Eps", this, &TsData::eps);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);
  new ClassInt<TsData>(ca, "Residual", this, &TsData::residual);
  new ClassDouble<TsData>(ca, "Cfl0", this, &TsData::cfl0);
  new ClassDouble<TsData>(ca, "Cfl1", this, &TsData::cflCoef1);
  new ClassDouble<TsData>(ca, "Cfl2", this, &TsData::cflCoef2);
  new ClassDouble<TsData>(ca, "CflMax", this, &TsData::cflMax);
  new ClassDouble<TsData>(ca, "Ser", this, &TsData::ser);
  new ClassStr<TsData>(ca, "Output", this, &TsData::output);  

  expl.setup("Explicit", ca);
  implicit.setup("Implicit", ca);

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
     reinterpret_cast<int DefoMeshMotionData::*>(&DefoMeshMotionData::element), 4,
     "LinearFiniteElement", 0, "NonLinearFiniteElement", 1, "TorsionalSprings", 2, "BallVertexSprings", 3);

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

  new ClassDouble<AeroelasticData>(ca, "InsidePressure", this, &AeroelasticData::pressure);
  new ClassDouble<AeroelasticData>(ca, "DisplacementScaling", this, &AeroelasticData::displacementScaling);
  new ClassDouble<AeroelasticData>(ca, "ForceScaling", this, &AeroelasticData::forceScaling);
  new ClassDouble<AeroelasticData>(ca, "PowerScaling", this, &AeroelasticData::powerScaling);

}

//------------------------------------------------------------------------------

ForcedData::ForcedData()
{

  type = FLEXIBLE;
  positions = "";
  amplification = 1.0;
  frequency = -1.0;
  timestep = -1.0;

}

//------------------------------------------------------------------------------

void ForcedData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);
  
  new ClassToken<ForcedData>
    (ca, "Type", this, 
     reinterpret_cast<int ForcedData::*>(&ForcedData::type), 2,
     "Rigid", 0, "Flexible", 1);

  new ClassStr<ForcedData>(ca, "Position", this, &ForcedData::positions);
  new ClassDouble<ForcedData>(ca, "Amplification", this, &ForcedData::amplification);
  new ClassDouble<ForcedData>(ca, "Frequency", this, &ForcedData::frequency);
  new ClassDouble<ForcedData>(ca, "TimeStep", this, &ForcedData::timestep);

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
  eps = 1e-4;
  eps2 = -1;
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

void LinearizedData::setup(char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 15, father);

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
                                                        
                                                        
void PadeData::setup(char *name, ClassAssigner *father)
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

  computeForces = UNSPECIFIED;
  forceResults = NO;

  rotationID = -1;
  velocity = 0.0;
}

//------------------------------------------------------------------------------

static RootClassAssigner nullAssigner;
Assigner *SurfaceData::getAssigner()  {

  ClassAssigner *ca = new ClassAssigner("normal", 8, &nullAssigner);

  new ClassDouble<SurfaceData>(ca, "Nx", this, &SurfaceData::nx);
  new ClassDouble<SurfaceData>(ca, "Ny", this, &SurfaceData::ny);
  new ClassDouble<SurfaceData>(ca, "Nz", this, &SurfaceData::nz);
  new ClassToken<SurfaceData> (ca, "ComputeForces", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::computeForces), 2, "False", 0, "True", 1);
  new ClassToken<SurfaceData> (ca, "SeparateForces", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::forceResults), 2, "False", 0, "True", 1);
  new ClassToken<SurfaceData> (ca, "SeparateFile", this, reinterpret_cast<int SurfaceData::*>(&SurfaceData::forceResults), 2, "False", 0, "True", 1);

  new ClassInt<SurfaceData>(ca, "VelocityID", this, &SurfaceData::rotationID);
  new ClassDouble<SurfaceData>(ca, "Velocity", this, &SurfaceData::velocity);

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

  ClassAssigner *ca = new ClassAssigner("normal", 8, &nullAssigner);

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

Assigner *VolumeData::getAssigner()  {

  ClassAssigner *ca = new ClassAssigner("normal", 17, &nullAssigner);

  new ClassDouble<VolumeData>(ca, "Ix", this, &VolumeData::iprimex);
  new ClassDouble<VolumeData>(ca, "Iy", this, &VolumeData::iprimey);
  new ClassDouble<VolumeData>(ca, "Iz", this, &VolumeData::iprimez);

  new ClassDouble<VolumeData>(ca, "Jx", this, &VolumeData::jprimex);
  new ClassDouble<VolumeData>(ca, "Jy", this, &VolumeData::jprimey);
  new ClassDouble<VolumeData>(ca, "Jz", this, &VolumeData::jprimez);

  new ClassDouble<VolumeData>(ca, "Kx", this, &VolumeData::kprimex);
  new ClassDouble<VolumeData>(ca, "Ky", this, &VolumeData::kprimey);
  new ClassDouble<VolumeData>(ca, "Kz", this, &VolumeData::kprimez);

  new ClassDouble<VolumeData>(ca, "Alphax", this, &VolumeData::alphax);
  new ClassDouble<VolumeData>(ca, "Alphay", this, &VolumeData::alphay);
  new ClassDouble<VolumeData>(ca, "Alphaz", this, &VolumeData::alphaz);

  new ClassDouble<VolumeData>(ca, "Betax", this, &VolumeData::betax);
  new ClassDouble<VolumeData>(ca, "Betay", this, &VolumeData::betay);
  new ClassDouble<VolumeData>(ca, "Betaz", this, &VolumeData::betaz);

  new ClassDouble<VolumeData>(ca, "Idr", this, &VolumeData::idr);

  new ClassDouble<VolumeData>(ca, "Ldr", this, &VolumeData::ldr);

  return ca;
}

//------------------------------------------------------------------------------

void PorousMedia::setup(const char *name)  {

  ClassAssigner *ca = new ClassAssigner(name, 0, 0);
  volumeMap.setup("VolumeData", 0);
}

//------------------------------------------------------------------------------

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
  mf.setup("MultiFluid");
  schemes.setup("Space");
  ts.setup("Time");
  dmesh.setup("MeshMotion");
  rmesh.setup("Accelerated");
  aero.setup("Aeroelastic");
  forced.setup("Forced");
  linearizedData.setup("Linearized");
  surfaces.setup("Surfaces");
  rotations.setup("Velocity");
  porousmedia.setup("PorousMedia");

}

//------------------------------------------------------------------------------

void IoData::readCmdFile() 
{

  extern int yyCmdfparse();

  setupCmdFileVariables();

  cmdFilePtr = freopen(cmdFileName, "r", stdin);
 
  if (!cmdFilePtr) {
    com->fprintf(stderr,"*** Error: could not open \'%s\'\n", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    com->fprintf(stderr,"*** Error: command file contained parsing errors\n");
    exit(error);
  }

  if (input.rstdata[0] != 0) {
    char *name = new char[strlen(input.prefix) + strlen(input.rstdata) + 1];
    if (strncmp(input.rstdata, "/", 1) == 0)
      sprintf(name, "%s", input.rstdata);
    else
      sprintf(name, "%s%s", input.prefix, input.rstdata);
    FILE *fp = freopen(name, "r", stdin);
    if (!fp) {
      com->fprintf(stderr, "*** Error: could not open \'%s\'\n", name);
      exit(-1);
    }
    error = yyCmdfparse();
    if (error) {
      com->fprintf(stderr, "*** Error: parameter file contained parsing errors\n");
      exit(error);
    }
  }

  resetInputValues();
  error = checkFileNames();
  error += checkInputValues();
  error += checkSolverValues();
  if (error) {
    com->fprintf(stderr, "*** Error: command file contained %d error%s\n", 
		 error, error>1? "s":"");
    exit(-1);
  }

  com->setMaxVerbose(problem.verbose);

}

//------------------------------------------------------------------------------

void IoData::resetInputValues()
{

  // part 1

  for (int i=0; i<ProblemData::SIZE; ++i)
    problem.type[i] = false;

  if (problem.alltype == ProblemData::_UNSTEADY_ ||
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
      problem.alltype == ProblemData::_POD_CONSTRUCTION_ ||
      problem.alltype == ProblemData::_ROM_AEROELASTIC_ ||
      problem.alltype == ProblemData::_ROM_ ||
      problem.alltype == ProblemData::_INTERPOLATION_)
    problem.type[ProblemData::LINEARIZED] = true;

  // part 2

  if (problem.type[ProblemData::AERO] || problem.type[ProblemData::THERMO] || 
      problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_ || 
      problem.alltype == ProblemData::_ROM_AEROELASTIC_)
    problem.mode = ProblemData::DIMENSIONAL;

  if (!problem.type[ProblemData::UNSTEADY])
    ts.implicit.type = ImplicitData::BACKWARD_EULER;

  if (ts.type == TsData::IMPLICIT &&
      (ts.implicit.newton.failsafe == NewtonData<KspFluidData>::YES ||
        ts.implicit.newton.failsafe == NewtonData<KspFluidData>::ALWAYS)) {
    if (schemes.ns.reconstruction == SchemeData::CONSTANT)
      ts.implicit.newton.failsafe = NewtonData<KspFluidData>::NO;
    //else
      //schemes.ns.advectiveOperator = SchemeData::FE_GALERKIN;
  }

  if (problem.type[ProblemData::AERO] && !problem.type[ProblemData::UNSTEADY]) {
    ts.implicit.normals = ImplicitData::LATEST_CFG;
    ts.implicit.velocities = ImplicitData::ZERO;
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

  if (eqs.numPhase == 1){
    if (eqs.fluidModel.fluid == FluidModelData::GAS)
      if(eqs.fluidModel.gasModel.type == GasModelData::IDEAL)
        com->fprintf(stderr, " ----- PERFECT GAS SIMULATION -----\n");
      else if(eqs.fluidModel.gasModel.type == GasModelData::STIFFENED)
        com->fprintf(stderr, " ----- STIFFENED GAS SIMULATION -----\n");
      else{
        com->fprintf(stderr, " ----- UNDEFINED SIMULATION -----\n -----> exiting program\n");
        exit(1);
      }
    else if (eqs.fluidModel.fluid == FluidModelData::LIQUID)
      com->fprintf(stderr, " ----- BAROTROPIC LIQUID SIMULATION -----\n");
    else {
      com->fprintf(stderr, " ----- UNDEFINED SIMULATION -----\n -----> exiting program\n");
      exit(1);
    }
  }
  else if (eqs.numPhase == 2) {
     com->fprintf(stderr, " ----- TWO-PHASE FLOW SIMULATION -----\n");
     if (eqs.fluidModel.fluid == FluidModelData::GAS &&
         eqs.fluidModel2.fluid == FluidModelData::GAS)
       if (eqs.fluidModel.gasModel.type == GasModelData::IDEAL && 
           eqs.fluidModel2.gasModel.type == GasModelData::IDEAL)
         com->fprintf(stderr, " ----- PERFECT GAS-PERFECT GAS SIMULATION -----\n");
       else if(eqs.fluidModel.gasModel.type == GasModelData::STIFFENED &&
               eqs.fluidModel2.gasModel.type == GasModelData::STIFFENED)
         com->fprintf(stderr, " ----- STIFFENED GAS-STIFFENED GAS SIMULATION -----\n");
       else if(eqs.fluidModel.gasModel.type == GasModelData::IDEAL && 
               eqs.fluidModel2.gasModel.type == GasModelData::STIFFENED)
         com->fprintf(stderr, " ----- PERFECT GAS-STIFFENED GAS SIMULATION -----\n");
       else if(eqs.fluidModel.gasModel.type == GasModelData::STIFFENED &&
               eqs.fluidModel2.gasModel.type == GasModelData::IDEAL)
         com->fprintf(stderr, " ----- STIFFENED GAS-PERFECT GAS SIMULATION -----\n");
       else{
         com->fprintf(stderr, " ----- UNDEFINED SIMULATION -----\n -----> exiting program\n");
         exit(1);
       }
     else if (eqs.fluidModel.fluid == FluidModelData::LIQUID &&
              eqs.fluidModel2.fluid == FluidModelData::GAS)
       if (eqs.fluidModel2.gasModel.type == GasModelData::IDEAL)
         com->fprintf(stderr, " ---- BAROTROPIC LIQUID-PERFECT GAS SIMULATION -----\n");
       else if (eqs.fluidModel2.gasModel.type == GasModelData::STIFFENED)
         com->fprintf(stderr, " ---- BAROTROPIC LIQUID-STIFFENED GAS SIMULATION -----\n");
       else{
         com->fprintf(stderr, " ----- UNDEFINED SIMULATION -----\n -----> exiting program");
         exit(1);
       }
     else if (eqs.fluidModel.fluid == FluidModelData::LIQUID &&
              eqs.fluidModel2.fluid == FluidModelData::LIQUID)
       com->fprintf(stderr, " ---- BAROTROPIC LIQUID-BAROTROPIC LIQUID SIMULATION ----\n");
     else{
       com->fprintf(stderr, " ----- GAS-LIQUID SIMULATIONS ARE NOT DEFINED -----\n -----> exiting program");
			 /* need to change the following to run a liquid in gas simulation:
			  *
			  */	
       //exit(1);
     }
  }
  else{
    com->fprintf(stderr, " ----- ONLY SINGLE AND TWO-PHASE FLOW SIMULATIONS ARE POSSIBLE ----\n -----> exiting program");
    exit(1);
  }
                                                                                                  
  if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
    if(schemes.bc.type == BoundarySchemeData::STEGER_WARMING &&
       bc.inlet.type != BcsFreeStreamData::INTERNAL){
      com->fprintf(stderr, "*** Error: for an hydrodynamic simulation, numerical treatment of boundary conditions needs to be an extrapolation method\n");
      exit(1);
    }
  }
  // to avoid having inlet nodes when computing Internal BCs
  if(bc.inlet.type == BcsFreeStreamData::INTERNAL)
    schemes.bc.type = BoundarySchemeData::STEGER_WARMING;

  int nCoarseFreq = 0;
  if (problem.alltype == ProblemData::_POD_CONSTRUCTION_) {
    // Assign the values for the coarse freq to an arra
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
}

//------------------------------------------------------------------------------

int IoData::checkFileNames()
{

  int error = 0;
  
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
       strcmp(input.match, "") == 0) {
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
/*
    double k1water = eqs.fluidModel.liquidModel.k1water;
    double k2water = eqs.fluidModel.liquidModel.k2water;
    double Prefwater = eqs.fluidModel.liquidModel.Prefwater;
    double RHOrefwater = eqs.fluidModel.liquidModel.RHOrefwater;
    double Pref, awater, bwater;                                                                                              
    if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
      Pref = -k1water/k2water;
      awater = (Prefwater + k1water/k2water)/pow(RHOrefwater, k2water);
      bwater = k2water;
    }                                                                                                
    double k1water2, k2water2, Prefwater2, RHOrefwater2;
    double Pref2, awater2, bwater2;
                                                                                                  
    if(eqs.numPhase == 2){
      k1water2 = eqs.fluidModel2.liquidModel.k1water;
      k2water2 = eqs.fluidModel2.liquidModel.k2water;
      Prefwater2 = eqs.fluidModel2.liquidModel.Prefwater;
      RHOrefwater2 = eqs.fluidModel2.liquidModel.RHOrefwater;
      if (eqs.fluidModel2.fluid == FluidModelData::LIQUID){
        Pref2 = -k1water2/k2water2;
        awater2 = (Prefwater2 + k1water2/k2water2)/pow(RHOrefwater2, k2water2);
        bwater2 = k2water2;
      }
                                                                                                  
    }
*/
  int error = 0;
  error += checkInputValuesEssentialBC();
  error += checkInputValuesStateEquation();

  error += checkInputValuesNonDimensional(); 
  error += checkInputValuesDimensional(); 

  checkInputValuesTurbulence();
                                                                                                  
  checkInputValuesDefaultOutlet();
                                                                                                  
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


  error += checkInputValuesInitializeMulti();
                                                                                                  
  return error;

}
//------------------------------------------------------------------------------
int IoData::checkInputValuesNonDimensional()
{
  int error = 0;
  if (problem.mode == ProblemData::NON_DIMENSIONAL) {
                                                                                                  
    if (eqs.type != EquationsData::EULER) {
      if (ref.reynolds_mu < 0.0) {
        com->fprintf(stderr, "*** Error: no valid Reynolds number (%d) given\n", ref.reynolds_mu);
        ++error;
      }
                                                                                                  
      if (ref.reynolds_lambda < 0.0 && eqs.fluidModel.fluid == FluidModelData::LIQUID){
        com->fprintf(stderr, "*** Error: no valid Reynolds number (%d) given for lambda\n", ref.reynolds_lambda);
        ++error;
      }
      else if (eqs.fluidModel.fluid == FluidModelData::GAS && eqs.fluidModel.gasModel.type == GasModelData::IDEAL)
        ref.reynolds_lambda = -3.0/2.0 * ref.reynolds_mu;
                                                                                                  
                                                                                                  
      if (eqs.viscosityModel.type == ViscosityModelData::SUTHERLAND && ref.temperature < 0.0) {
        com->fprintf(stderr, "*** Error: no valid reference temperature (%d) given\n", ref.temperature);
        ++error;
      }
    }

    double gamma = eqs.fluidModel.gasModel.specificHeatRatio;

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

    if (bc.inlet.density < 0.0)
      bc.inlet.density = 1.0;
    if (bc.inlet.pressure < 0.0)
      if (eqs.fluidModel.fluid == FluidModelData::GAS)
				if(ref.mach>0.0)
          bc.inlet.pressure = bc.inlet.pressure / (gamma * ref.mach * ref.mach * (bc.inlet.pressure + eqs.fluidModel.gasModel.pressureConstant));
		    else
					com->fprintf(stderr, "*** Error: no valid Mach number for non-dimensional simulation\n");
      else if(eqs.fluidModel.fluid == FluidModelData::LIQUID)
        bc.inlet.pressure = Prefwater/((Prefwater+k1water/k2water)*k2water*ref.mach*ref.mach);
    if (bc.inlet.temperature < 0.0 && eqs.fluidModel.fluid == FluidModelData::LIQUID){
      com->fprintf(stderr, "*** Error: no valid non-dimensionalized temperature (%d) given\n", bc.inlet.temperature);
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
int IoData::checkInputValuesDimensional()
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

  double Cv = eqs.fluidModel.liquidModel.Cv;
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
                                                                                                        
  double R2, gamma2;
  double Cv2, k1water2, k2water2, Prefwater2, RHOrefwater2;
  double Pref2, awater2, bwater2;
                                                                                                        
  if(eqs.numPhase == 2){
    R2 = eqs.fluidModel2.gasModel.idealGasConstant;
    gamma2 = eqs.fluidModel2.gasModel.specificHeatRatio;
                                                                                                      
    Cv2 = eqs.fluidModel2.liquidModel.Cv;
    k1water2 = eqs.fluidModel2.liquidModel.k1water;
    k2water2 = eqs.fluidModel2.liquidModel.k2water;
    Prefwater2 = eqs.fluidModel2.liquidModel.Prefwater;
    RHOrefwater2 = eqs.fluidModel2.liquidModel.RHOrefwater;
    if(eqs.fluidModel2.fluid == FluidModelData::LIQUID){
      Pref2 = -k1water2/k2water2;
      awater2 = (Prefwater2 + k1water2/k2water2)/pow(RHOrefwater2, k2water2);
      bwater2 = k2water2;
    }
  }
  if (problem.mode == ProblemData::DIMENSIONAL) {
    if (bc.inlet.pressure < 0.0)  {
      if(eqs.fluidModel.fluid == FluidModelData::GAS){
        com->fprintf(stderr, "*** Error: no valid inlet pressure (%f) given\n", bc.inlet.pressure);
        ++error;
      }else if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
        bc.inlet.pressure = Prefwater;
        com->fprintf(stderr, "*** Warning: inlet pressure set to reference pressure for Tait's EOS\n");
      }
    }
    if (bc.inlet.density < 0.0)
      if(eqs.fluidModel.fluid == FluidModelData::GAS){
        com->fprintf(stderr, "*** Error: no valid inlet density (%f) given\n", bc.inlet.density);
        ++error;
      }
                                                                                                        
    if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
			if (mf.problem == MultiFluidData::BUBBLE)
        bc.inlet.density = pow( (bc.inlet.pressure - Pref)/awater, 1.0/bwater);
    }

    if(bc.inlet.temperature < 0.0) {
      if(eqs.fluidModel.fluid == FluidModelData::GAS)
        bc.inlet.temperature = (bc.inlet.pressure + gamma*Pstiff)/(R*bc.inlet.density);
      else if(eqs.fluidModel.fluid == FluidModelData::LIQUID){
        com->fprintf(stderr, "*** Error: no valid inlet temperature (%f) given\n", bc.inlet.temperature);
        ++error;
      }
    }

    if(eqs.fluidModel.fluid == FluidModelData::GAS){
      if (ref.density < 0.0)
        ref.density = bc.inlet.density;
      if (ref.pressure < 0.0)
        ref.pressure = bc.inlet.pressure;
			if (ref.mach <= 0.0)
				ref.mach = ref.velocity /sqrt(gamma * (ref.pressure+Pstiff) / ref.density);
      double velocity = ref.mach * sqrt(gamma * (ref.pressure+Pstiff) / ref.density);
      ref.temperature = (ref.pressure + gamma*Pstiff)/ (ref.density * R);
      double viscosity = eqs.viscosityModel.sutherlandConstant * sqrt(ref.temperature) /
        (1.0 + eqs.viscosityModel.sutherlandReferenceTemperature/ref.temperature);
      ref.reynolds_mu = velocity * ref.length * ref.density / viscosity;
      ref.reynolds_lambda = -3.0 * ref.reynolds_mu/2.0;
        //as we are considering a gas whose Lame coefficients respect the Stokes relation (3*lambda+2*mu=0)
        // we have a relation between Re_mu and Re_lambda

      ref.rv.mode = RefVal::DIMENSIONAL;
      ref.rv.density = ref.density;
      ref.rv.velocity = velocity;
      ref.rv.pressure = ref.density * velocity*velocity;
      ref.rv.temperature = gamma*(gamma - 1.0) * ref.mach*ref.mach * (ref.pressure+Pstiff)/(R*ref.density);
//      ref.rv.temperature = gamma*(gamma - 1.0) * ref.mach*ref.mach * ref.temperature;
      ref.rv.viscosity_mu = viscosity;
      ref.rv.viscosity_lambda = -2.0/3.0 * viscosity;
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

      if (linearizedData.eps2 < 0)
        linearizedData.eps2 = ref.rv.time;
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
      double viscosity_mu = 0.000001; //NEEDS TO BE CHANGED ACCORDING TO THE VISCO LAW FOR WATER
      double viscosity_lambda = 0.00000000001;
      ref.reynolds_mu = velocity * ref.length * ref.density / viscosity_mu;
      ref.reynolds_lambda = velocity * ref.length * ref.length / viscosity_lambda;
      ref.rv.mode = RefVal::DIMENSIONAL;

      ref.rv.density = ref.density;
      ref.rv.velocity = velocity;
      ref.rv.pressure = ref.density * velocity*velocity;
      ref.rv.temperature = velocity * velocity / Cv;
                                                                                                        
                                                                                                        
      ref.rv.viscosity_mu = viscosity_mu;
      ref.rv.viscosity_lambda = viscosity_lambda;
      ref.rv.nutilde = viscosity_mu / ref.density;
      ref.rv.kenergy = velocity*velocity;
      ref.rv.epsilon = velocity*velocity*velocity / ref.length;
      ref.rv.time = ref.length / velocity;
      ref.rv.force = ref.density * velocity*velocity * ref.length*ref.length;
      ref.rv.energy = ref.density * velocity*velocity * ref.length*ref.length*ref.length;
      ref.rv.power = ref.density * velocity*velocity*velocity * ref.length*ref.length;
      ref.rv.tvelocity = velocity / aero.displacementScaling;
                                                                                                        
      ref.rv.tforce = ref.rv.force / aero.forceScaling;
      ref.rv.tpower = ref.rv.power / aero.powerScaling;
                                                                                                        
      eqs.fluidModel.liquidModel.Pref = Pref / (ref.rv.density * ref.rv.velocity *ref.rv.velocity);
      eqs.fluidModel.liquidModel.alpha = awater * pow(ref.rv.density, bwater - 1.0)/(ref.rv.velocity *ref.rv.velocity);
      eqs.fluidModel.liquidModel.beta  = bwater;

    }
    if (eqs.fluidModel2.fluid == FluidModelData::LIQUID) {
      eqs.fluidModel2.liquidModel.Pref = Pref2 / (ref.rv.density * ref.rv.velocity *ref.rv.velocity);
      eqs.fluidModel2.liquidModel.alpha = awater2 * pow(ref.rv.density, bwater2 - 1.0)/(ref.rv.velocity *ref.rv.velocity);
      eqs.fluidModel2.liquidModel.beta  = bwater2;
    }

    eqs.fluidModel.pmin /= ref.rv.pressure;
    eqs.fluidModel2.pmin /= ref.rv.pressure;

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
    linearizedData.stepsize = ts.timestep;
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
      rmesh.vpts[j]->time     /= ref.rv.time;
      rmesh.vpts[j]->velocityX /= ref.rv.velocity;
      rmesh.vpts[j]->velocityY /= ref.rv.velocity;
      rmesh.vpts[j]->velocityZ /= ref.rv.velocity;
    }
    aero.pressure /= ref.rv.pressure;
    forced.timestep /= ref.rv.time;
    forced.frequency *= ref.rv.time;

    bc.hydro.gravity /= ref.rv.velocity / ref.rv.time;
    bc.hydro.depth /= ref.length;
    bc.hydro.alpha *= acos(-1.0) / 180.0;
    bc.hydro.beta *= acos(-1.0) / 180.0;
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

  return error;

}
//------------------------------------------------------------------------------
int IoData::checkInputValuesStateEquation()
{
  int error = 0;

  if (eqs.fluidModel.fluid == FluidModelData::LIQUID){
    if (eqs.fluidModel.liquidModel.Prefwater < 0.0){
      com->fprintf(stderr, "*** Error: no valid reference pressure (%e) given for Tait's EOS\n", eqs.fluidModel.liquidModel.Prefwater);
      ++error;
    }
    if (eqs.fluidModel.liquidModel.RHOrefwater < 0.0){
      com->fprintf(stderr, "*** Error: no valid reference density (%e) given for Tait's EOS\n", eqs.fluidModel.liquidModel.RHOrefwater);
      ++error;
    }
    if ( eqs.fluidModel.liquidModel.Cv < 0.0 ) {
      com->fprintf(stderr, "*** Error: no valid reference specific heat coefficient (%e) given\n", eqs.fluidModel.liquidModel.Cv);
      ++error;
    }
  }
  else if(eqs.fluidModel.fluid         == FluidModelData::GAS &&
	  eqs.fluidModel.gasModel.type == GasModelData::IDEAL)
    eqs.fluidModel.gasModel.pressureConstant = 0.0;
                                                                                                  
  if (eqs.numPhase == 2){ 
    if (eqs.fluidModel2.fluid == FluidModelData::LIQUID){
      if (eqs.fluidModel2.liquidModel.Prefwater < 0.0){
        com->fprintf(stderr, "*** Error: no valid reference pressure (%e) given for 2nd Tait's EOS\n", eqs.fluidModel2.liquidModel.Prefwater);
        ++error;
      }
      if (eqs.fluidModel2.liquidModel.RHOrefwater < 0.0){
        com->fprintf(stderr, "*** Error: no valid reference density (%e) given for 2nd Tait's EOS\n", eqs.fluidModel2.liquidModel.RHOrefwater);
        ++error;
      }
      if ( eqs.fluidModel2.liquidModel.Cv < 0.0 ) {
        com->fprintf(stderr, "*** Error: no valid reference specific heat coefficient (%e) given\n", eqs.fluidModel2.liquidModel.Cv);
        ++error;
      }
    }
    else if(eqs.fluidModel2.fluid == FluidModelData::GAS && eqs.fluidModel2.gasModel.type == GasModelData::IDEAL)
      eqs.fluidModel2.gasModel.pressureConstant = 0.0;
  }
  return error;
}
//------------------------------------------------------------------------------
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

}

//------------------------------------------------------------------------------


int IoData::checkSolverValues()
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
    com->fprintf(stderr, "*** Error: no valid wall temperature (%d) given\n", bc.wall.temperature);
    ++error;
  }

  return error;

}  

//------------------------------------------------------------------------------

int IoData::checkInputValuesInitializeMulti()
{

  int error = 0;
  if(eqs.numPhase == 2){
    for(int i=0;i<10; i++){
      if(mf.icd.sphere[i]->r>0.0){
        mf.icd.nspheres +=1;
        if(eqs.fluidModel2.fluid == FluidModelData::GAS){
          if(mf.icd.sphere[i]->p < 0.0){
            ++error;
            com->fprintf(stderr, "*** Error: no valid pressure specified for initial sphere %d\n", i+1);
          }
          if(mf.icd.sphere[i]->rho < 0.0){
            ++error;
            com->fprintf(stderr, "*** Error: no valid density specified for initial sphere %d\n", i+1);
          }
        }
        if(eqs.fluidModel2.fluid == FluidModelData::LIQUID){
          if(mf.icd.sphere[i]->t < 0.0){
            ++error;
            com->fprintf(stderr, "*** Error: no valid temperature specified for initial sphere %d\n", i+1);
          }
          if(mf.icd.sphere[i]->p < 0.0 && mf.icd.sphere[i]->rho < 0.0){
            ++error;
            com->fprintf(stderr, "***Error: at least presssure or density must be specified for initial sphere %d\n", i+1);
          }else{
            double k1water2 = eqs.fluidModel2.liquidModel.k1water;
            double k2water2 = eqs.fluidModel2.liquidModel.k2water;
            double Prefwater2 = eqs.fluidModel2.liquidModel.Prefwater;
            double RHOrefwater2 = eqs.fluidModel2.liquidModel.RHOrefwater;
            double Pref2 = -k1water2/k2water2;
            double awater2 = (Prefwater2 + k1water2/k2water2)/pow(RHOrefwater2, k2water2);
            double bwater2 = k2water2;
            if(mf.icd.sphere[i]->p < 0.0){     
              mf.icd.sphere[i]->p = Pref2 + awater2*pow(mf.icd.sphere[i]->rho,bwater2);
            }
            else{
              mf.icd.sphere[i]->rho = pow( (mf.icd.sphere[i]->p - Pref2)/awater2, 1.0/bwater2);
              com->fprintf(stderr, "*** Warning: the pressure inside the bubble is set using the state equation and the input density\n");
            }
          }
        }
        mf.icd.sphere[i]->rho /= ref.rv.density;
        mf.icd.sphere[i]->p /= ref.rv.pressure;
        mf.icd.sphere[i]->t /= ref.rv.temperature;
        mf.icd.sphere[i]->vel /= ref.rv.velocity;
      }
    }
  }

  return error;
}
