#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <RefVal.h>
#include <cstdio>
#include <map>
#include "parser/ParseTree.h"
#include "parser/Dictionary.h"

using std::map;

class Assigner;
class ClassAssigner;
class Communicator;

//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:

  map<int, DataType *> dataMap;
  void setup(const char *name, ClassAssigner *);
  ~ObjectMap()
    {
      for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
	{
	  delete it->second;
	}
    }
};

//------------------------------------------------------------------------------
struct FluidRemapData {

  FluidRemapData();
  ~FluidRemapData() {}
  
  int oldID,newID;

  void setup(const char *, ClassAssigner * = 0);

  Assigner *getAssigner();
};

struct OneDimensionalInputData {

  OneDimensionalInputData();
  ~OneDimensionalInputData() {}
  const char* file;
  double x0,y0,z0;

  ObjectMap<FluidRemapData> fluidRemap;

  Assigner *getAssigner();

  void setup(const char *, ClassAssigner * = 0);
};
  

struct InputData {

  const char *prefix;
  const char *connectivity;
  const char *geometry;
  const char *decomposition;
  const char *cpumap;
  const char *match;
  const char *d2wall;
  const char *perturbed;
  const char *solutions;
  const char *positions;
  const char *levelsets;
  const char *rstdata;
  const char *podFile;
  const char *podFile2;
  const char *strModesFile;
  const char *embeddedSurface;

// Included (MB)
  const char *shapederivatives;

  ObjectMap< OneDimensionalInputData > oneDimensionalInput;

  InputData();
  ~InputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TransientData {

  const char *prefix;
  const char *solutions;
  const char *density;
  const char *tavdensity;
  const char *mach;
  const char *speed;
  const char *wtmach;
  const char *wtspeed;
  const char *absvelocity;
  const char *tavmach;
  const char *pressure;
  const char *diffpressure;
  const char *tavpressure;
  const char *hydrostaticpressure;
  const char *hydrodynamicpressure;
  const char *pressurecoefficient;
  const char *temperature;
  const char *tavtemperature;
  const char *totalpressure;
  const char *tavtotalpressure;
  const char *vorticity;
  const char *tavvorticity;
  const char *nutturb;
  const char *kturb;
  const char *epsturb;
  const char *eddyvis;
  const char *dplus;
  const char *sfric;
  const char *tavsfric;
  const char *psensor;
  const char *csdles;
  const char *tavcsdles;
  const char *csdvms;
  const char *tavcsdvms;
  const char *mutOmu;
  const char *velocity;
  const char *tavvelocity;
  const char *displacement;
  const char *tavdisplacement;
  const char *flightDisplacement;
  const char *localFlightDisplacement;
  const char *forces;
  const char *tavforces;
  const char *hydrostaticforces;
  const char *hydrodynamicforces;
  const char *generalizedforces;
  const char *lift;
  const char *tavlift;
  const char *hydrostaticlift;
  const char *hydrodynamiclift;
  const char *residuals;
  const char *materialVolumes;
  const char *conservation;
  const char *podFile;
  const char *romFile;
  const char *philevel;
  const char *controlvolume;
  const char* fluidid;
// Included (MB)
  const char *velocitynorm;
  const char *dSolutions;
  const char *dDensity;
  const char *dMach;
  const char *dPressure;
  const char *dTemperature;
  const char *dTotalpressure;
  const char *dNutturb;
  const char *dEddyvis;
  const char *dVelocityScalar;
  const char *dVelocityVector;
  const char *dDisplacement;
  const char *dForces;

  const char *tempnormalderivative;
  const char *surfaceheatflux;
  const char *heatfluxes;

  const char *sparseGrid;

  int frequency;
  double x0, y0, z0;
  double length;
  double surface;

  TransientData();
  ~TransientData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RestartData {

  enum Type {SINGLE = 0, DOUBLE = 1} type;

  const char *prefix;

  const char *solutions;
  const char *positions;
  const char *levelsets;
  const char *data;

  int frequency;

  RestartData();
  ~RestartData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct OutputData {

  TransientData transient;
  RestartData restart;

  OutputData();
  ~OutputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RestartParametersData {

  int iteration;

  double etime;
  double dt_nm1;
  double dt_nm2;
  double residual;
  double energy;

  RestartParametersData();
  ~RestartParametersData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ProblemData {

  enum Type {UNSTEADY = 0, ACCELERATED = 1, AERO = 2, THERMO = 3, FORCED = 4,
	     ROLL = 5, RBM = 6, LINEARIZED = 7, SIZE = 8};
  bool type[SIZE];

  enum AllType {_STEADY_ = 0, _UNSTEADY_ = 1, _ACC_UNSTEADY_ = 2, _STEADY_AEROELASTIC_ = 3,
		_UNSTEADY_AEROELASTIC_ = 4, _ACC_UNSTEADY_AEROELASTIC_ = 5,
		_STEADY_THERMO_ = 6, _UNSTEADY_THERMO_ = 7, _STEADY_AEROTHERMOELASTIC_ = 8,
		_UNSTEADY_AEROTHERMOELASTIC_ = 9, _FORCED_ = 10, _ACC_FORCED_ = 11,
		_ROLL_ = 12, _RBM_ = 13, _UNSTEADY_LINEARIZED_AEROELASTIC_ = 14,
		_UNSTEADY_LINEARIZED_ = 15, _POD_CONSTRUCTION_ = 16,
		_ROM_AEROELASTIC_ = 17, _ROM_ = 18, _FORCED_LINEARIZED_ = 19,
		_INTERPOLATION_ = 20, _STEADY_SENSITIVITY_ANALYSIS_ = 21,
		_SPARSEGRIDGEN_ = 22, _ONE_DIMENSIONAL_ } alltype;
  enum Mode {NON_DIMENSIONAL = 0, DIMENSIONAL = 1} mode;
  enum Test {REGULAR = 0} test;
  enum Prec {NON_PRECONDITIONED = 0, PRECONDITIONED = 1} prec;
  enum Framework {BODYFITTED = 0, EMBEDDED = 1} framework;
  int verbose;

  ProblemData();
  ~ProblemData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PreconditionData {

  double mach;
  double cmach;
  double k;
  double betav;
  double shockreducer;

  PreconditionData();
  ~PreconditionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ReferenceStateData {

  double mach;
  double velocity;
  double density;
  double pressure;
  double temperature;
  double reynolds_mu;
  double reynolds_lambda;
  double length;

// Included (MB)
  double dRe_mudMach;
  double dRe_lambdadMach;

  RefVal rv;

  ReferenceStateData();
  ~ReferenceStateData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsFreeStreamData {

  enum Type {EXTERNAL = 0, INTERNAL = 1} type;

  double mach;
  double velocity;
  double density;
  double pressure;
  double temperature;
  double nutilde;
  double kenergy;
  double eps;
  double alpha;
  double beta;

  BcsFreeStreamData();
  ~BcsFreeStreamData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsWallData {

  enum Type {ISOTHERMAL = 0, ADIABATIC = 1} type;
  enum Integration {AUTO = 0, WALL_FUNCTION = 1, FULL = 2} integration;
  enum Reconstruction {CONSTANT = 0, EXACT_RIEMANN = 1} reconstruction;

  double temperature;
  double delta;

  BcsWallData();
  ~BcsWallData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsHydroData {

  double depth;

  BcsHydroData();
  ~BcsHydroData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsData {

  BcsFreeStreamData inlet;
  BcsFreeStreamData outlet;
  BcsWallData wall;
  BcsHydroData hydro;

  BcsData();
  ~BcsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct GasModelData {

  enum Type {IDEAL = 0, STIFFENED = 1} type;

  double specificHeatRatio;
  double idealGasConstant;
  double pressureConstant;

  GasModelData();
  ~GasModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct JWLModelData {

  enum Type {IDEAL = 0, JWL = 1} type;

  double omega; // = specificHeatRatio-1.0
  double idealGasConstant;
  double A1,R1,rhoref,A2,R2;

  JWLModelData();
  ~JWLModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LiquidModelData {

  enum Type { COMPRESSIBLE = 0 } type;

  enum Check {YES = 0, NO = 1 } check;
  // the state equation is derived from a linearization of the bulk modulus wrt
  // pressure: K = k1 + k2 * P
  // the integration constant of the ODE is given by the couple (RHOrefwater,Prefwater)
  double specificHeatRatio;
  double Cv;
  double k1water;
  double k2water;
  double Prefwater;
  double RHOrefwater;

  //the state equation can be put in the form P=Pref+alpha*rho^beta
  // with Pref, alpha and beta derived from k1, k2 and the 'initial' couple
  double Pref;
  double alpha;
  double beta;

  LiquidModelData();
  ~LiquidModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct FluidModelData {

  enum Fluid { GAS = 0, LIQUID = 1, JWL = 2, UNDEFINED = 3} fluid;
  double pmin;

  GasModelData gasModel;
  JWLModelData jwlModel;
  LiquidModelData liquidModel;

  FluidModelData();
  ~FluidModelData() {}

  Assigner *getAssigner();
  void setup(const char *, ClassAssigner * = 0);
  //FluidModelData &operator=(const FluidModelData &);

};

//------------------------------------------------------------------------------

struct ViscosityModelData {

  enum Type {CONSTANT = 0, SUTHERLAND = 1, PRANDTL = 2, WATER = 3} type;

  double sutherlandReferenceTemperature;
  double sutherlandConstant;

  ViscosityModelData();
  ~ViscosityModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ThermalCondModelData {

  enum Type {CONSTANT_PRANDTL = 0, WATER = 1} type;

  double prandtl;

  ThermalCondModelData();
  ~ThermalCondModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PorousMedia  {

  double iprimex, iprimey, iprimez;
  double jprimex, jprimey, jprimez;
  double kprimex, kprimey, kprimez;
  double alphax,alphay,alphaz;
  double betax,betay,betaz;
  double idr,ldr;

  PorousMedia();
  //Assigner *getAssigner();
  void setup(const char*, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct InitialConditions {

  double mach;
  double velocity;
  double alpha, beta;
  double pressure;
  double density;
  double temperature;

  InitialConditions();
  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VolumeData  {

  enum Type {FLUID = 0, POROUS = 1} type;
  int fluidModelID;

  PorousMedia   porousMedia;
  InitialConditions initialConditions;

  VolumeData();
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SAModelData {

  double cb1;
  double cb2;
  double cw2;
  double cw3;
  double cv1;
  double cv2;
  double sigma;
  double vkcst;

  SAModelData();
  ~SAModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DESModelData {

  double cb1;
  double cb2;
  double cw2;
  double cw3;
  double cv1;
  double cv2;
  double cdes;
  double sigma;
  double vkcst;

  DESModelData();
  ~DESModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KEModelData {

  double sigma_k;
  double sigma_eps;
  double sigma_eps1;
  double sigma_eps2;
  double c_mu;

  KEModelData();
  ~KEModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TurbulenceModelData {

  enum Type {ONE_EQUATION_SPALART_ALLMARAS = 0, ONE_EQUATION_DES = 1, TWO_EQUATION_KE = 2} type;

  SAModelData sa;
  DESModelData des;
  KEModelData ke;

  TurbulenceModelData();
  ~TurbulenceModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SmagorinskyLESData {

  double c_s;

  SmagorinskyLESData();
  ~SmagorinskyLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct WaleLESData {

  double c_w;

  WaleLESData();
  ~WaleLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ClippingData {


  double cs_max;
  double pt_min;
  double pt_max;

  ClippingData();
  ~ClippingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DynamicLESData {

  ClippingData clip;

  DynamicLESData();
  ~DynamicLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------


struct DynamicVMSData {

  enum Type {D1VMSLES = 0, D2VMSLES = 1, D3VMSLES = 2} type;

  double c_s_prime;
  int agglomeration_width;
  int agglomeration_depth1;
  int agglomeration_depth2;

  ClippingData clip;

  DynamicVMSData();
  ~DynamicVMSData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VMSLESData {

  double c_s_prime;
  int agglomeration_width;
  int agglomeration_depth;

  VMSLESData();
  ~VMSLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LESModelData {

  enum Type {SMAGORINSKY = 0, DYNAMIC = 1, VMS = 2, DYNAMICVMS = 3, WALE = 4} type;
  enum Delta {VOLUME = 0, SIDE = 1} delta;

  SmagorinskyLESData sma;
  DynamicLESData dles;
  VMSLESData vms;
  DynamicVMSData dvms;
  WaleLESData wale;

  LESModelData();
  ~LESModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TBFixData {

  double x0;
  double y0;
  double z0;
  double x1;
  double y1;
  double z1;

  TBFixData();
  ~TBFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TripDomainData {

  TBFixData bfix;

  TripDomainData();
  ~TripDomainData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TurbulenceClosureData {

  enum Type {NONE = 0, EDDY_VISCOSITY = 1, LES = 2} type;

  double prandtlTurbulent;
  TurbulenceModelData tm;
  LESModelData les;
  TripDomainData tr;

  TurbulenceClosureData();
  ~TurbulenceClosureData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ProgrammedBurnData {

  int unburnedEOS,burnedEOS;
  double ignitionX0,ignitionY0,ignitionZ0;
  double e0;
  double cjDetonationVelocity;
  double cjPressure;
  double cjDensity;
  double cjEnergy;
  double ignitionTime;
  double factorB;
  double factorS;
  int ignited;
  int limitPeak;
  
  ProgrammedBurnData();
  ~ProgrammedBurnData();

  void setup(const char*, ClassAssigner* = 0);

};

struct SphereData {

  double cen_x, cen_y, cen_z, radius;
  int fluidModelID;
  InitialConditions initialConditions;

  ProgrammedBurnData programmedBurn;

  SphereData();
  ~SphereData() {}
  Assigner *getAssigner();

};
//------------------------------------------------------------------------------
struct PrismData {

  double cen_x, cen_y, cen_z, w_x,w_y,w_z;
  int fluidModelID;
  InitialConditions initialConditions;

  ProgrammedBurnData programmedBurn;

  bool inside(double x,double y,double z) const;

  PrismData();
  ~PrismData() {}
  Assigner *getAssigner();

};
//------------------------------------------------------------------------------

struct PlaneData {

  double cen_x, cen_y, cen_z, nx, ny, nz;
  int fluidModelID;
  InitialConditions initialConditions;

  PlaneData();
  ~PlaneData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct PointData {

  int fluidModelID;
  double x,y,z;
  InitialConditions initialConditions;

  PointData();
  ~PointData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MultiInitialConditionsData {

  ObjectMap<SphereData> sphereMap;
  ObjectMap<PrismData>  prismMap;
  ObjectMap<PlaneData>  planeMap;
  ObjectMap<PointData>  pointMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct SparseGridData {

  // to use already created sparse grids
  const char *tabulationFileName;
  int numberOfTabulations;

  // to generate sparse grids
  int verbose;

  int minPoints;
  int maxPoints;
  double relAccuracy;
  double absAccuracy;
  double dimAdaptDegree;

  double range1min, range1max, mapBaseValue1; int numDomainDim1;
  double range2min, range2max, mapBaseValue2; int numDomainDim2;
  double range3min, range3max, mapBaseValue3; int numDomainDim3;
  double range4min, range4max, mapBaseValue4; int numDomainDim4;
  double range5min, range5max, mapBaseValue5; int numDomainDim5;
  double range6min, range6max, mapBaseValue6; int numDomainDim6;
  typedef double Range[2];
  Range *range;
  double *mapBaseValue;
  int *numDomainDim;

  int numOutputs;
  int numInputs;

  SparseGridData();
  ~SparseGridData() {}
  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MultiFluidData {
  enum Method {NONE = 0, GHOSTFLUID_FOR_POOR = 1, GHOSTFLUID_WITH_RIEMANN} method;
  enum FictitiousTime {GLOBAL = 0, LOCAL = 1} localtime;
  enum InterfaceTracking {LINEAR = 0, GRADIENT = 1, HERMITE = 2} typeTracking;
  enum RiemannComputation {FE = 0, RK2 = 1, TABULATION2 = 2, TABULATION5 = 3} riemannComputation;
  int bandlevel;
  int subIt;
  double cfl;
  int frequency;
  double eps;
  int outputdiff;
  int numBlur;
  enum Problem {BUBBLE = 0, SHOCKTUBE = 1} problem;
  enum TypePhaseChange {ASIS = 0, RIEMANN_SOLUTION = 1, EXTRAPOLATION = 2} typePhaseChange;
  enum CopyCloseNodes {FALSE = 0, TRUE = 1} copy;
  enum LSInit {VOLUMES = 1, OLD = 0, GEOMETRIC = 2} lsInit;
  enum InterfaceType {FSF = 0, FF = 1, FSFandFF = 2} interfaceType;

  MultiInitialConditionsData multiInitialConditions;

  SparseGridData sparseGrid;

  MultiFluidData();
  ~MultiFluidData() {}
  void setup(const char *, ClassAssigner * = 0);
};


//------------------------------------------------------------------------------

struct Volumes  {

  ObjectMap<VolumeData> volumeMap;

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct EquationsData {

  int dimension;

  enum Type {EULER = 0, NAVIER_STOKES = 1} type;

  int numPhase;

// it is assumed that in a two-phase flow, fluidModel represents the surrounding fluid
// whereas fluidModel2 represents the fluid in the bubbles. This means that the
// uniform boundary conditions and the uniform initial conditions are computed
// with the values given to characterize the first fluid!

  double gravity_x, gravity_y, gravity_z;

  ObjectMap<FluidModelData> fluidModelMap;
  FluidModelData fluidModel;
  ViscosityModelData viscosityModel;
  ThermalCondModelData thermalCondModel;
  TurbulenceClosureData tc;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeData {

  enum AdvectiveOperator {FINITE_VOLUME = 0, FE_GALERKIN = 1} advectiveOperator;
  enum Flux {ROE = 0, VANLEER = 1, HLLE = 2, HLLC = 3} flux;

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruction;

  enum Limiter {NONE = 0, VANALBADA = 1, BARTH = 2, VENKAT = 3, P_SENSOR = 4} limiter;
  enum Gradient {LEAST_SQUARES = 0, GALERKIN = 1, NON_NODAL = 2} gradient;
  enum Dissipation {SECOND_ORDER = 0, SIXTH_ORDER = 1} dissipation;

  double beta;
  double gamma;
  double xiu;
  double xic;
  double eps;

  int allowsFlux;

  SchemeData(int allowsFlux = 1);
  ~SchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};
//------------------------------------------------------------------------------

struct CFixData {

  // nodal coordinates of first circle
  double x0;
  double y0;
  double z0;

  // nodal coordinates of second circle
  double x1;
  double y1;
  double z1;

  // radii of 1st and 2nd circle
  double r0;
  double r1;

  CFixData();
  ~CFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SFixData {

  double x0;
  double y0;
  double z0;
  double r;

  SFixData();
  ~SFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BFixData {

  double x0;
  double y0;
  double z0;
  double x1;
  double y1;
  double z1;

  BFixData();
  ~BFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeFixData {

  static const int num = 10;

  enum Symmetry {NONE = 0, X = 1, Y = 2, Z = 3} symmetry;
  double dihedralAngle;

  SFixData* spheres[num];
  BFixData* boxes[num];
  CFixData* cones[num];

  SFixData sfix1;
  SFixData sfix2;
  SFixData sfix3;
  SFixData sfix4;
  SFixData sfix5;
  SFixData sfix6;
  SFixData sfix7;
  SFixData sfix8;
  SFixData sfix9;
  SFixData sfix10;

  BFixData bfix1;
  BFixData bfix2;
  BFixData bfix3;
  BFixData bfix4;
  BFixData bfix5;
  BFixData bfix6;
  BFixData bfix7;
  BFixData bfix8;
  BFixData bfix9;
  BFixData bfix10;

  CFixData cfix1;
  CFixData cfix2;
  CFixData cfix3;
  CFixData cfix4;
  CFixData cfix5;
  CFixData cfix6;
  CFixData cfix7;
  CFixData cfix8;
  CFixData cfix9;
  CFixData cfix10;

  SchemeFixData();
  ~SchemeFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BoundarySchemeData {

  enum Type { STEGER_WARMING = 0,
	      CONSTANT_EXTRAPOLATION = 1,
	      LINEAR_EXTRAPOLATION = 2,
	      GHIDAGLIA = 3 } type;

  BoundarySchemeData();
  ~BoundarySchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemesData {

  SchemeData ns;
  SchemeData tm;
  SchemeData ls;
  SchemeFixData fixes;
  BoundarySchemeData bc;

  SchemesData();
  ~SchemesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ExplicitData {

//time-integration scheme used
  enum Type {RUNGE_KUTTA_4 = 0, RUNGE_KUTTA_2 = 1, FORWARD_EULER = 2, ONE_BLOCK_RK2 = 3, ONE_BLOCK_RK2bis = 4} type;

  ExplicitData();
  ~ExplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PcData {

  enum Type {IDENTITY = 0, JACOBI = 1, AS = 2, RAS = 3, ASH = 4, AAS = 5} type;
  enum Renumbering {NATURAL = 0, RCM = 1} renumbering;

  int fill;

  PcData();
  ~PcData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KspData {

  enum Type {RICHARDSON = 0, CG = 1, GMRES = 2, GCR = 3} type;
  enum EpsFormula {CONSTANT = 0, EISENSTADT = 1} epsFormula;
  enum CheckFinalRes {NO = 0, YES = 1} checkFinalRes;

  int maxIts;
  int numVectors;
  double eps;

  const char *output;

  PcData pc;

  KspData();
  ~KspData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KspFluidData {

  KspData ns;
  KspData tm;
  KspData lsi;

  KspFluidData();
  ~KspFluidData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

template<class GenericKrylov>
struct NewtonData {

  enum FailSafe {NO = 0, YES = 1, ALWAYS = 2} failsafe;
  int maxIts;
  double eps;
  double epsAbsRes, epsAbsInc;

  GenericKrylov ksp;

  NewtonData();
  ~NewtonData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ImplicitData {

  enum Type {BACKWARD_EULER = 0, CRANK_NICOLSON = 1, THREE_POINT_BDF = 2, FOUR_POINT_BDF = 3} type;
  enum Startup {REGULAR = 0, MODIFIED = 1} startup;
  enum TurbulenceModelCoupling {WEAK = 0, STRONG = 1} tmcoupling;
  enum Mvp {FD = 0, H1 = 1, H2 = 2, H1FD = 3} mvp;
  enum FiniteDifferenceOrder {FIRST_ORDER = 1, SECOND_ORDER = 2} fdOrder; 
  NewtonData<KspFluidData> newton;

  /// UH (09/10)
  /// This flag is not visible from the input file.
  /// It governs the computation of the Jacobian of the flux function,
  /// a component of the 'H' matrix (from the MatrixVectorProduct).
  enum FluxFcnJacobian {FINITE_DIFFERENCE = 0, APPROXIMATE = 1, EXACT = 2} ffjacobian;

  ImplicitData();
  ~ImplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData {

  enum Type {EXPLICIT = 0, IMPLICIT = 1} type;
  enum TypeTimeStep {AUTO = 0, LOCAL = 1, GLOBAL = 2} typeTimeStep;
  enum Clipping {NONE = 0, ABS_VALUE = 1, FREESTREAM = 2} typeClipping;

  enum Prec {NO_PREC = 0, PREC = 1} prec;
  double viscousCst;

  int maxIts;
  double eps;
  double timestep;
  double timestepinitial;
  double maxTime;

  int residual;
  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflMax;
  double ser;

  const char *output;

  ExplicitData expl;
  ImplicitData implicit;

  TsData();
  ~TsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DGCLData{

  enum Normals {AUTO = 0, IMPLICIT_FIRST_ORDER_GCL = 1, IMPLICIT_SECOND_ORDER_GCL = 2,
                IMPLICIT_FIRST_ORDER_EZGCL = 3, IMPLICIT_SECOND_ORDER_EZGCL = 4, IMPLICIT_THIRD_ORDER_EZGCL = 5,
                IMPLICIT_CURRENT_CFG = 6, IMPLICIT_LATEST_CFG = 7, EXPLICIT_RK2 = 8} normals;
  enum Velocities {AUTO_VEL = 0, IMPLICIT_BACKWARD_EULER_VEL = 1, IMPLICIT_THREE_POINT_BDF_VEL = 2,
                   IMPLICIT_IMPOSED_VEL = 3, IMPLICIT_IMPOSED_BACKWARD_EULER_VEL = 4,
                   IMPLICIT_IMPOSED_THREE_POINT_BDF_VEL = 5, IMPLICIT_ZERO = 6, EXPLICIT_RK2_VEL = 7} velocities;

  DGCLData();
  ~DGCLData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

// Included (MB)
struct SensitivityAnalysis {

  enum Method {DIRECT = 0, ADJOINT = 1} method;
  enum SensitivityComputation {ANALYTICAL = 0, SEMIANALYTICAL = 1,  FINITEDIFFERENCE = 2} scFlag;
  enum Compatible3D {OFF_COMPATIBLE3D = 0, ON_COMPATIBLE3D = 1} comp3d;
  enum AngleRadians {OFF_ANGLERAD = 0, ON_ANGLERAD = 1} angleRad;

  enum SensitivityMesh {OFF_SENSITIVITYMESH = 0, ON_SENSITIVITYMESH = 1} sensMesh;
  enum SensitivityMach {OFF_SENSITIVITYMACH = 0, ON_SENSITIVITYMACH = 1} sensMach;
  enum SensitivityAOA {OFF_SENSITIVITYALPHA = 0, ON_SENSITIVITYALPHA = 1} sensAlpha;
  enum SensitivityYAW {OFF_SENSITIVITYBETA = 0, ON_SENSITIVITYBETA = 1} sensBeta;

  // This flag repeats the linear solves until the number of iterations
  // is smaller than the maximum allowed.
  // Default Value = OFF_EXACTSOLUTION
  enum ExactSolution {OFF_EXACTSOLUTION = 0, ON_EXACTSOLUTION = 1} excsol;

  enum HomotopyComputation {OFF_HOMOTOPY = 0, ON_HOMOTOPY = 1} homotopy;
  enum FixSolution {NONEFIX = 0, PREVIOUSVALEUSFIX = 1} fixsol;

  double machref;
  double alpharef;
  double betaref;

  const char* meshderiv;
  const char* sensoutput;

  bool densFlag;
  bool pressFlag;
  bool apressFlag;

  int si;
  int sf;
  int avgsIt;

  double eps;
  double fres;

  KspData ksp;

  SensitivityAnalysis();
  ~SensitivityAnalysis() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SymmetryData {

  double nx, ny, nz;

  SymmetryData();
  ~SymmetryData() {};

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct BLMeshMotionData {
   int numLayers;
   int numIncrements;
   double power;
   double neiSelectionDistFactor;

   enum Type {PSEUDOSTRUCTURAL = 0, ALGEBRAIC = 1 } type;
   enum FractionalStrategy {Distance = 1, DotProduct = 2} fractionalStrategy;
   enum ClosestNeighbor {Fixed = 1, Variable = 0} bestNeiStrategy;
   int feedbackFrequency;

   BLMeshMotionData();
   ~BLMeshMotionData() {};

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct DefoMeshMotionData {

  enum Type {BASIC = 0, COROTATIONAL = 1} type;
  enum Element {LINEAR_FE = 0, NON_LINEAR_FE = 1, TORSIONAL_SPRINGS = 2, BALL_VERTEX = 3, NL_BALL_VERTEX = 4 } element;

  double volStiff;
  enum Mode {Recursive = 1, NonRecursive = 2} mode;
  int numIncrements;

  BLMeshMotionData blmeshmotion;
  NewtonData<KspData> newton;
  SymmetryData symmetry;

  DefoMeshMotionData();
  ~DefoMeshMotionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VelocityPoints {

  double time;
  double velocityX;
  double velocityY;
  double velocityZ;

  VelocityPoints();
  ~VelocityPoints() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ForcePoints {

  double time;
  double force;

  ForcePoints();
  ~ForcePoints() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RigidMeshMotionData {

  static const int num = 10;

  enum Tag {MACH = 0, TIME = 1, VELOCITY = 2} tag;
  enum LawType {VELOCITYPOINTS = 0, CONSTANTACCELERATION = 1} lawtype;

  double vx;
  double vy;
  double vz;

  double ax;
  double ay;
  double az;

  VelocityPoints* vpts[num];
  VelocityPoints vpts1;
  VelocityPoints vpts2;
  VelocityPoints vpts3;
  VelocityPoints vpts4;
  VelocityPoints vpts5;
  VelocityPoints vpts6;
  VelocityPoints vpts7;
  VelocityPoints vpts8;
  VelocityPoints vpts9;
  VelocityPoints vpts10;

  double timestep;

  RigidMeshMotionData();
  ~RigidMeshMotionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct AeroelasticData {

  enum Force {LAST = 0, AVERAGED = 1, LAST_KRIS = 2} force;
  double pressure;
  double displacementScaling;
  double forceScaling;
  double powerScaling;

  AeroelasticData();
  ~AeroelasticData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct HeavingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  double ax;
  double ay;
  double az;

  HeavingData();
  ~HeavingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct PitchingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  double alpha_in;
  double alpha_max;
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;

  PitchingData();
  ~PitchingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct DeformingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  const char *positions;

  double amplification;

  DeformingData();
  ~DeformingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct ForcedData {

  enum Type {HEAVING = 0, PITCHING = 1, DEFORMING = 2} type;

  double frequency;
  double timestep;

  HeavingData hv;
  PitchingData pt;
  DeformingData df;

  ForcedData();
  ~ForcedData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct PadeData {

  static const int num = 11;

  double freq[num];
  double freq1;
  double freq2;
  double freq3;
  double freq4;
  double freq5;
  double freq6;
  double freq7;
  double freq8;
  double freq9;
  double freq10;
  double freq11;

  int nPoints;
  int degNum;
  int degDen;

  PadeData();
  ~PadeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LinearizedData {

  enum PadeReconstruction {TRUE = 1, FALSE = 0} padeReconst;
  // Type is a defunct variable - it is not documented and only used for development and testing
  enum Type {DEFAULT = 0, ROM = 1, FORCED = 2} type;

  enum Domain {TIME = 0, FREQUENCY = 1} domain;
  enum InitialCondition {DISPLACEMENT = 0, VELOCITY = 1} initCond;

  double amplification;
  double frequency;
  double stepsize;
  double stepsizeinitial;
  double freqStep;
  double eps;
  double eps2;
  double tolerance;
  double refLength;
  const char *strModesFile;
  int modeNumber;
  int numSteps;
  int numPOD;
  int numStrModes;
  const char *romFile;

  PadeData pade;

  LinearizedData();
  ~LinearizedData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SurfaceData  {

  double nx, ny, nz;
  int sBit;
  static const int UNSPECIFIED = -1;
  enum ComputeForces {FALSE = 0, TRUE = 1 } computeForces;
  enum ForceResults {NO = 0, YES = 1} forceResults;
  int rotationID;
  double velocity;

  enum Type { ADIABATIC = 1, ISOTHERMAL = 2 } type;
  double temp;

  enum ComputeHeatPower {FALSE_HF = 0, TRUE_HF = 1 } computeHeatFluxes;
  enum HeatFluxResults {UNSPECIFIED_HF = -1, NO_HF = 0, YES_HF = 1} heatFluxResults;
  //the HF (Heat Flux) index ensures that there is no confusion with the force related data.

  SurfaceData();
  Assigner *getAssigner();
  void setBit(int b) { sBit = b; }
};

//------------------------------------------------------------------------------

struct Surfaces  {

  ObjectMap<SurfaceData> surfaceMap;
  void setup(const char *);
};

//------------------------------------------------------------------------------

struct PadeFreq  {

  ObjectMap<double> freqMap;
  void setup(const char *);
};

//------------------------------------------------------------------------------

struct RotationData  {

  double nx, ny, nz;
  double x0, y0, z0;
  double omega;
  enum InfRadius {FALSE = 0, TRUE = 1} infRadius;

  RotationData();
  Assigner *getAssigner();

};

//-----------------------------------------------------------------------------

struct Velocity  {

  ObjectMap<RotationData> rotationMap;

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct EmbeddedFramework { 

  enum IntersectorName {PHYSBAM = 0, FRG = 1} intersectorName;
  enum StructureNormal {ELEMENT_BASED = 0, NODE_BASED = 1} structNormal;
  enum EOSChange {NODAL_STATE = 0, RIEMANN_SOLUTION = 1} eosChange;
  enum ForceAlgorithm {RECONSTRUCTED_SURFACE = 0, CONTROL_VOLUME_BOUNDARY = 1} forceAlg;

  MultiInitialConditionsData embedIC;
  
  int nLevelset; //number of level-sets. Currently only consider bubbles.
  
  //Debug variables
  enum Coupling {TWOWAY = 0, ONEWAY = 1} coupling;
  enum Dim2Treatment {NO = 0, YES = 1} dim2Treatment;
  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruct;
  enum RiemannNormal {STRUCTURE = 0, FLUID = 1, AVERAGED_STRUCTURE = 2, AUTO = 3} riemannNormal;
  enum StructVelocity {COMPUTED_BY_STRUCTURE = 0, FINITE_DIFFERENCE = 1} structVelocity;
  
  EmbeddedFramework();
  ~EmbeddedFramework() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------
struct OneDimensionalInfo {
  enum CoordinateType {CARTESIAN = 0, CYLINDRICAL = 1, SPHERICAL = 2} coordType;
  enum VolumeType { CONSTANT_VOLUME = 0, REAL_VOLUME = 1} volumeType;
  double maxDistance; //mesh goes from 0 to maxDistance
  
  int numPoints; //mesh has numPoints elements
  int fluidId2;

  double interfacePosition;

  double density1, velocity1, pressure1,temperature1;
  double density2, velocity2, pressure2,temperature2;

  ProgrammedBurnData programmedBurn;

  OneDimensionalInfo();
  ~OneDimensionalInfo() {}

  void setup(const char *);
};
//------------------------------------------------------------------------------

struct ImplosionSetup {
  double Prate, Pinit;
  ImplosionSetup();
  ~ImplosionSetup() {}
  void setup(const char *);
};

struct MultigridInfo {
 
  MultigridInfo();
  ~MultigridInfo() {}
  void setup(const char *);

  const char* fineMesh;
  const char* coarseMesh;

  const char* fineDec;
  const char* coarseDec;

  const char* packageFile;
  const char* collectionFile;

  double radius0;
  double radiusf;
  int threshold;
};
//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

  Communicator *com;

public:

  InputData input;
  OutputData output;
  PreconditionData prec;
  RestartParametersData restart;
  ProblemData problem;
  ReferenceStateData ref;
  BcsData bc;
  EquationsData eqs;
  MultiFluidData mf;
  SchemesData schemes;

// Included (MB)
  SensitivityAnalysis sa;

  TsData ts;
  DGCLData dgcl;
  DefoMeshMotionData dmesh;
  RigidMeshMotionData rmesh;
  AeroelasticData aero;
  ForcedData forced;
  LinearizedData linearizedData;
  Surfaces surfaces;
  Velocity rotations;
  Volumes volumes;
  EmbeddedFramework embed;
  OneDimensionalInfo oneDimensionalInfo;
  ImplosionSetup implosion;

  MultigridInfo multigrid;

public:

  IoData(Communicator *);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void resetInputValues();
  int checkFileNames();
  int checkInputValues();
  int checkInputValuesAllEquationsOfState();
  int checkInputValuesProgrammedBurn();
  int checkProgrammedBurnLocal(ProgrammedBurnData& programmedBurn,
			       InitialConditions& IC);
  int checkInputValuesAllInitialConditions();
  void nonDimensionalizeAllEquationsOfState();
  void nonDimensionalizeAllInitialConditions();
  void nonDimensionalizeOneDimensionalProblem();
  int checkInputValuesNonDimensional();
  int checkInputValuesDimensional(map<int,SurfaceData*>& surfaceMap);
  int checkInputValuesEssentialBC();
  void checkInputValuesTurbulence();
  void checkInputValuesDefaultOutlet();
  int checkSolverValues(map<int,SurfaceData*>& surfaceMap);
  int checkInputValuesInitialConditions(InitialConditions &initialConditions,
                                        int fluidModelID);
  int checkInputValuesEquationOfState(FluidModelData &fluidModel, int fluidModelID);
  void nonDimensionalizeInitialConditions(InitialConditions &initialConditions);
  void nonDimensionalizeFluidModel(FluidModelData &fluidModel);
  int checkInputValuesSparseGrid(SparseGridData &sparseGrid);
  int checkInputValuesEmbeddedFramework();
  void printDebug();

  void setupOneDimensional();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <IoData.C>
#endif

#endif
