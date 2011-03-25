#ifndef _ONE_DIMENSIONAL_SOLVER_H_
#define _ONE_DIMENSIONAL_SOLVER_H_


//------------------------------------------------------------------------------
// Created March 2010 by A. Rallu @Stanford University
//------------------------------------------------------------------------------
// This class is used to solve the one-dimensional Euler equations.
// They can be written in cartesian, cylindrical or spherical coordinates.
// All these cases are considered in this class.
// This class is geared toward solving two-phase one-dimensional Euler equations
// with a possible burn process of explosives.
//
// In order to allow some flexibility and ease of use, this class
// uses only some of the local-level classes of the rest of the AERO-F code
// (like SVec, FluxFcn, LocalRiemann, VarFcn). Hence, it does not use
// classes like DistTimeState or SpaceOperator. One reason is that these
// classes then always refer to Edges and Faces and Elems which are not
// available for the simulations intended with this class. (It would be possible
// to go through Edges, Faces, Elems and therefore use SpaceOperator and 
// DistTimeState and others, but that requires more work/time than I can 
// afford now).
//
//
// Note that state vectors have 5 coordinates. Only three of them are used.
// The coordinates [2] and [3] are unused but must be present in order to 
// be able to use the functions of VarFcn, FluxFcn, LocalRiemann (as are).
//
// Note: the spherical one-D Euler equations can be written as
// d(r^2 U)/dt + d(r^2 F(U))/dr = S2(U,r) = {0, 2*r*pressure, 0}
// or
// dU/dt + dF(U)/dx = S(U,r) = {-a*density*velocity/r, -a*density*velocity^2/r, -a*(density*energy+pressure)*velocity/r}
// where a = 2 for spherical coordinates (if a = 1, these are the cylindrical
// one-dimensional Euler equations)
// The first should require that the conserved quantities are r^2*U.
// The second requires that the source be integrated in a specific manner at r=0
//
//------------------------------------------------------------------------------
#include "Vector.h"
#include "FluidSelector.h"
#include "IoData.h"
#include "RefVal.h"
#include "RecFcn.h"

#include "RKIntegrator.h"
#include "ExactRiemannSolver.h"
#include "ProgrammedBurn.h"

class FluxFcn;
class VarFcn;
class LocalRiemann;
class OneDimensionalSourceTermBase;

class Domain;
//------------------------------------------------------------------------------

class OneDimensional {
  const static int dim = 5, dimLS = 1;
  OneDimensionalInfo::CoordinateType coordType;
  OneDimensionalInfo::VolumeType volumeType;

  int numPoints;
  double maxDistance;
  SVec<double,1> ctrlVol;
  // for cartesian, cylindrical and spherical one-D simulation with control surfaces
  SVec<double,1> ctrlSurf;
  SVec<double,1> X; // center of control volumes
  SVec<double,1> Y; // control surface locations

  double finalTime;
  double cfl;

  double BC[2][5];
  double BCphi[2];
  
  SVec<double,5> U;
  SVec<double,5> V;
  SVec<double,5> R;
  SVec<double,5> gradV;

  SVec<double,1> Phi;
  SVec<double,1> Phin;
  SVec<double,1> Rphi;
  SVec<double,1> gradPhi;
  Vec<int> fluidId;
  Vec<int> fluidIdn;

  FluxFcn **fluxFcn;
  VarFcn *varFcn;
  ExactRiemannSolver<5>* riemann;
  FluidSelector fluidSelector;
  // for cartesian, cylindrical and spherical one-D simulation with source term
  OneDimensionalSourceTermBase *source; 

  int frequency; //postprocessing output frequency
  char *outfile;
  RefVal refVal;

  // Riemann solutions at the interface
  // Needed for appropriate handling of phase change updates.
  SVec<double,5> Wr;

  SVec<double,5> Vslope;
  SVec<double,1> Phislope;

  Vec<int> riemannStatus;

  RecFcn* recFcn, *recFcnLS;

  RecFcn* createRecFcn(IoData &ioData);
  RecFcn* createRecFcnLS(IoData &ioData);

  RKIntegrator< SVec<double, 5> >* Vintegrator;
  RKIntegrator< SVec<double, 1> >* Phiintegrator;

  double time;

  SVec<double,dim> rupdate;
  Vec<double> weight;
  SVec<double,dim-2> interfacialWi;
  SVec<double,dim-2> interfacialWj;

  void loadSparseGrid(IoData&);

  SparseGridCluster* tabulationC;

  ProgrammedBurn* programmedBurn;

 public:

  void EulerF(double t, SVec<double,5>& y,SVec<double,5>& k);
  void PhiF(double t, SVec<double,1>& y,SVec<double,1>& k);

public:
  OneDimensional(int, double*, IoData &ioData, Domain *domain);
  ~OneDimensional();

  void spatialSetup();
  void temporalSetup();
  void stateInitialization(OneDimensionalInfo &data);
  void totalTimeIntegration();
  void computeTimeSteps(SVec<double,1> &timeSteps);
  void singleTimeIntegration(double dt);
  void computeEulerFluxes(SVec<double,5>& y);
  void computeLevelSetFluxes(SVec<double,1>& y);
  void resultsOutput(double time, int iteration);

  static void load1DMesh(IoData& ioData,int& numPts,double* &meshPoints);

  template <int dim>
    void computeSlopes(SVec<double,dim>& VV, SVec<double,dim>& slopes,
		       Vec<int>& fid, bool crossInterface);

};

#endif

