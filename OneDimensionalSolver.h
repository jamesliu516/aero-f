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

#include "PostFcn.h"

#include <fstream>

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

  double sscale[PostFcn::SSIZE];
  double vscale[PostFcn::SSIZE];

  char *solutions;
  char *scalars[PostFcn::SSIZE];
  char *vectors[PostFcn::VSIZE];

  void setupOutputFiles(IoData& iod);
  void setupFixes(IoData& iod);

  int* loctag;

 public:

  void EulerF(double t, SVec<double,5>& y,SVec<double,5>& k);
  void PhiF(double t, SVec<double,1>& y,SVec<double,1>& k);

  enum ReadMode { ModeU, ModePhi };
  template <int dimp,int dimLS>
  static void read1DSolution(IoData& iod, DistSVec<double,dimp>& Up, 
			     DistSVec<double,dimLS>* Phi,
			     FluidSelector* fluidSelector,
			     VarFcn* varFcn,
			     DistSVec<double,3>& X,
			     ReadMode mode) {
    // read 1D solution
    for (map<int, OneDimensionalInputData *>::iterator itr = iod.input.oneDimensionalInput.dataMap.begin();
	 itr != iod.input.oneDimensionalInput.dataMap.end(); ++itr) {
      
      std::fstream input;
      int sp = strlen(iod.input.prefix) + 1;
      char *filename = new char[sp + strlen(itr->second->file)];
      sprintf(filename, "%s%s", iod.input.prefix, itr->second->file);
      input.open(filename, fstream::in);
      //cout << filename << endl;
      if (!input.is_open()) {
	cout<<"*** Error: could not open 1D solution file "<<filename<<endl;
	exit(1);
      }
      
      input.ignore(256,'\n');
      input.ignore(2,' ');
      int numPoints = 0;
      input >> numPoints;
      //cout << " Num 1d points = " << numPoints << endl;
      double* x_1D = new double[numPoints];
      double* v_1D = new double[numPoints*5];/* rho, u, p, phi, T*/
      int* fids = new int[numPoints];
      
      for(int i=0; i<numPoints; i++){
	input >> x_1D[i] >> v_1D[i*5] >> v_1D[i*5+1] >> v_1D[i*5+2] >> v_1D[i*5+3]>> fids[i] >> v_1D[i*5+4];
	x_1D[i]    /= iod.ref.rv.length;
	v_1D[i*5] /= iod.ref.rv.density;
	v_1D[i*5+1] /= iod.ref.rv.velocity;
	v_1D[i*5+2] /= iod.ref.rv.pressure;
	//v_1D[i][3] /= iod.ref.rv.length;
	v_1D[i*5+4] /= iod.ref.rv.temperature;
      }
      
      input.close();

      int lsdim = 0;
      if (fluidSelector) {
	int a = 0;
	int fid_new = fids[a];
	if (itr->second->fluidRemap.dataMap.find(fids[a]) != itr->second->fluidRemap.dataMap.end())
	  fid_new = itr->second->fluidRemap.dataMap.find(fids[a])->second->newID;
	lsdim = fluidSelector->getLevelSetDim(0,fid_new);
      }
      
      // interpolation assuming 1D solution is centered on bubble_coord0
      double bubble_x0 = itr->second->x0;
      double bubble_y0 = itr->second->y0;
      double bubble_z0 = itr->second->z0;
      double max_distance = x_1D[numPoints-1]; 
#pragma omp parallel for
      for (int iSub=0; iSub<Up.numLocSub(); ++iSub) {
	SVec<double,dimp> &u(Up(iSub));
	SVec<double, 3> &x(X(iSub));
	
	double localRadius; int np;
	double localAlpha, velocity_r;
	double localV[5]; double localPhi;
	for(int i=0; i<u.size(); i++) {
	  localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
	  
	  // If the node is inside the sphere, set its values accordingly
	  if (localRadius < max_distance) {
	    
	    // Do a binary search to find the closest point whose r
	    // is less than or equal to localRadius
	    int a = numPoints/2;
	    int rmin = 0,rmax=numPoints-2;
	    while (!(x_1D[a] <= localRadius && x_1D[a+1] > localRadius)  ) {
	      if (x_1D[a] < localRadius)
		rmin = a;
	      else
		rmax = a;

	      if (a == rmin)
		a = (rmax+rmin)/2 + 1;
	      else
		a = (rmax+rmin)/2 ;
	    }
	    int fid_new = fids[a];
	    if (itr->second->fluidRemap.dataMap.find(fids[a]) != itr->second->fluidRemap.dataMap.end())
	      fid_new = itr->second->fluidRemap.dataMap.find(fids[a])->second->newID;

	    //cout << "fid_new = " << fid_new << " a = " << a << " x[] = {" << x_1D[a] << " " << x_1D[a+1] << "}" << endl;
	    double alpha = (localRadius-x_1D[a])/(x_1D[a+1]-x_1D[a]);
	    //cout << "alpha = " << alpha << endl;
	    if (mode == ModeU) {
	      localV[0] = v_1D[a*5]*(1.0-alpha)+v_1D[(a+1)*5]*(alpha);
	      localV[1] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(x[i][0]-bubble_x0)/max(localRadius,1.0e-8);
	      localV[2] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(x[i][1]-bubble_y0)/max(localRadius,1.0e-8);
	      localV[3] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(x[i][2]-bubble_z0)/max(localRadius,1.0e-8);
	      if (varFcn->getType(fid_new) == VarFcnBase::TAIT)
		localV[4] = v_1D[a*5+4]*(1.0-alpha)+v_1D[(a+1)*5+4]*(alpha);
	      else
		localV[4] = v_1D[a*5+2]*(1.0-alpha)+v_1D[(a+1)*5+2]*(alpha);
	      varFcn->primitiveToConservative(localV,u[i],fid_new);
	    } else {

	      
	      int lsdim=0;
	      if (fluidSelector) {
		(*Phi)(iSub)[i][lsdim] = v_1D[(a+1)*5+3];
	      }
	    }
	    
	  }
	  
	}
      }
    }
  }

  void resultsOutput(double time, int iteration);
  void restartOutput(double time, int iteration);

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

  static void load1DMesh(IoData& ioData,int& numPts,double* &meshPoints);

  template <int dim>
    void computeSlopes(SVec<double,dim>& VV, SVec<double,dim>& slopes,
		       Vec<int>& fid, bool crossInterface);

};

#endif

