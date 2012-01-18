#include "OneDimensionalSolver.h"

#include <fstream>
#include <iostream>
using namespace std;
#include <cmath>

#include "FluxFcn.h"
#include "VarFcn.h"
#include "LocalRiemannDesc.h"
#include "FluidSelector.h"
#include "OneDimensionalSourceTerm.h"
#include "OneDimensionalInterpolator.h"
#include "RefVal.h"

#include "IoData.h"
#include "Domain.h"

//------------------------------------------------------------------------------
OneDimensional::OneDimensional(int np,double* mesh,IoData &ioData, Domain *domain) : 
  numPoints(np),
  ctrlVol(numPoints), ctrlSurf(numPoints+1), X(numPoints), Y(numPoints+1), cutCellStatus(numPoints+1),
  U(numPoints), V(numPoints), R(numPoints),
  gradV(numPoints), Phi(numPoints), Rphi(numPoints), gradPhi(numPoints),
  fluidId(numPoints),fluidIdn(numPoints),
  fluidSelector(ioData.eqs.numPhase,ioData), refVal(ioData.ref.rv),
  Wr(numPoints),Vslope(numPoints),Phislope(numPoints),
  rupdate(numPoints), weight(numPoints), interfacialWi(numPoints),
  interfacialWj(numPoints), riemannStatus(numPoints), Phin(numPoints),
  programmedBurn(NULL)
{
  // equation modelling
  coordType  = ioData.oneDimensionalInfo.coordType;
  volumeType = ioData.oneDimensionalInfo.volumeType;
  interfaceTreatment = 0;

  // time and space domain definition
  maxDistance = mesh[np-1];
  finalTime = ioData.ts.maxTime;
  cfl = ioData.ts.cfl0;

  // Copy 1D mesh to X
  X = 0.0;
  for (int i = 0; i < np; ++i)
    X[i][0] = mesh[i];
  delete [] mesh;

  bubbleRadiusFile = new char[256];
  

  // output
  frequency = ioData.output.transient.frequency;

  // necessary for computation: varFcn to compute different state quantities
  //                            fluxFcn to compute fluxes at interfaces
  //                            riemann to solve riemann problem between two fluids
  //                            source to compute source term(s) for spherical and cylindrical problems
  varFcn = new VarFcn(ioData);

  fluxFcn = new FluxFcn *[3];
  fluxFcn[0] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn);
  fluxFcn[1] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn);
  fluxFcn[2] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn);

  //riemann = new LocalRiemannGfmparGasJWL(varFcn,0,1,0,MultiFluidData::RK2);
  //riemann = new LocalRiemannGfmpGasJWL(varFcn,0,1);

  source = 0;
  if(volumeType == OneDimensionalInfo::CONSTANT_VOLUME){
    if(coordType == OneDimensionalInfo::CYLINDRICAL ||
       coordType == OneDimensionalInfo::SPHERICAL)
      source = new OneDimensionalSourceTerm();
  }else if(volumeType == OneDimensionalInfo::REAL_VOLUME){
  }

  strcpy(bubbleRadiusFile,"");

  recFcn = createRecFcn(ioData);
  recFcnLS = createRecFcnLS(ioData);

  if (ioData.ts.type != TsData::EXPLICIT) {
    fprintf(stderr,"Only explcit integration available for the 1D solver!\n");
    exit(1);
  }

  switch (ioData.ts.expl.type) {

  case ExplicitData::RUNGE_KUTTA_4:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK4, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK4, numPoints);
    break;
  case ExplicitData::RUNGE_KUTTA_2:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK2, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK2, numPoints);
    break;
  case ExplicitData::FORWARD_EULER:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::FE, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::FE, numPoints);
    break;

  default:
    fprintf(stderr,"Unavailable explicit integration for the 1D solver!\n");
    exit(1);
    break;
  }

  loadSparseGrid(ioData);

  riemann = new ExactRiemannSolver<5>(ioData,rupdate,weight, interfacialWi,
				      interfacialWj, varFcn,
				      tabulationC);
  
  if (ioData.oneDimensionalInfo.programmedBurn.unburnedEOS >= 0) {
    programmedBurn = new ProgrammedBurn(ioData,&this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }

  setupOutputFiles(ioData);

  setupFixes(ioData);

  cutCellStatus = 0;


  if (ioData.schemes.ns.dissipation == SchemeData::SIXTH_ORDER) {
    isSixthOrder = true;
  } else
    isSixthOrder = false;

  beta = ioData.schemes.ns.beta;

  if (ioData.mf.interfaceTreatment == MultiFluidData::SECONDORDER)
    interfaceTreatment = 1;

  levelSetMethod = 0;
 
  interfaceExtrapolation = 0;
  if (ioData.mf.interfaceExtrapolation == MultiFluidData::EXTRAPOLATIONSECONDORDER)
    interfaceExtrapolation = 1;
  
  
  if (ioData.mf.levelSetMethod == MultiFluidData::HJWENO)
    levelSetMethod = 1;
  else if (ioData.mf.levelSetMethod == MultiFluidData::SCALAR/* || ioData.mf.levelSetMethod == MultiFluidData::CONSERVATIVE*/)
    levelSetMethod = 2;
 
  int sto = ioData.oneDimensionalInfo.sourceTermOrder; 
  if(volumeType == OneDimensionalInfo::CONSTANT_VOLUME){
    if(coordType == OneDimensionalInfo::CYLINDRICAL)
      source->initialize(1.0,sto,X);
    else if (coordType == OneDimensionalInfo::SPHERICAL)
      source->initialize(2.0,sto,X);
  }
}
//------------------------------------------------------------------------------
OneDimensional::~OneDimensional(){

  delete varFcn;
  delete riemann;
  for(int i=0; i<3; i++) delete fluxFcn[i];
  delete [] fluxFcn;
  
  if (source)
    delete source;

  if (tabulationC)
    delete tabulationC;

}

void OneDimensional::setupFixes(IoData& ioData) {


  double spheres[20][4];
  double boxes[20][2][3];
  double cones[20][2][4];
  int j, nspheres = 0, nboxes = 0, ncones = 0;
  for (j=0; j<ioData.schemes.fixes.num; ++j) {
    if (ioData.schemes.fixes.spheres[j]->r > 0.0) {
      spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
      spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
      spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
      spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
      ++nspheres;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	spheres[nspheres][0] = - ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = - ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = - ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
    }
    if (ioData.schemes.fixes.boxes[j]->x0 < ioData.schemes.fixes.boxes[j]->x1) {
      boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
      boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
      boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
      boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
      boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
      boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
      ++nboxes;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	boxes[nboxes][0][0] = -ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = -ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = -ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = -ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
     if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = -ioData.schemes.fixes.boxes[j]->z1;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = -ioData.schemes.fixes.boxes[j]->z0;
	++nboxes;
      }
    }
    if (ioData.schemes.fixes.cones[j]->r0 >= 0.0 && ioData.schemes.fixes.cones[j]->r1 >= 0.0) {
      cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
      cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
      cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
      cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
      cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
      cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
      cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
      cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
      ++ncones;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
        cones[ncones][0][0] = -ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = -ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = -ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = -ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = -ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = -ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
    }
  }

  loctag = new int[numPoints];
  memset(loctag,0,sizeof(int)*numPoints);

  if (nspheres > 0 || nboxes > 0 || ncones > 0) {
    for (j=0; j<nspheres; ++j)
      printf( "*** Warning: set the gradients to zero in [(%g, %g, %g), %g]\n",
		  spheres[j][0], spheres[j][1], spheres[j][2], spheres[j][3]);
    for (j=0; j<nboxes; ++j)
      printf( "*** Warning: set the gradients to zero in [(%g, %g, %g), (%g, %g, %g)]\n",
		  boxes[j][0][0], boxes[j][0][1], boxes[j][0][2],
		  boxes[j][1][0], boxes[j][1][1], boxes[j][1][2]);

    for (j=0; j<ncones; ++j)
      printf( "*** Warning: set the gradients to zero in cone [(%g, %g, %g), %g; (%g, %g, %g), %g]\n",
                  cones[j][0][0], cones[j][0][1], cones[j][0][2], cones[j][0][3],
                  cones[j][1][0], cones[j][1][1], cones[j][1][2], cones[j][1][3]);

    for (int i=0; i<numPoints; ++i) {
      double x0[3] = {X[i][0],0.0,0.0};
      loctag[i] = 0;
      for (j=0; j<nspheres; ++j) {
	double r = sqrt( (x0[0] - spheres[j][0])*(x0[0] - spheres[j][0]) +
			 (x0[1] - spheres[j][1])*(x0[1] - spheres[j][1]) +
			 (x0[2] - spheres[j][2])*(x0[2] - spheres[j][2]) );
	if (r <= spheres[j][3])
	  loctag[i] = 1;
      }
      for (j=0; j<nboxes; ++j) {
	if ((x0[0] >= boxes[j][0][0]) && (x0[0] <= boxes[j][1][0]) &&
	    (x0[1] >= boxes[j][0][1]) && (x0[1] <= boxes[j][1][1]) &&
	    (x0[2] >= boxes[j][0][2]) && (x0[2] <= boxes[j][1][2]))
	    loctag[i] = 1;
      }
      for (j=0; j<ncones; ++j)  {
	Vec3D dr(cones[j][1][0]-cones[j][0][0], cones[j][1][1]-cones[j][0][1], cones[j][1][2]-cones[j][0][2]);
	double height = dr.norm();
	dr /= height;
	Vec3D xp;
	Vec3D pr0(x0[0]-cones[j][0][0], x0[1]-cones[j][0][1], x0[2]-cones[j][0][2]);
	double h = pr0*dr;
	if (h >= 0.0 && h <= height)  {
	  xp = pr0 - (h*dr);
	  double r = cones[j][0][3] + (cones[j][1][3]-cones[j][0][3]) * h / height;
	  if (xp.norm() < r)
	    loctag[i] = 1;
          }
      }
    }
  }
}

void OneDimensional::setupOutputFiles(IoData& iod) {

  int i;


  int sp = strlen(iod.output.transient.prefix) + 1;
  int spr = strlen(iod.output.restart.prefix) + 1;


  if (iod.output.restart.solutions[0] != 0) {
    outfile = new char[spr + strlen(iod.output.restart.solutions)];
    sprintf(outfile, "%s%s", 
	    iod.output.restart.prefix, iod.output.restart.solutions);
  }

  for (i=0; i<PostFcn::SSIZE; ++i) {
    sscale[i] = 1.0;
    scalars[i] = 0;
  }

  for (i=0; i<PostFcn::VSIZE; ++i) {
    vscale[i] = 1.0;
    vectors[i] = 0;
  }

  if (iod.output.transient.density[0] != 0) {
    sscale[PostFcn::DENSITY] = iod.ref.rv.density;
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.density);
  }
  /*if (iod.output.transient.tavdensity[0] != 0) {
    avsscale[PostFcn::DENSITYAVG] = iod.ref.rv.density;
    avscalars[PostFcn::DENSITYAVG] = new char[sp + strlen(iod.output.transient.tavdensity)];
    sprintf(avscalars[PostFcn::DENSITYAVG], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.tavdensity);
	    }*/
  if (iod.output.transient.mach[0] != 0) {
    scalars[PostFcn::MACH] = new char[sp + strlen(iod.output.transient.mach)];
    sprintf(scalars[PostFcn::MACH], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.mach);
  }
  /*if (iod.output.transient.wtmach[0] != 0) {
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
	    }*/
  if (iod.output.transient.pressure[0] != 0) {
    sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.pressure);
  }
  /*if (iod.output.transient.diffpressure[0] != 0) {
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
	    }*/
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
//    sscale[PostFcn::TEMPERATURE] = 1;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.temperature);
  }
  /*if (iod.output.transient.tavtemperature[0] != 0) {
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
	    }*/
  //  if (iod.output.transient.dplus[0] != 0) {
  //#if defined(HEAT_FLUX)
    /*
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    double dT = iod.bc.wall.temperature - 1.0 / (gam*(gam-1.0)*iod.bc.inlet.mach*iod.bc.inlet.mach);
    sscale[PostFcn::DELTA_PLUS] = iod.ref.reynolds * iod.eqs.thermalCondModel.prandtl / (gam * dT);
    */
  //    sscale[PostFcn::DELTA_PLUS] = iod.ref.rv.tpower / (iod.ref.length*iod.ref.length);
  //#endif
  //    scalars[PostFcn::DELTA_PLUS] = new char[sp + strlen(iod.output.transient.dplus)];
  //  sprintf(scalars[PostFcn::DELTA_PLUS], "%s%s", 
  //	    iod.output.transient.prefix, iod.output.transient.dplus);
  //}
  /*if (iod.output.transient.sfric[0] != 0) {
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
	    }*/
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
  /*if (iod.output.transient.controlvolume[0] != 0) {
    sscale[PostFcn::CONTROL_VOLUME] = iod.ref.rv.length * iod.ref.rv.length * iod.ref.rv.length;
    scalars[PostFcn::CONTROL_VOLUME] = new char[sp + strlen(iod.output.transient.controlvolume)];
    sprintf(scalars[PostFcn::CONTROL_VOLUME], "%s%s",
            iod.output.transient.prefix, iod.output.transient.controlvolume);
	    }*/
  if (iod.output.transient.velocity[0] != 0) {
    vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocity);
  }
  /*if (iod.output.transient.tavvelocity[0] != 0) {
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
  hydrodynamicforces = 0;*/

  if (iod.output.transient.bubbleRadius[0] != 0) {
    bubbleRadiusFile = new char[sp + strlen(iod.output.transient.bubbleRadius)];
    sprintf(bubbleRadiusFile, "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.bubbleRadius);
  }
}

//------------------------------------------------------------------------------
void OneDimensional::load1DMesh(IoData& ioData,int& numPts,double* &meshPoints) {

  if (ioData.input.geometry[0] != 0) {
    char mesh1d[256];
    
    sprintf(mesh1d,"%s%s",ioData.input.prefix,ioData.input.geometry);
    FILE* fin = fopen(mesh1d,"r");
    int n = fscanf(fin, "%i",&numPts);
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i) {

      int m = fscanf(fin,"%lf",&meshPoints[i]);
    }
  } else {

    numPts = ioData.oneDimensionalInfo.numPoints;
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i)
      meshPoints[i] = (double)i / (numPts-1)*ioData.oneDimensionalInfo.maxDistance;// + 0.25*((double)rand()/RAND_MAX-0.5)/(numPts-1)*ioData.oneDimensionalInfo.maxDistance;
  }
}
//------------------------------------------------------------------------------
void OneDimensional::spatialSetup(){

  //cout << "computing cell boundaries" <<endl;
  Y = 0.0;
  Y[0][0] = X[0][0];
  for(int i=1; i<numPoints; i++) Y[i][0] = 0.5*(X[i-1][0]+X[i][0]);
  Y[numPoints][0] = X[numPoints-1][0];

  //cout << "computing control volumes"<<endl;
  // computation of control volumes, for volumeType == CONSTANT_VOLUME
  ctrlVol[0][0] = 0.5*(X[1][0]-X[0][0]);
  for(int i=1; i<numPoints-1; i++)
    ctrlVol[i][0] = 0.5*(X[i+1][0]-X[i-1][0]);
  ctrlVol[numPoints-1][0] = 0.5*(X[numPoints-1][0]-X[numPoints-2][0]);

  // computation of control surfaces, for volumeType == CONSTANT_VOLUME
  //cout << "computing control surfaces"<<endl;
  for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = 1.0;

  // in case a real Finite Volume approach is considered
  // with real cylindrical/spherical volumes, the control volumes and
  // surfaces must be computed differently as follows:
  if(volumeType == OneDimensionalInfo::REAL_VOLUME){

    if(coordType == OneDimensionalInfo::SPHERICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = (Y[i+1][0]-Y[i][0])*
                        (Y[i+1][0]*Y[i+1][0]+Y[i+1][0]*Y[i][0]+Y[i][0]*Y[i][0])/3.0;
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0]*Y[i][0];

    }else if(coordType == OneDimensionalInfo::CYLINDRICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = 0.5*(Y[i+1][0]-Y[i][0])*(Y[i+1][0]+Y[i][0]);
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0];
    }
  }

   
}
//------------------------------------------------------------------------------
void OneDimensional::temporalSetup(){
}
//------------------------------------------------------------------------------
void OneDimensional::stateInitialization(OneDimensionalInfo &data){
  V = 0.0;
  // initialize V
  if (data.mode == OneDimensionalInfo::CONVTEST1) {
    data.pressure1 = data.pressure2;
  }
  std::cout << "p2 = " << data.pressure2 << std::endl;
  interfaceLocation = data.interfacePosition;
  for(int i=0; i<numPoints; i++){
    if(X[i][0]<data.interfacePosition){
      if (varFcn->getType(1) != VarFcnBase::TAIT) {
	//if (X[i][0] < 3.8) {
	  V[i][0] = data.density1;
	  V[i][1] = data.velocity1;
	  V[i][4] = data.pressure1;
	  /*} else {
	  V[i][0] = data.density1;
	  V[i][1] = data.velocity1;
	  V[i][4] = data.pressure1*10.0;
	  }*/
      } else {
	V[i][0] = data.density1;
	V[i][1] = data.velocity1;
	V[i][4] = data.temperature1;
      }
    }else{
      if (varFcn->getType(0) != VarFcnBase::TAIT) {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.pressure2;
      } else {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.temperature2;
      }
    }

    if (data.mode == OneDimensionalInfo::CONVTEST1) {
      V[i][0] = 1.0;
      if (fabs(X[i][0]-data.interfacePosition) < 0.2)
	V[i][4] = //data.pressure2+100.0*pow(cos(3.14159265358979323846/2.0*(X[i][0]-data.interfacePosition)/0.2),2.0);
	  100000000.0*pow((X[i][0]-data.interfacePosition)+0.2,4.0)*pow((X[i][0]-data.interfacePosition)-0.2,4.0)+data.pressure2;
      else
	V[i][4] = data.pressure2;
      V[i][1] = 0.0;
    }
  }

  isSinglePhase = varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1));

  //cout << data.density1 << " " << data.velocity1 << " " << data.pressure1 << " " << data.temperature1 << endl;
  //cout << data.density2 << " " << data.velocity2 << " " << data.pressure2 << " " << data.temperature2 << endl;
  

  // initialize Phi
  if (levelSetMethod == 0) {
    for(int i=0; i<numPoints; i++)
      Phi[i][0]  = -V[i][0]*(X[i][0]-data.interfacePosition);
  } else {
    for(int i=0; i<numPoints; i++)
      Phi[i][0]  = -(X[i][0]-data.interfacePosition);
  }

  // initialize fluidId and U
  fluidSelector.getFluidId(fluidId,Phi);
  varFcn->primitiveToConservative(V,U,&fluidId);

  fluidIdn = fluidId;

  // compute boundary states
  double temp0[5] = {data.density1, data.velocity1, 0.0, 0.0, data.pressure1};
  double temp1[5] = {data.density2, data.velocity2, 0.0, 0.0, data.pressure2};
  varFcn->primitiveToConservative(temp0,BC[0],fluidId[0]);
  varFcn->primitiveToConservative(temp1,BC[1],fluidId[numPoints-1]);
  BCphi[0] = Phi[0][0]/V[0][0];
  BCphi[1] = Phi[numPoints-1][0]/V[numPoints-1][0];

  // output
  cout<<"**primitive boundary conditions are:"<<endl;
  cout<<"*    "<<temp0[0]<<" "<<temp0[1]<<" "<<temp0[2]<<" "<<temp0[3]<<" "<<temp0[4]<<endl;
  cout<<"*    "<<temp1[0]<<" "<<temp1[1]<<" "<<temp1[2]<<" "<<temp1[3]<<" "<<temp1[4]<<endl;
  cout<<"**conservative boundary conditions are:"<<endl;
  cout<<"*    "<<BC[0][0]<<" "<<BC[0][1]<<" "<<BC[0][2]<<" "<<BC[0][3]<<" "<<BC[0][4]<<endl;
  cout<<"*    "<<BC[1][0]<<" "<<BC[1][1]<<" "<<BC[1][2]<<" "<<BC[1][3]<<" "<<BC[1][4]<<endl;
  
}
//------------------------------------------------------------------------------
void OneDimensional::totalTimeIntegration(){

  // messages
  cout<<"***************************************************************"<<endl;
  cout<<"***  ctrlVol[0] = "<<ctrlVol[0][0]<<" --- ctrlVol[end] = "<<ctrlVol[numPoints-1][0]<<endl;
  cout<<"***  ctrlSur[0] = "<<ctrlSurf[0][0]<<" --- ctrlSur[end] = "<<ctrlSurf[numPoints][0]<<endl;
  cout<<"***************************************************************"<<endl;
  SVec<double,1> timeSteps(numPoints);

  time = 0.0;
  double dt   = 0.0;
  int iteration = 0;

  resultsOutput(time,iteration);

  while(time<finalTime){
    // determine how much to advance in time
    computeTimeSteps(timeSteps);
    dt = cfl*timeSteps.min();

    if (programmedBurn)
      programmedBurn->setCurrentTime(time,varFcn, U,fluidId,fluidIdn);

    if(time+dt>finalTime) dt = finalTime-time;
    if(iteration % frequency == 0)
      cout <<"*** Iteration " << iteration <<": Time = "<<time*refVal.time<<", and dt = "<<dt*refVal.time<<endl;
    iteration++;

    // advance one iteration
    singleTimeIntegration(dt);
    time += dt;

    if(iteration % frequency == 0)
      resultsOutput(time,iteration);
  }
  resultsOutput(time,iteration);
  restartOutput(time,iteration);

}
//------------------------------------------------------------------------------
void OneDimensional::computeTimeSteps(SVec<double,1> &timeSteps){

  // very crude CFL law
  double c = 0.0;
  varFcn->conservativeToPrimitive(U,V,&fluidId);

  for(int i=0; i<numPoints; i++){
    c = varFcn->computeSoundSpeed(V[i],fluidId[i]);
    if (c == 0.0)
      std::cout << "c = " << c <<  " " << fluidId[i] <<  std::endl;
    timeSteps[i][0] = 0.5*(Y[i+1][0]-Y[i][0])/c;
  }

}

void OneDimensional::levelSetDerivative(double t0, Vec<double>& phi, Vec<double>& k) {
  
  double u;
  if (interfaceTreatment == 1) {
   
    /*for(int i=0; i<numPoints-1; i++){
      if (cutCellStatus[i] == 1)
	V[i][1] = 0.5*(V[i-1][1]+V[i+1][1]);
	}*/
    for(int i=0; i<numPoints-1; i++){
      if (cutCellStatus[i] == 1) {
	double xi = phi[0];
	double ul = V[i-1][1]*(xi-X[i-2][0])/(X[i-1][0]-X[i-2][0])+V[i-2][1]*(xi-X[i-1][0])/(X[i-2][0]-X[i-1][0]);
	double ur = V[i+1][1]*(xi-X[i+2][0])/(X[i+1][0]-X[i+2][0])+V[i+2][1]*(xi-X[i+1][0])/(X[i+2][0]-X[i+1][0]);
	k[0] = (ul+ur)*0.5;
      }
      
    }
  } else {
    OneDimensionalInterpolator::Interpolate<5>(V,u , 1,
					       X, phi[0] ,2);
    
    k[0] = u;
    //std::cout << u << std::endl;
  }
}

//------------------------------------------------------------------------------
void OneDimensional::singleTimeIntegration(double dt){
// for now, assume forward Euler

  double Vtemp[5];
  
  if (levelSetMethod == 1 || levelSetMethod == 0) {

    fluidSelector.getFluidId(fluidId,Phi);
  }

  //std::cout << interfaceTreatment << " " << interfaceExtrapolation << std::endl;
  if (interfaceTreatment == 1) {
    //varFcn->conservativeToPrimitive(U,V,&fluidId);
    for(int i=0; i<numPoints-1; i++){
      
      if (Y[i][0] < interfaceLocation && Y[i+1][0] >= interfaceLocation) {
	for (int k = 0; k < 5; ++k) {
	  if (interfaceLocation > X[i][0])
	    V[i][k] = (2.0*V[i-2][k]-V[i-1][k]);
	  else
	    V[i][k] = (2.0*V[i+1][k]-V[i+2][k]);
	}
	cutCellStatus[i] = 1;
      }
      else {
	if (cutCellStatus[i] == 1) {
	  if (interfaceExtrapolation == 0) {
	    memcpy(V[i],Wr[i],sizeof(double)*5);
	    for (int k = 0; k < 5; ++k) {
	      std::cout << V[i][k] << " ";
	    }
	  }
	  else if (interfaceExtrapolation == 1) {

	    // Linear extrapolation to populate the ghost state
	    int j = i+1;
	    if (fluidId[i-1] == fluidId[i])
	      j = i-1;
	    double xi = (X[j][0]+X[i][0])*0.5;
	    for (int k = 0; k < 5; ++k) {
	      V[i][k] = (X[i][0]-X[j][0])/(xi-X[j][0])*Wr[j][k]-
		(X[i][0]-xi)/(xi-X[j][0])*V[j][k];
	      std::cout << V[i][k] << " ";
	    }
	  }
	  std::cout << std::endl;	
	}
	cutCellStatus[i] = 0;
      }
    }
    varFcn->primitiveToConservative(V,U,&fluidId);
  }
  
  riemannStatus = 0;
  //std::cout << "U0 = " << U*U << std::endl;

  Vintegrator->integrate(this,&OneDimensional::EulerF,
			 U,time,dt);
  //R = 0.0;
  //computeEulerFluxes();

  //std::cout << "U1 = " << U*U << std::endl;

  if (levelSetMethod == 1 || levelSetMethod == 0) {
    Phin = Phi;
    Phiintegrator->integrate(this,&OneDimensional::PhiF,
			     Phi,time,dt);
    //Rphi = 0.0;
    //computeLevelSetFluxes();
    
    //preliminary update of U and Phi
    for(int i=0; i<numPoints; i++){
      if (Phi[i][0]*Phin[i][0] < 0.0 &&
	  !riemannStatus[i])
	Phi[i][0] = Phin[i][0];
    }

    // store previous primitive with old fluidId
    varFcn->conservativeToPrimitive(U,V,&fluidId);
    
    // update fluidId
    fluidIdn = fluidId;
    fluidSelector.getFluidId(fluidId,Phi);
  } else {

    fluidIdn = fluidId;
    Vec<double> phi(1);
    phi[0] = interfaceLocation;
    varFcn->conservativeToPrimitive(U,V,&fluidId);
    RKIntegrator<Vec<double> > phiI( RKIntegrator< Vec<double> >::RK4, 1);
    phiI.integrate(this, &OneDimensional::levelSetDerivative,
		   phi,time,dt);
    interfaceLocation = phi[0];
    
    //fprintf(stderr,"int. loc = %lf, x = %lf\n", interfaceLocation, X[75][0]);
    for(int i=0; i<numPoints; i++){
      if (X[i][0] < interfaceLocation)
	fluidId[i] = 1;
      else
	fluidId[i] = 0;
    }

  }

  for(int i=0; i<numPoints; i++){

    if (fluidId[i] != fluidIdn[i]&& !varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1))) { // Phase change
      
      if (!riemannStatus[i])
	std::cout << "Have a problem!" << std::endl;
      //fprintf(stderr,"Node %d changes phase!\n", i);
      //fprintf(stderr,"Wr[i] = %lf %lf %lf %lf %lf\n", Wr[i][0],Wr[i][1],Wr[i][2],Wr[i][3],Wr[i][4]);
      if (interfaceTreatment == 0) {
	if (interfaceExtrapolation == 0)
	  memcpy(V[i],Wr[i],sizeof(double)*5);
	else {

	  // Linear extrapolation to populate the ghost state
	  int j = i+1;
	  if (fluidId[i-1] == fluidId[i])
	    j = i-1;
	  double xi = (X[j][0]+X[i][0])*0.5;
	  for (int k = 0; k < 5; ++k) {
	    V[i][k] = (X[i][0]-X[j][0])/(xi-X[j][0])*Wr[i][k]-
	      (X[i][0]-xi)/(xi-X[j][0])*V[j][k];
	  }
	}
      }	
    }
  }

  double interfaceLocation;
  // initialize Phi
  for(int i=0; i<numPoints-1; i++) {

    if (Phi[i+1][0] < 0.0) {
      if (levelSetMethod == 0)
	interfaceLocation = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
      else if (levelSetMethod == 1)
	interfaceLocation = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
      break;
    }
  }

    
  if (levelSetMethod == 0) {
    /*for(int i=0; i<numPoints-1; i++) {
      Phi[i][0] = V[i][0]*(interfaceLocation-X[i][0]);
      }*/
  } else if (levelSetMethod == 1) {
    for(int i=0; i<numPoints-1; i++) {
      Phi[i][0] = (interfaceLocation-X[i][0]);
    }

  }
  

  //update PhaseChange with new fluidId
  varFcn->primitiveToConservative(V,U,&fluidId);

  // Check solution (clip pressure, that is), if necessary
  if (varFcn->doVerification()) {
    for(int i=0; i<numPoints; i++){
      varFcn->conservativeToPrimitiveVerification(i+1, U[i], Vtemp, fluidId[i]);
    }
  }
   
  if (programmedBurn) {

    programmedBurn->setFluidIds(time, fluidId,U);
  }

}   

void OneDimensional::EulerF(double t, SVec<double,5>& y,SVec<double,5>& k) {

  R = 0.0;
  computeEulerFluxes(y);
  for(int i=0; i<numPoints; i++){
    for(int idim=0; idim<dim; idim++)
      k[i][idim] = -R[i][idim] / ctrlVol[i][0];
  }
}

class EulerSource {

public:

  EulerSource(VarFcn* _vf,Vec<int>& _fid) :  vf(_vf), fid(_fid) { }

  ~EulerSource() { }

  void compute(double* U, double* f,int i) {

    double V[5];
    vf->conservativeToPrimitive(U,V,fid[i]);
    f[0] = U[1];
    f[1] = U[1]*U[1]/U[0];
    f[4] = (U[4]+vf->getPressure(V, fid[i]))*(U[1]/U[0]);
    f[2] = f[3] = 0.0;
  }

  VarFcn* vf;
  Vec<int>& fid;
};

//------------------------------------------------------------------------------
void OneDimensional::computeEulerFluxes(SVec<double,5>& y){

  double normal[3] = {1.0, 0.0, 0.0};
  double length = 1.0;
  double normalVel = 0.0; // no ALE
  double flux[dim];
  double Udummy[dim],Vtemp[dim];
  int i,j,k;
  
  for(int i=0; i<numPoints; i++){
    if (y[i][0] < 0.0)
      std::cout << "Error: node " << i << " has negative density " <<
            y[i][0] << "; fid = " << fluidId[i] << std::endl;
  }

  // Check solution (clip pressure, that is), if necessary
  if (varFcn->doVerification()) {
    for(int i=0; i<numPoints; i++){
      varFcn->conservativeToPrimitiveVerification(i+1, y[i], Vtemp, fluidId[i]);
    }
  } 

  varFcn->conservativeToPrimitive(y,V,&fluidId);
  computeSlopes(V,Vslope,fluidId,!varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1)));
  
  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    i = iEdge;
    j = iEdge+1;
  
    double Vi[dim*2],Vj[dim*2],Vsi[dim],Vsj[dim],VslopeI[dim],VslopeJ[dim];
    if (!isSixthOrder || !(i > 0 && i < numPoints-3 && (fluidId[i] == fluidId[i+1] && 
				     fluidId[i] == fluidId[i+2] && fluidId[i] == fluidId[i-1] || isSinglePhase))) {
      
      for (k = 0; k < dim; ++k) {
	VslopeI[k] = Vslope[i][k];
	VslopeJ[k] = Vslope[j][k];
      }
    } else {

      for (k = 0; k < dim; ++k) {
	VslopeI[k] = (V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0])+
	  (-1.0/30.0)*((V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])-2.0*(V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Vslope[i-1][k]-2.0*Vslope[i][k]+Vslope[j][k])/beta;
	
	  VslopeJ[k] = (V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])+
	  (-1.0/30.0)*((V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])-2.0*(V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Vslope[j+1][k]-2.0*Vslope[j][k]+Vslope[i][k])/beta;
	
	VslopeI[k] *= 0.5;
	VslopeJ[k] *= 0.5;
	
      }
      
    }

    for (k = 0; k < dim; ++k) {
      Vsi[k] = VslopeI[k]*(X[j][0]-X[i][0]);
      Vsj[k] = VslopeJ[k]*(X[j][0]-X[i][0]);
			  
    }
    
    // One dimensional (source) term
    /*if(volumeType == OneDimensionalInfo::REAL_VOLUME){
      
      if(coordType == OneDimensionalInfo::SPHERICAL) {

	double pi = varFcn->getPressure(Vi,fluidId[i]),pj = varFcn->getPressure(Vj,fluidId[j]);
	R[i][1] -= pi*ctrlSurf[iEdge+1][0];
	R[j][1] += pj*ctrlSurf[iEdge+1][0];
      }
      }*/

    length = (X[j][0]-X[i][0]);

    recFcn->compute(V[i], Vsi,V[j], Vsj, Vi, Vj);
    
    //std::cout << "Hello" << std::endl;
    varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vi);
    varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vj);
    for (k = 0; k < dim; ++k) {
      Vi[k+dim] = V[i][k];
      Vj[k+dim] = V[j][k];
    }

    if (interfaceTreatment == 0) {
      
      if(fluidId[i] == fluidId[j]|| isSinglePhase){

        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Vj, flux, fluidId[i]);
        for(int k=0; k<dim; ++k) {
          R[i][k] += ctrlSurf[iEdge+1][0]*flux[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*flux[k];
        }
      }else{
        double gradphi[3] = {1.0, 0.0, 0.0};
        double Wi[2*dim], Wj[2*dim];
        double Wir[2*dim], Wjr[2*dim];
        double Vir[dim],Vjr[dim];
        double wi, wj;
        double dx[3] = {X[j][0]-X[i][0], 0.0, 0.0};
        int iteration = 0;
        double fluxi[dim], fluxj[dim];

	memcpy(Vir, Vi, sizeof(double)*dim);
	memcpy(Vjr, Vj, sizeof(double)*dim);
	

        riemann->computeRiemannSolution(Vir,Vjr,fluidId[i],fluidId[j],gradphi,varFcn,
	  			      Wir,Wjr,i,j,i,dx);

	memcpy(Wi, Wir, sizeof(double)*dim);
	memcpy(Wj, Wjr, sizeof(double)*dim);
      
        memcpy(Wr[i], Wj, sizeof(double)*5);
        memcpy(Wr[j], Wi, sizeof(double)*5);
        riemannStatus[i] = riemannStatus[j] = 1;

        //fprintf(stderr,"Edge %d-%d crosses interface\n",i,j);

        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Wi, fluxi, fluidId[i]);
        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wj, Vj, fluxj, fluidId[j]);
        for (int k=0; k<dim; k++){
          R[i][k] += ctrlSurf[iEdge+1][0]*fluxi[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*fluxj[k];
        }

      }
    } else if (interfaceTreatment == 1) {

      if(cutCellStatus[i] == 0 && cutCellStatus[j] == 0){

        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Vj, flux, fluidId[i]);
        for(int k=0; k<dim; ++k) {
          R[i][k] += ctrlSurf[iEdge+1][0]*flux[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*flux[k];
        }
      }else{
        double gradphi[3] = {1.0, 0.0, 0.0};
        double Wi[2*dim], Wj[2*dim];
        double Wir[2*dim], Wjr[2*dim];
        double Vir[dim],Vjr[dim];
        double wi, wj;
        double dx[3] = {X[j][0]-X[i][0], 0.0, 0.0};
        int iteration = 0;
        double fluxi[dim], fluxj[dim];
        int I,J;
        if (cutCellStatus[j] == 1) {

          if (interfaceExtrapolation == 1) {
	    //std::cout << "\t" << X[i-1][0] << " " << X[i][0] << " " << X[j][0] << " " << X[j+1][0] << " " << X[j+2][0] << " " << interfaceLocation << std::endl;
	    for (int k = 0; k < dim; ++k) {

	      Vir[k] = (interfaceLocation-X[i-1][0])/(X[i][0]-X[i-1][0])*V[i][k] - (interfaceLocation-X[i][0])/(X[i][0]-X[i-1][0])*V[i-1][k];
	      Vjr[k] = (interfaceLocation-X[j+2][0])/(X[j+1][0]-X[j+2][0])*V[j+1][k] - (interfaceLocation-X[j+1][0])/(X[j+1][0]-X[j+2][0])*V[j+2][k];
	      //std::cout << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j+1][k] << " " << V[i-1][k] << " " << V[j+2][k] << std::endl;
	    }
          } else {
	    memcpy(Vir, V[i], sizeof(double)*dim);
	    memcpy(Vjr, V[j+1], sizeof(double)*dim);
	    for (int k = 0; k < dim; ++k) {
	      //std::cout << i << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
	    //std::cout << "p = " << varFcn->getPressure(Vjr,fluidId[j+1]) << std::endl;
	  }
        } else { // cutCellStatus[i] == 1
          
          if (interfaceExtrapolation == 1) {
	    for (int k = 0; k < dim; ++k) {
	      
	      Vir[k] = (interfaceLocation-X[i-2][0])/(X[i-1][0]-X[i-2][0])*V[i-1][k] - (interfaceLocation-X[i-1][0])/(X[i-1][0]-X[i-2][0])*V[i-2][k];
	      Vjr[k] = (interfaceLocation-X[j+1][0])/(X[j][0]-X[j+1][0])*V[j][k] - (interfaceLocation-X[j][0])/(X[j][0]-X[j+1][0])*V[j+1][k];
	      //std::cout << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
          } else {
	    memcpy(Vir, V[i-1], sizeof(double)*dim);
	    memcpy(Vjr, V[j], sizeof(double)*dim);
	    for (int k = 0; k < dim; ++k) {
	      //std::cout << i << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
          }
        }
 
        memset(fluxi,0,sizeof(double)*dim);
        memset(fluxj,0,sizeof(double)*dim);
	varFcn->getVarFcnBase(fluidId[i])->verification(0,Udummy,Vir);
	varFcn->getVarFcnBase(fluidId[j])->verification(0,Udummy,Vjr);

        riemann->computeRiemannSolution(Vir,Vjr,fluidId[i-1],fluidId[j+1],gradphi,varFcn,
	  			        Wir,Wjr,i,j,i,dx);

        if (interfaceExtrapolation == 1) {
          if (cutCellStatus[j] == 1) {
	    for (int k = 0; k < dim; ++k) {
	      Wi[k] = (Y[i+1][0]-X[i][0])/(interfaceLocation-X[i][0])*Wir[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[i][0])*V[i][k];
	    }
          } else {
	    for (int k = 0; k < dim; ++k) {
	      Wj[k] = (Y[i+1][0]-X[j][0])/(interfaceLocation-X[j][0])*Wjr[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[j][0])*V[j][k];
            }
          }
	} else {
          memcpy(Wi, Wir, sizeof(double)*dim);
	  memcpy(Wj, Wjr, sizeof(double)*dim);
        }

	if (interfaceTreatment == 0) {
	  memcpy(Wr[i], Wj, sizeof(double)*5);
	  memcpy(Wr[j], Wi, sizeof(double)*5);
	} else {
	  if (cutCellStatus[j] == 1)
	    memcpy(Wr[i], Wi, sizeof(double)*5);
	  else
	    memcpy(Wr[j], Wj, sizeof(double)*5);
	  
	}

	  
	riemannStatus[i] = riemannStatus[j] = 1;

        //fprintf(stderr,"Edge %d-%d crosses interface\n",i,j);

        if (cutCellStatus[j] == 1) {
          fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wi, Wi, fluxi, fluidId[i]);
        } else {

          fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wj, Wj, fluxj, fluidId[j]);
        }

        for (int k=0; k<dim; k++){
          R[i][k] += ctrlSurf[iEdge+1][0]*fluxi[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*fluxj[k];
        }
      }
    }
  }
	
  double dummy[dim];
  // flux at maxDistance
  fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[numPoints-1], BC[1], flux, fluidId[numPoints-1]);
  for (int k=0; k<dim; ++k)
    R[numPoints-1][k] += ctrlSurf[numPoints][0]*flux[k];

  //R[numPoints-1][1] -= varFcn->getPressure(V[numPoints-1],0)*ctrlSurf[numPoints][0];

  // flux at left (for cartesian) - use of non-reflecting BC
  // flux at center (for cylindrical and spherical) - use of wall=symmetry
  normal[0] = ctrlSurf[0][0];
  if(coordType == OneDimensionalInfo::CARTESIAN)
    fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[0], BC[0], flux, fluidId[0]);
  else if(coordType == OneDimensionalInfo::CYLINDRICAL)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);
  else if(coordType == OneDimensionalInfo::SPHERICAL)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);

  //std::cout << "fnew = " << flux[0] << " " << flux[1] << " " << flux[2] << " " << flux[3] << " " << flux[4] << std::endl;

  for (int k=0; k<dim; ++k)
    R[0][k] -= ctrlSurf[0][0]*flux[k];
  
  
  if (source) {
    EulerSource E(varFcn,fluidId);
    source->compute(E,y, R, X, Y,fluidId);
  }

  //std::cout << R*R << std::endl;

  // for debug
  //cout<<"flux[0] = "<<flux[0]<<" "<<flux[1]<<" "<<flux[2]<<" "<<flux[3]<<" "<<flux[4]<<endl;
  //cout<<"ctrlSurf[0] = "<<ctrlSurf[0][0]<<endl;
}
//------------------------------------------------------------------------------
void OneDimensional::PhiF(double t, SVec<double,1>& y,SVec<double,1>& k) {

  Rphi = 0.0;
  computeLevelSetFluxes(y);
  for(int i=0; i<numPoints; i++){
    k[i][0] = -Rphi[i][0] / ctrlVol[i][0];
  }
}

class LevelSetSource {

public:

  LevelSetSource(VarFcn* _vf,SVec<double,5>& _lu) :  vf(_vf), locU(_lu) { }

  ~LevelSetSource() { }

  void compute(double* U, double* f,int i) {

    f[0] = locU[i][1]*U[0]/locU[i][0];
  }

  VarFcn* vf;
  SVec<double,5>& locU;
};

double max(double e, double f) {
  return (e > f ? e : f);
}

double max(double d, double e, double f) {
  double q = max(e,f);
  return (d > q ? d : q);
}

double max(double c,double d, double e, double f) {
  double q = max(d,e,f);
  return (c > q ? c : q);
}

double max(double b,double c,double d, double e, double f) {
  double q = max(c,d,e,f);
  return (b > q ? b : q);
}

double max(double a,double b,double c,double d, double e, double f) {
  double q = max(b,c,d,e,f);
  return (a > q ? a : q);
}

void OneDimensional::computeLevelSetFluxes(SVec<double,1>& y){

  double uroe,flux;

  varFcn->conservativeToPrimitive(U,V,&fluidId);
  
  computeSlopes(y,Phislope,fluidId,false);

  if (levelSetMethod == 0) {
    double Phii[1],Phij[1],Phisi[1],Phisj[1];
    for(int iEdge=0; iEdge<numPoints-1; iEdge++){
      int i = iEdge;
      int j = iEdge+1;
      
      if (!isSixthOrder || !(i > 0)) {
	
	Phisi[0] = Phislope[i][0]*(X[j][0]-X[i][0]);
	Phisj[0] = Phislope[j][0]*(X[j][0]-X[i][0]);
      } else {
	
	Phisi[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0])+
	  (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Phislope[i-1][0]-2.0*Phislope[i][0]+Phislope[j][0])/beta;
	
	Phisj[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])+
	  (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Phislope[j+1][0]-2.0*Phislope[j][0]+Phislope[i][0])/beta;
	
	Phisi[0] *= 0.5*(X[j][0]-X[i][0]);
	Phisj[0] *= 0.5*(X[j][0]-X[i][0]);
	
      }
      
      recFcnLS->compute(y[i], Phisi,y[j], Phisj, Phii, Phij);
      
      double uav = 0.5*(V[i][1]+V[j][1]);
      
      if(uav > 0.0){
	flux = Phii[0]*uav;
      }
      else{
	flux = Phij[0]*uav;
      }
      Rphi[i][0] += flux*ctrlSurf[iEdge+1][0];
      Rphi[j][0] -= flux*ctrlSurf[iEdge+1][0];
    }
    
    // flux at maxDistance
    flux = (V[numPoints-1][1] > 0) ? Phi[numPoints-1][0]*V[numPoints-1][1] : V[numPoints-1][0]*BCphi[1]*V[numPoints-1][1];
    Rphi[numPoints-1][0] += flux*ctrlSurf[numPoints][0];
    
    // flux at center is zero since u=0 for radial and spherical
    // but for cartesian, flux needs to be computed
    if(coordType == OneDimensionalInfo::CARTESIAN){
      flux = (V[0][1] < 0) ? Phi[0][0]*V[0][1] : V[0][0]*BCphi[0]*V[0][1];
      Rphi[0][0] -= flux*ctrlSurf[0][0];
    }
    
    // source term
    
    if (source) {
      for (int i = 0; i < numPoints; ++i)
        Rphi[i][0] += 2.0*U[i][1]/U[i][0]*Phi[i][0]*ctrlVol[i][0]/(X[i][0] > 0 ? X[i][0] : 1e8);
        //LevelSetSource E(varFcn,U);
      
        //source->compute(E,y, Rphi, X, Y,fluidId);
    }
  } else if (levelSetMethod == 1) // H-J WENO 
    {
      for(int i=0; i<numPoints; i++){
	double u = V[i][1];
	
	double v1 = 0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;
	if (u > 0) {
	  if (i > 2)
	    v1 = (y[i-2][0]-y[i-3][0])/(X[i-2][0]-X[i-3][0]);
	  if (i > 1)
	    v2 = (y[i-1][0]-y[i-2][0])/(X[i-1][0]-X[i-2][0]);
	  if (i > 0)
	    v3 = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	  if (i < numPoints-1)
	    v4 = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	  if (i < numPoints-2)
	    v5 = (y[i+2][0]-y[i+1][0])/(X[i+2][0]-X[i+1][0]);
	} else {
	  if (i > 1)
	    v5 = (y[i-1][0]-y[i-2][0])/(X[i-1][0]-X[i-2][0]);
	  if (i > 0)
	    v4 = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	  if (i < numPoints-1)
	    v3 = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	  if (i < numPoints-2)
	    v2 = (y[i+2][0]-y[i+1][0])/(X[i+2][0]-X[i+1][0]);
	  if (i < numPoints-3)
	    v1 = (y[i+3][0]-y[i+2][0])/(X[i+3][0]-X[i+2][0]);
	}

	double S1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3)+0.25*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double S2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4)+0.25*(v2-v4)*(v2-v4);
	double S3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5)+0.25*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)+1.0e-30;
	double alpha1 = 0.1/((S1+eps)*(S1+eps));
	double alpha2 = 0.6/((S2+eps)*(S2+eps));
	double alpha3 = 0.3/((S3+eps)*(S3+eps));
	double asum = (alpha1+alpha2+alpha3);
	double om1 = alpha1/asum, om2 = alpha2/asum, om3 = alpha3/asum;
	
	double phi1 = v1/3.0-7.0/6.0*v2+11.0/6.0*v3;
	double phi2 = -v2/6.0+5.0/6.0*v3+1.0/3.0*v4;
	double phi3 = v3/3.0+5.0/6.0*v4-v5/6.0;
	double phix = om1*phi1+om2*phi2+om3*phi3;
	/*double phix = 0.0;
	if (u > 0) {
	  if (i > 1) {
	    double chi = (X[i][0]-X[i-1][0])/(X[i-1][0]-X[i-2][0]);
	    phix = (1.0+2.0*chi)/(1.0+chi)*y[i][0] + (-1.0-chi)*y[i-1][0]+(-(1.0+2.0*chi)/(1.0+chi)+1.0+chi)*y[i-2][0];
	    phix /= (X[i][0]-X[i-1][0]);
	  } else if (i > 0)
	    phix = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	} else if (u < 0) {
	  if (i < numPoints-2) {
	    double chi = (X[i][0]-X[i+1][0])/(X[i+1][0]-X[i+2][0]);
	    phix = (1.0+2.0*chi)/(1.0+chi)*y[i][0] + (-1.0-chi)*y[i+1][0]+(-(1.0+2.0*chi)/(1.0+chi)+1.0+chi)*y[i+2][0];
	    phix /= (X[i][0]-X[i+1][0]);
	  } else if (i < numPoints-1)
	    phix = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	    }*/
	Rphi[i][0] = u*phix * ctrlVol[i][0];
	 
      }   
      
      /*if (source) {
	LevelSetSource E(varFcn,U);
	
	source->compute(E,y, Rphi, X, Y,fluidId);
	}*/
    }  
}
//------------------------------------------------------------------------------
void OneDimensional::resultsOutput(double time, int iteration){

  fstream output;
  int i;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  for (i=0; i<PostFcn::SSIZE; ++i) {
    if (scalars[i]) {

      if (iteration == 0)
	output.open(scalars[i], fstream::out);
      else
	output.open(scalars[i], fstream::out | fstream::app);
      output << time*refVal.time  << endl;
      for (int j = 0; j < numPoints; ++j) {
	switch ( (PostFcn::ScalarType)i ) {
	case PostFcn::DENSITY:
	  output << V[j][0]*sscale[i]; break;
	case PostFcn::MACH:
	  output << fabs(V[j][1]) / varFcn->computeSoundSpeed(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::PRESSURE:
	  output << varFcn->getPressure(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::TEMPERATURE:
	  output << varFcn->computeTemperature(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::PHILEVEL:
	  output << Phi[j][0] / V[j][0]*sscale[i]; break;
	case PostFcn::VELOCITY_NORM:
	  output << fabs(V[j][1])*sscale[i]; break;
	case PostFcn::FLUIDID:
	  output << fluidId[j]*sscale[i]; break;
	default:
	  break;	
	}  
	output << endl;
      }
      output << endl;
      output.close();
    }
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    if (vectors[i]) {
      if (iteration == 0)
	output.open(vectors[i], fstream::out);	
      else
	output.open(vectors[i], fstream::out | fstream::app);

      output << time*refVal.time << endl;
      for (int j = 0; j < numPoints; ++j) {
	switch ( (PostFcn::VectorType)i ) {
	case PostFcn::VELOCITY:
	  output << V[j][1]*vscale[i] << " " << 0.0 << " " << 0.0; break;
	default:
	  break;
	} 
        output << endl; 
      }
      output << endl;
      output.close();
    }
  }

  if (bubbleRadiusFile[0] != 0) {
    
    if (iteration == 0)
      output.open(bubbleRadiusFile, fstream::out);
    else
      output.open(bubbleRadiusFile, fstream::out | fstream::app);
    double rad = 0.0;
    for (i=0; i<numPoints; ++i) {
      if (fluidId[i+1] == 0) {
	if (levelSetMethod == 0)
	  rad = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
	else if (levelSetMethod == 1)
	  rad = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
	else
	  rad = interfaceLocation;
	break;
      }
    }
    output << time*refVal.time << " " <<  rad*refVal.length << endl;
    output.close();
  }
}

void OneDimensional::restartOutput(double time, int iteration){

  //if(iteration==0) cout << "outputting results in file "<<outfile<<endl;
  fstream output;
  int sp = strlen(outfile)+1;
  //char str[10];
  //sprintf(str,"%d",iteration);
  //char *currentfile = new char[sp + strlen(str) ];
  output.open(outfile, fstream::out | fstream::trunc);
  output.setf(ios::scientific);
  output.precision(20);

  cout << "outputting last results in file "<<outfile<<endl;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  double rad = 0.0;
  for (int i=0; i<numPoints; ++i) {
    if (fluidId[i+1] == 0) {
      if (levelSetMethod == 0)
	rad = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
      else if (levelSetMethod == 1)
	rad = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
      else
	rad = interfaceLocation;
      break;
    }
  }
  
  output << "# time = " << time*refVal.time << endl;
  output << "# " << numPoints << endl;
  for(int i=0; i<numPoints; i++) {
    output << X[i][0]*refVal.length <<" "<< V[i][0]*refVal.density <<" "<<
      V[i][1]*refVal.velocity <<" "<<
      varFcn->getPressure(V[i],fluidId[i])*refVal.pressure <<" ";
    if (levelSetMethod == 0)
      output << Phi[i][0]/V[i][0];
    else if (levelSetMethod == 1)
      output << Phi[i][0];
    else if (levelSetMethod == 2)
      output << -(X[i][0]-rad)*refVal.length;
    output << " " << fluidId[i] << " " << 
      varFcn->computeTemperature(V[i],fluidId[i])*refVal.temperature << endl;
  }
  output.close();

}

//------------------------------------------------------------------------------

template <int neq>
void OneDimensional::computeSlopes(SVec<double,neq>& VV, SVec<double,neq>& slopes,
				   Vec<int>& fid,bool crossInterface) {
  
  int i,j;
  double r,sig;
  int stat = 0;
  slopes = 0.0;
  for(i=0; i<numPoints; ++i) {

    //double A[4] = { X[i+1][0]*X[i+1][0]+X[i-1][0]*X[i-1][0]+X[i][0]*X[i][0], X[i+1][0]+X[i-1][0]+X[i][0],
    //		   X[i+1][0]+X[i-1][0]+X[i][0], 3};
  //double det = A[0]*A[3]-A[1]*A[2];

    if (loctag[i])
      continue;
    
    //if (crossInterface) {
      if (i > 0 && i < numPoints-1 &&
	  (interfaceTreatment == 0 && fid[i] == fid[i+1] && fid[i] == fid[i-1] ||
	   interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i-1] && cutCellStatus[i] == cutCellStatus[i+1]
	   || !crossInterface))
	stat = 0;
      else if (i < numPoints-1 && 
	       (interfaceTreatment == 0 && fid[i] == fid[i+1] ||
		interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i+1] ||
		!crossInterface))
	stat = 1;
      else if (i > 0 && 
	       (interfaceTreatment == 0 && fid[i] == fid[i-1] ||
		interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i-1] ||
		!crossInterface) )
	stat = 2;
      else
	stat = 3;
      //}
    //std::cout << stat << std::endl;
    for (j = 0; j < neq; ++j) {

      if (stat == 0) {   
	double dx1 = X[i+1][0]-X[i][0], dx2 = X[i-1][0]-X[i][0];
	slopes[i][j] = ((V[i+1][j]-V[i][j])*dx1+(VV[i-1][j]-VV[i][j])*dx2)/(dx1*dx1+dx2*dx2); 
	//slopes[i][j] = (VV[i+1][j]-VV[i-1][j])/(X[i+1][0]-X[i-1][0]);
      }
      else if (stat == 1) {
	//double dx1 = X[i+1][0]-X[i][0], dx2 = X[i-1][0]-X[i][0];
	slopes[i][j] = (VV[i+1][j]-VV[i][j])/(X[i+1][0]-X[i][0]);
      }
      else if (stat == 2)
	slopes[i][j] = (VV[i][j]-VV[i-1][j])/(X[i][0]-X[i-1][0]);
      else
	slopes[i][j] = 0.0;

      //slopes[i][j] *= 0.5;
    }
  }

}

RecFcn* OneDimensional::createRecFcn(IoData &ioData)
{

  RecFcn *rf = 0;

  double beta = ioData.schemes.ns.beta;
  double eps = ioData.schemes.ns.eps;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    /*if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
      if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) {
	if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	  rf = new RecFcnConstant<6>;
      }
      else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
	if (ioData.schemes.ns.limiter == SchemeData::NONE) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLinearConstant<6>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnLinearVanAlbada<6>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnVanAlbadaConstant<6>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnVanAlbada<6>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLtdLinearConstant<6>(beta, eps);
	}
      }
    }
    else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
      if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) {
	if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	  rf = new RecFcnConstant<7>;
      }
      else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
	if (ioData.schemes.ns.limiter == SchemeData::NONE) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLinearConstant<7>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnLinearVanAlbada<7>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnVanAlbadaConstant<7>(beta, eps);
	  else if (ioData.schemes.tm.reconstruction == SchemeData::LINEAR) {
	    if (ioData.schemes.tm.limiter == SchemeData::VANALBADA)
	      rf = new RecFcnVanAlbada<7>(beta, eps);
	  }
	}
	else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR) {
	  if (ioData.schemes.tm.reconstruction == SchemeData::CONSTANT)
	    rf = new RecFcnLtdLinearConstant<7>(beta, eps);
	}
      }
      }*/
  } else {
    if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT)
      rf = new RecFcnConstant<dim>;
    else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
      if (ioData.schemes.ns.limiter == SchemeData::NONE)
	rf = new RecFcnLinear<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA)
	rf = new RecFcnVanAlbada<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::BARTH)
	rf = new RecFcnBarth<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VENKAT)
	rf = new RecFcnVenkat<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR)
	rf = new RecFcnLtdLinear<dim>(beta, eps);
    }
  }

  if (!rf) {
    fprintf(stderr, "*** Error: no valid choice for the reconstruction\n");
    exit(1);
  }

  return rf;

}

RecFcn* OneDimensional::createRecFcnLS(IoData &ioData)
{
  RecFcn *rf = 0;

  double beta = ioData.schemes.ls.beta;
  double eps = ioData.schemes.ls.eps;

  if (ioData.schemes.ls.reconstruction == SchemeData::CONSTANT)
    rf = new RecFcnConstant<1>;
  else if (ioData.schemes.ls.reconstruction == SchemeData::LINEAR) {
    if (ioData.schemes.ls.limiter == SchemeData::NONE)
      rf = new RecFcnLinear<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VANALBADA)
      rf = new RecFcnVanAlbada<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::BARTH)
      rf = new RecFcnBarth<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VENKAT)
      rf = new RecFcnVenkat<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::P_SENSOR)
      rf = new RecFcnLtdLinear<1>(beta, eps);
  }

  return rf;
}

void OneDimensional::loadSparseGrid(IoData& ioData) {

  if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
    // only the ioData.eqs.fluidModel is considered since only the Riemann invariant of one EOS is tabulated!
    // (no need to specify two different EOS)

    double *refIn  = new double[2];
    double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = pow(ioData.ref.rv.density,-ioData.eqs.fluidModel.jwlModel.omega)*ioData.ref.rv.velocity*ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.velocity;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);

  }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){

    double *refIn = new double[5]; double *refOut = new double[2];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = ioData.ref.rv.pressure;
    refIn[2] = ioData.ref.rv.density;
    refIn[3] = ioData.ref.rv.pressure;
    refIn[4] = ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.density;
    refOut[1] = ioData.ref.rv.density;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);
  }else{
    tabulationC = 0;
  }
}
