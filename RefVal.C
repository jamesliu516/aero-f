#include <RefVal.h>
#include <IoData.h>
#include <math.h>

//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	new way to define reference values for viscosity, velocity and 
//		temperature 
//------------------------------------------------------------------------------

RefVal::RefVal()
{

  mode = NON_DIMENSIONAL;
  length = 1.0;
  density = 1.0;
  velocity = 1.0;
  pressure = 1.0;
  temperature = 1.0;
  viscosity_mu = 1.0;
  viscosity_lambda = 1.0;
  nutilde = 1.0;
  kenergy = 1.0;
  epsilon = 1.0;
  time = 1.0;
  force = 1.0;
  energy = 1.0;
  power = 1.0;

  tlength = 1.0;
  tvelocity = 1.0;
  tforce = 1.0;
  tpower = 1.0;

}

//------------------------------------------------------------------------------
RefVal::RefVal(IoData &ioData)
{


  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {
    mode = NON_DIMENSIONAL;
    length = 1.0;
    density = 1.0;
    velocity = 1.0;
    pressure = 1.0;
    temperature = 1.0;
    viscosity_mu = 1.0;
    viscosity_lambda = 1.0;
    nutilde = 1.0;
    kenergy = 1.0;
    epsilon = 1.0;
    time = 1.0;
    force = 1.0;
    energy = 1.0;
    power = 1.0;

    tlength = 1.0;
    tvelocity = 1.0;
    tforce = 1.0;
    tpower = 1.0;
  }
  else if (ioData.problem.mode  == ProblemData::DIMENSIONAL) {
    mode = DIMENSIONAL;
    if(ioData.eqs.fluidModel.fluid == FluidModelData::GAS){
      double gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
      double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
      double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      velocity = mach * sqrt(gam * (ioData.ref.pressure + Pstiff) / density);
      pressure = density * velocity*velocity;
      temperature = gam*(gam - 1.0) * mach*mach * (ioData.ref.pressure + Pstiff)/(density*R);
    }
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID){
      double Cv = ioData.eqs.fluidModel.liquidModel.Cv;
      double awater = ioData.eqs.fluidModel.liquidModel.alpha;
      double bwater = ioData.eqs.fluidModel.liquidModel.beta;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      velocity = mach * sqrt(awater*bwater*pow(density, bwater - 1.0));
      pressure = density * velocity*velocity;
      temperature = velocity*velocity/Cv;
    }

    length = ioData.ref.length;
    //surface = ioData.ref.surface;

    force = density * velocity*velocity * length*length;
    moment = force * length;
    energy = force * length;

    //cf2force = force * surface / (2.0 * length*length);
    //cm2moment = moment * surface / (2.0 * length*length);
    cf2force = force  / (2.0 * length*length);
    cm2moment = moment / (2.0 * length*length);

    time = length / velocity;
  }

}

//------------------------------------------------------------------------------
