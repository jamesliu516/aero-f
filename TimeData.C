#include <TimeData.h>
#include <DistVector.h>

#include <cmath>

//------------------------------------------------------------------------------

TimeData::TimeData(IoData &ioData)
{ 

  if (ioData.ts.typeTimeStep == TsData::AUTO) {
    if (ioData.problem.type[ProblemData::LINEARIZED])
      typeTimeStep = TsData::GLOBAL;
    else if (ioData.problem.type[ProblemData::UNSTEADY])
      typeTimeStep = TsData::GLOBAL;
    else
      typeTimeStep = TsData::LOCAL;
  }
  else
    typeTimeStep = ioData.ts.typeTimeStep;

  if (ioData.ts.type == TsData::IMPLICIT)
    typeIntegrator = ioData.ts.implicit.type;
  else if (ioData.ts.type == TsData::EXPLICIT)
    typeIntegrator = ImplicitData::BACKWARD_EULER;
  typeStartup = ioData.ts.implicit.startup;

  dt_imposed = ioData.ts.timestep;
  dt_n = ioData.restart.dt_nm1;
  dt_nm1 = ioData.restart.dt_nm1;
  dt_nm2 = ioData.restart.dt_nm2;
  output_newton_step = ioData.restart.output_newton_step;

  errorTol = ioData.ts.errorTol;

  exist_nm1 = false;
  exist_nm2 = false;
  use_nm1 = false;
  if (typeIntegrator == ImplicitData::THREE_POINT_BDF || typeIntegrator == ImplicitData::FOUR_POINT_BDF) 
    use_nm1 = true;
  use_nm2 = false;
  if (typeIntegrator == ImplicitData::FOUR_POINT_BDF) 
    use_nm2 = true;

  if (ioData.problem.type[ProblemData::LINEARIZED])
    use_modal = true;
  else
    use_modal = false;


  if (ioData.ts.form == TsData::DESCRIPTOR) {
    descriptor_form = 1;
  } else if (ioData.ts.form == TsData::HYBRID) {
    descriptor_form = 2;
  } else {
    descriptor_form = 0;
  }


// Included (MB)
  if (ioData.sa.comp3d == SensitivityAnalysis::OFF_COMPATIBLE3D)
    use_modal = true;

  if (ioData.linearizedData.domain == LinearizedData::FREQUENCY) use_freq = true;
  else use_freq = false;

}

//------------------------------------------------------------------------------

void TimeData::update()
{

  dt_nm2 = dt_nm1;
  dt_nm1 = dt_n;

}

//------------------------------------------------------------------------------

void TimeData::computeCoefficients(DistVec<double> &dt, double dt_glob)
{

  dt_n = dt_glob;
  if (typeTimeStep == TsData::GLOBAL)
    dt = dt_n;

  tau_n = dt_n / dt_nm1;
  tau_nm1 = dt_nm1 / dt_nm2;

  if (typeIntegrator == ImplicitData::BACKWARD_EULER || 
      typeIntegrator == ImplicitData::CRANK_NICOLSON || 
      (typeIntegrator == ImplicitData::THREE_POINT_BDF && !exist_nm1) ||
      (typeIntegrator == ImplicitData::FOUR_POINT_BDF && !exist_nm1 && !exist_nm2)) {
    alpha_np1 = 1.0;
    alpha_n = -1.0;
    alpha_nm1 = 0.0;
    alpha_nm2 = 0.0;
  } 
  else if ((typeIntegrator == ImplicitData::THREE_POINT_BDF && exist_nm1) ||
	   (typeIntegrator == ImplicitData::FOUR_POINT_BDF && exist_nm1 && !exist_nm2)) {
    alpha_np1 = (1.0 + 2.0 * tau_n) / (1.0 + tau_n);
    alpha_n = -1.0 - tau_n;
    alpha_nm1 = - alpha_np1 - alpha_n;
    alpha_nm2 = 0.0;
  }
  else if (typeIntegrator == ImplicitData::FOUR_POINT_BDF && exist_nm1 && exist_nm2) {
    alpha_np1 = (3.0*tau_n*tau_n*tau_nm1 + 4.0*tau_n*tau_nm1 + 2.0*tau_n + tau_nm1 + 1.0) /
      ((tau_n + 1.0)*(tau_n*tau_nm1 + tau_nm1 + 1.0));
    alpha_n = - (tau_n + 1.0) * (tau_n*tau_nm1 + tau_nm1 + 1.0) / (tau_nm1 + 1.0);
    alpha_nm1 = tau_n*tau_n * (tau_n*tau_nm1 + tau_nm1 + 1.0) / (tau_n + 1.0);
    alpha_nm2 = - alpha_np1 - alpha_n - alpha_nm1;
  }

}

//------------------------------------------------------------------------------

void TimeData::computeVelocities(DGCLData::Velocities typeVelocities, 
				 DistSVec<double,3> &Xnm1, DistSVec<double,3> &Xn,
				 DistSVec<double,3> &X, DistSVec<double,3> &Xdot)
{

  if (typeVelocities == DGCLData::IMPLICIT_BACKWARD_EULER_VEL ||
      typeVelocities == DGCLData::IMPLICIT_IMPOSED_BACKWARD_EULER_VEL ||
      typeVelocities == DGCLData::EXPLICIT_RK2_VEL)  {
    Xdot = 1.0 / dt_n * (X - Xn);
  }
  else if (typeVelocities == DGCLData::IMPLICIT_THREE_POINT_BDF_VEL ||
	   typeVelocities == DGCLData::IMPLICIT_IMPOSED_THREE_POINT_BDF_VEL) {
    double c_np1 = alpha_np1 / dt_n;
    double c_n = alpha_n / dt_n;
    double c_nm1 = alpha_nm1 / dt_n;

    Xdot = c_np1 * X + c_n * Xn + c_nm1 * Xnm1;
  }
  else if (typeVelocities == DGCLData::IMPLICIT_ZERO)
    Xdot = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void TimeData::rstVar(IoData &ioData)
{

  dt_imposed = ioData.ts.timestep;
  dt_n = ioData.restart.dt_nm1;
  dt_nm1 = ioData.restart.dt_nm1;
  dt_nm2 = ioData.restart.dt_nm2;

}

//------------------------------------------------------------------------------
