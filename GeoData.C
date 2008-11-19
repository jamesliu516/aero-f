#include <GeoData.h>

//------------------------------------------------------------------------------

GeoData::GeoData(IoData &ioData)
{

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_STEADY_SENSITIVITY_ANALYSIS_) {

    if (ioData.ts.type != TsData::IMPLICIT)
      ioData.ts.type = TsData::IMPLICIT;

    if (ioData.ts.implicit.type != ImplicitData::BACKWARD_EULER)
      ioData.ts.implicit.type = ImplicitData::BACKWARD_EULER;

    if (ioData.dgcl.normals != DGCLData::IMPLICIT_FIRST_ORDER_GCL)
      ioData.dgcl.normals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;
      
    if (ioData.dgcl.velocities != DGCLData::IMPLICIT_ZERO)
      typeVelocities = DGCLData::IMPLICIT_ZERO;

  }

  use_n = false;
  use_nm1 = false;
  use_nm2 = false;
  use_save = false;
  typeNormals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;
  typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;
  typeVolumeChanges = DGCLData::AUTO_VOL;

  if (ioData.problem.type[ProblemData::ACCELERATED] ||
      ioData.problem.type[ProblemData::AERO] ||
      ioData.problem.type[ProblemData::FORCED] ||
      ioData.problem.type[ProblemData::ROLL]) {
    if (ioData.ts.type == TsData::IMPLICIT) {
      use_n = true;
      if (ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
	use_nm1 = true;
      if (ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF) {
	use_nm1 = true;
	use_nm2 = true;
      }

      if (ioData.dgcl.normals == DGCLData::AUTO) {
	if (ioData.ts.implicit.type == ImplicitData::BACKWARD_EULER ||
	    ioData.ts.implicit.type == ImplicitData::CRANK_NICOLSON)
	  typeNormals = DGCLData::IMPLICIT_FIRST_ORDER_GCL;
	else if (ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
	  typeNormals = DGCLData::IMPLICIT_SECOND_ORDER_GCL;
	else if (ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF)
	  typeNormals = DGCLData::IMPLICIT_THIRD_ORDER_EZGCL;
      }
      else
	typeNormals = ioData.dgcl.normals;    
      
      if (ioData.dgcl.velocities == DGCLData::AUTO_VEL)
	typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;
      else
	typeVelocities = ioData.dgcl.velocities;
    }
    else if (ioData.ts.type == TsData::EXPLICIT) {
      use_n = true;
      use_save = true;
      if (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2){
        typeNormals = DGCLData::EXPLICIT_RK2;
        typeVelocities = DGCLData::EXPLICIT_RK2_VEL;
        typeVolumeChanges = DGCLData::EXPLICIT_RK2_VOL;
      }
      else{ // RK4 does not have proper algorithm yet!
        if (ioData.dgcl.normals == DGCLData::AUTO)
	  typeNormals = DGCLData::IMPLICIT_LATEST_CFG;
        else
	  typeNormals = ioData.dgcl.normals;  // CBM  
  
        if (ioData.dgcl.velocities == DGCLData::AUTO_VEL)
          typeVelocities = DGCLData::IMPLICIT_BACKWARD_EULER_VEL;
        else
          typeVelocities = ioData.dgcl.velocities;  // CBM

        typeVolumeChanges = DGCLData::EXPLICIT_RK2_VOL;
      }
    }
  }
  if (ioData.problem.type[ProblemData::LINEARIZED])  {
    use_n = true;
    use_nm1 = true;
  }
}

//------------------------------------------------------------------------------
