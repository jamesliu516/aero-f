#include <GeoData.h>

//------------------------------------------------------------------------------

GeoData::GeoData(IoData &ioData)
{

  use_n = false;
  use_nm1 = false;
  use_nm2 = false;
  typeNormals = ImplicitData::FIRST_ORDER_GCL;
  typeVelocities = ImplicitData::BACKWARD_EULER_VEL;

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

      if (ioData.ts.implicit.normals == ImplicitData::AUTO) {
	if (ioData.ts.implicit.type == ImplicitData::BACKWARD_EULER ||
	    ioData.ts.implicit.type == ImplicitData::CRANK_NICOLSON)
	  typeNormals = ImplicitData::FIRST_ORDER_GCL;
	else if (ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
	  typeNormals = ImplicitData::SECOND_ORDER_GCL;
	else if (ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF)
	  typeNormals = ImplicitData::THIRD_ORDER_EZGCL;
      }
      else
	typeNormals = ioData.ts.implicit.normals;    
      
      if (ioData.ts.implicit.velocities == ImplicitData::AUTO_VEL)
	typeVelocities = ImplicitData::BACKWARD_EULER_VEL;
      else
	typeVelocities = ioData.ts.implicit.velocities;
    }
    else if (ioData.ts.type == TsData::EXPLICIT) {
      use_n = true;
      if (ioData.ts.implicit.normals == ImplicitData::AUTO)
	typeNormals = ImplicitData::LATEST_CFG;
      else
	typeNormals = ioData.ts.implicit.normals;  // CBM  

      if (ioData.ts.implicit.velocities == ImplicitData::AUTO_VEL)
        typeVelocities = ImplicitData::BACKWARD_EULER_VEL;
      else
        typeVelocities = ioData.ts.implicit.velocities;  // CBM
    }
  }
  if (ioData.problem.type[ProblemData::LINEARIZED])  {
    use_n = true;
    use_nm1 = true;
  }
}

//------------------------------------------------------------------------------
