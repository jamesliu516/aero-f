cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(MPI)

FIND_PATH(PHYSBAM_BASE PHYSBAM_INTERFACE.h PATHS ../../PhysBAM)
SET(INC_DIR .. . ../LevelSet ${PHYSBAM_BASE} ${PHYSBAM_BASE}/Public_Library)
FIND_LIBRARY(PhysBAM_Base_Lib NAMES physbam PATHS ${PHYSBAM_BASE})

IF(MPI_FOUND)
  MESSAGE("Will compile with MPI at: " ${MPI_INCLUDE_PATH})
  SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
  add_definitions(-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK)
ENDIF(MPI_FOUND)

INCLUDE_DIRECTORIES(${INC_DIR})

#same definitions as AERO-F main code
add_definitions(-DTEMPLATE_FIX -D__LINUX -DF_NEEDS_UNDSC -DTYPE_PREC=float)

add_library(IntersectorFRG SHARED PhysBAMIntersect.C)

TARGET_LINK_LIBRARIES(IntersectorFRG ${PhysBAM_Base_Lib})
