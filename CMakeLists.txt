cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
enable_language(Fortran)

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
  if(OpenMP_FLAG_DETECTED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_CXX_FLAGS}")
    #add_definitions(-D_OPENMP)
  endif(OpenMP_FLAG_DETECTED)
ENDIF(USE_OPENMP)

# Code to determine if we should add some link libraries
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES GNU)
  EXEC_PROGRAM(${CMAKE_Fortran_COMPILER}
      ARGS "-### CMake/Modules/CMakeTestFortran.f"
      OUTPUT_VARIABLE CMAKE_COMPILER_OUTPUT RETURN_VALUE CMAKE_COMPILER_RETURN)
  IF(NOT CMAKE_COMPILER_RETURN)
#      MESSAGE("Will add -lgfortran")
      IF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*-lgfortran.*" )
         SET(FORTRAN_HYBRID_LIBRARIES "-lgfortran")
      ENDIF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*-lgfortran.*" )
#  ELSE(NOT CMAKE_COMPILER_RETURN)
#     MESSAGE("Problem: " ${CMAKE_COMPILER_OUTPUT} )
  ENDIF(NOT CMAKE_COMPILER_RETURN)
ENDIF("${CMAKE_Fortran_COMPILER_ID}" MATCHES GNU)

# SET(CMAKE_VERBOSE_MAKEFILE true)

IF(NOT EXTRALIB)
SET(EXTRALIB "")
ENDIF(NOT EXTRALIB)

SET(EXTRALIB ${EXTRALIB} CACHE STRING "Extra link parameters")

SET(CMAKE_DEBUG_POSTFIX "_debug")

FIND_PACKAGE(LAPACK)
IF(BLAS_LIBRARIES)
  SET(EXTRALIBRARIES ${BLAS_LIBRARIES})
ELSE(BLAS_LIBRARIES)
  SET(EXTRALIBRARIES "")
ENDIF(BLAS_LIBRARIES)
IF(LAPACK_LIBRARIES)
  SET(EXTRALIBRARIES ${LAPACK_LIBRARIES} ${EXTRALIBRARIES})
ENDIF(LAPACK_LIBRARIES)

FIND_PACKAGE(MPI)

SET(TARGETS_INST fluid)

IF(UNIX)
  MESSAGE("-- Building for system type: " ${CMAKE_SYSTEM_NAME} ".")
  IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)

    FIND_LIBRARY(ARPACKLIB NAMES arpack 
       PATHS  . /usr/lib /usr/local/lib)

    SET(DO_MODAL 0)
    IF(ARPACKLIB AND LAPACK_FOUND)
      MESSAGE("-- ARPACK library " ${ARPACKLIB} " found.")
      add_definitions(-DDO_MODAL)
      SET(EXTRALIBRARIES 
          ${ARPACKLIB}
          ${EXTRALIBRARIES})
      SET(DO_MODAL 1)
    ENDIF(ARPACKLIB AND LAPACK_FOUND)

    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
    FIND_LIBRARY(SCALAPACKLIB NAMES scalapack scalapack-openmpi
                 PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})
    FIND_LIBRARY(BLACSLIB NAMES blacs blacs-openmpi
                 PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})
    FIND_LIBRARY(BLACSCLIB NAMES blacsC blacsCinit-openmpi
                 PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SCALAPACK ${SCALAPACKDIR})

    #SET(USE_SCALAPACK 0)
    IF(SCALAPACKLIB AND BLACSLIB AND BLACSCLIB)
      MESSAGE("-- SCALAPACK library " ${SCALAPACKLIB} " found.")
      MESSAGE("-- BLACS library " ${BLACSLIB},${BLACSCLIB} " found.")    
      add_definitions(-DDO_SCALAPACK)
      SET(EXTRALIBRARIES
          ${SCALAPACKLIB} ${BLACSCLIB} ${BLACSLIB} 
          ${EXTRALIBRARIES})
      SET(USE_SCALAPACK 1)
    ELSE(USE_SCALAPACK AND SCALAPACKLIB AND BLACSLIB)
      SET(USE_SCALAPACK 0)
    ENDIF(SCALAPACKLIB AND BLACSLIB AND BLACSCLIB)
    
  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
ENDIF(UNIX)

SET(BINDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)

FIND_PROGRAM(MERCURIAL_EXECUTABLE hg)

SET(MAIN_FILES Main.C)

ADD_LIBRARY(FluidBase
             FSI/DynamicNodalTransfer.cpp
             LevelSet/FluidTypeCriterion.cpp
             NavierStokesSolver.C 
             Solvers/LevelSetSolver5.C
             Solvers/NavierStokesCoupledSolver5.C
             Solvers/NavierStokesEmbedded5.C
             Solvers/NavierStokesEmbedded6.C
             Solvers/NavierStokesEmbedded7.C
             Solvers/NavierStokesMultiPhysicsEmbedded5.C
             Solvers/NavierStokesCoupledSolver6.C
             Solvers/NavierStokesCoupledSolver7.C
             Solvers/NavierStokesSegSolver651.C
             Solvers/NavierStokesSegSolver752.C
             BcFcnCore.C
             CommunicatorCore.C 
             CommunicatorMPI.C 
             ConnectivityCore.C 
             CorotSolver.C 
             CurvatureDetection.C 
             DistGeoState.C 
             DistMacroCellCore.C 
             DomainCore.C 
             EdgeCore.C 
             FaceCore.C 
             FaceTriaCore.C 
             ElemCore.C 
             ElemTetCore.C 
             FemEquationTermDesc.C 
             FluxFcnDesc.C 
             FluxFcnDescSG.C 
             FluxFcnDescTait.C 
             FluxFcnDescJwl.C
             GeoData.C 
             GeoSource.C 
             GeoState.C 
             HeatTransferHandlerCore.C 
             ImplicitRiemann.C
             InletNodeCore.C 
             IoDataCore.C 
             KspConvCriterion.C 
             MacroCellCore.C 
             MatchNodeCore.C 
             MemoryPool.C 
             MeshMotionHandlerCore.C 
             MeshMotionSolver.C 
             ModalSolver.C 
             Node.C 
             PostFcn.C 
             ProgrammedBurnCore.C
             RefVal.C 
             SmagorinskyLESTerm.C 
             WaleLESTerm.C 
             DynamicLESTerm.C 
             StructExc.C 
             SubDomainCore.C 
             SparseGridCore.C
             SparseGridClusterCore.C
             TimeData.C 
             Timer.C 
             TsInput.C 
             TsParameters.C 
             TsRestartCore.C 
             VMSLESTerm.C 
             DynamicVMSTerm.C 
             WallFcnCore.C 
             BCApplierCore.C 
	     BCond.C 
             BlockAlloc.C 
             TriangulatedSurface.C
             FluidSelectorCore.C
             OneDimensionalSolver.C
            )

ADD_LIBRARY(Utils utils/Predicate.C utils/Laguerre.C)

SET_TARGET_PROPERTIES(Utils
                PROPERTIES COMPILE_FLAGS -fPIC)

IF(MPI_FOUND)
  ADD_EXECUTABLE(loader tools/loader.C)
  TARGET_LINK_LIBRARIES(loader ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
  SET_TARGET_PROPERTIES(loader
                      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
ENDIF(MPI_FOUND)

# Parser
# search bison
MACRO(FIND_BISON)
    IF(NOT BISON_EXECUTABLE)
        FIND_PROGRAM(BISON_EXECUTABLE bison)
        IF (NOT BISON_EXECUTABLE)
          MESSAGE(FATAL_ERROR "bison not found - aborting")
        ENDIF (NOT BISON_EXECUTABLE)
    ENDIF(NOT BISON_EXECUTABLE)
ENDMACRO(FIND_BISON)

MACRO(ADD_BISON_FILES _sources )
    FIND_BISON()

    FOREACH (_current_FILE ${ARGN})
      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)

      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.tab.C)

      ADD_CUSTOM_COMMAND(
         OUTPUT ${_out}
         COMMAND ${BISON_EXECUTABLE}
         ARGS
         -d
         -v
         -p yyCmdf
         -o${_out}
         ${_in}
         DEPENDS ${_in}
      )

      SET(${_sources} ${${_sources}} ${_out} )
   ENDFOREACH (_current_FILE)
ENDMACRO(ADD_BISON_FILES)


# flex a .ll file

# search flex
MACRO(FIND_FLEX)
    IF(NOT FLEX_EXECUTABLE)
        FIND_PROGRAM(FLEX_EXECUTABLE flex)
        IF (NOT FLEX_EXECUTABLE)
          MESSAGE(FATAL_ERROR "flex not found - aborting")
        ENDIF (NOT FLEX_EXECUTABLE)
    ENDIF(NOT FLEX_EXECUTABLE)
ENDMACRO(FIND_FLEX)

MACRO(ADD_FLEX_FILES _sources )
    FIND_FLEX()

    FOREACH (_current_FILE ${ARGN})
      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)

      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/flex_${_basename}.cpp)

      ADD_CUSTOM_COMMAND(
         OUTPUT ${_out}
         COMMAND ${FLEX_EXECUTABLE}
         ARGS
         -PyyCmdf 
         -o${_out}
         ${_in}
         DEPENDS ${_in}
      )

      SET(${_sources} ${${_sources}} ${_out} )
   ENDFOREACH (_current_FILE)
ENDMACRO(ADD_FLEX_FILES)

ADD_FLEX_FILES(MyLib parser/Lexer.l)
ADD_BISON_FILES(MyLib parser/Parser.y)

ADD_LIBRARY(Parser parser/Dictionary.C parser/StringManager.C parser/AssignerCore.C parser/ParseTree.cpp ${MyLib})
SET_TARGET_PROPERTIES(Parser
                PROPERTIES COMPILE_FLAGS -fPIC)

install(TARGETS FluidBase DESTINATION lib EXPORT fluid-targets)
install(EXPORT fluid-targets DESTINATION lib/base)
export(TARGETS FluidBase FILE ${CMAKE_CURRENT_BINARY_DIR}/base-exports.cmake)

	
IF(MERCURIAL_EXECUTABLE)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"" `${MERCURIAL_EXECUTABLE} identify -in` "\\\"" "\;" > ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C
                 DEPENDS Main.C FluidBase
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ELSE (MERCURIAL_EXECUTABLE)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"\\\"\;" > ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C
                 DEPENDS Main.C FluidBase
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF (MERCURIAL_EXECUTABLE)

###################
# Local Variables #
###################
SITE_NAME(HOSTNAME)
FIND_PACKAGE(LOCAL NAMES ${HOSTNAME} PATHS . ./localconfig QUIET)
MARK_AS_ADVANCED(HOSTNAME)
MARK_AS_ADVANCED(LOCAL_DIR)
###################

IF(EXTRALIB)
  SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB})
ENDIF(EXTRALIB)

SET(INC_DIR .)

##############################
# For the embedded framework #
##############################
SET (DO_EMBEDDED 1)
if (DO_EMBEDDED)

  FIND_PACKAGE(Boost 1.19.0 REQUIRED)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  endif(Boost_FOUND)

  SET(INC_DIR PhysBAM_Interface PhysBAM_Interface/Public_Library)
  INCLUDE_DIRECTORIES(${INC_DIR})

  ADD_LIBRARY(physbam SHARED
  PhysBAM_Interface/PHYSBAM_INTERFACE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Arrays/ARRAY_COLLECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Vectors/VECTOR.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/EXCEPTIONS.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/TIMER.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Utilities/PROCESS_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Data_Structures/KD_TREE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Parsing/STRING_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/DEBUG_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/LOG.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/LOG_ENTRY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Log/DEBUG_PRINT.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_BASE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X2.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_3X3.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_MXN.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/ROTATION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Matrices/MATRIX_4X4.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Polynomials/CUBIC.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Random_Numbers/MT19937.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/POLYGON.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/RAY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/PLANE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/BOX.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/POINT_FACE_COLLISION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Continuous_Collision_Detection/EDGE_EDGE_COLLISION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_LINE_2D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_PLANE_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/BOX_POINT_SIMPLEX_1D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_PLANE_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_BOX_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_ORIENTED_BOX_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_SEGMENT_2D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/SIMPLEX_MESH.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/SEGMENT_MESH.cpp
  PhysBAM_Interface/Public_Library/PhysBAM_Geometry/Topology/TRIANGLE_MESH.cpp
  )
  
  SET_TARGET_PROPERTIES(physbam
                PROPERTIES COMPILE_FLAGS "-fPIC -DCOMPILE_WITHOUT_READ_WRITE_SUPPORT -IPublic_Library")

  SET(INC_DIR_TMP . IntersectorFRG IntersectorPhysBAM LevelSet ${PHYSBAM_BASE} ${PHYSBAM_BASE}/Public_Library)
  INCLUDE_DIRECTORIES(${INC_DIR_TMP})

  #same definitions as AERO-F main code
  add_definitions(-fPIC -DTEMPLATE_FIX -D__LINUX -DF_NEEDS_UNDSC -DTYPE_PREC=float)
  add_library(Intersector SHARED
              IntersectorFRG/IntersectorFRG.C
              IntersectorPhysBAM/IntersectorPhysBAM.C
              IntersectorPhysBAM/FloodFill.C
              IntersectorPhysBAM/Mpi_Utilities.C)
  TARGET_LINK_LIBRARIES(Intersector physbam)

  ADD_DEFINITIONS(-DDO_EMBEDDED)
  SET (EXTRALIBRARIES ${EXTRALIBRARIES} Intersector)
ENDIF (DO_EMBEDDED)
##############################

ADD_EXECUTABLE(fluid ${MAIN_FILES} ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C)

SET_TARGET_PROPERTIES(FluidBase
                PROPERTIES COMPILE_FLAGS -fPIC)

ADD_LIBRARY(FortranLib 
f77src/ballvertex.f
f77src/boundflux5.f
f77src/boundjac2.f
f77src/dsvdc.f
f77src/eriemann_gg.f
f77src/eriemann_gw.f
f77src/eriemann_ww.f
f77src/g_boundflux5.f
f77src/genbcfluxgas.f
f77src/genbcfluxtait.f
f77src/genbcfluxjwl.f
f77src/g_roeflux5.f
f77src/g_roeflux6.f
f77src/roeflux5jwl.f
f77src/hllcflux.f
f77src/hlleflux.f
f77src/hlleflux1.f
f77src/hllejac.f
f77src/printmask.f
f77src/roeflux5.f
f77src/roeflux5prim.f
f77src/roeflux5waterdissprec.f
f77src/roeflux6.f
f77src/roejac2.f
f77src/roejac2water.f
f77src/roejac5.f
f77src/roejac5waterdissprec.f
f77src/roejac6.f
f77src/roejac6jwl.f
f77src/roejacwaterdissprim.f
f77src/torsionspring.f)

SET_TARGET_PROPERTIES(FortranLib
                      PROPERTIES
                      COMPILE_FLAGS -fPIC
                      RELEASE_OUTPUT_NAME Fortran
                      DEBUG_OUTPUT_NAME Fortran)

IF (USE_SCALAPACK)
  ADD_LIBRARY(FortranScalapackLib
              f77src/thinsvd.f
              f77src/lworksize.f)
  SET_TARGET_PROPERTIES(FortranScalapackLib
                        PROPERTIES
                        COMPILE_FLAGS -fPIC
                        RELEASE_OUTPUT_NAME FortranScalapack
                        DEBUG_OUTPUT_NAME FortranScalapack)
ENDIF(USE_SCALAPACK)


IF(MPI_FOUND)
  MESSAGE("-- Will compile with MPI API " ${MPI_INCLUDE_PATH} ".")
  SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
  add_definitions(-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK)
  
  IF(NOT USE_SCALAPACK)
    ADD_LIBRARY(fluidlib SHARED Main.C ${CMAKE_CURRENT_BINARY_DIR}/hgidentity.C)
    SET_TARGET_PROPERTIES(fluidlib 
                          PROPERTIES COMPILE_FLAGS -DCREATE_DSO)
    TARGET_LINK_LIBRARIES(fluidlib FluidBase Parser FortranLib Utils 
                          ${EXTRALIBRARIES} ${FORTRAN_HYBRID_LIBRARIES})
    SET_TARGET_PROPERTIES(fluidlib 
                          PROPERTIES 
                          LIBRARY_OUTPUT_DIRECTORY ${BINDIR}
                          RELEASE_OUTPUT_NAME aerof
                          DEBUG_OUTPUT_NAME aerof.debug
                          PREFIX ""
                         )
    SET(TARGETS_INST ${TARGETS_INST} fluidlib)
  ENDIF(NOT USE_SCALAPACK)

ELSE(MPI_FOUND)
  MESSAGE("-- Only a non-MPI based executable will be compiled.")
  SET(MPI_LIBRARIES "")
ENDIF(MPI_FOUND)

# SET(CMAKE_VERBOSE_MAKEFILE true)

#MESSAGE("All include dirs" ${INC_DIR})
INCLUDE_DIRECTORIES( ${INC_DIR} )

#add_subdirectory(parser)
#add_subdirectory(f77src)
#add_subdirectory(utils)
#add_subdirectory(tools)


add_definitions(-DTEMPLATE_FIX -D__LINUX -DF_NEEDS_UNDSC -DTYPE_PREC=float)

IF(NOT USE_SCALAPACK)
  TARGET_LINK_LIBRARIES(fluid FluidBase Parser FortranLib Utils ${EXTRALIBRARIES} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
ELSE(NOT USE_SCALAPACK)
  TARGET_LINK_LIBRARIES(fluid FluidBase Parser FortranLib FortranScalapackLib Utils ${EXTRALIBRARIES} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} -lpthread)
ENDIF(NOT USE_SCALAPACK)
SET_TARGET_PROPERTIES(fluid 
                      PROPERTIES 
		        RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
			RELEASE_OUTPUT_NAME aerof
			DEBUG_OUTPUT_NAME aerof.debug
		      )


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

INSTALL(TARGETS ${TARGETS_INST}
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)

MESSAGE("")
MESSAGE("=================================================")
MESSAGE("           Summary of build options")
MESSAGE("-------------------------------------------------")
IF(MPI_FOUND)
 MESSAGE("Distributed Execution:    YES")
 MESSAGE("Aeroelastic:              YES")
ELSE(MPI_FOUND)
 MESSAGE("Distributed Execution:    NO *******")
 MESSAGE("Aeroelastic:              NO *******")
ENDIF(MPI_FOUND)
IF(DO_EMBEDDED)
 MESSAGE("Embedded framework:       YES")
ELSE(DO_EMBEDDED)
 MESSAGE("Embedded framework:       NO *******")
ENDIF(DO_EMBEDDED)
IF(DO_MODAL)
 MESSAGE("Modal capability:         YES")
ELSE(DO_MODAL)
 MESSAGE("Modal capability:         NO *******")
ENDIF(DO_MODAL)
IF(USE_SCALAPACK)
 MESSAGE("Parallel SVD capability:  YES")
ELSE(USE_SCALAPACK)
 MESSAGE("Parallel SVD capability:  NO *******")
ENDIF(USE_SCALAPACK)
MESSAGE("Build type:               " ${CMAKE_BUILD_TYPE})
MESSAGE("Extra link flags:         " ${EXTRALIB})
MESSAGE("=================================================")
MESSAGE("")

