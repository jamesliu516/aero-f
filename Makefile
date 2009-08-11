# This make file does not need to be changed if you are compiling on lancer or thunderbird or nivation

#uncoment this to use the intel compiler
#compiler = intel

PREFIX     = fluid
THESYS     = `/bin/uname -s`
EXE        = bin/$(PREFIX).opt.$(THESYS)
DSOEXE     = bin/$(PREFIX).so.$(THESYS)
THEDATE    = `date "+%b-%d-%Y"`
THEARCHIVE = $(PREFIX).$(THEDATE).tar
THEFILES   = `find . \( -name '*.C' -o -name '*.c' -o -name '*.f' -o -name '*.a' -o\
		        -name '*.h' -o -name '*.l' -o -name '*.y' -o \
		        -name Makefile -o -name CVS -o -name '*.defs' -o \
			-name README \)`
%.dep: %.C
	@set -e; rm -f $@; \
         $(CXX) -M $(CXXFLAGS) $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
         rm -f $@.$$$$

%.dep: %.f
	@set -e; touch $@

CXXOBJS    = BcFcnCore.o \
             CommunicatorCore.o \
             CommunicatorMPI.o \
             ConnectivityCore.o \
             CorotSolver.o \
             CurvatureDetection.o \
             DistGeoState.o \
             DistMacroCellCore.o \
             DomainCore.o \
             EdgeCore.o \
             FaceCore.o \
             FaceTriaCore.o \
             ElemCore.o \
             ElemTetCore.o \
             FemEquationTermDesc.o \
             FluxFcnDesc.o \
             FluxFcnDescPerfectGas.o \
             FluxFcnDescWaterCompressible.o \
             FluxFcnDescJWL.o \
             FluxFcnDescGasInGas.o \
             FluxFcnDescJWLInGas.o \
             FluxFcnDescLiquidInLiquid.o \
             FluxFcnDescGasInLiquid.o \
             GeoData.o \
             GeoSource.o \
             GeoState.o \
             HeatTransferHandlerCore.o \
             InletNodeCore.o \
             IoDataCore.o \
             KspConvCriterion.o \
             LevelSetCore.o \
             MacroCellCore.o \
             MatchNodeCore.o \
             MemoryPool.o \
             MeshMotionHandlerCore.o \
             MeshMotionSolver.o \
             ModalSolver.o \
             NavierStokesSolver.o \
             Node.o \
             PostFcn.o \
             RefVal.o \
             SmagorinskyLESTerm.o \
             WaleLESTerm.o \
             DynamicLESTerm.o \
             StructExc.o \
             SubDomainCore.o \
             TimeData.o \
             Timer.o \
             TsInput.o \
             TsParameters.o \
             TsRestartCore.o \
             VarFcn.o \
             VMSLESTerm.o \
             DynamicVMSTerm.o \
             WallFcnCore.o \
             BCApplierCore.o \
	     BCond.o \
	     SparseGridCore.o \
	     BlockAlloc.o


# INTEL
# -----
#CXX = icc -Kc++
#CXX = mpiCC -Kc++
#F77 = ifc -w95
# SUN
# ---
#CXX = g++
# ASCI RED
# --------
#CC  = cicc
#CXX = ciCC
#CCC = ciCC
#F77 = cif77
#FC  = cif77

# default defs - (defaults assume Linux OS)
ARPACKDIR = /lib
OMPFLAGS  = # OpenMP flags (defined in machine-specific def file) to be added to the OMPFLAG variable
FORTLIBS  = # Fortran libraries needed for linking when compiled with the Intel C/C++ compiler
CXX       = g++
#include default.defs
#include nivation.defs

#usescalapack = false
usescalapack = true
ifeq ($(usescalapack),true)
  SFLAGS       = -DDO_SCALAPACK
else
  SFLAGS       =
endif


host := $(shell hostname -s)
ifeq ($(host), frontend-0)
  host = nivation
  CXX = g++
endif

ifeq ($(host), nivation)
  host = nivation
  CXX = g++
endif

ifeq ($(host), regelation)
  host = 64bit
  CXX = g++
endif

ifeq ($(findstring nas-0-, $(host)), nas-0-)
  host = nivation
  CXX = g++
endif

ifeq ($(findstring compute-2-, $(host)), compute-2-)
  host = 64bit
  CXX = g++
endif

ifeq ($(host), thunderbird)
  host = thunderbird
  CXX = CC
endif


# Look for machine-specific definitions
#finddefs = $(shell ls $(host).defs)
#ifeq ($(findstring $(host).defs, $(finddefs)), $(host).defs)
#  include $(host).defs
#endif

deffile = $(host).defs
finddefs = $(shell ls $(deffile))

ifeq ($(findstring $(host), $(finddefs)), $(host))
  include $(host).defs
else
  include default.defs
endif

# search for arpack directory and make appropriate definitons
findarpack = $(shell ls $(ARPACKLIB))

ifeq ($(findstring $(ARPACKLIB), $(findarpack)), $(ARPACKLIB))
  DFLAGS = -DF_NEEDS_UNDSC -DTYPE_PREC=float -DDO_MODAL
else
  ARPACKLIB = 
  DFLAGS = -DF_NEEDS_UNDSC -DTYPE_PREC=float
endif

DYNLINKFLAGS = -rdynamic --export-dynamic

MPIFLAG = -DUSE_MPI -DMPI_NO_CPPBIND \
          -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 \
          -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 -DHAVE_MPI_CPP -fexceptions \
          #-I/opt/mpich/myrinet/gnu/include/mpi2c++ -fexceptions -I/opt/mpich/myrinet/gnu/include

OMPFLAG   = $(OMPFLAGS) #-mp
CCFLAGS   = $(CXXFLAGS)

default: $(CXXOBJS) f77src parser utils
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -c Main.C
	$(CXX) $(CXXFLAGS) -o $(EXE) Main.o $(CXXOBJS) -Lf77src -Lparser -lf77src -lparser -Lutils -lutils $(LIBS) $(MPLIBS) $(FORTLIBS)

fluid.so: $(CXXOBJS) f77src parser utils
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -DCREATE_DSO -c Main.C
	$(LD) $(CXXFLAGS) $(SOFLAGS) -o $(DSOEXE) Main.o $(CXXOBJS) -Lf77src -Lparser -lf77src -lparser -Lutils -lutils $(LIBS) $(MPLIBS) $(FORTLIBS)

f77src::
	(cd f77src; $(MAKE) "FC=$(F77)" "F77=$(F77)" "FFLAGS=$(OFLAGS)" "usescalapack=$(usescalapack)")

parser::
	(cd parser; $(MAKE) "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "SAR=$(SAR)")

tools::
	(cd tools; $(MAKE) "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "MPLIBS=$(MPLIBS)" "DLLIBS=$(DLLIBS)" "DYNLINKFLAGS=$(DYNLINKFLAGS)")

utils::
	 (cd utils; $(MAKE) "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) -ffloat-store")

all:
	make
	make fluid.so
	make tools

clean: 
	rm -f $(CXXOBJS) $(EXE) $(DSOEXE) so_locations

allclean: 
	make clean
	(cd f77src; $(MAKE) clean)
	(cd parser; $(MAKE) clean)
	(cd tools; $(MAKE) clean)
	(cd utils; $(MAKE) clean)

prune:
	rm -f */*~

tar:
	@if (test -f $(THEARCHIVE)) ;\
	then mv -f $(THEARCHIVE) $(THEARCHIVE)~ ;\
	fi
	@if (test -f $(THEARCHIVE).gz) ;\
	then mv -f $(THEARCHIVE).gz $(THEARCHIVE).gz~ ;\
	fi
	@tar cvf $(THEARCHIVE) $(THEFILES)
	@chmod 600 $(THEARCHIVE)
	@gzip $(THEARCHIVE)

include $(CXXOBJS:.o=.dep)
