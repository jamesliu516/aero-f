# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lustre/home/ibuenros/Fluid

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lustre/home/ibuenros/Fluid

# Include any dependencies generated for this target.
include CMakeFiles/loader.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/loader.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/loader.dir/flags.make

CMakeFiles/loader.dir/tools/loader.C.o: CMakeFiles/loader.dir/flags.make
CMakeFiles/loader.dir/tools/loader.C.o: tools/loader.C
	$(CMAKE_COMMAND) -E cmake_progress_report /lustre/home/ibuenros/Fluid/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/loader.dir/tools/loader.C.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/loader.dir/tools/loader.C.o -c /lustre/home/ibuenros/Fluid/tools/loader.C

CMakeFiles/loader.dir/tools/loader.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/loader.dir/tools/loader.C.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lustre/home/ibuenros/Fluid/tools/loader.C > CMakeFiles/loader.dir/tools/loader.C.i

CMakeFiles/loader.dir/tools/loader.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/loader.dir/tools/loader.C.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lustre/home/ibuenros/Fluid/tools/loader.C -o CMakeFiles/loader.dir/tools/loader.C.s

CMakeFiles/loader.dir/tools/loader.C.o.requires:
.PHONY : CMakeFiles/loader.dir/tools/loader.C.o.requires

CMakeFiles/loader.dir/tools/loader.C.o.provides: CMakeFiles/loader.dir/tools/loader.C.o.requires
	$(MAKE) -f CMakeFiles/loader.dir/build.make CMakeFiles/loader.dir/tools/loader.C.o.provides.build
.PHONY : CMakeFiles/loader.dir/tools/loader.C.o.provides

CMakeFiles/loader.dir/tools/loader.C.o.provides.build: CMakeFiles/loader.dir/tools/loader.C.o
.PHONY : CMakeFiles/loader.dir/tools/loader.C.o.provides.build

# Object files for target loader
loader_OBJECTS = \
"CMakeFiles/loader.dir/tools/loader.C.o"

# External object files for target loader
loader_EXTERNAL_OBJECTS =

bin/loader: CMakeFiles/loader.dir/tools/loader.C.o
bin/loader: /usr/mpi/gcc/openmpi-1.2.6/lib64/libmpi_cxx.so
bin/loader: /usr/mpi/gcc/openmpi-1.2.6/lib64/libmpi.so
bin/loader: /usr/mpi/gcc/openmpi-1.2.6/lib64/libopen-rte.so
bin/loader: /usr/mpi/gcc/openmpi-1.2.6/lib64/libopen-pal.so
bin/loader: /usr/lib64/libdl.so
bin/loader: /usr/lib64/libnsl.so
bin/loader: /usr/lib64/libutil.so
bin/loader: /usr/lib64/libm.so
bin/loader: /usr/lib64/libdl.so
bin/loader: /usr/lib64/libnsl.so
bin/loader: /usr/lib64/libutil.so
bin/loader: /usr/lib64/libm.so
bin/loader: CMakeFiles/loader.dir/build.make
bin/loader: CMakeFiles/loader.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable bin/loader"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/loader.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/loader.dir/build: bin/loader
.PHONY : CMakeFiles/loader.dir/build

CMakeFiles/loader.dir/requires: CMakeFiles/loader.dir/tools/loader.C.o.requires
.PHONY : CMakeFiles/loader.dir/requires

CMakeFiles/loader.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/loader.dir/cmake_clean.cmake
.PHONY : CMakeFiles/loader.dir/clean

CMakeFiles/loader.dir/depend:
	cd /lustre/home/ibuenros/Fluid && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lustre/home/ibuenros/Fluid /lustre/home/ibuenros/Fluid /lustre/home/ibuenros/Fluid /lustre/home/ibuenros/Fluid /lustre/home/ibuenros/Fluid/CMakeFiles/loader.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/loader.dir/depend

