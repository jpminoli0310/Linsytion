# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jpminoli/Documentos/linsytion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jpminoli/Documentos/linsytion

# Include any dependencies generated for this target.
include CMakeFiles/sal.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/sal.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sal.out.dir/flags.make

CMakeFiles/sal.out.dir/mainProgram.f90.o: CMakeFiles/sal.out.dir/flags.make
CMakeFiles/sal.out.dir/mainProgram.f90.o: mainProgram.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jpminoli/Documentos/linsytion/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/sal.out.dir/mainProgram.f90.o"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/jpminoli/Documentos/linsytion/mainProgram.f90 -o CMakeFiles/sal.out.dir/mainProgram.f90.o

CMakeFiles/sal.out.dir/mainProgram.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/sal.out.dir/mainProgram.f90.i"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/jpminoli/Documentos/linsytion/mainProgram.f90 > CMakeFiles/sal.out.dir/mainProgram.f90.i

CMakeFiles/sal.out.dir/mainProgram.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/sal.out.dir/mainProgram.f90.s"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/jpminoli/Documentos/linsytion/mainProgram.f90 -o CMakeFiles/sal.out.dir/mainProgram.f90.s

CMakeFiles/sal.out.dir/mainProgram.f90.o.requires:

.PHONY : CMakeFiles/sal.out.dir/mainProgram.f90.o.requires

CMakeFiles/sal.out.dir/mainProgram.f90.o.provides: CMakeFiles/sal.out.dir/mainProgram.f90.o.requires
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/mainProgram.f90.o.provides.build
.PHONY : CMakeFiles/sal.out.dir/mainProgram.f90.o.provides

CMakeFiles/sal.out.dir/mainProgram.f90.o.provides.build: CMakeFiles/sal.out.dir/mainProgram.f90.o


CMakeFiles/sal.out.dir/Linsytion.f90.o: CMakeFiles/sal.out.dir/flags.make
CMakeFiles/sal.out.dir/Linsytion.f90.o: Linsytion.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jpminoli/Documentos/linsytion/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/sal.out.dir/Linsytion.f90.o"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/jpminoli/Documentos/linsytion/Linsytion.f90 -o CMakeFiles/sal.out.dir/Linsytion.f90.o

CMakeFiles/sal.out.dir/Linsytion.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/sal.out.dir/Linsytion.f90.i"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/jpminoli/Documentos/linsytion/Linsytion.f90 > CMakeFiles/sal.out.dir/Linsytion.f90.i

CMakeFiles/sal.out.dir/Linsytion.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/sal.out.dir/Linsytion.f90.s"
	/usr/bin/f95  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/jpminoli/Documentos/linsytion/Linsytion.f90 -o CMakeFiles/sal.out.dir/Linsytion.f90.s

CMakeFiles/sal.out.dir/Linsytion.f90.o.requires:

.PHONY : CMakeFiles/sal.out.dir/Linsytion.f90.o.requires

CMakeFiles/sal.out.dir/Linsytion.f90.o.provides: CMakeFiles/sal.out.dir/Linsytion.f90.o.requires
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/Linsytion.f90.o.provides.build
.PHONY : CMakeFiles/sal.out.dir/Linsytion.f90.o.provides

CMakeFiles/sal.out.dir/Linsytion.f90.o.provides.build: CMakeFiles/sal.out.dir/Linsytion.f90.o


# Object files for target sal.out
sal_out_OBJECTS = \
"CMakeFiles/sal.out.dir/mainProgram.f90.o" \
"CMakeFiles/sal.out.dir/Linsytion.f90.o"

# External object files for target sal.out
sal_out_EXTERNAL_OBJECTS =

sal.out: CMakeFiles/sal.out.dir/mainProgram.f90.o
sal.out: CMakeFiles/sal.out.dir/Linsytion.f90.o
sal.out: CMakeFiles/sal.out.dir/build.make
sal.out: CMakeFiles/sal.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jpminoli/Documentos/linsytion/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran executable sal.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sal.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sal.out.dir/build: sal.out

.PHONY : CMakeFiles/sal.out.dir/build

CMakeFiles/sal.out.dir/requires: CMakeFiles/sal.out.dir/mainProgram.f90.o.requires
CMakeFiles/sal.out.dir/requires: CMakeFiles/sal.out.dir/Linsytion.f90.o.requires

.PHONY : CMakeFiles/sal.out.dir/requires

CMakeFiles/sal.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sal.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sal.out.dir/clean

CMakeFiles/sal.out.dir/depend:
	cd /home/jpminoli/Documentos/linsytion && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jpminoli/Documentos/linsytion /home/jpminoli/Documentos/linsytion /home/jpminoli/Documentos/linsytion /home/jpminoli/Documentos/linsytion /home/jpminoli/Documentos/linsytion/CMakeFiles/sal.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sal.out.dir/depend

