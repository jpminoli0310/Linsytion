# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jpminoli/Documentos/linsytion/CMakeFiles /home/jpminoli/Documentos/linsytion/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jpminoli/Documentos/linsytion/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named sal.out

# Build rule for target.
sal.out: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sal.out
.PHONY : sal.out

# fast build rule for target.
sal.out/fast:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/build
.PHONY : sal.out/fast

Linsytion.o: Linsytion.f90.o

.PHONY : Linsytion.o

# target to build an object file
Linsytion.f90.o:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/Linsytion.f90.o
.PHONY : Linsytion.f90.o

Linsytion.i: Linsytion.f90.i

.PHONY : Linsytion.i

# target to preprocess a source file
Linsytion.f90.i:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/Linsytion.f90.i
.PHONY : Linsytion.f90.i

Linsytion.s: Linsytion.f90.s

.PHONY : Linsytion.s

# target to generate assembly for a file
Linsytion.f90.s:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/Linsytion.f90.s
.PHONY : Linsytion.f90.s

mainProgram.o: mainProgram.f90.o

.PHONY : mainProgram.o

# target to build an object file
mainProgram.f90.o:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/mainProgram.f90.o
.PHONY : mainProgram.f90.o

mainProgram.i: mainProgram.f90.i

.PHONY : mainProgram.i

# target to preprocess a source file
mainProgram.f90.i:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/mainProgram.f90.i
.PHONY : mainProgram.f90.i

mainProgram.s: mainProgram.f90.s

.PHONY : mainProgram.s

# target to generate assembly for a file
mainProgram.f90.s:
	$(MAKE) -f CMakeFiles/sal.out.dir/build.make CMakeFiles/sal.out.dir/mainProgram.f90.s
.PHONY : mainProgram.f90.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... sal.out"
	@echo "... Linsytion.o"
	@echo "... Linsytion.i"
	@echo "... Linsytion.s"
	@echo "... mainProgram.o"
	@echo "... mainProgram.i"
	@echo "... mainProgram.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

