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
CMAKE_SOURCE_DIR = /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build

# Include any dependencies generated for this target.
include CMakeFiles/process.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/process.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/process.dir/flags.make

CMakeFiles/process.dir/process.cpp.o: CMakeFiles/process.dir/flags.make
CMakeFiles/process.dir/process.cpp.o: ../process.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/process.dir/process.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/process.dir/process.cpp.o -c /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/process.cpp

CMakeFiles/process.dir/process.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/process.dir/process.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/process.cpp > CMakeFiles/process.dir/process.cpp.i

CMakeFiles/process.dir/process.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/process.dir/process.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/process.cpp -o CMakeFiles/process.dir/process.cpp.s

CMakeFiles/process.dir/process.cpp.o.requires:

.PHONY : CMakeFiles/process.dir/process.cpp.o.requires

CMakeFiles/process.dir/process.cpp.o.provides: CMakeFiles/process.dir/process.cpp.o.requires
	$(MAKE) -f CMakeFiles/process.dir/build.make CMakeFiles/process.dir/process.cpp.o.provides.build
.PHONY : CMakeFiles/process.dir/process.cpp.o.provides

CMakeFiles/process.dir/process.cpp.o.provides.build: CMakeFiles/process.dir/process.cpp.o


# Object files for target process
process_OBJECTS = \
"CMakeFiles/process.dir/process.cpp.o"

# External object files for target process
process_EXTERNAL_OBJECTS =

process: CMakeFiles/process.dir/process.cpp.o
process: CMakeFiles/process.dir/build.make
process: CMakeFiles/process.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable process"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/process.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/process.dir/build: process

.PHONY : CMakeFiles/process.dir/build

CMakeFiles/process.dir/requires: CMakeFiles/process.dir/process.cpp.o.requires

.PHONY : CMakeFiles/process.dir/requires

CMakeFiles/process.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/process.dir/cmake_clean.cmake
.PHONY : CMakeFiles/process.dir/clean

CMakeFiles/process.dir/depend:
	cd /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build /home/john/Dropbox/Psychology/Neuroscience/Ermentrout2009/processes/sim/build/CMakeFiles/process.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/process.dir/depend

