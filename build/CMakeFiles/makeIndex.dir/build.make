# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/cmake

# The command to remove a file.
RM = /dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build

# Include any dependencies generated for this target.
include CMakeFiles/makeIndex.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/makeIndex.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/makeIndex.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/makeIndex.dir/flags.make

CMakeFiles/makeIndex.dir/src/makeIndex.o: CMakeFiles/makeIndex.dir/flags.make
CMakeFiles/makeIndex.dir/src/makeIndex.o: ../src/makeIndex.cpp
CMakeFiles/makeIndex.dir/src/makeIndex.o: CMakeFiles/makeIndex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/makeIndex.dir/src/makeIndex.o"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/makeIndex.dir/src/makeIndex.o -MF CMakeFiles/makeIndex.dir/src/makeIndex.o.d -o CMakeFiles/makeIndex.dir/src/makeIndex.o -c /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/makeIndex.cpp

CMakeFiles/makeIndex.dir/src/makeIndex.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/makeIndex.dir/src/makeIndex.i"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/makeIndex.cpp > CMakeFiles/makeIndex.dir/src/makeIndex.i

CMakeFiles/makeIndex.dir/src/makeIndex.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/makeIndex.dir/src/makeIndex.s"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/makeIndex.cpp -o CMakeFiles/makeIndex.dir/src/makeIndex.s

# Object files for target makeIndex
makeIndex_OBJECTS = \
"CMakeFiles/makeIndex.dir/src/makeIndex.o"

# External object files for target makeIndex
makeIndex_EXTERNAL_OBJECTS =

makeIndex: CMakeFiles/makeIndex.dir/src/makeIndex.o
makeIndex: CMakeFiles/makeIndex.dir/build.make
makeIndex: libreadfasta_library.so
makeIndex: /dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env/lib/libboost_serialization.so
makeIndex: CMakeFiles/makeIndex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable makeIndex"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/makeIndex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/makeIndex.dir/build: makeIndex
.PHONY : CMakeFiles/makeIndex.dir/build

CMakeFiles/makeIndex.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/makeIndex.dir/cmake_clean.cmake
.PHONY : CMakeFiles/makeIndex.dir/clean

CMakeFiles/makeIndex.dir/depend:
	cd /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles/makeIndex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/makeIndex.dir/depend

