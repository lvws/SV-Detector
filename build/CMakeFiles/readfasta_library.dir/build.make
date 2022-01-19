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
include CMakeFiles/readfasta_library.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/readfasta_library.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/readfasta_library.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/readfasta_library.dir/flags.make

CMakeFiles/readfasta_library.dir/src/reference.o: CMakeFiles/readfasta_library.dir/flags.make
CMakeFiles/readfasta_library.dir/src/reference.o: ../src/reference.cpp
CMakeFiles/readfasta_library.dir/src/reference.o: CMakeFiles/readfasta_library.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/readfasta_library.dir/src/reference.o"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/readfasta_library.dir/src/reference.o -MF CMakeFiles/readfasta_library.dir/src/reference.o.d -o CMakeFiles/readfasta_library.dir/src/reference.o -c /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/reference.cpp

CMakeFiles/readfasta_library.dir/src/reference.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/readfasta_library.dir/src/reference.i"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/reference.cpp > CMakeFiles/readfasta_library.dir/src/reference.i

CMakeFiles/readfasta_library.dir/src/reference.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/readfasta_library.dir/src/reference.s"
	/dssg/home/lvws/programs/miniconda3/envs/conda_gcc_env/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/src/reference.cpp -o CMakeFiles/readfasta_library.dir/src/reference.s

# Object files for target readfasta_library
readfasta_library_OBJECTS = \
"CMakeFiles/readfasta_library.dir/src/reference.o"

# External object files for target readfasta_library
readfasta_library_EXTERNAL_OBJECTS =

libreadfasta_library.so: CMakeFiles/readfasta_library.dir/src/reference.o
libreadfasta_library.so: CMakeFiles/readfasta_library.dir/build.make
libreadfasta_library.so: CMakeFiles/readfasta_library.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libreadfasta_library.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/readfasta_library.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/readfasta_library.dir/build: libreadfasta_library.so
.PHONY : CMakeFiles/readfasta_library.dir/build

CMakeFiles/readfasta_library.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/readfasta_library.dir/cmake_clean.cmake
.PHONY : CMakeFiles/readfasta_library.dir/clean

CMakeFiles/readfasta_library.dir/depend:
	cd /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build /dssg/home/lvws/GS-WanShengLv/gscap/Tools/SV-Detector/build/CMakeFiles/readfasta_library.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/readfasta_library.dir/depend
