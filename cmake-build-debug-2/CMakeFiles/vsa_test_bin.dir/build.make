# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/vector_cat/gits/OTMapping

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/vector_cat/gits/OTMapping/cmake-build-debug-2

# Include any dependencies generated for this target.
include CMakeFiles/vsa_test_bin.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vsa_test_bin.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vsa_test_bin.dir/flags.make

CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o: CMakeFiles/vsa_test_bin.dir/flags.make
CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o: ../src/vsa_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/vector_cat/gits/OTMapping/cmake-build-debug-2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o -c /Users/vector_cat/gits/OTMapping/src/vsa_test.cpp

CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/vector_cat/gits/OTMapping/src/vsa_test.cpp > CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.i

CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/vector_cat/gits/OTMapping/src/vsa_test.cpp -o CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.s

# Object files for target vsa_test_bin
vsa_test_bin_OBJECTS = \
"CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o"

# External object files for target vsa_test_bin
vsa_test_bin_EXTERNAL_OBJECTS =

vsa_test_bin: CMakeFiles/vsa_test_bin.dir/src/vsa_test.cpp.o
vsa_test_bin: CMakeFiles/vsa_test_bin.dir/build.make
vsa_test_bin: libotmapping.a
vsa_test_bin: /usr/local/lib/libboost_thread-mt.a
vsa_test_bin: /usr/local/lib/libboost_system-mt.a
vsa_test_bin: /usr/local/lib/libboost_chrono-mt.a
vsa_test_bin: /usr/local/lib/libboost_date_time-mt.a
vsa_test_bin: /usr/local/lib/libboost_atomic-mt.a
vsa_test_bin: /usr/local/lib/libmpfr.dylib
vsa_test_bin: /usr/local/lib/libgmp.dylib
vsa_test_bin: /usr/local/lib/libboost_thread-mt.a
vsa_test_bin: /usr/local/lib/libboost_chrono-mt.a
vsa_test_bin: /usr/local/lib/libboost_system-mt.a
vsa_test_bin: /usr/local/lib/libboost_date_time-mt.a
vsa_test_bin: /usr/local/lib/libboost_atomic-mt.a
vsa_test_bin: /System/Library/Frameworks/OpenGL.framework/OpenGL
vsa_test_bin: imgui/libimgui.a
vsa_test_bin: glfw/src/libglfw3.a
vsa_test_bin: glad/libglad.a
vsa_test_bin: libtinyexpr.a
vsa_test_bin: CMakeFiles/vsa_test_bin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/vector_cat/gits/OTMapping/cmake-build-debug-2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vsa_test_bin"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vsa_test_bin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vsa_test_bin.dir/build: vsa_test_bin

.PHONY : CMakeFiles/vsa_test_bin.dir/build

CMakeFiles/vsa_test_bin.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vsa_test_bin.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vsa_test_bin.dir/clean

CMakeFiles/vsa_test_bin.dir/depend:
	cd /Users/vector_cat/gits/OTMapping/cmake-build-debug-2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vector_cat/gits/OTMapping /Users/vector_cat/gits/OTMapping /Users/vector_cat/gits/OTMapping/cmake-build-debug-2 /Users/vector_cat/gits/OTMapping/cmake-build-debug-2 /Users/vector_cat/gits/OTMapping/cmake-build-debug-2/CMakeFiles/vsa_test_bin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vsa_test_bin.dir/depend

