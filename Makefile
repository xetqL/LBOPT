# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/xetql/lb_exhaustive_search

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xetql/lb_exhaustive_search

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/xetql/lb_exhaustive_search/CMakeFiles /home/xetql/lb_exhaustive_search/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/xetql/lb_exhaustive_search/CMakeFiles 0
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
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named lb_exhaustive_search

# Build rule for target.
lb_exhaustive_search: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 lb_exhaustive_search
.PHONY : lb_exhaustive_search

# fast build rule for target.
lb_exhaustive_search/fast:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/build
.PHONY : lb_exhaustive_search/fast

lbnode.o: lbnode.cpp.o

.PHONY : lbnode.o

# target to build an object file
lbnode.cpp.o:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/lbnode.cpp.o
.PHONY : lbnode.cpp.o

lbnode.i: lbnode.cpp.i

.PHONY : lbnode.i

# target to preprocess a source file
lbnode.cpp.i:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/lbnode.cpp.i
.PHONY : lbnode.cpp.i

lbnode.s: lbnode.cpp.s

.PHONY : lbnode.s

# target to generate assembly for a file
lbnode.cpp.s:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/lbnode.cpp.s
.PHONY : lbnode.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/main.cpp.s
.PHONY : main.cpp.s

utils.o: utils.cpp.o

.PHONY : utils.o

# target to build an object file
utils.cpp.o:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/utils.cpp.o
.PHONY : utils.cpp.o

utils.i: utils.cpp.i

.PHONY : utils.i

# target to preprocess a source file
utils.cpp.i:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/utils.cpp.i
.PHONY : utils.cpp.i

utils.s: utils.cpp.s

.PHONY : utils.s

# target to generate assembly for a file
utils.cpp.s:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/utils.cpp.s
.PHONY : utils.cpp.s

zupply/src/zupply.o: zupply/src/zupply.cpp.o

.PHONY : zupply/src/zupply.o

# target to build an object file
zupply/src/zupply.cpp.o:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/zupply/src/zupply.cpp.o
.PHONY : zupply/src/zupply.cpp.o

zupply/src/zupply.i: zupply/src/zupply.cpp.i

.PHONY : zupply/src/zupply.i

# target to preprocess a source file
zupply/src/zupply.cpp.i:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/zupply/src/zupply.cpp.i
.PHONY : zupply/src/zupply.cpp.i

zupply/src/zupply.s: zupply/src/zupply.cpp.s

.PHONY : zupply/src/zupply.s

# target to generate assembly for a file
zupply/src/zupply.cpp.s:
	$(MAKE) -f CMakeFiles/lb_exhaustive_search.dir/build.make CMakeFiles/lb_exhaustive_search.dir/zupply/src/zupply.cpp.s
.PHONY : zupply/src/zupply.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... lb_exhaustive_search"
	@echo "... edit_cache"
	@echo "... lbnode.o"
	@echo "... lbnode.i"
	@echo "... lbnode.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... utils.o"
	@echo "... utils.i"
	@echo "... utils.s"
	@echo "... zupply/src/zupply.o"
	@echo "... zupply/src/zupply.i"
	@echo "... zupply/src/zupply.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

