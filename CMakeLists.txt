cmake_minimum_required(VERSION 3.25...3.31)
project(Tree2Secondaries LANGUAGES CXX)

# support folders in IDEs #

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# C++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# includes and macros #

include(CheckIPOSupported)
include(FetchContent)
include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...
include(CMakePrintHelpers)

# enable LTO #

check_ipo_supported(RESULT LTO_SUPPORTED)
if(LTO_SUPPORTED)
  message(STATUS "LTO/IPO is supported")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "LTO/IPO is not supported")
endif()

# profiling #

option(ENABLE_PROFILING "Enable profiling with gprof" OFF)
if(ENABLE_PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

# find dependencies #

message(STATUS "Looking for ROOT")
find_package(ROOT 6.36 CONFIG REQUIRED COMPONENTS Core RIO Tree EG)
if(ROOT_FOUND)
  message(STATUS "Looking for ROOT -- found")
  message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")
endif()

message(STATUS "Looking for KFParticle")
find_package(KFParticle)
if(KFParticle_FOUND)
  message(STATUS "Looking for KFParticle -- found")
  message(STATUS "KFParticle found at ${KFParticle_LIB_DIR}")
endif()

message(STATUS "Fetching https://github.com/CLIUtils/CLI11")
FetchContent_Declare(
  cli11_proj
  QUIET
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.5.0)
FetchContent_MakeAvailable(cli11_proj)
message(STATUS "Fetching https://github.com/CLIUtils/CLI11 -- fetched")

# main targets of the project in various subdirectories. order matters. #

add_subdirectory(src)
