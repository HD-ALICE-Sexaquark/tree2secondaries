cmake_minimum_required(VERSION 3.15...3.31)
project(Tree2Secondaries LANGUAGES CXX)

# support folders in IDEs #

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

include(GNUInstallDirs) # GNU-like installation paths, e.g. lib/, include/, ...

# C++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -O3 -DNDEBUG -march=native -mtune=native"
)

# options

set(CMAKE_EXPORT_COMPILE_COMMANDS True)

option(T2S_DEBUG "Print debug messages." OFF)
if(T2S_DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DT2S_DEBUG")
endif(T2S_DEBUG)

option(T2S_USE_ALICE "Use ALICE Offline V0 Finder instead of KFParticle." OFF)
if(T2S_USE_ALICE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DT2S_USE_ALICE")
endif(T2S_USE_ALICE)

# enable LTO #

include(CheckIPOSupported)
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
find_package(ROOT 6.32 CONFIG REQUIRED COMPONENTS Core RIO Tree EG)
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

include(FetchContent)

message(STATUS "Fetching https://github.com/CLIUtils/CLI11")
FetchContent_Declare(
  cli11_proj
  QUIET
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.5.0)
FetchContent_MakeAvailable(cli11_proj)
message(STATUS "Fetching https://github.com/CLIUtils/CLI11 -- fetched")

message(STATUS "Fetching https://github.com/fmtlib/fmt")
FetchContent_Declare(
  fmt
  QUIET
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 11.2.0)
FetchContent_MakeAvailable(fmt)
message(STATUS "Fetching https://github.com/fmtlib/fmt -- fetched")

# get ready to build target #

message(STATUS "Using compiler flags: ${CMAKE_CXX_FLAGS}")

# main targets of the project in various subdirectories. order matters. #

add_subdirectory(src)
