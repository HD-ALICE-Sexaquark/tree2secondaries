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

# get ready to build target #

message(STATUS "Using compiler flags: ${CMAKE_CXX_FLAGS}")

# main targets of the project in various subdirectories. order matters. #

add_subdirectory(src)
