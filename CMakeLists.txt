cmake_minimum_required(VERSION 3.25...3.31)
project(
  RNTuple2DoubleStrangeness
  VERSION 26.05.24
  LANGUAGES CXX)

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

# enable LTO #

check_ipo_supported(RESULT LTO_SUPPORTED)
if(LTO_SUPPORTED)
  message(STATUS "LTO/IPO is supported")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(STATUS "LTO/IPO is not supported")
endif()

# option: profiling #

option(ENABLE_PROFILING "Enable profiling with gprof" OFF)
if(ENABLE_PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

option(LEGACY_KF "Enable original KFParticle functions" OFF)

# find dependencies #

message(STATUS "Fetching https://github.com/CLIUtils/CLI11")
FetchContent_Declare(
  cli11_proj
  QUIET
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.6.2)
FetchContent_MakeAvailable(cli11_proj)
message(STATUS "Fetching https://github.com/CLIUtils/CLI11 -- fetched")

message(STATUS "Fetching https://gitlab.com/libeigen/eigen")
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 5.0.1)
FetchContent_MakeAvailable(eigen)
message(STATUS "Fetching https://gitlab.com/libeigen/eigen -- fetched")

message(STATUS "Looking for ROOT")
find_package(ROOT 6.36 REQUIRED CONFIG COMPONENTS Core Hist RIO EG ROOTNTuple)
if(ROOT_FOUND)
  message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS}")
endif()

# main target #

add_executable(
  ${PROJECT_NAME}
  $<$<BOOL:${LEGACY_KF}>:
  src/Legacy/LegacyAddDaughter.cxx
  src/Legacy/LegacyAddDaughterWithEnergyFit.cxx
  src/Legacy/LegacyGetMeasurement.cxx
  src/Legacy/LegacyHelixHelix.cxx
  src/Legacy/LegacyHelixPoint.cxx
  src/Legacy/LegacyInvertCholesky.cxx
  src/Legacy/LegacyLineLine.cxx
  src/Legacy/LegacyLinePoint.cxx
  src/Legacy/LegacyMultQSQt.cxx
  src/Legacy/LegacyParticle.cxx
  src/Legacy/LegacyTransportBz.cxx
  src/Legacy/LegacyTransportLine.cxx
  >
  src/Seeder/SeederHelixHelix.cxx
  src/Seeder/SeederHelixLine.cxx
  src/Seeder/SeederHelixVertex.cxx
  src/Seeder/SeederLineLine.cxx
  src/Seeder/SeederLineVertex.cxx
  src/KalmanFitter/KalmanFitter.cxx
  src/KalmanFitter/KalmanFitterParticle.cxx
  src/Packager.cxx
  src/Finder.cxx
  src/Verifier.cxx
  src/App.cxx)

target_link_libraries(
  ${PROJECT_NAME} PRIVATE CLI11::CLI11 Eigen3::Eigen ROOT::Core ROOT::RIO
                          ROOT::Hist ROOT::ROOTNTuple)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_compile_options(
  ${PROJECT_NAME}
  PRIVATE -Wall
          -Wextra
          -Wshadow
          -Wnon-virtual-dtor
          -Wpedantic
          -Wconversion
          $<$<CONFIG:Release>:-mtune=native
          -march=native
          -ffp-contract=fast
          -fno-math-errno
          -mfma>)

target_compile_definitions(
  ${PROJECT_NAME}
  PUBLIC EIGEN_DONT_PARALLELIZE $<$<BOOL:${LEGACY_KF}>:R2DS_LEGACY_KF>
         $<$<CONFIG:Debug>:R2DS_DEBUG> $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)
