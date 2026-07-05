cmake_minimum_required(VERSION 3.28...3.31)
project(
  t2ds
  VERSION 26.07.05
  LANGUAGES CXX)

# set "Release" as default build type #

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# prevent in-source builds #

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# C++ standard #

set(CMAKE_CXX_EXTENSIONS OFF) # ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# includes and macros #

include(CheckIPOSupported)
include(FetchContent)

# enable LTO #

check_ipo_supported(RESULT LTO_SUPPORTED)
if(LTO_SUPPORTED)
  message(STATUS "LTO/IPO is supported")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
else()
  message(STATUS "LTO/IPO is not supported")
endif()

option(LEGACY_KF "Enable original KFParticle functions" OFF)

# find dependencies #

message(STATUS "Fetching https://github.com/CLIUtils/CLI11")
FetchContent_Declare(
  cli11_proj
  SYSTEM
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.6.2
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(cli11_proj)
message(STATUS "Fetching https://github.com/CLIUtils/CLI11 -- fetched")

message(STATUS "Fetching https://gitlab.com/libeigen/eigen")
FetchContent_Declare(
  eigen
  SYSTEM
  EXCLUDE_FROM_ALL
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 5.0.1
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)
message(STATUS "Fetching https://gitlab.com/libeigen/eigen -- fetched")

message(STATUS "Looking for ROOT")
find_package(ROOT 6.40 REQUIRED CONFIG COMPONENTS Core Tree MathCore GenVector Hist RIO ROOTNTuple)
if(ROOT_FOUND)
  message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS}")
endif()

# flags #

set(T2DS_RELEASE_FLAGS -march=x86-64-v3 -mtune=native -ffp-contract=fast -fno-math-errno -fno-trapping-math)

# option: profiling #

option(ENABLE_PROFILING "Enable profiling with gprof" OFF)
if(ENABLE_PROFILING)
  # these affect only targets created after the call
  add_compile_options(-pg)
  add_link_options(-pg)
endif()

# root dictionary #

add_library(POD SHARED)
root_generate_dictionary(G__POD common/Constants.hpp
                                common/POD_Event.hpp
                                common/POD_InjectedHdib.hpp
                                common/POD_InjectedSexa.hpp
                                common/POD_LambdaPair.hpp
                                common/POD_McParticle.hpp
                                common/POD_PreFoundLambda.hpp
                                common/POD_Sexaquark.hpp
                                common/POD_Track.hpp
                                common/POD_V0.hpp
                                LINKDEF common/POD_LinkDef.h
                                MODULE POD)

target_link_libraries(POD PUBLIC ROOT::Core ROOT::RIO)
target_link_options(POD PRIVATE -Wno-alloc-size-larger-than)
target_include_directories(POD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                      ${CMAKE_CURRENT_SOURCE_DIR}/common)

# main target #

add_executable(${PROJECT_NAME} src/Seeder/SeederHelixHelix.cxx
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

if(LEGACY_KF)
  # add sources to a target
  target_sources(${PROJECT_NAME} PRIVATE src/Legacy/LegacyAddDaughter.cxx
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
                                         src/Legacy/LegacyTransportLine.cxx)
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE POD
                                              CLI11::CLI11
                                              Eigen3::Eigen
                                              ROOT::Core ROOT::Tree ROOT::MathCore ROOT::GenVector
                                              ROOT::RIO ROOT::Hist ROOT::ROOTNTuple)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                                   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall
                                               -Wextra
                                               -Wshadow
                                               -Wnon-virtual-dtor
                                               -Wpedantic
                                               -Wconversion
                                               -Wdouble-promotion
                                               $<$<CONFIG:Release>:${T2DS_OPT_FLAGS}>)

target_compile_definitions(${PROJECT_NAME} PRIVATE EIGEN_DONT_PARALLELIZE
                                                   $<$<BOOL:${LEGACY_KF}>:T2DS_LEGACY_KF>
                                                   $<$<CONFIG:Debug>:T2DS_DEBUG>
                                                   $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)
