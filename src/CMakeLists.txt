add_executable(App Packager.cxx Finder.cxx App.cxx)

target_link_libraries(
  App
  PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist KFParticle::KFParticle
  PRIVATE CLI11::CLI11)

target_include_directories(
  App PUBLIC ${ROOT_INCLUDE_DIRS} ${KFParticle_INCLUDE_DIR}
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

target_compile_features(App PUBLIC cxx_std_23)
target_compile_options(App PRIVATE -Wall -Wextra -Wshadow -Wnon-virtual-dtor)
target_compile_options(App PRIVATE -O3 -DNDEBUG -mtune=native -march=native
                                   -ffp-contract=fast)

option(T2S_DEBUG "Print debug messages." OFF)
if(T2S_DEBUG)
  target_compile_definitions(App PUBLIC T2S_DEBUG)
endif(T2S_DEBUG)

# print target compile options #

cmake_print_properties(TARGETS App PROPERTIES COMPILE_FEATURES COMPILE_OPTIONS
                                              COMPILE_DEFINITIONS)
