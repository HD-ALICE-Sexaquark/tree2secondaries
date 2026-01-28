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

# apply option: optimize #

if(OPTIMIZE)
  target_compile_options(App PRIVATE -O3 -DNDEBUG -mtune=native -march=native
                                     -ffp-contract=fast)
endif(OPTIMIZE)

# apply option: debug #

if(DEBUG)
  target_compile_definitions(App PUBLIC T2S_DEBUG)
endif(DEBUG)

# apply option: legacy KFParticle #

if(LEGACY_KF)
  target_compile_definitions(App PUBLIC T2S_LEGACY_KF HomogeneousField)
endif(LEGACY_KF)
