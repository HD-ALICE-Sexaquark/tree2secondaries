add_executable(App Packager.cxx Finder.cxx App.cxx)

target_link_libraries(
  App
  PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist KFParticle::KFParticle
  PRIVATE fmt::fmt CLI11::CLI11)

target_include_directories(
  App PUBLIC ${ROOT_INCLUDE_DIRS} ${KFParticle_INCLUDE_DIR}
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
