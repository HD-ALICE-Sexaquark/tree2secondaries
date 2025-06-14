add_executable(App Packager.cxx Finder.cxx App.cxx)

target_link_libraries(
  App
  PUBLIC ROOT::Core ROOT::RIO ROOT::Tree KFParticle::KFParticle
  PRIVATE CLI11::CLI11)

target_include_directories(
  App PUBLIC ${ROOT_INCLUDE_DIRS} ${KFParticle_INCLUDE_DIR}
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
