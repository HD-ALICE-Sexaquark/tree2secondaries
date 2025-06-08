add_executable(App Packager.cxx Finder.cxx App.cxx)

target_link_libraries(App PUBLIC ROOT::Core ROOT::RIO ROOT::Tree
                                 KFParticle::KFParticle)

target_include_directories(
  App
  PUBLIC ${ROOT_INCLUDE_DIRS} ${KFParticle_INCLUDE_DIR}
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
