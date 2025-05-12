add_executable(App Manager.cxx App.cxx)

target_link_libraries(App PUBLIC ROOT::Core ROOT::RIO ROOT::Tree)

target_include_directories(
  App
  PUBLIC ${ROOT_INCLUDE_DIRS}
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)
