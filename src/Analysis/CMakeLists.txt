add_library(Analysis STATIC Manager.cxx)

target_link_libraries(Analysis PUBLIC ROOT::Core ROOT::Tree)

target_include_directories(
  Analysis PUBLIC ${ROOT_INCLUDE_DIRS}
                  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>)
