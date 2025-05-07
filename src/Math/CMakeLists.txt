add_library(Math STATIC Common.cxx Vertexer.cxx)

target_link_libraries(Math PUBLIC ROOT::Core)

target_include_directories(
  Math PUBLIC ${ROOT_INCLUDE_DIRS}
              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>)
