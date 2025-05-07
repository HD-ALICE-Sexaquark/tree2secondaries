add_executable(App App.cxx Parser.cxx)

target_link_libraries(App PUBLIC ROOT::EG ROOT::Tree ROOT::RIO ROOT::Core
                                 Analysis Math)

target_include_directories(
  App PUBLIC ${ROOT_INCLUDE_DIRS}
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>)
