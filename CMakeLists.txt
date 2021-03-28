cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)

set(CMAKE_BUILD_TYPE Debug)

project(sequential_test)

include_directories(sequential)

add_executable(main main.cpp)

add_subdirectory(sequential)

target_link_libraries(main sequential)
