cmake_minimum_required(VERSION 3.12)

project(sequential)

add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)

set(
  HEADERS
  fragment.hpp
  plane.hpp
  optimizer.hpp
  one_dim_point_tree.hpp
  point_tree.hpp
  map_plane.hpp
)
set(
  SOURCE
  fragment.cpp
  plane.cpp
  optimizer.cpp
  one_dim_point_tree.cpp
  point_tree.cpp
  map_plane.cpp
)

add_library(sequential ${SOURCE} ${HEADERS})

target_include_directories(sequential PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
