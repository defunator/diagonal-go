cmake_minimum_required(VERSION 3.12)

project(tests)

add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)

set(CMAKE_BUILD_TYPE Debug)

set(HEADERS test.hpp tests.hpp)
set(SOURCE test.cpp tests.cpp)

add_library(tests ${SOURCE} ${HEADERS})

target_include_directories(tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

