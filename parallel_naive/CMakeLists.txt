cmake_minimum_required(VERSION 3.15)

project(parallel_naive)

add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_BUILD_TYPE Debug)

set(HEADERS parallel_naive.hpp)
set(SOURCE parallel_naive.cpp)

add_library(parallel_naive ${SOURCE} ${HEADERS})

# include_directories(sequential)
# link_directories(sequential)

# target_link_libraries(parallel_naive sequential)

target_include_directories(parallel_naive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(parallel_naive PRIVATE OpenMP::OpenMP_CXX)
