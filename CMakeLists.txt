cmake_minimum_required(VERSION 3.15)

project(sequential_test)

add_definitions(-std=c++17)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fno-omit-frame-pointer" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_BUILD_TYPE Debug)

include_directories(sequential)
# include_directories(parallel_naive)
include_directories(tests)

add_executable(main main.cpp)

add_subdirectory(sequential)
# add_subdirectory(parallel_naive)
add_subdirectory(tests)

target_link_libraries(main PUBLIC sequential)
# target_link_libraries(main PUBLIC parallel_naive)
target_link_libraries(main PUBLIC tests)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)

