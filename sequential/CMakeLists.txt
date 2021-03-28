cmake_minimum_required(VERSION 2.8)

project(sequential)

set(HEADERS fragment.hpp plane.hpp optimizer.hpp)
set(SOURCE fragment.cpp plane.cpp optimizer.cpp)

add_library(sequential ${SOURCE} ${HEADERS})

target_include_directories(sequential PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
