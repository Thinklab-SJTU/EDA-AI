project(RSMTests)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

file(GLOB SRCS eval.c)

include_directories("geosteiner")
include_directories("geosteiner/lp_solve_2.3")

add_library(eval SHARED ${SRCS})


target_link_libraries(eval ${CMAKE_CURRENT_SOURCE_DIR}/geosteiner/.libs/libgeosteiner.a)
target_link_libraries(eval ${CMAKE_CURRENT_SOURCE_DIR}/geosteiner/lp_solve_2.3/libLPS.a)
message(STATUS "I am fine.")
