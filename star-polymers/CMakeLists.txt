cmake_minimum_required(VERSION 2.6)
project(STAR_POLYMERS)

list( APPEND CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS} -Wall -g -fopenmp -O3 ")

 
include_directories ("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_subdirectory (lib)

include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")
add_subdirectory (lib/Analysis)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
add_subdirectory (lib/Hydrodynamics)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")
add_subdirectory (lib/Thermostats)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")
add_subdirectory (lib/Molecules)

