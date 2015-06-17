include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")


add_library(Hydrodynamics Hydrodynamics.cpp) 
target_link_libraries(Hydrodynamics Molecule Box)

add_library(MPC MPC.cpp)
target_link_libraries(MPC Particle Rand Box Hydrodynamics)
