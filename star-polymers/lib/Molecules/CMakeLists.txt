include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")


include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")

include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")



add_library(Particle Particle.cpp)
target_link_libraries(Particle MatVec)

add_library(Molecule Molecule.cpp)
target_link_libraries(Molecule Particle MatVec Rand)