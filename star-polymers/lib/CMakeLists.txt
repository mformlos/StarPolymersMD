include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")


include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")

include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")

include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")





add_library(MatVec MatVec.cpp)



add_library(Box Box.cpp)
target_link_libraries(Box MatVec Particle Molecule)

add_library(Rand Rand.cpp) 



