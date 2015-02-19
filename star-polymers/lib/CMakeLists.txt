add_library(MatVec MatVec.cpp)

add_library(Particle Particle.cpp)
target_link_libraries(Particle MatVec)

add_library(Molecule Molecule.cpp)
target_link_libraries(Molecule Particle MatVec)

add_library(Box Box.cpp)
target_link_libraries(Box MatVec Particle Molecule)
