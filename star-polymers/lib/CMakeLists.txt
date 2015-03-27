add_library(MatVec MatVec.cpp)

add_library(Particle Particle.cpp)
target_link_libraries(Particle MatVec)

add_library(Molecule Molecule.cpp)
target_link_libraries(Molecule Particle MatVec Rand)

add_library(Box Box.cpp)
target_link_libraries(Box MatVec Particle Molecule)

add_library(Rand Rand.cpp) 

add_library(Thermostat Thermostat.cpp)
target_link_libraries(Thermostat Box)

add_library(Thermostat_None Thermostat_None.cpp)
target_link_libraries(Thermostat_None Box Thermostat) 

add_library(Lowe_Andersen Lowe_Andersen.cpp)
target_link_libraries(Lowe_Andersen Box Rand Thermostat) 

