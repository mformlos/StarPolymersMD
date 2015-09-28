include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")


include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")

include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")




add_executable(star-polymers star-polymers.cpp)
target_link_libraries (star-polymers MatVec Particle Molecule Box Thermostat Hydrodynamics Thermostat_None Andersen MPC stdc++)

#add_executable(continue_star_polymers continue_star_polymers.cpp)
#target_link_libraries (continue_star_polymers MatVec Particle Molecule Box Thermostat Hydrodynamics Thermostat_None Andersen MPC stdc++)

add_executable(mpc-test mpc_test.cpp) 
target_link_libraries(mpc-test Particle Molecule Box MPC Function_Output stdc++)

#add_executable(shear shear.cpp)
#target_link_libraries(shear Particle Molecule Box Thermostat Hydrodynamics Thermostat_None Andersen MPC stdc++)
