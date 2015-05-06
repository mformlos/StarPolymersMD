add_executable(star-polymers star-polymers.cpp)
target_link_libraries (star-polymers MatVec Particle Molecule Box Thermostat Thermostat_None Andersen MPC)

add_executable(mpc-test mpc_test.cpp) 
target_link_libraries(mpc-test Particle Molecule Box MPC Function_Output)



