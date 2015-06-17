include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")




add_library(Thermostat Thermostat.cpp)
target_link_libraries(Thermostat Box)

add_library(Thermostat_None Thermostat_None.cpp)
target_link_libraries(Thermostat_None Box Thermostat) 

add_library(Lowe_Andersen Lowe_Andersen.cpp)
target_link_libraries(Lowe_Andersen Box Rand Thermostat) 

add_library(Andersen Andersen.cpp)
target_link_libraries(Andersen Box Rand Thermostat) 

add_library(Nose_Hoover Nose_Hoover.cpp)
target_link_libraries(Nose_Hoover Box Thermostat)
