include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Thermostats")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Molecules")


add_library(Function_Output Function_Output.cpp)
add_library(Streamlines Streamlines.cpp)
