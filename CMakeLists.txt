cmake_minimum_required(VERSION 3.12)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)

# Set compiler flags
add_compile_options(-Wall -O3 -flto -g)


# Add paths to external libraries
find_package(Eigen3 REQUIRED)

# List of source files
set(SRCS 
    src/main.cpp 
    src/Orbit/Orbit.cpp 
    src/Coordinate/Coordinate.cpp 
    src/Time/Time.cpp
)

# List of header files
set(HDRS 
    src/Orbit/Orbit.h 
    src/Coordinate/Coordinate.h 
    src/Time/Time.h
)

# List of object files
add_executable(main ${SRCS} ${HDRS})

# Include directories
target_include_directories(main PRIVATE ${EIGEN3_INCLUDE_DIR})

# Set the output directory for the executable
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add a custom target for running the test case
add_custom_target(test
    COMMAND $<TARGET_FILE:main>
    COMMAND ${CMAKE_COMMAND} -E echo "Running Gnuplot..."
    COMMAND gnuplot ${CMAKE_SOURCE_DIR}/plot/3d.gnu
    COMMAND open {CMAKE_BINARY_DIR}/bin/3d.png
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running Test and Gnuplot"
)


# Add the test case for ctest
add_test(NAME test COMMAND $<TARGET_FILE:main> ${TEST_PARAMS})