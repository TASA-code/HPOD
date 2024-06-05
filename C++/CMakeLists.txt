cmake_minimum_required(VERSION 3.12)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)

# Set compiler flags
add_compile_options(-Wall -O3 -flto -g)

# Add include directory
include_directories(include)  # This line adds -Iinclude

# Add paths to external libraries
find_package(Eigen3 REQUIRED)
find_package(NLOPT REQUIRED)

# List of source files
set(SRCS 
    src/main.cpp 
    src/propagator.cpp
    src/Coordinate/coordinate.cpp 
    src/Time/time.cpp
    src/accel/accel.cpp
    src/utils/rk45.cpp
    src/utils/input.cpp
)

# List of header files
set(HDRS 
    include/propagator.h
    include/coordinate/coordinate.h 
    include/time/time.h
    include/accel/accel.h
    include/utils/rk45.h
    include/utils/input.h
)

# List of object files
add_executable(main ${SRCS} ${HDRS})

# Include directories
target_include_directories(main PRIVATE ${EIGEN3_INCLUDE_DIR})

# Link the NLOPT library
target_link_libraries(main PRIVATE ${NLOPT_LIBRARIES})

# Set the output directory for the executable
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add a custom target for running the test case
add_custom_target( sim
    COMMAND $<TARGET_FILE:main>
    COMMAND ${CMAKE_COMMAND} -E echo "Running Gnuplot..."
    COMMAND gnuplot ${CMAKE_SOURCE_DIR}/plot/3d.gnu
    COMMAND open ${CMAKE_SOURCE_DIR}/plot/3d.png
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running Simulation and Gnuplot"
)


# Add the test case for ctest
add_test(NAME test COMMAND $<TARGET_FILE:main> ${TEST_PARAMS})