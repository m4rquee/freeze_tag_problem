cmake_minimum_required(VERSION 3.21)

project(Freeze_Tag_Problem_ILP
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Set a default build type if none was specified
set(default_build_type "Release")

# Headers directories
include_directories("include")
include_directories("lemon/lemon-1.3.1/include")
include_directories($ENV{GUROBI_HOME}/include)

# Lemon lib directory
link_directories(${CMAKE_SOURCE_DIR}/lemon/lemon-1.3.1/build/lemon)

# Gurobi lib directory
link_directories($ENV{GUROBI_HOME}/lib)

# Utils used in all targets
add_library(graph_utils src/mygraphlib.cpp src/geompack.cpp src/myutils.cpp
        src/mycolor.cpp src/problem_utils.cpp)
link_libraries(graph_utils -l:libemon.a) # link the lemon binary

# Compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wunused-local-typedefs")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE("Activating -Ofast optimizations")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
endif ()

# Target executable and linking
add_executable(solver_ilp src/solver_ilp.cpp)
target_link_libraries(solver_ilp m pthread)
target_link_libraries(solver_ilp graph_utils)
target_link_libraries(solver_ilp gurobi_c++ gurobi100)

#[[add_executable(ftp_qntd src/solver_ilp_ftp_qntd.cpp)
target_link_libraries(ftp_qntd m pthread)
target_link_libraries(ftp_qntd graph_utils)
target_link_libraries(ftp_qntd gurobi_c++ gurobi100)]]
