cmake_minimum_required(VERSION 3.21)

project(Freeze_Tag_Problem_ILP
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Set a build type if none was specified
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
add_library(graph_utils
        include/mygraphlib.hpp src/mygraphlib.cpp
        include/myutils.hpp src/myutils.cpp
        include/mycolor.hpp src/mycolor.cpp
        include/problem_utils.hpp src/problem_utils.cpp)
link_libraries(graph_utils -l:libemon.a) # link the lemon binary

# Compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wunused-local-typedefs")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE("Activating -Ofast optimizations")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math")
endif ()

# Target executable and linking (Freeze-Tag Problem):
add_executable(ftp_solver src/solver_ilp.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver m pthread)
target_link_libraries(ftp_solver graph_utils)
target_link_libraries(ftp_solver gurobi_c++ gurobi110)

# Target executable and linking (Degree-Bounded Minimum Height Spanning Tree):
add_executable(bdhst_solver src/solver_ilp.cpp src/problem_utils.cpp)
set_target_properties(bdhst_solver PROPERTIES COMPILE_DEFINITIONS "BDHST")
target_link_libraries(bdhst_solver m pthread)
target_link_libraries(bdhst_solver graph_utils)
target_link_libraries(bdhst_solver gurobi_c++ gurobi110)

# Target executable and linking (Freeze-Tag Problem with Shortest Path Tree formulation):
add_executable(ftp_solver_spt src/solver_ilp_spt.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver_spt m pthread)
target_link_libraries(ftp_solver_spt graph_utils)
target_link_libraries(ftp_solver_spt gurobi_c++ gurobi110)
