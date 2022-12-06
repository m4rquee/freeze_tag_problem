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
include_directories("brkga/src")
include_directories($ENV{GUROBI_HOME}/include)

# Utils used in all targets
add_library(graph_utils src/mygraphlib.cpp src/geompack.cpp src/myutils.cpp
        src/mycolor.cpp src/ftp_utils.cpp)

# Lemon lib directory
link_directories(${CMAKE_SOURCE_DIR}/lemon/lemon-1.3.1/lib)

# Gurobi lib directory
link_directories($ENV{GUROBI_HOME}/lib)

# Compiler flags
set(CMAKE_CXX_FLAGS "-m64 -Ofast -Wall -Wunused-local-typedefs -lemon -fopenmp")

# Target executable and linking
add_executable(ftp_ilp src/ftp_ilp.cpp)
target_link_libraries(ftp_ilp m)
target_link_libraries(ftp_ilp pthread)
target_link_libraries(ftp_ilp graph_utils)
target_link_libraries(ftp_ilp gurobi_c++ gurobi95)
