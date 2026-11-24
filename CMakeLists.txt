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

# Keep the bundled dependency's generated files in the main build tree.
set(LEMON_ENABLE_GLPK NO CACHE STRING "Enable Lemon GLPK backend")
set(LEMON_ENABLE_ILOG NO CACHE STRING "Enable Lemon CPLEX backend")
set(LEMON_ENABLE_COIN NO CACHE STRING "Enable Lemon COIN backend")
set(LEMON_ENABLE_SOPLEX NO CACHE STRING "Enable Lemon SoPlex backend")
add_subdirectory(lemon/lemon-1.3.1 "${CMAKE_BINARY_DIR}/lemon" EXCLUDE_FROM_ALL)

# Gurobi lib directory
link_directories($ENV{GUROBI_HOME}/lib)
find_library(GUROBI_LIBRARY NAMES gurobi130 gurobi120 gurobi110 gurobi100 gurobi95
        HINTS "$ENV{GUROBI_HOME}/lib" REQUIRED)

# Utils used in all targets
add_library(graph_utils
        include/mygraphlib.hpp src/mygraphlib.cpp
        include/myutils.hpp src/myutils.cpp
        include/mycolor.hpp src/mycolor.cpp
        include/problem_utils.hpp src/problem_utils.cpp)
target_link_libraries(graph_utils PUBLIC lemon)
link_libraries(graph_utils)

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
target_link_libraries(ftp_solver gurobi_c++ ${GUROBI_LIBRARY})

# Target executable and linking (Degree-Bounded Minimum Height Spanning Tree):
add_executable(bdhst_solver src/solver_ilp.cpp src/problem_utils.cpp)
set_target_properties(bdhst_solver PROPERTIES COMPILE_DEFINITIONS "BDHST")
target_link_libraries(bdhst_solver m pthread)
target_link_libraries(bdhst_solver graph_utils)
target_link_libraries(bdhst_solver gurobi_c++ ${GUROBI_LIBRARY})

# Target executable and linking (Freeze-Tag Problem with Shortest Path Tree formulation):
add_executable(ftp_solver_spt src/solver_ilp_spt.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver_spt m pthread)
target_link_libraries(ftp_solver_spt graph_utils)
target_link_libraries(ftp_solver_spt gurobi_c++ ${GUROBI_LIBRARY})

# Target executable and linking (Freeze-Tag Problem with Flow formulation):
add_executable(ftp_solver_flow src/solver_ilp_flow.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver_flow m pthread)
target_link_libraries(ftp_solver_flow graph_utils)
target_link_libraries(ftp_solver_flow gurobi_c++ ${GUROBI_LIBRARY})

# Target executable and linking (Freeze-Tag Problem with Normal formulation + leaves restriction):
add_executable(ftp_solver_leaf src/solver_ilp_leaf.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver_leaf m pthread)
target_link_libraries(ftp_solver_leaf graph_utils)
target_link_libraries(ftp_solver_leaf gurobi_c++ ${GUROBI_LIBRARY})

# Target executable and linking (Freeze-Tag Problem with Normal formulation + Lagrangian relaxation):
add_executable(ftp_solver_lagran src/solver_ilp_lagran.cpp src/problem_utils.cpp)
target_link_libraries(ftp_solver_lagran m pthread)
target_link_libraries(ftp_solver_lagran graph_utils)
target_link_libraries(ftp_solver_lagran gurobi_c++ ${GUROBI_LIBRARY})

include(CTest)
if(BUILD_TESTING)
    add_executable(ftp_pdf_viewer_test tests/pdf_viewer.cpp)
    target_link_libraries(ftp_pdf_viewer_test graph_utils)
    add_test(NAME pdf_viewer COMMAND ftp_pdf_viewer_test)
endif()
