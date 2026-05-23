# cmake_minimum_required(VERSION 3.16)
# project(FORCE_AI VERSION 0.1 LANGUAGES CXX)

# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_MODULE_PATH "/Users/sebastianmyrback/Documents/doktorandprojekt/cutfem_library/cmake/")

# message(CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH})

# find_package(MPI REQUIRED)
# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
# include_directories(${MPI_CXX_INCLUDE_PATH})

# find_package(CUTFEM REQUIRED)

# set(name_files  
#     fictitious
#     bulk_dg
#     bulk_dg_algoim
#     eriks_maxwell
#     eriks_maxwell2
#     maxwell3D)

# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LDFLAGS}  -std=c++2a -fPIC -O3")

# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

# foreach(name_file ${name_files})
#     add_executable(${name_file} "src/${name_file}.cpp")
#     target_link_libraries("${name_file}" PUBLIC ${CUTFEM_LIBRARIES}  )
#     target_include_directories("${name_file}" PUBLIC ${CMAKE_BINARY_DIR} 
#                                                      ${CUTFEM_INCLUDE_DIR})
# endforeach(name_file)

cmake_minimum_required(VERSION 3.16)
project(LIBCUTFEM VERSION 0.1 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the module path
set(CMAKE_MODULE_PATH "/mnt/hot/User/Work/repositories/SF250X/")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

# Find MPI
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
include_directories(${MPI_CXX_INCLUDE_PATH})

# Find CUTFEM
find_package(CUTFEM REQUIRED)
# Gustaf
find_library(CUTFEM_PARALLEL_LIBRARY
    NAMES parallel libparallel
    PATHS
        /home/gbjurstam/Work/repositories/SF250X/CutFEM-Library/build/lib
    REQUIRED
)

message(STATUS "CUTFEM_PARALLEL_LIBRARY: ${CUTFEM_PARALLEL_LIBRARY}")

# Ensure CXX_FLAGS are set correctly
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LDFLAGS} -fopenmp -fPIC -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LDFLAGS} -fPIC -O3 -g")

# Set output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# List of source files
# set(name_files  
#     fictitious
#     bulk_dg
#     bulk_dg_algoim
#     eriks_maxwell2
#     maxwell3D
#     3field
#     3field_old
#     kikuchi
#     wave)

# Include Boost if needed for Bessel functions
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


find_path(GMSH_INCLUDE_DIR NAMES gmsh.h PATHS /opt/homebrew/include REQUIRED)
find_library(GMSH_LIBRARY   NAMES gmsh libgmsh PATHS /opt/homebrew/lib REQUIRED)

add_library(gmsh::gmsh SHARED IMPORTED)
set_target_properties(gmsh::gmsh PROPERTIES
  IMPORTED_LOCATION "${GMSH_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
)


function(add_example FULL_NAME)
    get_filename_component(EXE_NAME ${FULL_NAME} NAME_WLE )
    add_executable(${EXE_NAME} "src/${FULL_NAME}")
    
    # ${ARGN} will expand to any extra arguments passed to the function
    target_link_libraries("${EXE_NAME}" PUBLIC 
        ${CUTFEM_LIBRARIES} 
        ${CUTFEM_PARALLEL_LIBRARY} # Gustaf
        MPI::MPI_CXX # Gustaf
        ${Boost_LIBRARIES} 
        gmsh::gmsh 
        ${ARGN}
    )
    
    target_include_directories("${EXE_NAME}" PUBLIC 
        ${CMAKE_BINARY_DIR} 
        ${CUTFEM_INCLUDE_DIR}
    )
endfunction()

add_example(other/poissonFEM.cpp)

add_example(convection_diffusion/merging_circles.cpp)
add_example(convection_diffusion/bulk_dg.cpp)
add_example(convection_diffusion/bulk_dg_algoim.cpp)
add_example(convection_diffusion/langmuir.cpp)

add_example(transport/high_order_transport.cpp)
add_example(transport/transport.cpp)

add_example(stokes/fictitious.cpp)
add_example(stokes/fictitious_nitsche.cpp)
add_example(stokes/fictitious_lagr.cpp)
add_example(stokes/fictitious_lagr_nit.cpp)
add_example(stokes/fictitious_vel_lagr_nit.cpp)
add_example(stokes/fictitious_double_lagr.cpp)
add_example(stokes/fictitious_split_lagr.cpp)
add_example(stokes/normal_lagr_nit.cpp)
add_example(stokes/lagr_fictitious.cpp)
add_example(stokes/split_lagr_fictitious.cpp)
add_example(stokes/lagr_nit_fictitious.cpp)
add_example(stokes/block_fictitious.cpp)
add_example(stokes/block_fitted.cpp)
add_example(stokes/test_block.cpp)
add_example(stokes/stokes_fitted.cpp)
add_example(stokes/no_flow.cpp)
add_example(stokes/cartesian_nitsche.cpp)
add_example(stokes/scott_vogelius_fitted.cpp)
add_example(stokes/dock_intf.cpp)
add_example(stokes/interface.cpp)
add_example(stokes/interface_lagr.cpp)
add_example(stokes/interface_lagr_nit.cpp)
add_example(stokes/interface_split_lagr.cpp)
add_example(stokes/time_nitsche.cpp)
add_example(stokes/time_lagr_nit.cpp)

add_example(darcy/darcy.cpp)

add_example(maxwell/3field.cpp)
add_example(maxwell/3field_old.cpp)
add_example(maxwell/3field_new.cpp)
add_example(maxwell/3field_maxwell.cpp)
add_example(maxwell/3field_23april25.cpp)
add_example(maxwell/3field_7may25.cpp)
add_example(maxwell/3field_9may25.cpp)
add_example(maxwell/kikuchi.cpp)
add_example(maxwell/wave.cpp)
add_example(maxwell/wave_new.cpp)
add_example(maxwell/kikuchi_new.cpp)

add_example(navier_stokes/taylor_green_vortex.cpp)
add_example(navier_stokes/nitsche.cpp)
add_example(navier_stokes/scott_vogelius.cpp)

add_example(other/test_barycentric.cpp)
add_example(other/test_BDM2_seb.cpp)

add_example(active_surfaces/surface_problem.cpp)
add_example(active_surfaces/test.cpp)
add_example(active_surfaces/generate_mesh.cpp)


# add_example("curvature/test_curvature.cpp")
add_example("curvature/curvature_2d_algoim.cpp" lapacke lapack)
add_example("curvature/curvature_2d.cpp")
add_example("curvature/curvature_3d.cpp")
add_example("curvature/curvature_3d_hexa.cpp")
