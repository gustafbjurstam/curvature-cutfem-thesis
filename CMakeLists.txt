cmake_minimum_required(VERSION 3.16)
project(curvature_cutfem_thesis VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LDFLAGS} -fPIC -O3 -g")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MPI REQUIRED)
find_package(CUTFEM REQUIRED)

find_library(CUTFEM_PARALLEL_LIBRARY
    NAMES parallel libparallel
    PATHS ${CUTFEM_LIBRARY_DIR}
    REQUIRED
)

function(add_curvature_executable FULL_NAME)
    get_filename_component(EXE_NAME ${FULL_NAME} NAME_WLE)

    add_executable(${EXE_NAME} "src/${FULL_NAME}")

    target_link_libraries(${EXE_NAME} PUBLIC
        ${CUTFEM_LIBRARIES}
        ${CUTFEM_PARALLEL_LIBRARY}
        MPI::MPI_CXX
        ${ARGN}
    )

    target_include_directories(${EXE_NAME} PUBLIC
        ${CUTFEM_INCLUDE_DIR}
    )
endfunction()

add_curvature_executable("curvature/curvature_2d.cpp")
add_curvature_executable("curvature/curvature_2d_algoim.cpp" lapacke lapack)
add_curvature_executable("curvature/curvature_3d.cpp")
add_curvature_executable("curvature/curvature_3d_hexa.cpp")
add_curvature_executable("curvature/curvature_3d_algoim.cpp" lapacke lapack)