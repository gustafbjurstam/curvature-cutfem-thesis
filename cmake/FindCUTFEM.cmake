set(CUTFEM_FOUND "NO")

set(CUTFEM_ROOT "" CACHE PATH "Path to the CutFEM-Library source directory")
set(CUTFEM_BUILD_DIR "" CACHE PATH "Path to the CutFEM-Library build directory")

if(NOT CUTFEM_ROOT)
    message(FATAL_ERROR "Set CUTFEM_ROOT, e.g. -DCUTFEM_ROOT=/path/to/CutFEM-Library")
endif()

if(NOT CUTFEM_BUILD_DIR)
    set(CUTFEM_BUILD_DIR "${CUTFEM_ROOT}/build")
endif()

set(CUTFEM_COMPONENT
    libcommon.so
    libFESpace.so
    libsolver.so
    libproblem.so
)

find_path(CUTFEM_INCLUDE_DIR
    NAMES cutfem.hpp
    PATHS "${CUTFEM_ROOT}/cpp"
    REQUIRED
)

find_path(CUTFEM_LIBRARY_DIR
    NAMES libcommon.so
    PATHS "${CUTFEM_BUILD_DIR}/lib"
    REQUIRED
)

foreach(lib_name ${CUTFEM_COMPONENT})
    find_library(COMP_LIBRARY
        NAMES ${lib_name}
        PATHS ${CUTFEM_LIBRARY_DIR}
        NO_DEFAULT_PATH
        REQUIRED
    )
    list(APPEND CUTFEM_LIBRARIES ${COMP_LIBRARY})
    unset(COMP_LIBRARY CACHE)
endforeach()

set(CUTFEM_FOUND "YES")

message(STATUS "CutFEM include directory: ${CUTFEM_INCLUDE_DIR}")
message(STATUS "CutFEM library directory: ${CUTFEM_LIBRARY_DIR}")