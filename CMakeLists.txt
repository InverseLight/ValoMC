cmake_minimum_required (VERSION 3.5)
project (ValoMC)
find_package(OpenMP)
include(FindMatlab)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ./)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DUSE_OMP)
endif()

add_executable(MC2D cpp/2d/MC2D.cpp)
add_executable(MC3D cpp/3d/MC3D.cpp)

if (Matlab_FOUND)
   matlab_add_mex(NAME MC2Dmex OUTPUT_NAME MC2Dmex ${CMAKE_CURRENT_SOURCE_DIR}/cpp/2d/MC2Dmex.cpp)
   matlab_add_mex(NAME MC3Dmex OUTPUT_NAME MC3Dmex ${CMAKE_CURRENT_SOURCE_DIR}/cpp/3d/MC3Dmex.cpp)
   matlab_add_mex(NAME createBH3 OUTPUT_NAME createBH3 ${CMAKE_CURRENT_SOURCE_DIR}/cpp/3d/createBH3mex.cpp)
endif()
