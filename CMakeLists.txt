cmake_minimum_required (VERSION 3.5)
project (ValoMC)
find_package(OpenMP)

set (CMAKE_BUILD_TYPE "Release")

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message("64 bit target platform detected.")
else()
    message(FATAL_ERROR "Did not detect a 64 bit target platform. Use \"cmake --help\" to list all available target platforms and select a 64 bit one using e.g. \n cmake . -G \"Visual Studio 15 2017 Win64\" ")
endif()

find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DUSE_OMP)
endif()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR})

add_executable(MC2D cpp/2d/MC2D.cpp)
add_executable(MC3D cpp/3d/MC3D.cpp)


if (Matlab_FOUND)
   matlab_add_mex(NAME MC2Dmex OUTPUT_NAME MC2Dmex ${CMAKE_CURRENT_SOURCE_DIR}/cpp/2d/MC2Dmex.cpp)
   matlab_add_mex(NAME MC3Dmex OUTPUT_NAME MC3Dmex ${CMAKE_CURRENT_SOURCE_DIR}/cpp/3d/MC3Dmex.cpp)
   matlab_add_mex(NAME createBH3 OUTPUT_NAME createBH3 ${CMAKE_CURRENT_SOURCE_DIR}/cpp/3d/createBH3mex.cpp)
   get_filename_component(Matlab_MX_LIBRARY_PATH ${Matlab_MX_LIBRARY} DIRECTORY)
   find_library(Matlab_UT_LIBRARY NAME libut.lib libut.so HINTS ${Matlab_MX_LIBRARY_PATH})
   if(Matlab_UT_LIBRARY)
           message("UT library found, adding support for CTRL-C ")
           target_link_libraries(MC2Dmex ${Matlab_UT_LIBRARY} )
           target_link_libraries(MC3Dmex ${Matlab_UT_LIBRARY} )
           target_link_libraries(createBH3 ${Matlab_UT_LIBRARY} )
           add_definitions(-D_MAKE_CTRL_C_POSSIBLE_)
   else()
           message("No UT library found, no support for CTRL-C ")
   endif()
endif()
