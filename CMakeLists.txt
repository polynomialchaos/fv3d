####################################################################################################################################
# FV3D - Finite volume solver
# (c) 2020 | Florian Eigentler
####################################################################################################################################
cmake_minimum_required( VERSION 3.15 )
enable_language( C )
project( "FV3D" )
include( CTest )

# hide following contents from ccmake panel
mark_as_advanced( CMAKE_OSX_ARCHITECTURES )
mark_as_advanced( CMAKE_OSX_DEPLOYMENT_TARGET )
mark_as_advanced( CMAKE_OSX_SYSROOT )

####################################################################################################################################
# Compiler settings
#-----------------------------------------------------------------------------------------------------------------------------------
if( CMAKE_C_COMPILER_ID MATCHES GNU )
    set( ADDITIONAL_C_COMPILE_FLAGS "-std=c99" )
    set( CMAKE_C_FLAGS_RELEASE "${ADDITIONAL_C_COMPILE_FLAGS} -O3" )
    set( CMAKE_C_FLAGS_DEBUG "${ADDITIONAL_C_COMPILE_FLAGS} -O0 -g -Wall -Wextra" )
    set( CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -pg" )
elseif( CMAKE_C_COMPILER_ID MATCHES AppleClang )
    set( ADDITIONAL_C_COMPILE_FLAGS "-std=c99" )
    set( CMAKE_C_FLAGS_RELEASE "${ADDITIONAL_C_COMPILE_FLAGS} -O3" )
    set( CMAKE_C_FLAGS_DEBUG "${ADDITIONAL_C_COMPILE_FLAGS} -O0 -g -Wall -Wextra" )
    set( CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -pg" )
else()
   MESSAGE( FATAL_ERROR "Provided unsupported C compiler ${CMAKE_C_COMPILER_ID}!" )
endif()

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release" FORCE )
    set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile )
endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
    add_definitions( "-DDEBUG" )
endif()

####################################################################################################################################
# System directories
#-----------------------------------------------------------------------------------------------------------------------------------
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_MODULE_PATH              ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )

####################################################################################################################################
# User options
#-----------------------------------------------------------------------------------------------------------------------------------
set( BUILD_TESTING ON CACHE BOOL "Build the testing tree" )
set( MEMORYCHECK ON CACHE BOOL "Build the testing tree with memory check" )
set( SHARED OFF CACHE BOOL "Provide a shared library" )
set( C_DIR "${CMAKE_SOURCE_DIR}/../00_libs/GNU-MPI/LIBC" CACHE PATH "C library path" )
set( METIS_DIR "/usr" CACHE PATH "HDF5 library path" )

if( BUILD_TESTING )
    enable_testing()
    if( MEMORYCHECK )
        find_program( MEMORYCHECK_COMMAND NAMES valgrind )
        set( MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --show-reachable=yes --track-origins=yes" )
    endif()
endif()

if( SHARED )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" )
endif()

if( CMAKE_SYSTEM_NAME MATCHES Darwin )
    set( suffix "dylib" )
else()
    set( suffix "so" )
endif()

####################################################################################################################################
# External libraries
#-----------------------------------------------------------------------------------------------------------------------------------
set( external_libraries "" )
set( external_library_compile_flags "" )
set( external_library_link_flags "" )

include( "${C_DIR}/share/cmake/C.cmake" )

####################################################################################################################################
# Add subdirecotries
#-----------------------------------------------------------------------------------------------------------------------------------
add_subdirectory( src )

####################################################################################################################################
# Install target
#-----------------------------------------------------------------------------------------------------------------------------------
install( DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h" )
install( EXPORT FV3D DESTINATION share/cmake )
