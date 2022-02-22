#*******************************************************************************
# @file CMakeLists.txt
# @author Florian Eigentler
# @brief
# @version 0.1
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
#*******************************************************************************
cmake_minimum_required(VERSION 3.16.3)
project(FV3D VERSION 1.0.0 DESCRIPTION "Finite volume solver")

#*******************************************************************************
# CMake settings
#******************************************************************************/
mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
mark_as_advanced(CMAKE_OSX_SYSROOT)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#*******************************************************************************
# Compiler settings
#******************************************************************************/
if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_STANDARD 99)
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -pg")
elseif(CMAKE_C_COMPILER_ID MATCHES AppleClang)
    set(CMAKE_C_STANDARD 99)
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wextra")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -pg")
else()
    MESSAGE(FATAL_ERROR
        "Provided unsupported C compiler ${CMAKE_C_COMPILER_ID}!")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS Debug Release Profile)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions("-DDEBUG")
endif()

#*******************************************************************************
# User options
#******************************************************************************/
set(build_unit_tests ON CACHE BOOL "Build unit tests")
set(memory_check ON CACHE BOOL "Build unit tests together with memory check")
set(memory_check_sup OFF CACHE BOOL "Use supression file with memory check")
set(basec_dir "../basec/install" CACHE PATH "Base C library path" )

if(build_unit_tests)
    enable_testing()
    if(memory_check)
        set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full \
            --show-leak-kinds=all --show-reachable=yes --track-origins=yes")
        if(memory_check_sup)
            set(MEMORYCHECK_SUPPRESSIONS_FILE
                "${PROJECT_BINARY_DIR}/basec.sup")
        endif()
        include(CTest)
        find_program(MEMORYCHECK_COMMAND NAMES valgrind)
    endif()
endif()

if(use_mpi)
    add_definitions("-DMPI")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install"
         CACHE PATH "default install path" FORCE)
endif()

include("${basec_dir}/share/cmake/basec.cmake")

#*******************************************************************************
# Add subdirecotries
#******************************************************************************/
add_subdirectory(src)