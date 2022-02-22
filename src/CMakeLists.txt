#*******************************************************************************
# @file CMakeLists.txt
# @author Florian Eigentler
# @brief
# @version 0.1
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
#*******************************************************************************
set(target ${CMAKE_PROJECT_NAME})
set(install_target ${target})

#*******************************************************************************
# Add subdirecotries
#******************************************************************************/
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(Solver)
add_subdirectory(FV3D)

#*******************************************************************************
# Install target
#******************************************************************************/
install(EXPORT ${install_target} DESTINATION share/cmake)