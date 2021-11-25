#*******************************************************************************
# @file unit_test.cmake
# @author Florian Eigentler
# @brief
# @version 0.1
# @date 2021-11-08
# @copyright Copyright (c) 2021
#******************************************************************************/
function(unit_test target ut_target_source)
    get_filename_component(test_name ${ut_target_source} NAME_WE)
    set(test_target unit_test_${target}_${test_name})

    add_executable(${test_target} ${ut_target_source})
    add_dependencies(${test_target} ${target})
    target_link_libraries(${test_target} ${CMAKE_PROJECT_NAME})

    set_target_properties(${test_target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

    add_test(
        NAME ${test_target} COMMAND ${test_target}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
    )
endfunction()