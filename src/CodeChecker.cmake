function(add_prefix PREFIX LIST)
    foreach(ITEM ${LIST})
        list(APPEND RESULT "${PREFIX}${ITEM}")
    endforeach()
    set(RETURN_VALUE ${RESULT} PARENT_SCOPE)
endfunction()


# include all source files to get checked
set(FILES_TO_CHECK ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/LocalIntegrationRecorder.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/NeighbIntegrationRecorder.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/PlasticityRecorder.cpp)

# handle all include directories. Note, cannot use target_include_directories b/c a custom target doesn't build anything
set(INCLUDE_DIRS src
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer
                 ${CMAKE_CURRENT_SOURCE_DIR}/src
                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/${EQUATIONS}
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/async
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/easi/include
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen3
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/glm
                 ${CMAKE_CURRENT_SOURCE_DIR}/submodules/Device)

add_prefix("-I" "${INCLUDE_DIRS}")
set(INCLUDE_DIRS ${RETURN_VALUE})

# handle all include directories. See the reason above
set(DEFINITION_LIST "${HARDWARE_DEFINITIONS}"
                    CONVERGENCE_ORDER=${ORDER}
                    NUMBER_OF_QUANTITIES=${NUMBER_OF_QUANTITIES}
                    NUMBER_OF_RELAXATION_MECHANISMS=${NUMBER_OF_MECHANISMS})
add_prefix("-D" "${DEFINITION_LIST}")
set(DEFINITION_LIST ${RETURN_VALUE})

set(CLANG_MINIMUM_VERSION_REQUIRED 6.0)
if (FILES_TO_CHECK)
    #find_package(ClangFormat REQUIRED)
    #if (ClangFormat_FOUND)
    #    add_custom_target(formating ALL COMMAND clang-format --style=file ${FILES_TO_CHECK} > ./format-log)
    #endif()


    find_package(ClangTidy REQUIRED)
    if (ClangTidy_FOUND)
        add_custom_target(static_analyze ALL COMMAND clang-tidy -extra-arg="-std=c++11" ${FILES_TO_CHECK} --config=""
                          --
                          ${INCLUDE_DIRS} ${DEFINITION_LIST}
                          DEPENDS
                            src/generated_code/tensor.h
                            src/generated_code/subroutine.h
                            src/generated_code/init.h
                            src/generated_code/kernel.h)

        # add target-dependencies
        if (ClangFormat_FOUND)
            add_dependencies(static_analyze formating)
        endif()
        add_dependencies(SeisSol-lib static_analyze)
    endif()
endif()