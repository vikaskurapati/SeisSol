set(FILES_TO_CHECK ${FILES_TO_CHECK}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/LocalIntegrationRecorder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/NeighbIntegrationRecorder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording/PlasticityRecorder.cpp)

set(CHECK_INCLUDE_DIRS {CHECK_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/recording
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer)

set(CLANG_MINIMUM_VERSION_REQUIRED 6.0)
if (FILES_TO_CHECK)
    #find_package(ClangFormat REQUIRED)
    #if (ClangFormat_FOUND)
    #    add_custom_target(formating ALL COMMAND clang-format --style=file ${FILES_TO_CHECK} > ./format-log)
    #endif()


    find_package(ClangTidy REQUIRED)
    if (ClangTidy_FOUND)
        add_custom_target(static_analyze ALL COMMAND clang-tidy ${FILES_TO_CHECK} --config="" -- -I${CHECK_INCLUDE_DIRS})

        # add target-dependencies
        if (ClangFormat_FOUND)
            add_dependencies(static_analyze formating)
        endif()
        add_dependencies(SeisSol-lib static_analyze)
    endif()
endif()