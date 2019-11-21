
string(REPLACE ":" ";" default_bin_search_paths $ENV{PATH})
string(REPLACE ":" ";" default_lib_search_paths $ENV{LIBRARY_PATH})
string(REPLACE ":" ";" default_includes_search_paths $ENV{C_INCLUDE_PATH} ":" $ENV{CPLUS_INCLUDE_PATH} ":" $ENV{CPATH})


string(REPLACE "," ";" component_list ${GemmTools_FIND_COMPONENTS})
set(GemmTools_FOUND True)
foreach(component ${component_list})

    if ("${component}" STREQUAL "LIBXSMM")
        set(component_found True)


    elseif ("${component}" STREQUAL "PSpaMM")
        set(component_found True)


    elseif ("${component}" STREQUAL "MKL")
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
        set(GemmTools_DEFINITIONS ${GemmTools_DEFINITIONS} -DMKL_ILP64)

        if (DEFINED ENV{MKLROOT})
            set(GemmTools_INCLUDE_DIR $ENV{MKLROOT}/include)
        endif()

        set(component_found True)


    elseif ("${component}" STREQUAL "OpenBLAS")
        # check wether we can find the openblas lib within user's environment
        find_file(OpenBLAS_lib
                  NAMES libopenblas.so
                  HINTS ${default_lib_search_paths})


        # check wether we can find an openblas header within user's environment
        find_file(OpenBLAS_include
                NAMES cblas.h
                HINTS ${default_includes_search_paths})

        if (OpenBLAS_lib STREQUAL OpenBLAS_lib-NOTFOUND  OR  OpenBLAS_include EQUAL OpenBLAS_include-NOTFOUND)
            message(STATUS "OpenBLAS INFO:")
            message(STATUS "\tlib: ${OpenBLAS_lib}")
            message(STATUS "\tincl.: ${OpenBLAS_include}")
            set(component_found False)
        else()
            set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} -lopenblas)
            set(component_found True)
        endif()


    elseif ("${component}" STREQUAL "BLIS")
        set(component_found True)


    elseif ("${component}" STREQUAL "ACL_DEVICE_BLAS")
        add_subdirectory(submodules/Device/cuda)
        set(GemmTools_LIBRARIES ${GemmTools_LIBRARIES} custom_blas)
        set(GemmTools_INCLUDE_DIR ${GemmTools_INCLUDE_DIR} submodules/Device/cuda/src)
        set(GemmTools_DEFINITIONS ${GemmTools_DEFINITIONS} ACL_DEVICE)

        set(component_found True)


    else()
        message(FATAL_ERROR "Gemm Tools do not have a requested component, i.e. ${component}. \
                Please, see the documentation")
    endif()


    # GemmTools is not found if at least one requested component is missing
    if (NOT component_found)
        set(GemmTools_FOUND False)
    endif()

endforeach()


if (NOT GemmTools_FOUND)
    message(STATUS "\tNOTE: clean cmake-cache after changing your enviroment variables")
endif()