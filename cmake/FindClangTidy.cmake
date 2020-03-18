#  Clang-Tidy - C++ “linter” tool
#  website: https://clang.llvm.org/extra/clang-tidy/
#
#  NOTE: it is not an official cmake search file
#  author: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de 
#
#  ClangTidy_FOUND        - system has Clang-Tidy
#
#  Additional env. variables that may be used by this module. 
#  They can change the default behaviour and
#  could to be set before calling find_package:
#
#  CLANG_TIDY_DIR         - the root directory of the Clang-Tidy installation
#
#  Variables that must be set in the outer scope:
#  CLANG_MINIMUM_VERSION_REQUIRED  - a version of the required package
#       example:  set(CLANG_MINIMUM_VERSION_REQUIRED 6.0)

include(FindPackageHandleStandardArgs)

find_program(CLANG_TIDY_EXE clang-tidy
             HINTS ${CLANG_TIDY_DIR})

if (CLANG_TIDY_EXE)
    execute_process(COMMAND clang-tidy --version
                    COMMAND grep "version ${CLANG_MINIMUM_VERSION_REQUIRED}"
                    RESULT_VARIABLE TERMINATION_CODE
                    OUTPUT_QUIET)
    set(SUCCESSFUL_TERMINATION 0)
    if (${TERMINATION_CODE} STREQUAL ${SUCCESSFUL_TERMINATION})
        set(CLANG_CORRECT_VERSION True)
        set(CLANG_TIDY_DIR ${CLANG_TIDY_EXE} CACHE STRING "path to the correct clang-tidy")
    else()
        message(WARNING "the executable has been found but it doesn't match"
                "the requested version: ${CLANG_MINIMUM_VERSION_REQUIRED}")
    endif()
endif()

find_package_handle_standard_args(ClangTidy CLANG_TIDY_EXE
                                  CLANG_CORRECT_VERSION)