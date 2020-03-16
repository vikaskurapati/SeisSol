#  Clang-Format - C++ formating tool
#  website: https://clang.llvm.org/docs/ClangFormat.html
#
#  NOTE: it is not an official cmake search file
#  author: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de
#
#  ClangFormat_FOUND        - system has Clang-Format
#
#  Additional env. variables that may be used by this module.
#  They can change the default behaviour and
#  could to be set before calling find_package:
#
#  CLANG_FORMAT_DIR         - the root directory of the Clang-Format installation
#
#  Variables that must be set in the outer scope:
#  CLANG_MINIMUM_VERSION_REQUIRED  - a version of the required package
#       example:  set(CLANG_MINIMUM_VERSION_REQUIRED 6.0)

include(FindPackageHandleStandardArgs)

find_program(CLANG_FORMAT_EXE clang-format
             HINTS ENV CLANG_FORMAT_DIR)

if (CLANG_FORMAT_EXE)
    execute_process(COMMAND clang-format --version
                    COMMAND grep "version ${CLANG_MINIMUM_VERSION_REQUIRED}"
                    RESULT_VARIABLE TERMINATION_CODE
                    OUTPUT_QUIET)
    set(SUCCESSFUL_TERMINATION 0)
    if (${TERMINATION_CODE} STREQUAL ${SUCCESSFUL_TERMINATION})
        set(CLANG_CORRECT_VERSION True)
    else()
        message(WARNING "the executable has been found but it doesn't match"
                "the requested version: ${CLANG_MINIMUM_VERSION_REQUIRED}")
    endif()
endif()

find_package_handle_standard_args(ClangFormat CLANG_FORMAT_EXE
                                              CLANG_CORRECT_VERSION)