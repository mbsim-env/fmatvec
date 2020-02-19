cmake_minimum_required (VERSION 3.10)
project(fmatvec)

include(CheckSymbolExists)

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

################################################################################
### The version number etc.
# Name of package 
SET(PACKAGE ${CMAKE_PROJECT_NAME})
SET (fmatvec_MAJOR_VERSION 2)
SET (fmatvec_MINOR_VERSION 0)
SET (fmatvec_PATCH_VERSION 0)
# Version number of package 
SET(VERSION ${fmatvec_MAJOR_VERSION}.${fmatvec_MINOR_VERSION}.${fmatvec_PATCH_VERSION})
# Define to the address where bug reports for this package should be sent. 
SET(PACKAGE_BUGREPORT "martin.o.foerg@googlemail.com")
# Define to the full name of this package. 
SET(PACKAGE_NAME ${CMAKE_PROJECT_NAME})
# Define to the full name and version of this package. 
SET(PACKAGE_STRING "${CMAKE_PROJECT_NAME} ${VERSION}")
# Define to the one symbol short name of this package. 
SET(PACKAGE_TARNAME ${CMAKE_PROJECT_NAME})
# Define to the home page for this package. 
SET(PACKAGE_URL "")
# Define to the version of this package. 
SET(PACKAGE_VERSION ${VERSION})
# Define to 1 if you have the ANSI C header files. 
SET(STDC_HEADERS 1)

#############################################
### Options #################################
#############################################

SET(INSTALL_PKGCONFIG_DIR "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

################################################################################
add_compile_options(-std=c++11)
add_compile_options(-Wno-unused-function)
add_compile_options(-Werror)

SET(HAVE_BOOST_PREPROCESSOR_ITERATION_LOCAL_HPP ${Boost_preprocessor_FOUND})
FIND_FILE(HAVE_BOOST_PREPROCESSOR_ITERATION_LOCAL_HPP "boost/preprocessor/iteration/local.hpp" DOC "Does the boost/preprocessor/iteration/local.hpp exist?")
FIND_FILE(HAVE_DLFCN_H     "dlfcn.h"     DOC "Does the dlfcn.h exist?")
FIND_FILE(HAVE_INTTYPES_H  "inttypes.h"  DOC "Does the inttypes.h exist?")
FIND_FILE(HAVE_MEMORY_H    "memory.h"    DOC "Does the memory.h exist?")
FIND_FILE(HAVE_STDBOOL_H   "stdbool.h"   DOC "Does the stdbool.h exist?")
FIND_FILE(HAVE_STDINT_H    "stdint.h"    DOC "Does the stdint.h exist?")
FIND_FILE(HAVE_STDLIB_H    "stdlib.h"    DOC "Does the stdlib.h exist?")
FIND_FILE(HAVE_STRINGS_H   "strings.h"   DOC "Does the strings.h exist?")
FIND_FILE(HAVE_STRING_H    "string.h"    DOC "Does the string.h exist?")
FIND_FILE(HAVE_SYS_STAT_H  "sys/stat.h"  DOC "Does the sys/stat.h exist?")
FIND_FILE(HAVE_SYS_TYPES_H "sys/types.h" DOC "Does the sys/types.h exist?")
FIND_FILE(HAVE_UNISTD_H    "unistd.h"    DOC "Does the unistd.h exist?")

SET(HAVE_LIBATLAS  ${ATLAS_FOUND})
SET(HAVE_LIBBLAS   ${BLAS_FOUND})
SET(HAVE_LIBLAPACK ${LAPACK_FOUND})
check_symbol_exists(sqrt "math.h" HAVE_SQRT)
check_symbol_exists(pow "math.h" HAVE_POW)


CONFIGURE_FILE("config.h.in-cmake" "config.h")

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")
################################################################################

#############################################
### Search for dependencies #################
#############################################
# Checking for the blas and lapack dependency
FIND_PACKAGE(LAPACK REQUIRED) # checks for blas and lapack!

# Check for the boost dependency
find_package(Boost REQUIRED) 

#############################################
### Include the subdirectories ##############
#############################################
add_subdirectory(fmatvec)
#add_subdirectory(doc)



################################################################################
## Config.cmake generation and installation

set(fmatvec_needed_libraries ${fmatvec_needed_libraries} ${LAPACK_LIBRARIES} )


set(ConfigPackageLocation "${CMAKE_INSTALL_LIBDIR}/cmake/fmatvec")
install(EXPORT fmatvec
     NAMESPACE fmatvec::
     DESTINATION ${ConfigPackageLocation})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fmatvecConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfig.cmake" @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfigVersion.cmake"
    VERSION ${VERSION}
    COMPATIBILITY AnyNewerVersion
    )

install(FILES 
             "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfig.cmake" 
             "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfigVersion.cmake" 
        DESTINATION ${ConfigPackageLocation})


#############################################
###  Installing #############################
#############################################
###  Install the pkg-config file
configure_file("fmatvec.pc.in" "fmatvec.pc" @ONLY)
install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/fmatvec.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
