cmake_minimum_required(VERSION 3.10)

# read config from "configure.sc"
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/configure.ac CONF_AC)
string(REGEX MATCH "AC_INIT\\(\\ *([a-zA-Z0-9]+)\\ *,\\ *([0-9\.]*)\\ *,(.*)\\)" _ ${CONF_AC}) # regex matching needs double escaping, e.g. "\\(", for special characters
  set(PROJECT_NAME "${CMAKE_MATCH_1}")
  set(VERSION      "${CMAKE_MATCH_2}")
# set(PROJECT_MAINTAINER "${CMAKE_MATCH_3}")
project(${PROJECT_NAME} VERSION ${VERSION})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(MSVC)
  # https://devblogs.microsoft.com/cppblog/permissive-switch/
  add_compile_options(/permissive- /wd4996)
else()
# string(APPEND CMAKE_CXX_FLAGS_INIT "-Wno-format-security")
  add_compile_options(-Wno-format-security) # this is required by types.h stating "fprintf(stderr, out_c_str);" which is marked by gcc -Wformat-security which is enabled by default on Ubuntu distributions
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

################################################################################
### Search for dependencies ####################################################
################################################################################

# cmake 3.12.18081601-MSVC_2 fails to find c:\local\boost_1_72_0 (it ONLY detects the coexisting c:\local\boost_1_60_0)
#     With boost 1.67 it works as expected
# cmake 3.17 finds c:\local\boost_1_72_0 extracted from boost_<version>-bin-msvc-all-32-64.7z
# Alternative: set BOOST_INCLUDEDIR to boost installation
if( WIN32 )
  set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.67 REQUIRED ) # COMPONENTS spirit uuid ... => all header only ???

# Checking for the blas and lapack dependency
IF( WIN32 )
  # Path with blas and lapack lib/dll needs to be in environment variable LIB
  # find_package(BLAS) called separately for slightly better error message
  find_package(BLAS)
  if(BLAS)
    find_package(LAPACK)
  endif()

  if(NOT LAPACK)
    MESSAGE("   ... trying to find Intel MKL runtime libraries instead of LAPACK/BLAS")
    # fall back to Intel MKL runtime libs

    # get link library
    find_library(LAPACK_LIBRARIES name "mkl_rt.lib" PATHS $ENV{LIB} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} REQUIRED)
    # ... and runtime libraries needed to be installed
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
    set(mkl_needed_libs "mkl_rt;mkl_intel_thread")
    foreach(lib ${mkl_needed_libs})
      find_library(mkl_${lib} name "${lib}" PATHS $ENV{LIB} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
      install(FILES ${mkl_${lib}} DESTINATION bin)
      MESSAGE("    install mkl_needed_libs: mkl_${lib}= ${mkl_${lib}}")
    endforeach()
  endif()
ELSE()
  FIND_PACKAGE(LAPACK REQUIRED) # checks for blas and lapack!
ENDIF()

################################################################################
set(FMATVEC_CONFIG_INCLUDE ${CMAKE_BINARY_DIR}/include)
file(WRITE ${FMATVEC_CONFIG_INCLUDE}/config.h "/* DUMMY FILE*/")


################################################################################
### Include the subdirectories #################################################
################################################################################
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


################################################################################
###  Installing ################################################################
################################################################################
###  Install the pkg-config file
set(prefix ${CMAKE_INSTALL_PREFIX})
set(LIBDIRS "-L.")
set(LIBS "-lpthread")
foreach(LIB ${fmatvec_needed_libraries})
  string(REGEX MATCH "(.*)/lib(.*)\\.(.*)" LIBTMP "${LIB}") # "." must have double excaping "\\"
  set(LIBS "${LIBS} -L${CMAKE_MATCH_1} -l${CMAKE_MATCH_2}")
endforeach()
configure_file("fmatvec.pc.in" "fmatvec.pc" @ONLY)
install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/fmatvec.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
