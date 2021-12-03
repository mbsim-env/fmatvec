cmake_minimum_required(VERSION 3.10)

project(fmatvec VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(MSVC)
  # https://devblogs.microsoft.com/cppblog/permissive-switch/
  add_compile_options(/permissive- /wd4996)
else()
  add_compile_options(-Wno-format-security) # this is required by types.h stating "fprintf(stderr, out_c_str);" which is marked by gcc -Wformat-security which is enabled by default on Ubuntu distributions
  add_compile_options(-Wno-unused-function)
  add_compile_options(-Werror)
  add_compile_options(-fvisibility=hidden)
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

################################################################################
### Search for dependencies ####################################################
################################################################################

# Checking for the blas and lapack dependency
if( WIN32 )
  # Path with blas and lapack lib/dll needs to be in environment variable LIB
  # find_package(BLAS) called separately for slightly better error message
  find_package(BLAS)
  if(BLAS)
    find_package(LAPACK)
    set(BLAS_LAPACK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
  endif()

  if(NOT LAPACK)
    # fall back to Intel MKL runtime libs
    message("   ... trying to find Intel MKL runtime libraries instead of LAPACK/BLAS")
    # get link library
    find_library(LAPACK_LIBRARIES name "mkl_rt.lib" PATHS $ENV{LIB} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} REQUIRED)
    # ... and runtime libraries needed to be installed
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
    set(mkl_needed_libs "mkl_rt;mkl_intel_thread")
    foreach(lib ${mkl_needed_libs})
      find_library(mkl_${lib} name "${lib}" PATHS $ENV{LIB} ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
      install(FILES ${mkl_${lib}} DESTINATION bin)
      message("    install mkl_needed_libs: mkl_${lib}= ${mkl_${lib}}")
    endforeach()
    set(BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} )
  endif()
  set(Boost_USE_STATIC_LIBS ON)
else()
  find_package(LAPACK REQUIRED) # checks for blas and lapack!
  set(BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} )
endif()
# Notes for cmake limitations under Windows:
#     cmake 3.12.18081601-MSVC_2 fails to find c:\local\boost_1_72_0 (it ONLY detects the coexisting c:\local\boost_1_60_0)
#         With boost 1.67 it works as expected
#     cmake 3.17 finds c:\local\boost_1_72_0 extracted from boost_<version>-bin-msvc-all-32-64.7z
#     Alternative: set BOOST_INCLUDEDIR to boost installation
find_package(Boost 1.67 REQUIRED COMPONENTS system) # COMPONENTS spirit uuid dll

# C++ compilation will fail due to missing config.h -> write an empty config.h for the MSVC build
set(FMATVEC_CONFIG_INCLUDE ${CMAKE_BINARY_DIR}/include)
file(WRITE ${FMATVEC_CONFIG_INCLUDE}/config.h "/* DUMMY FILE*/")

################################################################################
### Include the subdirectories #################################################
################################################################################
add_subdirectory(fmatvec)
add_subdirectory(doc)

################################################################################
### Config.cmake generation and installation
################################################################################


set(ConfigPackageLocation "${CMAKE_INSTALL_LIBDIR}/cmake/fmatvec")
install(EXPORT fmatvec
     NAMESPACE fmatvec::
     DESTINATION ${ConfigPackageLocation})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fmatvecConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfig.cmake" @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfigVersion.cmake"
    COMPATIBILITY AnyNewerVersion
    )

install(FILES 
             "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfig.cmake" 
             "${CMAKE_CURRENT_BINARY_DIR}/config/fmatvecConfigVersion.cmake" 
        DESTINATION ${ConfigPackageLocation})

# install license
install(FILES "COPYING" "NOTICE" DESTINATION ${CMAKE_INSTALL_PREFIX}/share/fmatvec)

################################################################################
###  Installing ################################################################
################################################################################
###  Install the pkg-config file
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix \${prefix})
set(libdir      \${prefix}/lib)
set(includedir  \${prefix}/include)
set(LIBDIRS "-L.")
set(LIBS "-lpthread")
set(INCDIRS "-I${Boost_INCLUDE_DIRS}")
foreach(LIB ${BLAS_LAPACK_LIBRARIES})
  string(REGEX MATCH "(.*)/lib([^\\.]*)\\.(.*)" LIBTMP "${LIB}") # "." must have double excaping "\\"
  set(LIBS "${LIBS} -L${CMAKE_MATCH_1} -l${CMAKE_MATCH_2}")
endforeach()
foreach(LIB Boost::system)
  get_property(LIBLOC TARGET ${LIB} PROPERTY LOCATION)
  string(REGEX MATCH "(.*)/lib([^\\.]*)\\.(.*)" LIBTMP "${LIBLOC}") # "." must have double excaping "\\"
  set(LIBS "${LIBS} -L${CMAKE_MATCH_1} -l${CMAKE_MATCH_2}")
endforeach()
configure_file("fmatvec.pc.in" "fmatvec.pc" @ONLY)
install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/fmatvec.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
