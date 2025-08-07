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
  add_compile_options(-Wnon-virtual-dtor)
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# clang compile DB
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_BINARY_DIR}/compile_commands.json ${CMAKE_SOURCE_DIR}/compile_commands.json)

################################################################################
### Search for dependencies ####################################################
################################################################################

# Checking for the blas and lapack dependency
if( WIN32 )
  set(Boost_USE_STATIC_LIBS ON)
endif()

if( USE_STATIC_INTEL_MKL )
  # use Intel MKL on the windows side
  message("Enforce usage of static Intel MKL instead of arbitrary LAPACK/BLAS")
  # works with IntelOne integration providing Intel(R) Math Kernel Library 2020 Update 1
  set(BLA_STATIC ON)
  set(BLA_VENDOR Intel10_64lp_seq)
endif()
find_package(LAPACK REQUIRED) # checks for blas and lapack!
set(BLAS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} )

# Notes for cmake limitations under Windows:
#     cmake 3.12.18081601-MSVC_2 fails to find c:\local\boost_1_72_0 (it ONLY detects the coexisting c:\local\boost_1_60_0)
#         With boost 1.67 it works as expected
#     cmake 3.17 finds c:\local\boost_1_72_0 extracted from boost_<version>-bin-msvc-all-32-64.7z
#     Alternative: set BOOST_INCLUDEDIR to boost installation
find_package(Boost 1.67 REQUIRED COMPONENTS system) # COMPONENTS spirit uuid dll

if(NOT DEFINED ARPACK_INCLUDE_DIRS OR NOT DEFINED ARPACK_LIBRARIES)
  find_package(arpack-ng QUIET)
endif()
if(DEFINED ARPACK_INCLUDE_DIRS OR DEFINED ARPACK_LIBRARIES)
  add_definitions(-DHAVE_ARPACK)
  message(INFO " arpack library found.")
else()
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(ARPACK arpack)
    if(ARPACK_FOUND)
      add_definitions(-DHAVE_ARPACK)
      message(INFO " arpack library found.")
    else()
      message(WARNING "arpack library not found: large sparse eigenvalue solve is done using memory/CPU expensive dense routines.")
    endif()
  else()
    message(WARNING "arpack library not found: large sparse eigenvalue solve is done using memory/CPU expensive dense routines.")
  endif()
endif()
if(NOT DEFINED SPOOLES_INCLUDE_DIRS OR NOT DEFINED SPOOLES_LIBRARIES)
  find_package(spooles QUIET)
endif()
if(DEFINED SPOOLES_INCLUDE_DIRS OR DEFINED SPOOLES_LIBRARIES)
  add_definitions(-DHAVE_SPOOLES)
  message(INFO " spooles library found.")
else()
  message(WARNING "spooles library not found: large sparse linear system solve is done using memory/CPU expensive dense routines.")
endif()

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
set(DEFINES "")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEFINES "${DEFINES} -DFMATVEC_DEBUG")
endif()
set(LIBDIRS "-L.")
set(LIBS "")
if(NOT WIN32)
  set(LIBS "${LIBS} -lpthread")
else()
  set(LIBS "${LIBS} -lbcrypt")
endif()
set(INCDIRS "-I${Boost_INCLUDE_DIRS}")
foreach(LIB ${BLAS_LAPACK_LIBRARIES})
  string(REGEX MATCH "(.*)/lib([^\\.]*)\\.(.*)" LIBTMP "${LIB}") # "." must have double excaping "\\"
  if(NOT "${CMAKE_MATCH_1}" STREQUAL "" AND NOT "${CMAKE_MATCH_2}" STREQUAL "")
    set(LIBS "${LIBS} -L${CMAKE_MATCH_1} -l${CMAKE_MATCH_2}")
  endif()
endforeach()
foreach(LIB Boost::system)
  get_property(LIBLOC TARGET ${LIB} PROPERTY LOCATION)
  string(REGEX MATCH "(.*)/lib([^\\.]*)\\.(.*)" LIBTMP "${LIBLOC}") # "." must have double excaping "\\"
  if(NOT "${CMAKE_MATCH_1}" STREQUAL "" AND NOT "${CMAKE_MATCH_2}" STREQUAL "")
    set(LIBS "${LIBS} -L${CMAKE_MATCH_1} -l${CMAKE_MATCH_2}")
  endif()
endforeach()
configure_file("fmatvec.pc.in" "fmatvec.pc" @ONLY)
install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/fmatvec.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
