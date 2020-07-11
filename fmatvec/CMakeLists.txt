cmake_minimum_required(VERSION 3.10)

#############################################
# add the project library with the sources
set( fmatvecSrc
   _memory.cc
   ast.cc
   atom.cc
   linear_algebra_complex.cc
   linear_algebra_double.cc
   stream.cc 
   wrapper.cc 
)

file(GLOB_RECURSE fmatvecHdr "*.h")

add_library(fmatvec SHARED ${fmatvecSrc} ${fmatvecHdr})

target_include_directories( fmatvec 
   PUBLIC 
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
   PRIVATE
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

if( MSVC )
  set_source_files_properties(stream.cc PROPERTIES COMPILE_FLAGS /bigobj)
endif()

target_link_libraries(fmatvec ${LAPACK_LIBRARIES})

#############################################
### Install
#############################################

# Install the library fmatvec
if( WIN32 ) # WIN32 diffs in destination for the library (bin on Windows?!?) and must include lapack stuff for portability
  include(InstallRequiredSystemLibraries)
  install(TARGETS fmatvec EXPORT fmatvec
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION bin
    ARCHIVE DESTINATION lib
    CONFIGURATIONS ${CMAKE_CONFIGURATION_TYPES}
    )
  install(FILES ${LAPACK_LIBRARIES} DESTINATION lib)
else()
  install(TARGETS fmatvec EXPORT fmatvec
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    CONFIGURATIONS ${CMAKE_CONFIGURATION_TYPES}
    )
endif()

# Install the library-headers
install(FILES  ${fmatvecHdr}
  DESTINATION include/fmatvec)

# Make build takes that info from config.h
target_compile_definitions(fmatvec
  PUBLIC
  HAVE_LIBBLAS
)

target_compile_definitions(fmatvec
  PUBLIC
  HAVE_LIBLAPACK
)

#############################################
### Include the subdirectories
#############################################
add_subdirectory(check)
add_subdirectory(example)
