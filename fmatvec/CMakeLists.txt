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

target_link_libraries(fmatvec ${LAPACK_LIBRARIES})

#############################################
### Install
#############################################

# Install the library fmatvec
if( WIN32 ) # WIN32 diffs in destination for the library (bin on Windows?!?) and must include lapack stuff for portability
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

if( WIN32 )
  target_compile_definitions(fmatvec
      PRIVATE
      "MSVC_DLL_DECLSPEC=__declspec(dllexport)"
  )
endif()


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
