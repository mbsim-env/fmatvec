cmake_minimum_required(VERSION 3.10)

#############################################
# add the project library with the sources
set( fmatvecSrc
   _memory.cc
   ast.cc
   atom.cc
   linear_algebra_complex.cc
   linear_algebra_double.cc
   sparse_linear_algebra_double.cc
   stream.cc 
   stream100.cc 
   stream200.cc 
   stream300.cc 
   wrapper.cc 
   spooles.cc
)

set(fmatvecHdr
   var_fixed_general_matrix.h
   ast.h
   fixed_var_general_matrix.h
   _memory.h
   function.h
   symbolic_function1_temp.h
   range.h
   vector.h
   var_square_matrix.h
   band_matrix.h
   wrapper.h
   atom.h
   stream.h
   square_matrix.h
   row_vector.h
   matrix.h
   fixed_square_matrix.h
   linear_algebra_complex.h
   diagonal_matrix.h
   stream_impl.h
   types.h
   var_row_vector.h
   symbolic_function.h
   fixed_general_matrix.h
   toString.h
   fixed_row_vector.h
   linear_algebra.h
   var_general_matrix.h
   var_vector.h
   fmatvec.h
   general_matrix.h
   symbolic_function2_temp.h
   sparse_matrix.h
   symmetric_sparse_matrix.h
   fixed_symmetric_matrix.h
   symmetric_matrix.h
   var_symmetric_matrix.h
   linear_algebra_double.h
   symbolic.h
   fixed_vector.h
   indices.h
   spooles.h
   sparse_linear_algebra_double.h
)

add_library(fmatvec SHARED ${fmatvecSrc} ${fmatvecHdr})

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(fmatvec
    PUBLIC
    FMATVEC_DEBUG
  )
endif()

target_include_directories( fmatvec 
   PUBLIC 
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
   $<BUILD_INTERFACE:${ARPACK_INCLUDE_DIRS}>
   $<BUILD_INTERFACE:${SPOOLES_INCLUDE_DIRS}>
   PRIVATE
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

if( MSVC )
  set_source_files_properties(
     stream.cc
     stream100.cc
     stream200.cc
     stream300.cc
     PROPERTIES COMPILE_FLAGS /bigobj
  )
endif()

target_link_libraries(fmatvec ${BLAS_LAPACK_LIBRARIES} Boost::system ${ARPACK_LIBRARIES} ${SPOOLES_LIBRARIES})
if(NOT WIN32)
  target_link_libraries(fmatvec pthread)
else()
  target_link_libraries(fmatvec bcrypt)
endif()

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
