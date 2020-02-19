#############################################
### Include the subdirectories
#############################################
# add the project library with the sources
file(GLOB src *.cc)

    
add_library(fmatvec SHARED ${src})
target_link_libraries(fmatvec ${LAPACK_LIBRARIES})

#############################################
### Install
#############################################

# Install the library fmatvec
install(TARGETS fmatvec EXPORT fmatvec
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)

# Install the library-headers
install(DIRECTORY ./
  DESTINATION include/fmatvec
  FILES_MATCHING PATTERN "*.h")

