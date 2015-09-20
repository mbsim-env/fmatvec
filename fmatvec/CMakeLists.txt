#############################################
### Include the subdirectories
#############################################

# Checking for the blas dependency
FIND_PACKAGE(BLAS REQUIRED)
SET(HAVE_LIBBLAS 1)

# Check for the boost dependency
find_package(Boost REQUIRED) 

# add the project library with the sources
add_library(fmatvec SHARED 
atom.cc
CMakeLists.txt
linear_algebra_complex.cc
linear_algebra_double.cc
_memory.cc
wrapper.cc
)
    
target_link_libraries(fmatvec)

# Create a config.h due to the checks made before
configure_file ("${PROJECT_SOURCE_DIR}/fmatvec/config.h.in"
                "${PROJECT_BINARY_DIR}/config.h" )
# Include the path when compiling
include_directories(${PROJECT_BINARY_DIR})

#############################################
### Install
#############################################

# Install the library fmatvec
install(TARGETS fmatvec
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)

# Install the library-headers
install(DIRECTORY ./
  DESTINATION include/fmatvec
  FILES_MATCHING PATTERN "*.h")

