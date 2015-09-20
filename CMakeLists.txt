CMAKE_MINIMUM_REQUIRED (VERSION 2.8.12)

PROJECT(fmatvec)

### The version number
SET (fmatvec_MAJOR_VERSION 2)
SET (fmatvec_MINOR_VERSION 0)
SET (fmatvec_PATCH_VERSION 0)
set(fmatvec_VERSION ${fmatvec_MAJOR_VERSION}.${fmatvec_MINOR_VERSION}.${fmatvec_PATCH_VERSION})

#############################################
### Options #################################
#############################################

#############################################
### Search for dependencies #################
#############################################
# Checking for the blas and lapack dependency
FIND_PACKAGE(LAPACK REQUIRED) # checks for blas and lapack!

# Check for the boost dependency
find_package(Boost REQUIRED) 

# Add all libs together in one
SET(LIBS "-lblas -llapack")


#############################################
### Include the subdirectories ##############
#############################################
add_subdirectory(fmatvec)

#############################################
###  Installing #############################
#############################################
###  Install the pkg-config file
configure_file("fmatvec.pc.in" "fmatvec.pc" @ONLY)
install(FILES       "${CMAKE_CURRENT_BINARY_DIR}/fmatvec.pc"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
