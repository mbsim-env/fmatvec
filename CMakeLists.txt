CMAKE_MINIMUM_REQUIRED (VERSION 2.8.12)

PROJECT(fmatvec)

### The version number
SET (fmatvec_VERSION_MAJOR 0)
SET (fmatvec_VERSION_MINOR 1)

#############################################
### Options #################################
#############################################

#############################################
### Search for dependencies #################
#############################################

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
