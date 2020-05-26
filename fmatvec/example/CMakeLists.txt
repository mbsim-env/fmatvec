cmake_minimum_required(VERSION 3.10)

add_executable(eigval
               "eigval.cc")

target_include_directories(eigval
   PRIVATE
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

target_link_libraries(eigval fmatvec)

