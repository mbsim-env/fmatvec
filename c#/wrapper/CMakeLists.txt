cmake_minimum_required(VERSION 3.10)

add_library(wrapper SHARED "wrapper.cc" "wrapper.h" )

target_include_directories(wrapper
   PRIVATE
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

target_link_libraries(wrapper fmatvec)

