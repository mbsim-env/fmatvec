cmake_minimum_required(VERSION 3.10)

add_executable(testfunction
               "testfunction.cc")

target_include_directories(testfunction
   PRIVATE
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

target_link_libraries(testfunction fmatvec)

add_executable(testast
               "testast.cc")

target_include_directories(testast
   PRIVATE
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

target_link_libraries(testast fmatvec)

add_executable(testsymfunction
               "testsymfunction.cc")

target_include_directories(testsymfunction
   PRIVATE
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
   $<BUILD_INTERFACE:${FMATVEC_CONFIG_INCLUDE}>
)

target_link_libraries(testsymfunction fmatvec)

add_custom_command(OUTPUT testsymfunction.out
    COMMAND testsymfunction > testsymfunction.out # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    DEPENDS testsymfunction
)

add_custom_command(OUTPUT testast.out
    COMMAND testast > testast.out # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    DEPENDS testast
)

add_custom_target(check
    testfunction # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    COMMAND ${CMAKE_COMMAND} -E compare_files "${CMAKE_CURRENT_SOURCE_DIR}/testsymfunction.ref" testsymfunction.out
    COMMAND ${CMAKE_COMMAND} -E compare_files "${CMAKE_CURRENT_SOURCE_DIR}/testast.ref" testast.out
    DEPENDS testfunction testast testsymfunction.out testast.out
    COMMENT "Run test binaries"
)
