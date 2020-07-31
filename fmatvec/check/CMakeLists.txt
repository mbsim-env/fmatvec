cmake_minimum_required(VERSION 3.10)

set(VALGRINDARGS "")
if( WIN32 )
  find_program(DIFF NAMES fc diff)
  set(VALGRIND "")
else()
  find_program(DIFF NAMES diff)
  find_program(VALGRIND NAMES valgrind)
  if(VALGRIND)
    set(VALGRINDARGS --gen-suppressions=all --suppressions=/home/zander-ste/Development/src/mbsim/fmatvec/fmatvec/check/valgrindtestast.supp --error-exitcode=200 --num-callers=150 --leak-check=full --show-reachable=yes)
  else()
    set(VALGRIND "")
  endif()
endif()

################################################################################
# generate test executables
################################################################################
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

################################################################################
# define tests, each taking steps run and diff as individual parts
################################################################################
add_custom_target(testsymfunction_run
    COMMAND ${VALGRIND} ${VALGRINDARGS} $<TARGET_FILE_DIR:testsymfunction>/$<TARGET_FILE_NAME:testsymfunction>  > testsymfunction.out # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    DEPENDS testsymfunction
    COMMENT "Run testsymfunction"
)
add_custom_target(testsymfunction_diff
    COMMAND ${DIFF} ${CMAKE_CURRENT_SOURCE_DIR}/testsymfunction.ref testsymfunction.out
    DEPENDS testsymfunction_run
    COMMENT "Diff reulsts of testsymfunction"
)
add_dependencies(testsymfunction_diff testsymfunction_run)

add_custom_target(testast_run
    COMMAND ${VALGRIND} ${VALGRINDARGS} $<TARGET_FILE_DIR:testast>/$<TARGET_FILE_NAME:testast> > testast.out # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    DEPENDS testast
    COMMENT "Run testast"
)
add_custom_target(testast_diff
    COMMAND ${DIFF} ${CMAKE_CURRENT_SOURCE_DIR}/testast.ref testast.out
    COMMENT "Diff reulsts of testast"
)
add_dependencies(testast_diff testast_run)

add_custom_target(testfunction_run
    COMMAND ${VALGRIND} ${VALGRINDARGS} $<TARGET_FILE_DIR:testfunction>/$<TARGET_FILE_NAME:testfunction> # add_custom_command can take target names and expands to regular platform-specific paths/executable names
    COMMENT "Run testfunction"
    DEPENDS testfunction
)

add_custom_target(check
    DEPENDS testfunction_run testsymfunction_diff testast_diff
    COMMENT "Run test binaries"
)

