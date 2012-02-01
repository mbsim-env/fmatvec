SUBDIRS =  .
  
lib_LTLIBRARIES = libfmatvecTestbench.la
libfmatvecTestbench_la_LDFLAGS = -version-info @LIBFMATVECTESTBENCH_VERSION_INFO@
libfmatvecTestbench_la_SOURCES = testbench.cc\
                                 randomgenerator.cc\
                                 matgenerator.cc

libfmatvecTestbench_la_LIBADD = $(DEPS_LIBS)
INCLUDES = $(DEPS_CFLAGS)

fmatvecTestbenchincludedir = $(includedir)/fmatvecTestbench
fmatvecTestbenchinclude_HEADERS = testbench.h\
                                  randomgenerator.h\
                                  matgenerator.h
