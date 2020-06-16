#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(fmatvec, 2.0.0, martin.o.foerg@googlemail.com)
AC_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([fmatvec/fmatvec.h])
LT_INIT([win32-dll])
AC_PROG_CXX
AC_PROG_CXXCPP
AC_LANG([C++])

# enable C++17
CXXFLAGS="$CXXFLAGS -std=c++17 -fvisibility=hidden"

dnl By autoscan
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([pow])
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

atlas_libs=yes
atlas_headers=yes

intelmkl_libs=yes
intelmkl_headers=yes

AC_ARG_WITH(atlas-inc-prefix,[  --with-atlas-inc-prefix=PFX  Prefix where atlas includes are installed], [atlas_inc_prefix="$withval"], [atlas_inc_prefix="."])

AC_ARG_WITH(atlas-lib-prefix,[  --with-atlas-lib-prefix=PFX  Prefix where atlas libraries are installed], [atlas_lib_prefix="$withval"], [atlas_lib_prefix="."])

AC_ARG_WITH(intelmkl-inc-prefix,[  --with-intelmkl-inc-prefix=PFX  Prefix where IntelMKL includes are installed], [intelmkl_inc_prefix="$withval"], [intelmkl_inc_prefix="."])

AC_ARG_WITH(intelmkl-lib-prefix,[  --with-intelmkl-lib-prefix=PFX  Prefix where IntelMKL libraries are installed], [intelmkl_lib_prefix="$withval"], [intelmkl_lib_prefix="."])

AC_ARG_WITH(lapack-lib-prefix,[  --with-lapack-lib-prefix=PFX  Prefix where lapack library is installed], [lapack_lib_prefix="$withval"], [lapack_lib_prefix="."])

AC_ARG_WITH(blas-lib-prefix,[  --with-blas-lib-prefix=PFX  Prefix where blas library is installed], [blas_lib_prefix="$withval"], [blas_lib_prefix="."])

AC_ARG_WITH([f2c],[  --with-f2c=PFX  the path and name to the libf2c.a library(optional)],[f2c="$withval"],[f2c=""])

dnl Define the memory allocator
AC_ARG_WITH([allocator-header],[  --with-allocator-header=<header file>  Include <header file> defining the STL-style allocator class to use (Default: memory)],[allocator_header="$withval"],[allocator_header="NONE"])
AC_ARG_WITH([allocator-class],[  --with-allocator-class=<class>  Use <class> for STL-style allocator class name. (fmatvec::PoolAllocator (default), ...)],[allocator_class="$withval"],[allocator_class="DEFAULT"])
if test "_$allocator_header" != "_NONE"; then
  CXXFLAGS="$CXXFLAGS -DALLOCATORHEADER=$allocator_header"
  DEFINES="$DEFINES -DALLOCATORHEADER=$allocator_header"
fi
if test "_$allocator_class" != "_DEFAULT"; then
  CXXFLAGS="$CXXFLAGS -DALLOCATORCLASS=$allocator_class"
  DEFINES="$DEFINES -DALLOCATORCLASS=$allocator_class"
fi
AC_SUBST([DEFINES])

dnl Include path for Boost
AC_ARG_WITH([boost-inc],
  [  --with-boost-inc=PATH  The path Boost include directory.],
  [CPPFLAGS="$CPPFLAGS -I$withval"; DEFINES="$DEFINES -I$withval"])

AC_CHECK_HEADERS([boost/preprocessor/iteration/local.hpp],[],[AC_MSG_ERROR([Can not find boost headers (boost/preprocessor/iteration/local.hpp)])])

dnl Define the Memory class
AC_ARG_WITH([memory-class],[  --with-memory-class=<class>  Use the defined memory management: (INTRUSIVECOUNT (default); STDSHAREDPTR; FMATVEC); STDSHAREDPTR uses str::shared_ptr<AT[]> using new operator; INTRUSIVECOUNT uses a fmatvec builtin intrusive reference counting using new operator],[memory_class="$withval"],[memory_class="INTRUSIVECOUNT"])
if test "_$memory_class" == "_STDSHAREDPTR"; then
  CXXFLAGS="$CXXFLAGS -DMEMORYCLASS_STDSHAREDPTR"
  DEFINES="$DEFINES -DMEMORYCLASS_STDSHAREDPTR"
elif test "_$memory_class" == "_FMATVEC"; then
  CXXFLAGS="$CXXFLAGS -DMEMORYCLASS_FMATVEC"
  DEFINES="$DEFINES -DMEMORYCLASS_FMATVEC"
fi

AC_PROG_FC
if test -n "$FC"; then
 	AC_FC_LIBRARY_LDFLAGS
 	if test -n "$FFLAGS"; then
 		AC_MSG_RESULT([defining FFLAGS to be $FFLAGS])
 	fi
 	have_fortran_compiler=true
fi

LIBS=""

AC_ARG_ENABLE([pthread],[  --disable-pthread   Disable linking against pthread],[USEPTHREAD=$enableval],[USEPTHREAD=yes])
if test "_$USEPTHREAD" = "_yes"; then
  AC_CHECK_LIB(pthread, pthread_mutex_init,PTHREADLIB="-lpthread",PTHREADLIB="","-lpthread")
fi

dnl ----------------------------------------------------------

AC_ARG_ENABLE([atlas],[  --enable-atlas  Use ATLAS],[atlas="yes"],[atlas="no"])
AC_ARG_ENABLE([parallelatlas],[  --enable-parallelatlas  Use parallel ATLAS],[parallelatlas="yes"],[parallelatlas="no"])
AC_ARG_ENABLE([intelmkl],[  --enable-intelmkl  Use IntelMKL],[intelmkl="yes"],[intelmkl="no"])
AC_ARG_ENABLE([parallelIntelMKL],[  --enable-parallelIntelMKL  Use parallel IntelMKL],[parallelIntelMKL="yes"],[parallelIntelMKL="no"])

if test "$atlas" = "yes"; then
  AC_DEFINE([COSIM],[1],[Use Cosimulation Interface])

if test "$atlas_inc_prefix" != "" ; then
	CFLAGS="$CFLAGS  -I$atlas_inc_prefix"
	CPPFLAGS="$CPPFLAGS  -I$atlas_inc_prefix"
fi

INCDIRS="-I$atlas_inc_prefix"

AC_MSG_CHECKING(for ATLAS headers)
AC_TRY_COMPILE([#include <clapack.h>],[],[AC_MSG_RESULT(found)], [AC_MSG_ERROR(not present)])
  
if test "$atlas_lib_prefix" != "" ; then
	LDFLAGS="$LDFLAGS -L$atlas_lib_prefix"
fi
LIBDIRS="-L$atlas_lib_prefix"

LIBS=""

if test "x$parallelatlas" = xyes; then
  ATLAS_LIBS="-llapack -lptcblas -latlas $PTHREADLIB $FCLIBS"
else
  ATLAS_LIBS="-llapack -lcblas -latlas $PTHREADLIB $FCLIBS"
fi

AC_CHECK_LIB(atlas,cblas_dcopy,,atlas_libs=no,$ATLAS_LIBS)
if test "x$atlas_libs" = xno; then
	AC_MSG_ERROR(ATLAS libraries were not found on the system. You can get the ATLAS package from http://math-atlas.sourceforge.net.)
fi

if test "x$parallelatlas" = xyes; then
  LIBS="-llapack -lptcblas -latlas -lblas $PTHREADLIB $FCLIBS"
else
  LIBS="-llapack -lcblas -latlas -lblas $PTHREADLIB $FCLIBS"
fi

else
  if test "$intelmkl" = "yes"; then
    if test "$intelmkl_inc_prefix" != "" ; then
	    CFLAGS="$CFLAGS  -I$intelmkl_inc_prefix"
	    CPPFLAGS="$CPPFLAGS  -I$intelmkl_inc_prefix"
    fi

    INCDIRS="-I$atlas_inc_prefix"
    
    AC_MSG_CHECKING(for IntelMKL headers)
    AC_TRY_COMPILE([#include <mkl_lapack.h>],[],[AC_MSG_RESULT(found)], [AC_MSG_ERROR(not present)])
  
    if test "$intelmkl_lib_prefix" != "" ; then
    	LDFLAGS="$LDFLAGS -L$intelmkl_lib_prefix"
    fi
    LIBDIRS="-L$intelmkl_lib_prefix"
    
    LIBS=""
    
    if test "x$parallelIntelMKL" = xyes; then
      AC_CHECK_LIB([mkl_intel_lp64],DftiComputeForward,, AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0]),[-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread])
      LIBS="-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp $LIBS"
    else
      AC_CHECK_LIB([mkl_intel_lp64],DftiComputeForward,, AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0]),[-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread])
      LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgomp $LIBS"
    fi
    
    LIBS="$LIBS $PTHREADLIB $FCLIBS"

  else
  	if test "$blas_lib_prefix" != "" ; then
  		LDFLAGS="$LDFLAGS -L$blas_lib_prefix"
  	fi
  	if test "$lapack_lib_prefix" != "" ; then
  		LDFLAGS="$LDFLAGS -L$lapack_lib_prefix"
  	fi
  	INCDIRS=""
  	LIBDIRS="-L$blas_lib_prefix -L$lapack_lib_prefix"
  	
  	LIBS=""
  	LAPACK_LIBS="-llapack -lblas $PTHREADLIB $FCLIBS"
  	
  	AC_CHECK_LIB(blas,dcopy_,,blas_libs=no,$LAPACK_LIBS)
  	if test "x$blas_libs" = xno; then
  		AC_MSG_ERROR(LAPACK library were not found on the system. You can get the
  				LAPACK package from http://.)
  	fi
  	AC_CHECK_LIB(lapack,dgemm_,,lapack_libs=no,$LAPACK_LIBS)
  	if test "x$lapack_libs" = xno; then
  		AC_MSG_ERROR(LAPACK library were not found on the system. You can get the
  				LAPACK package from http://.)
  	fi
  	
  	LIBS="-llapack -lblas $PTHREADLIB $FCLIBS"
  fi
fi
AM_CONDITIONAL([COND_ATLAS],[test "$atlas" = "yes"])

CPPFLAGS="$CPPFLAGS -Wall -Werror -Wno-unknown-pragmas"

AC_ARG_WITH([doxygenpath],[  --with-doxygenpath=PATH  The path to the 'doxygen' program.],[doxygenpath="$withval"],[doxygenpath=""])
AC_PATH_PROG([doxygen],[doxygen],[no],[$doxygenpath:$PATH])
if test $doxygen = "no" ; then
  AC_MSG_RESULT(["can't generate documentation. 'doxygen' required."])
fi

# thread-safe?
AC_ARG_ENABLE([threadsafety],[  --enable-threadsafety[=TYPE]   enables thread safety for memory managment: Posible TYPE: pthread; gccbuiltin; openmp (default)],[threadsafety="$enableval" ],[threadsafety=no])
if test $threadsafety = pthread; then
  DEFINES="$DEFINES -DFMATVEC_THREADSAFE_PTHREAD"
  CXXFLAGS="$CXXFLAGS -DFMATVEC_THREADSAFE_PTHREAD"
  LDFLAGS="$LDFLAGS -lpthread"
  LIBS="$LIBS -lpthread"
  AC_MSG_RESULT([Use pthread for thread safety.])
elif test $threadsafety = gccbuiltin; then
  DEFINES="$DEFINES -DFMATVEC_THREADSAFE_GCCBUILTIN"
  CXXFLAGS="$CXXFLAGS -DFMATVEC_THREADSAFE_GCCBUILTIN"
  AC_MSG_RESULT([Use gcc build in atomic operations for thread safety.])
elif test $threadsafety = yes -o $threadsafety = openmp; then
  DEFINES="$DEFINES -DFMATVEC_THREADSAFE_OPENMP"
  CXXFLAGS="$CXXFLAGS -fopenmp -DFMATVEC_THREADSAFE_OPENMP"
  LDFLAGS="$LDFLAGS -fopenmp -lgomp" # remove -lgomp if libtool is fixed
  LIBS="$LIBS -fopenmp -lgomp" # remove -lgomp if libtool is fixed
  AC_MSG_RESULT([Use OpenMP for thread safety.])
elif test $threadsafety = no; then
  AC_MSG_RESULT([Thread safety is disabled.])
else
  AC_MSG_ERROR([Unknown thread safety type specified by --enable-threadsafety.])
fi

# from now on only libtool is used for linking
LDFLAGS="$LDFLAGS $LDFLAGS_LIBTOOL -no-undefined -Wl,--no-undefined"
test "_$host_os" != "_mingw32" && LDFLAGS="$LDFLAGS -Wl,--disable-new-dtags,-rpath,\\\$\$ORIGIN/../lib"

if test "$cross_compiling" == "yes"; then
  AC_SUBST([XC_EXEC_PREFIX],[wine])
fi

dnl Check for MinGW
AC_CANONICAL_HOST
if test "_$host_os" = "_mingw32" ; then
  WIN32=yes
fi
AM_CONDITIONAL([COND_WIN32], [test "_$WIN32" = "_yes"])

# check for debug build
DEBUGBUILD=1
echo $CFLAGS $CXXFLAGS $CPPFLAGS | grep "\-DNDEBUG" > /dev/null && DEBUGBUILD=0
AC_SUBST([DEBUGBUILD])

AC_SUBST(LIBS)
AC_SUBST([INCDIRS])
AC_SUBST([LIBDIRS])
AC_CONFIG_FILES([Makefile fmatvec/Makefile fmatvec/check/Makefile doc/checkHtml/Makefile doc/doxy_template/Makefile doc/doxy_template/doxyfile_common doc/doxyfile doc/Makefile fmatvec.pc])
AC_CONFIG_COMMANDS([doc/doxy_template/doxyfile_common_fmatvec],
  [sed -re "s|$prefix/share/fmatvec/doxy_template/|$abs_srcdir/doc/doxy_template/|g" doc/doxy_template/doxyfile_common > doc/doxy_template/doxyfile_common_fmatvec],
  [prefix=$prefix; abs_srcdir=$(realpath $srcdir)]
  )
AC_CONFIG_FILES([fmatvec/check/testfunction.sh], [chmod +x fmatvec/check/testfunction.sh])
AC_CONFIG_FILES([fmatvec/check/testast.sh], [chmod +x fmatvec/check/testast.sh])
AC_CONFIG_FILES([fmatvec/check/testsymfunction.sh], [chmod +x fmatvec/check/testsymfunction.sh])
AC_CONFIG_FILES([fmatvec/check/valgrindtestast.sh], [chmod +x fmatvec/check/valgrindtestast.sh])
AC_CONFIG_FILES([fmatvec/check/valgrindtestfunction.sh], [chmod +x fmatvec/check/valgrindtestfunction.sh])
hardcode_into_libs=no # do not add hardcoded libdirs to ltlibraries
hardcode_libdir_flag_spec_CXX= # do not add hardcodeed libdirs to ltbinaries
AC_OUTPUT

echo "
Configuration for $PACKAGE $VERSION :
--------------------------------

 Source code location ....: $srcdir
 Host System Type ........: $host
 Preprocessor ............: $CC $CPPFLAGS
 Compiler ................: $CC $CXXFLAGS $LIBFLAC_CFLAGS
 Linker ..................: $CC $LDFLAGS $LIBS $LIBFLAC_LIBS
 Install path ............: $prefix

 Now type 'make' to build $PACKAGE $VERSION,
 and then 'make install' for installation.
"
