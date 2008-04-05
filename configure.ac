#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(fmatvec, 0.9.3, mfoerg@users.berlios.de)
AC_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/fmatvec.h])
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_PROG_CXXCPP
AC_LANG([C++])

dnl versioning (see http://www.gnu.org/software/libtool/manual.html)
dnl current: 	The most recent interface number that this library implements. 
dnl revision: 	The implementation number of the current interface. 
dnl age: 	The difference between the newest and oldest interfaces that
dnl		this library implements. In other words, the library implements all the
dnl 		interface numbers in the range from number current (age<=current).
LIBFMATVEC_CURRENT=1
LIBFMATVEC_REVISION=0
LIBFMATVEC_AGE=0
AC_SUBST([LIBFMATVEC_CURRENT])
AC_SUBST([LIBFMATVEC_REVISION])
AC_SUBST([LIBFMATVEC_AGE])
LIBFMATVEC_VERSION_INFO=$LIBFMATVEC_CURRENT:$LIBFMATVEC_REVISION:$LIBFMATVEC_AGE
AC_SUBST([LIBFMATVEC_VERSION_INFO])

dnl By autoscan
AC_CHECK_FUNCS([sqrt])
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

atlas_libs=yes
atlas_headers=yes

AC_ARG_WITH(atlas-inc-prefix,[  --with-atlas-inc-prefix=PFX  Prefix where atlas includes are installed], [atlas_inc_prefix="$withval"], [atlas_inc_prefix="."])

AC_ARG_WITH(atlas-lib-prefix,[  --with-atlas-lib-prefix=PFX  Prefix where atlas libraries are installed], [atlas_lib_prefix="$withval"], [atlas_lib_prefix="."])

AC_ARG_WITH(lapack-lib-prefix,[  --with-lapack-lib-prefix=PFX  Prefix where lapack library is installed], [lapack_lib_prefix="$withval"], [lapack_lib_prefix="."])

AC_ARG_WITH(blas-lib-prefix,[  --with-blas-lib-prefix=PFX  Prefix where blas library is installed], [blas_lib_prefix="$withval"], [blas_lib_prefix="."])

AC_ARG_WITH([f2c],[  --with-f2c=PFX  the path and name to the libf2c.a library(optional)],[f2c="$withval"],[f2c=""])

CXXFLAGS="-O2"
if test -n "$F77"; then
 	AC_FC_LIBRARY_LDFLAGS
 	if test -n "$FFLAGS"; then
 		AC_MSG_RESULT([defining FFLAGS to be $FFLAGS])
 	fi
 	have_fortran_compiler=true
fi

LIBS=""

AC_CHECK_LIB(pthread, pthread_mutex_init,PTHREADLIB="-lpthread",PTHREADLIB="","-lpthread")

dnl ----------------------------------------------------------

AC_ARG_ENABLE([atlas],[  --enable-atlas  Use ATLAS],[atlas="yes"],[atlas="no"])

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
ATLAS_LIBS="-llapack -lcblas -lblas -latlas $PTHREADLIB $FCLIBS"

AC_CHECK_LIB(atlas,cblas_dcopy,,atlas_libs=no,$ATLAS_LIBS)
if test "x$atlas_libs" = xno; then
	AC_MSG_ERROR(ATLAS libraries were not found on the system. You can get the ATLAS package from http://math-atlas.sourceforge.net.)
fi

LIBS="-llapack -lcblas -latlas -lblas $PTHREADLIB $FCLIBS"

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
AM_CONDITIONAL([COND_ATLAS],[test "$atlas" = "yes"])

AC_ARG_WITH([doxygenpath],[  --with-doxygenpath=PATH  The path to the 'doxygen' program.],[doxygenpath="$withval"],[doxygenpath=""])
AC_PATH_PROG([doxygen],[doxygen],[no],[$doxygenpath:$PATH])
if test $doxygen = "no" ; then
  AC_MSG_RESULT([can't generate documentation. 'doxygen' required.])
fi

AC_SUBST(LIBS)
AC_SUBST([INCDIRS])
AC_SUBST([LIBDIRS])
AC_CONFIG_FILES([Makefile src/Makefile doc/doxyfile doc/Makefile fmatvec.pc])
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
