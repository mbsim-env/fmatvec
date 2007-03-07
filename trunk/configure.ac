#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(fmatvec, 0.9.0, mfoerg@users.berlios.de)
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
LIBFMATVEC_CURRENT=0
LIBFMATVEC_REVISION=0
LIBFMATVEC_AGE=0
AC_SUBST([LIBFMATVEC_CURRENT])
AC_SUBST([LIBFMATVEC_REVISION])
AC_SUBST([LIBFMATVEC_AGE])
LIBFMATVEC_VERSION_INFO=$LIBFMATVEC_CURRENT:$LIBFMATVEC_REVISION:$LIBFMATVEC_AGE
AC_SUBST([LIBFMATVEC_VERSION_INFO])

atlas_libs=yes
atlas_headers=yes
pthread_libs=yes

AC_ARG_WITH(atlas-inc-prefix,[  --with-atlas-inc-prefix=PFX  Prefix where atlas includes are installed], [atlas_inc_prefix="$withval"], [atlas_inc_prefix="."])

AC_ARG_WITH(atlas-lib-prefix,[  --with-atlas-lib-prefix=PFX  Prefix where atlas libraries are installed], [atlas_lib_prefix="$withval"], [atlas_lib_prefix="."])

AC_ARG_WITH([f2c],[  --with-f2c=PFX  the path and name to the libf2c.a library(optional)],[f2c="$withval"],[f2c=""])

CXXFLAGS="-O2"
if test -n "$F77"; then
 	AC_F77_LIBRARY_LDFLAGS
 	if test -n "$FFLAGS"; then
 		AC_MSG_RESULT([defining FFLAGS to be $FFLAGS])
 	fi
 	have_fortran_compiler=true
fi

LIBS=""

AC_CHECK_LIB(pthread, pthread_mutex_init,,pthread_libs=no,"-lpthread")
if test "x$pthread_libs" = xno; then
	AC_MSG_ERROR(pthread library needed!)
fi

dnl Add any special include directories
#AC_MSG_CHECKING(for ATLAS CPPFLAGS)
if test "$atlas_inc_prefix" != "" ; then
	CFLAGS="$CFLAGS  -I$atlas_inc_prefix"
	CPPFLAGS="$CPPFLAGS  -I$atlas_inc_prefix"
fi

AC_MSG_CHECKING(for ATLAS headers)
AC_TRY_COMPILE([#include <clapack.h>],[],[AC_MSG_RESULT(found)], [AC_MSG_ERROR(not present)])
  
#AC_MSG_CHECKING(for ATLAS LDFLAGS)
if test "$atlas_lib_prefix" != "" ; then
	LDFLAGS="$LDFLAGS -L$atlas_lib_prefix"
fi

LIBS=""

dnl add the alsa library
ATLAS_LIBS="-llapack -lcblas -lblas -latlas -lpthread $FLIBS"

dnl AC_MSG_RESULT($LDFLAGS)

AC_CHECK_LIB(atlas,cblas_dcopy,,atlas_libs=no,$ATLAS_LIBS)
if test "x$atlas_libs" = xno; then
	AC_MSG_ERROR(ATLAS libraries were not found on the system. You can get the ATLAS package from http://math-atlas.sourceforge.net.)
fi

AC_SUBST([atlas_inc_prefix])
AC_SUBST([atlas_lib_prefix])

dnl Checking for doxygen

AC_ARG_WITH([doxygenpath],[  --with-doxygenpath=PATH  The path to the 'doxygen' program.],[doxygenpath="$withval"],[doxygenpath=""])
AC_PATH_PROG([doxygen],[doxygen],[no],[$doxygenpath:$PATH])
if test $doxygen = "no" ; then
  AC_MSG_RESULT([can't generate documentation. 'doxygen' required.])
fi

LIBS="-llapack -lcblas -latlas -lblas -lpthread $FLIBS"
AC_SUBST(LIBS)
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
