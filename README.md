# Dashboard

Current Daily Build Status of the MBSim-Environment

| Build Type | Variant | Failed |
|------------|---------|--------|
| linux64-dailydebug | build | [![image](https://www.mbsim-env.de/service/builds/current/linux64-dailydebug/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/builds/current/linux64-dailydebug/nrAll.svg)](https://www.mbsim-env.de/builds/run/current/linux64-dailydebug/master/master/master/master/) |
| linux64-dailydebug | examples | [![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug/nrAll.svg)](https://www.mbsim-env.de/runexamples/run/current/linux64-dailydebug/master/master/master/master/) |
| linux64-dailydebug | coverage | [![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug/coverageRate.svg)](https://www.mbsim-env.de/runexamples/run/current/linux64-dailydebug/master/master/master/master/#coverage) |
| linux64-dailydebug | examples-valgrind | [![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug-valgrind/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug-valgrind/nrAll.svg)](https://www.mbsim-env.de/runexamples/run/current/linux64-dailydebug-valgrind/master/master/master/master/) |
| linux64-dailydebug | coverage-valgrind | [![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailydebug-valgrind/coverageRate.svg)](https://www.mbsim-env.de/runexamples/run/current/linux64-dailydebug-valgrind/master/master/master/master/#coverage) |
| linux64-dailyrelease | build | [![image](https://www.mbsim-env.de/service/builds/current/linux64-dailyrelease/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/builds/current/linux64-dailyrelease/nrAll.svg)](https://www.mbsim-env.de/builds/run/current/linux64-dailyrelease/master/master/master/master/) |
| linux64-dailyrelease | examples | [![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailyrelease/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/runexamples/current/linux64-dailyrelease/nrAll.svg)](https://www.mbsim-env.de/runexamples/run/current/linux64-dailyrelease/master/master/master/master/) |
| msys2win64-dailyrelease | build | [![image](https://www.mbsim-env.de/service/builds/current/msys2win64-dailyrelease/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/builds/current/msys2win64-dailyrelease/nrAll.svg)](https://www.mbsim-env.de/builds/run/current/msys2win64-dailyrelease/master/master/master/master/) |
| msys2win64-dailyrelease | examples | [![image](https://www.mbsim-env.de/service/runexamples/current/msys2win64-dailyrelease/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/runexamples/current/msys2win64-dailyrelease/nrAll.svg)](https://www.mbsim-env.de/runexamples/run/current/msys2win64-dailyrelease/master/master/master/master/) |
| - | manuals | [![image](https://www.mbsim-env.de/service/manuals/nrFailed.svg) / ![image](https://www.mbsim-env.de/service/manuals/nrAll.svg)](https://www.mbsim-env.de/service/home/#manuals) |



# fmatvec (Fast Matrix-Vector Library)

The purpose of this project is to provide a C++-library for high performance
matrix-vector computations. The software has an object oriented design allowing
the use of matrices and vectors in a comfortable way.  One of the main tasks
ist storing and managing of data types. For most of the computations ATLAS
(Automatically Tuned Linear Algebra Software) is used providing C and Fortran77
interfaces to a portably efficient BLAS implementation, as well as a few
routines from LAPACK. For a full linear algebra package the original LAPACK
library is required as well.

The library supports vectors, general matrices, symmetric matrices, banded
matrices as well as sparse matrices. Computations with dimensionless matrices
and vectors are possible. It makes use of the template mechanism in order to
provide matrices and vectors of any type like int, double, complex, etc..

# Requirements

- ATLAS from http://math-atlas.sourceforge.net (BSD license)
- LAPACK from http://www.netlib.org/lapack (BSD license)

# Installation

1. Install ATLAS with full LAPACK support, if necessary. For more
   details see http://math-atlas.sourceforge.net/faq.html.
   Alternatively, use fmatvec without ATLAS and install the reference
   libraries of BLAS and LAPACK from http://www.netlib.org/lapack.

2. Run aclocal, autoheader, autoconf, libtoolize und automake, if
   necessary:
   ```
   aclocal
   autoheader
   autoconf
   libtoolize
   automake -a
   ```
		                
3. Run configure

   `./configure`

   Hereby, the following options are important:

   `--prefix=PFX` (prefix, where the library will be installed).

   `--enable-atlas` (use ATLAS)

   `--with-blas-lib-prefix=PFX` (prefix, where the BLAS lib is
       installed, when ATLAS is not used)

   `--with-lapack-lib-prefix=PFX` (prefix, where the LAPACK lib is
       installed, when ATLAS is not used)

   `--with-atlas-inc-prefix=PFX` (prefix, where the ATLAS includes 
       are installed, when ATLAS is used)

   `--with-atlas-lib-prefix=PFX` (prefix, where the ATLAS libs are
       installed, when ATLAS is used)

4. Run make in order to build the shared library libfmatvec.so.
          
   `make`

5. Run make install in order to install the library libfmatvec.so.

   `make install`
				  
# Usage

1. To compile a program that includes the fmatvec headers use

   `pkg-config --cflags fmatvec`

2. To link a program against the fmatvec library use

   `pkg-config --libs fmatvec`

# Documentation

See the doxygen documentation. Use

`make doc`

to create this documentation.
