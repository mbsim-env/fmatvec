# Dashboard

Current Daily Build Status of the MBSim-Environment

| Build Type | Variant | Failed |
|------------|---------|--------|
| linux64-dailydebug | build | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-build.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-build.nrAll.svg)](https://www.mbsim-env.de/mbsim/linux64-dailydebug/report/result_current/) |
| linux64-dailydebug | examples | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-examples.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-examples.nrAll.svg)](https://www.mbsim-env.de/mbsim/linux64-dailydebug/report/result_current/runexamples_report/result_current/) |
| linux64-dailydebug | coverage | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-coverage.svg)](https://www.mbsim-env.de/mbsim/linux64-dailydebug/report/result_current/runexamples_report/result_current/coverage/) |
| linux64-dailydebug | examples-valgrind | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-valgrind-examples.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-valgrind-examples.nrAll.svg)](https://www.mbsim-env.de/mbsim/linux64-dailydebug/report/runexamples_valgrind_report/result_current/) |
| linux64-dailydebug | coverage-valgrind | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailydebug-valgrind-coverage.svg)](https://www.mbsim-env.de/mbsim/linux64-dailydebug/report/runexamples_valgrind_report/result_current/coverage/) |
| linux64-dailyrelease | build | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailyrelease-build.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailyrelease-build.nrAll.svg)](https://www.mbsim-env.de/mbsim/linux64-dailyrelease/report/result_current/) |
| linux64-dailyrelease | examples | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailyrelease-examples.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/linux64-dailyrelease-examples.nrAll.svg)](https://www.mbsim-env.de/mbsim/linux64-dailyrelease/report/result_current/runexamples_report/result_current/) |
| win64-dailyrelease | build | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/win64-dailyrelease-build.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/win64-dailyrelease-build.nrAll.svg)](https://www.mbsim-env.de/mbsim/win64-dailyrelease/report/result_current/) |
| win64-dailyrelease | examples | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/win64-dailyrelease-examples.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/win64-dailyrelease-examples.nrAll.svg)](https://www.mbsim-env.de/mbsim/win64-dailyrelease/report/result_current/runexamples_report/result_current/) |
| - | manuals | [![image](https://www.mbsim-env.de/mbsim/buildsystemstate/build-manuals.nrFailed.svg) / ![image](https://www.mbsim-env.de/mbsim/buildsystemstate/build-manuals.nrAll.svg)](https://www.mbsim-env.de/mbsim/doc_manualsbuild.log) |



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
