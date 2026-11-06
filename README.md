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

The fmatvec (Fast Matrix-Vector Library) project aims to provide a high-performance C++ library for matrix-vector computations. Its object-oriented design enables convenient and intuitive handling of matrices and vectors, with a particular focus on efficient storage and management of different data types.

Most computational operations rely on ATLAS (Automatically Tuned Linear Algebra Software), which provides a portable and efficient implementation of BLAS through C and Fortran 77 interfaces, as well as selected LAPACK routines. For access to the full range of linear algebra functionality, the original LAPACK library is also required.

The library supports a wide range of data structures, including vectors, general matrices, symmetric matrices, banded matrices, and sparse matrices. It also supports computations involving dimensionless matrices and vectors. By leveraging C++ templates, fmatvec allows matrices and vectors to be instantiated with various data types, including integers, floating-point numbers (such as double), and complex numbers.



# The MBSim Environment

This repository is embed in the GitHub organization [MBSim-Env](https://github.com/mbsim-env/)

The MBSim Environment (MBSim-Env) is a scientific software project with
a long-standing history and a dedicated development and release infrastructure.

For further information about the project, its scientific background,
and its development, please visit the official project homepage:

[www.mbsim-env.de](https://www.mbsim-env.de/)

## Project Resources

The project maintains its own infrastructure for continuous integration,
automated testing, documentation generation, and official software releases.

The following resources are available on the official project website:

* [Documentation](https://www.mbsim-env.de/service/docu/) - Technical documentation and user guides.
* [Releases](https://www.mbsim-env.de/service/releases/) - Official software releases.
* [Nightly and CI Builds](https://www.mbsim-env.de/service/home/) - Development builds and continuous integration results.

## Project History and Maintenance

The official homepage also provides an overview of the scientific history and development of MBSim-Env.

It is primarily maintained by [Markus Friedrich](https://github.com/friedrichatgc).
