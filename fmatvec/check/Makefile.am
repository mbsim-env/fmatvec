if COND_WIN32
else
  MAYBE_VALGRIND_TESTS = valgrindtestfunction.sh valgrindtestast.sh
endif

check_PROGRAMS = testfunction testast testsymfunction

TESTS = testfunction.sh testast.sh testsymfunction.sh $(MAYBE_VALGRIND_TESTS)

testfunction_SOURCES = testfunction.cc
testfunction_CPPFLAGS = -I$(top_srcdir)
testfunction_LDADD = ../libfmatvec.la

testast_SOURCES = testast.cc
testast_CPPFLAGS = -I$(top_srcdir)
testast_LDADD = ../libfmatvec.la

testsymfunction_SOURCES = testsymfunction.cc
testsymfunction_CPPFLAGS = -I$(top_srcdir)
testsymfunction_LDADD = ../libfmatvec.la
