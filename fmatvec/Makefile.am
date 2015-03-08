SUBDIRS = . check

lib_LTLIBRARIES = libfmatvec.la
libfmatvec_la_LDFLAGS = -version-info @LIBFMATVEC_VERSION_INFO@ -Wl,--disable-new-dtags,-rpath,\$$ORIGIN/../lib
libfmatvec_la_SOURCES = linear_algebra_double.cc \
			linear_algebra_complex.cc \
			atom.cc \
			_memory.cc \
			wrapper.cc

libfmatvec_la_CPPFLAGS = -I$(top_srcdir)
fmatvecincludedir=$(includedir)/fmatvec
fmatvecinclude_HEADERS = diagonal_matrix.h \
			 atom.h \
			 fmatvec.h \
			 general_matrix.h \
			 linear_algebra_double.h \
			 linear_algebra_complex.h \
			 linear_algebra.h \
			 matrix.h \
			 _memory.h \
			 row_vector.h \
			 square_matrix.h \
			 symmetric_matrix.h \
			 types.h \
			 vector.h \
			 band_matrix.h\
			 sparse_matrix.h\
			 range.h\
			 fixed_general_matrix.h\
			 fixed_square_matrix.h\
			 fixed_symmetric_matrix.h\
			 fixed_vector.h\
			 fixed_row_vector.h\
			 var_general_matrix.h\
			 var_square_matrix.h\
			 var_symmetric_matrix.h\
			 var_vector.h\
			 var_row_vector.h\
			 fixed_var_general_matrix.h\
			 var_fixed_general_matrix.h\
			 wrapper.h\
			 function.h
