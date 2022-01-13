# GNU Automake file

lib_LTLIBRARIES = libcca.la

libcca_la_SOURCES = cca.c

pkginclude_HEADERS = cca.h

libcca_la_LDFLAGS = $(AM_LDFLAGS)

if BUILD_SHARED_LIBRARY
  AM_CFLAGS = -DBUILD_SHARED_LIBRARY
endif


# End of file
