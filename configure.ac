dnl -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(CCA, 1.1.0, software@scec.org)
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE(foreign tar-pax)

dnl
dnl Dependencies
dnl

dnl Etree (Euclid)
AC_ARG_ENABLE([etree],
    [AC_HELP_STRING([--enable-etree],
        [install Etree library @<:@default=yes@:>@])],
	[if test "$enableval" = yes ; then install_etree=yes; else install_etree=no; fi],
	[install_etree=yes])
AM_CONDITIONAL([INSTALL_ETREE], [test "$install_etree" = yes])
AC_ARG_WITH([etree-incdir],
  [AS_HELP_STRING([--with-etree-incdir], [directory containing Etree header files])],
  [with_etree_incdir="$withval"],
  [with_etree_incdir=no])
AC_SUBST(with_etree_incdir)
AC_ARG_WITH([etree-libdir],
  [AS_HELP_STRING([--with-etree-libdir], [directory containing Etree library])],
  [with_etree_libdir="$withval"],
  [with_etree_libdir=no])

dnl PROJ
AC_ARG_ENABLE([proj],
    [AC_HELP_STRING([--enable-proj],
        [install Proj @<:@default=yes@:>@])],
	[if test "$enableval" = "yes" ; then install_proj=yes; else install_proj=no; fi],
	[install_proj=yes])
AM_CONDITIONAL([INSTALL_PROJ], [test "$install_proj" = yes])
AC_ARG_WITH([proj-incdir],
    [AC_HELP_STRING([--with-proj-incdir], [directory containing Proj header files @<:@default=no@:>@])],
	[with_proj_incdir="$withval"],
	[with_proj_incdir=no])
AC_SUBST(with_proj_incdir)
AC_ARG_WITH([proj-libdir],
    [AC_HELP_STRING([--with-proj-libdir], [directory containing Proj library @<:@default=no@:>@])],
	[with_proj_libdir="$withval"],
	[with_proj_libdir=no])
AC_SUBST(with_proj_libdir)

dnl
dnl Command line tools
dnl

AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MKDIR_P

dnl CURL (used to download model files)
AC_PATH_PROG(CURL, curl)
if test -z "$CURL" ; then
  AC_MSG_FAILURE([cannot find 'curl' program.])
  CURL=`echo "Error: curl is not installed." ; false`
fi

dnl TAR
AC_PATH_PROG(TAR, tar)
if test -z "$TAR" ; then
  AC_MSG_FAILURE([cannot find 'tar' program.])
  TAR=`echo "Error: tar is not installed." ; false`
fi

dnl
dnl Setup environment so dependencies are used in build
dnl
CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib -L$prefix/lib64 $LDFLAGS"
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

dnl
dnl Verify configuration
dnl

dnl ETREE (required)
if test "$install_etree" = yes ; then
  ETREE_INCLUDES=""
  ETREE_LDFLAGS="-letree"
else
  if test "$with_etree_incdir" != no ; then
    ETREE_INCLUDES="-I$with_etree_incdir"
  fi
  if test "$with_etree_libdir" != no ; then
    ETREE_LDFLAGS="-L$with_etree_libdir -letree"
  fi
  SCEC_ETREE_HEADER
  SCEC_ETREE_LIB
fi
AC_SUBST(ETREE_INCLUDES)
AC_SUBST(ETREE_LDFLAGS)

dnl PROJ (required)
if test "$install_proj" = yes ; then
  PROJ_INCLUDES=""
  PROJ_LDFLAGS="-lproj"
else
  if test "$with_proj_incdir" != no ; then
    PROJ_INCLUDES="-I$with_proj_incdir"
  fi
  if test "$with_proj_libdir" != no ; then
    PROJ_LDFLAGS="-L$with_proj_libdir -lproj"
  fi
  SCEC_PROJ_HEADER
  SCEC_PROJ_LIB
fi
AC_SUBST(PROJ_INCLUDES)
AC_SUBST(PROJ_LDFLAGS)


AC_CONFIG_FILES([
  Makefile
  src/Makefile
  dependencies/Makefile
  data/Makefile
  ])

AC_OUTPUT

AM_CONDITIONAL([BUILD_SHARED_LIBRARY], [test "$enable_shared" = yes])
if test "$enable_shared" = yes; then
  AC_MSG_NOTICE([Building shared libraries. UCVM must be built with shared libraries for plugins to work.])
else
  AC_MSG_WARN([Building static libraries. Shared libraries are not enabled.])
fi


dnl End of file
