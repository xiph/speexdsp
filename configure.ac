dnl Process this file with autoconf to produce a configure script. -*-m4-*-

AC_INIT(libspeex/speex.c)

AM_CONFIG_HEADER([config.h])

SPEEX_MAJOR_VERSION=1
SPEEX_MINOR_VERSION=1
SPEEX_MICRO_VERSION=8
SPEEX_EXTRA_VERSION=
#SPEEX_VERSION=1.1.7
SPEEX_VERSION=$SPEEX_MAJOR_VERSION.$SPEEX_MINOR_VERSION.$SPEEX_MICRO_VERSION$SPEEX_EXTRA_VERSION

SPEEX_LT_CURRENT=3
SPEEX_LT_REVISION=0
SPEEX_LT_AGE=2

AC_SUBST(SPEEX_LT_CURRENT)
AC_SUBST(SPEEX_LT_REVISION)
AC_SUBST(SPEEX_LT_AGE)

# For automake.
VERSION=$SPEEX_VERSION
PACKAGE=speex

AC_SUBST(SPEEX_VERSION)

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define)
AM_MAINTAINER_MODE

AC_CANONICAL_HOST
AM_PROG_LIBTOOL

AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT

AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)

if test x$has_char16 = "xyes" ; then
        case 1 in
                $ac_cv_sizeof_short) SIZE16="short";;
                $ac_cv_sizeof_int) SIZE16="int";;
        esac
else
        case 2 in
                $ac_cv_sizeof_short) SIZE16="short";;
                $ac_cv_sizeof_int) SIZE16="int";;
        esac
fi

if test x$has_char16 = "xyes" ; then
        case 2 in
                $ac_cv_sizeof_int) SIZE32="int";;
                $ac_cv_sizeof_long) SIZE32="long";;
                $ac_cv_sizeof_short) SIZE32="short";;
        esac
else
        case 4 in
                $ac_cv_sizeof_int) SIZE32="int";;
                $ac_cv_sizeof_long) SIZE32="long";;
                $ac_cv_sizeof_short) SIZE32="short";;
        esac
fi

AC_SUBST(SIZE16)
AC_SUBST(SIZE32)


AC_CHECK_HEADERS(sys/soundcard.h sys/audioio.h)

XIPH_PATH_OGG([src="src"], [src=""])
AC_SUBST(src)

AC_CHECK_LIB(m, sin)

# Check for getopt_long; if not found, use included source.
AC_CHECK_FUNCS([getopt_long],,
[# FreeBSD has a gnugetopt library.
  AC_CHECK_LIB([gnugetopt],[getopt_long],
[AC_DEFINE([HAVE_GETOPT_LONG])],
[# Use the GNU replacement.
AC_LIBOBJ(getopt)
AC_LIBOBJ(getopt1)])])

AC_CHECK_LIB(winmm, main)

AC_DEFINE_UNQUOTED(SPEEX_VERSION, "${SPEEX_VERSION}", [Complete version string])
AC_DEFINE_UNQUOTED(SPEEX_MAJOR_VERSION, ${SPEEX_MAJOR_VERSION}, [Version major])
AC_DEFINE_UNQUOTED(SPEEX_MINOR_VERSION, ${SPEEX_MINOR_VERSION}, [Version minor])
AC_DEFINE_UNQUOTED(SPEEX_MICRO_VERSION, ${SPEEX_MICRO_VERSION}, [Version micro])
AC_DEFINE_UNQUOTED(SPEEX_EXTRA_VERSION, "${SPEEX_EXTRA_VERSION}", [Version extra])

AC_ARG_ENABLE(wideband, [  --disable-wideband       disable wideband codec],
[if test "$enableval" = no; then
  AC_DEFINE([DISABLE_WIDEBAND], , [Disable wideband codec])
fi])

AC_ARG_ENABLE(valgrind, [  --enable-valgrind       enable valgrind extra checks],
[if test "$enableval" = yes; then
  AC_DEFINE([ENABLE_VALGRIND], , [Enable valgrind extra checks])
fi])

AC_ARG_ENABLE(sse, [  --enable-sse            enable SSE support], [if test "$enableval" = yes; then
AC_DEFINE([_USE_SSE], , [Enable SSE support])
CFLAGS="$CFLAGS -O3 -msse"
fi
])

AC_ARG_ENABLE(fixed-point, [  --enable-fixed-point    compile as fixed-point],
[if test "$enableval" = yes; then
  AC_DEFINE([FIXED_POINT], , [Compile as fixed-point])
fi])

AC_ARG_ENABLE(arm4-asm, [  --enable-arm4-asm        make use of ARM4 assembly instructions],
[if test "$enableval" = yes; then
  AC_DEFINE([ARM4_ASM], , [Make use of ARM4 assembly instructions])
fi])

AC_ARG_ENABLE(arm5e-asm, [  --enable-arm5e-asm        make use of ARM5E assembly instructions],
[if test "$enableval" = yes; then
  AC_DEFINE([ARM5_ASM], , [Make use of ARM5E assembly instructions])
fi])

AC_ARG_ENABLE(fixed-point-debug, [  --enable-fixed-point-debug  Debug fixed-point implementation],
[if test "$enableval" = yes; then
  AC_DEFINE([FIXED_DEBUG], , [Debug fixed-point implementation])
fi])

AC_ARG_ENABLE(epic-48k, [  --enable-epic-48k       enable support for Epic 4.8 kbps mode],
[if test "$enableval" = yes; then
  AC_DEFINE([EPIC_48K], , [Enable support for Epic 4.8 kbps mode])
fi])

AC_ARG_ENABLE(ti-c55x, [  --enable-ti-c55x         enable support for TI C55X DSP],
[if test "$enableval" = yes; then
  AC_DEFINE([TI_C55X], , [Enable support for TI C55X DSP])
fi])

dnl Output the makefiles and version.h.

AC_OUTPUT([Makefile libspeex/Makefile src/Makefile doc/Makefile Speex.spec
           include/Makefile include/speex/Makefile speex.pc
           win32/Makefile win32/libspeex/Makefile win32/speexenc/Makefile
           win32/speexdec/Makefile include/speex/speex_types.h ])

if test "x$src" = "x"; then 
echo "**IMPORTANT**"
echo "You don't seem to have the development package for libogg (libogg-devel) installed. Only the Speex library (libspeex) will be built (no encoder/decoder executable)"
echo "You can download libogg from http://www.vorbis.com/download.psp"
fi

echo "Type \"make; make install\" to compile and install Speex";
