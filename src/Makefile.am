## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.7 2002/08/23 05:47:07 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

INCLUDES = -I$(top_srcdir)/libspeex

EXTRA_DIST = getopt_win.h getopt.c getopt1.c wave_out.c wave_out.h

include_HEADERS = 
noinst_HEADERS = wav_io.h

bin_PROGRAMS = speexenc speexdec

speexenc_SOURCES = speexenc.c wav_io.c
speexenc_LDADD = $(top_srcdir)/libspeex/libspeex.la $(OGG_LDFLAGS) $(OGG_LIBS)

speexdec_SOURCES = speexdec.c wav_io.c
speexdec_LDADD = $(top_srcdir)/libspeex/libspeex.la $(OGG_LDFLAGS) $(OGG_LIBS)
