## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.6 2002/08/23 04:38:22 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

INCLUDES = -I$(top_srcdir)/libspeex

include_HEADERS = 
noinst_HEADERS = wav_io.h

bin_PROGRAMS = speexenc speexdec

speexenc_SOURCES = speexenc.c wav_io.c
speexenc_LDADD = $(top_srcdir)/libspeex/libspeex.la $(OGG_LDFLAGS) $(OGG_LIBS)

speexdec_SOURCES = speexdec.c wav_io.c
speexdec_LDADD = $(top_srcdir)/libspeex/libspeex.la $(OGG_LDFLAGS) $(OGG_LIBS)
