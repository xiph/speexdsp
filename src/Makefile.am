## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.2 2002/03/27 06:30:37 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

INCLUDES = -I$(top_srcdir)/libspeex

include_HEADERS = 
noinst_HEADERS = 

bin_PROGRAMS = speexenc speexdec

speexenc_SOURCES = speexenc.c
speexenc_LDADD = $(top_srcdir)/libspeex/libspeex.la

speexdec_SOURCES = speexdec.c
speexdec_LDADD = $(top_srcdir)/libspeex/libspeex.la
