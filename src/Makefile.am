## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.3 2002/05/15 21:47:51 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

INCLUDES = -I$(top_srcdir)/libspeex

include_HEADERS = 
noinst_HEADERS = 

bin_PROGRAMS = speexenc speexdec

speexenc_SOURCES = speexenc.c
speexenc_LDADD = $(top_srcdir)/libspeex/libspeex.la $(LIB_OGG)

speexdec_SOURCES = speexdec.c
speexdec_LDADD = $(top_srcdir)/libspeex/libspeex.la $(LIB_OGG)
