## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.1 2002/03/26 23:24:47 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

INCLUDES = -I$(top_srcdir)/libspeex

include_HEADERS = 
noinst_HEADERS = 

bin_PROGRAMS = speexenc

speexenc_SOURCES = speexenc.c
speexenc_LDADD = $(top_srcdir)/libspeex/libspeex.la
