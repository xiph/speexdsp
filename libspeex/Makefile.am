## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.3 2002/02/19 22:44:28 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

lib_LTLIBRARIES = libspeex.la

# Sources for compilation in the library
libspeex_la_SOURCES = speex.c \
	lpc.c \
	ltp.c \
	lsp.c \
	quant_lsp.c \
	lsp_tables_nb.c \
	gain_table.c \
	cb_search.c \
	stoc.c

include_HEADERS = speex.h

noinst_HEADERS = lsp.h \
	lpc.h \
	ltp.h \
	quant_lsp.h \
	cb_search.h
	
	
libspeex_la_LDFLAGS = -release $(LT_RELEASE)

bin_PROGRAMS = testenc

testenc_SOURCES = testenc.c

testenc_LDADD = libspeex.la
