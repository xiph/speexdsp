## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# $Id: Makefile.am,v 1.26 2002/05/03 19:27:52 jmvalin Exp $

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

lib_LTLIBRARIES = libspeex.la

# Sources for compilation in the library
libspeex_la_SOURCES = nb_celp.c \
	sb_celp.c \
	lpc.c \
	ltp.c \
	lsp.c \
	quant_lsp.c \
	lsp_tables_nb.c \
	lsp_tables_wb.c \
	gain_table.c \
	cb_search.c \
	stoc.c \
	filters.c \
	bits.c \
	modes.c \
	exc_table.c \
	exc_wb_table.c \
	exc_gains_table.c \
	exc_gains_wb_table.c \
	exc_gains_wb2_table.c \
	vq.c \
	matrix.c \
	mpulse.c \
	exc_sb_table.c \
	high_lsp_tables.c \
	vbr.c \
	hexc_table.c


include_HEADERS =  speex.h \
	speex_bits.h

noinst_HEADERS = lsp.h \
	nb_celp.h \
	lpc.h \
	ltp.h \
	quant_lsp.h \
	cb_search.h \
	filters.h \
	stack_alloc.h \
	vq.h \
	matrix.h \
	mpulse.h \
	modes.h \
	sb_celp.h \
	vbr.h
	
	
libspeex_la_LDFLAGS = -release $(LT_RELEASE)

bin_PROGRAMS = testenc testenc_wb

testenc_SOURCES = testenc.c
testenc_LDADD = libspeex.la

testenc_wb_SOURCES = testenc_wb.c
testenc_wb_LDADD = libspeex.la
