# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies


EXTRA_DIST=testenc.c testenc_wb.c testenc_uwb.c testdenoise.c testecho.c

INCLUDES = -I$(top_srcdir)/include -I$(top_builddir)/include -I$(top_builddir) @OGG_CFLAGS@

lib_LTLIBRARIES = libspeex.la

# Sources for compilation in the library
libspeex_la_SOURCES = nb_celp.c \
	sb_celp.c \
	lpc.c \
	ltp.c \
	lsp.c \
	quant_lsp.c \
	lsp_tables_nb.c \
	gain_table.c \
	gain_table_lbr.c \
	cb_search.c \
	filters.c \
	bits.c \
	modes.c \
	speex.c \
	vq.c \
	high_lsp_tables.c \
	vbr.c \
	hexc_table.c \
	exc_5_256_table.c \
	exc_5_64_table.c \
	exc_8_128_table.c \
	exc_10_32_table.c \
	exc_10_16_table.c \
	exc_20_32_table.c \
	hexc_10_32_table.c \
	misc.c \
	speex_header.c \
	speex_callbacks.c \
	math_approx.c \
	stereo.c \
	preprocess.c \
	smallft.c \
	lbr_48k_tables.c \
	jitter.c \
	mdf.c

noinst_HEADERS = lsp.h \
	nb_celp.h \
	lpc.h \
	lpc_bfin.h \
	ltp.h \
	quant_lsp.h \
	cb_search.h \
	filters.h \
	stack_alloc.h \
	vq.h \
	vq_sse.h \
	vq_arm4.h \
	vq_bfin.h \
	modes.h \
	sb_celp.h \
	vbr.h \
	misc.h \
	misc_bfin.h \
	ltp_sse.h \
	ltp_arm4.h \
	ltp_bfin.h \
	filters_sse.h \
	filters_arm4.h \
	filters_bfin.h \
	math_approx.h \
	smallft.h \
	arch.h \
	fixed_arm4.h \
	fixed_arm5e.h \
	fixed_bfin.h \
	fixed_debug.h \
	fixed_generic.h \
	cb_search_sse.h \
	cb_search_arm4.h \
	cb_search_bfin.h


libspeex_la_LDFLAGS = -version-info @SPEEX_LT_CURRENT@:@SPEEX_LT_REVISION@:@SPEEX_LT_AGE@

