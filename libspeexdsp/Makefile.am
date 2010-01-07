# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

EXTRA_DIST=echo_diagnostic.m

INCLUDES = -I$(top_builddir)/include -I$(top_builddir) @FFT_CFLAGS@

lib_LTLIBRARIES = libspeexdsp.la

# Sources for compilation in the library
if BUILD_KISS_FFT
  FFTSRC=kiss_fft.c _kiss_fft_guts.h kiss_fft.h kiss_fftr.c kiss_fftr.h 
else
if BUILD_SMALLFT
  FFTSRC=smallft.c
else
  FFTSRC=
endif
endif

libspeexdsp_la_SOURCES = preprocess.c jitter.c mdf.c fftwrap.c filterbank.c resample.c buffer.c scal.c $(FFTSRC)

noinst_HEADERS = 	arch.h 	bfin.h 
		fixed_arm4.h \
		fixed_arm5e.h 	fixed_bfin.h 	fixed_debug.h 	\
		math_approx.h 		misc_bfin.h 	\
		stack_alloc.h 	fftwrap.h \
	filterbank.h fixed_generic.h os_support.h \
	pseudofloat.h smallft.h vorbis_psy.h resample_sse.h

libspeexdsp_la_LDFLAGS = -no-undefined -version-info @SPEEX_LT_CURRENT@:@SPEEX_LT_REVISION@:@SPEEX_LT_AGE@

noinst_PROGRAMS = testdenoise testecho testjitter
testdenoise_SOURCES = testdenoise.c
testdenoise_LDADD = libspeexdsp.la @FFT_LIBS@
testecho_SOURCES = testecho.c
testecho_LDADD = libspeexdsp.la @FFT_LIBS@
testjitter_SOURCES = testjitter.c
testjitter_LDADD = libspeexdsp.la @FFT_LIBS@
