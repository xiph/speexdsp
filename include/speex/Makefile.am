# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

install-data: installdirs
	$(INSTALL_DATA) $(top_builddir)/include/speex/speex_types.h $(includedir)/speex/speex_types.h

pkginclude_HEADERS =  speex.h \
	speex_bits.h \
	speex_header.h \
	speex_callbacks.h \
	speex_stereo.h \
	speex_preprocess.h \
	speex_jitter.h \
	speex_echo.h

