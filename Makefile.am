## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

m4datadir = $(datadir)/aclocal
m4data_DATA = speex.m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = speex.pc

EXTRA_DIST = Speex.spec Speex.spec.in speex.m4 speex.pc.in

SUBDIRS = libspeex include @src@ doc win32

DIST_SUBDIRS = libspeex include src doc win32

rpm: dist
	rpmbuild -ta ${PACKAGE}-${VERSION}.tar.gz
