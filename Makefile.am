## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = speex.pc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = speex.pc

EXTRA_DIST = Speex.spec Speex.spec.in Speex.kdevelop

#Fools KDevelop into including all files
SUBDIRS = libspeex include @src@ doc win32

DIST_SUBDIRS = libspeex include src doc win32

rpm: dist
	rpmbuild -ta ${PACKAGE}-${VERSION}.tar.gz
