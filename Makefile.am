## Process this file with automake to produce Makefile.in. -*-Makefile-*-

# Disable automatic dependency tracking if using other tools than gcc and gmake
#AUTOMAKE_OPTIONS = no-dependencies

EXTRA_DIST = Speex.spec Speex.spec.in Speex.kdevelop

#Fools KDevelop into including all files
if TRUE
SUBDIRS = libspeex @src@ doc win32
else
SUBDIRS = libspeex src doc win32
endif

DIST_SUBDIRS = libspeex src doc win32

rpm: dist
	rpmbuild -ta ${PACKAGE}-${VERSION}.tar.gz
