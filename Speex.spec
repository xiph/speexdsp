%define name     Speex
%define ver      0.2.0
%define rel      1

Summary: An open-source, patent-free speech codec
Name: %name
Version: %ver
Release: %rel
Copyright: LGPL
Group: Application/Devel
Source: http://prdownloads.sourceforge.net/speex/%{name}-%{ver}.tar.gz
URL: http://speex.sourceforge.net/
Vendor: Speex
Packager: Jean-Marc Valin (jean-marc.valin@hermes.usherb.ca)

%description
Speex is a patent-free audio codec designed especially for voice (unlike 
Vorbis which targets general audio) signals and providing good narrowband 
and wideband quality. This project aims to be complementary to the Vorbis
codec.

%prep
%setup

%build
export CFLAGS='-O3'
./configure --prefix=/usr --enable-shared --enable-static
make

%install
make install

%files
/usr/lib/libspeex*
/usr/bin/speexenc
/usr/bin/speexdec
/usr/include/speex.h
/usr/include/speex_bits.h
