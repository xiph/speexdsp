%define name     Speex
%define ver      0.6.0
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
BuildRoot: /var/tmp/%{name}-build-root
Docdir: /usr/share/doc

%description
Speex is a patent-free audio codec designed especially for voice (unlike 
Vorbis which targets general audio) signals and providing good narrowband 
and wideband quality. This project aims to be complementary to the Vorbis
codec.

%changelog
* Tue Jul 30 2002 Fredrik Rambris <boost@users.sourceforge.net> 0.5.2
- Added buildroot and docdir and ldconfig. Makes it builadble by non-roots
  and also doesn't write to actual library paths when building.

%prep
%setup

%build
export CFLAGS='-O3'
./configure --prefix=/usr --enable-shared --enable-static
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%defattr(-, root, root)
%doc AUTHORS COPYING NEWS TODO README
/usr/lib/libspeex*
/usr/bin/speexenc
/usr/bin/speexdec
/usr/include/speex.h
/usr/include/speex_bits.h
/usr/include/speex_header.h
