AC_DEFUN(AC_FIND_FILE,
[
$3=NONE
for i in $2;
do
  for j in $1;
  do
    if test -r "$i/$j"; then
      $3=$i
      break 2
    fi
  done
done
])

AC_DEFUN(AC_PATH_LIBOGG,
[
OGG_LIBS="-logg"

AC_MSG_CHECKING([for libogg])

ac_ogg_includes=NONE ac_ogg_libraries=NONE ac_ogg_bindir=NONE
ogg_libraries=""
ogg_includes=""
AC_ARG_WITH(ogg-dir,
    [  --with-ogg-dir=DIR       where the root of OGG is installed ],
    [  ac_ogg_includes="$withval"/include
       ac_ogg_libraries="$withval"/lib
    ])

AC_ARG_WITH(ogg-includes,
    [  --with-ogg-includes=DIR  where the OGG includes are. ],
    [  
       ac_ogg_includes="$withval"
    ])
    
ogg_libs_given=no

AC_ARG_WITH(ogg-libraries,
    [  --with-ogg-libraries=DIR where the OGG library is installed.],
    [  ac_ogg_libraries="$withval"
       ogg_libs_given=yes
    ])

ogg_incdirs="/usr/include /usr/lib/ogg/include /opt/include /usr/local/ogg/include /usr/include/ogg /usr/include /usr/local/include"
if test ! "$ac_ogg_includes" = "NONE"; then
  ogg_incdirs="$ac_ogg_includes $ac_ogg_includes/.. $ogg_incdirs"
fi
AC_FIND_FILE(ogg/ogg.h, $ogg_incdirs, ogg_incdir)
echo "Ogg includes in $ogg_incdir"


ogg_libdirs="$ac_ogg_libraries /usr/lib/ogg/lib /usr/lib /opt/lib /usr/local/ogg/lib /usr/local/lib /usr/lib/ogg /usr/local/lib"
test -n "$OGGDIR" && ogg_libdirs="$OGGDIR/lib $OGGDIR $ogg_libdirs"
if test ! "$ac_ogg_libraries" = "NONE"; then
  ogg_libdirs="$ac_ogg_libraries $ogg_libdirs"
fi

test=NONE
ogg_libdir=NONE
for dir in $ogg_libdirs; do
  try="ls -1 $dir/libogg*"
  if test=`eval $try 2> /dev/null`; then ogg_libdir=$dir; break; else echo "tried $dir" >&AC_FD_CC ; fi
done

echo "Ogg libraries in $ogg_libdir"

if test "$ogg_libdir" = "NONE" || test "$ogg_incdir" = "NONE"; then
   have_libogg=no
else
   have_libogg=yes
   AC_DEFINE(HAVE_LIBOGG)
fi

OGG_INCLUDES="-I$ogg_incdir"
OGG_LDFLAGS="-L$ogg_libdir"


AC_SUBST(OGG_LIBS)
AC_SUBST(OGG_INCLUDES)
AC_SUBST(OGG_LDFLAGS)

])

AC_DEFUN(AC_PATH_OGG,
[
LIBOGG="-logg"

AC_MSG_CHECKING([for libogg])

LIBOGG="$LIBOGG"
ac_ogg_includes=NO ac_ogg_libraries=NO ac_ogg_bindir=NO
ogg_libraries=""
ogg_includes=""
AC_ARG_WITH(ogg-dir,
    [  --with-ogg-dir=DIR       where the root of OGG is installed ],
    [  ac_ogg_includes="$withval"/include
       ac_ogg_libraries="$withval"/lib
       ac_ogg_bindir="$withval"/bin
    ])

AC_ARG_WITH(ogg-includes,
    [  --with-ogg-includes=DIR  where the OGG includes are. ],
    [  
       ac_ogg_includes="$withval"
    ])
    
ogg_libs_given=no

AC_ARG_WITH(ogg-libraries,
    [  --with-ogg-libraries=DIR where the OGG library is installed.],
    [  ac_ogg_libraries="$withval"
       ogg_libs_given=yes
    ])
AC_CACHE_VAL(ac_cv_have_ogg,
[#try to guess OGG locations

ogg_incdirs="/usr/lib/ogg/include /opt/include /usr/local/ogg/include /usr/include/ogg /usr/include /usr/local/include $OGGINC"
test -n "$OGGDIR" && ogg_incdirs="$OGGDIR/include $OGGDIR $ogg_incdirs"
ogg_incdirs="$ac_ogg_includes $ogg_incdirs"
AC_FIND_FILE(ogg/ogg.h, $ogg_incdirs, ogg_incdir)
echo $ogg_incdir
ac_ogg_includes="$ogg_incdir"

ogg_libdirs="/usr/lib/ogg/lib /usr/lib /opt/lib /usr/local/ogg/lib /usr/local/lib /usr/lib/ogg /usr/local/lib $OGGLIB"
test -n "$OGGDIR" && ogg_libdirs="$OGGDIR/lib $OGGDIR $ogg_libdirs"
if test ! "$ac_ogg_libraries" = "NO"; then
  ogg_libdirs="$ac_ogg_libraries $ogg_libdirs"
fi

test=NONE
ogg_libdir=NONE
for dir in $ogg_libdirs; do
  try="ls -1 $dir/libogg*"
  if test=`eval $try 2> /dev/null`; then ogg_libdir=$dir; break; else echo "tried $dir" >&AC_FD_CC ; fi
done

ac_ogg_libraries="$ogg_libdir"

ac_cxxflags_safe="$CXXFLAGS"
ac_ldflags_safe="$LDFLAGS"
ac_libs_safe="$LIBS"

INCLUDE="$INCLUDE -I$ogg_incdir $all_includes"
LDFLAGS="-L$ogg_libdir $all_libraries"
#LIBS="$LIBS $LIBOGG"

CXXFLAGS="$ac_cxxflags_safe"
LDFLAGS="$ac_ldflags_safe"
#LIBS="$ac_libs_safe"

if test "$ac_ogg_includes" = NO || test "$ac_ogg_libraries" = NO; then
  ac_cv_have_ogg="have_ogg=no"
  ac_ogg_notfound=""
  if test "$ac_ogg_includes" = NO; then
    if test "$ac_ogg_libraries" = NO; then
      ac_ogg_notfound="(headers and libraries)";
    else
      ac_ogg_notfound="(headers)";
    fi
  else
    ac_ogg_notfound="(libraries)";
  fi

else
  have_ogg="yes"
fi
])

eval "$ac_cv_have_ogg"

if test "$have_ogg" != yes; then
  AC_MSG_RESULT([$have_ogg]);
else
  ac_cv_have_ogg="have_ogg=yes \
    ac_ogg_includes=$ac_ogg_includes ac_ogg_libraries=$ac_ogg_libraries"
  AC_MSG_RESULT([libraries $ac_ogg_libraries, headers $ac_ogg_includes])
  
  ogg_libraries="$ac_ogg_libraries"
  ogg_includes="$ac_ogg_includes"
  AC_DEFINE(HAVE_LIBOGG)
fi

dnl if test ! "$ogg_libs_given" = "yes"; then
dnl CHECK_OGG_DIRECT(ogg_libraries= ,[])
dnl fi

AC_SUBST(ogg_libraries)
AC_SUBST(ogg_includes)

if test "$ogg_includes" = "/usr/include" || "$ogg_includes" = "$x_includes" || test -z "$ogg_includes"; then
 OGG_INCLUDES="";
else
 OGG_INCLUDES="-I$ogg_includes"
 all_includes="$OGG_INCLUDES $all_includes"
fi

if test "$ogg_libraries" = "$x_libraries" || test -z "$ogg_libraries"; then
 OGG_LDFLAGS=""
LIBOGG=""
else
 OGG_LDFLAGS="-L$ogg_libraries"
LIBOGG='-logg'
 all_libraries="$OGG_LDFLAGS $all_libraries"
fi

AC_SUBST(OGG_INCLUDES)
AC_SUBST(OGG_LDFLAGS)

AC_SUBST(LIBOGG)

])
