

dnl CHECK CGAL BEGIN
dnl This script takes too arguments, the action on success and the action on failure.
dnl It first checks if CGAL_MAKEFILE is defined (or a --with-cgalmakefile) value is passed.
dnl If that fails, it seaches for CGAL in the standard places.
dnl CGAL_CXXFLAGS, CGAL_CPPFLAGS, CGAL_LDFLAGS and CGAL_LIBS are all defined.
AC_DEFUN([ACX_CGAL],
[
acx_cgal_found=no
AC_ARG_WITH(cgalmakefile,
         [AC_HELP_STRING([--with-cgalmakefile=makefile], [Use the following CGAL makefile])])
case $with_cgalmakefile in
        yes | "") ;;
        no) acx_cgal_found=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) CGAL_MAKEFILE="$with_cgalmakefile" ;;
        *) CGAL_MAKEFILE="$with_cgalmakefile" ;;
esac

if test "$acx_cgal_found" == no; then
    AC_MSG_CHECKING(CGAL_MAKEFILE)

    if test \! -z "$CGAL_MAKEFILE"; then

	if test -e "$CGAL_MAKEFILE"; then
	    tname=`mktemp /tmp/cgal_makefile_dsrXXXXXX`
    
cat > $tname << _ACEOF
include $CGAL_MAKEFILE

cppflags:
	@echo \$(CGAL_CXXFLAGS)

cxxflags:
	@echo
ldflags:
	@echo \$(CGAL_LDFLAGS)
_ACEOF
	    CGAL_CPPFLAGS="`make -s -f $tname cppflags`"
	    CGAL_CXXFLAGS="`make -s -f $tname cxxflags`"
	    CGAL_LDFLAGST="`make -s -f $tname ldflags`"
	    for i in $CGAL_LDFLAGST; do
		if test `echo $i| grep -c ^-l`; then 
		    CGAL_LIBS="$CGAL_LIBS $i"
		else
		    CGAL_LDFLAGS="$CGAL_LDFLAGS $i"
		fi
	    done
	    rm -f $tname
	    AC_MSG_RESULT(yes)
	    acx_cgal_found=yes
	   dnl  echo CGAL_CPPFLAGS are $CGAL_CPPFLAGS
	   dnl  echo CGAL_LDFLAGS are $CGAL_LDFLAGS
	   dnl  echo CGAL_LIBS are $CGAL_LIBS
        else
	    AC_MSG_RESULT(invalid)
 	    AC_ERROR([CGAL_MAKEFILE defined, but the makefile does not exist.])
	fi
    else 
	AC_MSG_RESULT(not defined)
    fi
fi

if test "$acx_cgal_found" == no; then
	AC_CHECK_HEADER(CGAL/Exact_predicates_inexact_constructions_kernel.h, cgal_have_header=yes, cgal_have_header=no)
	if test "$cgal_have_header" == yes; then
		AC_CHECK_LIB(CGAL, main, cgal_have_lib=yes, cgal_have_lib=no)
		if test "$cgal_have_lib" == no; then
			save_LIBS="$LIBS"; LIBS="$LIBS -lgmp -lmpfr -lm"
        		AC_CHECK_LIB(CGAL, main, [CGAL_LIBS="-lCGAL -lgmp -lmpfr"
						  cgal_have_lib=yes], cgal_have_lib=no)
        		LIBS="$save_LIBS"
		else
			CGAL_LIBS="-lCGAL"
			AC_CHECK_LIB(mpfr, main, [CGAL_LIBS="$CGAL_LIBS -lmpfr"])
			AC_CHECK_LIB(gmp, main, [CGAL_LIBS="$CGAL_LIBS -lgmp"])
			AC_CHECK_LIB(gmpxx, main, [CGAL_LIBS="$CGAL_LIBS -lgmpxx"])
		fi

		if test "$cgal_have_lib" == yes; then 
			acx_cgal_found=yes
		fi
	fi 
	if test "$acx_cgal_found" == yes; then 
		AC_CHECK_LIB(Core, main, [CGAL_LIBS="$CGAL_LIBS -lCore"])
	fi
fi



AC_MSG_CHECKING(CGAL)
AC_SUBST(CGAL_MAKEFILE)
AC_SUBST(CGAL_CXXFLAGS)
AC_SUBST(CGAL_CPPFLAGS)
AC_SUBST(CGAL_LDFLAGS)
AC_SUBST(CGAL_LIBS)
if test "$acx_cgal_found" == yes; then 
	AC_MSG_RESULT(yes)
	$1
else
	AC_MSG_RESULT(no)
	$2
fi])


dnl CHECK CGAL END








dnl CHECK CGALQt BEGIN
dnl This script checks for Qt support in CGAL.
dnl It takes two arguments, the action on success and the action on failure. 
dnl It first checks for a CGAL_MAKEFILE.
dnl If not, it looks for CGALQt in the usual places.
dnl CGALQT_CXXFLAGS, CGALQT_CPPFLAGS, CGALQT_LDFLAGS and CGALQT_LIBS are all defined. 
dnl If no CGAL_MAKEFILE is defined, it expects the variable QT_LIBS to be defined.
dnl This variable can be defined by ACX_QT (which is below)
AC_DEFUN([ACX_CGALQT],
[
acx_cgalqt_found=no

if test "$acx_cgalqt_found" == no; then

    if test \! -z "$CGAL_MAKEFILE" -a -e "$CGAL_MAKEFILE"; then
	tname=`mktemp /tmp/cgal_makefile_dsrXXXXXX`
    
cat > $tname << _ACEOF
include $CGAL_MAKEFILE
cppflags:
	@echo \$(CGAL_CXXFLAGS)

ldflags:
	@echo \$(CGAL_QT_LDFLAGS)
_ACEOF
	CGALQT_CPPFLAGST="`make -s -f $tname cppflags`"
	CGALQT_LDFLAGST="`make -s -f $tname ldflags`"
	if test `echo $CGALQT_CPPFLAGST | grep -c -e -DCGAL_USE_QT` == 1; then
    	    for i in $CGALQT_LDFLAGST; do
	        if test `echo $i| grep -c ^-l`; then 
		    CGALQT_LIBS="$CGALQT_LIBS $i"
	        else
		    CGALQT_LDFLAGS="$CGALQT_LDFLAGS $i"
	        fi
	    done
	    acx_cgalqt_found=yes
	fi
        rm -f $tname
    fi
fi

if test "$acx_cgalqt_found" == no; then 
    save_LIBS="$LIBS"; LIBS="$LIBS $CGAL_LIBS $CGAL_LDFLAGS $QT_LIBS"
    AC_CHECK_LIB(CGALQt, main, [acx_cgalqt_found=yes; CGALQT_LIBS="-lCGALQt"])
    LIBS="$save_LIBS"
fi
AC_MSG_CHECKING(Qt support for CGAL)
if test "$acx_cgalqt_found" == yes; then
    AC_MSG_RESULT(yes);
    $1
else
    AC_MSG_RESULT(no);
    $2
fi
])

dnl CHECK CGALQt END


dnl This script checks for Qt support in CGAL.
dnl It takes two arguments, the action on success and the action on failure. 
dnl It uses 
AC_DEFUN([ACX_QT],
[
gw_CHECK_QT

if test "$QT_MAJOR" == 3; then
   if test "$QT_IS_MT" == "yes"; then
     	QT_LIBS="-lqt-mt"
   else
      	QT_LIBS="-lqt"
   fi
   QT_LDFLAGS="$QT_LDADD"
   $1
else
   $2
fi

])




# This is taken from http://autoqt.sourceforge.net/
# Copyright (c) 2002, Geoffrey Wossum
# All rights reserved.
 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#  - Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.

#  - Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.

#  - Neither the name of Geoffrey Wossum nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.


# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Check for Qt compiler flags, linker flags, and binary packages
AC_DEFUN([gw_CHECK_QT],
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_PATH_X])

AC_MSG_CHECKING([QTDIR])
AC_ARG_WITH([qtdir], [  --with-qtdir=DIR        Qt installation directory [default=$QTDIR]], QTDIR=$withval)
# Check that QTDIR is defined or that --with-qtdir given
if test x"$QTDIR" = x ; then
    QT_SEARCH="/usr/lib/qt31 /usr/local/qt31 /usr/lib/qt3 /usr/local/qt3 /usr/lib/qt2 /usr/local/qt2 /usr/lib/qt /usr/local/qt"
    for i in $QT_SEARCH; do
        if test -f $i/include/qglobal.h -a x$QTDIR = x; then QTDIR=$i; fi
    done
fi
if test x"$QTDIR" = x ; then
    AC_MSG_ERROR([*** QTDIR must be defined, or --with-qtdir option given])
fi
AC_MSG_RESULT([$QTDIR])

# Change backslashes in QTDIR to forward slashes to prevent escaping
# problems later on in the build process, mainly for Cygwin build
# environment using MSVC as the compiler
# TODO: Use sed instead of perl
QTDIR=`echo $QTDIR | perl -p -e 's/\\\\/\\//g'`

# Figure out which version of Qt we are using
AC_MSG_CHECKING([Qt version])
QT_VER=`grep 'define.*QT_VERSION_STR\W' $QTDIR/include/qglobal.h | perl -p -e 's/\D//g'`
case "${QT_VER}" in
    2*)
        QT_MAJOR="2"
    ;;
    3*)
        QT_MAJOR="3"
    ;;
    *)
        AC_MSG_ERROR([*** Don't know how to handle this Qt major version])
    ;;
esac
AC_MSG_RESULT([$QT_VER ($QT_MAJOR)])

if test \! -z "$MOC"; then
	if test \! -e "$MOC"; then
		AC_MSG_ERROR([MOC defined but is not a valid file.])
	fi
fi

if test -z "$MOC"; then 
	AC_MSG_CHECKING(QTDIR/bin/moc)
	if test -e "$QTDIR/bin/moc"; then
		MOC=$QTDIR/bin/moc
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
	fi
fi

# Check that moc is in path
if test -z "$MOC"; then
	AC_CHECK_PROG(MOC, moc, moc, $PATH)
	if test x$MOC = x ; then
        	AC_MSG_ERROR([*** moc must be in path])
	fi
fi

# uic is the Qt user interface compiler
AC_CHECK_PROG(UIC, uic, uic)
if test x$UIC = x ; then
        AC_MSG_ERROR([*** uic must be in path])
fi

# qembed is the Qt data embedding utility.
# It is located in $QTDIR/tools/qembed, and must be compiled and installed
# manually, we'll let it slide if it isn't present
AC_CHECK_PROG(QEMBED, qembed, qembed)


# Calculate Qt include path
QT_CXXFLAGS="-I$QTDIR/include"

QT_IS_EMBEDDED="no"
# On unix, figure out if we're doing a static or dynamic link
case "${host}" in
    *-cygwin)
	AC_DEFINE_UNQUOTED(WIN32, "", Defined if on Win32 platform)
        if test -f "$QTDIR/lib/qt.lib" ; then
            QT_LIB="qt.lib"
            QT_IS_STATIC="yes"
            QT_IS_MT="no"
        elif test -f "$QTDIR/lib/qt-mt.lib" ; then
            QT_LIB="qt-mt.lib" 
            QT_IS_STATIC="yes"
            QT_IS_MT="yes"
        elif test -f "$QTDIR/lib/qt$QT_VER.lib" ; then
            QT_LIB="qt$QT_VER.lib"
            QT_IS_STATIC="no"
            QT_IS_MT="no"
        elif test -f "$QTDIR/lib/qt-mt$QT_VER.lib" ; then
            QT_LIB="qt-mt$QT_VER.lib"
            QT_IS_STATIC="no"
            QT_IS_MT="yes"
        fi
        ;;

    *)
        QT_IS_STATIC=`ls $QTDIR/lib/*.a 2> /dev/null`
        if test "x$QT_IS_STATIC" = x; then
            QT_IS_STATIC="no"
        else
            QT_IS_STATIC="yes"
        fi
        if test x$QT_IS_STATIC = xno ; then
            QT_IS_DYNAMIC=`ls $QTDIR/lib/*.so 2> /dev/null` 
            if test "x$QT_IS_DYNAMIC" = x;  then
                AC_MSG_ERROR([*** Couldn't find any Qt libraries])
            fi
        fi

        if test "x`ls $QTDIR/lib/libqt.* 2> /dev/null`" != x ; then
            QT_LIB="-lqt"
            QT_IS_MT="no"
        elif test "x`ls $QTDIR/lib/libqt-mt.* 2> /dev/null`" != x ; then
            QT_LIB="-lqt-mt"
            QT_IS_MT="yes"
        elif test "x`ls $QTDIR/lib/libqte.* 2> /dev/null`" != x ; then
            QT_LIB="-lqte"
            QT_IS_MT="no"
            QT_IS_EMBEDDED="yes"
        elif test "x`ls $QTDIR/lib/libqte-mt.* 2> /dev/null`" != x ; then
            QT_LIB="-lqte-mt"
            QT_IS_MT="yes"
            QT_IS_EMBEDDED="yes"
        fi
        ;;
esac
AC_MSG_CHECKING([if Qt is static])
AC_MSG_RESULT([$QT_IS_STATIC])
AC_MSG_CHECKING([if Qt is multithreaded])
AC_MSG_RESULT([$QT_IS_MT])
AC_MSG_CHECKING([if Qt is embedded])
AC_MSG_RESULT([$QT_IS_EMBEDDED])

QT_GUILINK=""
QASSISTANTCLIENT_LDADD="-lqassistantclient"
case "${host}" in
    *irix*)
        QT_LIBS="$QT_LIB"
        if test $QT_IS_STATIC = yes ; then
            QT_LIBS="$QT_LIBS -L$x_libraries -lXext -lX11 -lm -lSM -lICE"
        fi
        ;;

    *linux*)
        QT_LIBS="$QT_LIB"
        if test $QT_IS_STATIC = yes && test $QT_IS_EMBEDDED = no; then
            QT_LIBS="$QT_LIBS -L$x_libraries -lXext -lX11 -lm -lSM -lICE -ldl -ljpeg"
        fi
        ;;


    *osf*) 
        # Digital Unix (aka DGUX aka Tru64)
        QT_LIBS="$QT_LIB"
        if test $QT_IS_STATIC = yes ; then
            QT_LIBS="$QT_LIBS -L$x_libraries -lXext -lX11 -lm -lSM -lICE"
        fi
        ;;

    *solaris*)
        QT_LIBS="$QT_LIB"
        if test $QT_IS_STATIC = yes ; then
            QT_LIBS="$QT_LIBS -L$x_libraries -lXext -lX11 -lm -lSM -lICE -lresolv -lsocket -lnsl"
        fi
        ;;


    *win*)
        # linker flag to suppress console when linking a GUI app on Win32
        QT_GUILINK="/subsystem:windows"

	if test $QT_MAJOR = "3" ; then
	    if test $QT_IS_MT = yes ; then
        	QT_LIBS="/nodefaultlib:libcmt"
            else
            	QT_LIBS="/nodefaultlib:libc"
            fi
        fi

        if test $QT_IS_STATIC = yes ; then
            QT_LIBS="$QT_LIBS $QT_LIB kernel32.lib user32.lib gdi32.lib comdlg32.lib ole32.lib shell32.lib imm32.lib advapi32.lib wsock32.lib winspool.lib winmm.lib netapi32.lib"
            if test $QT_MAJOR = "3" ; then
                QT_LIBS="$QT_LIBS qtmain.lib"
            fi
        else
            QT_LIBS="$QT_LIBS $QT_LIB"        
            if test $QT_MAJOR = "3" ; then
                QT_CXXFLAGS="$QT_CXXFLAGS -DQT_DLL"
                QT_LIBS="$QT_LIBS qtmain.lib qui.lib user32.lib netapi32.lib"
            fi
        fi
        QASSISTANTCLIENT_LDADD="qassistantclient.lib"
        ;;

esac


if test x"$QT_IS_EMBEDDED" = "xyes" ; then
        QT_CXXFLAGS="-DQWS $QT_CXXFLAGS"
fi

if test x"$QT_IS_MT" = "xyes" ; then
        QT_CXXFLAGS="$QT_CXXFLAGS -D_REENTRANT -DQT_THREAD_SUPPORT"
fi

QT_LDADD="-L$QTDIR/lib $QT_LIBS"

if test x$QT_IS_STATIC = xyes ; then
    OLDLIBS="$LIBS"
    LIBS="$QT_LDADD"
    AC_CHECK_LIB(Xft, XftFontOpen, QT_LDADD="$QT_LDADD -lXft")
    LIBS="$LIBS"
fi

AC_MSG_CHECKING([QT_CXXFLAGS])
AC_MSG_RESULT([$QT_CXXFLAGS])
AC_MSG_CHECKING([QT_LDADD])
AC_MSG_RESULT([$QT_LDADD])

AC_SUBST(QT_CXXFLAGS)
AC_SUBST(QT_LDADD)
AC_SUBST(QT_GUILINK)
AC_SUBST(QASSISTANTCLIENT_LDADD)

])
