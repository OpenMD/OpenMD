dnl We need a function similar to AC_CHECK_LIB to check for C++ libraries.
dnl AC_CHECK_CXX_LIB provides a similar interface like AC_CHECK_LIB
dnl and uses AC_TRY_LINK.
dnl
dnl $1 library name (without "-l")
dnl $2 object name to check for
dnl $3 neccessary include directive(s)
dnl $4 command to create object $2
dnl $5 yes-action
dnl $6 no-action
dnl $7 include dir for $3 (with -I)
dnl $8 additional libraries to link with

AC_DEFUN(AC_CHECK_CXX_LIB, AC_MSG_CHECKING([for $2 in -l$1])
save_CXXFLAGS_CHECK_CXX_LIB="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS $7"
save_LIBS_CHECK_CXX_LIB="$LIBS"
LIBS="-l$1 $8 $LIBS"
[AC_TRY_LINK([$3], [$4], [AC_MSG_RESULT([yes])
$5], [AC_MSG_RESULT([no])
$6])]
CXXFLAGS="$save_CXXFLAGS_CHECK_CXX_LIB"
LIBS="$save_LIBS_CHECK_CXX_LIB")dnl

dnl
dnl AC_CHECK_MODSUFFIX
dnl
AC_DEFUN([AC_CHECK_MODSUFFIX],[

AC_MSG_CHECKING(for module suffix)
rm -f conftest*
# Intel ifc compiler generates files by the name of work.pc and work.pcl (!)
rm -f work*
cat >conftest.$ac_ext <<EOF
        module conftest
        integer n
        parameter (n=1)
        end module conftest
EOF
# SGI and absoft compilers generates module name in upper case!
testname="conftest"
modcase="lower"
if (eval $ac_compile) 2>/dev/null ; then
    MOD=`ls conftest* | grep -v conftest.$ac_ext | grep -v conftest.o`
    MOD=`echo "$MOD" | sed -e 's/conftest\.//g'`
    if test -z "$MOD" ; then
        MOD=`ls CONFTEST* 2>/dev/null \
                | grep -v CONFTEST.$ac_ext | grep -v CONFTEST.o`
        MOD=`echo "$MOD" | sed -e 's/CONFTEST\.//g'`
        if test -n "$MOD" ; then
            testname="CONFTEST"
            modcase="upper"
        fi
    fi
    if test -z "$MOD" ; then
        AC_MSG_RESULT(unknown)
        # Use mod if we can't figure it out
        MOD="mod"
    else
        AC_MSG_RESULT($MOD)
    fi
    if test -s work.pcl ; then
        AC_MSG_WARN([Compiler generates auxillery files!])
    fi
else
    AC_MSG_RESULT(unknown)
fi
AC_SUBST(MOD)

])

AC_DEFUN(ACX_CHECK_CC_FLAGS,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether ${CC} accepts $1, ac_cv_$2,
[echo 'void f(){}' > conftest.c
if test -z "`${CC} $1 -c conftest.c 2>&1`"; then
        ac_cv_$2=yes
else
        ac_cv_$2=no
fi
rm -f conftest*
])
if test "$ac_cv_$2" = yes; then
        :
        $3
else
        :
        $4
fi
])

AC_DEFUN(ACX_CHECK_CXX_FLAGS,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether ${CXX} accepts $1, ac_cv_$2,
[echo 'void f(){}' > conftest.cpp
if test -z "`${CXX} $1 -c conftest.cpp 2>&1`"; then
        ac_cv_$2=yes
else
        ac_cv_$2=no
fi
rm -f conftest*
])
if test "$ac_cv_$2" = yes; then
        :
        $3
else
        :
        $4
fi
])

AC_DEFUN(ACX_PROG_GCC_VERSION,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether we are using gcc $1.$2 or later, ac_cv_prog_gcc_$1_$2,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUC__ && !defined (__INTEL_COMPILER)
#  if (__GNUC__ > $1) || (__GNUC__ == $1 && __GNUC_MINOR__ >= $2)
     yes;
#  endif
#endif
EOF
if AC_TRY_COMMAND(${CC-cc} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gcc_$1_$2=yes
else
  ac_cv_prog_gcc_$1_$2=no
fi
])
if test "$ac_cv_prog_gcc_$1_$2" = yes; then
        :
        $3
else
        :
        $4
fi
])

AC_DEFUN(ACX_PROG_GXX_VERSION,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether we are using g++ $1.$2 or later, ac_cv_prog_gxx_$1_$2,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.cpp <<EOF
#ifdef __GNUC__ && !defined (__INTEL_COMPILER)
#  if (__GNUC__ > $1) || (__GNUC__ == $1 && __GNUC_MINOR__ >= $2)
     yes;
#  endif
#endif
EOF
if AC_TRY_COMMAND(${CXX-c++} -E conftest.cpp) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gxx_$1_$2=yes
else
  ac_cv_prog_gxx_$1_$2=no
fi
])
if test "$ac_cv_prog_gxx_$1_$2" = yes; then
        :
        $3
else
        :
        $4
fi
])

AC_DEFUN(ACX_PROG_REALLY_GCC,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether we are *really* using GNU cc, ac_cv_prog_really_gcc,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUC__
  #if defined(__INTEL_COMPILER) || defined(__PATHCC__)
     no;
  #else
     yes;
  #endif
#endif
EOF
if AC_TRY_COMMAND(${CC-cc} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_really_gcc=yes
else
  ac_cv_prog_really_gcc=no
fi
])
if test "$ac_cv_prog_really_gcc" = yes; then
        :
        $1
else
        :
        $2
fi
])

AC_DEFUN(ACX_PROG_REALLY_GXX,
[
AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK(whether we are *really* using GNU c++, ac_cv_prog_really_gxx,
[
dnl The semicolon after "yes" below is to pacify NeXT's syntax-checking cpp.
cat > conftest.cpp <<EOF
#ifdef __GNUC__
  #if defined(__INTEL_COMPILER) || defined(__PATHCC__)
     no;
  #else
     yes;
  #endif
#endif
EOF
if AC_TRY_COMMAND(${CXX-c++} -E conftest.cpp) | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_really_gxx=yes
else
  ac_cv_prog_really_gxx=no
fi
])
if test "$ac_cv_prog_really_gxx" = yes; then
        :
        $1
else
        :
        $2
fi
])


AC_DEFUN(ACX_PROG_CC_MAXOPT,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_CANONICAL_HOST])

ACX_PROG_REALLY_GCC

# Try to determine "good" native compiler flags if none specified on command
# line
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$CC" = icc; then
                    CFLAGS="-O2"
                fi;;
  sparc-solaris2*) if test "$CC" = cc; then
                    CFLAGS="-O -dalign"
                 fi;;

  alpha*-osf*)  if test "$CC" = cc; then
                    CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host -arch host -std1"
                fi;;

  hppa*-hpux*)  if test "$CC" = cc; then
                    CFLAGS="-Ae +O3 +Oall"
                fi;;

   rs6000*-aix*)  if test "$CC" = cc -o "$CC" = xlc; then
                    CFLAGS="-O3 -qtune=auto -qansialias -w"
                fi;;
   powerpc*-aix*)
	if test "$CC" = cc -o "$CC" = xlc; then
        	CFLAGS="-O3 -qtune=auto -qansialias -w"
		echo "*******************************************************"
		echo "*  You have AIX on an unknown powerpc system.  It is  *"
		echo "*  recommended that you use                           *"
		echo "*                                                     *"
		echo "*    CFLAGS=-O3 -qarch=ppc -qtune=xxx -qansialias -w  *"
		echo "*                                 ^^^                 *"
		echo "*  where xxx is 601, 603, 604, or whatever kind of    *"
                echo "*  PowerPC CPU you have.   For more info, man cc.     *"
		echo "*******************************************************"
        fi;;
   *darwin*)
	if test "$CC" = xlc; then
        	CFLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto -qaltivec"
        fi
	if test "$CC" = icc; then
        	CFLAGS="-O2"
        fi
        if test $ac_cv_prog_really_gcc = yes; then
                CFLAGS="-Os"
        fi;;
  esac

  # use default flags for gcc on all systems
  if test $ac_cv_prog_really_gcc = yes -a -z "$CFLAGS"; then
     CFLAGS="-O2"
  fi

  if test -z "$CFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use  make CFLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of CFLAGS=-O will be used)     *"
	echo "********************************************************"
	echo ""
        CFLAGS="-O"
  fi

  ACX_CHECK_CC_FLAGS(${CFLAGS}, ac_cv_guessed_cflags, , [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use  make CFLAGS=..., or edit the top level Makefile *"
        echo "********************************************************"
        echo ""
        CFLAGS=""
  ])

fi
])

AC_DEFUN(ACX_PROG_CXX_MAXOPT,
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_CANONICAL_HOST])

ACX_PROG_REALLY_GXX

# Try to determine "good" native compiler flags if none specified on command
# line
if test "$ac_test_CXXFLAGS" != "set"; then
  CXXFLAGS=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$CXX" = icc -o "$CXX" = icpc; then
                    CXXFLAGS="-O2"
                fi;;
  sparc-solaris2*) if test "$CXX" = CC; then
                    CXXFLAGS="-features=extensions -O -dalign"
                 fi;;
   rs6000*-aix*)  if test "$CXX" = xlC; then
                    CXXFLAGS="-O3 -qarch=pwrx -qtune=pwrx -qansialias -w"
                fi;;
   powerpc*-aix*)
	if test "$CXX" = xlC; then
        	CXXFLAGS="-O3 -qarch=ppc -qansialias -w"
		echo "*******************************************************"
		echo "*  You have AIX on an unknown powerpc system.  It is  *"
		echo "*  recommended that you use                           *"
		echo "*                                                     *"
		echo "*  CXXFLAGS=-O3 -qarch=ppc -qtune=xxx -qansialias -w  *"
		echo "*                                 ^^^                 *"
		echo "*  where xxx is 601, 603, 604, or whatever kind of    *"
                echo "*  PowerPC CPU you have.   For more info, man cc.     *"
		echo "*******************************************************"
        fi;;
   *darwin*)
	if test "$CXX" = xlc++ -o "$CXX" = xlC ; then
        	CXXFLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto -qaltivec"
        fi
	if test "$CXX" = icpc; then
        	CXXFLAGS="-O2"
        fi
        if test $ac_cv_prog_really_gxx = yes; then
                CXXFLAGS="-Os"
        fi;;
  esac

  # use default flags for gcc on all systems
  if test $ac_cv_prog_really_gxx = yes -a -z "$CXXFLAGS"; then
     CXXFLAGS="-O2"
  fi

  if test -z "$CXXFLAGS"; then
	echo ""
	echo "**********************************************************"
        echo "* WARNING: Don't know the best CXXFLAGS for this system  *"
        echo "* Use  make CXXFLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of CXXFLAGS=-O will be used)     *"
	echo "**********************************************************"
	echo ""
        CXXFLAGS="-O"
  fi

  ACX_CHECK_CXX_FLAGS(${CXXFLAGS}, ac_cv_guessed_cxxflags, , [
	echo ""
        echo "**********************************************************"
        echo "* WARNING: The guessed CXXFLAGS don't seem to work with  *"
        echo "* your compiler.                                         *"
        echo "* Use  make CXXFLAGS=..., or edit the top level Makefile *"
        echo "**********************************************************"
        echo ""
        CXXFLAGS=""
  ])

fi
])

AC_DEFUN([adl_FUNC_GETOPT_LONG],
 [AC_PREREQ(2.49)dnl
  # clean out junk possibly left behind by a previous configuration
  rm -f src/getopt.h
  # Check for getopt_long support
  AC_CHECK_HEADERS([getopt.h])
  AC_CHECK_FUNCS([getopt_long],,
   [# FreeBSD has a gnugetopt library for this
    AC_CHECK_LIB([gnugetopt],[getopt_long],[AC_DEFINE([HAVE_GETOPT_LONG])],
     [# use the GNU replacement
      AC_LIBOBJ(getopt)
      AC_LIBOBJ(getopt1)
      AC_CONFIG_LINKS([src/getopt.h:src/utils/gnugetopt.h])])])])


AC_DEFUN([ACX_CONFIG_HOME], [
 myDir=${0%/*}
 if [ "$myDir" = "$0" ]; then
    # Ran from local directory
     myDir=$PWD
 fi
 # Resolve symlinks.
 myProgram="$0"
 while [ -L "$myProgram" ]; do
    ls=`/bin/ls -ld "$myProgram"`
    link=`/usr/bin/expr "$ls" : '.*-> \(.*\)$'`
    if /usr/bin/expr "$link" : '.*/.*' > /dev/null; then
       myProgram="$link"
    else
       myProgram="`AS_DIRNAME([$myProgram])`/$link"
    fi
 done
 myDir=`AS_DIRNAME([$myProgram])`
fi
CONFIG_HOME=$myDir
])

AC_DEFUN(BB_ENABLE_DOXYGEN,
[
AC_ARG_ENABLE(doxygen, [  --enable-doxygen        enable documentation generation with doxygen (auto)])
AC_ARG_ENABLE(dot, [  --enable-dot            use 'dot' to generate graphs in doxygen (auto)])              
AC_ARG_ENABLE(html-docs, [  --enable-html-docs      enable HTML generation with doxygen (yes)], [], [ enable_html_docs=yes])              
AC_ARG_ENABLE(latex-docs, [  --enable-latex-docs     enable LaTeX documentation generation with doxygen (no)], [], [ enable_latex_docs=no])              
if test "x$enable_doxygen" = xno; then
        enable_doc=no
else 
        AC_PATH_PROG(DOXYGEN, doxygen, , $PATH)
        if test "x$DOXYGEN" = x; then
                if test "x$enable_doxygen" = xyes; then
                        AC_MSG_ERROR([could not find doxygen])
                fi
                enable_doc=no
        else
                enable_doc=yes
                AC_PATH_PROG(DOT, dot, , $PATH)
        fi
fi

if test "x$enable_doc" = xyes; then
  DOC=1
else
  DOC=0
fi
AC_SUBST(DOC)

if test x$DOT = x; then
        if test "x$enable_dot" = xyes; then
                AC_MSG_ERROR([could not find dot])
        fi
        enable_dot=no
else
        enable_dot=yes
fi
AC_SUBST(enable_dot)
AC_SUBST(enable_html_docs)
AC_SUBST(enable_latex_docs)
])
#
#
#
AC_DEFUN([AX_SYS_PERLSHARPBANG],[dnl

   AC_PATH_PROG(PERLINTERP,perl,perl)
   ac_cv_path_perlinterp="$PERLINTERP"
   _sHpB='#!'

   AC_ARG_WITH(perl-shebang,
                AC_HELP_STRING([--with-perl-shebang],
           [override what perl thinks is the way for the kernel to start it (seldom needed)]dnl
		           ),
		[opt_perl_shebang="$withval"]dnl
		            ,dnl
		[opt_perl_shebang="not_set"]dnl
    )dnl

   AC_CACHE_CHECK([whether explicit instead of detected sharpbang is to be used],
		   ax_cv_opt_perl_shebang,
		  [ case "$opt_perl_shebang" in
		      not_set  ) ax_cv_opt_perl_shebang=''
		               ;;
		         *     )
	ax_cv_opt_perl_shebang=`echo "$opt_perl_shebang" | sed -e's|^#!\s*\(.*\)$|\1|'`
		    esac
		  ]dnl
    )dnl

   if test "A$ax_cv_opt_perl_shebang" != "A"
     then
       ac_cv_sys_kernshrpbang_perl="$ax_cv_opt_perl_shebang"
       PERL_SHEBANG="$ac_cv_sys_kernshrpbang_perl"
       AC_SUBST(PERL_SHEBANG)dnl
       AC_MSG_NOTICE([OK - PERL_SHEBANG is $_sHpB$PERL_SHEBANG.])

# Automatic detection of sharpbang formula starts here
     else dnl

   _somian_shbangperl=`$PERLINTERP -V:startperl`
   negclass="[[^']]"; dnl
# must leave this comment:  m4 will remove the outer brackets for us, heheh
   AC_CACHE_CHECK([for kernel sharpbang invocation to start perl],
	          ac_cv_sys_kernshrpbang_perl,
	[_somian_kspb_perl=`echo "$_somian_shbangperl" | sed -ne"s|.*='\($negclass*\)';$|\1|p"`
	if test "x$_somian_kspb_perl" == x
	  then _somian_ksbp_warn_empty='durnit'
	  else
	  case "A$_somian_kspb_perl" in
	         A#!*perl* )
           ac_cv_sys_kernshrpbang_perl=`echo "$_somian_kspb_perl" | sed -e's|#!\(.*\)$|\1|'`
			;;
	             A*    )  _somian_ksbp_warn_defau='trouble'
		              ac_cv_sys_kernshrpbang_perl="$PERLINTERP"
	  esac
	fi
])dnl Done with testing sharpbang

# The above prints Checking ... result message to user.
   PERL_SHEBANG="$ac_cv_sys_kernshrpbang_perl"
   AC_SUBST(PERL_SHEBANG)
    if test A${_somian_ksbp_warn_empty+set} == Aset
      then   AC_MSG_WARN([dnl
In last check, doing $PERLINTERP -V:startperl yielded empty result! That should not happen.])
    fi
# Inform user after printing result value
    if test A${_somian_ksbp_warn_defau+set} == Aset
      then AC_MSG_NOTICE([Maybe Not good -])
	   AC_MSG_WARN([dnl
In last check perl's Config query did not work so we bunted: $_sHpB$PERLINTERP])
      else AC_MSG_NOTICE([OK Good result - ])
	   AC_MSG_NOTICE([dnl
In last check we got a proper-looking answer from perl's Config: $_somian_shbangperl])
dnl Done with user info messages
    fi
dnl Outer loop checked for user override term here
  fi dnl

])dnl EOMACRO DEF

AC_DEFUN([ACX_CHECK_ZLIB],
#
# Handle user hints
#
[AC_ARG_WITH(zlib,
                AC_HELP_STRING([--with-zlib=DIR],
           [root directory path of zlib installation (defaults to /usr/local or /usr if not found in /usr/local)]dnl
		           ),
		[zlib_dir="$withval"]dnl
		            ,dnl
		[zlib_dir="not_set"]dnl
    )dnl

if test "$zlib_dir" != "no"; then

if test "$zlib_dir" != "not_set" ; then
  if test -d "$zlib_dir"
  then
    ZLIB_HOME="$zlib_dir"
  else
    AC_MSG_WARN([Sorry, $zlib_dir does not exist, checking usual places])
	ZLIB_HOME=/usr/local
	if test ! -f "${ZLIB_HOME}/include/zlib.h"
	then
        	ZLIB_HOME=/usr
	fi
  fi
fi
#
# Locate zlib, if wanted
#
if test -n "${ZLIB_HOME}"
then
        ZLIB_OLD_LDFLAGS=$LDFLAGS
        ZLIB_OLD_CFLAGS=$CFLAGS
        LDFLAGS="$LDFLAGS -L${ZLIB_HOME}/lib"
        CFLAGS="$CFLAGS -I${ZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [zlib_cv_libz=yes], [zlib_cv_libz=no])
        AC_CHECK_HEADER(zlib.h, [zlib_cv_zlib_h=yes], [zlib_cv_zlib_h=no])
        AC_LANG_RESTORE

        if test "$zlib_cv_libz" = "yes" -a "$zlib_cv_zlib_h" = "yes"; then
                AC_DEFINE(HAVE_ZLIB_H, 1, [have zlib.h])
                AC_DEFINE(HAVE_LIBZ, 1, [have libz.a])
                ZLIB_INC_DIR="${ZLIB_HOME}/include"
                ZLIB_LIB_DIR="${ZLIB_HOME}/lib"
                ZLIB="-lz"
        else
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                ZLIB_INC_DIR=
                ZLIB_LIB_DIR=
                ZLIB=
                LDFLAGS="$ZLIB_OLD_LDFLAGS"
                CFLAGS="$ZLIB_OLD_CFLAGS"
                AC_MSG_RESULT(failed)
	        echo ""
	        echo "*********************************************************"
                echo "* WARNING: Could not find a working zlib installation   *"
                echo "* If you need OpenMD to be able to deal with compressed *"
                echo "* trajectory dump files be sure to specify a valid zlib *"
	        echo "* installation with --with-zlib=DIR                     *"
                echo "*                                                       *"
                echo "* OpenMD will still work without zlib installed.        *"
	        echo "*********************************************************"
	        echo ""
        fi
        AC_SUBST(ZLIB_INC_DIR)
        AC_SUBST(ZLIB_LIB_DIR)
        AC_SUBST(ZLIB)
fi
fi
])

AC_DEFUN([ACX_CHECK_QHULL],
#       
# Handle user hints
#
[AC_ARG_WITH(qhull,
             AC_HELP_STRING([--with-qhull=DIR],
                            [root directory path of qhull installation (defaults to /usr/local or /usr if not found in /usr/local)]dnl
		           ),
             [qhull_dir="$withval"],
             [qhull_dir="not_set"]dnl
            )dnl

QHULL_INC_DIR=
QHULL_LIB_DIR=
QHULL=
USE_QHULL=no

if test "$qhull_dir" != "no"; then
   if test "$qhull_dir" != "not_set" ; then
     if test -d "$qhull_dir"; then
       QHULL_HOME="$qhull_dir"
     else
       AC_MSG_WARN([Sorry, $qhull_dir does not exist, checking usual places])
       QHULL_HOME=/usr/local
       if test ! -f "${QHULL_HOME}/include/qhull/libqhull.h"; then
          QHULL_HOME=/usr
       fi
     fi
   fi
   #
   # Locate qhull, if wanted
   #
   if test -n "${QHULL_HOME}"; then
        AC_MSG_NOTICE([Checking for qhull in ${QHULL_HOME}])
        AC_LANG_SAVE
        AC_LANG_C
        QHULL_OLD_LDFLAGS=$LDFLAGS
        QHULL_OLD_CFLAGS=$CFLAGS
        QHULL_OLD_CPPFLAGS=$CPPFLAGS
        LDFLAGS="$LDFLAGS -L${QHULL_HOME}/lib"
        CFLAGS="$CFLAGS -I${QHULL_HOME}/include"
        CPPFLAGS="$CPPFLAGS -I${QHULL_HOME}/include"
        AC_CHECK_HEADER(qhull/libqhull.h, [qhull_cv_libqhull_h=yes], [qhull_cv_libqhull_h=no])
        AC_CHECK_LIB(qhull, qh_qhull, [qhull_cv_libqhull=yes], [qhull_cv_libqhull=no])        
        AC_CHECK_LIB(qhull6, qh_qhull, [qhull_cv_libqhull6=yes], [qhull_cv_libqhull6=no])
        LDFLAGS="$QHULL_OLD_LDFLAGS"
        CFLAGS="$QHULL_OLD_CFLAGS"
        CPPFLAGS="$QHULL_OLD_CPPFLAGS"
        AC_LANG_RESTORE

        if test "$qhull_cv_libqhull_h" = "yes" -a "$qhull_cv_libqhull" = "yes" -o "$qhull_cv_libqhull6" = "yes"; then
           AC_DEFINE(HAVE_LIBQHULL_H, 1, [have libqhull.h])
           if test "$qhull_cv_libqhull" = "yes"; then
                AC_DEFINE(HAVE_QHULL, 1, [have libqhull.a])
                QHULL="-lqhull"
           fi
           if test "$qhull_cv_libqhull6" = "yes"; then
                AC_DEFINE(HAVE_QHULL, 1, [have libqhull6.a])
                QHULL="-lqhull6"
           fi
           USE_QHULL=yes                
           QHULL_INC_DIR="${QHULL_HOME}/include"
           QHULL_LIB_DIR="${QHULL_HOME}/lib"
           AC_MSG_RESULT([Working qhull found, will proceed.])
        else
	   AC_MSG_WARN([])
           AC_MSG_WARN([Could not find a working qhull installation])
           AC_MSG_WARN([If you need OpenMD to be able to deal with convex    ])
           AC_MSG_WARN([hulls be sure to specify a valid qhull installation ])
	   AC_MSG_WARN([with --with-qhull=DIR                               ])
           AC_MSG_WARN([])
           AC_MSG_WARN([OpenMD will still work without qhull installed.      ])
	   AC_MSG_WARN([])
        fi
    fi
fi
AC_SUBST(QHULL_INC_DIR)
AC_SUBST(QHULL_LIB_DIR)
AC_SUBST(QHULL)
AC_SUBST(USE_QHULL)
])

AC_DEFUN([ACX_CHECK_OPENBABEL],
#
# Handle user hints
#
[AC_ARG_WITH(openbabel,
                AC_HELP_STRING([--with-openbabel=DIR],
           [root directory path of openbabel-2.x installation (defaults to /usr/local or /usr if not found in /usr/local)]dnl
		           ),
		[openbabel_dir="$withval"]dnl
		            ,dnl
		[openbabel_dir="not_set"]dnl
    )dnl

if test "$openbabel_dir" != "no"; then

if test "$openbabel_dir" != "not_set" ; then
  if test -d "$openbabel_dir"
  then
    OPENBABEL_HOME="$openbabel_dir"
  else
    AC_MSG_WARN([Sorry, $openbabel_dir does not exist, checking usual places])
	OPENBABEL_HOME=/usr/local
	if test ! -f "${OPENBABEL_HOME}/include/openbabel-2.0/openbabel/babelconfig.h" -a -f "${OPENBABEL_HOME}/include/openbabel-2.0/openbabel/obconversion.h"
	then
		OPENBABEL_HOME=/usr
	fi
  fi
fi
#
# Locate openbabel, if wanted
#
if test -n "${OPENBABEL_HOME}"
then
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        OPENBABEL_OLD_LDFLAGS=$LDFLAGS
        OPENBABEL_OLD_CPPFLAGS=$CPPFLAGS
        LDFLAGS="$LDFLAGS -L${OPENBABEL_HOME}/lib -lopenbabel"
        CPPFLAGS="$CPPFLAGS -I${OPENBABEL_HOME}/include/openbabel-2.0"
        AC_CHECK_HEADER(openbabel/babelconfig.h, [openbabel_cv_openbabel_h=yes], [openbabel_cv_openbabel_h=no])
        AC_CHECK_HEADER(openbabel/obconversion.h, [openbabel_cv_obconversion_h=yes], [openbabel_cv_obconversion_h=no])
        AC_LINK_IFELSE([
             AC_LANG_PROGRAM(
                    [[
@%:@include <openbabel/babelconfig.h>
@%:@include <openbabel/obconversion.h>
using namespace std;
using namespace OpenBabel;
                    ]],
                    [[
OBConversion Conv(&cin, &cout);
                    ]]
                )],
                [
                openbabel_lib_found="yes"
                AC_MSG_RESULT([found])
                ],
                [
                openbabel_lib_found="no"
                AC_MSG_RESULT([not found])
                ]
            )
        AC_LANG_RESTORE
        LDFLAGS="$OPENBABEL_OLD_LDFLAGS"
        CPPFLAGS="$OPENBABEL_OLD_CPPFLAGS"

        if test "$openbabel_lib_found" = "yes" -a "$openbabel_cv_openbabel_h" = "yes" -a "$openbabel_cv_obconversion_h" = "yes"; then
                USE_OPENBABEL=yes                
                OPENBABEL_INC_DIR="${OPENBABEL_HOME}/include/openbabel-2.0"
                OPENBABEL_LIB_DIR="${OPENBABEL_HOME}/lib"
                OPENBABEL_LIB="-lopenbabel"
        else
                AC_MSG_CHECKING(openbabel in ${OPENBABEL_HOME})
                OPENBABEL_INC_DIR=
                OPENBABEL_LIB_DIR=
                OPENBABEL_LIB=
                USE_OPENBABEL=no
                AC_MSG_RESULT(failed)
	        echo ""
	        echo "*********************************************************"
                echo "* WARNING: Could not find a working openbabel-2.x       *"
                echo "* installation If you need OpenMD to be able to convert *"
                echo "* xyz or pdb files you need to specify a valid          *"
	        echo "* openbabel-2.x installation with --with-openbabel=DIR  *"
                echo "*                                                       *"
                echo "* OpenMD will still work without openbabel installed.   *"
	        echo "*********************************************************"
	        echo ""
        fi
        AC_SUBST(OPENBABEL_INC_DIR)
        AC_SUBST(OPENBABEL_LIB_DIR)
        AC_SUBST(OPENBABEL_LIB)
	AC_SUBST(USE_OPENBABEL)
fi
fi
])


AC_DEFUN([ACX_CHECK_FFTW],
#
# Handle user hints
#
[AC_ARG_WITH(fftw,
             AC_HELP_STRING([--with-fftw=DIR],
             [root directory path of fftw installation (defaults to /usr/local or /usr if not found in /usr/local)]dnl
		           ),
             [fftw_dir="$withval"]dnl
		            ,dnl
	     [fftw_dir="not_set"]dnl
    )dnl

if test "$fftw_dir" != "no"; then
  if test "$fftw_dir" != "not_set" ; then
    if test -d "$fftw_dir"; then
      FFTW_HOME="$fftw_dir"
    else
      AC_MSG_WARN([Sorry, $fftw_dir does not exist, checking usual places])
      FFTW_HOME=/usr/local
      if test ! -f "${FFTW_HOME}/include/fftw3.h" -o -f "${FFTW_HOME}/include/fftw.h" -o  -f "${FFTW_HOME}/include/dfftw.h"; then
        FFTW_HOME=/usr
      fi
    fi
    #
    # Locate fftw, if wanted
    #
    if test -n "${FFTW_HOME}"; then
      FFTW_OLD_LDFLAGS=$LDFLAGS
      FFTW_OLD_CFLAGS=$CFLAGS
      LDFLAGS="$LDFLAGS -L${FFTW_HOME}/lib"
      CFLAGS="$CFLAGS -I${FFTW_HOME}/include"
      AC_LANG_SAVE
      AC_LANG_C
      AC_CHECK_LIB(fftw3, fftw_execute, [fftw_cv_libfftw3=yes], [fftw_cv_libfftw3=no])
      AC_CHECK_HEADER(fftw3.h, [fftw_cv_fftw3_h=yes], [fftw_cv_fftw3_h=no])
      if test "$fftw_cv_libfftw3" = "no" -o "$fftw_cv_fftw3_h" = "no"; then
        AC_CHECK_LIB(fftw, fftwnd_one, [fftw_cv_libfftw=yes], [fftw_cv_libfftw=no])
        AC_CHECK_HEADER(fftw.h, [fftw_cv_fftw_h=yes], [fftw_cv_fftw_h=no])
        if test "$fftw_cv_libfftw" = "no" -o "$fftw_cv_fftw_h" = "no"; then
          AC_CHECK_LIB(dfftw, fftwnd_one, [fftw_cv_libdfftw=yes], [fftw_cv_libdfftw=no])
          AC_CHECK_HEADER(dfftw.h, [fftw_cv_dfftw_h=yes], [fftw_cv_dfftw_h=no])
        fi
      fi                      
      AC_LANG_RESTORE
      if test "$fftw_cv_libfftw3" = "yes" -a "$fftw_cv_fftw3_h" = "yes"; then
        AC_DEFINE(HAVE_FFTW3_H, 1, [have fftw3.h])
        FFTW_INC_DIR="${FFTW_HOME}/include"
        FFTW_LIB_DIR="${FFTW_HOME}/lib"
        FFTW_LIBS="-lfftw3"
      else
        if test "$fftw_cv_libfftw" = "yes" -a "$fftw_cv_fftw_h" = "yes"; then
          AC_DEFINE(HAVE_FFTW_H, 1, [have fftw.h])
          FFTW_INC_DIR="${FFTW_HOME}/include"
          FFTW_LIB_DIR="${FFTW_HOME}/lib"
          FFTW_LIBS="-lfftw"
        else
          if test "$fftw_cv_libdfftw" = "yes" -a "$fftw_cv_dfftw_h" = "yes"; then
            AC_DEFINE(HAVE_DFFTW_H, 1, [have dfftw.h])
            FFTW_INC_DIR="${FFTW_HOME}/include"
            FFTW_LIB_DIR="${FFTW_HOME}/lib"
            FFTW_LIBS="-ldfftw"
          else
            AC_MSG_CHECKING(fftw in ${FFTW_HOME})
            FFTW_INC_DIR=
            FFTW_LIB_DIR=
            FFTW_LIBS=
            LDFLAGS="$FFTW_OLD_LDFLAGS"
            CFLAGS="$FFTW_OLD_CFLAGS"
            AC_MSG_RESULT(failed)
	    echo ""
	    echo "*********************************************************"
            echo "* WARNING: Could not find a working FFTW installation   *"
            echo "* If you need the staticProps program to be able to     *"
            echo "* compute undulation spectra, be sure to specify a      *"
	    echo "* valid fftw installation with --with-fftw=DIR          *"
            echo "*                                                       *"
            echo "* OpenMD will still work without fftw installed.        *"
	    echo "*********************************************************"
	    echo ""
          fi
        fi
      fi
      AC_SUBST(FFTW_INC_DIR)
      AC_SUBST(FFTW_LIB_DIR)
      AC_SUBST(FFTW_LIBS)
    fi
  fi
fi
])


dnl
dnl These functions were taken from the GRASS GIS toolkit:
dnl 
dnl autoconf undefines "eval", so use "builtin([eval], ...)"
dnl

AC_DEFUN(LOC_CHECK_USE,[
AC_MSG_CHECKING(whether to use $2)
AC_MSG_RESULT("$with_$1")
case "$with_$1" in
	"no")	$3=	;;
	"yes")	$3="1"	;;
	*)	AC_MSG_ERROR([*** You must answer yes or no.])	;;
esac

])

AC_DEFUN(LOC_CHECK_INC_PATH,[
AC_MSG_CHECKING(for location of $2 includes)
case "$with_$1_includes" in
y | ye | yes | n | no)
	AC_MSG_ERROR([*** You must supply a directory to --with-$1-includes.])
	;;
esac
AC_MSG_RESULT($with_$1_includes)

if test -n "$with_$1_includes" ; then
    for dir in $with_$1_includes; do
        if test -d "$dir"; then
            $3="$$3 -I$dir"
        else
            AC_MSG_ERROR([*** $2 includes directory $dir does not exist.])
        fi
    done
fi
])

AC_DEFUN(LOC_CHECK_LIB_PATH,[
AC_MSG_CHECKING(for location of $2 library)
case "$with_$1_libs" in
y | ye | yes | n | no)
	AC_MSG_ERROR([*** You must supply a directory to --with-$1-libs.])
	;;
esac
AC_MSG_RESULT($with_$1_libs)

if test -n "$with_$1_libs"; then
    for dir in $with_$1_libs; do
        if test -d "$dir"; then
            $3="$$3 -L$dir"
        else
            AC_MSG_ERROR([*** $2 library directory $dir does not exist.])
        fi
    done
fi
])

AC_DEFUN(LOC_CHECK_SHARE_PATH,[
AC_MSG_CHECKING(for location of $2 data files)
case "$with_$1_share" in
y | ye | yes | n | no)
        AC_MSG_ERROR([*** You must supply a directory to --with-$1-share.])
        ;;
esac
AC_MSG_RESULT($with_$1_share)

if test -n "$with_$1_share" ; then
    if test -d "$with_$1_share"; then
        $3="$with_$1_share"
    else
        AC_MSG_ERROR([*** $2 data directory $dir does not exist.])
    fi
fi
])

AC_DEFUN(LOC_CHECK_INCLUDES,[
ac_save_cppflags="$CPPFLAGS"
CPPFLAGS="$3 $CPPFLAGS"
AC_CHECK_HEADERS($1, [], ifelse($4,[],[
    AC_MSG_ERROR([*** Unable to locate $2 includes.])
], $4))
CPPFLAGS=$ac_save_cppflags
])

dnl autoconf undefines "shift", so use "builtin([shift], ...)"

define(LOC_SHIFT1,[builtin([shift],$*)])
define(LOC_SHIFT2,[LOC_SHIFT1(LOC_SHIFT1($*))])
define(LOC_SHIFT4,[LOC_SHIFT2(LOC_SHIFT2($*))])
define(LOC_SHIFT8,[LOC_SHIFT4(LOC_SHIFT4($*))])
define(LOC_SHIFT9,[LOC_SHIFT1(LOC_SHIFT8($*))])

dnl $1  = library
dnl $2  = function
dnl $3  = descriptive name
dnl $4  = LDFLAGS initialiser
dnl $5  = result variable
dnl $6  = mandatory dependencies (not added to $5)
dnl $7  = mandatory dependencies (added to $5)
dnl $8  = ACTION-IF-NOT-FOUND
dnl $9+ = optional dependencies

define(LOC_CHECK_LIBS_0,[
AC_CHECK_LIB($1, $2, $5="$$5 -l$1 $7",[
[$8]
],$6 $7)
])

define(LOC_CHECK_LIBS_1,[
ifelse($9,[],
LOC_CHECK_LIBS_0($1,$2,,,$5,$6,$7,$8),
[
LOC_CHECK_LIBS_1($1,$2,,,$5,$6,$7,
LOC_CHECK_LIBS_1($1,$2,,,$5,$6,$7 $9,$8,LOC_SHIFT9($*)),
LOC_SHIFT9($*))
]
)
])

define(LOC_CHECK_LIBS,[
ac_save_ldflags="$LDFLAGS"
LDFLAGS="$4 $LDFLAGS"
LOC_CHECK_LIBS_1($1,$2,,,$5,$6,$7,
LDFLAGS=${ac_save_ldflags}
ifelse($8,[],[
    AC_MSG_ERROR([*** Unable to locate $3 library.])
],$8),LOC_SHIFT8($*))
LDFLAGS=${ac_save_ldflags}
])

AC_DEFUN(LOC_CHECK_VERSION_STRING,[
AC_MSG_CHECKING($3 version)
ac_save_cppflags="$CPPFLAGS"
CPPFLAGS="$5 $CPPFLAGS"
AC_TRY_RUN([
#include <stdio.h> 
#include <$1>
int main(void) {
 FILE *fp = fopen("conftestdata","w");
 fputs($2, fp);
 return 0;
}
],
[   $4=`cat conftestdata`
    AC_MSG_RESULT($$4)],
[   AC_MSG_ERROR([*** Could not determine $3 version.]) ],
[   $4=$6
    AC_MSG_RESULT([unknown (cross-compiling)]) ])
CPPFLAGS=$ac_save_cppflags
])

AC_DEFUN(LOC_CHECK_SHARE,[
AC_CHECK_FILE($3/$1, [], ifelse($4,[],[
    AC_MSG_ERROR([*** Unable to locate $2 data files.])
], $4))
])

AC_DEFUN(LOC_CHECK_VERSION_INT,[
AC_MSG_CHECKING($3 version)
ac_save_cppflags="$CPPFLAGS"
CPPFLAGS="$5 $CPPFLAGS"
AC_TRY_RUN([
#include <stdio.h>
#include <$1>
int main(void) {
 FILE *fp = fopen("conftestdata","w");
 fprintf(fp, "%d", $2);
 return 0;
}
    ],
    [   $4=`cat conftestdata`
        AC_MSG_RESULT($$4)],
    [   AC_MSG_ERROR([*** Could not determine $3 version.]) ],
    [   $4=$6
        AC_MSG_RESULT([unknown (cross-compiling)]) ])
CPPFLAGS=$ac_save_cppflags
])

dnl autoconf undefines "eval", so use "builtin([eval], ...)"

AC_DEFUN(LOC_PAD,[$1[]ifelse(builtin([eval],len($1) > 23),1,[
                          ],substr([                        ],len($1)))])

AC_DEFUN(LOC_ARG_WITH,[
AC_ARG_WITH($1,
LOC_PAD([  --with-$1])[support $2 functionality (default: ]ifelse([$3],,yes,[$3])[)],,
[with_$1=]ifelse([$3],,yes,[$3]))
])

AC_DEFUN(LOC_ARG_WITH_INC,[
AC_ARG_WITH($1-includes,
LOC_PAD([  --with-$1-includes=DIRS])[$2 include files are in DIRS])
])

AC_DEFUN(LOC_ARG_WITH_LIB,[
AC_ARG_WITH($1-libs,
LOC_PAD([  --with-$1-libs=DIRS])[$2 library files are in DIRS])
])

AC_DEFUN(LOC_ARG_WITH_SHARE,[
AC_ARG_WITH($1-share,
LOC_PAD([  --with-$1-share=DIR])[$2 data files are in DIR])
])

AC_DEFUN(LOC_OPTIONAL,[
AC_MSG_CHECKING(whether to build $1)
if test -n "$USE_$2" ; then
	AC_MSG_RESULT(yes)
	BUILD_$3="$4"
else
	AC_MSG_RESULT(no)
	BUILD_$3=
fi
AC_SUBST(BUILD_$3)
])

AC_DEFUN(LOC_MSG,[
echo "$1"
])

AC_DEFUN(LOC_PAD_26,[substr([                           ],len($1))])

AC_DEFUN(LOC_YES_NO,[if test -n "${$1}" ; then echo yes ; else echo no ; fi])

AC_DEFUN(LOC_MSG_USE,[
[echo "  $1:]LOC_PAD_26($1)`LOC_YES_NO($2)`"])
