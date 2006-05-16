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


dnl 
dnl AC_CHECK_MODDIRFLAG
dnl
AC_DEFUN([AC_CHECK_MODDIRFLAG],[

# Check for module include path (some use -I, some (Solaris) use -M, some
# (absoft) use -p).
# Intel compilers use a wierd system: -cl,filename.pcl .  If no file is
# specified, work.pcl and work.pc are created.  However, if you specify
# a file, it must contain a the name of a file ending in .pc .  Ugh!
# Use the module made above
AC_MSG_CHECKING(for module directory path flag)
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
if (eval $ac_compile) 2>/dev/null ; then
   mod=`ls CONFTEST* 2>/dev/null | grep -v CONFTEST.$ac_ext | grep -v CONFTEST.o`
   mod=`echo "$mod" | sed -e 's/CONFTEST\.//g'`
   if test -n "$mod" ; then
      testname="CONFTEST"
   fi
   madedir=0
   if test ! -d conf ; then mkdir conf ; madedir=1; fi
   cp $testname.$MOD conf
   rm -f conftest* CONFTEST*
   cat >conftest1.$ac_ext <<EOF
        program main
        use conftest
        print *, n
        end
EOF
   F90_WORK_FILES_ARG=""
   F90MODINCSPEC=""
   if $FC -c -Iconf $FCFLAGS $FCFLAGS_SRCEXT conftest1.$ac_ext > conftest.out 2>&1 ; then
       MODDIRFLAG="-I"
       F90MODINCSPEC="-I<dir>"
       AC_MSG_RESULT(-I)
   elif $FC -c -Mconf $FCFLAGS $FCFLAGS_SRCEXT conftest1.$ac_ext >> conftest.out 2>&1 ; then
       MODDIRFLAG="-M"
       F90MODINCSPEC="-M<dir>"
       AC_MSG_RESULT(-M)
   elif $FC -c -pconf $FCFLAGS $FCFLAGS_SRCEXT conftest1.$ac_ext >> conftest.out 2>&1 ; then
       MODDIRFLAG="-p"
       F90MODINCSPEC="-p<dir>"
       AC_MSG_RESULT(-p)
   elif test -s work.pc ; then
        cp work.pc conf/mpimod.pc
        echo "mpimod.pc" > conf/mpimod.pcl
        echo "`pwd`/conf/mpimod.pc" >> conf/mpimod.pcl
        if $FC -c -cl,conf/mpimod.pcl $FCFLAGS $FCFLAGS_SRCEXT conftest1.$ac_ext >>conftest.out 2>&1 ; then
            MODDIRFLAG='-cl,mpimod.pcl'
            AC_MSG_RESULT([-cl,filename where filename contains a list of files and directories])
            F90_WORK_FILES_ARG="-cl,mpimod.pcl"
            F90MODINCSPEC="-cl,<dir>/<file>mod.pcl"
            AC_SUBST(F90_WORK_FILES_ARG)
        else
            # The version of the Intel compiler that I have refuses to let
            # you put the "work catalog" list anywhere but the current
            # directory. For example, you cannot in
         :
        fi
   fi
   if test -z "MODDIRFLAG" ; then
       AC_MSG_RESULT(unknown)
   fi
   AC_SUBST(MODDIRFLAG)
   AC_SUBST(F90MODINCSPEC)
   rm -f conftest* conf/conftest* conf/CONFTEST* CONFTEST* conf/mpimod*
   if test $madedir = 1 ; then rmdir conf ; fi
fi

])

AC_DEFUN(ACX_CHECK_CC_FLAGS,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether ${CC} accepts $1, ac_$2,
[echo 'void f(){}' > conftest.c
if test -z "`${CC} $1 -c conftest.c 2>&1`"; then
        ac_$2=yes
else
        ac_$2=no
fi
rm -f conftest*
])
if test "$ac_$2" = yes; then
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
AC_CACHE_CHECK(whether ${CXX} accepts $1, ac_$2,
[echo 'void f(){}' > conftest.cpp
if test -z "`${CXX} $1 -c conftest.cpp 2>&1`"; then
        ac_$2=yes
else
        ac_$2=no
fi
rm -f conftest*
])
if test "$ac_$2" = yes; then
        :
        $3
else
        :
        $4
fi
])

dnl -------------------------------------------------------------------------
dnl ACX_CHECK_FC_FLAGS()
dnl
dnl     Check for optimizer flags the Fortran compiler can use.
dnl
AC_DEFUN(ACX_CHECK_FC_FLAGS,
[
AC_CACHE_CHECK(whether ${FC} accepts $1, ac_$2,
[
AC_LANG_SAVE
AC_LANG(Fortran)
echo 'program main' > conftest.$ac_ext
echo 'end program main' >> conftest.$ac_ext
ac_compile='${FC} -c $1 $FCFLAGS $FCFLAGS_SRCEXT conftest.$ac_ext 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_compile); then
        ac_$2=yes
else
        ac_$2=no
fi
rm -f conftest*
AC_LANG_RESTORE()
])
if test "$ac_$2" = yes; then
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
                    CFLAGS="-O"
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
	if test "$CC" = xlc -o "$CC" = cc; then
        	CFLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto -qaltivec"
        fi
	if test "$CC" = icc; then
        	CFLAGS="-O3 -ip -no-prec-div -mdynamic-no-pic"
        fi;;
  esac

  # use default flags for gcc on all systems
  if test $ac_cv_prog_really_gcc = yes; then
     CFLAGS="-O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
  fi

  # test for gcc-specific flags:
  if test $ac_cv_prog_really_gcc = yes; then
    # -malign-double for x86 systems
    ACX_CHECK_CC_FLAGS(-malign-double,align_double, CFLAGS="$CFLAGS -malign-double")
    # -fstrict-aliasing for gcc-2.95+
    ACX_CHECK_CC_FLAGS(-fstrict-aliasing,fstrict_aliasing, CFLAGS="$CFLAGS -fstrict-aliasing")
  fi

  CPU_FLAGS=""
  if test $ac_cv_prog_really_gcc = yes; then
	  dnl try to guess correct CPU flags, at least for linux
	  case "${host_cpu}" in
	  i586*)  ACX_CHECK_CC_FLAGS(-mcpu=pentium,cpu_pentium,
			[CPU_FLAGS=-mcpu=pentium],
			[ACX_CHECK_CC_FLAGS(-mpentium,pentium,
				[CPU_FLAGS=-mpentium])])
		  ;;
	  i686*)  ACX_CHECK_CC_FLAGS(-mcpu=pentiumpro,cpu_pentiumpro,
			[CPU_FLAGS=-mcpu=pentiumpro],
			[ACX_CHECK_CC_FLAGS(-mpentiumpro,pentiumpro,
				[CPU_FLAGS=-mpentiumpro])])
		  ;;
	  powerpc*)
		cputype=`(grep cpu /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/ //g') 2> /dev/null`
		is60x=`echo $cputype | egrep "^60[0-9]e?$"`
		if test -n "$is60x"; then
			ACX_CHECK_CC_FLAGS(-mcpu=$cputype,m_cpu_60x,
				CPU_FLAGS=-mcpu=$cputype)
		elif test "$cputype" = 750; then
                        ACX_PROG_GCC_VERSION(2,95,
                                ACX_CHECK_CC_FLAGS(-mcpu=750,m_cpu_750,
					CPU_FLAGS=-mcpu=750))
		fi
		if test -z "$CPU_FLAGS"; then
		        ACX_CHECK_CC_FLAGS(-mcpu=powerpc,m_cpu_powerpc,
				CPU_FLAGS=-mcpu=powerpc)
		fi
		if test -z "$CPU_FLAGS"; then
			ACX_CHECK_CC_FLAGS(-mpowerpc,m_powerpc,
				CPU_FLAGS=-mpowerpc)
		fi
	  esac
  fi

  if test -n "$CPU_FLAGS"; then
        CFLAGS="$CFLAGS $CPU_FLAGS"
  fi

  if test -z "$CFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use  make CFLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        CFLAGS="-O3"
  fi

  ACX_CHECK_CC_FLAGS(${CFLAGS}, guessed_cflags, , [
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
                    CXXFLAGS="-O"
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
        	CXXFLAGS="-O3 -ip -no-prec-div -mdynamic-no-pic"
        fi;;
  esac

  # use default flags for gcc on all systems
  if test $ac_cv_prog_really_gxx = yes; then
     CXXFLAGS="-O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
  fi

  # test for gcc-specific flags:
  if test $ac_cv_prog_really_gxx = yes; then
    # -malign-double for x86 systems
    ACX_CHECK_CXX_FLAGS(-malign-double,align_double, CXXFLAGS="$CXXFLAGS -malign-double")
    # -fstrict-aliasing for gcc-2.95+
    ACX_CHECK_CXX_FLAGS(-fstrict-aliasing,fstrict_aliasing, CXXFLAGS="$CXXFLAGS -fstrict-aliasing")
  fi

  CPU_FLAGS=""
  if test $ac_cv_prog_really_gxx = yes; then
	  dnl try to guess correct CPU flags, at least for linux
	  case "${host_cpu}" in
	  i586*)  ACX_CHECK_CXX_FLAGS(-mcpu=pentium,cpu_pentium,
			[CPU_FLAGS=-mcpu=pentium],
			[ACX_CHECK_CXX_FLAGS(-mpentium,pentium,
				[CPU_FLAGS=-mpentium])])
		  ;;
	  i686*)  ACX_CHECK_CXX_FLAGS(-mcpu=pentiumpro,cpu_pentiumpro,
			[CPU_FLAGS=-mcpu=pentiumpro],
			[ACX_CHECK_CXX_FLAGS(-mpentiumpro,pentiumpro,
				[CPU_FLAGS=-mpentiumpro])])
		  ;;
	  powerpc*)
		cputype=`(grep cpu /proc/cpuinfo | head -1 | cut -d: -f2 | sed 's/ //g') 2> /dev/null`
		is60x=`echo $cputype | egrep "^60[0-9]e?$"`
		if test -n "$is60x"; then
			ACX_CHECK_CXX_FLAGS(-mcpu=$cputype,m_cpu_60x,
				CPU_FLAGS=-mcpu=$cputype)
		elif test "$cputype" = 750; then
                        ACX_PROG_GXX_VERSION(2,95,
                                ACX_CHECK_CXX_FLAGS(-mcpu=750,m_cpu_750,
					CPU_FLAGS=-mcpu=750))
		fi
		if test -z "$CPU_FLAGS"; then
		        ACX_CHECK_CXX_FLAGS(-mcpu=powerpc,m_cpu_powerpc,
				CPU_FLAGS=-mcpu=powerpc)
		fi
		if test -z "$CPU_FLAGS"; then
			ACX_CHECK_CXX_FLAGS(-mpowerpc,m_powerpc,
				CPU_FLAGS=-mpowerpc)
		fi
	  esac
  fi

  if test -n "$CPU_FLAGS"; then
        CXXFLAGS="$CXXFLAGS $CPU_FLAGS"
  fi

  if test -z "$CXXFLAGS"; then
	echo ""
	echo "**********************************************************"
        echo "* WARNING: Don't know the best CXXFLAGS for this system  *"
        echo "* Use  make CXXFLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of CXXFLAGS=-O3 will be used)    *"
	echo "**********************************************************"
	echo ""
        CXXFLAGS="-O3"
  fi

  ACX_CHECK_CXX_FLAGS(${CXXFLAGS}, guessed_cxxflags, , [
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

AC_DEFUN(ACX_PROG_FC_MAXOPT,
[
AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified on command
# line

if test "$ac_test_FFLAGS" != "set"; then
  FCFLAGS=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$FC" = ifc -o "$FC" = ifort; then
                    FCFLAGS="-O3 -ip -no-prec-div -cxxlib-icc"
                fi;;
   rs6000*-aix*)  if test "$FC" = xlf90 -o "$FC" = f90 -o "$FC" = xlf95; then
                    FCFLAGS="-O3 -qarch=pwrx -qtune=pwrx -qansialias -w"
                fi;;
   powerpc*-aix*)
	if test "$FC" = f90 -o "$FC" = xlf90 -o "$FC" = xlf95; then
        	FCFLAGS="-O3 -qarch=ppc -qansialias -w"
		echo "*******************************************************"
		echo "*  You have AIX on an unknown powerpc system.  It is  *"
		echo "*  recommended that you use                           *"
		echo "*                                                     *"
		echo "*   FCFLAGS=-O3 -qarch=ppc -qtune=xxx -qansialias -w  *"
		echo "*                                 ^^^                 *"
		echo "*  where xxx is 601, 603, 604, or whatever kind of    *"
                echo "*  PowerPC CPU you have.   For more info, man xlf.    *"
		echo "*******************************************************"
        fi;;
   *darwin*)
	if test "$FC" = f90 -o "$FC" = xlf90 -o "$FC" = xlf95; then
        	FCFLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto"
        fi
	if test "$FC" = ifort; then
        	FCFLAGS="-O3 -ip -no-prec-dev -mdynamic-no-pic"
        fi;;
  esac

  if test -n "$CPU_FLAGS"; then
        FCFLAGS="$FCFLAGS $CPU_FLAGS"
  fi

  if test -z "$FCFLAGS"; then
	echo ""
	echo "*********************************************************"
        echo "* WARNING: Don't know the best FCFLAGS for this system  *"
        echo "* Use  make FCFLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of FCFLAGS=-O3 will be used)    *"
	echo "*********************************************************"
	echo ""
        FCFLAGS="-O3"
  fi

  ACX_CHECK_FC_FLAGS(${FCFLAGS}, guessed_f90flags, , [
	echo ""
        echo "**********************************************************"
        echo "* WARNING: The guessed FCFLAGS don't seem to work with  *"
        echo "* your compiler.                                        *"
        echo "* Use  make FCFLAGS=..., or edit the top level Makefile *"
        echo "*********************************************************"
        echo ""
        FCFLAGS=""
  ])

fi
])

AC_DEFUN(ACX_PROG_F90_PREPFLAG,
[
AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine native compiler flags that allow us to use F90 suffix
# for preprocessed f90 source. 

if test "$ac_test_PREPFLAG" != "set"; then
  PREPFLAG=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$FC" = ifc -o "$FC" = ifort; then
                    PREPFLAG="-fpp1 "
                fi;;
  *aix*)  if test "$FC" = xlf90 -o "$FC" = f90 -o "$FC" = xlf95; then
                    PREPFLAG="-qsuffix=cpp=F90 "
                fi;;
  *darwin*)
	if test "$FC" = f90 -o "$FC" = xlf90 -o "$FC" = xlf95; then
        	PREPFLAG="-qsuffix=cpp=F90 "
        fi;;
  esac

  if test -z "$PREPFLAG"; then
        AC_MSG_WARN("Using empty PREPFLAG")
        PREPFLAG=""
  fi

  AC_MSG_CHECKING(to make sure F90 preprocessor flag works)
  AC_LANG_SAVE()
  AC_LANG(Fortran)
  ac_save_ext=$ac_ext
  ac_ext=F90
  ac_save_FCFLAGS_SRCEXT=$FCFLAGS_SRCEXT

  AS_IF([test "$PREPFLAG"], [FCFLAGS_SRCEXT="${PREPFLAG}"])
    _AC_COMPILE_IFELSE([
      AC_LANG_SOURCE([
program conftest
  integer :: i
  i = 1
end program conftest
])], [prepflagworks=1], [prepflagworks=0])

  FCFLAGS_SRCEXT=$ac_save_FCFLAGS_SRCEXT
  ac_ext=$ac_save_ext
  AC_LANG_RESTORE()

  if test "$prepflagworks" = 1; then
    AC_MSG_RESULT(yes)
    FCFLAGS_SRCEXT="${PREPFLAG}"
    AC_SUBST(FCFLAGS_SRCEXT)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Can't figure out working Fortran90 preprocessor flag])
  fi
fi
])


AC_DEFUN(ACX_PROG_F90_PREPDEFFLAG,
[
AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine native compiler flags that allow us to use F90 suffix
# for preprocessed f90 source with -D type defines

if test "$ac_test_PREPDEFFLAG" != "set"; then
  PREPDEFFLAG=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$FC" = ifc -o "$FC" = ifort; then
                    PREPDEFFLAG=" "
                fi;;
  *aix*)  if test "$FC" = xlf90 -o "$FC" = f90 -o "$FC" = xlf95; then
                    PREPDEFFLAG="-WF,"
                fi;;
  *darwin*)
	if test "$FC" = f90 -o "$FC" = xlf90 -o "$FC" = xlf95; then
        	PREPDEFFLAG="-WF,"
        fi;;
  esac

  if test -z "$PREPDEFFLAG"; then
        AC_MSG_WARN("Using empty PREPDEFFLAG")
        PREPDEFFLAG=" "
  fi

  AC_MSG_CHECKING(to make sure F90 preprocessor define flag works)
  AC_LANG_SAVE()
  AC_LANG(Fortran)
  ac_save_ext=$ac_ext
  ac_ext=F90
  ac_save_FCFLAGS=$FCFLAGS

  AS_IF([test "$PREPDEFFLAG"], [FCFLAGS="${FCFLAGS} ${PREPDEFFLAG}-DTEST"])
    _AC_COMPILE_IFELSE([
      AC_LANG_SOURCE([
program conftest
  integer :: i
#ifdef TEST
  i = 1
#else
  choke me
#endif
end program conftest
])], [prepdefflagworks=1], [prepdefflagworks=0])

  FCFLAGS=$ac_save_FCFLAGS 
  ac_ext=$ac_save_ext
  AC_LANG_RESTORE()

  if test "$prepdefflagworks" = 1; then
    AC_MSG_RESULT(yes)
    AC_SUBST(PREPDEFFLAG)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Can't figure out working Fortran90 preprocessor define flag])
  fi
fi
])

dnl check for the required MPI library
AC_DEFUN([ACX_MPI], [

# Set variables...
MPI_LIB_DIR="$MPI/lib"
MPI_INC_DIR="$MPI/include"
AC_SUBST([MPI_LIB_DIR])
AC_SUBST([MPI_INC_DIR])

AC_MSG_CHECKING([for mpi.h])
have_mpi_h=0
rm -f conftest*
echo '#include <mpi.h>' > conftest.cc
if ${CXX} -I${MPI_INC_DIR} -c conftest.cc 2>&1 ; then
        AC_MSG_RESULT(yes)
        have_mpi_h=1
else
	if test -s conftest.out ; then
		cat conftest.out >> config.log
	fi
        AC_MSG_RESULT(no! Check MPI include paths)
        USE_MPI="no"
fi
rm -f conftest*
if test "$have_mpi_h" = 1; then
    AC_DEFINE(HAVE_MPI_H, 1, [have mpi.h])
fi

AC_MSG_CHECKING([whether mpif.h is usable])
have_mpif_h=0
rm -f conftest*
cat >conftest.$ac_ext <<EOF
program main
   include 'mpif.h'
end
EOF
if $FC -I$MPI_INC_DIR -c $FCFLAGS $FCFLAGS_SRCEXT conftest.$ac_ext > conftest.out 2>&1 ; then
        AC_MSG_RESULT(yes)
        MPI_F90_INC="$MPI_INC_DIR"
        have_mpif_h=1
else
        if test -s conftest.out ; then 
                cat conftest.out >> config.log
        fi
        AC_MSG_RESULT([no! Check MPI include paths])
        USE_MPI="no"
fi
rm -f conftest*
AC_SUBST(MPI_F90_INC)
if test "$have_mpif_h" = 1; then
    AC_DEFINE(HAVE_MPIF_H, 1, [have mpif.h])
fi

AC_LANG_PUSH(C)
ac_save_LDFLAGS=$LDFLAGS
LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR} "

if test x = x"$MPI_LIB"; then
        AC_CHECK_LIB(mpich, MPI_Init, [MPI_LIB="-lmpich"])
fi
if test x = x"$MPI_LIB"; then
        AC_CHECK_LIB(mpi, MPI_Init, [MPI_LIB="-lmpi"])
fi
$as_unset ac_cv_lib_mpi_MPI_Init
if test x = x"$MPI_LIB"; then
        AC_CHECK_LIB(mpi, MPI_Init, [MPI_LIB="-lmpi -llam"], [],
                     "-llam")
fi
$as_unset ac_cv_lib_mpich_MPI_Init
if test x = x"$MPI_LIB"; then
        AC_CHECK_LIB(mpich, MPI_Init, [MPI_LIB="-lpmpich -lmpich"], [],
                     "-lpmpich")
fi

$as_unset ac_cv_lib_mpi_MPI_Init
if test x = x"$MPI_LIB"; then
AC_CHECK_LIB(mpi, MPI_Init, [MPI_LIB="-lmpi -llam -lpthread"],[
             AC_MSG_ERROR([Didn't find liblam, libmpi, or libmpich; check path for MPI package first...])
             USE_MPI="no"
             ],
             [-llam -lpthread])
fi

AC_SUBST(MPI_LIB)

AC_MSG_CHECKING([for MPI Fortran library])
MPI_F90_LIB=""
if test -f "$MPI_LIB_DIR/libfmpich.a" ; then
        MPI_F90_LIB="-lfmpich"
elif test -f "$MPI_LIB_DIR/liblamf77mpi.a" ; then
        MPI_F90_LIB="-llamf77mpi"
else
        dnl nothing special found, we'll assume that the C 
        dnl library is all we need
        MPI_F90_LIB="  "
fi
AC_MSG_RESULT([found $MPI_F90_LIB])
AC_SUBST(MPI_F90_LIB)
])dnl ACX_MPI

# ACX_CHECK_FFTW()
# ----------------
# This macro checks for fftw header files and libraries,
# including the possible prefixing with s or d to determine precision.
# Arg 1 is the fftw header/library name to check for, without 
# prefix or anything else (e.g. rfftw_mpi for real MPI transforms)
# Arg 2 is the size of the real variable used.
AC_DEFUN(ACX_CHECK_FFTW,
[
if test -z "$ac_fftw_firstname"; then

sizeof_real=$2
if test $sizeof_real = 8; then
  prec="double"
  fftwcheckprefix=d
else
  prec="single"
  fftwcheckprefix=s
fi

xfftwname=${fftwcheckprefix}$1

ok="no"
# check header doesn't work, since we must use mpicc to get includes, 
# we cant trust cpp.
AC_MSG_CHECKING([for $xfftwname.h])
AC_TRY_COMPILE([#include <$xfftwname.h>],,
[
fftwname=$xfftwname 
AC_MSG_RESULT(yes)
],
AC_MSG_RESULT(no))

# fftwname was set if we found a header

if test -n "$fftwname"; then
# we cannot run the code since an MPI program might not be allowed 
# on a login node of a supercomputer
AC_TRY_COMPILE([#include <$fftwname.h>], 
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $sizeof_real)]; ],
[
ok=yes
usedprefix=$fftwcheckprefix
],[ok=no])
fi

if test "$ok" != "yes"; then
  AC_MSG_CHECKING([for $1.h])
  AC_TRY_COMPILE([#include <$1.h>],,AC_MSG_RESULT(yes),
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cannot find any $prec precision $xfftwname.h or $1.h]
[Do you have $prec precision FFTW installed? If you are using packages,]
[note that you also need fftw-devel to use FFTW with OOPSE. You can find the ]
[software at www.fftw.org.]
[If you compiled FFTW yourself:                                        ]
[Note that the default FFTW setup is double precision.  If you want MPI support,]
[use --with-mpi. It is a good idea to install both single & double.] 
[If you have installed FFTW in a non-standard location, you should ]
[provide the correct paths in the CPPFLAGS and LDFLAGS environment ]
[variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.])
])
AC_TRY_COMPILE([#include <$1.h>],
[int _array_ [1 - 2 * !((sizeof(fftw_real)) == $sizeof_real)];],
[
usedprefix=""
fftwname=$1
],
[
AC_MSG_ERROR([Cannot find any $prec precision $xfftwname.h or $1.h]
[Do you have $prec precision FFTW installed? If you are using packages,]
[note that you also need fftw-devel to use FFTW with OOPSE. You can find the ]
[software at www.fftw.org.]
[If you compiled FFTW yourself:                                        ]
[Note that the default FFTW setup is double precision.  If you want MPI support,]
[use --with-mpi. It is a good idea to install both single & double.] 
[If you have installed FFTW in a non-standard location, you should ]
[provide the correct paths in the CPPFLAGS and LDFLAGS environment ]
[variables before running configure.]
[That is also necessary to do if your compiler doesn't search]
[/usr/local/include and /usr/local/lib by default.])])
fi

AC_CHECK_LIB($fftwname,main,,
AC_MSG_ERROR([Can't find a library to match the $fftwname header]))
ac_fftw_savedprefix=$usedprefix
ac_fftw_firstname=$fftwname

else

fftwname=${ac_fftw_savedprefix}$1
AC_MSG_CHECKING([for $fftwname.h])
AC_TRY_COMPILE(
[#include <$fftwname.h>],,
[AC_MSG_RESULT(yes)
LIBS="-l$fftwname $LIBS"
AC_TRY_LINK_FUNC([main],,,
AC_MSG_ERROR([Can't find a library to match the $fftwname header]))],
[
AC_MSG_RESULT(no)
AC_MSG_ERROR([Cant find $fftwname.h header. Make sure all your 
fftw prefixes match - we already use $ac_fftw_firstname.h])
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
[AC_MSG_CHECKING(if zlib is wanted)
AC_ARG_WITH(zlib,
[  --with-zlib=DIR root directory path of zlib installation [defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-zlib to disable zlib usage completely],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test -d "$withval"
  then
    ZLIB_HOME="$withval"
  else
    AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
  fi
else
  AC_MSG_RESULT(no)
fi])

ZLIB_HOME=/usr/local
if test ! -f "${ZLIB_HOME}/include/zlib.h"
then
        ZLIB_HOME=/usr
fi

#
# Locate zlib, if wanted
#
if test -n "${ZLIB_HOME}"
then
        ZLIB_OLD_LDFLAGS=$LDFLAGS
        ZLIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${ZLIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${ZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [zlib_cv_libz=yes], [zlib_cv_libz=no])
        AC_CHECK_HEADER(zlib.h, [zlib_cv_zlib_h=yes], [zlib_cv_zlib_h=no])
        AC_LANG_RESTORE
        if test "$zlib_cv_libz" = "yes" -a "$zlib_cv_zlib_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(z, inflateEnd)
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                LDFLAGS="$ZLIB_OLD_LDFLAGS"
                CPPFLAGS="$ZLIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid zlib installation with --with-zlib=DIR or disable zlib usage with --without-zlib)
        fi
fi

])
# AC_F90_MODULE_NAMES
# -------------------
#
# Figure out how the Fortran 90 compiler constructs module file names
# 
AC_DEFUN([AC_F90_MODULE_NAMES],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_CACHE_CHECK([for Fortran 90 module file names],
               ac_cv_f90_module_names,
[AC_LANG_PUSH(Fortran)
# carry out the test in a new directory, so that we don't miss anything
mkdir conftest
cd conftest
AC_COMPILE_IFELSE(
[MODULE Bar
END MODULE Bar],
ac_cv_f90_module_names=
[ac_file_list=*
for ac_file in $ac_file_list; do
   case $ac_file in
      # don't care for original source and object files 
      conftest.$ac_ext | conftest.$ac_objext | conftest.err )
          :
          ;;
      # look for new files derived from the file name 
      *conftest*)
          ac_pat=`echo $ac_file | sed s/conftest/%FILE%/`
	  _AC_LIST_MEMBER_IF($ac_pat, $ac_cv_f90_module_names,,
              ac_cv_f90_module_names="$ac_cv_f90_module_names $ac_pat")
          ;;
      # look for new files derived from the module name,
      # with different case translation schemes 
      *Bar*)
          ac_pat=`echo $ac_file | sed s/Bar/%Module%/`
	  _AC_LIST_MEMBER_IF($ac_pat, $ac_cv_f90_module_names,,
              ac_cv_f90_module_names="$ac_cv_f90_module_names $ac_pat")
          ;;
      *bar*)
          ac_pat=`echo $ac_file | sed s/bar/%module%/`
	  _AC_LIST_MEMBER_IF($ac_pat, $ac_cv_f90_module_names,,
              ac_cv_f90_module_names="$ac_cv_f90_module_names $ac_pat")
          ;;
      *BAR*)
          ac_pat=`echo $ac_file | sed s/BAR/%MODULE%/`
	  _AC_LIST_MEMBER_IF($ac_pat, $ac_cv_f90_module_names,,
              ac_cv_f90_module_names="$ac_cv_f90_module_names $ac_pat")
          ;;
       # Other files - we have no idea how they are generated
       *)
          AC_MSG_WARN([Bogus file found: $ac_file])
          ;;
   esac
done
if test "x$ac_cv_f90_module_names" = "x"; then
  AC_MSG_WARN([Couldn't determine module file names])
fi
],
[ac_cv_f90_module_names=
AC_MSG_WARN([Couldn't determine module file names])])
cd ..
# cleanup
rm -rf conftest
AC_LANG_POP()dnl
]) # AC_CACHE_CHECK

# We now generate a shell script that will help us to figure out the correct 
# module file names, using the value of ac_cv_f90_module_names

echo "Generating shell script modnam"

cat > scripts/modnam << EOF
#! /bin/sh
# This script is auto-generated by configure
#
usage="\\
Usage: \$[0] [[FILES]]

[[FILES]] are Fortran 90 source files. 
The output is a list of module file names that the Fortran 90 compiler 
generates when compiling [[FILES]]."

list=
empty=

if test \$[@%:@] -eq 0; then
   echo "\$usage"; exit 0
fi

while test \$[@%:@] != 0; do

  file=\$[1]
  shift

# strip suffix
  base=\`echo \$file | sed 's/[[.]][[^.]]*$//'\`

  test ! -f \$file && continue

# Look for module definitions and transform them to upper / lower case
  mods=\`cat \$file | sed '/^ *[[mM][oO][dD][uU][lL][eE]]/!d;s/^ *[[mM][oO][dD][uU][lL][eE]] *\([[A-Za-z_][A-Za-z0-9_]]*\).*\$/\1/'\`
  upper=\`echo \$mods | tr a-z A-Z\`
  lower=\`echo \$mods | tr A-Z a-z\`

# Here, the patterns for generating module file names were inserted by configure
  for trans in $ac_cv_f90_module_names; do

    pat=\`echo \$trans | sed 's/.*\(%.*%\).*/\1/'\`
    var=empty
    case \$pat in 
       %MODULE%) 
          var=upper ;;
       %Module%)
          var=mods ;;
       %module%) 
          var=lower ;;
       %FILE%)
          test -n "\$mods" && var=base ;;
    esac
    new=\`eval '(for i in \$'\$var '; do echo \$trans | sed s/\$pat/\$i/; done)'\`
    list="\$list \$new"
  done
done

echo \$list
# end of configure-generated script
EOF
chmod 755 scripts/modnam
]) # AC_F90_MODULE_NAMES
