dnl
dnl AC_CHECK_MODSUFFIX
dnl
AC_DEFUN([AC_CHECK_MODSUFFIX],[

AC_MSG_CHECKING(for module suffix)
rm -f conftest*
# Intel ifc compiler generates files by the name of work.pc and work.pcl (!)
rm -f work*
cat >conftest.f90 <<EOF
        module conftest
        integer n
        parameter (n=1)
        end module conftest
EOF
# SGI and absoft compilers generates module name in upper case!
testname="conftest"
modcase="lower"
if $F90 -c conftest.f90 > conftest.out 2>&1 ; then
    MOD=`ls conftest* | grep -v conftest.f | grep -v conftest.o`
    MOD=`echo "$MOD" | sed -e 's/conftest\.//g'`
    if test -z "$MOD" ; then
        MOD=`ls CONFTEST* 2>/dev/null \
                | grep -v CONFTEST.f | grep -v CONFTEST.o`
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
cat >conftest.f90 <<EOF
        module conftest
        integer n
        parameter (n=1)
        end module conftest
EOF
# SGI and absoft compilers generates module name in upper case!
testname="conftest"
if $F90 -c conftest.f90 > conftest.out 2>&1 ; then
   mod=`ls CONFTEST* 2>/dev/null | grep -v CONFTEST.f | grep -v CONFTEST.o`
   mod=`echo "$mod" | sed -e 's/CONFTEST\.//g'`
   if test -n "$mod" ; then
      testname="CONFTEST"
   fi
   madedir=0
   if test ! -d conf ; then mkdir conf ; madedir=1; fi
   cp $testname.$MOD conf
   rm -f conftest* CONFTEST*
   cat >conftest1.f90 <<EOF
        program main
        use conftest
        print *, n
        end
EOF
   F90_WORK_FILES_ARG=""
   F90MODINCSPEC=""
   if $F90 -c -Iconf conftest1.f90 > conftest.out 2>&1 ; then
       MODDIRFLAG="-I"
       F90MODINCSPEC="-I<dir>"
       AC_MSG_RESULT(-I)
   elif $F90 -c -Mconf conftest1.f90 >> conftest.out 2>&1 ; then
       MODDIRFLAG="-M"
       F90MODINCSPEC="-M<dir>"
       AC_MSG_RESULT(-M)
   elif $F90 -c -pconf conftest1.f90 >> conftest.out 2>&1 ; then
       MODDIRFLAG="-p"
       F90MODINCSPEC="-p<dir>"
       AC_MSG_RESULT(-p)
   elif test -s work.pc ; then
        cp work.pc conf/mpimod.pc
        echo "mpimod.pc" > conf/mpimod.pcl
        echo "`pwd`/conf/mpimod.pc" >> conf/mpimod.pcl
        if $F90 -c -cl,conf/mpimod.pcl conftest1.f >>conftest.out 2>&1 ; then
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


dnl
dnl AC_CHECK_MPI_F90MOD
dnl
AC_DEFUN([AC_CHECK_MPI_F90MOD],[

  AC_ARG_WITH(mpi_f90_mods, [  --with-mpi_f90_mods=<dir>
                          Location where MPI f90 modules are installed ],
	                   mpi_f90_mods="$withval", 
	                   mpi_f90_mods="/usr/local/include/f90choice")

  AC_MSG_CHECKING(for MPI F90 modules)
  AC_LANG_SAVE()
  AC_LANG([Fortran 90])
  ac_save_F90FLAGS=$F90FLAGS

  AS_IF([test "$mpi_f90_mods"], [F90FLAGS="${F90FLAGS} ${MODDIRFLAG}${mpi_f90_mods}"])
    _AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([
use mpi
integer :: ierr
call MPI_Init(ierr)
])], [HAVE_MPI_MOD=1], [HAVE_MPI_MOD=0])
 
  F90FLAGS=$ac_save_F90FLAGS 
  AC_LANG_RESTORE()

  if test "$HAVE_MPI_MOD" = 1; then
    AC_MSG_RESULT(yes)
    AC_DEFINE(MPI_MOD, 1, [have mpi module])
    MPI_F90_MODS="${mpi_f90_mods}"
    AC_SUBST(MPI_F90_MODS)        
    # The library name:
    if test -z "$MPI_LIB" ; then 
       MPI_LIBNAME=-lmpich
    else 
       MPI_LIBNAME="$MPI_LIB" 
    fi
    if test -z "$MPIMODLIBNAME" ; then
       MPIMODLIBNAME="${MPI_LIBNAME}f90"
    fi
    AC_SUBST(MPIMODLIBNAME)
  else
    AC_MSG_RESULT(no)
    AC_MSG_WARN([Couldn't locate MPI F90 Modules])
  fi

])




dnl
dnl AM_PATH_SPRNG
dnl
AC_DEFUN([AM_PATH_SPRNG],[

  AC_ARG_WITH(sprng_prefix, [  --with-sprng_prefix=PREFIX
                          Prefix where SPRNG is installed ],
                           sprng_prefix="$withval",
                           sprng_prefix="/usr/local")
  AC_ARG_WITH(sprng-libdir, [  --with-sprng-libdir=PREFIX  SPRNG library directory],
                           sprng_libdir="$withval",
                           sprng_libdir="/usr/local/lib")
  AC_ARG_WITH(sprng-include, [  --with-sprng-include=PREFIX
                          SPRNG include directory],
                           sprng_include="$withval",
                           sprng_include="/usr/local/include/sprng")

  if test x$sprng_libdir = x ; then
    sprng_libdir=${sprng_prefix}/lib
  fi

  if test x$sprng_include = x ; then
    sprng_include=${sprng_prefix}/include
  fi

  AC_MSG_CHECKING(for SPRNG include files in $sprng_include)
  if test -f ${sprng_include}/sprng.h; then
    have_sprngincl=yes
    AC_MSG_RESULT(yes)
  else
    have_sprngincl=no
    AC_MSG_RESULT(no)
    AC_MSG_ERROR(Could not locate the SPRNG include files)
  fi

  AC_MSG_CHECKING(for SPRNG libraries in $sprng_libdir)
  if test -f ${sprng_libdir}/libsprng.a; then
    have_sprnglib=yes
    AC_MSG_RESULT(yes)
  else
    have_sprnglib=no
    AC_MSG_RESULT(no)
    AC_MSG_ERROR(Could not locate the SPRNG libraries)
  fi

  AC_LANG_SAVE()
  AC_LANG([C])
  ac_save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="${CPPFLAGS} -I${sprng_include}"
  ac_save_LDFLAGS=$LDFLAGS
  LDFLAGS="${LDFLAGS} -L${sprng_libdir} -lsprng"
  AC_CHECK_HEADER(sprng.h, [
    AC_CHECK_LIB(sprng, 
                 init_rng, 
                    [SPRNG_LIBDIR="${sprng_libdir}"
                     SPRNG_LIB="-lsprng" SPRNG_INC="-I${sprng_include}"
                     HAVE_SPRNG="yes"])])
  CPPFLAGS=$ac_save_CPPFLAGS 
  LDFLAGS=$ac_save_LDFLAGS 
  AC_LANG_RESTORE()
  
  if test x_$HAVE_SPRNG != x_yes; then
        AC_MSG_ERROR(Can't build with SPRNG)
  fi
  AC_SUBST(SPRNG_LIBDIR)
  AC_SUBST(SPRNG_LIB)
  AC_SUBST(SPRNG_INC)
  AC_SUBST(HAVE_SPRNG)
])

dnl
dnl CHECK_MKL
dnl
AC_DEFUN([CHECK_MKL],
[AC_MSG_CHECKING(if MKL is wanted)
AC_ARG_WITH(mkl,
	      [  --with-mkl              Do we want MKL [will check /usr/local/intel/mkl61 /opt/intel/mkl61]],
[   AC_MSG_RESULT(yes)
    for dir in $withval /usr/local/intel/mkl61 /opt/intel/mkl61; do
        mkldir="$dir"
        if test -f "$dir/include/mkl.h"; then
            found_mkl="yes";
            break;
        fi
    done
    if test x_$found_mkl != x_yes; then
        AC_MSG_ERROR(Cannot find MKL includes)
    else
        printf "MKL includes found in $mkldir/include\n";
    fi

  AC_LANG_SAVE()
  AC_LANG([C])
  ac_save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="${CPPFLAGS} -I${mkldir}/include"
  ac_save_LDFLAGS=$LDFLAGS
  LDFLAGS="${LDFLAGS} -L${mkldir}/lib/32 -lmkl -lvml -lguide"
  AC_CHECK_HEADER(mkl.h, [
    AC_CHECK_LIB(mkl, 
                 vslNewStream, 
                    [MKL_LIBDIR="${mkldir}/lib/32", 
                     MKL_LIB="-lmkl -lvml -lguide",
                     HAVE_MKL="yes"]) 
    ], [MKL_INC="-I${mkldir}/include"])
  CPPFLAGS=$ac_save_CPPFLAGS 
  LDFLAGS=$ac_save_LDFLAGS 
  AC_LANG_RESTORE()
  
  if test x_$HAVE_MKL != x_yes; then
        AC_MSG_ERROR(Can't build with MKL)
  fi
  AC_SUBST(MKL_LIBDIR)
  AC_SUBST(MKL_LIB)
  AC_SUBST(MKL_INC)
  AC_SUBST(HAVE_MKL)
],
[
    AC_MSG_RESULT(no)
])
])
dnl


AC_DEFUN(ACX_CHECK_CC_FLAGS,
[
AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK(whether ${CC-cc} accepts $1, ac_$2,
[echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} $1 -c conftest.c 2>&1`"; then
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
AC_CACHE_CHECK(whether ${CXX-c++} accepts $1, ac_$2,
[echo 'void f(){}' > conftest.cpp
if test -z "`${CXX-c++} $1 -c conftest.cpp 2>&1`"; then
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
dnl ACX_CHECK_F90_FLAGS()
dnl
dnl     Check for optimizer flags the Fortran compiler can use.
dnl
AC_DEFUN(ACX_CHECK_F90_FLAGS,
[
AC_CACHE_CHECK(whether ${F90-f90} accepts $1, ac_$2,
[
AC_LANG_SAVE
AC_LANG([Fortran 90])
echo 'program main' > conftest.$ac_ext
echo 'end program main' >> conftest.$ac_ext
ac_compile='${F90-f90} -c $1 $F90FLAGS conftest.$ac_ext 1>&AC_FD_CC'
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
  #ifndef __INTEL_COMPILER
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
  #ifndef __INTEL_COMPILER
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
                    CFLAGS="-native -fast -xO5 -dalign"
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
                    CXXFLAGS="-native -fast -xO5 -dalign"
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
	if test "$CXX" = xlc++; then
        	CXXFLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto -qaltivec"
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

AC_DEFUN(ACX_PROG_F90_MAXOPT,
[
AC_REQUIRE([AC_PROG_F90])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine "good" native compiler flags if none specified on command
# line

if test x"$F90FLAGS" = x ; then
  F90FLAGS=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$F90" = ifc -o "$F90" = ifort; then
                    F90FLAGS="-O"
                fi;;
   rs6000*-aix*)  if test "$F90" = xlf90 -o "$F90" = f90; then
                    F90FLAGS="-O3 -qarch=pwrx -qtune=pwrx -qansialias -w"
                fi;;
   powerpc*-aix*)
	if test "$F90" = f90 -o "$F90" = xlf90; then
        	F90FLAGS="-O3 -qarch=ppc -qansialias -w"
		echo "*******************************************************"
		echo "*  You have AIX on an unknown powerpc system.  It is  *"
		echo "*  recommended that you use                           *"
		echo "*                                                     *"
		echo "*  F90FLAGS=-O3 -qarch=ppc -qtune=xxx -qansialias -w  *"
		echo "*                                 ^^^                 *"
		echo "*  where xxx is 601, 603, 604, or whatever kind of    *"
                echo "*  PowerPC CPU you have.   For more info, man xlf.    *"
		echo "*******************************************************"
        fi;;
   *darwin*)
	if test "$F90" = f90 -o "$F90" = xlf90 -o "$F90" = xlf95; then
        	F90FLAGS="-qthreaded -O -qtune=auto -qarch=auto -qunroll=auto"
        fi;;
  esac

  if test -n "$CPU_FLAGS"; then
        F90FLAGS="$F90FLAGS $CPU_FLAGS"
  fi

  if test -z "$F90FLAGS"; then
	echo ""
	echo "**********************************************************"
        echo "* WARNING: Don't know the best F90FLAGS for this system  *"
        echo "* Use  make F90FLAGS=..., or edit the top level Makefile *"
	echo "* (otherwise, a default of F90FLAGS=-O3 will be used)    *"
	echo "**********************************************************"
	echo ""
        F90FLAGS="-O3"
  fi

  ACX_CHECK_F90_FLAGS(${F90FLAGS}, guessed_f90flags, , [
	echo ""
        echo "**********************************************************"
        echo "* WARNING: The guessed F90FLAGS don't seem to work with  *"
        echo "* your compiler.                                         *"
        echo "* Use  make F90FLAGS=..., or edit the top level Makefile *"
        echo "**********************************************************"
        echo ""
        F90FLAGS=""
  ])

fi
])

AC_DEFUN(ACX_PROG_F90_PREPFLAG,
[
AC_REQUIRE([AC_PROG_F90])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine native compiler flags that allow us to use F90 suffix
# for preprocessed f90 source. 

if test "$ac_test_PREPFLAG" != "set"; then
  PREPFLAG=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$F90" = ifc -o "$F90" = ifort; then
                    PREPFLAG="-fpp1 "
                fi;;
  *aix*)  if test "$F90" = xlf90 -o "$F90" = f90; then
                    PREPFLAG="-qsuffix=cpp=F90 "
                fi;;
  *darwin*)
	if test "$F90" = f90 -o "$F90" = xlf90; then
        	PREPFLAG="-qsuffix=cpp=F90 "
        fi;;
  esac

  if test -z "$PREPFLAG"; then
        AC_MSG_WARN("Using empty PREPFLAG")
        PREPFLAG=""
  fi

  AC_MSG_CHECKING(to make sure F90 preprocessor flag works)
  AC_LANG_SAVE()
  AC_LANG([Fortran 90])
  ac_save_ext=$ac_ext
  ac_ext=F90
  ac_save_F90FLAGS=$F90FLAGS

  AS_IF([test "$PREPFLAG"], [F90FLAGS="${F90FLAGS} ${PREPFLAG}-DTEST"])
    _AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([
  integer :: i
  i = 1
])], [prepflagworks=1], [prepflagworks=0])

  F90FLAGS=$ac_save_F90FLAGS 
  ac_ext=$ac_save_ext
  AC_LANG_RESTORE()

  if test "$prepflagworks" = 1; then
    AC_MSG_RESULT(yes)
    AC_SUBST(PREPFLAG)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Can't figure out working Fortran90 preprocessor flag])
  fi
fi
])


AC_DEFUN(ACX_PROG_F90_PREPDEFFLAG,
[
AC_REQUIRE([AC_PROG_F90])
AC_REQUIRE([AC_CANONICAL_HOST])

# Try to determine native compiler flags that allow us to use F90 suffix
# for preprocessed f90 source with -D type defines

if test "$ac_test_PREPDEFFLAG" != "set"; then
  PREPDEFFLAG=""
  case "${host_cpu}-${host_os}" in

  *linux*) if test "$F90" = ifc -o "$F90" = ifort; then
                    PREPDEFFLAG=" "
                fi;;
  *aix*)  if test "$F90" = xlf90 -o "$F90" = f90; then
                    PREPDEFFLAG="-WF,"
                fi;;
  *darwin*)
	if test "$F90" = f90 -o "$F90" = xlf90; then
        	PREPDEFFLAG="-WF,"
        fi;;
  esac

  if test -z "$PREPDEFFLAG"; then
        AC_MSG_WARN("Using empty PREPDEFFLAG")
        PREPDEFFLAG=" "
  fi

  AC_MSG_CHECKING(to make sure F90 preprocessor define flag works)
  AC_LANG_SAVE()
  AC_LANG([Fortran 90])
  ac_save_ext=$ac_ext
  ac_ext=F90
  ac_save_F90FLAGS=$F90FLAGS

  AS_IF([test "$PREPDEFFLAG"], [F90FLAGS="${F90FLAGS} ${PREPFLAG} ${PREPDEFFLAG}-DTEST"])
    _AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([
  integer :: i
#ifdef TEST
  i = 1
#else
  choke me
#endif
])], [prepdefflagworks=1], [prepdefflagworks=0])

  F90FLAGS=$ac_save_F90FLAGS 
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
echo '#include <mpi.h>' > conftest.cc
if test -z "`${CXX} -I${MPI_INC_DIR} -c conftest.cc 2>&1`"; then
        AC_MSG_RESULT(yes)
        have_mpi_h=1
else
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
cat >conftest.f90 <<EOF
program main
   include 'mpif.h'
end
EOF
if $F90 -I$MPI_INC_DIR -c conftest.f90 > conftest.out 2>&1 ; then
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


dnl check for the required SPRNG library
AC_DEFUN([ACX_SPRNG], [

# Set variables...
SPRNG_LIB_DIR="$SPRNG/lib"
SPRNG_INC_DIR="$SPRNG/include"
AC_SUBST([SPRNG_LIB_DIR])
AC_SUBST([SPRNG_INC_DIR])

AC_MSG_CHECKING([for sprng.h])
have_sprng_h=0
echo '#include <sprng.h>' > conftest.cc
if test -z "`${CXX} -I${SPRNG_INC_DIR} -c conftest.cc 2>&1`"; then
        AC_MSG_RESULT(yes)
        have_sprng_h=1
else
        AC_MSG_RESULT([no! Check SPRNG include path!])
        USE_SPRNG="no"
fi
rm -f conftest*
if test "$have_sprng_h" = 1; then
    AC_DEFINE(HAVE_SPRNG_H, 1, [have sprng.h])
fi

AC_LANG_PUSH(C)
ac_save_LDFLAGS=$LDFLAGS
LDFLAGS="${LDFLAGS} -L${SPRNG_LIB_DIR} "

AC_CHECK_LIB(sprng, init_rng, [SPRNG_LIB="-lsprng"], [
             AC_MSG_ERROR([Didn't find libsprng; check path for SPRNG package first...])
             USE_SPRNG="no"
             ])

if test "$USE_SPRNG" = "no"; then
  AC_MSG_ERROR(No working SPRNG library found)
fi
AC_SUBST(SPRNG_LIB)
])dnl ACX_SPRNG

AC_DEFUN([adl_FUNC_GETOPT_LONG],
 [AC_PREREQ(2.49)dnl
  # clean out junk possibly left behind by a previous configuration
  rm -f src/utils/getopt.h
  # Check for getopt_long support
  AC_CHECK_HEADERS([getopt.h])
  AC_CHECK_FUNCS([getopt_long],,
   [# FreeBSD has a gnugetopt library for this
    AC_CHECK_LIB([gnugetopt],[getopt_long],[AC_DEFINE([HAVE_GETOPT_LONG])],
     [# use the GNU replacement
      AC_LIBOBJ(getopt)
      AC_LIBOBJ(getopt1)
      AC_CONFIG_LINKS([src/utils/getopt.h:src/utils/gnugetopt.h])])])])


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

