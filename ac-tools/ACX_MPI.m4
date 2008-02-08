dnl ACX_MPI_GET_PROG(IF-FOUND, IF-NOT-FOUND)
AC_DEFUN([ACX_MPI_GET_PROG],[
dnl
dnl check for MPI C-compiler first
dnl
AC_CHECK_PROGS(acx_mpi_mpicc,[$CC openmpicc mpicc],[no])
AS_IF([test AS_VAR_GET([acx_mpi_mpicc]) != no],[
  AC_PATH_PROG(acx_mpi_mpicc_path,[$acx_mpi_mpicc],[no])])
dnl
dnl then check for MPI f90 compiler (f77 will not be sufficient)
dnl
AC_CHECK_PROGS(acx_mpi_mpif90,[$FC openmpif90 mpif90],[no])
AS_IF([test AS_VAR_GET([acx_mpi_mpif90]) != no],[
  AC_PATH_PROG(acx_mpi_mpif90_path,[$acx_mpi_mpif90],[no])])
dnl
dnl last, check for MPI c++ compiler
dnl
AC_CHECK_PROGS(acx_mpi_mpicxx,[$CXX openmpicxx openmpiCC openmpic++ mpicxx mpiCC mpic++],[no])
AS_IF([test AS_VAR_GET([acx_mpi_mpicxx]) != no],[
  AC_PATH_PROG(acx_mpi_mpicxx_path,[$acx_mpi_mpicxx],[no])])
dnl
dnl should we substitute MPI c++ compiler for C code?  Or vice versa?
dnl
AS_IF([test AS_VAR_GET([acx_mpi_mpicc]) = no],[
  AS_IF([test AS_VAR_GET([acx_mpi_mpicxx]) = no],[
    AS_VAR_SET([acx_mpi_compiler],[no])
  ],[
    AS_VAR_SET([acx_mpi_compiler],[AS_VAR_GET([acx_mpi_mpicxx])])
    AS_VAR_SET([acx_mpi_compiler_path],[AS_VAR_GET([acx_mpi_mpicxx_path])])])
],[
  AS_IF([test AS_VAR_GET([acx_mpi_mpicxx]) = no],[
    AS_VAR_SET([acx_mpi_compiler],[AS_VAR_GET([acx_mpi_mpicc])])
    AS_VAR_SET([acx_mpi_compiler_path],[AS_VAR_GET([acx_mpi_mpicc_path])])
  ],[
    AC_MSG_CHECKING([whether to use AS_VAR_GET([acx_mpi_mpicc]) or AS_VAR_GET([acx_mpi_mpicxx])])
    AC_LANG_PUSH(C++)
    AC_LANG_CONFTEST([AC_LANG_PROGRAM([
#include <mpi.h>
],[
MPI_Finalize();
])])
    AS_IF([_AC_DO_STDERR($acx_mpi_mpicxx_path -c conftest.$ac_ext >&AS_MESSAGE_LOG_FD) && {
         test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag" ||
         test ! -s conftest.err
       } && test -s conftest.$ac_objext],[
      AS_VAR_SET([acx_mpi_compiler],[AS_VAR_GET([acx_mpi_mpicxx])])
      AS_VAR_SET([acx_mpi_compiler_path],[AS_VAR_GET([acx_mpi_mpicxx_path])])
    ],[
      AS_VAR_SET([acx_mpi_compiler],[AS_VAR_GET([acx_mpi_mpicc])])
      AS_VAR_SET([acx_mpi_compiler_path],[AS_VAR_GET([acx_mpi_mpicc_path])])])
    AC_LANG_POP
    AC_MSG_RESULT([AS_VAR_GET([acx_mpi_compiler])])])])
dnl
dnl let's check to see if the fortran compiler works
dnl
AC_MSG_CHECKING([whether mpif.h is usable])
AC_LANG_PUSH(Fortran)
ac_save_ext=$ac_ext
ac_ext=F90
rm -f conftest*
cat >conftest.$ac_ext <<EOF
program main
   include 'mpif.h'
   call MPI_Init
end
EOF
if $acx_mpi_mpif90_path -c conftest.$ac_ext > conftest.out 2>&1 ; then
        AS_VAR_SET([acx_mpi_f90_compiler],[yes])
        AC_MSG_RESULT([yes])
else
        if test -s conftest.out ; then                     
                cat conftest.out >> config.log
        fi
        AS_VAR_SET([acx_mpi_f90_compiler],[no])
        AC_MSG_RESULT([no! Check MPI fortran include paths])
fi
rm -f conftest*
ac_ext=$ac_save_ext
AC_LANG_POP
AS_IF([test AS_VAR_GET([acx_mpi_compiler]) = no || test AS_VAR_GET([acx_mpi_f90_compiler]) = no],[$2],[$1])])# ACX_MPI_GET_PROG

dnl ACX_DETECT_MPI_IMPLEMENTATION(DEFAULT-IMPL, IF-FOUND, IF-NOT-FOUND)
AC_DEFUN([ACX_DETECT_MPI_IMPLEMENTATION],[
ACX_MPI_GET_PROG([dnl
AS_VAR_PUSHDEF([ac_var],[acx_cv_mpi_implementation])dnl
AC_CACHE_CHECK([for the MPI implementation flavor], [ac_var],
[dnl
AS_VAR_SET([ac_var],[no])
AS_IF([_AC_DO_STDERR(grep -q LAM $acx_mpi_compiler_path)],[
  AS_VAR_SET([ac_var],[lammpi])
],[
  AS_IF([_AC_DO_STDERR(grep -q showme $acx_mpi_compiler_path)],[
    AS_VAR_SET([ac_var],[openmpi])
  ],[
    AS_IF([_AC_DO_STDERR(grep -q MPICH $acx_mpi_compiler_path)],[
      AS_VAR_SET([ac_var],[mpich])
    ],[AS_VAR_SET([ac_var],[$1])])])])])
AS_VAR_POPDEF([ac_var])],[
AS_VAR_SET([acx_cv_mpi_implementation],[no])])dnl
case AS_VAR_GET([acx_cv_mpi_implementation]) in
(lammpi)
  AC_DEFINE(MPI_IS_LAMMPI,[1],[Indicate that MPI implementation is LAMMPI])
  ;;
(openmpi)
  AC_DEFINE(MPI_IS_OPENMPI,[1],[Indicate that MPI implementation is OPEN-MPI])
  ;;
(mpich)
  AC_DEFINE(MPI_IS_MPICH,[1],[Indicate that MPI implementation is MPICH])
  ;;
esac
AS_IF([test AS_VAR_GET([acx_cv_mpi_implementation]) = no],[$3],[$2])
])# ACX_DETECT_MPI_IMPLEMENTATION


AC_DEFUN([ACX_MPI_FILTER_CFLAGS],
  acx_mpi_tmp_mode="normal"
  acx_mpi_tmp_prefix=""
  AS_VAR_SET([$1],[""])
  for acx_mpi_tmp in $2; do
    case "$acx_mpi_tmp_mode" in
      (normal)
        case "$acx_mpi_tmp" in
          (-I|-D)
            acx_mpi_tmp_prefix="$acx_mpi_tmp"
            acx_mpi_tmp_mode="accept"
            ;;
          (-I*|-D*)
            AS_VAR_SET([$1],["AS_VAR_GET([$1]) $acx_mpi_tmp"])
            ;;
          (-L|-l)
            acx_mpi_tmp_mode="skip"
            ;;
          (*)
            ;;
        esac
        ;;
      (accept)
        AS_VAR_SET([$1],["AS_VAR_GET([$1]) $acx_mpi_tmp_prefix $acx_mpi_tmp"])
        ;;
      (skip)
        ;;
    esac
  done
)


dnl ACX_MPI(DEFAULT-IMPL, IF-FOUND, IF-NOT-FOUND)
AC_DEFUN([ACX_MPI],[
AC_ARG_WITH(mpi,
  AS_HELP_STRING([--with-mpi=auto|mpich|lam|openmpi|no],
                 [Indicates which kind of MPI implementation to use @<:@default=auto@:>@]),
  [],
  [with_mpi="auto"])
if test "x$with_mpi" != "xno"; then
ACX_DETECT_MPI_IMPLEMENTATION(AS_VAR_GET([with_mpi]),[
dnl
AS_VAR_PUSHDEF([ac_var],[acx_mpi_cflags])
AC_CACHE_CHECK([how to compile MPI-C code],[ac_var],[
case AS_VAR_GET([acx_cv_mpi_implementation]) in
  (lammpi)
    cfo="--showme:compile"
    ;;
  (openmpi)
    cfo="--showme:compile"
    ;;
  (mpich)
    case "$acx_mpi_compiler_path" in
      (*mpiCC)
        sc_cv_cn="-CC="
	;;
      (*mpicxx)
        sc_cv_cn="-cxx="
        ;;
      (*mpicc)
        sc_cv_cn="-cc="
	;;
      (*)
        sc_cv_cn=""
	;;
    esac
    cfo="-compile_info $sc_cv_cn"
    ;;
esac
_AS_ECHO_LOG([mpi_pre_cflags="`$acx_mpi_compiler_path $cfo || echo ' no'`"])
mpi_pre_cflags="`$acx_mpi_compiler_path $cfo 2>conftest.er1 || echo ' no'`"
grep -v '^ *+' conftest.er1 >conftest.err
rm -f conftest.er1
cat conftest.err >&AS_MESSAGE_LOG_FD
_AS_ECHO_LOG([mpi_pre_cflags = $mpi_pre_cflags])
case "$mpi_pre_cflags" in
  (*no)
    AS_VAR_SET([ac_var],[no])
    ac_var="no"
    ;;
  (*)
    ACX_MPI_FILTER_CFLAGS([ac_var],[$mpi_pre_cflags])
    ;;
esac
])
AS_VAR_POPDEF([ac_var])
dnl
AS_VAR_PUSHDEF([ac_var],[acx_mpi_libs])
AC_CACHE_CHECK([how to link MPI-C code],[ac_var],[
case AS_VAR_GET([acx_cv_mpi_implementation]) in
  (lammpi)
    lfo="--showme:compile --showme:link"
    ;;
  (openmpi)
    lfo="--showme:link"
    ;;
  (mpich)
    case "$acx_mpi_compiler_path" in
      (*mpiCC)
        sc_cv_cn="-CC="
	;;
      (*mpicxx)
        sc_cv_cn="-cxx="
        ;;
      (*mpicc)
        sc_cv_cn="-cc="
	;;
      (*)
        sc_cv_cn=""
	;;
    esac
    lfo="-link_info $sc_cv_cn"
    ;;
esac
dnl
_AS_ECHO_LOG([mpi_pre_libs="`$acx_mpi_compiler_path $lfo || echo ' no'`"])
mpi_pre_libs="`$acx_mpi_compiler_path $lfo 2>conftest.er1 || echo ' no'`"
grep -v '^ *+' conftest.er1 >conftest.err
rm -f conftest.er1
cat conftest.err >&AS_MESSAGE_LOG_FD
_AS_ECHO_LOG([mpi_pre_libs = $mpi_pre_libs])
case "$mpi_pre_libs" in
  (*no)
    AS_VAR_SET([ac_var],[no])
    ;;
  (*)
    AS_VAR_SET([ac_var],["$mpi_pre_libs"])
    ;;
esac
])
AS_VAR_POPDEF([ac_var])
dnl
AS_IF([test AS_VAR_GET([acx_mpi_mpif90_path]) != no],[
  AS_VAR_PUSHDEF([ac_var],[acx_mpi90_libs])
  AC_CACHE_CHECK([how to link MPI-Fortran code],[ac_var],[
    _AS_ECHO_LOG([mpi_pre_libs="`$acx_mpi_mpif90_path $lfo || echo " no"`"])
    mpi_pre_libs="`$acx_mpi_mpif90_path $lfo 2>conftest.er1 || echo " no"`"
    grep -v '^ *+' conftest.er1 >conftest.err
    rm -f conftest.er1
    cat conftest.err >&AS_MESSAGE_LOG_FD
    _AS_ECHO_LOG([mpi_pre_libs = $mpi_pre_libs])
    case "$mpi_pre_libs" in
      (*no)
        AS_VAR_SET([ac_var],[no])
        ;;
      (*)
        AS_VAR_SET([ac_var],["$mpi_pre_libs"])
        ;;
    esac
  ])
  AS_VAR_POPDEF([ac_var])
])
dnl
AS_IF([test "AS_VAR_GET([acx_mpi_cflags])" = no -o "AS_VAR_GET([acx_mpi_libs])" = no],[
AS_VAR_SET([acx_mpi_cflags],[""])
AS_VAR_SET([acx_mpi_libs],[""])
$2],[
AC_SUBST(MPI_CPPFLAGS,AS_VAR_GET([acx_mpi_cflags]))
AC_SUBST(MPI_LIBS,AS_VAR_GET([acx_mpi_libs]))
AS_IF([test "AS_VAR_GET([acx_mpi90_libs])" != no],[AC_SUBST([MPI90_LIBS],[AS_VAR_GET([acx_mpi90_libs])])])
$1])
],[
AS_VAR_SET([acx_mpi_cflags],[""])
AS_VAR_SET([acx_mpi_libs],[""])
$2])
fi
])# ACX_MPI

dnl ACX_MPI2([IF-FOUND],[IF-NOT-FOUND])
AC_DEFUN([ACX_MPI2],[dnl
AC_REQUIRE([ACX_MPI])
AS_IF([test AS_VAR_GET([acx_cv_mpi_implementation]) = no],[$2],[dnl
sc_cv_oldCPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
AC_CHECK_DECLS([MPI_Comm_spawn],[
   sc_cv_oldLIBS="$LIBS"
   LIBS="$LIBS $MPI_LIBS"
   AC_CHECK_FUNCS([MPI_Init_thread MPI_Comm_disconnect MPI_Comm_get_parent MPI_Finalized MPI_Comm_spawn],
                  [acx_cv_mpi_is_mpi2=yes],
                  [acx_cv_mpi_is_mpi2=no])
   LIBS="$sc_cv_oldLIBS"
],[acx_cv_mpi_is_mpi2=no],
     [
#include <mpi.h>
     ])
CPPFLAGS="$sc_cv_oldCPPFLAGS"
AS_IF([test AS_VAR_GET([acx_cv_mpi_is_mpi2]) = yes],[$1],[$2])
])])

dnl ACX_MPI_SEEK([IF-OK],[IF-NOT-OK])
AC_DEFUN([ACX_MPI_SEEK],[dnl
AC_REQUIRE([ACX_MPI])
AS_IF([test AS_VAR_GET([acx_cv_mpi_implementation]) = no],[$2],[dnl
AS_VAR_PUSHDEF([ac_var],[acx_mpi_seek_ok])
AC_CACHE_CHECK([whether MPI accepts stdio SEEK_* to be defined],[ac_var],[dnl
acx_mpi_old_cppflags="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
AC_PREPROC_IFELSE(AC_LANG_SOURCE([[
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2

#include <mpi.h>
]]),[AS_VAR_SET([ac_var],[yes])],[AS_VAR_SET([ac_var],[no])])
CPPFLAGS="$acx_mpi_old_cppflags"
])
AS_IF([test AS_VAR_GET([ac_var]) = yes],[
  AC_DEFINE([MPI_ACCEPTS_SEEK],[1],[Indicates whether MPI accepts stdio SEEK_* to be defined])
  $1
],[
  AC_DEFINE([MPI_REFUSES_SEEK],[1],[Indicates whether MPI accepts stdio SEEK_* to be defined])
  $2])
AS_VAR_POPDEF([ac_var])
])])
