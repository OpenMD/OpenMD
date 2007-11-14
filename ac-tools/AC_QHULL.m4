AC_DEFUN(AC_CHECK_QHULL_VERSION,
[AC_MSG_CHECKING([for qh_qhull in -lqhull with qh_version])
AC_CACHE_VAL(ac_cv_lib_qhull_version,
[changequote(, )dnl
cat > conftest.c <<EOF
#include <stdio.h>
char qh_version[] = "version";
char qh_qhull();
int
main(argc, argv)
  int argc;
  char *argv[];
{
  qh_qhull();
  return 0;
}
EOF
changequote([, ])dnl
ac_try="${CC-cc} $CFLAGS $CPPFLAGS $LDFLAGS conftest.c -o conftest -lqhull $LIBS"
if AC_TRY_EVAL(ac_try) && test -s conftest ; then
    ac_cv_lib_qhull_version=yes
else
    ac_cv_lib_qhull_version=no
fi
rm -f conftest.c conftest.o conftest
])dnl
if test "$ac_cv_lib_qhull_version" = "yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$1], , , [$1])
else
  AC_MSG_RESULT(no)
  ifelse([$2], , , [$2])
fi
])

dnl Do we have libqhull?
AC_SUBST(have_qhull)
AC_SUBST(need_qhull_version)
AC_CHECK_HEADER(qhull/qhull_a.h, have_qhull=yes, have_qhull=no)
if test $have_qhull = yes ; then
   OF_CHECK_LIB(qhull, qh_qhull, have_qhull=yes, have_qhull=no)
    if test $have_qhull = yes ; then
        need_qhull_version=no
    else
	AC_CHECK_QHULL_VERSION(have_qhull=yes, have_qhull=no)
	need_qhull_version=yes
    fi
fi
if test $have_qhull = yes ; then
    QHULLSTATUS=yes
else
    QHULLSTATUS="Qhull not found --- see docs/README.geometry"
fi
