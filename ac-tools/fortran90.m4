# This file is part of Autoconf.                       -*- Autoconf -*-
# Programming languages support.
# Copyright 2000, 2001
# Free Software Foundation, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
#
# Written by Akim Demaille, Christian Marquardt, Martin Wilks (and probably
# many others). 


# Table of Contents:
#
# 1. Language selection
#    and routines to produce programs in a given language.
#  a. Fortran 77 (to be moved from aclang.m4)
#  b. Fortran 90
#  c. Fortran 95
#
# 2. Producing programs in a given language.
#  a. Fortran 77 (to be moved from aclang.m4)
#  b. Fortran 90
#  c. Fortran 95
#
# 3. Looking for a compiler
#    And possibly the associated preprocessor.
#  a. Fortran 77 (to be moved from aclang.m4)
#  b. Fortran 90
#  c. Fortran 95
#
# 4. Compilers' characteristics.
#  a. Fortran 77 (to be moved from aclang.m4)
#  b. Fortran 90
#  c. Fortran 95

# _AC_LIST_MEMBER_IF(ELEMENT, LIST, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ---------------------------------------------------------------------------
#
# Processing the elements of a list is tedious in shell programming,
# as lists tend to be implemented as space delimited strings.
#
# This macro searches LIST for ELEMENT, and executes ACTION-IF-FOUND
# if ELEMENT is a member of LIST, otherwise it executes
# ACTION-IF-NOT-FOUND.
AC_DEFUN([_AC_LIST_MEMBER_IF],
[dnl Do some sanity checking of the arguments.
m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  ac_exists=false
  for ac_i in $2; do
    if test x"$1" = x"$ac_i"; then
      ac_exists=true
      break
    fi
  done

  AS_IF([test x"$ac_exists" = xtrue], [$3], [$4])[]dnl
])# _AC_LIST_MEMBER_IF



# _AC_LINKER_OPTION(LINKER-OPTIONS, SHELL-VARIABLE)
# -------------------------------------------------
#
# Specifying options to the compiler (whether it be the C, C++ or
# Fortran 77 compiler) that are meant for the linker is compiler
# dependent.  This macro lets you give options to the compiler that
# are meant for the linker in a portable, compiler-independent way.
#
# This macro take two arguments, a list of linker options that the
# compiler should pass to the linker (LINKER-OPTIONS) and the name of
# a shell variable (SHELL-VARIABLE).  The list of linker options are
# appended to the shell variable in a compiler-dependent way.
#
# For example, if the selected language is C, then this:
#
#   _AC_LINKER_OPTION([-R /usr/local/lib/foo], foo_LDFLAGS)
#
# will expand into this if the selected C compiler is gcc:
#
#   foo_LDFLAGS="-Xlinker -R -Xlinker /usr/local/lib/foo"
#
# otherwise, it will expand into this:
#
#   foo_LDFLAGS"-R /usr/local/lib/foo"
#
# You are encouraged to add support for compilers that this macro
# doesn't currently support.
# FIXME: Get rid of this macro.
AC_DEFUN([_AC_LINKER_OPTION],
[if test "$ac_compiler_gnu" = yes; then
  for ac_link_opt in $1; do
    $2="[$]$2 -Xlinker $ac_link_opt"
  done
else
  $2="[$]$2 $1"
fi[]dnl
])# _AC_LINKER_OPTION



## ----------------------- ##
## 1. Language selection.  ##
## ----------------------- ##

# ----------------------------- #
# 1b. The Fortran 90 language.  #
# ----------------------------- #

# AC_LANG(Fortran 90)
# -------------------
m4_define([AC_LANG(Fortran 90)],
[ac_ext=f90
ac_compile='$F90 -c $F90FLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$F90 -o conftest$ac_exeext $F90FLAGS $LD90FLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_f90_compiler_gnu
])


# _AC_LANG_ABBREV(Fortran 90)
# ---------------------------
m4_define([_AC_LANG_ABBREV(Fortran 90)], [f90])


# ----------------------------- #
# 1c. The Fortran 95 language.  #
# ----------------------------- #

# AC_LANG(Fortran 95)
# -------------------
m4_define([AC_LANG(Fortran 95)],
[ac_ext=f95
ac_compile='$F95 -c $F95FLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$F95 -o conftest$ac_exeext $F95FLAGS $LD95FLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$ac_cv_f95_compiler_gnu
])


# _AC_LANG_ABBREV(Fortran 95)
# ---------------------------
m4_define([_AC_LANG_ABBREV(Fortran 95)], [f95])


## ---------------------- ##
## 2.Producing programs.  ##
## ---------------------- ##

# ------------------------ #
# 2b. Fortran 90 sources.  #
# ------------------------ #

# AC_LANG_SOURCE(Fortran 90)(BODY)
# --------------------------------
m4_copy([AC_LANG_SOURCE(Fortran 77)], [AC_LANG_SOURCE(Fortran 90)])


# AC_LANG_PROGRAM(Fortran 90)([PROLOGUE], [BODY])
# -----------------------------------------------
m4_define([AC_LANG_PROGRAM(Fortran 90)], [
program main
$1
$2
end program main
])

# AC_LANG_CALL(Fortran 90)(PROLOGUE, FUNCTION)
# --------------------------------------------
m4_define([AC_LANG_CALL(Fortran 90)],
[AC_LANG_PROGRAM([$1],
[call $2])])


# ------------------------ #
# 2c. Fortran 95 sources.  #
# ------------------------ #

# AC_LANG_SOURCE(Fortran 95)(BODY)
# --------------------------------
m4_copy([AC_LANG_SOURCE(Fortran 90)], [AC_LANG_SOURCE(Fortran 95)])

# AC_LANG_PROGRAM(Fortran 95)([PROLOGUE], [BODY])
# -----------------------------------------------
m4_copy([AC_LANG_PROGRAM(Fortran 90)], [AC_LANG_PROGRAM(Fortran 95)])

# AC_LANG_CALL(Fortran 95)(PROLOGUE, FUNCTION)
# --------------------------------------------
m4_copy([AC_LANG_CALL(Fortran 90)], [AC_LANG_CALL(Fortran 90)])


## -------------------------------------------- ##
## 3. Looking for Compilers and Preprocessors.  ##
## -------------------------------------------- ##

# ----------------------------- #
# 3b. The Fortran 90 compiler.  #
# ----------------------------- #


# AC_LANG_PREPROC(Fortran 90)
# ---------------------------
# Find the Fortran 90 preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran 90)],
[m4_warn([syntax],
         [$0: No preprocessor defined for ]_AC_LANG)])


# AC_LANG_COMPILER(Fortran 90)
# ----------------------------
# Find the Fortran 90 compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran 90)],
[AC_REQUIRE([AC_PROG_F90])])


# AC_PROG_F90([COMPILERS...])
# ---------------------------
# COMPILERS is a space separated list of Fortran 90 compilers to search
# for.
#
# Compilers are ordered by
#  1. F90, F95
#  2. Good/tested native compilers, bad/untested native compilers
#
# pgf90 is the Portland Group F90 compilers.
# xlf90/xlf95 are IBM (AIX) F90/F95 compilers.
# lf95 is the Lahey-Fujitsu compiler.
# epcf90 is the "Edinburgh Portable Compiler" F90.
# fort is the Compaq Fortran 90 (now 95) compiler for Tru64 and Linux/Alpha.
AC_DEFUN([AC_PROG_F90],
[AC_LANG_PUSH(Fortran 90)dnl
AC_ARG_VAR([F90],      [Fortran 90 compiler command])dnl
AC_ARG_VAR([F90FLAGS], [Fortran 90 compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
AC_CHECK_TOOLS(F90,
      [m4_default([$1],
                  [f90 xlf90 pgf90 epcf90 f95 xlf95 lf95 fort g95])])

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
# If we don't use `.F90' as extension, the preprocessor is not run on the
# input file.
ac_save_ext=$ac_ext
ac_ext=F90
_AC_LANG_COMPILER_GNU
ac_ext=$ac_save_ext
G95=`test $ac_compiler_gnu = yes && echo yes`
AC_LANG_POP(Fortran 90)dnl
])# AC_PROG_F90


# ----------------------------- #
# 3c. The Fortran 95 compiler.  #
# ----------------------------- #


# AC_LANG_PREPROC(Fortran 95)
# ---------------------------
# Find the Fortran 95 preprocessor.  Must be AC_DEFUN'd to be AC_REQUIRE'able.
AC_DEFUN([AC_LANG_PREPROC(Fortran 95)],
[m4_warn([syntax],
         [$0: No preprocessor defined for ]_AC_LANG)])


# AC_LANG_COMPILER(Fortran 95)
# ----------------------------
# Find the Fortran 95 compiler.  Must be AC_DEFUN'd to be
# AC_REQUIRE'able.
AC_DEFUN([AC_LANG_COMPILER(Fortran 95)],
[AC_REQUIRE([AC_PROG_F95])])


# AC_PROG_F95([COMPILERS...])
# ---------------------------
# COMPILERS is a space separated list of Fortran 95 compilers to search
# for.
#
# Compilers are ordered by
#  1. Good/tested native compilers, bad/untested native compilers
#
# xlf95 is the IBM (AIX) F95 compiler.
# lf95 is the Lahey-Fujitsu compiler.
# fort is the Compaq Fortran 90 (now 95) compiler for Tru64 and Linux/Alpha.
AC_DEFUN([AC_PROG_F95],
[AC_LANG_PUSH(Fortran 95)dnl
AC_ARG_VAR([F95],      [Fortran 95 compiler command])dnl
AC_ARG_VAR([F95FLAGS], [Fortran 95 compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
AC_CHECK_TOOLS(F95,
      [m4_default([$1],
                  [f95 xlf95 lf95 fort g95])])

m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
# If we don't use `.F95' as extension, the preprocessor is not run on the
# input file.
ac_save_ext=$ac_ext
ac_ext=F95
_AC_LANG_COMPILER_GNU
ac_ext=$ac_save_ext
G95=`test $ac_compiler_gnu = yes && echo yes`
AC_LANG_POP(Fortran 95)dnl
])# AC_PROG_F95


## ------------------------------- ##
## 4. Compilers' characteristics.  ##
## ------------------------------- ##


# ---------------------------------------- #
# 4b. Fortan 90 compiler characteristics.  #
# ---------------------------------------- #


# _AC_PROG_F90_V_OUTPUT([FLAG = $ac_cv_prog_f90_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (which default value, $ac_cv_prog_f90_v, is computed by
# _AC_PROG_F90_V), and return the output in $ac_f90_v_output.  This
# output is processed in the way expected by AC_F90_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_F90_V_OUTPUT],
[AC_REQUIRE([AC_PROG_F90])dnl
AC_LANG_PUSH(Fortran 90)dnl

AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 90 compiler in order to get
# "verbose" output that we can then parse for the Fortran 90 linker
# flags.
ac_save_F90FLAGS=$F90FLAGS
F90FLAGS="$F90FLAGS m4_default([$1], [$ac_cv_prog_f90_v])"

(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_f90_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f90_v_output" >&AS_MESSAGE_LOG_FD
F90FLAGS=$ac_save_F90FLAGS

rm -f conftest.*
AC_LANG_POP(Fortran 90)dnl

# If we are using xlf then replace all the commas with spaces.
if echo $ac_f90_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f90_v_output=`echo $ac_f90_v_output | sed 's/,/ /g'`
fi

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_f90_v_output | grep cft90 >/dev/null 2>&1; then
  ac_f90_v_output=`echo $ac_f90_v_output | sed "s/\"//g"`
fi[]dnl
])# _AC_PROG_F90_V_OUTPUT


# _AC_PROG_F90_V
# --------------
#
# Determine the flag that causes the Fortran 90 compiler to print
# information of library and object files (normally -v)
# Needed for AC_F90_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_F90_V],
[AC_CACHE_CHECK([how to get verbose linking output from $F90],
                [ac_cv_prog_f90_v],
[AC_LANG_ASSERT(Fortran 90)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_f90_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_F90_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_f90_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_f90_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_f90_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $F90])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_F90_V


# AC_F90_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran 90
# intrinsic and run-time libraries that are required to successfully
# link a Fortran 90 program or shared library.  The output variable
# F90LIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran 90, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran 90 compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran 90 intrinsic and run-time libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran 90 libraries.  Hence, the macro
# "AC_F90_LIBRARY_LDFLAGS" was created to determine these Fortran 90
# libraries.
#
# This macro was copied from the Fortran 77 version by Matthew D. Langston.
AC_DEFUN([AC_F90_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran 90)dnl
_AC_PROG_F90_V
AC_CACHE_CHECK([for Fortran 90 libraries], ac_cv_flibs,
[if test "x$F90LIBS" != "x"; then
  ac_cv_f90libs="$F90LIBS" # Let the user override the test.
else

_AC_PROG_F90_V_OUTPUT

ac_cv_f90libs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_f90_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
              ac_cv_f90libs="$ac_cv_f90libs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_f90libs)])
          ;;
          # Ignore these flags.
        -lang* | [-lcrt[012].o] | -lc | -lgcc | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_f90libs="$ac_cv_f90libs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
             "" | -*);;
             *)
		ac_arg="$ac_arg$[2]"

		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_f90libs, ,
                            [ac_arg="$ac_arg $ac_j"
                             ac_cv_f90libs="$ac_cv_f90libs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f90libs, ,
                          ac_cv_f90libs="$ac_cv_f90libs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f90_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&

        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_f90libs)
      ;;
esac
fi # test "x$F90LIBS" = "x"
])
F90LIBS="$ac_cv_f90libs"
AC_SUBST(F90LIBS)
AC_LANG_POP(Fortran 90)dnl
])# AC_F90_LIBRARY_LDFLAGS


# _AC_F90_NAME_MANGLING
# ---------------------
# Test for the name mangling scheme used by the Fortran 90 compiler.
#
# Sets ac_cv_f90_mangling. The value contains three fields, separated
# by commas:
#
# lower case / upper case:
#    case translation of the Fortan 90 symbols
# underscore / no underscore:
#    whether the compiler appends "_" to symbol names
# extra underscore / no extra underscore:
#    whether the compiler appends an extra "_" to symbol names already
#    containing at least one underscore
#
AC_DEFUN([_AC_F90_NAME_MANGLING],
[AC_REQUIRE([AC_F90_LIBRARY_LDFLAGS])dnl
AC_CACHE_CHECK([for Fortran 90 name-mangling scheme],
               ac_cv_f90_mangling,
[AC_LANG_PUSH(Fortran 90)dnl
AC_COMPILE_IFELSE(
[subroutine foobar()
return
end
subroutine foo_bar()
return
end],
[mv conftest.$ac_objext cf90_test.$ac_objext

  AC_LANG_PUSH(C)dnl

  ac_save_LIBS=$LIBS
  LIBS="cf90_test.$ac_objext $F90LIBS $LIBS"

  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      AC_TRY_LINK_FUNC($ac_func,
         [ac_success=yes; break 2])
    done
  done

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
        foobar)
           ac_case=lower
           ac_foo_bar=foo_bar
           ;;
        FOOBAR)
           ac_case=upper
           ac_foo_bar=FOO_BAR
           ;;
     esac

     ac_success_extra=no
     for ac_extra in "" "_"; do
        ac_func="$ac_foo_bar$ac_underscore$ac_extra"

        AC_TRY_LINK_FUNC($ac_func,
        [ac_success_extra=yes; break])
     done

     if test "$ac_success_extra" = "yes"; then
	ac_cv_f90_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_f90_mangling="$ac_cv_f90_mangling, no underscore"
	else
           ac_cv_f90_mangling="$ac_cv_f90_mangling, underscore"

        fi
        if test -z "$ac_extra"; then
           ac_cv_f90_mangling="$ac_cv_f90_mangling, no extra underscore"
	else
           ac_cv_f90_mangling="$ac_cv_f90_mangling, extra underscore"
        fi
      else
	ac_cv_f90_mangling="unknown"
      fi
  else
     ac_cv_f90_mangling="unknown"

  fi

  LIBS=$ac_save_LIBS
  AC_LANG_POP(C)dnl
  rm -f cf90_test* conftest*])
AC_LANG_POP(Fortran 90)dnl
])
])# _AC_F90_NAME_MANGLING

# The replacement is empty.
AU_DEFUN([AC_F90_NAME_MANGLING], [])


# AC_F90_WRAPPERS
# ---------------
# Defines C macros F90_FUNC(name,NAME) and F90_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran 90 compiler.
AC_DEFUN([AC_F90_WRAPPERS],
[AC_REQUIRE([_AC_F90_NAME_MANGLING])dnl
AH_TEMPLATE([F90_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran 90.])dnl
AH_TEMPLATE([F90_FUNC_],
    [As F90_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_f90_mangling in
  "lower case, no underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name]) ;;
  "lower case, no underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [name ## __]) ;;
  "upper case, no underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME]) ;;
  "upper case, no underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, no extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, extra underscore")
          AC_DEFINE([F90_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F90_FUNC_(name,NAME)], [NAME ## __]) ;;
  *)
          AC_MSG_WARN([unknown Fortran 90 name-mangling scheme])
          ;;
esac
])# AC_F90_WRAPPERS


# AC_F90_FUNC(NAME, [SHELLVAR = NAME])
# ------------------------------------
# For a Fortran subroutine of given NAME, define a shell variable
# $SHELLVAR to the Fortran 90 mangled name.  If the SHELLVAR
# argument is not supplied, it defaults to NAME.
AC_DEFUN([AC_F90_FUNC],
[AC_REQUIRE([_AC_F90_NAME_MANGLING])dnl
case $ac_cv_f90_mangling in
  upper*) ac_val="m4_toupper([$1])" ;;
  lower*) ac_val="m4_tolower([$1])" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f90_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac
m4_if(m4_index([$1],[_]),-1,[],
[case $ac_cv_f90_mangling in *," extra underscore"*) ac_val="$ac_val"_ ;; esac
])
m4_default([$2],[$1])="$ac_val"

])# AC_F90_FUNC


# ---------------------------------------- #
# 4c. Fortan 95 compiler characteristics.  #
# ---------------------------------------- #


# _AC_PROG_F95_V_OUTPUT([FLAG = $ac_cv_prog_f95_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (which default value, $ac_cv_prog_f95_v, is computed by
# _AC_PROG_F95_V), and return the output in $ac_f95_v_output.  This
# output is processed in the way expected by AC_F95_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_F95_V_OUTPUT],
[AC_REQUIRE([AC_PROG_F95])dnl
AC_LANG_PUSH(Fortran 95)dnl

AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran 95 compiler in order to get
# "verbose" output that we can then parse for the Fortran 95 linker
# flags.
ac_save_F95FLAGS=$F95FLAGS
F95FLAGS="$F95FLAGS m4_default([$1], [$ac_cv_prog_f95_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_f95_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_f95_v_output" >&AS_MESSAGE_LOG_FD
F95FLAGS=$ac_save_F95FLAGS

rm -f conftest.*
AC_LANG_POP(Fortran 95)dnl

# If we are using xlf then replace all the commas with spaces.
if echo $ac_f95_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_f95_v_output=`echo $ac_f95_v_output | sed 's/,/ /g'`
fi

# If we are using Cray Fortran then delete quotes.
# Use "\"" instead of '"' for font-lock-mode.
# FIXME: a more general fix for quoted arguments with spaces?
if echo $ac_f95_v_output | grep cft95 >/dev/null 2>&1; then
  ac_f95_v_output=`echo $ac_f95_v_output | sed "s/\"//g"`
fi[]dnl
])# _AC_PROG_F95_V_OUTPUT


# _AC_PROG_F95_V
# --------------
#
# Determine the flag that causes the Fortran 95 compiler to print
# information of library and object files (normally -v)
# Needed for AC_F95_LIBRARY_FLAGS
# Some compilers don't accept -v (Lahey: -verbose, xlf: -V, Fujitsu: -###)
AC_DEFUN([_AC_PROG_F95_V],
[AC_CACHE_CHECK([how to get verbose linking output from $F95],
                [ac_cv_prog_f95_v],
[AC_LANG_ASSERT(Fortran 95)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_f95_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_F95_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_f95_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_f95_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_f95_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $F95])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_F95_V


# AC_F95_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran 95
# intrinsic and run-time libraries that are required to successfully
# link a Fortran 95 program or shared library.  The output variable
# F95LIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran 95, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran 95 compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran 95 intrinsic and run-time libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran 95 libraries.  Hence, the macro
# "AC_F95_LIBRARY_LDFLAGS" was created to determine these Fortran 95
# libraries.
#
# This macro was copied from the Fortran 77 version by Matthew D. Langston.
AC_DEFUN([AC_F95_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(Fortran 95)dnl
_AC_PROG_F95_V
AC_CACHE_CHECK([for Fortran 95 libraries], ac_cv_flibs,
[if test "x$F95LIBS" != "x"; then
  ac_cv_f95libs="$F95LIBS" # Let the user override the test.
else

_AC_PROG_F95_V_OUTPUT

ac_cv_f95libs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_f95_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f95libs, ,
              ac_cv_f95libs="$ac_cv_f95libs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f95libs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_f95libs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -LANG:=*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_f95libs="$ac_cv_f95libs $ac_arg"
          ;;
        -[[LRuY]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
             "" | -*);;
             *)
		ac_arg="$ac_arg$[2]"

		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_f95libs, ,
                            [ac_arg="$ac_arg $ac_j"
                             ac_cv_f95libs="$ac_cv_f95libs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_f95libs, ,
                          ac_cv_f95libs="$ac_cv_f95libs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_f95_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&

        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_f95libs)
      ;;
esac
fi # test "x$F95LIBS" = "x"
])
F95LIBS="$ac_cv_f95libs"
AC_SUBST(F95LIBS)
AC_LANG_POP(Fortran 95)dnl
])# AC_F95_LIBRARY_LDFLAGS


# _AC_F95_NAME_MANGLING
# ---------------------
# Test for the name mangling scheme used by the Fortran 95 compiler.
#
# Sets ac_cv_f95_mangling. The value contains three fields, separated
# by commas:
#
# lower case / upper case:
#    case translation of the Fortan 95 symbols
# underscore / no underscore:
#    whether the compiler appends "_" to symbol names
# extra underscore / no extra underscore:
#    whether the compiler appends an extra "_" to symbol names already
#    containing at least one underscore
#
AC_DEFUN([_AC_F95_NAME_MANGLING],
[AC_REQUIRE([AC_F95_LIBRARY_LDFLAGS])dnl
AC_CACHE_CHECK([for Fortran 95 name-mangling scheme],
               ac_cv_f95_mangling,
[AC_LANG_PUSH(Fortran 95)dnl
AC_COMPILE_IFELSE(
[subroutine foobar()
return
end
subroutine foo_bar()
return
end],
[mv conftest.$ac_objext cf95_test.$ac_objext

  AC_LANG_PUSH(C)dnl

  ac_save_LIBS=$LIBS
  LIBS="cf95_test.$ac_objext $F95LIBS $LIBS"

  ac_success=no
  for ac_foobar in foobar FOOBAR; do
    for ac_underscore in "" "_"; do
      ac_func="$ac_foobar$ac_underscore"
      AC_TRY_LINK_FUNC($ac_func,
         [ac_success=yes; break 2])
    done
  done

  if test "$ac_success" = "yes"; then
     case $ac_foobar in
        foobar)
           ac_case=lower
           ac_foo_bar=foo_bar
           ;;
        FOOBAR)
           ac_case=upper
           ac_foo_bar=FOO_BAR
           ;;
     esac

     ac_success_extra=no
     for ac_extra in "" "_"; do
        ac_func="$ac_foo_bar$ac_underscore$ac_extra"

        AC_TRY_LINK_FUNC($ac_func,
        [ac_success_extra=yes; break])
     done

     if test "$ac_success_extra" = "yes"; then
	ac_cv_f95_mangling="$ac_case case"
        if test -z "$ac_underscore"; then
           ac_cv_f95_mangling="$ac_cv_f95_mangling, no underscore"
	else
           ac_cv_f95_mangling="$ac_cv_f95_mangling, underscore"

        fi
        if test -z "$ac_extra"; then
           ac_cv_f95_mangling="$ac_cv_f95_mangling, no extra underscore"
	else
           ac_cv_f95_mangling="$ac_cv_f95_mangling, extra underscore"
        fi
      else
	ac_cv_f95_mangling="unknown"
      fi
  else
     ac_cv_f95_mangling="unknown"

  fi

  LIBS=$ac_save_LIBS
  AC_LANG_POP(C)dnl
  rm -f cf95_test* conftest*])
AC_LANG_POP(Fortran 95)dnl
])
])# _AC_F95_NAME_MANGLING

# The replacement is empty.
AU_DEFUN([AC_F95_NAME_MANGLING], [])


# AC_F95_WRAPPERS
# ---------------
# Defines C macros F95_FUNC(name,NAME) and F95_FUNC_(name,NAME) to
# properly mangle the names of C identifiers, and C identifiers with
# underscores, respectively, so that they match the name mangling
# scheme used by the Fortran 95 compiler.
AC_DEFUN([AC_F95_WRAPPERS],
[AC_REQUIRE([_AC_F95_NAME_MANGLING])dnl
AH_TEMPLATE([F95_FUNC],
    [Define to a macro mangling the given C identifier (in lower and upper
     case), which must not contain underscores, for linking with Fortran 95.])dnl
AH_TEMPLATE([F95_FUNC_],
    [As F95_FUNC, but for C identifiers containing underscores.])dnl
case $ac_cv_f95_mangling in
  "lower case, no underscore, no extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [name])
          AC_DEFINE([F95_FUNC_(name,NAME)], [name]) ;;
  "lower case, no underscore, extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [name])
          AC_DEFINE([F95_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, no extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F95_FUNC_(name,NAME)], [name ## _]) ;;
  "lower case, underscore, extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [name ## _])
          AC_DEFINE([F95_FUNC_(name,NAME)], [name ## __]) ;;
  "upper case, no underscore, no extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F95_FUNC_(name,NAME)], [NAME]) ;;
  "upper case, no underscore, extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [NAME])
          AC_DEFINE([F95_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, no extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F95_FUNC_(name,NAME)], [NAME ## _]) ;;
  "upper case, underscore, extra underscore")
          AC_DEFINE([F95_FUNC(name,NAME)],  [NAME ## _])
          AC_DEFINE([F95_FUNC_(name,NAME)], [NAME ## __]) ;;
  *)
          AC_MSG_WARN([unknown Fortran 95 name-mangling scheme])
          ;;
esac
])# AC_F95_WRAPPERS


# AC_F95_FUNC(NAME, [SHELLVAR = NAME])
# ------------------------------------
# For a Fortran subroutine of given NAME, define a shell variable
# $SHELLVAR to the Fortran 95 mangled name.  If the SHELLVAR
# argument is not supplied, it defaults to NAME.
AC_DEFUN([AC_F95_FUNC],
[AC_REQUIRE([_AC_F95_NAME_MANGLING])dnl
case $ac_cv_f95_mangling in
  upper*) ac_val="m4_toupper([$1])" ;;
  lower*) ac_val="m4_tolower([$1])" ;;
  *)      ac_val="unknown" ;;
esac
case $ac_cv_f95_mangling in *," underscore"*) ac_val="$ac_val"_ ;; esac
m4_if(m4_index([$1],[_]),-1,[],
[case $ac_cv_f95_mangling in *," extra underscore"*) ac_val="$ac_val"_ ;; esac
])
m4_default([$2],[$1])="$ac_val"

])# AC_F95_FUNC


