# Determine whether ps is posix or bsd.
# AC_PATH_PROG(PS, ps)
#AC_CACHE_CHECK([for POSIX or BSD ps syntax], ac_cv_prog_ps_syntax, [
        if $PS ax -o rss > /dev/null 2>&1; then
                ac_cv_prog_ps_ax=yes
        else
                ac_cv_prog_ps_ax=no
        fi
        if $PS -ef -o rss > /dev/null 2>&1; then
                ac_cv_prog_ps_ef=yes
        else
                ac_cv_prog_ps_ef=no
        fi
        if test "$ac_cv_prog_ps_ax" = yes; then
                ac_cv_prog_ps_syntax=BSD
        else
                if test "$ac_cv_prog_ps_ef" = yes; then
                        ac_cv_prog_ps_syntax=POSIX
                else
                        AC_MSG_ERROR(Could not determine ps syntax)
                fi
        fi
])

if (UNIX)
   if(NOT ps_cmd)
   	  find_program(ps_cmd "ps")
   endif(NOT ps_cmd)
 

   if(ps_cmd)
	execute_process(COMMAND "${ps_cmd}")