dnl 
dnl AC_CXX_STD
dnl
dnl Description
dnl 
dnl Tests for various things that we might expect to find in
dnl C++ namespace std.
dnl
dnl Looks for <iostream> or else <iostream.h>.
dnl Looks for <iomanip> or else <iomanip.h>.
dnl Looks for <cmath>.
dnl Defines macro HAVE_FSTREAM_ATTACH if we have:
dnl  void ofstream::attach (int FILE)
dnl Defines macro HAVE_FSTREAM_OPEN if we have:
dnl  void ofstream::open (const char* FNAME, int MODE)
dnl Defines macro FSTREAM_OPEN_PROT if we have:
dnl  void ofstream::open (const char* FNAME, int MODE, int PROT)
dnl Defines macro HAVE_STD_IOSTREAM if this works:
dnl  #include <iostream>
dnl  std::cout<<"Hello World"<<std::endl;
dnl Defines macro HAVE_STD_STL if this works:
dnl  #include <list>
dnl  std::list<int> foo;
dnl 
dnl Copyright (C) 2003, Alex Tingle <alex.autoconf@firetree.net>
dnl 
dnl License:
dnl GNU General Public License
dnl [http://www.gnu.org/software/ac-archive/htmldoc/COPYING.html]
dnl with this special exception
dnl [http://www.gnu.org/software/ac-archive/htmldoc/COPYING-Exception.html]. 
dnl 

AC_DEFUN([AC_CXX_STD],[
  AC_CXX_HAVE_IOSTREAM
  AC_CXX_HAVE_IOMANIP
  AC_CHECK_HEADERS([cmath])
  AC_CXX_HAVE_STD_IOSTREAM
  AC_CXX_HAVE_STD_STL
  AC_CXX_HAVE_FSTREAM_ATTACH
  AC_CXX_HAVE_FSTREAM_OPEN
])


dnl 
dnl AC_CXX_HAVE_IOSTREAM
dnl Looks for <iostream> or else <iostream.h>.
dnl 

AC_DEFUN([AC_CXX_HAVE_IOSTREAM],[
  AC_CHECK_HEADER([iostream],[
    AC_DEFINE([HAVE_IOSTREAM],[1],[Define to 1 if you have the <iostream> header file.])
  ],[
    AC_CHECK_HEADERS([iostream.h])
  ])
])


dnl 
dnl AC_CXX_HAVE_FSTREAM_ATTACH
dnl Defines macro HAVE_FSTREAM_ATTACH if we have:
dnl  void ofstream::attach (int FILE)
dnl 

AC_DEFUN([AC_CXX_HAVE_FSTREAM_ATTACH],[
  AC_REQUIRE([AC_CXX_HAVE_STD_IOSTREAM])
  AC_CACHE_CHECK([for fstream::attach()],[ac_cv_cxx_have_fstream_attach],[
    ac_cv_cxx_have_fstream_attach=no
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_TRY_COMPILE([
#ifdef HAVE_IOSTREAM
#include <fstream>
#else
#include <fstream.h>
#endif
#ifdef HAVE_STD_IOSTREAM
using namespace std;
#endif
],[int fd=0;ofstream s;s.attach(fd);],
      [ac_cv_cxx_have_fstream_attach=yes])
    AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_have_fstream_attach" = yes; then
    AC_DEFINE([HAVE_FSTREAM_ATTACH],[1],[define if we have fstream::attach().])
  fi
])


dnl 
dnl AC_CXX_HAVE_FSTREAM_OPEN
dnl Defines macro HAVE_FSTREAM_OPEN if we have:
dnl  void ofstream::open (const char* FNAME, int MODE)
dnl Defines macro FSTREAM_OPEN_PROT if we have:
dnl  void ofstream::open (const char* FNAME, int MODE, int PROT)
dnl 

AC_DEFUN([AC_CXX_HAVE_FSTREAM_OPEN],[
  AC_REQUIRE([AC_CXX_HAVE_STD_IOSTREAM])
  AC_CACHE_CHECK([for fstream::open()],[ac_cv_cxx_have_fstream_open],[
    ac_cv_cxx_have_fstream_open=no
    ac_cv_cxx_fstream_open_prot=no
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    # Try with 2 parameters
    AC_TRY_COMPILE([
#ifdef HAVE_IOSTREAM
#include <fstream>
#else
#include <fstream.h>
#endif
#ifdef HAVE_STD_IOSTREAM
using namespace std;
#endif
],[ofstream s;s.open("conftest.txt",ios::out|ios::trunc);],
      [ac_cv_cxx_have_fstream_open=yes])
    # Try with mode parameter
    AC_TRY_COMPILE([
#ifdef HAVE_IOSTREAM
#include <fstream>
#else
#include <fstream.h>
#endif
#ifdef HAVE_STD_IOSTREAM
using namespace std;
#endif
],[ofstream s;s.open("conftest.txt",ios::out|ios::trunc,0666);],
      [ac_cv_cxx_fstream_open_prot=yes])
    AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_have_fstream_open" = yes; then
    AC_DEFINE([HAVE_FSTREAM_OPEN],[1],
      [define if we have fstream::open().])
  fi
  if test "$ac_cv_cxx_fstream_open_prot" = yes; then
    AC_DEFINE([FSTREAM_OPEN_PROT],[1],
      [define if fstream::open() accepts third parameter.])
  fi
])


dnl 
dnl AC_CXX_HAVE_IOMANIP
dnl Looks for <iomanip> or else <iomanip.h>.
dnl 

AC_DEFUN([AC_CXX_HAVE_IOMANIP],[
  AC_CHECK_HEADER([iomanip],[
    AC_DEFINE([HAVE_IOMANIP],[1],[Define to 1 if you have the <iomanip> header file.])
  ],[
    AC_CHECK_HEADERS([iomanip.h])
  ])
])


dnl 
dnl AC_CXX_HAVE_STD_IOSTREAM
dnl Defines macro HAVE_STD_IOSTREAM if this works:
dnl  #include <iostream>
dnl  std::cout<<"Hello World"<<std::endl;
dnl 

AC_DEFUN([AC_CXX_HAVE_STD_IOSTREAM],[
  AC_REQUIRE([AC_CXX_NAMESPACES])
  AC_REQUIRE([AC_CXX_HAVE_IOSTREAM])
  AC_CACHE_CHECK([for C++ iostream in namespace std],
    ac_cv_cxx_have_std_iostream,[
      ac_cv_cxx_have_std_iostream=no
      ac_cv_cxx_need_use_std_iostream=no
      if test "x$ac_cv_cxx_namespaces" = xyes; then
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_TRY_COMPILE([
#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif
],[std::cout<<"Hello World"<<std::endl;return 0;],
          [ac_cv_cxx_have_std_iostream=yes])
        AC_TRY_COMPILE([
#define __USE_STD_IOSTREAM 1
#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif
],[std::cout<<"Hello World"<<std::endl;return 0;],
          [ac_cv_cxx_have_std_iostream=yes;ac_cv_cxx_need_use_std_iostream=yes])
        AC_LANG_RESTORE
      fi
  ])
  if test "$ac_cv_cxx_have_std_iostream" = yes; then
    AC_DEFINE([HAVE_STD_IOSTREAM],[1],[define if C++ iostream is in namespace std.])
  fi
  if test "$ac_cv_cxx_need_use_std_iostream" = yes; then
    AC_DEFINE([__USE_STD_IOSTREAM],[1],[needed by DEC/Compaq/HP cxx to activate ANSI standard iostream.])
  fi
])


dnl 
dnl AC_CXX_HAVE_STD_STL
dnl Defines macro HAVE_STD_STL if this works:
dnl  #include <list>
dnl  std::list<int> foo;
dnl 

AC_DEFUN([AC_CXX_HAVE_STD_STL],[
  AC_REQUIRE([AC_CXX_NAMESPACES])
  AC_REQUIRE([AC_CXX_HAVE_STL])
  AC_CACHE_CHECK([for C++ Standard Template Library in namespace std.],
    ac_cv_cxx_have_std_stl,[
      ac_cv_cxx_have_std_stl=no
      if test "x$ac_cv_cxx_namespaces" = xyes; then
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_TRY_COMPILE([#include <list>
          ],[std::list<int> foo;return 0;],
          [ac_cv_cxx_have_std_stl=yes])
        AC_LANG_RESTORE
      fi
  ])
  if test "$ac_cv_cxx_have_std_stl" = yes; then
    AC_DEFINE([HAVE_STD_STL],[1],[define if C++ Standard Template Library is in namespace std])
  fi
])

