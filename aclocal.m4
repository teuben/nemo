# generated automatically by aclocal 1.11.3 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004,
# 2005, 2006, 2007, 2008, 2009, 2010, 2011 Free Software Foundation,
# Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_f77_cmain_fflags.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_F77_CMAIN_FFLAGS([ACTION-IF-SUCCEED], [ACTION-IF-FAIL])
#
# DESCRIPTION
#
#   This macro figures out if extra Fortran compiler flags are required in
#   order to use the Fortran linker to link programs where the main()
#   function is defined via C (or other language). On some systems, notably
#   the Alpha with Compaq compilers, the Fortran libraries have their own
#   main() function which must be disabled.
#
#   Runs ACTION-IF-SUCCEED if successful, and ACTION-IF-FAIL if not. Defines
#   the output variable F77_CMAIN_FFLAGS to any discovered flags. (If
#   ACTION-IF-FAIL is not specified, defaults to halting with an error.)
#
#   This macro is especially useful in conjunction with automake, since by
#   default automake uses $F77 to link programs mixing C and Fortran,
#   leading to a link error on some systems. In this case, you should set
#   the FFLAGS for that program to include F77_CMAIN_FFLAGS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 5

AU_ALIAS([ACX_F77_CMAIN_FFLAGS], [AX_F77_CMAIN_FFLAGS])
AC_DEFUN([AX_F77_CMAIN_FFLAGS],
[AC_CACHE_CHECK([for f77 flags to use C main function], ax_cv_f77_cmain_fflags,
[ax_cv_f77_cmain_fflags="unknown"
AC_LANG_PUSH(C)
AC_COMPILE_IFELSE([[int main(void) { return 0; }]],
		  [mv conftest.$ac_objext conftest_cmain.$ac_objext],
		  [ax_cv_f77_cmain_fflags=error])
AC_LANG_POP(C)
if test "x$ax_cv_f77_cmain_fflags" != xerror; then
    AC_LANG_PUSH(Fortran 77)
    ax_save_LIBS=$LIBS
    LIBS="conftest_cmain.$ac_objext $LIBS"
    ax_save_FFLAGS=$FFLAGS
    for ax_flag in none -nofor_main; do
	case $ax_flag in
	    none) FFLAGS=$ax_save_FFLAGS ;;
	    *)    FFLAGS="$ax_save_FFLAGS $ax_flag" ;;
	esac
	AC_LINK_IFELSE([
      subroutine foobar()
      return
      end
], [ax_cv_f77_cmain_fflags=$ax_flag; break]);
    done
    FFLAGS=$ax_save_FFLAGS
    LIBS=$ax_save_LIBS
    AC_LANG_POP(Fortran 77)
fi])
    case $ax_cv_f77_cmain_fflags in
	error|unknown)
	    F77_CMAIN_FFLAGS=""
	    ifelse([$2],,[AC_MSG_ERROR([cannot link C main with Fortran])],[$2])
	    ;;
	*)
	    if test "x$ax_cv_f77_cmain_fflags" = xnone; then
		F77_CMAIN_FFLAGS=""
	    else
		F77_CMAIN_FFLAGS="$ax_cv_f77_cmain_fflags"
	    fi
	    $1
	    ;;
    esac
    AC_SUBST(F77_CMAIN_FFLAGS)
])

# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000


AC_DEFUN([AX_PATH_GSL],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no)
  min_gsl_version=ifelse([$1], ,0.2.5,$1)
  AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
  no_gsl=""
  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
  else
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
       gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
       gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
       gsl_micro_version=0
    fi

    if test "x$enable_gsltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSL_CFLAGS"
      LIBS="$LIBS $GSL_LIBS"

      rm -f conf.gsltest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gsltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gsl_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gsl_version");
     exit(1);
   }

   if (($gsl_major_version > major) ||
      (($gsl_major_version == major) && ($gsl_minor_version > minor)) ||
      (($gsl_major_version == major) && ($gsl_minor_version == minor) && ($gsl_micro_version >= micro)))
     { 
       exit(0);
     }   
   else
     {
       exit(1);
     }
}

],, no_gsl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gsl" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GSL_CONFIG" = "no" ; then
       echo "*** The gsl-config script installed by GSL could not be found"
       echo "*** If GSL was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GSL_CONFIG environment variable to the"
       echo "*** full path to gsl-config."
     else
       if test -f conf.gsltest ; then
        :
       else
          echo "*** Could not run GSL test program, checking why..."
          CFLAGS="$CFLAGS $GSL_CFLAGS"
          LIBS="$LIBS $GSL_LIBS"
          AC_TRY_LINK([
#include <stdio.h>
],      [ return 0; ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GSL or finding the wrong"
          echo "*** version of GSL. If it is not finding GSL, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH"],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GSL was incorrectly installed"
          echo "*** or that you have moved GSL since it was installed. In the latter case, you"
          echo "*** may want to edit the gsl-config script: $GSL_CONFIG" ])
          CFLAGS="$ac_save_CFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
#     GSL_CFLAGS=""
#     GSL_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)
  rm -f conf.gsltest
])

AU_ALIAS([AM_PATH_GSL], [AX_PATH_GSL])

m4_include([acinclude.m4])
