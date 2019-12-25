dnl
dnl				acx_plplot.m4
dnl
dnl Figure out if the PLPlot library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
dnl
dnl	License:		GNU General Public License
dnl
dnl	AstrOmatic software is free software: you can redistribute it and/or
dnl	modify it under the terms of the GNU General Public License as
dnl	published by the Free Software Foundation, either version 3 of the
dnl	License, or (at your option) any later version.
dnl	AstrOmatic software is distributed in the hope that it will be useful,
dnl	but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl	GNU General Public License for more details.
dnl	You should have received a copy of the GNU General Public License
dnl	along with AstrOmatic software.
dnl	If not, see <http://www.gnu.org/licenses/>.
dnl
dnl	Last modified:		10/10/2010
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_PLPLOT([PLPLOT_LIBDIR, PLPLOT_INCDIR,
dnl                      [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the PlPlot library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$PLPLOT_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $PLPLOT_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if PlPlot
dnl is found (HAVE_PLPLOT is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_PLPLOT], [
AC_REQUIRE([AC_CANONICAL_HOST])

PLPLOT_LIBS=""
OLIBS="$LIBS"
LIBS=""

acx_plplot_ok=yes
acx_plplotpkg_ok=no
if test x$2 = x && test x$1 = x; then
  AC_MSG_CHECKING([for PLPlot pkg-config info])
  if pkg-config --exists plplotd; then
    AC_MSG_RESULT([yes])
    [PLPLOT_CFLAGS=`pkg-config --cflags plplotd`]
    [PLPLOT_LIBS=`pkg-config --libs plplotd`]
    AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
    AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
    acx_plplotpkg_ok=yes
  else
    AC_MSG_RESULT([no])
  fi
fi
if test x$acx_plplotpkg_ok = xno; then
  if test x$2 = x; then
    AC_CHECK_HEADER(plplot.h, [acx_plplothead_ok=yes], [acx_plplothead_ok=no])
    if test x$acx_plplothead_ok = xyes; then
      AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
      AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
    else
      AC_CHECK_HEADER(plplot/plplot.h,
		[acx_plplothead_ok=yes], [acx_plplothead_ok=no])
      if test x$acx_plplothead_ok = xyes; then
        AC_DEFINE(PLPLOT_H, "plplot/plplot.h", [PLPlot header filename.])
        AC_DEFINE(PLPLOTP_H, "plplot/plplotP.h",
		[PLPlot private header filename.])
      else
        acx_plplot_ok=no
      fi
    fi
  else
    AC_CHECK_HEADER($2/plplot.h,
		[acx_plplothead_ok=yes], [acx_plplothead_ok=no])
    if test x$acx_plplothead_ok = xyes; then
      AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
      AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
     [PLPLOT_CFLAGS="-I$2"]
    else
      acx_plplot_ok=no
    fi
  fi
  if test x$1 = x; then
    AC_CHECK_LIB(plplotd, c_plinit,, [acx_plplot_ok=no])
    [PLPLOT_LIBS="-lplplotd"]
  else
    AC_CHECK_LIB(plplotd, c_plinit,, [acx_plplot_ok=no], [-L$1])
    [PLPLOT_LIBS="-L$1 -lplplotd"]
  fi
fi

LIBS="$OLIBS"
if test x$acx_plplot_ok = xyes; then
  AC_SUBST(PLPLOT_CFLAGS)
  AC_SUBST(PLPLOT_LIBS)
  AC_DEFINE(HAVE_PLPLOT,1,
	[Define if you have the PLPlot libraries and header files.])
  $3
else
  $4
fi

])dnl ACX_PLPLOT

