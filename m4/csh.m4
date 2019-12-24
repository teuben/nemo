#
# NEMO still needs csh for many scripts
#      but with the flag --without-csh you can bypass this test
#      but ... YMMV
#

AC_DEFUN([AC_PROG_CSH],
[
  AC_ARG_WITH(csh, [  --with-csh              bypass test for csh requirement], with_csh=$withval, with_csh=yes)
  AC_MSG_CHECKING([for the (t)csh shell $with_csh])
  if test "$with_csh" = "yes"; then
     if test -f /bin/csh; then
        have_csh="yes"
     else
        have_csh="no"
     fi
     AC_MSG_RESULT([$have_csh])
     if test "$have_csh" = "no"; then
        AC_MSG_ERROR([/bin/csh is required - you can bypass with --without-csh (YMMV)])
     fi
  fi
  ]
)

