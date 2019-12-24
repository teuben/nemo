#
# NEMO still needs csh for many scripts
#
# @todo   add a --with-csh flag so it overrides this error, just in case....
#

AC_DEFUN([AC_PROG_CSH],
  [AC_MSG_CHECKING([for the (t)csh shell])
  if test -f /bin/csh; then
     have_csh="yes"
  else
     have_csh="no"
  fi
  AC_MSG_RESULT([$have_csh])
  if test "$have_csh" = "no"; then
     AC_MSG_ERROR([/bin/csh is required])
  fi
  ]
)

