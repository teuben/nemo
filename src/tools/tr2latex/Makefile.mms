# MMS Description file for for tr2LaTeX
#  Converted from the Unix Makefile
#   by Richard L. Dyson 26-NOV-1991 (dyson@iowasp.physics.uiowa.edu)
#
#   Updated  4-MAR-1992 for tr2latex v2.0   RLD
#

DEBUG  = /NoDebug
OPTIMIZE = /Optimize
CFLAGS = $(CFLAGS) $(DEBUG) $(OPTIMIZE) /Define = ("HAVE_SGTTY=0","NO_SGTTY")

SRCS = tr2latex.c tr.c subs.c version.c
HFILES = setups.h simil.h greek.h macros.h maths.h flip.h forbid.h

all :	tr2latex.exe tr2latex.hlb
	@ Continue

help :	tr2latex.hlb
	@ Continue

tr2latex.exe :	tr2latex.obj tr.obj subs.obj version.obj vaxcrtl.opt
	$(LINK) $(LINKFLAGS) $(DEBUG) tr2latex,tr,subs,version,vaxcrtl.opt /Option

tr2latex.hlb :	tr2latex.hlp
tr2latex.obj :	tr2latex.c setups.h
tr.obj :	tr.c setups.h
subs.obj :	subs.c $(HFILES)

clean :
	@- Set Protection = Owner:RWED *.*;*
	@- Purge /NoLog *.*
	@- Set Protection = Owner:RWE *.*;*
