$ CC /NOLIST/OBJECT=TR2LATEX.OBJ /NoDebug /Optimize /Define = ("HAVE_SGTTY=0","NO_SGTTY") TR2LATEX.C
$ CC /NOLIST/OBJECT=TR.OBJ /NoDebug /Optimize /Define = ("HAVE_SGTTY=0","NO_SGTTY") TR.C
$ CC /NOLIST/OBJECT=SUBS.OBJ /NoDebug /Optimize /Define = ("HAVE_SGTTY=0","NO_SGTTY") SUBS.C
$ CC /NOLIST/OBJECT=VERSION.OBJ /NoDebug /Optimize /Define = ("HAVE_SGTTY=0","NO_SGTTY") VERSION.C
$ LINK /TRACE/NOMAP/EXEC=TR2LATEX.EXE /NoDebug tr2latex,tr,subs,version,vaxcrtl.opt /Option
$ continue
