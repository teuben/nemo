.TH TKRUN 1 "24 December 2009"
.SH NAME
tkrun \- automatic GUI builder and frontend to run shell scripts or programs
.SH SYNOPSIS
\fBtkrun\fP [-option value...] \fIscript\fP
.SH DESCRIPTION
\fItkrun\fP presents an interactive graphical user interface to
a shell \fIscript\fP or program that conforms to a ``key=val'' interface.
It does this by extracting
named keywords from the shell script and building a
Tcl/Tk run script that provideds the interface to running the
script itself (or any program for that matter, see the \fB-r\fP option
below).
.SH OPTIONS
The following options parameters are recognized (the ones with a * still
have to be implemented):
.TP 15
\fB-r \fP\fIrun_path\fP
Override default name of executable to run in the tk-script
Should be in your $PATH.
(Default: \fIscript\fP).
.TP
\fB-t \fP\fItk_path\fP
Override default name of the tk script to run at the end of \fBtkrun\fP
(Default: \fIscript\fP\fB.tk\fP in the local directory).
.TP
\fB-w \fP\fIwish_path\fP
Override location and name of the \fIwish(1)\fP executable. 
Should be in your $PATH.
(Default: \fBwish\fP).
.TP
\fB-x\fP
Do not run the script, merely create the Tcl/Tk script.
.TP
\fB-e\fP
Show an example/template. Serves as a reminder how to create a new
script from scratch. 
.TP
\fB-c\fP
* Cleanup the script after done (maybe reverse the logic here, and cleanup
by default; currently no cleanup is done)
.TP
\fB-a\fP
* Autorun the script after each change. (toggle button on menu?)
.TP
\fB-n\fP
* Do not use keywords that have no value in running the program.
.TP
\fB-h\fP
* A bit of inline help
.SH SCRIPT
The following tags are recongized GUI builders, using
\fB#>\fP at the start of a line to mark a keyword
.nf

    #> ENTRY    \fIkey=val\fP
    #> SCALE    \fIkey=val\fP     \fIfrom:to:step\fP
    #> RADIO    \fIkey=val\fP     \fIopt1,opt2,opt3,...\fP
    #> CHECK    \fIkey=val\fP     \fIopt1,opt2,opt3,...\fP
    #> IFILE    \fIkey=val\fP     
    #> OFILE    \fIkey=val\fP
    #<      \fI-tkopt1 tkval1 -tkopt2 -tkval2 ... \fP
    #
    #> HELP  \fIkey\fP    One-line helpers for this particular keyword
    #
    #> -EXEC
    #           any code inbetween here is not executed by tkrun
    #> +EXEC
    #
.fi
It is also possible to put these lines on the same line as the shell variable
in e.g. c-shell scripts, e.g.:
.nf
    set a=1    #> SCALE  a=1     0:10:0.1
.fi
but this is a new feature and will possibly be modified to not have to specify
the variable setting "a=1"  twice.
.SH CAVEATS
Resetting (right click on the keyword) a checkbutton doesn't work properly.
.PP
checkbuttons and radiobuttons do not layout their buttons efficiently, and
the overal layout is limited to fit your (vertical) screensize. Planned
fix is to use the \fIgrid(3)\fP geometry manager, instead of
\fIpack(3)\fP.
.PP
The additional Tcl/Tk options (on a \fB#<\fP line immediately following the
\fB#>\fP line) are not working.
.PP
The code is horribly non-defensive, and does not handle syntax errors on the
\fB#>\fP tag lines very well.
.SH TODO
Help is not implemented yet. Neither is the -EXEC/+EXEC section to block
out sections of the script.
.SH SEE ALSO
wish(1), zenity(1)
.PP
.nf
.ta +2i
http://sunscript.sun.com/	Home Page of Tcl/Tk (now gone)
http://www.neosoft.com/		Archive of Tcl/Tk tools 
http://tcl.sourceforge.net/	Current home page and source code repo
.fi
.SH EXAMPLES
.nf
    % \fBtkrun design\fP
.fi
Would create a Tcl/Tk script called ``\fBdesign.tk\fP'', extract keyword info
and run this script, which in turn will call ``\fBdesign\fP'' with your
selection of ``\fIkeyword=value\fP'' pairs.
.nf
    % \fBtkrun -r tsf itemize.csh\fP
.fi
parses the GUI elements from \fBitemize.csh\fP, and will call the
\fBtsf\fP program.
.PP
Here is a very simple script, called \fBtestscript\fP in the source code
distribution:
.nf
#! /bin/csh -f
#
#       here are some basic GUI elements
#
#>  IFILE   in=
#>  OFILE   out=
#>  ENTRY   eps=0.01
#>  RADIO   mode=gauss              gauss,newton,leibniz
#>  CHECK   options=mean,sigma      sum,mean,sigma,skewness,kurtosis
#>  SCALE   n=3.141592              0:10:0.01
#<          -digits 2
#<  
#		parse named arguments
foreach a ($*)
   set $a
end

#>  -EXEC
echo If this is executed, tkrun does not work properly yet
#>  +EXEC

echo ARGS: in=$in out=$out eps=$eps mode=$mode options=$options n=$n

echo TESTSCRIPT
echo 0 : $0
set count=0

again:
  if ($#argv == 0) exit 0
  @ count++
  echo $count : \"$1\"
  shift argv
  goto again
.fi
.SH AUTHOR
Peter Teuben (teuben@astro.umd.edu) Astronomy Department,University of Maryland,College Park,MD,USA)
.SH FILES
.nf
.ta +2.0i
~/.tkrunrc	(optional) global configuration file
$run.key	keyword file in which key=val pairs are maintained
.fi
.SH "UPDATE HISTORY"
.nf
.ta +1.0i +4.5i
18-aug-97	V0.1 released for ADASS97 paper (tcl/tk 7.6/4.2) 	PJT
30-jan-98	V0.2 added -e, and some minor features	PJT
22-apr-04	V1.2 added option to share the set and #> comment on one line	PJT
12-jun-08	V1.4 updated doc and minor code cleanup	PJT
24-dec-09	V1.7 implemented the load/save of key=val files    	PJT
.fi
