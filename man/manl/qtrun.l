.TH QTRUN 1 "8 December 2023"

.SH "NAME"
qtrun \- automatic GUI builder and frontend to run shell scripts or programs

.SH "SYNOPSIS"
\fBqtrun\fP [-option value...] \fIscript\fP

.SH "DESCRIPTION"
\fIqtrun\fP presents an interactive graphical user interface to
a shell \fIscript\fP or program that conforms to a ``key=val'' interface.
.PP
It is based off \fItkrun(lNEMO)\fP, and supposed to be mostly compatible
with it. A few format enhancements were proposed for qtrun. See "TKRUN V2"
below.

.SH "OPTIONS"
The following options parameters are currently recognized. 
.TP 15
\fB-h, --help\fP
help
.TP
\fB-d, --debug\fP
debug
.TP
\fB-v, --version\fP
show version

.SH "TKRUN V2"
In version 2 of this format we have implemented this via this new python based qtrun
Here the format can be run via bash, csh and python in the
.nf

    key=val  # help    #> GUI

.fi
specification, but the client software (bash/csh/python/....) is responsible to properly parse the \fIkey=val\fP
command line arguments.  Here are the proposed GUI tags with generic key=val/help sections

.EX

   in=foo        # input file      #> IFILE
   out=bar       # output file     #> OFILE
   text=hello    # some text       #> ENTRY
   mode=gauss    # the mode        #> RADIO gauss,newton,leibniz
   stats=mean    # what to show    #> CHECK sum,mean,sigma,skewness,kurtosis
   n=3.141592    # the n value     #> SCALE 0:10:0.01
   k=3.141592    # k, use n        #> LINK  n

.EE

.SH "SEE ALSO"
wish(1), zenity(1), qtun(1NEMO)
.nf
https://github.com/teuben/agui - development version 
.fi

.SH "AUTHOR"
Kylie Gong, Peter Teuben

.SH "FILES"
.nf
.ta +2.5i
~/.tkrunrc	(optional) global configuration file
$run.key	keyword file in which key=val pairs are maintained
$NEMO/scripts/templates 	example templates of key=val parsing
.fi

.SH "HISTORY"
.nf
.ta +1.5i +5.5i
aug-2023	experimental version released	KG/PJT
.fi
