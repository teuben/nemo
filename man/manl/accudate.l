.TH CPROTO 1 "..."
.SH NAME
\fBaccudate\fP \- 
.SH DESCRIPTION
An alternative to \fIdate(1)\fP, with more accurate output.
.SH OPTIONS
.nf
 
   [no args]     write time as YYYY-MM-DD:hh:mm:ss.lll  to stdout
   -t0           write time as "seconds microseconds" to stdout
   -df seconds microseconds
                 write elapsed time as "DDDD-hh:mm.ss.lll" to stdout
   -ds seconds microseconds
                 write elapsed time as "ssssss.lll" to stdout
   -h            print this help message (also -help --h --help -? --?)
   -i            emit version, copyright, license and contact information

.fi
.SH AUTHOR
David Mathog (Caltech)
.SH FILES
src/tools/misc/accudate.c
.SH HISTORY
.nf
.ta +1i +2i
4-jan-2003    imported into NEMO  PJT
8-jan-2005    man page written    PJT
.fi
