.TH WAISSEARCH 1 "Thu Feb 20 1992" "Thinking Machines"
.SH NAME
waissearch \- Simple shell user interface for making WAIS requests

.SH SYNOPSIS
.PP
waissearch [\-h host-machine] [\-p [service-or-port]] [\-d database]
[\-m max-results] word word...

.SH DESCRIPTION
Part of the Wide Area Information Server system.
.PP
waissearch makes a WAIS request over TCP/IP.

.SH OPTIONS
.TP 10
\fB\-h\fI host-machine
Server machine name as a IP name.  Defaults to local machine.
.TP 10
\fB\-p\fI [ service-or-port ]
Use this port to communicate.  If it is not supplied, then z39_50 is used
from the services file.
.TP 10
\fB\-m\fI max-results
Set the maximum number of results returned for a query.  Default is 40.
.TP 10
\fB\-d\fI database
The database to use in the query.  By convention, many servers support the
INFO database.
.TP 10
\fB\-v
Print the version number and date of this command.

.SH EXAMPLE
.TP
to query the directory of servers about tcp servers:
.br
waissearch \-h quake.think.com \-d directory\-of\-servers \\ \-p 210 tcp
.TP
to ask what tcp databases are supported on a particular server:
.br
waissearch \-h quake.think.com \-d INFO \-p 210 tcp
.TP
to ask for information on a database:
.br
waissearch \-h quake.think.com \-d directory\-of\-servers \\ \-p 210 help

.SH SEE ALSO
waisserver(1), waisindex(1), waisq(1), xwais(1), xwaisq(1)
.PP
Wide Area Information Servers Concepts by Brewster Kahle.
.br
Brewster@think.com

.SH DIAGNOSTICS
The diagnostics produced by the waissearch are meant to be self\-explanatory.

.SH BUGS

