.TH WAISSERVER 1 "Tue Apr 28 1992" "Thinking Machines"
.SH NAME
waisserver \- serves WAIS requests

.SH SYNOPSIS
.B waisserver
[ \-p [ port_number ] ] [ \-s ] [ \-d directory ] [\ \-e [ pathname ] ] 
[\ \-l log_level ] [\-u user ] [ \-v ]
.TP
.B waisserver.d
[ same arguments ]

.SH DESCRIPTION
.PP
Part of the Wide Area Information Server system.
.PP
waisserver will take WAIS requests from a TCP port or standard-io and
return the appropriate response.  If the name of the command is
waisserver.d, then it is assumed it is running from inetd, and it uses
stdio for its I/O sockets.  See the examples below for inetd.conf.

In addition, waisserver can act like an anonymous FTP server.  When the
server is given a document of type "FTP" as a relevant document, it will
build a result list from the directory of the file.  Subdirectories
may be listed using by adding them to the relevant document list.

Note that a minimal level of security is present in two forms:

1.  The server will never present directories above the default server
directory (-d option, described below).

2.  The server will only build a directory listing from a file of type FTP,
and that file must be in the specified database.

.SH OPTIONS
.TP 10
\fB\-p\fI [ port ] 
Listen to the port.  If the port is supplied, then that port number is
used.  If it is not supplied then the Z39.50 port (210) is used.
.TP 10
\fB\-s
listen to standard I/O for queries.
.TP 10
\fB\-d\fI directory
Use this directory as the default location of the indexes.  Therefore if
the directory were /usr/local, then the database foo would be found in
/usr/local/foo (see waisindex for how to create an index)
.TP 10
\fB\-e\fI [ filename ]
Redirect error output to pathname, if supplied, or to /dev/null.
Error output defaults to stderr, unless -s is selected,
in which case it defaults to /dev/null.
.TP 10
\fB\-l\fI log_level 
set logging level.  Currently only levels 0, 1, 5 and 10 are meaningful:
Level 0 means log nothing (silent).  Level 1 logs only errors and warnings
(messages of HIGH priority), level 5 logs messages of MEDIUM priority (like
client init info).  Level 10 logs everything.
.TP 10
\fB\-u\fI user
Set the server's user id to the user specified after attaching the tcp-port. 
This is only used if the server is started as root.
.TP 10
\fB\-v
Print the current version and date of the server.
.SH EXAMPLES
The following are examples of \fIwaisserver\fR usage:
.TP
waisserver\ -p\ 8000\ -d\ wais-sources\ -e\ server.log
Runs waisserver as a standalone server, using tcp port 8000 on directory
wais-sources writing messages to server.log
.PP
Some example inetd.conf entries (note, these must be on one line in
inted.conf):
.PP
hpux 7.0/800, Interactive/386 2.2.1:
.TP 1
z3950\ stream\ tcp\ nowait\ root\ /etc/waisserver\ waisserver.d
-d\ /wais-sources\ -e\ /server.log
.PP
Ultrix 4.1:
.TP 1
z3950\ stream\ tcp\ nowait\ /etc/waisserver\ waisserver.d
-d\ /wais-sources\ -e\ /server.log
.PP
Also, add the next line to /etc/services, and tickle your YP server:
.PP
z3950\ 210/tcp\ #\ wide\ area\ information\ server\ (wais)
.SH SEE ALSO
.PP
waissearch(1), waisindex(1), waissearch-gmacs(1), xwais(1), xwaisq(1),
inetd(8C), inetd.conf(5)
.PP
Wide Area Information Servers Concepts by Brewster Kahle.
.br
Brewster@think.com
.SH DIAGNOSTICS
.PP
The diagnostics produced by the waisserver are meant to be self-explanatory.
.SH BUGS
Malformed protocol packets can cause the server to dump core (segmentation
violation).  These are logged in the server's log file.

