.TH SHELL 1l
.\" $Header$
.SH NAME
shell \- Interpreter for shell archives
.SH SYNOPSIS
.B shell
[ file...  ]
.SH DESCRIPTION
This program interprets enough UNIX shell syntax, and command usage,
to enable it to unpack many different types of UNIX shell archives,
or ``shar's.''
It is primarily intended to be used on non-UNIX systems that need to
unpack such archives.
.PP
.I Shell
does
.B not
check for security holes, and will blithely execute commands like
.RS
cp /dev/null /etc/passwd
.RE
which, if executed by the super-user, can be disastrous.
.PP
The
.I shell
parser is line-based, where lines are then split into tokens; it is not a
true token-based parser.
In general, it is best if all
.I sh
keywords that can appear alone on a line do so, and that compound
commands (i.e., using a semi-colon) be avoided.
For more details on the syntax, see the source (sorry...).
.SH BUGS
It is probably easier to write a true portable replacement for /bin/sh
than it is to write something which understands all shar formats.
.SH SEE ALSO
shar(1l).
