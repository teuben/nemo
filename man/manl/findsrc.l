.TH FINDSRC 1l
.\" $Header$
.SH NAME
findsrc \- walk directories, trying to find source files
.SH SYNOPSIS
.B findsrc
[
.B \-.
] [
.BI \-d\| y_or_n
] [
.BI \-o\| output_file
] [
.B \-R
] [
.B \-S
] [
.B \-v
] [ file... ]
.SH DESCRIPTION
.I Findsrc
recursively examines all directories and files specified on the command
line, and determines, based on the file name, whether the file contains
source code or not.
All directories are listed first, followed by all regular files,
with one item per line.
.PP
If
.I findsrc
is unable to make a decision, it invokes the
.IR file (1)
command, and prompts the user for a decision.
In reply to the prompt, type the letter ``y'' or ``n'' (either case);
RETURN means yes.
If the reply starts with an exclamation point, the rest of the line
is passed off to a sub-shell and the question is repeated.
The ``\-d'' option may be used with an argument of ``y'' or ``n''
to by-pass the interaction, and provide a default answer.
.PP
The ``\-o'' option may be used to specify an output filename.
This is designed to prevent confusion if a command like the following
is executed:
.RS
findsrc . * >LIST
.RE
.PP
By default,
.I findsrc
ignores files whose names begin with a period, like ``.newsrc'' or
``.tags''; such files may be included by using the ``\-.'' option.
.I Findsrc
also normally ignores
RCS and SCCS files and directories; using either the ``\-R'' or ``\-S''
option causes both to be included.
.PP
.I Findsrc
normally lists only the names of those files that have been selected.
If the ``\-v'' option is used, rejected files are also listed preceeded
by the word ``PUNTED.''
.PP
If no files are specified on the command line, the program operates as
a filter, reading a list of files and directories from the standard
input, one per line.
.SH "SEE ALSO"
makekit(1l).
