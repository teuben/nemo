.TH UNSHAR 1l
.\" $Header$
.SH NAME
unshar \- unpack shell archives from news, mail, notes, etc.
.SH SYNOPSIS
.B unshar
[
.BI \-c\| directory
] [
.BI \-d\| directory
] [
.BI \-h\| file
] [
.B \-f
] [
.B \-n
] [
.B \-s
] [ file... ]
.SH DESCRIPTION
.I Unshar
removes mail and news header lines from its input, and feeds the remainder
to
.IR /bin/sh (1)
so that a shell archive can be properly unpacked.
If no files are specified,
.I unshar
reads from standard input.
The program is designed to be useful when unpacking archives directly
from the news or mail systems (e.g., s | unshar).
.PP
.I Unshar
normally unpacks its files in the current directory.
Use the ``\-c'' option to have the program change to a new directory
before invoking the shell.
If the directory does not exist, it will try to create it.
If the directory name starts with a question mark, then
.I unshar
will ask for the directory name before doing anything; this is most useful
with the environment variable UNSHAREDIR.
If the directory name starts with a tilde, then the value of the HOME
environment variable is inserted in place of that character.
For convenience, the ``\-d'' option is a synonym for the ``\-c'' option.
.PP
.I Unshar
normally complains if the input looks like something other than a shar file.
(Among other things, it checks for files that resemble C, and Pascal code).
It can be fooled, however, by nonstandard versions of news, notes, etc.
The ``\-f'' option forces
.I unshar
to try unpacking files, even if they look like something else.
.PP
Depending on how the program is installed,
.I unshar
may or may not try to preserve the header part of file ``foo''
into the name ``foo.hdr'' (if the file is standard input, the name
will be ``UNSHAR.HDR'').
Using the ``\-s'' option forces the program to save the headers, while
using the ``\-n'' option forces it to discard the headers.
The file is appended to, if it already exists, so all headers can be easily
saved in one file.
The name of the file may be given by using the ``\-h'' option; this is
particularly useful when processing more than one file at a time.
.SH ENVIRONMENT
.ta \w'UNSHAREDIR  'u
HOME	Value used if a leading tilde is given in directory name.
.br
UNSHAREDIR	Default value for ``\-c'' option.
.SH SEE ALSO
shar(1).
