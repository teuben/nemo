.\"
.\" Copyright (c) 1992 David I. Bell and Landon Curt Noll
.\" Permission is granted to use, distribute, or modify this source,
.\" provided that this copyright notice remains intact.
.\"
.\" calculator by David I. Bell
.\" man page by Landon Noll
.TH calc 1 "^..^" "22jun91"
.SH NAME
\fBcalc\fP \- arbitrary precision calculator
.SH SYNOPSIS
\fBcalc\fP
[
\fB\-h\fP
] [
\fB\-q\fP
] [
.I calc_cmd
\&.\|.\|.
]
.SH DESCRIPTION
\&
.br
CALC COMMAND LINE
.PP
.TP
\fB \-h\fP
Print a help message.
This option implies \fB \-q\fP.
This is equivalent to the calc command \fBhelp help\fP.
.TP
\fB \-q\fP
Disable the use of the \fB$CALCRC\fP startup library scripts.
.PP
Without \fBcalc_cmd\fPs, \fBcalc\fP operates interactively.
If one or more \fBcalc_cmd\fPs are given on the command line,
\fBcalc\fP will execute them and exit.
.PP
Normally on startup, \fBcalc\fP attempts to execute a collection 
of library scripts.
The environment variable \fB$CALCRC\fP (if non-existent then
a compiled in value) contains a \fB:\fP separated list of
startup library scripts.
No error conditions are produced if these startup library scripts
are not found.
.PP
Filenames are subject to ``~'' expansion (see below).
The environment variable \fB$CALCPATH\fP (if non-existent then
a compiled in value) contains a \fB:\fP separated list of search
directories.
If a file does not begin with \fB/\fP, \fB~\fP or \fB./\fP,
then it is searched for under each directory listed in the \fB$CALCPATH\fP.
It is an error if no such readable file is found.
.PP
For more information use the following calc commands:
.PP
.in 1.0i
help usage
.br
help help
.br
help environment
.in -1.0i
.PP
OVERVIEW
.PP
\fBCalc\fP is arbitrary precision arithmetic system that uses 
a C-like language.
\fBCalc\fP is useful as a calculator, an algorithm prototyped
and as a mathematical research tool.
More importantly, \fBcalc\fP provides one with a machine
independent means of computation.
.PP
A rich set of builtin functions is provided.
A number of library scripts are also provided because they are
useful and to serve as examples of the \fBcalc\fP language.
.PP
One may further extend \fBcalc\fP permits further thru to
use of calc library scripts.
Written in the same C-like language, library scripts may be
read in and executed during a \fBcalc\fP session.
.PP
Internally calc represents numeric values as fractions reduced to their
lowest terms.
The numerators and denominators of these factions may grow to
arbitrarily large values.
Numeric values read in are automatically converted into rationals.
The user need not be aware of this internal representation.
.PP
For more information use the following calc commands:
.PP
.in 1.0i
help intro
.br
help builtin
.br
help stdlib
.br
help define
.br
show builtins
.br
show functions
.in -1.0i
.PP
DATA TYPES
.PP
Fundamental builtin data types include integers, real numbers, 
rational numbers, complex numbers and strings.
.PP
By use of an object, one may define an arbitrarily complex
data types.
One may define how such objects behave a wide range of
operations such as addition, subtraction,
multiplication, division, negation, squaring, modulus,
rounding, exponentiation, equality, comparison, printing
and so on.
.PP
For more information use the following calc commands:
.PP
.in 1.0i
help types
.br
help obj
.br
show objfuncs
.in -1.0i
.PP
VARIABLES
.PP
Variables in \fBcalc\fP are typeless.
In other words, the fundamental type of a variable is determined by its content.
Before variable is assigned a value is of type ``null''.
.PP
The scope of a variable may be global, or only a local to a procedure.
Values may be grouped together in a matrix, or into a
a list that permits stack and queue style operations.
.PP
For more information use the following calc commands:
.PP
.in 1.0i
help variable
.br
help mat
.br
help list
.br
show globals
.in -1.0i
.PP
INPUT/OUTPUT
.PP
A leading ``0x'' implies a hexadecimal value,
a leading ``0b'' implies a binary value,
and a ``0'' followed by a digit implies an octal value.
Complex numbers are indicated by a trailing ``i'' such as in ``3+4i''.
Strings may be delimited by either a pair of single or double quotes.
By default, \fBcalc\fP prints values as if they were floating point numbers.
One may change the default to print values in a number of modes
including fractions, integers and exponentials.
.PP
A number of stdio-like file I/O operations are provided.
One may open, read, write, seek and close files.
Filenames are subject to ``\~'' expansion to home directories
in a way similar to that of the Korn or C-Shell.
.PP
For example:
.PP
.in 1.0i
~/.calcrc
.br
~chongo/lib/fft_multiply.cal
.in -1.0i
.PP
For more information use the following calc command:
.PP
.in 1.0i
help file
.in -1.0i
.PP
CALC LANGUAGE
.PP
The \fBcalc\fP language is a C-like language.
The language includes commands such as variable declarations, 
expressions, tests, labels, loops, file operations, function calls.
These commands are very similar to their counterparts in C.
.PP
The language also include a number of commands particular
to \fBcalc\fP itself.
These include commands such as function definition, help, 
reading in library scripts, dump files to a file, error notification, 
configuration control and status.
.PP
For more information use the following calc command:
.PP
.in 1.0i
help command
.br
help statement
.br
help expression
.br
help operator
.br
help config
.in -1.0i
.PP
.SH FILES
.PD 0
.TP 20
${LIBDIR}
library scripts shipped with calc
.br
.sp
.TP 20
${LIBDIR}/help
help files
.br
.sp
Typically ${LIBDIR} is /usr/local/lib/calc
.sp
.SH CREDIT
Written by David I. Bell.
.sp
Thanks for suggestions and encouragement from Peter Miller,
Neil Justusson, and Landon Noll.
.sp
Portions of this program are derived from an earlier set of
public domain arbitrarily precision routines which was posted
to the net around 1984.  By now, there is almost no recognizable 
code left from that original source.
.sp
Most of this source and binary is:
.sp
.PP
.in 1.0i
Copyright (c) 1992 David I. Bell
.sp
.in -1.0i
.PP
Some files are a copyrighted David I. Bell and Landon Noll.
.sp
Permission is granted to use, distribute, or modify this source,
provided that this copyright notice remains intact.
.sp
Send calc comments, suggestions, bug fixes, enhancements
and interesting calc scripts that you would like you see included 
in future distributions to:
.sp
.PP
.in 1.0i
dbell@pdact.pd.necisa.oz.au\ \ and\ \ chongo@toad.com
.sp
.in -1.0i
.PP
.sp
Enjoy!
